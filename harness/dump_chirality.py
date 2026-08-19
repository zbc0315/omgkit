"""导**手性中心的真值**,给 `omgkit-conf` 的手性抽取当外部判官。

# 判什么

`chiral.rs` 里有两件事,难度完全不同:

- **几何**(有符号四点体积、镜像变号):本地判据就能钉死。
- **抽中心**(哪些原子是手性中心、四个配体按什么槽位排):**必须外部验**。
  槽位排错的话错法是**整批一致**的,于是"符号正确率"会是 0% 或 100% ——
  两个数看起来都像"约定定死了",实际一个全对一个全错。

所以这里导的是**真值**:每个手性中心,四个配体按 RDKit 的键序,以及
有符号体积在**真实三维构象**上的实际符号。omgkit 那边按同一张连接表建分子,
自己抽一遍中心、自己预测符号,与这个真值逐个比。

# 为什么另起一个文件,不并进 rdkit_bounds.jsonl

那份基准是 **RDKit 2025.09.2** 导的,而手边能跑的是 2022.09.5 ——
实测两版的界矩阵逐元素差到 **0.35 Å**,拿 2022 重导会把三条判官的参照悄悄换掉。
手性与界矩阵无关(它只依赖 SMILES 的立体标记与真实坐标),所以单独一个文件、
单独一个版本,互不牵连。

顺带实测过:2022 的立体标记与 2025 导出的坐标**逐个吻合**(39/39),
所以标记这一层在两版之间是稳的。

约定(实测得来,不是猜的):`@`(CCW)→ 有符号体积**为负**;`@@`(CW)→ **为正**,
配体取 `GetNeighbors()` 的顺序。这与 `omgkit-depict` 的 `read_chirality` 一致。

用法:

    .venv/bin/python harness/dump_chirality.py harness/corpus/large.smi out.jsonl 400
"""

import json
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdDistGeom

RDLogger.DisableLog("rdApp.*")

SEED = 0xF00D
TETRA = (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW)


def signed_volume(p0, p1, p2, p3):
    """det[p1−p0, p2−p0, p3−p0]。与 `chiral.rs::signed_volume` 同一个式子。"""
    def d(p):
        return [p[k] - p0[k] for k in range(3)]

    a, b, c = d(p1), d(p2), d(p3)
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def main() -> int:
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    corpus, out = sys.argv[1], sys.argv[2]
    limit = int(sys.argv[3]) if len(sys.argv) > 3 else 10**9

    n_ok = n_skip = n_centers = 0
    n_disagree = 0
    with open(out, "w", encoding="utf-8") as fh:
        for line in open(corpus, encoding="utf-8"):
            smi = line.split("\t")[0].strip()
            if not smi or n_ok >= limit:
                continue
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                n_skip += 1
                continue
            # 只收**带四配位手性中心**的分子 —— 没有中心的分子对这个判据没信息
            if not any(a.GetChiralTag() in TETRA for a in mol.GetAtoms()):
                continue
            mol = Chem.AddHs(mol)
            try:
                p = rdDistGeom.ETKDGv3()
                p.randomSeed = SEED
                if rdDistGeom.EmbedMolecule(mol, p) < 0:
                    n_skip += 1
                    continue
                AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
            except Exception:  # noqa: BLE001
                n_skip += 1
                continue
            conf = mol.GetConformer()
            coords = [
                [conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z]
                for i in range(mol.GetNumAtoms())
            ]

            centers = []
            for a in mol.GetAtoms():
                t = a.GetChiralTag()
                if t not in TETRA:
                    continue
                nb = [x.GetIdx() for x in a.GetNeighbors()]
                if len(nb) != 4:
                    continue
                v = signed_volume(*[coords[i] for i in nb])
                # 真值取**真实构象算出来的号**,不取标记推的号 ——
                # 标记推的那个正是待验的东西,拿它当真值就成了自证。
                actual = -1 if v < 0 else 1
                expected = -1 if t == Chem.ChiralType.CHI_TETRAHEDRAL_CCW else 1
                if actual != expected:
                    n_disagree += 1
                centers.append(
                    {
                        "atom": a.GetIdx(),
                        "nbrs": nb,
                        "tag": 2 if t == Chem.ChiralType.CHI_TETRAHEDRAL_CCW else 1,
                        "sign": actual,
                        "vol": v,
                    }
                )
            if not centers:
                continue
            n_centers += len(centers)

            bonds = []
            for b in mol.GetBonds():
                order = {
                    Chem.BondType.SINGLE: 1,
                    Chem.BondType.DOUBLE: 2,
                    Chem.BondType.TRIPLE: 3,
                    Chem.BondType.AROMATIC: 4,
                }.get(b.GetBondType(), 1)
                bonds.append([b.GetBeginAtomIdx(), b.GetEndAtomIdx(), order])
            fh.write(
                json.dumps(
                    {
                        "smiles": smi,
                        "n": mol.GetNumAtoms(),
                        "z": [a.GetAtomicNum() for a in mol.GetAtoms()],
                        "charge": [a.GetFormalCharge() for a in mol.GetAtoms()],
                        "bonds": bonds,
                        "coords": coords,
                        "centers": centers,
                    }
                )
                + "\n"
            )
            n_ok += 1

    print(f"导出 {n_ok} 个分子、{n_centers} 个四配位手性中心 → {out}(跳过 {n_skip})")
    # 这一行是**版本自检**:标记推的号与真实构象算的号必须处处一致。
    # 不一致说明立体信息在嵌入过程中丢了,那这份真值就不能用。
    print(f"  标记与真实构象不一致的中心:{n_disagree}(必须是 0)")
    return 1 if n_disagree else 0


if __name__ == "__main__":
    sys.exit(main())
