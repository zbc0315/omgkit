"""把 RDKit 的界矩阵导出成 JSONL,给 `omgkit-conf` 的三角光滑化当**外部判官**。

每个分子导:`raw`(未光滑界矩阵)、`smoothed`(RDKit 光滑过的)、
连接表(`z` + `bonds`)、以及一个 **MMFF 优化过的真实构象**(`coords`)。
判据的用法是:把 `raw` 喂给我们的 `triangle_smooth`,结果与 `smoothed` 逐元素比。

**为什么要外部判官**:自己写的光滑化拿自己的性质去验(比如"再光滑一次不变"),
两边同错就同绿。RDKit 不知道我们怎么写的,它只认数。

约定与 RDKit 一致:**上三角是上限、下三角是下限**。

用法:

    python3 harness/dump_bounds.py harness/corpus/large.smi out.jsonl [分子数上限]

注意:装在系统 python 里的 rdkit 可能是旧版且与 numpy 2 冲突。
本仓用的是 2025.09.2(见 harness/embed_ablation.py 顶部那段版本说明)。
"""

import json
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdDistGeom


def _order_tag(b) -> int:
    """键级:1/2/3,芳香记 4。"""
    from rdkit.Chem import BondType

    return {BondType.SINGLE: 1, BondType.DOUBLE: 2, BondType.TRIPLE: 3,
            BondType.AROMATIC: 4}.get(b.GetBondType(), 1)

RDLogger.DisableLog("rdApp.*")


def main() -> int:
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    corpus, out = sys.argv[1], sys.argv[2]
    limit = int(sys.argv[3]) if len(sys.argv) > 3 else 10**9

    n_ok = 0
    n_skip = 0
    with open(out, "w", encoding="utf-8") as fh:
        for line in open(corpus, encoding="utf-8"):
            smi = line.split("\t")[0].strip()
            if not smi:
                continue
            if n_ok >= limit:
                break
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                n_skip += 1
                continue
            mol = Chem.AddHs(mol)
            try:
                raw = rdDistGeom.GetMoleculeBoundsMatrix(mol, doTriangleSmoothing=False)
                sm = rdDistGeom.GetMoleculeBoundsMatrix(mol, doTriangleSmoothing=True)
            except Exception:  # noqa: BLE001 —— RDKit 会抛 C++ 异常
                n_skip += 1
                continue
            # **连接表也要导。** RDKit 的 AddHs 顺序与 omgkit 补氢的顺序不保证一致,
            # 两边各自解析 SMILES 就会错位。判官按这张表直接建分子,下标天生对齐。
            znum = [a.GetAtomicNum() for a in mol.GetAtoms()]
            # **形式电荷必须一起导。** 少了它,[NH3+] 到了对面就是个带四根键的
            # 中性氮,价键检查当场判死 —— 实测 400 个分子里 201 个因此建不出来。
            chg = [a.GetFormalCharge() for a in mol.GetAtoms()]
            rad = [a.GetNumRadicalElectrons() for a in mol.GetAtoms()]
            # **立体标记也要导。** 界矩阵靠 stereo + stereo_atoms 把 1-4 的
            # 顺反析取解掉;不导出来,判官那边这一支永远走不到,而它看不出区别 ——
            # 只是界更松。
            bonds = []
            for b in mol.GetBonds():
                st = int(b.GetStereo())
                sa = list(b.GetStereoAtoms())
                bonds.append(
                    [
                        b.GetBeginAtomIdx(),
                        b.GetEndAtomIdx(),
                        _order_tag(b),
                        st,
                        sa[0] if len(sa) == 2 else -1,
                        sa[1] if len(sa) == 2 else -1,
                    ]
                )
            # 一个真实构象:ETKDG 嵌入 + MMFF 优化。判据一("界必须包住真实构象")用它。
            coords = None
            probe = Chem.Mol(mol)
            p = rdDistGeom.ETKDGv3()
            p.randomSeed = 0xF00D
            if rdDistGeom.EmbedMolecule(probe, p) >= 0:
                try:
                    AllChem.MMFFOptimizeMolecule(probe, maxIters=2000)
                    conf = probe.GetConformer()
                    coords = [
                        [conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
                         conf.GetAtomPosition(i).z]
                        for i in range(probe.GetNumAtoms())
                    ]
                except Exception:  # noqa: BLE001
                    coords = None
            fh.write(
                json.dumps(
                    {
                        "smiles": smi,
                        "n": int(raw.shape[0]),
                        "z": znum,
                        "charge": chg,
                        "radical": rad,
                        "bonds": bonds,
                        "coords": coords,
                        "raw": [float(x) for x in raw.ravel()],
                        "smoothed": [float(x) for x in sm.ravel()],
                    }
                )
                + "\n"
            )
            n_ok += 1
    print(f"导出 {n_ok} 个分子,跳过 {n_skip} 个 → {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
