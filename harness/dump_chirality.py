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

# 真值用的是哪个体积 —— 换过一次,原因是老的那个在结构上抓不到东西

有符号体积有两种写法,**不是同一个量**:

- 四配体:`det[l₁−l₀, l₂−l₀, l₃−l₀]`,**完全不看中心原子在哪**。
- 中心基点:`det[l₀−c, l₁−c, l₂−c]`(RDKit 的 `assignChiralTypesFrom3D` 用这个)。

头一版这里用的是**前者**,注释还写着"与 `chiral.rs::signed_volume` 同一个式子" ——
真值与待验实现同一条式子,于是这个判官**在结构上不可能**抓到"中心原子被挤到
配体四面体外面"(伞形翻转)那一档:那时四配体行列式一点变化都没有,而真实
立体化学已经翻了。

这一档**当前没有在发生**(实测交付坐标上 484 个中心,0 个在四面体外、0 个号不符),
换过来是把洞堵上,不是修一个正在漏的洞 —— 但真值用一条对该失效模式天生失明的
式子,这个判官在结构上就守不住那一档,所以还是要换。

现在用**后者**。约定(实测 247 个中心:127/127、120/120):
`@`(CCW)→ 中心基点体积**为正**;`@@`(CW)→ **为负**,配体取 `GetNeighbors()` 的顺序。
注意这与旧的四配体口径**正好反号**(正四面体上 `V_配体 = −4·V_中心`)。

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

# **三配位立体中心**:三根键 + 一对孤对(亚砜/亚磺酰胺的 S、膦的 P……)。
# 与 `omgkit_core::element::has_stereogenic_lone_pair` 同一张表。
LONE_PAIR = {15, 16, 33, 34, 52}  # P S As Se Te

# 三配位那一档的真值要**跨 seed 稳定**才算数。
#
# RDKit 的 `AssignStereochemistryFrom3D` 不读三配位 P,但它的**嵌入器**认 ——
# 所以真值只能从嵌出来的构象上算。可嵌入器不是每次都摆对:实测
# `Cc1ccc(cc1)[S@@](=O)C2=CCCO2` 的 S 在 5 个 seed 下号是 [-1,+1,-1,-1,-1],
# 有一个 seed 摆反了。号不稳的中心**不进基准** —— 拿一个自己都摇摆的数当真值,
# 判据红了也不知道该信谁。
STABILITY_SEEDS = (0xF00D, 0xC0FFEE, 0xBEEF, 1, 7)


def det3(a, b, c):
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def center_volume(coords, center, nbrs):
    """det[l₀−c, l₁−c, l₂−c] —— **以中心原子为基点**。

    先前这里用的是四个配体的行列式 `det[l₁−l₀, l₂−l₀, l₃−l₀]`,注释还写着
    "与 `chiral.rs::signed_volume` 同一个式子" —— 那正是问题所在:
    真值与待验的实现用同一条式子,于是**这个判官在结构上不可能抓到
    "中心原子翻到配体四面体外面"那一档**(四配体行列式对它完全不敏感)。

    实测交付坐标上 484 个中心里 0 个在四面体外、0 个号不符 —— 这一档当前
    没有在发生,换基点是把洞堵上。
    """
    o = coords[center]

    def d(i):
        return [coords[i][k] - o[k] for k in range(3)]

    return det3(d(nbrs[0]), d(nbrs[1]), d(nbrs[2]))


def ligand_volume(coords, nbrs):
    """旧口径,只用来报"两者反号"这件事,不再当真值。"""
    p0 = coords[nbrs[0]]

    def d(i):
        return [coords[i][k] - p0[k] for k in range(3)]

    return det3(d(nbrs[1]), d(nbrs[2]), d(nbrs[3]))


def stable_sign(smi, atom_idx):
    """这个三配位中心的号,在几个 seed 下是不是稳定的?

    不稳就说明 RDKit 的嵌入器自己都没把这个标记摆稳,拿它当真值没有意义。
    """
    seen = set()
    for seed in STABILITY_SEEDS:
        m = Chem.AddHs(Chem.MolFromSmiles(smi))
        p = rdDistGeom.ETKDGv3()
        p.randomSeed = seed
        if rdDistGeom.EmbedMolecule(m, p) < 0:
            continue
        AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
        cf = m.GetConformer()
        pts = [
            [cf.GetAtomPosition(i).x, cf.GetAtomPosition(i).y, cf.GetAtomPosition(i).z]
            for i in range(m.GetNumAtoms())
        ]
        nb = [x.GetIdx() for x in m.GetAtomWithIdx(atom_idx).GetNeighbors()]
        if len(nb) != 3:
            return False
        seen.add(1 if center_volume(pts, atom_idx, nb) > 0 else -1)
    return len(seen) == 1


def main() -> int:
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    corpus, out = sys.argv[1], sys.argv[2]
    limit = int(sys.argv[3]) if len(sys.argv) > 3 else 10**9

    n_ok = n_skip = n_centers = 0
    n_disagree = n_same_sign = n_zero = n_unstable = n_three = 0
    with open(out, "w", encoding="utf-8") as fh:
        for line in open(corpus, encoding="utf-8"):
            # 语料格式:`SMILES<TAB>名字`,**`#` 开头是注释、空行忽略**。
            # 先前这里没认注释,把它们当 SMILES 解析 —— "跳过"那个数里因此
            # 混着注释行,报出来的数是错的(同一个坑在 `feasibility.rs` 里修过一次)。
            line = line.strip()
            if not line or line.startswith("#"):
                continue
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
                three = (
                    len(nb) == 3
                    and a.GetAtomicNum() in LONE_PAIR
                    and a.GetFormalCharge() <= 0
                )
                if len(nb) != 4 and not three:
                    continue
                if three and not stable_sign(smi, a.GetIdx()):
                    n_unstable += 1
                    continue
                v = center_volume(coords, a.GetIdx(), nb)
                # 四配体行列式要四个配体 —— 三配位那一档没有,记 0 并跳过反号自检
                vl = ligand_volume(coords, nb) if len(nb) == 4 else 0.0
                # 真值取**真实构象算出来的号**,不取标记推的号 ——
                # 标记推的那个正是待验的东西,拿它当真值就成了自证。
                # `v == 0` 时 `actual` 会静默取 +1 —— 那是掷硬币,必须计数
                if v == 0.0:
                    n_zero += 1
                actual = -1 if v < 0 else 1
                # 约定:`@`(CCW)→ 中心基点体积**为正**;`@@`(CW)→ **为负**。
                # 这与旧的四配体口径**正好反过来**(V_配体 = −4·V_中心),
                # 因为量的根本不是同一个体积。
                expected = 1 if t == Chem.ChiralType.CHI_TETRAHEDRAL_CCW else -1
                if actual != expected:
                    n_disagree += 1
                if len(nb) == 4 and v * vl > 0:
                    n_same_sign += 1
                centers.append(
                    {
                        "atom": a.GetIdx(),
                        "nbrs": nb,
                        "three_coordinate": three,
                        "tag": 2 if t == Chem.ChiralType.CHI_TETRAHEDRAL_CCW else 1,
                        "sign": actual,
                        "vol": v,
                        "vol_ligand": vl,
                    }
                )
            if not centers:
                continue
            n_centers += len(centers)
            n_three += sum(1 for c in centers if c["three_coordinate"])

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
    # 这一行说明"换了基点"确实换了量:真实构象上两者应当**处处反号**。
    # 同号的个数不是 0,就说明有中心已经翻伞(或者构象本身是废的)。
    #
    # **它必须进退出码。** 只 print 的自检等于没有:坏基准照样落盘,
    # 而下游所有手性判据都以这份基准为真值。
    print(f"  中心基点与四配体行列式**同号**的中心:{n_same_sign}(真实构象上应当是 0)")
    print(f"  中心基点体积恰好为 0 的中心:{n_zero}(必须是 0 —— 号是掷硬币)")
    print(f"  其中**三配位**(孤对)中心:{n_three};号跨 seed 不稳而被剔除的:{n_unstable}")
    return 1 if (n_disagree or n_same_sign or n_zero) else 0


if __name__ == "__main__":
    sys.exit(main())
