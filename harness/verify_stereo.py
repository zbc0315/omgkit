"""把我们交付的坐标交给 RDKit 读回立体化学,逐分子验它满不满足输入 SMILES。

**完全绕开我们自己的任何公式** —— 这是唯一真正外部的立体化学判据。
配套的导出工具是 `cargo run -p omgkit-conf --release --example dump_conformers`。

# 为什么不能直接比规范 SMILES

三维结构**必然**给桥头碳之类定一个构型,而输入 SMILES 往往不写(它不是独立的
立体中心)。头一版直接比规范 SMILES,301 个分子里报 41 个"不一致",逐个看下去
绝大多数是**读回的立体信息比输入更多** —— 那不是错。

要问的是"我们的结构满不满足输入**指定**的每一处立体"。RDKit 的子结构匹配开
`useChirality=True` 正是这个语义:查询里没指定的原子匹配任意,指定了的必须一致。
换成这个口径之后 290/301,而那 41 里多出来的 30 个正是桥头碳那一类。

# 判官够不着的那一档要**单独报**,不能混进失配

`elif not rdkit_can_read_itself(smi)`:失配时先拿 RDKit **自己的**构象走同一条
检验,它自己都读不回来的分子归"判官够不着"。实测这一档全是**三配位磷**
(`C[P@H]C…`、`[P@@]2CCC…`)—— RDKit 2022.09.5 的 `AssignStereochemistryFrom3D`
不给三配位 P 赋手性,连它自己嵌出来的构象都读不回。不分开的话这两条会被读成
"我们把 P 摆错了",而那是没有依据的结论。

与 `threading_oracle` 先校准检测器再量自己是同一条规矩:**判官先证明自己看得见,
它报的 0 才有意义。**

# 实测(`large.smi`,2026-08-20)

| | 覆盖分子 | 一致 | 判官够不着 |
|---|---|---|---|
| 只收带 `@` 的、中心基点手性项落地前 | 301 | 288(95.68%) | — |
| 同上,落地后 | 301 | 290(96.35%) | — |
| 补上双键顺反折算 + 放宽输入(也收 `/` `\`) | 632 | 631(99.84%) | — |
| **补上三配位立体中心 + 判官自校准** | **642** | **640 / 640(100%)** | **2** |

那 2 个是三配位 P。三配位 S(亚砜、亚磺酰胺)**13 个中心**全部通过
(11 是**分子**数,不是中心数 —— 头一版这里写混了)。

# P 那一档不是"验不了",只是这条判据看不见

先前这里写着"P 的绝对约定没有外部验证",那**过于悲观**:
`AssignStereochemistryFrom3D` 读不回三配位 P,但 RDKit 的**嵌入器**认它 ——
所以真值可以从**嵌出来的构象**上算。`harness/dump_chirality.py` 现在就是这么造
`harness/baseline/smoke.lonepair.jsonl` 的(17 个三配位中心,号跨 seed 不稳的剔除),
而那份基准**进了 CI**(`conformer_oracle` 跑 `smoke.lonepair.jsonl` 那一步)。
实测那条闸上 21/21 全对,把三配位的槽位前两个对调当场红 17 个。

(这里原先写着"第 14 道闸"。闸门序号是**手抄的数**,加一道闸就掉队 ——
`gates.sh` 现在连自己的总数都是算出来 + 末尾自查的,散文里更不该写序号。)

**这条判据本身仍然看不见 P** —— 所以那 2 个继续记在"够不着"里,别把它读成"我们错了"。

用法:

    cargo run -p omgkit-conf --release --example dump_conformers -- harness/corpus/large.smi > /tmp/ours.jsonl
    .venv/bin/python harness/verify_stereo.py /tmp/ours.jsonl
"""
import json
import sys

import rdkit
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

BT = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE,
      3: Chem.BondType.TRIPLE, 4: Chem.BondType.AROMATIC}

# **"判官够不着"的分子数上限。**
#
# 自校准是一个**单向过滤器**:它只会把"失配"变成"不计数"。只让判据变绿的东西
# 必须配一道上限闸,否则没人拦得住它 —— 独立审核实测过:把两个分子的交付坐标
# 整体镜像(= 交付对映体),没有这条闸时判官打印 `0/0 一致、判官够不着 2`
# 并**退出 0**,因为那两个分子恰好 `EmbedMolecule` 失败。
#
# 现值 2(都是三配位 P,RDKit 的 `AssignStereochemistryFrom3D` 不给它赋手性)。
# 设 5 是给多试几个 seed 之后的余量;涨上去要当场查,不是调大它。
MAX_BLIND = 5

# 自校准要试几个 seed。**一个 seed 不够** —— 嵌入失败与"读不回"是两回事,
# 前者说明这次没试出来,后者才说明判官够不着。审核实测:642 个分子里有 2 个
# 在单一 seed 下嵌入失败,于是真错也能被吞进"够不着"。
CALIBRATION_SEEDS = (0xF00D, 0xC0FFEE, 0xBEEF, 1, 7)


def rdkit_can_read_itself(smi):
    """RDKit 从**它自己**生成的构象里,能不能把这条 SMILES 的立体读回来?

    答 `False` 就说明这个分子落在判官的能力之外,拿它去判我们的产物没有依据。

    **多试几个 seed 才算数**:任何一个 seed 上读得回来,就说明判官看得见这一档。
    全部 seed 都嵌不出来的分子返回 `False`(判官确实用不上),但那一档同样计进
    `MAX_BLIND` —— 它是"判官够不着",不是"我们对了"。
    """
    ref = Chem.MolFromSmiles(smi)
    for seed in CALIBRATION_SEEDS:
        m = Chem.AddHs(Chem.MolFromSmiles(smi))
        ps = AllChem.ETKDGv3()
        ps.randomSeed = seed
        if AllChem.EmbedMolecule(m, ps) != 0:
            continue
        AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
        Chem.AssignStereochemistryFrom3D(m)
        if Chem.RemoveHs(m).HasSubstructMatch(Chem.RemoveHs(ref), useChirality=True):
            return True
    return False


# **把版本打出来。** CI 钉的是 `harness/requirements.lock`(2025.09.2),而开发机的
# `.venv` 未必是同一个。这条判据两边喂的是同一个 RDKit,版本差异会对消 ——
# 但"哪个版本给出的这个数"不该靠记。
print(f"外部实现:RDKit {rdkit.__version__}")

n_ok = n_bad = n_skip = n_blind = 0
bad = []
blind = []
for line in open(sys.argv[1], encoding="utf-8"):
    r = json.loads(line)
    smi = r["smiles"]
    ref = Chem.MolFromSmiles(smi)
    if ref is None:
        n_skip += 1
        continue
    want = Chem.MolToSmiles(Chem.RemoveHs(ref))

    rw = Chem.RWMol()
    for z, c in zip(r["z"], r["charge"]):
        a = Chem.Atom(int(z))
        a.SetFormalCharge(int(c))
        a.SetNoImplicit(True)
        rw.AddAtom(a)
    for i, j, o in r["bonds"]:
        rw.AddBond(int(i), int(j), BT.get(int(o), Chem.BondType.SINGLE))
    m = rw.GetMol()
    conf = Chem.Conformer(m.GetNumAtoms())
    for k, p in enumerate(r["xyz"]):
        conf.SetAtomPosition(k, tuple(float(v) for v in p))
    m.AddConformer(conf)
    try:
        Chem.SanitizeMol(m)
        Chem.AssignStereochemistryFrom3D(m)
        got = Chem.MolToSmiles(Chem.RemoveHs(m))
    except Exception as e:  # noqa: BLE001
        n_skip += 1
        continue
    # **不能直接比规范 SMILES**:三维结构必然给桥头碳之类定一个构型,而输入
    # SMILES 往往不写(它不是独立的立体中心)。多出来的立体信息不是错。
    #
    # 要问的是"我们的结构满不满足输入**指定**的每一处立体" ——
    # 子结构匹配开 `useChirality` 正是这个语义:查询里没指定的原子匹配任意,
    # 指定了的必须一致。
    ours = Chem.RemoveHs(m)
    if ours.HasSubstructMatch(Chem.RemoveHs(ref), useChirality=True):
        n_ok += 1
    elif not rdkit_can_read_itself(smi):
        # **判官够不着这一档,不是我们摆错。** 先拿 RDKit **自己的**构象走同一条
        # 检验:它自己都读不回来,那这条失配说明的是判官的能力,不是产物的对错。
        #
        # 实测这一档全是**三配位磷**(`C[P@H]C…`、`[P@@]2CCC…`):
        # RDKit 2022.09.5 的 `AssignStereochemistryFrom3D` 不给三配位 P 赋手性。
        # 不分开的话,这两条会被读成"我们把 P 摆错了" —— 而那是没有依据的结论。
        #
        # 这与 `threading_oracle` 先校准检测器再量自己是同一条规矩:
        # **判官先证明自己看得见,它报的 0 才有意义。**
        n_blind += 1
        if len(blind) < 10:
            blind.append(f"{smi}\n      我们读回 {got}\n      RDKit 自己也读不回")
    else:
        n_bad += 1
        if len(bad) < 10:
            bad.append(f"{smi}\n      读回 {got}\n      期望 {want}")

tot = n_ok + n_bad
print(f"RDKit 从我们的坐标读回立体化学:{n_ok}/{tot} 一致"
      f"({100.0 * n_ok / max(tot, 1):.2f}%);跳过 {n_skip}")
for b in bad:
    print("   ", b)
# **判官够不着的那一档单独报,不混进失配。** 混进去的话它会被读成"我们摆错了",
# 而那是没有依据的结论;单独报出来,数量一涨也看得见判官在退化。
print(f"  判官够不着(RDKit 自己也读不回)的分子:{n_blind}(上限 {MAX_BLIND})")
for b in blind:
    print("   ", b)
if n_blind > MAX_BLIND:
    print(f"\n判官够不着的分子涨到 {n_blind} > {MAX_BLIND} —— "
          "自校准是单向过滤器,它只会把失配变成不计数。"
          "涨上去要当场查是哪一档看不见了,不是调大这个数")
sys.exit(1 if (n_bad or n_blind > MAX_BLIND) else 0)
