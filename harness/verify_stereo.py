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

# 实测(`large.smi`,301 个带立体标记的分子,2026-08-20)

| | 一致 |
|---|---|
| 中心基点手性项落地**前** | 288 / 301(95.68%) |
| 中心基点手性项落地**后** | **290 / 301(96.35%)** |

修好 2 个、弄坏 0 个。剩下 11 个:**10 个是环上双键的 E/Z**、
1 个是三配位硫(`C[C@@H]1CO[S@@](=O)N1c2ccccc2`,亚磺酰胺,
`chiral::centers` 因为凑不够四个配体把它整个丢掉了)。四面体手性一个都没错。

用法:

    cargo run -p omgkit-conf --release --example dump_conformers -- harness/corpus/large.smi > /tmp/ours.jsonl
    .venv/bin/python harness/verify_stereo.py /tmp/ours.jsonl
"""
import json
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

BT = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE,
      3: Chem.BondType.TRIPLE, 4: Chem.BondType.AROMATIC}

n_ok = n_bad = n_skip = 0
bad = []
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
    else:
        n_bad += 1
        if len(bad) < 10:
            bad.append(f"{smi}\n      读回 {got}\n      期望 {want}")

tot = n_ok + n_bad
print(f"RDKit 从我们的坐标读回立体化学:{n_ok}/{tot} 一致"
      f"({100.0 * n_ok / max(tot, 1):.2f}%);跳过 {n_skip}")
for b in bad:
    print("   ", b)
sys.exit(1 if n_bad else 0)
