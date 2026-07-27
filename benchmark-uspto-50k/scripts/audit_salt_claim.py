"""把"盐有多普遍"这个数字拆开:小离子里有多少是**真反应物**。

# 为什么要拆

`measure_intra_salt.py` 数的是"反应物侧含小离子或带电片段"的记录数(3625 条、
7.2%),拿它支撑"盐不是边角情形"。但 USPTO 的原始串把所有东西都用 `.` 连起来,
文本上分不清

    真反应物      [NH4+] 当氮源、Wittig 的鏻盐、格氏试剂
    盐的另一半    Cl⁻、Na⁺ 这类只是配平电荷的反离子

两者对"契约缺口"这个论证的分量完全不同:前者是模板本来就该匹配的东西,后者
才是"模板碰不到、又不知道该归谁"的那一档。混在一起数,论证就虚了。

# 判据

`templates.jsonl` 的 `reactants` 是**参与反应**的分子(有原子映射号进产物的)。
原始串里的某个小离子/带电片段,在不在这个集合里,就是判据:

    在      真反应物
    不在    旁观,即真正的反离子那一档
"""

import argparse
import csv
import json
import os
from collections import Counter

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def stripped(smi):
    """去掉原子映射号,给出**规范** SMILES。

    最后那次"写出→读回→再写出"不能省。`MolToSmiles` 在"带映射号解析、事后清掉
    映射号"的分子上给出的**不是**规范式 —— 映射号参与了解析时的定序,清号并不会
    让它重算。实测 `[CH2:10][C@H:11]2…` 这类串清号后直接写出得到 `[C@@H]…`,
    重新读一遍再写才收敛到 `[C@H]…`。

    这一步漏掉,同一个分子会拿到两串,membership 判据就会把真反应物误判成旁观。
    """
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None, None
    for a in m.GetAtoms():
        a.SetAtomMapNum(0)
    again = Chem.MolFromSmiles(Chem.MolToSmiles(m))
    if again is None:
        return None, None
    return Chem.MolToSmiles(again), m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default=os.path.join(ROOT, "data", "uspto50k_raw.csv"))
    ap.add_argument("--tpl", default=os.path.join(ROOT, "data", "templates.jsonl"))
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "salt_claim.json"))
    args = ap.parse_args()

    tpl = {}
    with open(args.tpl) as fh:
        for line in fh:
            r = json.loads(line)
            tpl[r["row"]] = r

    c = Counter()
    examples = {"参与": [], "旁观": []}
    with open(args.raw) as fh:
        for i, row in enumerate(csv.DictReader(fh)):
            rec = tpl.get(i)
            if rec is None:
                continue
            rxn = row["rxn_smiles"]
            if rxn.count(">") != 2:
                continue
            c["记录"] += 1
            want = set()
            for s in rec["reactants"]:
                try:
                    want.add(Chem.CanonSmiles(s))
                except Exception:
                    pass
            saw_part = saw_spec = False
            for frag in rxn.split(">")[0].split("."):
                s, m = stripped(frag)
                if m is None:
                    continue
                small = m.GetNumHeavyAtoms() <= 2
                charged = sum(a.GetFormalCharge() for a in m.GetAtoms()) != 0
                if not (small or charged):
                    continue
                if s in want:
                    saw_part = True
                    if len(examples["参与"]) < 12:
                        examples["参与"].append({"row": i, "id": rec["id"], "frag": s})
                else:
                    saw_spec = True
                    if len(examples["旁观"]) < 12:
                        examples["旁观"].append({"row": i, "id": rec["id"], "frag": s})
            if saw_part:
                c["含小离子/带电片段·是参与反应的"] += 1
            if saw_spec:
                c["含小离子/带电片段·旁观(真反离子那一档)"] += 1
            if saw_part or saw_spec:
                c["含小离子/带电片段·合计"] += 1

    n = c["记录"]
    print(f"记录 {n} 条\n")
    for k, v in c.items():
        if k == "记录":
            continue
        print(f"  {k:40s} {v:6d}  ({100 * v / n:.1f}%)")
    print("\n== 举例:被判为参与反应的 ==")
    for e in examples["参与"][:6]:
        print(f"  行 {e['row']:6d} {e['id']:18s} {e['frag']}")
    print("\n== 举例:被判为旁观的 ==")
    for e in examples["旁观"][:6]:
        print(f"  行 {e['row']:6d} {e['id']:18s} {e['frag']}")

    with open(args.out, "w") as fh:
        json.dump({"tally": dict(c), "examples": examples}, fh, ensure_ascii=False, indent=2)
    print(f"\n已写出 {args.out}")


if __name__ == "__main__":
    main()
