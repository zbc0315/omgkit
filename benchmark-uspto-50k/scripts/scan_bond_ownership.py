"""语料上"两端都被匹配、模板却没匹配到"的键有多普遍,以及哪些真模板能当判据。

# 为什么要扫

`carry_over` 原来的判据是"两端都被匹配就不搬"。子结构匹配只要求模板的每根键
在底物里找得到,**不**要求这些原子之间没有别的键 —— 模板把环写成开链路径时,
环闭合的那根键两端都被匹配,模板却从没看见它。按旧判据它被删掉,环撕开。

判据必须建在**真模板**上:自己造的模板可能根本不是它看上去的意思。所以这里
从语料里挑,并把出现频次一并数出来。

# 分三档

  路径写法   底物里有一根键,两端都被反应物模板匹配到,而模板没匹配它
             —— 旧判据会删掉它
  产物补写   产物模板在两个映射原子之间写了一根键,而反应物模板没匹配它,
             且底物里**已经有**这根键 —— 只按"模板没匹配就搬"会加两遍
  真断键     反应物模板匹配到、产物模板不写 —— 这根键必须删,新判据不能放过

三档都要有真模板,判据才算不缺角。
"""

import argparse
import json
import os
from collections import Counter

from rdkit import Chem, RDLogger
from rdkit.Chem import rdChemReactions

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def pair(a, b):
    return (a, b) if a <= b else (b, a)


def template_bond_map(rt):
    """反应物模板里"映射号对 → 有键"。"""
    out = set()
    for b in rt.GetBonds():
        m1 = rt.GetAtomWithIdx(b.GetBeginAtomIdx()).GetAtomMapNum()
        m2 = rt.GetAtomWithIdx(b.GetEndAtomIdx()).GetAtomMapNum()
        if m1 and m2:
            out.add(pair(m1, m2))
    return out


def analyse(smarts, substrates):
    """对一条模板 + 一组底物,数出三档各命中多少。"""
    try:
        rx = rdChemReactions.ReactionFromSmarts(smarts)
        rx.Initialize()
    except Exception:
        return None
    n = rx.GetNumReactantTemplates()
    mols = [Chem.MolFromSmiles(s) for s in substrates]
    if any(m is None for m in mols):
        return None

    rt_bonds = set()
    for ti in range(n):
        rt_bonds |= template_bond_map(rx.GetReactantTemplate(ti))
    pt_bonds = set()
    kept_maps = set()
    for ti in range(rx.GetNumProductTemplates()):
        pt = rx.GetProductTemplate(ti)
        pt_bonds |= template_bond_map(pt)
        kept_maps |= {a.GetAtomMapNum() for a in pt.GetAtoms() if a.GetAtomMapNum()}

    res = Counter()
    # 真断键:反应物模板有、产物模板没有的映射号对
    res["broken"] = len(rt_bonds - pt_bonds)
    # 产物补写:产物模板有、反应物模板没有的映射号对
    added = pt_bonds - rt_bonds

    # 路径写法要看底物:枚举匹配,找"两端都匹配、模板没匹配"的键。
    #
    # 必须看**每一处**匹配,不能只看第一处 —— 同一条模板在同一个底物上往往
    # 有多处匹配,带额外键的那处不一定排在前面(实测嘌呤那条就排在后面)。
    # 每条反应最多记一次,记的是"最多有几根这样的键"。
    worst_extra = 0
    for ti in range(n):
        rt = rx.GetReactantTemplate(ti)
        by_map_q = {
            rt.GetAtomWithIdx(i).GetAtomMapNum(): i
            for i in range(rt.GetNumAtoms())
            if rt.GetAtomWithIdx(i).GetAtomMapNum()
        }
        for mol in mols:
            for m in mol.GetSubstructMatches(rt, uniquify=False, maxMatches=500):
                # 只算**两端都活到产物里**的那些原子。模板里没有映射号的原子
                # 是要删掉的,它们之间多出来的键无论如何都不该搬 —— 那一档由
                # 另一条判据("匹配到却不在产物里,遍历到此为止")挡下,新旧
                # 两版行为一致,算进来会把影响范围虚报好几倍。
                alive = {
                    m[i]
                    for i in range(rt.GetNumAtoms())
                    if rt.GetAtomWithIdx(i).GetAtomMapNum() in kept_maps
                }
                tb = {pair(m[b.GetBeginAtomIdx()], m[b.GetEndAtomIdx()]) for b in rt.GetBonds()}
                extra = sum(
                    1
                    for b in mol.GetBonds()
                    for p in [pair(b.GetBeginAtomIdx(), b.GetEndAtomIdx())]
                    if p[0] in alive and p[1] in alive and p not in tb
                )
                worst_extra = max(worst_extra, extra)
                # 产物补写,且底物**已经有**这根键
                for m1, m2 in added:
                    if m1 in by_map_q and m2 in by_map_q:
                        a, b2 = m[by_map_q[m1]], m[by_map_q[m2]]
                        if mol.GetBondBetweenAtoms(a, b2) is not None:
                            res["added_already_present"] = 1
    if worst_extra:
        res["path_written"] = 1
        res["path_written_bonds"] = worst_extra
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tpl", default=os.path.join(ROOT, "data", "templates.jsonl"))
    ap.add_argument("--limit", type=int, default=0, help="0 = 全量")
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "bond_ownership.json"))
    args = ap.parse_args()

    tally = Counter()
    picks = {"path_written": [], "added_already_present": [], "broken": []}
    n = 0
    with open(args.tpl) as fh:
        for line in fh:
            r = json.loads(line)
            if args.limit and n >= args.limit:
                break
            n += 1
            for d in ("retro", "fwd"):
                subs = [r["prod"]] if d == "retro" else r["reactants"]
                res = analyse(r[d], subs)
                if res is None:
                    tally[f"{d}_unparsable"] += 1
                    continue
                tally[f"{d}_ok"] += 1
                for k in ("path_written", "added_already_present", "broken"):
                    if res[k]:
                        tally[f"{d}_{k}"] += 1
                        if len(picks[k]) < 1000:
                            picks[k].append(
                                {
                                    "id": r["id"],
                                    "row": r["row"],
                                    "dir": d,
                                    "smarts": r[d],
                                    "substrates": subs,
                                    "count": res[k],
                                    "extra_bonds": res.get("path_written_bonds", 0),
                                }
                            )
                tally[f"{d}_path_bonds"] += res["path_written_bonds"]

    print(f"扫了 {n} 条记录,正逆各一次\n")
    for k, v in sorted(tally.items()):
        print(f"  {k:32s} {v}")
    with open(args.out, "w") as fh:
        json.dump({"n_records": n, "tally": dict(tally), "picks": picks}, fh,
                  ensure_ascii=False, indent=2)
    print(f"\n已写出 {args.out}")


if __name__ == "__main__":
    main()
