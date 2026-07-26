"""逐个原子/键地比对真值与 omgkit 的输出,看立体到底差在哪一档。

两个分子的**构成完全相同**(骨架已经比过),所以这里 CIP 可比 —— 跨反应比 CIP
不可靠是因为取代基换了,同一个分子的两种立体写法没有这个问题。

三档:
  lost         真值写了,omgkit 没写   → 信息丢了,是引擎的问题
  extra        omgkit 写了,真值没写   → 记录欠指定,omgkit 反而更全
  contradict   两边都写了但不一样      → 真·矛盾,必须查

原子与双键分开统计 —— 两者的成因完全不同,混在一起会把线索抹平。
"""

import argparse
import json
import os
from collections import Counter

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def prepared(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    return m


def atom_labels(m):
    """原子下标 -> 'R'/'S';未指定的不收。"""
    return {
        a.GetIdx(): a.GetPropsAsDict()["_CIPCode"]
        for a in m.GetAtoms()
        if a.HasProp("_CIPCode")
    }


def bond_labels(m):
    """(端点下标对) -> STEREOE/STEREOZ;未指定的不收。"""
    out = {}
    for b in m.GetBonds():
        if b.GetStereo() != Chem.BondStereo.STEREONONE:
            out[frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))] = str(b.GetStereo())
    return out


def correspondence(a, b):
    """b 的原子下标 -> a 的原子下标。构成相同才有解。"""
    fa = Chem.MolFromSmiles(Chem.MolToSmiles(a, isomericSmiles=False))
    fb = Chem.MolFromSmiles(Chem.MolToSmiles(b, isomericSmiles=False))
    if fa is None or fb is None:
        return None
    # 直接在原分子之间配,忽略手性
    match = b.GetSubstructMatch(a, useChirality=False)
    if not match or len(match) != a.GetNumAtoms():
        return None
    # match[i] 是 a 的第 i 个原子在 b 里的下标
    return {match[i]: i for i in range(len(match))}


def compare(truth_smi, pred_smi):
    t, p = prepared(truth_smi), prepared(pred_smi)
    if t is None or p is None:
        return None
    p2t = correspondence(t, p)
    if p2t is None:
        return None
    ta, pa = atom_labels(t), atom_labels(p)
    tb, pb = bond_labels(t), bond_labels(p)
    pa = {p2t[i]: v for i, v in pa.items() if i in p2t}
    pb = {
        frozenset(p2t[i] for i in k): v
        for k, v in pb.items()
        if all(i in p2t for i in k)
    }
    out = Counter()
    detail = []
    for key, kind in ((set(ta) | set(pa), "atom"), (set(tb) | set(pb), "bond")):
        tsrc, psrc = (ta, pa) if kind == "atom" else (tb, pb)
        for k in key:
            tv, pv = tsrc.get(k), psrc.get(k)
            if tv == pv:
                continue
            if tv is not None and pv is None:
                out[f"{kind}-lost"] += 1
                detail.append((kind, "lost", str(k), tv, pv))
            elif tv is None and pv is not None:
                out[f"{kind}-extra"] += 1
                detail.append((kind, "extra", str(k), tv, pv))
            else:
                out[f"{kind}-contradict"] += 1
                detail.append((kind, "contradict", str(k), tv, pv))
    return out, detail


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default=os.path.join(ROOT, "results", "misses_attributed.jsonl"))
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "stereo_diff.jsonl"))
    ap.add_argument("--show", type=int, default=2)
    args = ap.parse_args()

    tally = Counter()
    combo = Counter()
    samples = {}
    with open(args.src) as fh, open(args.out, "w") as out:
        for line in fh:
            rec = json.loads(line)
            closest = rec.get("closest")
            if closest is None:
                tally["无可比预测"] += 1
                continue
            res = compare(rec["truth"], closest)
            if res is None:
                tally["对应不上"] += 1
                continue
            counts, detail = res
            if not counts:
                tally["查不出差异"] += 1
                continue
            sig = "+".join(f"{k}x{v}" for k, v in sorted(counts.items()))
            kinds = frozenset(k.split("-")[1] for k in counts)
            label = "/".join(sorted(kinds))
            tally[label] += 1
            combo[(rec["direction"], sig)] += 1
            samples.setdefault(label, []).append((rec, detail))
            out.write(json.dumps({**rec, "stereo_diff": dict(counts), "detail": detail}) + "\n")

    print("== 立体差异的性质 ==")
    for k, v in tally.most_common():
        print(f"  {k:22s} {v}")
    print("\n== 最常见的具体组合(前 12)==")
    for (d, sig), v in combo.most_common(12):
        print(f"  {d:6s} {sig:44s} {v}")
    for label, rows in samples.items():
        print(f"\n---- {label} 举例 ----")
        for rec, detail in rows[: args.show]:
            print(f"  {rec['id']} [{rec['direction']}] {rec['bucket']}")
            print(f"    真值   {rec['truth']}")
            print(f"    omgkit {rec['closest']}")
            print(f"    差异   {detail[:4]}")
    print(f"\n已写入 {args.out}")


if __name__ == "__main__":
    main()
