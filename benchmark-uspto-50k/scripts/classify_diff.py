"""给两个引擎判得不一样的那些反应分类。

分三档,判据都用 RDKit(裁判统一):

A 记录欠指定  骨架相同,只是**赢家的真值比输家的预测少写了立体**。
              USPTO 的记录里常有产物不写顺反的,这时保住底物几何的那一方
              反而对不上真值 —— 是记录的问题,不是引擎的。
B 立体相矛盾  骨架相同、两边写的立体一样多,但内容不同。这一档才是缺陷线索。
C 骨架就不同  连去掉立体都对不上,与立体无关。

用法:python classify_diff.py [results/bench.diff.jsonl]
"""

import json
import os
import sys
from collections import Counter

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def flat_skeleton(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    return Chem.MolToSmiles(m, isomericSmiles=False)


def stereo_count(smi):
    """写死了的立体标注个数:双键顺反 + 四面体手性。"""
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return -1
    nb = sum(1 for b in m.GetBonds() if b.GetStereo() != Chem.BondStereo.STEREONONE)
    na = sum(
        1
        for a in m.GetAtoms()
        if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
    )
    return nb + na


def classify(rec):
    """返回 (档次, 最接近真值的那个预测)。"""
    truth = rec["truth"]
    loser = "rdkit" if rec["winner"] == "omgkit" else "omgkit"
    preds = rec[f"{loser}_preds"]
    if not preds:
        return "C", None, loser
    ts = flat_skeleton(truth)
    same = [p for p in preds if flat_skeleton(p) == ts]
    if not same:
        return "C", preds[0], loser
    # 骨架相同的预测里,挑立体标注最多的那个来判
    best = max(same, key=stereo_count)
    if stereo_count(truth) < stereo_count(best):
        return "A", best, loser
    return "B", best, loser


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "results", "bench.diff.jsonl")
    recs = [json.loads(l) for l in open(path)]
    tally = Counter()
    buckets = {"A": [], "B": [], "C": []}
    for r in recs:
        cls, pred, loser = classify(r)
        tally[(r["direction"], r["winner"], cls)] += 1
        buckets[cls].append((r, pred, loser))

    print("== 分档 ==  (方向, 命中方, 档次) -> 条数")
    for k in sorted(tally):
        print(f"  {k}  {tally[k]}")

    for cls in ("B", "C", "A"):
        rows = buckets[cls]
        print(f"\n== {cls} 档,共 {len(rows)} 条,列前 6 ==")
        for r, pred, loser in rows[:6]:
            print(f"  {r['id']} [{r['direction']}] 赢家={r['winner']}")
            print(f"    真值      {r['truth']}")
            print(f"    {loser:6s}最近 {pred}")

    out = os.path.join(ROOT, "results", "diff_classified.jsonl")
    with open(out, "w") as fh:
        for cls in ("A", "B", "C"):
            for r, pred, loser in buckets[cls]:
                fh.write(
                    json.dumps(
                        {
                            **r,
                            "bucket": cls,
                            "loser": loser,
                            "loser_closest": pred,
                        }
                    )
                    + "\n"
                )
    print(f"\n已写入 {out}")


if __name__ == "__main__":
    main()
