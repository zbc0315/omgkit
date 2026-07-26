"""给 omgkit 未命中的每一条定责:是引擎的问题,还是模板/记录的问题。

# 为什么能定责

抽模板用的就是这条反应本身,所以"把模板用回原反应"应当还原出原记录 ——
这是 rdchiral 抽取的设计意图。还原不出来时只有两种可能:

1. **模板不够**。抽取只保留反应中心加一圈环境,真实反应里有些改变没被写进
   模板(试剂参与、多处反应、互变异构写法不同…)。这时任何忠实执行模板的
   实现都还原不出来,不是实现的问题。
2. **执行有误**。模板说得清清楚楚,实现没做到。

分开的判据:**另一个独立实现给出什么**。两个实现给出同一个错答案,几乎必然是
模板不够;给出不同答案,才轮到查实现。

# 分档

engine-stereo   骨架对、立体不同,而且两个实现给的不一样 → 查实现
engine-skeleton 骨架就不同,两个实现给的不一样      → 查实现
record-stereo   骨架对、只差立体,两个实现给的一样  → 记录/模板的立体信息不足
record-skeleton 骨架不同,两个实现给的一样          → 模板还原不出这条反应
truth-partial   真值是预测的**子集/超集**(分子数不同)→ 记录里混了不参与的分子
"""

import argparse
import json
import os
from collections import Counter

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def skeleton(smi):
    m = Chem.MolFromSmiles(smi)
    return None if m is None else Chem.MolToSmiles(m, isomericSmiles=False)


def frags(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    return frozenset(
        Chem.MolToSmiles(f, isomericSmiles=False) for f in Chem.GetMolFrags(m, asMols=True)
    )


def attribute(rec):
    truth = rec["truth"]
    og = rec["omgkit_preds"]
    rd = rec["rdkit_preds"]
    agree = set(og) == set(rd)
    ts = skeleton(truth)
    tf = frags(truth)

    if not og:
        return ("no-output-agree" if agree else "no-output-engine"), None

    # 骨架层面对得上吗
    hit_skel = [p for p in og if skeleton(p) == ts]
    if hit_skel:
        kind = "stereo"
        closest = hit_skel[0]
    else:
        # 分子拆分不同?比较片段集合
        part = [p for p in og if tf is not None and frags(p) is not None and (frags(p) < tf or frags(p) > tf)]
        if part:
            return ("truth-partial-agree" if agree else "truth-partial-engine"), part[0]
        kind = "skeleton"
        closest = og[0]
    return (f"{'record' if agree else 'engine'}-{kind}"), closest


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default=os.path.join(ROOT, "results", "bench.diff.jsonl"))
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "misses_attributed.jsonl"))
    ap.add_argument("--show", type=int, default=3)
    args = ap.parse_args()

    tally = Counter()
    samples = {}
    with open(args.src) as fh, open(args.out, "w") as out:
        for line in fh:
            rec = json.loads(line)
            if rec["winner"] == "omgkit":
                continue  # omgkit 中了,不在本次归因范围
            bucket, closest = attribute(rec)
            key = (rec["direction"], rec["winner"], bucket)
            tally[key] += 1
            samples.setdefault(key, []).append((rec, closest))
            out.write(json.dumps({**rec, "bucket": bucket, "closest": closest}) + "\n")

    print("== 归因 ==  (方向, 谁命中, 档次) -> 条数")
    for k in sorted(tally, key=lambda x: -tally[x]):
        print(f"  {k[0]:6s} {k[1]:6s} {k[2]:22s} {tally[k]}")
    print(f"\n合计 {sum(tally.values())}")

    for k in sorted(tally, key=lambda x: -tally[x]):
        print(f"\n---- {k} ({tally[k]} 条) ----")
        for rec, closest in samples[k][: args.show]:
            print(f"  {rec['id']}")
            print(f"    真值   {rec['truth']}")
            print(f"    omgkit {closest}")
            print(f"    rdkit  {(rec['rdkit_preds'] or ['<空>'])[0]}")
            print(f"    模板   {(rec['retro'] if rec['direction'] == 'retro' else rec['fwd'])[:170]}")
    print(f"\n已写入 {args.out}")


if __name__ == "__main__":
    main()
