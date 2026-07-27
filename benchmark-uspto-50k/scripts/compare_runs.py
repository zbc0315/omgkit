"""比两次基准跑的逐条结果,看一处改动到底动了什么。

# 为什么不能只看总命中率

命中率是个标量,升了降了都盖着细节:一处修复可能同时救回 60 条、又碰坏 3 条,
净值 +57 看着漂亮,那 3 条却是新缺陷。所以逐条比,并把四种迁移分开数:

    未中 → 命中   救回来的
    命中 → 未中   **碰坏的,必须逐条解释**
    净化失败数变化 产物可用性,与命中率无关但同样要紧

# 与"预期影响范围"对账

改动若有可预先算出的影响范围(比如"只动某个条件成立的那些反应"),就把实际
变化的集合与预期集合对一遍。落在预期之外的每一条都是没想到的副作用。
"""

import argparse
import json
import os
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

DIRS = ("fwd", "retro")
ENGINES = ("omgkit", "rdkit")


def load(path):
    """按 (row, engine, direction) 收 hit / n_bad。"""
    out = {}
    with open(path) as fh:
        for line in fh:
            r = json.loads(line)
            if "err" in r:
                continue
            for e in ENGINES:
                for d in DIRS:
                    v = r.get(f"{e}_{d}")
                    if v and "err" not in v:
                        out[(r["row"], e, d)] = (v["hit"], v.get("n_bad", 0), v.get("n_uniq", 0))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--before", default=os.path.join(ROOT, "results", "bench.jsonl.before-bondfix"))
    ap.add_argument("--after", default=os.path.join(ROOT, "results", "bench.jsonl"))
    ap.add_argument(
        "--expect",
        default=os.path.join(ROOT, "results", "bond_ownership.json"),
        help="预期影响范围;给了就对账",
    )
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "compare_runs.json"))
    args = ap.parse_args()

    a, b = load(args.before), load(args.after)
    common = a.keys() & b.keys()
    print(f"前 {len(a)} 组,后 {len(b)} 组,可比 {len(common)} 组\n")

    moves = Counter()
    changed = {"gained": [], "lost": [], "bad_changed": []}
    for k in sorted(common):
        (h0, nb0, u0), (h1, nb1, u1) = a[k], b[k]
        eng, d = k[1], k[2]
        if h0 != h1:
            tag = "gained" if h1 else "lost"
            moves[f"{eng}_{d}_{tag}"] += 1
            changed[tag].append({"row": k[0], "engine": eng, "dir": d})
        if nb0 != nb1:
            moves[f"{eng}_{d}_bad_delta"] += nb1 - nb0
            changed["bad_changed"].append(
                {"row": k[0], "engine": eng, "dir": d, "before": nb0, "after": nb1}
            )

    print("== 命中迁移 ==")
    if not moves:
        print("  没有任何变化")
    for k, v in sorted(moves.items()):
        print(f"  {k:26s} {v:+d}" if "delta" in k else f"  {k:26s} {v}")

    # 净化失败的反应条数(不是产物数),两次各多少
    for tag, src in (("前", a), ("后", b)):
        c = Counter()
        for (row, e, d), (h, nb, u) in src.items():
            if nb:
                c[f"{e}_{d}"] += 1
        print(f"\n== {tag}:产物净化不过的反应条数 ==")
        for k, v in sorted(c.items()):
            print(f"  {k:16s} {v}")

    if args.expect and os.path.exists(args.expect):
        exp = json.load(open(args.expect))
        want = {
            (p["row"], p["dir"]) for p in exp["picks"]["path_written"]
        }
        n_expected = (
            exp["tally"].get("fwd_path_written", 0) + exp["tally"].get("retro_path_written", 0)
        )
        actual = {(c["row"], c["dir"]) for c in changed["gained"] + changed["lost"]}
        actual |= {(c["row"], c["dir"]) for c in changed["bad_changed"]}
        print(f"\n== 与预期影响范围对账 ==")
        print(f"  预期受影响(路径写法):{n_expected} 条(样本里留了 {len(want)} 条)")
        print(f"  实际有变化:          {len(actual)} 条")
        outside = actual - want
        if want:
            print(f"  实际变化里落在样本外的:{len(outside)} 条")
            if len(outside) <= 20:
                for x in sorted(outside):
                    print(f"      {x}")

    with open(args.out, "w") as fh:
        json.dump({"moves": dict(moves), "changed": changed}, fh, ensure_ascii=False, indent=2)
    print(f"\n已写出 {args.out}")


if __name__ == "__main__":
    main()
