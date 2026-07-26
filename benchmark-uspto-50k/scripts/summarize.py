"""汇总:命中率、耗时分位、归因分布。输出 results/summary.json 与一段 Markdown。"""

import argparse
import json
import os
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

KEYS = ["omgkit_fwd", "rdkit_fwd", "omgkit_retro", "rdkit_retro"]


def pct(v, n):
    return f"{100.0 * v / n:.2f}%" if n else "—"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench", default=os.path.join(ROOT, "results", "bench.jsonl"))
    ap.add_argument("--verdicts", default=os.path.join(ROOT, "results", "miss_verdicts_all.jsonl"))
    ap.add_argument("--rdchiral", default=os.path.join(ROOT, "results", "rdchiral_verdicts.jsonl"))
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "summary.json"))
    args = ap.parse_args()

    rows = [json.loads(l) for l in open(args.bench)]
    out = {"n_records": len(rows)}

    for d in ("fwd", "retro"):
        ok = [
            r
            for r in rows
            if "hit" in r.get(f"omgkit_{d}", {}) and "hit" in r.get(f"rdkit_{d}", {})
        ]
        a = sum(r[f"omgkit_{d}"]["hit"] for r in ok)
        b = sum(r[f"rdkit_{d}"]["hit"] for r in ok)
        out[d] = {
            "comparable": len(ok),
            "omgkit_hit": a,
            "rdkit_hit": b,
            "only_omgkit": sum(
                1 for r in ok if r[f"omgkit_{d}"]["hit"] and not r[f"rdkit_{d}"]["hit"]
            ),
            "only_rdkit": sum(
                1 for r in ok if r[f"rdkit_{d}"]["hit"] and not r[f"omgkit_{d}"]["hit"]
            ),
            "neither": sum(
                1 for r in ok if not r[f"omgkit_{d}"]["hit"] and not r[f"rdkit_{d}"]["hit"]
            ),
            "shape_skipped": sum(
                1 for r in rows if r.get(f"omgkit_{d}", {}).get("err") == "shape"
            ),
        }

    for k in KEYS:
        ts = sorted(r[k]["t_min"] * 1e6 for r in rows if "t_min" in r.get(k, {}))
        n = len(ts)
        out.setdefault("timing", {})[k] = {
            "n": n,
            "median_us": round(ts[n // 2], 2),
            "p90_us": round(ts[int(0.90 * n)], 2),
            "p99_us": round(ts[int(0.99 * n)], 2),
            "max_us": round(ts[-1], 2),
            "total_s": round(sum(ts) / 1e6, 3),
        }

    if os.path.exists(args.verdicts):
        c = Counter(json.loads(l)["verdict"] for l in open(args.verdicts))
        out["miss_attribution"] = dict(c.most_common())
    if os.path.exists(args.rdchiral):
        c = Counter(json.loads(l)["final"] for l in open(args.rdchiral))
        out["engine_adjudication"] = dict(c.most_common())

    with open(args.out, "w") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=2)

    print(f"记录 {out['n_records']} 条\n")
    print("| 方向 | 可比 | omgkit 命中 | RDKit 命中 | 仅 omgkit | 仅 RDKit | 都未中 | 契约跑不了 |")
    print("|---|---|---|---|---|---|---|---|")
    for d, name in (("fwd", "正向"), ("retro", "逆向")):
        s = out[d]
        print(
            f"| {name} | {s['comparable']} | {s['omgkit_hit']} ({pct(s['omgkit_hit'], s['comparable'])}) "
            f"| {s['rdkit_hit']} ({pct(s['rdkit_hit'], s['comparable'])}) | {s['only_omgkit']} "
            f"| {s['only_rdkit']} | {s['neither']} | {s['shape_skipped']} |"
        )
    print()
    print("| 组 | 中位 | p90 | p99 | 最大 | 合计 |")
    print("|---|---|---|---|---|---|")
    for k in KEYS:
        t = out["timing"][k]
        print(
            f"| {k} | {t['median_us']:.1f} µs | {t['p90_us']:.1f} µs | {t['p99_us']:.1f} µs "
            f"| {t['max_us']:.0f} µs | {t['total_s']:.2f} s |"
        )
    if "miss_attribution" in out:
        print("\n未命中归因(按反应计):")
        for k, v in out["miss_attribution"].items():
            print(f"  {k:34s} {v}")
    if "engine_adjudication" in out:
        print("\n疑似引擎缺陷的 rdchiral 裁定:")
        for k, v in out["engine_adjudication"].items():
            print(f"  {k:34s} {v}")
    print(f"\n已写出 {args.out}")


if __name__ == "__main__":
    main()
