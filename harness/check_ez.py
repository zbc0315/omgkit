#!/usr/bin/env python3
"""用外部实现裁判:写出 SMILES 时双键立体有没有守恒。

`dump_written` 的输出交给这里,每行 `原文<TAB>本实现写出的串`。两边都用外部
实现读进来,比每根双键的 E/Z。

# 为什么要单独一份判据

往返测试(解析 → 写出 → 再解析)只能保证**本实现自己**读得回来。方向键的
参照系换算若整体错了一个符号,自己读自己仍然自洽,顺反却全反了 —— 这类错误
拓扑完全正确,只有分子是镜像的,肉眼极难发现。所以判官必须是外部实现。

# 两种失败不是一回事

- **丢了立体**:少写一条方向。丢信息,不好,但下游拿到的是"未指定",
  不会把它当成确定的顺式或反式。
- **写错/多写**:写出了原文没有的、或与原文相反的立体。这更严重 ——
  下游会当真。

所以两者分开计数,不合并成一个"分歧数"。

用法:

    cargo run --release -p omgkit-io --example dump_written -- \\
        harness/corpus/large.smi > /tmp/written.tsv
    python3 harness/check_ez.py /tmp/written.tsv
"""

import argparse
import collections
import pathlib
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


def stereo_multiset(smi: str):
    """分子里每根双键的立体标注,多重集。读不了返回 None。"""
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    Chem.SetBondStereoFromDirections(m)
    return collections.Counter(
        str(b.GetStereo()) for b in m.GetBonds() if str(b.GetStereo()) != "STEREONONE"
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("tsv", type=pathlib.Path, help="dump_written 的输出")
    ap.add_argument("--limit", type=int, default=8, help="最多打印几条失败")
    args = ap.parse_args()

    stat: collections.Counter = collections.Counter()
    bad = []

    for line in args.tsv.read_text().splitlines():
        if not line.strip():
            continue
        src, got = line.split("\t")
        if got == "<parse-error>":
            stat["本实现解析失败"] += 1
            continue
        want = stereo_multiset(src)
        if want is None:
            # 原文外部实现就读不了 —— 无从比对,不算本实现的账
            stat["原文外部实现读不了"] += 1
            continue
        have = stereo_multiset(got)
        if have is None:
            stat["写出的串外部实现读不了"] += 1
            bad.append((src, got, "外部实现读不了"))
        elif want == have:
            stat["E/Z 守恒"] += 1
        elif len(have) < len(want):
            stat["丢了立体"] += 1
            bad.append((src, got, f"{dict(want)} -> {dict(have)}"))
        else:
            stat["立体写错或多写"] += 1
            bad.append((src, got, f"{dict(want)} -> {dict(have)}"))

    for k, v in stat.most_common():
        print(f"  {k:<22} {v}")
    for s, g, why in bad[: args.limit]:
        print(f"\n  {s}\n    -> {g}\n     ({why})")
    if len(bad) > args.limit:
        print(f"\n  ...(另有 {len(bad) - args.limit} 条)")

    # 写错比丢失严重得多,但两者都不该有
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
