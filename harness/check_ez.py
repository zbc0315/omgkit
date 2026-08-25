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
    python3 harness/check_ez.py /tmp/written.tsv harness/corpus/large.smi
"""

import argparse
import collections
import pathlib
import sys

import denominator
import rdkit
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


#: 语料里允许有多少条**没真正进比对**。见 `denominator.py`。
#:
#: 实测:`large.smi`(8839 行)没比到 **6**(2022.09.5)/ **8**(2025.09.2,CI 装的),
#: 全是外部实现读不了原串;
#: `smoke.smi`(149 行)没比到 **12**(8 条本实现故意解析不了 + 4 条判官读不了)。
#: 现值 15 = 实测最大加一点余量。这是分母闸,不是宽容度。
MAX_UNCOMPARED = 15


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("tsv", type=pathlib.Path, help="dump_written 的输出")
    ap.add_argument("corpus", type=pathlib.Path, help="喂给上游 dump 的那份语料(核分母用)")
    ap.add_argument("--limit", type=int, default=8, help="最多打印几条失败")
    args = ap.parse_args()

    stat: collections.Counter = collections.Counter()
    bad = []
    rows = 0

    for line in args.tsv.read_text().splitlines():
        if not line.strip():
            continue
        rows += 1
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

    # **分母闸。** 见 `MAX_UNCOMPARED` 与 `denominator.py`:少比几个分子,
    # 上面每一档都会变好看,而"分歧 0"是它退 0 的唯一依据。
    n_corpus = denominator.corpus_size(args.corpus)
    compared = stat["E/Z 守恒"] + len(bad)

    print(f"外部实现:RDKit {rdkit.__version__}")
    for k, v in stat.most_common():
        print(f"  {k:<22} {v}")
    print(denominator.line(n_corpus, rows, compared, MAX_UNCOMPARED))
    for s, g, why in bad[: args.limit]:
        print(f"\n  {s}\n    -> {g}\n     ({why})")
    if len(bad) > args.limit:
        print(f"\n  ...(另有 {len(bad) - args.limit} 条)")

    if bad:
        return 1
    why = denominator.verdict(n_corpus, rows, compared, MAX_UNCOMPARED)
    if why:
        print(f"\n{why}")
        return 1
    print("零分歧")
    return 0


if __name__ == "__main__":
    sys.exit(main())
