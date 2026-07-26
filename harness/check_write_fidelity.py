#!/usr/bin/env python3
"""用外部实现裁判:**净化之后**写出的 SMILES 忠不忠实。

`dump_sanitized` 的输出交给这里,每行 `原文<TAB>本实现写出的串`。

# 判据:两边交给同一个读者,得到同一个分子

比的是 `规范(外部实现读原文)` 与 `规范(外部实现读本实现写出的串)`。两边都由
外部实现读、外部实现规范化,所以比的**只是写出这一步忠不忠实** —— 本实现的
净化与外部实现有出入也不会混进来,那由 L2 的差分单独守。

不能只比分子式:方向键写反、手性写反都不改分子式,却是实打实的另一个分子。
规范串相同这条判据把它们一并盖住,而且严格更强。

# 为什么净化那一档要单独测

不净化的往返测试(解析 → 写出 → 再解析)守不住这个:净化会**重排氢的存放
位置** —— 第 12 步把一部分隐式氢挪进 `num_explicit_hs`,同时清掉 NO_IMPLICIT
标志。写出器若按标志决定要不要写方括号,吡咯型氮的 `[nH]` 就会被写成裸 `n`,
氢凭空消失。

这个缺口在两处都是看不见的:

- 不净化的往返测试里,标志还在,方括号照写,一切正常
- L2 的差分测试比的是**分子对象**的字段,不比写出的字符串

这类串坏得很彻底:氢丢了之后连凯库勒化都做不成。反应产物尤其受害 ——
产物必然是净化过的。

用法:

    cargo run --release -p omgkit-chem --example dump_sanitized -- \\
        harness/corpus/large.smi > /tmp/san.tsv
    python3 harness/check_write_fidelity.py /tmp/san.tsv
"""

import argparse
import collections
import pathlib
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors

RDLogger.DisableLog("rdApp.*")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("tsv", type=pathlib.Path, help="dump_sanitized 的输出")
    ap.add_argument("--limit", type=int, default=8, help="最多打印几条失败")
    args = ap.parse_args()

    stat: collections.Counter = collections.Counter()
    bad = []

    for line in args.tsv.read_text().splitlines():
        if not line.strip():
            continue
        src, got = line.split("\t")
        if got.startswith("<"):
            # <parse-error> / <sanitize-error> —— 本实现自己报的,不在本判据范围内
            stat[got] += 1
            continue
        want = Chem.MolFromSmiles(src)
        if want is None:
            stat["原文外部实现读不了"] += 1
            continue
        have = Chem.MolFromSmiles(got)
        if have is None:
            # 写出的串外部实现读不了 —— 通常是丢了氢导致凯库勒化失败
            stat["写出的串外部实现读不了"] += 1
            bad.append((src, got, "读不了"))
            continue
        cw = Chem.MolToSmiles(want)
        ch = Chem.MolToSmiles(have)
        if cw == ch:
            stat["写出忠实"] += 1
            continue
        stat["写出不忠实"] += 1
        # 分子式只当诊断线索:变了就是丢了原子,没变则多半是立体或电荷
        fw = rdMolDescriptors.CalcMolFormula(want)
        fh = rdMolDescriptors.CalcMolFormula(have)
        hint = f"分子式 {fw} -> {fh}" if fw != fh else "分子式没变,多半是立体或电荷"
        bad.append((src, got, f"{cw}\n        vs {ch}\n     ({hint})"))

    for k, v in stat.most_common():
        print(f"  {k:<24} {v}")
    for s, g, why in bad[: args.limit]:
        print(f"\n  {s}\n    -> {g}\n     ({why})")
    if len(bad) > args.limit:
        print(f"\n  ...(另有 {len(bad) - args.limit} 条)")

    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
