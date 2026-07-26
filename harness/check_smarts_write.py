#!/usr/bin/env python3
"""用外部实现裁判:写出的 SMARTS 与原 SMARTS **语义相同**吗。

`dump_smarts_written` 的输出交给这里,每行 `原文<TAB>本实现写出的串`。

# 判据是"匹配到的东西一样",不是字符串一样

同一个查询有很多写法(`[CH1]`、`[C&H1]`、`[H1&C]` 都一样),逐字节比对比的是
写法而不是语义。所以两个 SMARTS 都交给外部实现,去匹配同一批分子,比匹配
到的东西。

## 比的是**集合**,不是匹配元组

匹配元组的第 k 位对应"查询里第 k 个原子",而写出会**重排查询原子的编号** ——
同一处匹配在两边会给出顺序不同的元组。拿元组直接比,报出来的全是编号差异。

判据只要碰到"两边的编号可能不同",就得先想清楚比的是不是编号无关的量 ——
SMILES 往返、SMARTS 表达式树、这里,三处都要各自处理。

往返幂等(在 Rust 侧的 roundtrip_smarts 测试里)只保证"解析→写出→解析这一趟
没丢信息",保证不了语义没变 —— 一个系统性写错的运算符可以既幂等又是错的。
这一档补的正是那个缺口。

# 空匹配不算通过

一条 SMARTS 若两边都匹配不到任何东西,比出来当然"一致",但什么也没验证。
所以单独统计**确实匹配到东西**的条数,并在它太少时报错 —— 语料换了、探针
分子集换了都可能让判据悄悄变空。

用法:

    cargo run --release -p omgkit-io --example dump_smarts_written -- \\
        harness/corpus/smarts.txt > /tmp/sw.tsv
    python3 harness/check_smarts_write.py /tmp/sw.tsv \\
        --mols harness/corpus/large.smi
"""

import argparse
import collections
import pathlib
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


def load_probes(path: pathlib.Path, limit: int):
    """探针分子。解析不了的跳过 —— 那是分子语料的事,不在本判据范围内。"""
    out = []
    for line in path.read_text(errors="replace").splitlines():
        tok = line.split()[0] if line.split() else ""
        if not tok or tok.startswith("#"):
            continue
        m = Chem.MolFromSmiles(tok)
        if m is not None:
            out.append(m)
        if len(out) >= limit:
            break
    return out


def matches(patt, probes):
    """该模式在探针集上匹配到的东西:{分子序号: 匹配元组的集合}。"""
    out = {}
    for i, m in enumerate(probes):
        got = m.GetSubstructMatches(patt, uniquify=True, maxMatches=200)
        if got:
            # frozenset:去掉查询原子编号的影响,只留"匹配到了分子里的哪些原子"
            out[i] = {frozenset(t) for t in got}
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("tsv", type=pathlib.Path, help="dump_smarts_written 的输出")
    ap.add_argument("--mols", type=pathlib.Path, required=True, help="探针分子语料")
    ap.add_argument("--limit-mols", type=int, default=3000, help="用前几个分子当探针")
    ap.add_argument(
        "--min-hits",
        type=int,
        default=150,
        help="至少要有这么多条真的匹配到东西,低于此值说明判据快空过了",
    )
    ap.add_argument("--limit", type=int, default=8, help="最多打印几条失败")
    args = ap.parse_args()

    probes = load_probes(args.mols, args.limit_mols)
    if not probes:
        sys.exit("探针分子一个都没读到")

    stat: collections.Counter = collections.Counter()
    bad = []
    hits = 0

    for line in args.tsv.read_text().splitlines():
        if not line.strip():
            continue
        src, _, got = line.partition("\t")
        if got.startswith("<"):
            stat["本实现解析失败(语料里的非法输入)"] += 1
            continue
        pa = Chem.MolFromSmarts(src)
        if pa is None:
            # 外部实现读不了原文 —— 无从比对,不算本实现的账
            stat["原文外部实现读不了"] += 1
            continue
        pb = Chem.MolFromSmarts(got)
        if pb is None:
            stat["写出的串外部实现读不了"] += 1
            bad.append((src, got, "外部实现读不了"))
            continue
        ma, mb = matches(pa, probes), matches(pb, probes)
        if ma:
            hits += 1
        if ma == mb:
            stat["语义相同"] += 1
            continue
        stat["语义不同"] += 1
        only_a = {k: v for k, v in ma.items() if mb.get(k) != v}
        only_b = {k: v for k, v in mb.items() if ma.get(k) != v}
        bad.append((src, got, f"原文独有 {list(only_a)[:3]} / 写出独有 {list(only_b)[:3]}"))

    for k, v in stat.most_common():
        print(f"  {k:<28} {v}")
    print(f"  {'其中确实匹配到东西':<28} {hits}")
    for s, g, why in bad[: args.limit]:
        print(f"\n  {s}\n    -> {g}\n     ({why})")
    if len(bad) > args.limit:
        print(f"\n  ...(另有 {len(bad) - args.limit} 条)")

    if hits < args.min_hits:
        print(f"\n判据几乎是空过的:只有 {hits} 条真的匹配到东西(要求 ≥ {args.min_hits})")
        return 1
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
