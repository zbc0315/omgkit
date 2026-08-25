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
        --mols harness/corpus/large.smi --corpus harness/corpus/smarts.txt
"""

import argparse
import collections
import pathlib
import sys

import denominator
import rdkit
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


#: SMARTS 语料里允许有多少条**没真正进比对**。见 `denominator.py`。
#:
#: 这条判官原本只有 `--min-hits`(区分力闸:"真的匹配到东西的模式数"),
#: 那守的是"判据别空过",守不住"上游少喂了几条模式" —— 少喂几条,分歧数
#: 跟着变好看,而 `hits` 只要还够就照样退 0。两条闸各守一半。
#:
#: 实测 `smarts.txt`(776 条)没比到 **20** —— 全是语料里**故意写坏**的模式,
#: 本实现解析失败。现值 25 = 实测值加一点余量:解析器要是在合法模式上退化,
#: 这个数会往上走,当场红。
MAX_UNCOMPARED = 25


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("tsv", type=pathlib.Path, help="dump_smarts_written 的输出")
    ap.add_argument(
        "--corpus", type=pathlib.Path, required=True, help="喂给 dump 的那份 SMARTS 语料(核分母用)"
    )
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
    rows = 0
    n_mapped = 0

    for line in args.tsv.read_text().splitlines():
        if not line.strip():
            continue
        rows += 1
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
        # **原子映射号不参与匹配** —— 上面那条语义判据对它是**瞎的**:
        # 写出器把 `:n` 整批丢掉,匹配到的东西一个不差,这里照样"语义相同"。
        # 所以单独数一遍。
        #
        # 注意这一档在**当前语料上是空过的**:`smarts.txt` 的 776 条模式里
        # 一条映射号都没有(下面 `n_mapped` 会把这个数打出来,别让它静默)。
        # 真正带映射号的语料是 `reactions.txt`,由 Rust 侧的
        # `roundtrip_smarts.rs::reaction_templates_keep_their_atom_maps` 守着。
        maps_a = sorted(a.GetAtomMapNum() for a in pa.GetAtoms() if a.GetAtomMapNum())
        maps_b = sorted(a.GetAtomMapNum() for a in pb.GetAtoms() if a.GetAtomMapNum())
        if maps_a:
            n_mapped += 1
        if maps_a != maps_b:
            stat["原子映射号变了"] += 1
            bad.append((src, got, f"映射号 {maps_a} -> {maps_b}"))
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

    # **分母闸。** 见 `MAX_UNCOMPARED`:`--min-hits` 守的是"判据别空过",
    # 这一条守的是"上游别少喂"—— 两回事。
    n_corpus = denominator.corpus_size(args.corpus)
    compared = stat["语义相同"] + len(bad)

    print(f"外部实现:RDKit {rdkit.__version__}")
    for k, v in stat.most_common():
        print(f"  {k:<28} {v}")
    print(f"  {'其中确实匹配到东西':<28} {hits}")
    # **把 0 打出来。** 这份语料一条映射号都没有,上面那条映射号判据因此是
    # 空过的 —— 空过要看得见,不能只在源码注释里写着。
    print(f"  {'其中带原子映射号':<28} {n_mapped}(为 0 说明这份语料验不了映射号)")
    print(denominator.line(n_corpus, rows, compared, MAX_UNCOMPARED))
    for s, g, why in bad[: args.limit]:
        print(f"\n  {s}\n    -> {g}\n     ({why})")
    if len(bad) > args.limit:
        print(f"\n  ...(另有 {len(bad) - args.limit} 条)")

    if hits < args.min_hits:
        print(f"\n判据几乎是空过的:只有 {hits} 条真的匹配到东西(要求 ≥ {args.min_hits})")
        return 1
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
