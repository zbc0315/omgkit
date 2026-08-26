#!/usr/bin/env python3
"""用外部实现裁判 omgkit 的 SMILES 写出。

往返测试(`cargo test --test roundtrip_smiles`)只用自家解析器验证写出,
原理上留了一个漏网可能:**解析与写出共享同一个误解**,互为逆运算却都偏离
SMILES 的语义。手性尤其危险 —— 标记写反了,原子数、键集合、连通性全都对,
只有分子是镜像的,纯拓扑比对永远发现不了。

本脚本把两边各自规范化再比字符串,判官是外部实现,与 omgkit 无共谋。

用法:

    cargo run --release -p omgkit-io --example write_smiles -- \\
        harness/corpus/large.smi > /tmp/written.tsv
    python3 harness/check_write.py /tmp/written.tsv harness/corpus/large.smi

输入是两列 TSV:原始 SMILES、omgkit 写出的 SMILES。第二个参数是**同一份语料**,
用来核分母。

尚未写出的立体信息会被分桶而不是算作失败 —— 见下面 `--strict` 的说明。

# 分母也要有闸

这条判据的分母是上游喂的,而它只数分歧、不数"该数到多少"。实测:

- **喂一个空文件进去,它打印"零分歧"并退 0**;
- 把 19 行 TSV 换成垃圾,`逐条规范形式相同` 从 8833 悄悄掉到 8814,照样退 0。

所以现在要第二个参数(语料),而且数的是**真正比对过的分子数**,不是 TSV 行数
—— 行数掉不下来。少比到的超过 `MAX_UNCOMPARED` 就判失败,反向(TSV 比语料还长)
也判失败,否则传错语料时分母闸会静默失效。

同一条在 `harness/check_wedge_readback.py` 上栽过,那边是"dump 少喂几个分子,
每一档都会变好看"。
"""

import argparse
import pathlib
import sys

import rdkit
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

# 关掉 removeHs:它会删掉显式 [H],而**带立体信息的 [H] 会被保留**。
# 于是"没写出双键立体"就间接表现为原子数不同,把一个已知缺口伪装成
# 结构错误。关掉之后两边的图才对得齐。
PARAMS = Chem.SmilesParserParams()
PARAMS.removeHs = False


# 尚未写出的立体类别 —— **现在一个都没有,所以这里是空的。**
#
# 差别仅限于名单里某一类时不算失败,而是分桶计数(已登记的缺口,不是结构错误)。
# 先前这里登记着"双键立体"与"配位几何立体"两项;两者都写得出来之后名单清空,
# 于是加不加 `--strict` 是一回事。
#
# **名单留着,但它只能往绿的方向拨,所以必须配一道闸。** 那道闸就是
# `--strict`:`harness/gates.sh` 与 CI 的三处调用全都带着它,往名单里加一项
# 换不来绿。真要登记新的一类,先把闸摘掉才行 —— 那是个显式动作,看得见。
NOT_YET_WRITTEN = {}


#: 语料里允许有几行**没真正被比对到**。
#:
#: 数的不是 TSV 行数,是 `exact + 各档豁免 + 分歧` —— 也就是真的被判过的分子。
#: 头一版数的是行数,而行数掉不下来:`a is None`(外部实现读不了原串)那条
#: `continue` 既不计数也没上限。独立审核实测:把 19 行 TSV 换成垃圾,
#: `逐条规范形式相同` 从 8833 悄悄掉到 8814,判据照样打印"零分歧"退 0。
#:
#: 少比到的两条路:
#:
#: - **没写出**:本实现解析不了,`write_smiles` 直接跳过
#: - **外部实现读不了原串**:判官自己解析不了输入,`a is None`
#:
#: 后一条是**版本相关**的,所以上限按更严的那个版本定。逐版实测:
#:
#: | 语料 | RDKit 2022.09.5 | RDKit 2025.09.2(CI 钉的) |
#: |---|---|---|
#: | `large.smi`(8839 行) | 没比到 **6** | 没比到 **8** |
#: | `smoke.smi`(149 行) | 没比到 **12**(8 条故意解析不了 + 4 条判官读不了) | 没比到 **12** |
#:
#: 现值 15 = 实测最大 12 加一点余量。变异实测:把 19 行 TSV 换成垃圾,
#: 没比到变 25,当场退 1。
#:
#: 这是**分母闸**,不是宽容度:它管的是"少比几个分子,分歧数跟着变好看"。
#: 涨上去说明有一类分子进不了比对,要当场查,不是调大它。
MAX_UNCOMPARED = 15

#: **写出器写不出来的那些分子**,按 SMILES 逐条钉死,不是一个上限。
#:
#: 现在是**空的** —— `@SP` / `@TB` / `@OH` 三类都写得出来了。机制留着,因为它是
#: 双向的:写出器哪天在某个分子上退化,那个分子不在名单里,照旧算分歧;
#: 而要往名单里加一条,得先写清楚为什么。
#:
#: 它一度有 6 条。这件事先前**一条判据都没守**:CI 只拿 `large.smi` 跑这条判官,
#: 而那份语料里一条非四面体立体都没有。冒烟语料里有 6 条,`--strict` 下当场
#: 6 条分歧,只是从来没人在 CI 里跑过(`check_write_fidelity.py` 也一样)。
#: 补写出的时候,正是这条判据的双向钉子把它们逐条点名报出来的。
#:
NON_TETRAHEDRAL_GAP = frozenset()


def canonical_without(mol, drop):
    """抹掉若干类立体信息后的规范 SMILES。`mol` 会被就地修改,故传副本。"""
    for name in drop:
        NOT_YET_WRITTEN[name](mol)
    return Chem.MolToSmiles(mol)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tsv", type=pathlib.Path, help="write_smiles 的输出")
    ap.add_argument("corpus", type=pathlib.Path, help="喂给 write_smiles 的那份语料(核分母用)")
    ap.add_argument("--limit", type=int, default=10, help="最多打印几条分歧")
    ap.add_argument(
        "--strict",
        action="store_true",
        help="不给尚未实现的立体信息开豁免。相应写出做完之后应当打开",
    )
    args = ap.parse_args()

    exact = 0
    chiral = 0
    excused = {name: 0 for name in NOT_YET_WRITTEN}
    bad = []
    gap = set()
    rows = 0
    unreadable = 0

    for line in args.tsv.read_text().splitlines():
        parts = line.split("\t")
        if len(parts) != 2:
            continue
        rows += 1
        orig, got = parts
        a = Chem.MolFromSmiles(orig, PARAMS)
        b = Chem.MolFromSmiles(got, PARAMS)
        if a is None:
            # 原始就读不了,不是写出的问题 —— 但**要计数**:它是"少比一个分子"的
            # 两条路之一,而少比的分子会让下面每一个数变好看。见 MAX_UNCOMPARED。
            unreadable += 1
            continue
        if b is None:
            bad.append((orig, got, "写出结果无法解析"))
            continue

        if any(at.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for at in a.GetAtoms()):
            chiral += 1

        if Chem.MolToSmiles(a) == Chem.MolToSmiles(b):
            exact += 1
            continue

        # 逐类豁免:抹掉某一类尚未写出的立体信息后若一致,说明差别**仅限于**
        # 那一类,不是结构错误。归到对应的桶里而不是算失败。
        matched = None
        if not args.strict:
            for name in NOT_YET_WRITTEN:
                if canonical_without(Chem.Mol(a), [name]) == canonical_without(
                    Chem.Mol(b), [name]
                ):
                    matched = name
                    break
        if matched:
            excused[matched] += 1
        elif orig in NON_TETRAHEDRAL_GAP:
            # 已知未实现,**逐条钉死**的那一档 —— 见 NON_TETRAHEDRAL_GAP。
            # 收进集合而不是计数:待会儿要与名单**逐条**核,少一条也红。
            gap.add(orig)
        else:
            bad.append((orig, got, f"{Chem.MolToSmiles(a)}\n       != {Chem.MolToSmiles(b)}"))

    # **分母闸。** 见 MAX_UNCOMPARED:少比几个分子,上面每一个计数都会变好看。
    n_corpus = sum(
        1
        for l in args.corpus.read_text(encoding="utf-8").splitlines()
        if l.strip() and not l.lstrip().startswith("#")
    )
    compared = exact + sum(excused.values()) + len(bad) + len(gap)
    uncompared = n_corpus - compared

    # **已知未实现那一档:与名单逐条核,两个方向都红。** 见 NON_TETRAHEDRAL_GAP。
    # 语料里根本没有这些分子时(`large.smi`)期望集合是空的,这一段自然不作声。
    corpus_smis = {
        l.split()[0]
        for l in args.corpus.read_text(encoding="utf-8").splitlines()
        if l.strip() and not l.lstrip().startswith("#") and l.split()
    }
    want_gap = NON_TETRAHEDRAL_GAP & corpus_smis

    print(f"外部实现:RDKit {rdkit.__version__}")
    print(f"逐条规范形式相同 {exact} 条;含立体中心 {chiral} 条")
    print(
        f"语料 {n_corpus} 行,写出 {rows} 行,真正比对 {compared} 条,"
        f"没比到 {uncompared} 条(上限 {MAX_UNCOMPARED})"
    )
    print(f"  没比到的两条路:没写出 {n_corpus - rows} 条,外部实现读不了原串 {unreadable} 条")
    if want_gap or gap:
        print(f"  已知未实现(非四面体立体写出){len(gap)} 条,名单里 {len(want_gap)} 条")
    for name, count in excused.items():
        if count:
            print(f"仅 {name} 不同 {count} 条({name}写出尚未实现,已登记)")
    if gap != want_gap:
        print(
            f"\n已知未实现那一档对不上名单 —— 多出来的 {sorted(gap - want_gap)},"
            f"名单里却没红的 {sorted(want_gap - gap)}。\n"
            "写出器补上了就把它从 `NON_TETRAHEDRAL_GAP` 里划掉;多出来的是新的回归。"
        )
        return 1
    if bad:
        print(f"\n分歧 {len(bad)} 条(最多列 {args.limit} 条):")
        for orig, got, why in bad[: args.limit]:
            print(f"  原: {orig}\n  写: {got}\n  {why}\n")
        return 1
    # 传错语料时 `uncompared` 会变成负数,而负数当然 <= 上限 —— 那样分母闸就
    # 静默失效了。实测:拿 large 的 TSV 配 smoke 的语料,它打印"没写出 -8690 行"
    # 然后退 0。所以两个方向都要闸。
    if rows > n_corpus:
        print(
            f"\n写出的行数({rows})比语料还多({n_corpus})—— 十有八九是 TSV 与语料"
            "对不上(传错文件了)。分母核不了,这条判据算出来的数没有意义"
        )
        return 1
    if uncompared > MAX_UNCOMPARED:
        print(
            f"\n语料里有 {uncompared} 条没真正被比对,超过上限 {MAX_UNCOMPARED} —— "
            "分歧数是分子,覆盖面是分母。别调大这个数,先查是哪一类分子进不了比对"
        )
        return 1
    print("零分歧" + ("(严格模式)" if args.strict else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
