#!/usr/bin/env python3
"""生成子结构匹配的差分基准:「分子语料 × SMARTS 语料」的全部命中。

每行 `分子序号<TAB>模式序号<TAB>命中的原子集合`,与
`cargo run -p omgkit-match --example dump_matches` 的输出格式一致。
只写有命中的组合。

**序号必须与语料的行号一一对应**,解析失败的条目也要占位 —— 序号一旦错位,
比对会把全部结果报成分歧,而且看不出根因。

用法:

    python3 harness/oracle_matches.py \\
        --mols harness/corpus/large.smi \\
        --pats harness/corpus/smarts.txt \\
        --limit-mols 2000 \\
        --out harness/baseline/matches.tsv
"""

import argparse
import pathlib
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

# 与 Rust 侧保持一致。RDKit 默认 1000,超过会静默截断 —— 两边必须用同一个值,
# 否则高度对称的分子上会出现"只是截断点不同"的假分歧。
MAX_MATCHES = 1000

# **必须关掉 removeHs**。`MolFromSmiles` 默认会删掉显式 `[H]` 原子并把它折成
# 氢计数,而 omgkit 把 removeHs 划在净化之外、不做这一步。两边不一致时,
# 带显式氢的分子上原子编号会整体错位,氢计数也对不上 —— 表现为"匹配莫名其妙
# 地不命中",根因却在语料对齐上,极难往这个方向想。
PARAMS = Chem.SmilesParserParams()
PARAMS.removeHs = False


def read_lines(path: pathlib.Path, limit: int) -> list[str]:
    out = []
    for line in path.read_text(errors="replace").splitlines():
        tok = line.split()[0] if line.split() else ""
        if not tok or tok.startswith("#"):
            continue
        out.append(tok)
        if len(out) >= limit:
            break
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mols", required=True, type=pathlib.Path)
    ap.add_argument("--pats", required=True, type=pathlib.Path)
    ap.add_argument("--out", required=True, type=pathlib.Path)
    ap.add_argument("--limit-mols", type=int, default=1 << 30)
    ap.add_argument("--limit-pats", type=int, default=1 << 30)
    args = ap.parse_args()

    smis = read_lines(args.mols, args.limit_mols)
    pats_raw = read_lines(args.pats, args.limit_pats)

    mols = [Chem.MolFromSmiles(s, PARAMS) for s in smis]
    pats = [Chem.MolFromSmarts(s) for s in pats_raw]
    n_mol_ok = sum(m is not None for m in mols)
    n_pat_ok = sum(p is not None for p in pats)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    pairs = 0
    with_hits = 0
    with args.out.open("w") as f:
        # 覆盖范围必须**随数据走**:这份基准可能只跑了前 N 个分子,而比对方读的是
        # 整份语料。两个数对不上时,第 N 个之后的分子在比对方眼里全是"只有本实现
        # 有命中",凭空多出成千上万条假分歧。所以范围写进首行,读不到就报错 ——
        # 靠文档提醒对齐守不住,靠数据自己携带才守得住。
        f.write(f"#mols\t{len(smis)}\n")
        for mi, m in enumerate(mols):
            if m is None:
                continue
            for pi, p in enumerate(pats):
                if p is None:
                    continue
                pairs += 1
                # useChirality=True:查询里写的手性与顺反算数。
                # 外部实现**默认是关的**,而本库默认判 —— 口径不一致的话,
                # 2000 分子 × 776 模式里有 23% 的组合会对不上,报出来全是噪声。
                hits = m.GetSubstructMatches(
                    p, uniquify=True, maxMatches=MAX_MATCHES, useChirality=True
                )
                if not hits:
                    continue
                with_hits += 1
                sets = sorted(",".join(str(i) for i in sorted(h)) for h in hits)
                f.write(f"{mi}\t{pi}\t{'|'.join(sets)}\n")

    print(f"已写出 {args.out}")
    print(f"  分子 {len(smis)}(可解析 {n_mol_ok}),模式 {len(pats_raw)}(可解析 {n_pat_ok})")
    print(f"  组合 {pairs},有命中 {with_hits}")
    print(f"  RDKit: {Chem.rdBase.rdkitVersion}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
