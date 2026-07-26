#!/usr/bin/env python3
"""用外部实现裁判产物生成。

`dump_reactions` 的输出交给这里,两边的产物**都由外部实现规范化**再比。

# 为什么不能各自规范化后比字符串

两边的规范化是两套不同的算法,同一个分子写出来的字符串本来就不一样。拿本
实现的规范 SMILES 去比外部实现的规范 SMILES,比出来的"分歧"全是噪声。

判官必须是同一个:两边的产物都交给外部实现读进来再规范化。

# 判据是多重集

同一条反应在一个底物上可能有多处反应位点,产物之间可能重复。比集合会把
"两条路径"和"一条路径"混为一谈。

用法:

    cargo run --release -p omgkit-match --example dump_reactions -- \\
        harness/corpus/reactions.txt harness/corpus/large.smi 300 > /tmp/rx.tsv
    python3 harness/check_reactions.py /tmp/rx.tsv \\
        --rxns harness/corpus/reactions.txt \\
        --mols harness/corpus/large.smi
"""

import argparse
import collections
import pathlib
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

MAX_PRODUCT_SETS = 100
INVALID = "<invalid>"

PARAMS = Chem.SmilesParserParams()
PARAMS.removeHs = False


def read_corpus(path: pathlib.Path) -> list[str]:
    out = []
    for line in path.read_text(errors="replace").splitlines():
        tok = line.split()[0] if line.split() else ""
        if tok and not tok.startswith("#"):
            out.append(tok)
    return out


def canonical(smi: str) -> str:
    """把一个产物 SMILES 规范化。读不了的记成 invalid —— 与两侧的约定一致。"""
    if smi == INVALID:
        return INVALID
    m = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(m) if m is not None else INVALID


def canonical_group(group: str) -> str:
    return ".".join(sorted(canonical(s) for s in group.split(".")))


def rdkit_products(rxn, mol) -> collections.Counter:
    try:
        sets = rxn.RunReactants((mol,))
    except Exception:
        return collections.Counter()
    out = collections.Counter()
    for tup in sets[:MAX_PRODUCT_SETS]:
        parts = []
        for p in tup:
            try:
                q = Chem.Mol(p)
                Chem.SanitizeMol(q)
                parts.append(Chem.MolToSmiles(q))
            except Exception:
                parts.append(INVALID)
        out[".".join(sorted(parts))] += 1
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tsv", type=pathlib.Path, help="dump_reactions 的输出")
    ap.add_argument("--rxns", required=True, type=pathlib.Path)
    ap.add_argument("--mols", required=True, type=pathlib.Path)
    ap.add_argument("--limit", type=int, default=8, help="最多打印几条分歧")
    args = ap.parse_args()

    rxn_src = read_corpus(args.rxns)
    smis = read_corpus(args.mols)
    rxns = []
    for r in rxn_src:
        try:
            rxns.append(AllChem.ReactionFromSmarts(r))
        except Exception:
            rxns.append(None)

    # 覆盖范围由数据自己携带 —— dump 可以只跑前 N 个分子,而这里读的是整份
    # 语料。两个数对不上时,第 N 个之后的分子在这里全变成"只有基准有产物",
    # 凭空多出成千上万条假分歧。所以宁可硬报错,也不默默按整份语料比。
    n_mols = None
    ours: dict[tuple[int, int], collections.Counter] = {}
    for line in args.tsv.read_text().splitlines():
        if not line.strip():
            continue
        if line.startswith("#mols\t"):
            n_mols = int(line.split("\t")[1])
            continue
        ri, mi, groups = line.split("\t")
        ours[(int(ri), int(mi))] = collections.Counter(
            canonical_group(g) for g in groups.split("|")
        )

    if n_mols is None:
        sys.exit(
            "产物文件缺少 `#mols<TAB>N` 首行,无从知道它覆盖了多少个分子。\n"
            "用当前版本的 dump_reactions 重新生成:\n"
            "  cargo run --release -p omgkit-match --example dump_reactions -- \\\n"
            "      <反应.txt> <分子.smi> [分子数上限] > out.tsv"
        )
    if n_mols > len(smis):
        sys.exit(f"产物文件覆盖 {n_mols} 个分子,而 --mols 只有 {len(smis)} 个 —— 语料对不上")
    smis = smis[:n_mols]

    agreed = only_ours = only_rdkit = different = 0
    with_invalid = 0
    bad = []

    # 以两侧的并集为准 —— 只遍历本实现的话,"本实现少了一条"就查不出来
    keys = set(ours)
    for mi, s in enumerate(smis):
        m = Chem.MolFromSmiles(s, PARAMS)
        if m is None:
            continue
        for ri, rxn in enumerate(rxns):
            if rxn is None or rxn.GetNumReactantTemplates() != 1:
                continue
            if rdkit_products(rxn, m):
                keys.add((ri, mi))

    for ri, mi in sorted(keys):
        m = Chem.MolFromSmiles(smis[mi], PARAMS)
        theirs = rdkit_products(rxns[ri], m) if m is not None else collections.Counter()
        mine = ours.get((ri, mi), collections.Counter())
        if any(INVALID in k for k in theirs) or any(INVALID in k for k in mine):
            with_invalid += 1
        if theirs == mine:
            agreed += 1
            continue
        if not mine:
            only_rdkit += 1
        elif not theirs:
            only_ours += 1
        else:
            different += 1
        if len(bad) < args.limit:
            bad.append((rxn_src[ri], smis[mi], theirs, mine))

    print(f"一致 {agreed};只有基准 {only_rdkit},只有本实现 {only_ours},产物不同 {different}")
    print(f"其中含无法净化产物的组合 {with_invalid}")
    for r, s, t, o in bad:
        print(f"\n  反应 {r}\n  底物 {s}")
        print(f"    基准   {sorted(t.elements())}")
        print(f"    本实现 {sorted(o.elements())}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
