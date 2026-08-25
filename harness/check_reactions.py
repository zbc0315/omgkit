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

#: dump 至少要覆盖多少个分子。**分母闸。**
#:
#: 覆盖面缩水时,`DELIBERATE` 的期望集合会跟着缩(名单里够不着的条目自动
#: 不计入),于是"名单对得上"照样成立 —— 只是在更少的分子上成立。
#: CI 与 `gates.sh` 都跑 300,这里钉死 300。
MIN_MOLS = 300

#: **刻意分歧,按(反应模板, 底物)逐条钉死。**
#:
#: 这些不是缺陷,是本实现刻意选择的语义:产物模板描述的是反应中心的**片段**,
#: 不是"一个片段一个分子"。所有产物模板建进同一张图、模板之外的原子只搬一次,
#: 最后按**连通分量**切开。环状底物 + 断环模板于是给出**一个**开环产物,
#: 而外部实现逐产物各搬一次,共享的那一段被复制进每个产物,原子凭空变多。
#: 详见 `harness/README.md` 的"与外部实现的一处刻意分歧"。
#:
#: # 为什么钉集合而不是设一个上限
#:
#: README 里原本写着"这些的数目变了要重新查",而**没有任何判据在看这个数** ——
#: 这条判官是零容差的,当场就红,所以它一直没进 CI,谁想起来谁跑一次。于是
#: 文档里的 717 / 24 悄悄变成了 716 / 25 没人知道,而多出来的那一条根本不是
#: 这个故事(是双键顺反的参照原子被反应删掉之后整个作废,已修,见
#: `omgkit-match/tests/reaction.rs` 的
#: `bond_stereo_rebases_to_the_other_side_when_nothing_fills_the_slot`)。
#:
#: 上限只能变松;钉死的集合**两个方向都红**:少一条说明改动被撤销或产物语义
#: 变了,多一条说明有新的形状撞上来。两种都要当场查,不是调数。
#:
#: 语料覆盖不到的条目自动不计入期望(比如 `--mols` 只跑前 N 个分子)。
DELIBERATE = frozenset(
    {
        ('[C:1][O:2][C:3]>>[C:1][O:2].[C:3]', '[C@]1([C@H]([C@@]2(C)O[C@@H]1CC2)C(=O)OC)(C(=O)OC)C'),
        ('[C:1][O:2][C:3]>>[C:1][O:2].[C:3]', '[C@]12([C@@H](C3(CO)OC1CC3)COC2=O)C'),
        ('[C:1][O:2][C:3]>>[C:1][O:2].[C:3]', '[C@]12(C(OCC1C(=C)CC[C@@H]2O)=O)C'),
        ('[C:1][O:2][C:3]>>[C:1][O:2].[C:3]', '[C@]12(C3CC([C@@H]1COC2=O)CC3)C'),
        ('[C:1][O:2][C:3]>>[C:1][O:2].[C:3]', '[C@H]12[C@@H](C(CC3[C@]1(C4C(C5CCC([C@](CC4)(C)5)=O)CC3)C)=O)O2'),
        ('[C:1][N:2]>>[C:1].[N:2]', '[C@@H]12[C@@H](NC(=N1)O)N=C(N2)O'),
        ('[C:1][N:2]>>[C:1].[N:2]', '[C@@H]12[C@@H](NS(=O)(=O)N1)[N-]/C(=N\\\\[N+](=O)[O-])/N2'),
        ('[C:1][N:2]>>[C:1].[N:2]', '[H]/[O+]=C\\\\1/C=CC=C/C1=C\\\\2/[NH2+][C@H](CS2)C(=O)[O-]'),
        ('[C:1][N:2]>>[C:1].[N:2]', '[H]/N=C(\\\\C#N)/[C@@](C)([NH+](C)C)SC'),
        ('[C:1][N:2]>>[C:1].[N:2]', '[H]/N=C/1\\\\C(=C(\\\\C#N)/N)\\\\N=CN1N=C(C)C'),
        ('[C:1][N:2]>>[C:1].[N:2]', '[H]/N=C/1\\\\CC2(CC[NH+](CC2)C)c3c(nc([nH]3)COC)O1'),
        ('[C:1][N:2]>>[C:1].[N:2]', '[H]/N=C/1\\\\N[C@]2(CSC(=[NH+]2)N)CS1'),
        ('[C:1][N:2]>>[C:1].[N:2]', '[H]/N=C/1\\\\N=C([C@H](S1)CC(=O)[O-])O'),
        ('[C:1][N:2]>>[C:1].[N:2]', '[H]/N=C/1\\\\NN=C(CS1)CC(=O)OCC'),
        ('[C:1][N:2]>>[C:1].[N:2]', '[H]/N=C\\\\1/C(=O)C=C(NC1=O)ONCCO'),
        ('[C:1][N:2]>>[C:1].[N:2]', '[N+](/C=C/N1CCCCC1)([O-])=O'),
        ('[C:1][N:2]>>[C:1].[N:2]', '[N+](C1C(N2CCCCC2)=CC=NC=1)([O-])=O'),
        ('[C:1][N:2]>>[C:1].[N:2]', '[N+]1(=N/O)/C2=C(C=CC=C2)CCC1'),
        ('[C:1][N:2]>>[C:1].[N:2]', '[N+]1(C(N2CCCCC2)=CC=CC=1)[O-]'),
        ('[C:1][N:2]>>[C:1].[N:2]', '[N+]1(C=C(N2CCCCC2)C=CC=1)[O-]'),
        ('[C:1](=[O:2])[N:3]>>[C:1](=[O:2])[OH].[N:3]', '[H]/N=C\\\\1/C(=O)C=C(NC1=O)ONCCO'),
        ('[C:1][C:2]=[O:3]>>[C:1][C:2][OH:3]', '[H]/[O+]=C\\\\1/C=CC=C/C1=C\\\\2/[NH2+][C@H](CS2)C(=O)[O-]'),
    }
)
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
    if n_mols < MIN_MOLS:
        sys.exit(
            f"产物文件只覆盖 {n_mols} 个分子(至少要 {MIN_MOLS})—— 覆盖面一缩水,\n"
            "刻意分歧那份名单的期望也跟着缩,'名单对得上'就成了在更少的分子上成立"
        )
    smis = smis[:n_mols]

    agreed = only_ours = only_rdkit = different = 0
    with_invalid = 0
    bad = []
    gap = set()
    unpinned = 0

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
        if (rxn_src[ri], smis[mi]) in DELIBERATE:
            # 刻意分歧那一档 —— 收进集合而不是计数,待会儿与名单**逐条**核。
            gap.add((rxn_src[ri], smis[mi]))
            continue
        if len(bad) < args.limit:
            bad.append((rxn_src[ri], smis[mi], theirs, mine))
        unpinned += 1

    # **刻意分歧那一档:与名单逐条核,两个方向都红。** 见 DELIBERATE。
    # 期望只取语料真的覆盖到的那些 —— `--mols` 只跑前 N 个分子时,
    # 名单里够不着的条目不该被当成"少了一条"。
    covered = set(smis)
    want_gap = {(r, m) for (r, m) in DELIBERATE if m in covered}
    # **名单里够不着的条目要报出来。** 不报的话,语料一改(或者 `--mols`
    # 调小),对不上的条目会**静默**退出期望集合 —— 名单照样"对得上",
    # 只是它管的东西少了。变异实测:往名单里加一个语料里没有的分子,
    # 不报的那版退 0。
    orphan = sorted(DELIBERATE - want_gap)

    print(f"一致 {agreed};只有基准 {only_rdkit},只有本实现 {only_ours},产物不同 {different}")
    print(f"其中含无法净化产物的组合 {with_invalid}")
    print(f"刻意分歧(已逐条钉死){len(gap)} 条,名单里 {len(want_gap)} 条;名单外的分歧 {unpinned} 条")
    for r, s, t, o in bad:
        print(f"\n  反应 {r}\n  底物 {s}")
        print(f"    基准   {sorted(t.elements())}")
        print(f"    本实现 {sorted(o.elements())}")
    if bad:
        print(
            f"\n有 {unpinned} 条分歧不在 `DELIBERATE` 名单里 —— 要么是新的回归,\n"
            "要么是一处新的刻意选择。**先查清是哪一种**,别直接往名单里加。"
        )
        return 1
    if orphan:
        print(
            f"\n`DELIBERATE` 名单里有 {len(orphan)} 条语料够不着:{orphan[:3]}\n"
            "语料改了就把名单跟着改 —— 够不着的条目会静默退出期望集合,"
            "名单照样'对得上',只是它管的东西少了。"
        )
        return 1
    if gap != want_gap:
        print(
            f"\n刻意分歧那一档对不上名单 —— 名单里却没红的 {sorted(want_gap - gap)[:5]}。\n"
            "产物语义改了(或者某处改动被撤销了),把名单跟着改,并说明为什么。"
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
