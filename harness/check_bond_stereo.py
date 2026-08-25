#!/usr/bin/env python3
"""用外部实现裁判:双键顺反的感知结果。

`dump_bond_stereo` 的输出交给这里,每行
`SMILES<TAB>begin,end,顺反,参照a,参照b;...`。

# 判据是"哪根双键、什么顺反、相对谁"

只比顺反值是不够的 —— 顺反离开参照原子没有意义。四取代双键上参照挑得不同,
同一个几何会得出相反的顺反值。所以三样一起比,并且**把参照归一**:
参照对 `(a, b)` 换成双键另一侧的邻居时顺反要跟着翻,归一之后才能比。

# 不比 E/Z

E/Z 要 CIP 优先级,那是另一件事。这里比的是"相对记录下来的参照原子"的
顺反,不涉及 CIP。外部实现里对应的是 `SetBondStereoFromDirections` 给出的
STEREOCIS/STEREOTRANS,而不是 `AssignStereochemistry` 给出的 STEREOE/STEREOZ。

用法:

    cargo run --release -p omgkit-io --example dump_bond_stereo -- \\
        harness/corpus/large.smi > /tmp/bs.tsv
    python3 harness/check_bond_stereo.py /tmp/bs.tsv harness/corpus/large.smi
"""

import argparse
import collections
import pathlib
import sys

import denominator
import rdkit
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


# 高层解析:全套净化 + 立体感知,而且**留着显式氢**,所以下标与
# `MolFromSmiles(sanitize=False)` 那一份逐个对齐(默认的 `removeHs=True` 会
# 把 `[H]/N=C…` 开头的氢删掉,整串下标错位)。
_KEEP_HS = Chem.SmilesParserParams()
_KEEP_HS.removeHs = False


def stereogenic(smi: str, ref):
    """RDKit **自己**认哪些双键是立体元素 —— 用来给方向标记那一层过筛。

    `SetBondStereoFromDirections` 是低层原语:它把 `/` `\\` 机械地折算成
    CIS/TRANS,**不问这根双键有没有资格带顺反**。RDKit 的高层解析问 ——
    小环里的双键(最小环 < 8)被环锁死,没有 E/Z 可言,反式环辛烯是最小的
    能分离出来的反式环烯烃。

    实测大语料 8833 个分子,两者只在**一根**键上分歧:
    `CN1CCC\\2=C1/C(=N\\O)/S/C2=N\\c3ccc(cc3)F` 的 4=5,它在一个**五元环**里。
    方向标记那一层照折算,高层解析给 STEREONONE。本实现跟高层走。

    这不是把判据改成"跟我们一样"。过筛用的是 RDKit 的**另一条**解析路径,
    筛完之后"哪根键、什么顺反、相对谁"照旧逐个比 —— 我们要是给小环双键
    标了顺反,这里就是**多标**,照样红。
    """
    hi = Chem.MolFromSmiles(smi, _KEEP_HS)
    if hi is None or hi.GetNumAtoms() != ref.GetNumAtoms():
        return None
    ok = set()
    for b in hi.GetBonds():
        if b.GetStereo() != Chem.BondStereo.STEREONONE:
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            ok.add((min(i, j), max(i, j)))
    return ok


def reference(smi: str):
    """外部实现的感知结果:{(begin, end): (顺反, 参照a, 参照b)},端点已归一。"""
    m = Chem.MolFromSmiles(smi, sanitize=False)
    if m is None:
        return None
    try:
        Chem.SanitizeMol(m)
        Chem.SetBondStereoFromDirections(m)
    except Exception:
        return None
    ok = stereogenic(smi, m)
    if ok is None:
        return None
    out = {}
    for b in m.GetBonds():
        st = str(b.GetStereo())
        if st not in ("STEREOCIS", "STEREOTRANS"):
            continue
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if (min(i, j), max(i, j)) not in ok:
            continue
        sa = list(b.GetStereoAtoms())
        if len(sa) != 2:
            continue
        out[(min(i, j), max(i, j))] = normalise(m, i, j, st, sa)
    return out


def normalise(mol, i, j, stereo, refs):
    """把 (顺反, 参照对) 归一,使不同的参照选择可以直接比。

    做法:两端各挑**下标最小**的取代基当参照;换掉一端的参照就把顺反翻一次。
    """
    flips = 0
    picked = []
    for end, other, ref in ((i, j, refs[0]), (j, i, refs[1])):
        subs = sorted(
            n.GetIdx() for n in mol.GetAtomWithIdx(end).GetNeighbors() if n.GetIdx() != other
        )
        if ref not in subs:
            return None
        want = subs[0]
        if want != ref:
            flips += 1
        picked.append(want)
    flipped = flips % 2 == 1
    if flipped:
        stereo = "STEREOTRANS" if stereo == "STEREOCIS" else "STEREOCIS"
    # 端点顺序也归一
    if i > j:
        picked.reverse()
    return (stereo, picked[0], picked[1])


def ours(smi: str, cell: str):
    m = Chem.MolFromSmiles(smi, sanitize=False)
    if m is None:
        return None
    try:
        Chem.SanitizeMol(m)
    except Exception:
        return None
    out = {}
    for part in filter(None, cell.split(";")):
        b, e, st, ra, rb = part.split(",")
        b, e, ra, rb = int(b), int(e), int(ra), int(rb)
        st = "STEREOCIS" if st == "CIS" else "STEREOTRANS"
        got = normalise(m, b, e, st, [ra, rb])
        if got is None:
            return "参照原子不是该双键的邻居"
        out[(min(b, e), max(b, e))] = got
    return out


#: 语料里允许有多少条**没真正进比对**。
#:
#: 这条判据先前只数分歧、不数"该数到多少" —— 上游的 `dump_bond_stereo` 少喂
#: 几个分子,每一档都跟着变好看,而它退 0。喂**空文件**进去更彻底:一条分歧
#: 都没有,打印一片空白然后"全部通过"。
#:
#: 实测:`large.smi`(8839 行)没比到 **6**(2022.09.5)/ **8**(2025.09.2,CI 装的),
#: 全是外部实现读不了原串;
#: `smoke.smi`(149 行)没比到 **12**(8 条故意解析不了 + 4 条判官读不了)。
#: 现值 15 = 实测最大加一点余量。
#:
#: 这是分母闸,不是宽容度。涨上去说明有一类分子进不了比对,要当场查。
MAX_UNCOMPARED = 15


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("tsv", type=pathlib.Path)
    ap.add_argument("corpus", type=pathlib.Path, help="喂给 dump_bond_stereo 的那份语料(核分母用)")
    ap.add_argument("--limit", type=int, default=8)
    args = ap.parse_args()

    stat: collections.Counter = collections.Counter()
    bad = []
    rows = 0

    for line in args.tsv.read_text().splitlines():
        if not line.strip():
            continue
        rows += 1
        smi, _, cell = line.partition("\t")
        if cell.startswith("<"):
            stat[cell] += 1
            continue
        want = reference(smi)
        if want is None:
            stat["外部实现读不了"] += 1
            continue
        have = ours(smi, cell)
        if isinstance(have, str):
            stat["参照原子无效"] += 1
            bad.append((smi, have))
            continue
        if have is None:
            stat["外部实现读不了"] += 1
            continue
        if want == have:
            stat["一致"] += 1
            if want:
                stat["  └ 其中确实有标注"] += 1
        else:
            missing = {k: v for k, v in want.items() if k not in have}
            extra = {k: v for k, v in have.items() if k not in want}
            differ = {k: (want[k], have[k]) for k in want.keys() & have.keys() if want[k] != have[k]}
            if missing and not extra and not differ:
                stat["本实现漏标"] += 1
            elif extra and not missing and not differ:
                stat["本实现多标"] += 1
            else:
                stat["顺反判定不同"] += 1
            bad.append((smi, f"缺{missing} 多{extra} 不同{differ}"))

    # **分母闸。** 见 `MAX_UNCOMPARED` 与 `denominator.py`:少比几个分子,
    # 上面每一档都会变好看,而"分歧 0"是它退 0 的唯一依据 ——
    # 喂空文件进去分歧当然是 0。
    n_corpus = denominator.corpus_size(args.corpus)
    compared = stat["一致"] + len(bad)

    print(f"外部实现:RDKit {rdkit.__version__}")
    for k, v in stat.most_common():
        print(f"  {k:<24} {v}")
    print(denominator.line(n_corpus, rows, compared, MAX_UNCOMPARED))
    for s, why in bad[: args.limit]:
        print(f"\n  {s}\n     {why}")
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
