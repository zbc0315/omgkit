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
    python3 harness/check_bond_stereo.py /tmp/bs.tsv
"""

import argparse
import collections
import pathlib
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


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
    out = {}
    for b in m.GetBonds():
        st = str(b.GetStereo())
        if st not in ("STEREOCIS", "STEREOTRANS"):
            continue
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
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


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("tsv", type=pathlib.Path)
    ap.add_argument("--limit", type=int, default=8)
    args = ap.parse_args()

    stat: collections.Counter = collections.Counter()
    bad = []

    for line in args.tsv.read_text().splitlines():
        if not line.strip():
            continue
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

    for k, v in stat.most_common():
        print(f"  {k:<24} {v}")
    for s, why in bad[: args.limit]:
        print(f"\n  {s}\n     {why}")
    if len(bad) > args.limit:
        print(f"\n  ...(另有 {len(bad) - args.limit} 条)")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
