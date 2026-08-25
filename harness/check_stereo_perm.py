#!/usr/bin/env python3
"""配位几何(`@SP`/`@TB`/`@OH`)的排列序号:与外部实现**分组一致**。

# 判据的形状

序号的含义相对"配体按什么顺序列出"。对每一类几何,取**互不相同**的配体,
穷举"每个序号 × 每种列出顺序"的全部写法:

| 类别 | 配体数 | 序号数 | 写法数 |
|---|---|---|---|
| `@SP` | 4 | 3 | 3 × 24 = 72 |
| `@TB` | 5 | 20 | 20 × 120 = 2400 |
| `@OH` | 6 | 30 | 30 × 720 = 21600 |

两侧各自规范化,再比**分组**(哪些写法落到同一个分子)。比分组而不是比字符串:
两个实现的规范串本来就不一样,一致的应当是"谁和谁是同一个分子"这件事。

# 为什么必须用互不相同的配体

配体两两相同的话,不同序号会指向同一个分子(那是对的:`[Pt@SP1](Cl)(Cl)(N)N`
与 `[Pt@SP3](Cl)(Cl)(N)N` 就是同一个异构体),分组会塌,规则本身就被盖住了。
顺带一提,那种塌陷正是我方与参照的一处**已知差别**:我方把它收敛成一串,
RDKit 不收(把同一个分子重排原子再规范化,它会在两个序号之间跳)。
所以这条判据只跑互不相同的配体 —— 那一档没有自同构,两侧应当逐组吻合。

# 这条判据守的是什么

`crates/omgkit-core/src/polyhedron.rs` 里那三张表就是这么量出来的。表进了源码
之后,**没有任何东西再核对它** —— 这条判据每次 CI 都重新量一遍。

用法:python3 harness/check_stereo_perm.py
"""
import itertools
import sys

import omgkit
from rdkit import Chem

CASES = [
    ("SP", "Pt", ["F", "Cl", "Br", "I"], 3),
    ("TB", "P", ["F", "Cl", "Br", "I", "S"], 20),
    ("OH", "Co", ["F", "Cl", "Br", "I", "S", "P"], 30),
]


def partition(by_key):
    """{写法: 规范串} → {frozenset(写法)},即"谁和谁是同一个分子"。"""
    groups = {}
    for key, value in by_key.items():
        groups.setdefault(value, set()).add(key)
    return {frozenset(v) for v in groups.values()}


def main() -> int:
    print(f"  外部实现:RDKit {Chem.rdBase.rdkitVersion};omgkit wheel:{omgkit.__file__}")
    bad = 0
    for tag, elem, ligands, n_perm in CASES:
        ours, theirs = {}, {}
        for idx in range(1, n_perm + 1):
            for order in itertools.permutations(range(len(ligands))):
                lig = [ligands[i] for i in order]
                smi = f"[{elem}@{tag}{idx}](" + ")(".join(lig[:-1]) + ")" + lig[-1]
                m = Chem.MolFromSmiles(smi)
                if m is None:
                    print(f"✗ {tag}:RDKit 读不了 {smi}")
                    bad += 1
                    continue
                theirs[(idx, order)] = Chem.MolToSmiles(m)
                ours[(idx, order)] = omgkit.parse_smiles(smi).to_canonical_smiles()

        n = len(ours)
        want = n_perm * len(list(itertools.permutations(range(len(ligands)))))
        if n != want:
            print(f"✗ {tag}:只比了 {n} 种写法,应当是 {want} —— 判据被喂空了一部分")
            bad += 1
            continue
        pa, pb = partition(ours), partition(theirs)
        if len(pb) != n_perm:
            print(f"✗ {tag}:参照把 {n} 种写法归成 {len(pb)} 组,应当是 {n_perm} 组")
            bad += 1
            continue
        if pa != pb:
            print(f"✗ {tag}:分组不一致 —— 我方 {len(pa)} 组,参照 {len(pb)} 组")
            only_ours = pa - pb
            for g in list(only_ours)[:2]:
                print(f"    我方有而参照没有的一组({len(g)} 种写法):{sorted(g)[:4]} …")
            bad += 1
            continue
        print(f"✓ {tag}:{n} 种写法,两侧都归成 {n_perm} 组,逐组吻合")

    if bad:
        print(f"\n{bad} 类的排列序号与外部实现对不上 —— "
              "`polyhedron.rs` 的表或者读写两侧的换算有问题。")
        return 1
    print("\n三类配位几何的排列分组与外部实现完全一致。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
