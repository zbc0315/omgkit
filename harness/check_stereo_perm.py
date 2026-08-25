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

# 缺一个顶点的那几族

方括号里的氢、或者一个空的配位位置,也占多面体的一个顶点,而它不在键序列里。
它在列出顺序里落在哪一位是量出来的:"自身位置" —— 紧跟前驱原子之后、
环闭合之前。下面按上下文分族压住这条规则:手性原子居首 / 前面有原子 /
中心上带环闭合,空位 / 方括号里的氢。放错位置时解析与写出两侧会一起错、
正好抵消,往返判据看不出来 —— 只有这里的分组比对能看出来。

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

def branches(elem, tag, extra=""):
    """`[X@TAG{i}](a)(b)…z` —— 手性原子居首,配体全是分支。"""
    return lambda idx, lig: (
        f"[{elem}@{tag}{idx}{extra}](" + ")(".join(lig[:-1]) + ")" + lig[-1]
    )


def with_parent(elem, tag, extra=""):
    """`a[X@TAG{i}](b)…z` —— 第一个配体写在手性原子**前面**。"""
    return lambda idx, lig: (
        f"{lig[0]}[{elem}@{tag}{idx}{extra}](" + ")(".join(lig[1:-1]) + ")" + lig[-1]
    )


def ring(elem, tag):
    """`[X@TAG{i}]1(a)…CO1` —— 中心上带一个环闭合(缺的顶点排在它**前面**)。"""
    return lambda idx, lig: f"[{elem}@{tag}{idx}]1(" + ")(".join(lig) + ")CO1"


def ring_with_parent(elem, tag):
    """`a[X@TAG{i}]1(b)…CO1` —— 同一个分子,换个起笔位置。"""
    return lambda idx, lig: (
        f"{lig[0]}[{elem}@{tag}{idx}]1(" + ")(".join(lig[1:]) + ")CO1"
    )


def ring_rev(elem, tag):
    """`[X@TAG{i}]1(a)…OC1` —— 同一个三元环**反着写**,环闭合的那一端从 O 换成 C。"""
    return lambda idx, lig: f"[{elem}@{tag}{idx}]1(" + ")(".join(lig) + ")OC1"


def ring_rev_with_parent(elem, tag):
    """`a[X@TAG{i}]1(b)…OC1` —— 反着写,而且换个起笔位置。"""
    return lambda idx, lig: (
        f"{lig[0]}[{elem}@{tag}{idx}]1(" + ")(".join(lig[1:]) + ")OC1"
    )


# 每一族里的几种"形状"必须写的是**同一批分子**,只是写法不同。
#
# 一族里只放一种形状是不够的。分组判据看不见"给配体全局换个名字"这件事:
# 缺的那个顶点与它旁边那个配体如果在全族里都固定,把两者对调,每一条写法都
# 挨同样的一下,分组一点不变。变异实测(两侧一起把缺的顶点挪到环闭合**之后**):
# 只有 `ring` 一种形状时全绿;加上 `ring_with_parent` **仍然全绿**(幻影还是
# 紧挨着同一个环氧);要把三元环**反着写**一遍、让环闭合的那一端在 O 与 C
# 之间变,这条变异才红。
#
# 不带环的那几族没有这个问题:幻影旁边是前驱原子,而哪个配体当前驱原子
# 随列出顺序在变。
CASES = [
    # (族名, 序号数, 配体, 几种写法形状)
    ("SP·满配位", 3, ["F", "Cl", "Br", "I"], [branches("Pt", "SP")]),
    ("TB·满配位", 20, ["F", "Cl", "Br", "I", "S"], [branches("P", "TB")]),
    ("OH·满配位", 30, ["F", "Cl", "Br", "I", "S", "P"], [branches("Co", "OH")]),
    # 缺一个顶点:一个空的配位位置
    ("SP·空位", 3, ["F", "Cl", "Br"], [branches("Pt", "SP"), with_parent("Pt", "SP")]),
    (
        "TB·空位",
        20,
        ["F", "Cl", "Br", "I"],
        [branches("P", "TB"), with_parent("P", "TB")],
    ),
    (
        "OH·空位",
        30,
        ["F", "Cl", "Br", "I", "S"],
        [branches("Co", "OH"), with_parent("Co", "OH")],
    ),
    # 缺一个顶点:方括号里写了氢
    (
        "SP·方括号里的氢",
        3,
        ["F", "Cl", "Br"],
        [branches("Pt", "SP", "H"), with_parent("Pt", "SP", "H")],
    ),
    (
        "OH·方括号里的氢",
        30,
        ["F", "Cl", "Br", "I", "S"],
        [branches("Co", "OH", "H"), with_parent("Co", "OH", "H")],
    ),
    # 缺一个顶点,而且中心上有环闭合 —— 那个顶点排在环闭合**之前**。
    # 环里两个原子取 C 与 O:都是碳的话交换它们是这个分子的自同构,分组会塌。
    # 平面四方那一档没放进来:配体只剩一个分支,两种形状之间排不出差别。
    (
        "OH·环闭合",
        30,
        ["F", "Cl", "Br"],
        [
            ring("Co", "OH"),
            ring_with_parent("Co", "OH"),
            ring_rev("Co", "OH"),
            ring_rev_with_parent("Co", "OH"),
        ],
    ),
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
    for family, n_perm, ligands, shapes in CASES:
        ours, theirs = {}, {}
        for idx in range(1, n_perm + 1):
            for order in itertools.permutations(range(len(ligands))):
                for shape, build in enumerate(shapes):
                    smi = build(idx, [ligands[i] for i in order])
                    m = Chem.MolFromSmiles(smi)
                    if m is None:
                        print(f"✗ {family}:RDKit 读不了 {smi}")
                        bad += 1
                        continue
                    theirs[(idx, order, shape)] = Chem.MolToSmiles(m)
                    ours[(idx, order, shape)] = omgkit.parse_smiles(
                        smi
                    ).to_canonical_smiles()

        n = len(ours)
        want = (
            n_perm
            * len(list(itertools.permutations(range(len(ligands)))))
            * len(shapes)
        )
        if n != want:
            print(f"✗ {family}:只比了 {n} 种写法,应当是 {want} —— 判据被喂空了一部分")
            bad += 1
            continue
        pa, pb = partition(ours), partition(theirs)
        if len(pb) != n_perm:
            print(f"✗ {family}:参照把 {n} 种写法归成 {len(pb)} 组,应当是 {n_perm} 组")
            bad += 1
            continue
        if pa != pb:
            print(f"✗ {family}:分组不一致 —— 我方 {len(pa)} 组,参照 {len(pb)} 组")
            only_ours = pa - pb
            for g in list(only_ours)[:2]:
                print(f"    我方有而参照没有的一组({len(g)} 种写法):{sorted(g)[:4]} …")
            bad += 1
            continue
        print(f"✓ {family}:{n} 种写法,两侧都归成 {n_perm} 组,逐组吻合")

    if bad:
        print(f"\n{bad} 族的排列序号与外部实现对不上 —— "
              "`polyhedron.rs` 的表或者读写两侧的换算有问题。")
        return 1
    print(f"\n{len(CASES)} 族配位几何写法的排列分组与外部实现完全一致。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
