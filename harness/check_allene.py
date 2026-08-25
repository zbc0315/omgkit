#!/usr/bin/env python3
"""丙二烯型轴手性(`@AL`):与外部实现**分组一致**。

# 为什么这条判据不用 RDKit

**仓库钉的 RDKit 2025.09.2 完全不支持丙二烯立体。** 逐条实测:SMILES 读写、
`FindPotentialStereo`、`rdCIPLabeler`、从 3D 坐标反推、molblock 往返、带手性的
子结构匹配 —— 六条路全都把 `@AL1` 与 `@AL2` 当成同一个东西,连
`F/C=C=C=C/F` 的方向键也丢。所以这一档它当不了裁判。

Indigo 支持:`exactMatch(..., "ALL")` 认这个区别,而且自洽 ——
`@AL1` ≡ `@`、`@AL2` ≡ `@@`,从另一端起笔序号不变,读→写链上每一步都还是
同一个分子。它的**规范串**倒是把立体丢了,所以这里比的是 `exactMatch` 给出的
分组,不是字符串。

# 判据的形状

同一个分子的全部写法 × 两个序号,两侧各自分组,再比分组。自由度:
先写哪一端、每端内部两个取代基谁在前、中心是不是写在最前、环反着写。

**一族里必须有好几种写法。** 比分组的判据看不见"给配体全局换个名字":
若某两个配体的角色在全族里绑死,把它俩对调,每条写法都挨同样一下,分组一点
不变 —— 上一块(配位几何)实测过这个洞。所以这里每族都换起笔位置、
换环的书写方向。

用法:python3 harness/check_allene.py
"""
import itertools
import sys

import omgkit
from indigo import Indigo

IND = Indigo()

ENDS = [("N", "Br"), ("O", "F")]


def substituted():
    """两端各两个取代基,四个配体互不相同。"""
    out = []
    for idx, first, s0, s1, centre_first in itertools.product(
        (1, 2), (0, 1), (0, 1), (0, 1), (0, 1)
    ):
        e1 = list(ENDS[first])
        e2 = list(ENDS[1 - first])
        if s0:
            e1.reverse()
        if s1:
            e2.reverse()
        if centre_first:
            out.append(f"[C@AL{idx}](=C({e1[0]}){e1[1]})=C({e2[0]}){e2[1]}")
        else:
            out.append(f"{e1[0]}C({e1[1]})=[C@AL{idx}]=C({e2[0]}){e2[1]}")
    return out


def with_hydrogen():
    """一端只写一个取代基,另一个配体是那个原子上的氢。

    氢落在"自身位置"(紧跟前驱原子之后),裸写与方括号两种形式都要走一遍。
    """
    out = []
    for idx, s1, centre_first, bracket in itertools.product(
        (1, 2), (0, 1), (0, 1), (0, 1)
    ):
        e2 = list(ENDS[1])
        if s1:
            e2.reverse()
        end = "[CH]" if bracket else "C"
        if centre_first:
            out.append(f"[C@AL{idx}](={end}N)=C({e2[0]}){e2[1]}")
        else:
            out.append(f"N{end}=[C@AL{idx}]=C({e2[0]}){e2[1]}")
    return out


def with_ring():
    """一端的两个配体都在一个三元环里 —— 环还要**反着写**一遍。

    反着写会把环闭合的那一端从 O 换成 C。少了这一手,"缺的那个配体"
    与它旁边那个在全族里绑死,分组判据就看不见两者对调。
    """
    out = []
    for idx, (x, y) in itertools.product((1, 2), (("N", "Br"), ("Br", "N"))):
        tail = f"=[C@AL{idx}]=C({x}){y}"
        head = f"(=[C@AL{idx}]=C({x}){y})"
        out += [
            f"C1{head}OC1",  # 环闭合在自身位置,链上的邻居是 O
            f"C1{head}CO1",  # 环反着写:链上的邻居换成 C
            f"O1CC1{tail}",  # 换起笔位置:前驱是 C,环闭合到 O
            f"C1OC1{tail}",  # 再反着写一遍
        ]
    return out


FAMILIES = [
    ("两端各两个取代基", substituted, 2),
    ("一端带一个氢", with_hydrogen, 2),
    ("一端的配体在环里", with_ring, 2),
]


def strip_stereo(smi):
    return smi.replace("@AL1", "").replace("@AL2", "")


def indigo_partition(smis):
    """按 `exactMatch(..., "ALL")` 分组。"""
    groups = []
    for smi in smis:
        m = IND.loadMolecule(smi)
        for g in groups:
            if IND.exactMatch(m, IND.loadMolecule(g[0]), "ALL") is not None:
                g.append(smi)
                break
        else:
            groups.append([smi])
    return {frozenset(g) for g in groups}


def ours_partition(smis):
    by_canonical = {}
    for smi in smis:
        key = omgkit.parse_smiles(smi).to_canonical_smiles()
        by_canonical.setdefault(key, set()).add(smi)
    return {frozenset(v) for v in by_canonical.values()}


def main() -> int:
    print(f"  外部实现:Indigo {IND.version()};omgkit wheel:{omgkit.__file__}")
    bad = 0
    for name, build, want_groups in FAMILIES:
        smis = build()
        if len(set(smis)) != len(smis):
            print(f"✗ {name}:{len(smis)} 条写法里有重复,判据被喂窄了")
            bad += 1
            continue

        # 这一族必须是**同一个分子**的不同写法 —— 否则比的不是同一件事
        plain = {IND.loadMolecule(strip_stereo(s)).canonicalSmiles() for s in smis}
        if len(plain) != 1:
            print(f"✗ {name}:去掉立体之后不是同一个分子,{sorted(plain)}")
            bad += 1
            continue

        theirs = indigo_partition(smis)
        if len(theirs) != want_groups:
            print(f"✗ {name}:外部实现把 {len(smis)} 条归成 {len(theirs)} 组,应当是 {want_groups} 组")
            bad += 1
            continue

        ours = ours_partition(smis)
        if ours != theirs:
            print(f"✗ {name}:分组不一致 —— 我方 {len(ours)} 组,外部 {len(theirs)} 组")
            for g in list(ours - theirs)[:2]:
                print(f"    我方有而外部没有的一组({len(g)} 条):{sorted(g)[:3]} …")
            bad += 1
            continue

        print(f"✓ {name}:{len(smis)} 条写法,两侧都归成 {want_groups} 组,逐组吻合")

    if bad:
        print(f"\n{bad} 族的丙二烯轴手性与外部实现对不上。")
        return 1
    print(f"\n{len(FAMILIES)} 族丙二烯写法的分组与外部实现完全一致。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
