#!/usr/bin/env python3
"""**我方写的 molblock,我方自己读得回来吗** —— 真值由外部实现读同一份字节给出。

# 为什么这条路以前没人走过

已有的两条判据各守一半:

* `check_molblock.py` 比的是"**外部实现**读我方写的块",证明交出去的字节别人
  读得对;
* `check_molblock_read.py` 比的是"我方读**外部实现**写的块",证明别人的字节
  我方读得对。

两条都绿,**我方写的字节我方自己读**这条路仍然一次都没走过 —— 而那是最常用的
一条:写一个 `.mol` 出去,过一会儿再读回来。实测它当时是断的:写出侧给方括号
原子写了价键字段(`vvv=4`),读出侧把这个字段当成"这个原子上没有氢",于是
`C[C@H](N)C(=O)O` 读回来是个三配位的碳 —— 中心少一个配体,手性当场消失,
而分子式看着一点毛病没有。

# 判什么

同一份块,两侧各读一遍,规范串必须相同:

* **我方**:`omgkit.parse_molblock(block)`(Python 绑定,也就是用的人真正走的那条路)
* **外部实现**:`Chem.MolFromMolBlock(block)`

两侧的串都交给外部实现规范化再比 —— 跨实现不能直接比规范串。

# 还比第二件事:**写出去的字节忠不忠于原分子**

上面那一条只保证"两个读者对同一份字节的读法一致" —— 两边一起读出同一句**假话**
时它照样全绿。实测栽过一次:作者没写顺反的双键,我方写图时没标交叉双键,于是
图上那个由布局随手摆出来的几何被两个读者都读成了化学信息。**8831 个分子里
551 个(6.2%)**,而当时 40 道闸一条没红。

`check_molblock.py` 也够不着它:那条判据的口径是"读回来必须**满足**输入指定的
每一处立体",明写着「多出来的立体信息不是错」—— 对三维结构那是对的(桥头碳的
构型确实被坐标定死),对顺反不对:molblock 有表达"未知顺反"的标准写法
(交叉双键),不用它就是造信息。

所以这里加一档:**外部实现从我方写的块里读出的顺反,不得多于它从原串里读出的**。
上限 0。

# 分档

`--max-diff` 那一档是"骨架一样、只有立体不同"。它不是允许错,是给已知的感知
边界留位置,逐条查过才准往上抬。"读成了别的分子"没有上限可言。

**带手性的条数配一条下限**:两侧一起把手性丢光时这条判据同样全绿 —— 那正是
它要抓的故障。

用法:

    python3 harness/check_molblock_roundtrip.py harness/corpus/large.smi
"""
import argparse
import sys

import omgkit
import rdkit
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

# 骨架相同、只有立体不同的条数上限。**8831 条里 1 条,查过:**
#
# 第 987 行 `C[P@H]C[P@@H]CS(=O)(=O)[O-]` —— 我方读回来与**输入 SMILES 逐字符
# 一致**,外部实现读成 `CPCPCS(=O)(=O)[O-]`,把两个磷的构型都丢了。这是三价磷
# 的感知边界:我方认膦的 P 有孤对、四个方向定得下来,外部实现的二维感知不认。
# 分歧的方向是**我方多读出信息**,不是读错 —— 迁就它等于主动丢掉正确的东西。
#
# 与 `check_molblock3d_read.py` 那 16 条里的 4 条三价磷是同一个边界,只是那边
# 走三维坐标、这边走楔形。
MAX_STEREO_DIFF = 1

# 真正比过的条数下限。语料 8831 条,我方画得出二维图的实测远多于此。
MIN_CHECKED = 7000

# **参照侧读得出四面体的条数下限。** 上限为 0 的那一档单看是空断言:
# 两侧一起丢光手性时"逐条一致"照样成立。
MIN_WITH_CHIRAL = 200

# 参照侧读得出顺反的条数下限。手性与顺反是两段独立的代码,各配一条。
MIN_WITH_EZ = 100

# **原串里本来就带顺反的分子条数下限。** "凭空多出顺反"那一档上限是 0,
# 而 0 在语料里一个顺反都没有时同样成立 —— 那时判据什么也没守。实测 366 条。
MIN_EZ_FROM_SOURCE = 200


def canon(smiles):
    m = Chem.MolFromSmiles(smiles)
    return None if m is None else Chem.MolToSmiles(m)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("corpus")
    ap.add_argument("--max-diff", type=int, default=MAX_STEREO_DIFF)
    ap.add_argument("--min-checked", type=int, default=MIN_CHECKED)
    args = ap.parse_args()

    print(f"外部实现:RDKit {rdkit.__version__}")
    print(f"  omgkit wheel:{omgkit.__file__}")

    same = stereo_diff = with_chiral = with_ez = 0
    ours_cannot_write = ref_cannot_read = 0
    invented = 0          # 写出去之后凭空多出来的顺反
    ez_from_source = 0    # 原串里本来就有顺反的分子数(给上面那档配的下限)
    failures = []
    for lineno, line in enumerate(open(args.corpus, encoding="utf-8")):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        smi = line.split("\t")[0].strip()
        try:
            block = omgkit.parse_smiles(smi).to_molblock_2d()
        except ValueError:
            # 画不出二维图的分子不在这条判据的射程内,`check_molblock.py` 管那一档
            ours_cannot_write += 1
            continue

        # **凭空多出的顺反。** 数的是外部实现自己的两次读取:一次读原串、
        # 一次读我方写的块。同一个实现、同一套口径,差出来的就是我方在写出
        # 那一步添的。
        src = Chem.MolFromSmiles(smi)
        if src is not None:
            def n_ez(mol):
                return sum(
                    1 for b in mol.GetBonds()
                    if b.GetStereo() not in (Chem.BondStereo.STEREONONE,
                                             Chem.BondStereo.STEREOANY)
                )
            n_src = n_ez(src)
            if n_src:
                ez_from_source += 1
            written = Chem.MolFromMolBlock(block)
            if written is not None and n_ez(written) > n_src:
                invented += 1
                if len(failures) < 20:
                    failures.append(
                        f"第 {lineno} 行 {smi}:写出去之后多出 "
                        f"{n_ez(written) - n_src} 处顺反 —— "
                        "没写顺反的双键要标成交叉双键"
                    )

        ref = Chem.MolFromMolBlock(block)
        if ref is None:
            # 判官读不了我方写的块 —— 那是 `check_molblock.py` 该报的事,
            # 这里只是没有真值可比
            ref_cannot_read += 1
            continue
        if any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in ref.GetAtoms()):
            with_chiral += 1
        if any(
            b.GetStereo() not in (Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY)
            for b in ref.GetBonds()
        ):
            with_ez += 1
        want = Chem.MolToSmiles(ref)

        try:
            got = omgkit.parse_molblock(block).mol.to_smiles()
        except ValueError as e:
            failures.append(f"第 {lineno} 行 {smi}:外部实现读得出,我方读不了({e})")
            continue
        got_canon = canon(got)
        if got_canon is None:
            failures.append(f"第 {lineno} 行 {smi}:我方读出的 `{got}` 外部实现读不了")
            continue
        if got_canon == want:
            same += 1
            continue
        flat = Chem.MolToSmiles(Chem.MolFromSmiles(got), isomericSmiles=False)
        if flat == Chem.MolToSmiles(ref, isomericSmiles=False):
            stereo_diff += 1
            failures.append(f"第 {lineno} 行 {smi}:骨架对、立体不同 —— 我方 {got_canon},外部 {want}")
            continue
        failures.append(f"第 {lineno} 行 {smi}:我方读成 {got_canon},外部实现读成 {want}")

    hard = [f for f in failures if "立体不同" not in f and "多出" not in f]
    print(f"逐条一致 {same};骨架对但立体不同 {stereo_diff}(上限 {args.max_diff});"
          f"读成别的分子 {len(hard)}")
    print(f"  我方画不出二维图 {ours_cannot_write};外部实现读不了我方的块 {ref_cannot_read}")
    print(f"  参照侧带四面体的 {with_chiral} 条(下限 {MIN_WITH_CHIRAL});"
          f"带顺反的 {with_ez} 条(下限 {MIN_WITH_EZ})")
    print(f"  **凭空多出顺反的 {invented} 条(上限 0)**;原串本来就带顺反的 "
          f"{ez_from_source} 条(下限 {MIN_EZ_FROM_SOURCE})")
    for f in failures[:8]:
        print(f"  ✗ {f}")

    if hard:
        print("\n我方写出去的 molblock,我方自己读回来不是同一个分子。")
        return 1
    if invented:
        print(f"\n{invented} 条写出去之后多出了顺反 —— 图上那个几何是布局随手摆的,"
              "不是化学信息。没写顺反的双键要标成交叉双键(键块第四列 3)。")
        return 1
    if ez_from_source < MIN_EZ_FROM_SOURCE:
        print(f"\n原串带顺反的只有 {ez_from_source} 条,低于下限 {MIN_EZ_FROM_SOURCE} —— "
              "上面那条上限为 0 的判据被喂空了")
        return 1
    if stereo_diff > args.max_diff:
        print(f"\n只有立体不同的涨到 {stereo_diff} 条,超过上限 {args.max_diff}")
        return 1
    if same + stereo_diff < args.min_checked:
        print(f"\n只比过 {same + stereo_diff} 条,低于下限 {args.min_checked} —— 判据被喂空了")
        return 1
    if with_chiral < MIN_WITH_CHIRAL:
        print(f"\n参照侧只有 {with_chiral} 条带四面体,低于下限 {MIN_WITH_CHIRAL} —— "
              "手性那一档被喂空了,两侧一起丢光也是这个样子")
        return 1
    if with_ez < MIN_WITH_EZ:
        print(f"\n参照侧只有 {with_ez} 条带顺反,低于下限 {MIN_WITH_EZ} —— 顺反那一档被喂空了")
        return 1
    print("\n写出去再读回来是同一个分子。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
