#!/usr/bin/env python3
"""读 molblock:外部实现写的文件,我们读出来是不是同一个分子。

# 两边读的是同一批字节

外部实现把语料里的每条 SMILES 写成 V2000 molblock(二维坐标),我们读那个文件、
写回 SMILES;外部实现也读同一个文件、写回 SMILES。两串都交给它规范化再比 ——
**跨实现不能直接比规范串**(两套规范化算法给出的字符串本来就不一样),
所以统一由外部实现来判"是不是同一个分子"。

差别因此只可能来自**读**:元素符号、键类型、电荷(`M CHG` 与原子块旧字段的
优先级)、同位素(质量差 → 质量数)、自由基(`M RAD` 的编码)、价键字段
(不读的话 `[CH]` 会被按默认价补成 `[CH3]`)。

# 立体化学分三档比

一档"完全一致",一档"只差双键顺反",一档"骨架对、四面体不同"。分开是因为
"顺反读错了"与"手性读错了"是两件事,混成一档的话前者会把后者盖住。

两档各有上限(`MAX_EZ_ONLY` / `MAX_CHIRAL_DIFF`),而且**两个方向都卡** ——
少读一个、多读一个、读反一个,一律落进对应那一档。

上限为 0 的那一档单看是空断言,所以参照侧"有多少条带立体"也一起打出来并配
下限:参照侧一条都没有时,零分歧说明不了任何事。

用法:

    python3 harness/check_molblock_read.py --write  <out.sdf> <语料.smi>
    cargo run -q -p omgkit-io --release --example read_molblock -- <out.sdf> > <ours.txt>
    python3 harness/check_molblock_read.py --compare <out.sdf> <ours.txt>
"""
import argparse
import sys

import rdkit
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

# 外部实现写成 V3000、我方拒收的条数上限。
#
# 实测 `large.smi` 上是 1 条:二茂铁 `…[Fe++]23456789…`,铁上十根键,
# 外部实现自己就换成了 V3000。**贴着现值留了余量,不是随手写个大数** ——
# 这一档一旦变多,要么语料变了,要么该把 V3000 也读起来了,两种都得有人看见。
MAX_V3000 = 3

# **只差双键顺反**的条数上限。**两个方向都卡在这里**:我方少读一根、多读一根、
# 或者读反一根,三种都落进这一档(`no_ez` 把两侧的顺反一起抹掉再比,剩下的
# 差别只可能来自顺反本身)。
#
# 先前是 365 —— 那时我方根本不从坐标反读顺反。现在读了,实测 0。
MAX_EZ_ONLY = 0

# 参照侧**认得出顺反**的分子条数下限。
#
# `MAX_EZ_ONLY = 0` 单看是个空断言:参照侧一根顺反都没有的话,这一档永远是 0,
# 而判据会照常打印"逐条一致"。实测 366 条 —— 贴着现值留了余量。
MIN_EZ_MOLECULES = 300

# **四面体也不一样**的条数上限。
#
# 实测 5 条,全是桥环/稠环,而且**两个方向都有**:
#
# * 3 条我方少读一个中心 —— 二维布局退化时三个画出来的邻居张不出足够的体积,
#   我方按设计判"这张图定不出构型"(见 `omgkit_io::wedge` 的 `ZERO_VOLUME_TOL`);
# * 2 条我方多读一个 —— 桥头碳的构型由环系本身定死,外部实现的立体感知把这类
#   "不是真立体中心"的标记摘掉,我方读的是纯几何。
#
# 所以这不是"我方保守",是**立体感知的边界不一样**。钉住是为了让它不悄悄变多;
# 真要收敛得先把两边的"什么算真立体中心"对齐,那是另一块活。
MAX_CHIRAL_DIFF = 5


def write_blocks(out_path, corpus, limit):
    n = 0
    with open(out_path, "w", encoding="utf-8") as f:
        for lineno, line in enumerate(open(corpus, encoding="utf-8")):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            smi = line.split("\t")[0].strip()
            m = Chem.MolFromSmiles(smi)
            if m is None:
                continue
            AllChem.Compute2DCoords(m)
            f.write(f">>> {lineno}\t{smi}\n")
            f.write(Chem.MolToMolBlock(m))
            f.write("$$$$\n")
            n += 1
            if limit and n >= limit:
                break
    print(f"写了 {n} 条 molblock 到 {out_path}")
    return n


def blocks(path):
    smi = None
    lineno = None
    buf = []
    for line in open(path, encoding="utf-8"):
        if line.startswith(">>> "):
            lineno, smi = line[4:].rstrip("\n").split("\t", 1)
            buf = []
        elif line.rstrip("\n") == "$$$$":
            yield lineno, smi, "".join(buf)
        else:
            buf.append(line)


def no_ez(smiles):
    """只抹掉**双键顺反**、留着四面体的规范串。读不了时给 None。

    分档要分得开:"顺反没读"与"四面体读错了"是两件事,混在一档里的话前者会把
    后者盖住 —— 而后者才是真的读错。
    """
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    for b in m.GetBonds():
        b.SetStereo(Chem.BondStereo.STEREONONE)
        b.SetBondDir(Chem.BondDir.NONE)
    try:
        m = Chem.RemoveHs(m)
    except Exception:  # noqa: BLE001
        return None
    return Chem.MolToSmiles(m)


def canon(smiles, drop_stereo):
    """规范串。`drop_stereo` 为真时先抹掉立体。读不了时给 None。

    **去显式氢那一步不能省。** 外部实现的 SMILES 解析器会保留"承载方向键的
    显式氢"(`[H]/[O+]=c1…` 里那个 H 是一个原子),而两侧未必都写方向键 ——
    同一个氢一边是原子、一边被合并进邻居,原子数一个 9 一个 8,比出来是
    "不同的分子",而那是**比法**造成的,不是读法。实测 8831 条里 26 条栽在这上面。
    """
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    if drop_stereo:
        Chem.RemoveStereochemistry(m)
    try:
        m = Chem.RemoveHs(m)
    except Exception:  # noqa: BLE001
        return None
    return Chem.MolToSmiles(m)


def compare(blocks_path, ours_path, min_checked):
    print(f"外部实现:RDKit {rdkit.__version__}")
    ours = {}
    for line in open(ours_path, encoding="utf-8"):
        lineno, _smi, got = line.rstrip("\n").split("\t", 2)
        ours[lineno] = got

    same = diff = unreadable = skipped = with_stereo = v3000 = 0
    ez_only = chiral_diff = with_ez = 0
    failures = []
    for lineno, smi, block in blocks(blocks_path):
        got = ours.get(lineno)
        if got is None:
            failures.append(f"第 {lineno} 行 {smi}:我方一行输出都没有")
            unreadable += 1
            continue
        if got.startswith("<"):
            # **V3000 是一档有意的限制,不是失败。** 本模块只读 V2000;认出来
            # 明确拒收,好过把它当 V2000 硬读出一个错分子。单独计数、配上限,
            # 免得这一档悄悄长大。
            if "V3000" in got:
                v3000 += 1
            else:
                failures.append(f"第 {lineno} 行 {smi}:{got}")
                unreadable += 1
            continue
        ref = Chem.MolFromMolBlock(block)
        if ref is None:
            skipped += 1
            continue
        if any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in ref.GetAtoms()):
            with_stereo += 1
        # **`STEREOANY` 不算。** 交叉双键(键块第四列 `3`)在参照侧读成
        # `STEREOANY` —— 那是"作者说不知道",两侧都不该有构型。把它算进来的话,
        # 只含交叉双键的分子会把下限撑起来,而顺反那一档照样是空的。
        if any(
            b.GetStereo() not in (Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY)
            for b in ref.GetBonds()
        ):
            with_ez += 1
        raw = Chem.MolToSmiles(ref)
        want = canon(raw, drop_stereo=False)
        mine = canon(got, drop_stereo=False)
        if mine is None:
            unreadable += 1
            failures.append(f"第 {lineno} 行 {smi}:我方写出的 `{got}` 外部实现读不了")
        elif mine == want:
            same += 1
        elif no_ez(got) == no_ez(raw):
            # 只差双键顺反 —— 少读、多读、读反,三种都在这里
            ez_only += 1
            failures.append(
                f"第 {lineno} 行 {smi}:骨架对、顺反不同 —— 我方 {mine},外部 {want}"
            )
        elif canon(got, drop_stereo=True) == canon(raw, drop_stereo=True):
            # 骨架一样,四面体也不同 —— 桥环那一档
            chiral_diff += 1
            failures.append(
                f"第 {lineno} 行 {smi}:骨架对、四面体不同 —— 我方 {mine},外部 {want}"
            )
        else:
            diff += 1
            failures.append(f"第 {lineno} 行 {smi}:我方读成 {mine},外部实现读成 {want}")

    print(f"读回来一致 {same};不一致 {diff};读不了/写不出 {unreadable};"
          f"外部实现自己读不了 {skipped}")
    print(f"  外部实现写成了 V3000、我方明确拒收的 {v3000} 条(上限 {MAX_V3000})")
    print(f"  只差双键顺反 {ez_only} 条(上限 {MAX_EZ_ONLY});"
          f"四面体也不同 {chiral_diff} 条(上限 {MAX_CHIRAL_DIFF})")
    print(f"  参照侧带四面体的 {with_stereo} 条;带顺反的 {with_ez} 条"
          f"(下限 {MIN_EZ_MOLECULES})")
    if failures:
        for f in failures[:8]:
            print(f"  ✗ {f}")
    if diff or unreadable:
        print("\n读出来不是同一个分子。")
        return 1
    for got_n, cap, what in [
        (ez_only, MAX_EZ_ONLY, "只差双键顺反"),
        (chiral_diff, MAX_CHIRAL_DIFF, "四面体也不同"),
    ]:
        if got_n > cap:
            print(f"\n{what}的涨到 {got_n} 条,超过上限 {cap}")
            return 1
    if v3000 > MAX_V3000:
        print(f"\nV3000 那一档涨到 {v3000} 条,超过上限 {MAX_V3000} —— "
              "要么语料变了,要么该把 V3000 也读起来了")
        return 1
    if with_ez < MIN_EZ_MOLECULES:
        print(f"\n参照侧只有 {with_ez} 条带顺反,低于下限 {MIN_EZ_MOLECULES} —— "
              "顺反那一档被喂空了,它的零分歧说明不了任何事")
        return 1
    if same < min_checked:
        print(f"\n只比过 {same} 条,低于下限 {min_checked} —— 判据被喂空了")
        return 1
    print("\n读回来逐条一致(立体除外)。")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", nargs=2, metavar=("OUT", "CORPUS"))
    ap.add_argument("--compare", nargs=2, metavar=("BLOCKS", "OURS"))
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--min-checked", type=int, default=8000)
    args = ap.parse_args()
    if args.write:
        write_blocks(args.write[0], args.write[1], args.limit)
        return 0
    if args.compare:
        return compare(args.compare[0], args.compare[1], args.min_checked)
    ap.print_help()
    return 2


if __name__ == "__main__":
    sys.exit(main())
