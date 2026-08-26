#!/usr/bin/env python3
"""读**三维** molblock:外部实现写的三维结构,我们从坐标读出来是不是同一个分子。

# 与二维那条判据是两条路

`check_molblock_read.py` 喂的是二维图:手性靠**楔形**、顺反靠**平面投影**。
三维文件里楔形一般是空的,立体全在坐标里:手性靠**有符号体积**、顺反靠
**二面角**。两条路的代码没有一行是共用的,一条绿说明不了另一条。

# 真值从哪来

外部实现把每条 SMILES 嵌成三维构象(固定随机种子)、写成三维 molblock,
**再自己从那个文件读回来**(`AssignStereochemistryFrom3D`)。我方读同一份字节。
两侧各写回 SMILES,统一交给外部实现规范化再比 —— 跨实现不能直接比规范串。

真值不是原始 SMILES:嵌出来的构象**可能与原串的立体不符**(嵌入失败、
或者构象本身就把某个中心摆成了对映体)。拿原串当真值等于把嵌入的毛病算到
读取头上。所以真值是"外部实现自己从这份坐标读出来的东西"。

# 分档

`MAX_STEREO_DIFF` 那一档单独计:骨架一样、只有立体不同。实测 16 条,逐条查过,
只有两类,都是**能力/感知边界**而不是读错 —— 明细写在那个常数上面。

"读成了别的分子"那一档是 0,而且没有上限可言:骨架读错就是读错。

用法:

    python3 harness/check_molblock3d_read.py --write  <out.sdf> <语料.smi>
    cargo run -q -p omgkit-io --release --example read_sdf -- <out.sdf> > <ours.txt>
    python3 harness/check_molblock3d_read.py --compare <out.sdf> <ours.txt>
"""
import argparse
import sys

import rdkit
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

# 骨架相同、只有立体不同的条数上限。**8795 条里 16 条,只有两类,逐条查过:**
#
# * **12 条:外部实现读了非四面体立体**(`@TB` / `@OH` / `@SP`)—— 六氨合钴那一族、
#   八面体的镍、三角双锥的钒与磷。我方眼下只从三维坐标读四面体,这一档是**能力
#   空白**,不是读错;
# * **4 条:三价磷**。我方认它是立体中心(膦、亚磷酸酯的 P 有孤对,构型确定),
#   外部实现的三维感知不认。其中两条 P 上挂着两条组成相同的支路(靠支路里各自的
#   中心才成立),另两条是货真价实的膦中心 —— 两个方向都有,是**感知边界不一样**。
#
# 钉住是为了让它不悄悄变多。要收敛得先把非四面体那一档做出来,那是另一块活。
MAX_STEREO_DIFF = 16

# 真正比过的条数下限。语料 8831 条里外部实现嵌得出三维构象的实测 8795 条,
# 其余是它自己的界矩阵过不去。**这一档要有下限**:嵌入那一步半路断掉、
# 或者换个版本嵌得更少,判据都会在一份悄悄变短的文件上报"逐条一致"。
MIN_RECORDS = 7000

# 参照侧**读得出四面体**的分子条数下限。上限为 0 的那一档单看是空断言。
MIN_WITH_CHIRAL = 200

# 参照侧**读得出顺反**的分子条数下限。手性与顺反是两段独立的代码,
# 各配一条 —— 合成一条的话,顺反整档丢光会被手性那一档撑住。
MIN_WITH_EZ = 100


def write_sdf(out_path, corpus, limit):
    """嵌三维构象、写三维 molblock。

    **补氢再嵌。** 不补的话嵌出来的构型在手性中心上是欠定的,外部实现自己都
    读不回原来的标记 —— 那时判据量的是嵌入的毛病,不是读取的。
    """
    n = failed = 0
    with Chem.SDWriter(out_path) as w:
        for lineno, line in enumerate(open(corpus, encoding="utf-8")):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            smi = line.split("\t")[0].strip()
            m = Chem.MolFromSmiles(smi)
            if m is None:
                continue
            m = Chem.AddHs(m)
            ps = AllChem.ETKDGv3()
            ps.randomSeed = 0xF00D
            # **嵌不出来有两种表现:返回 −1,以及直接抛。**
            #
            # 后者是 `Invariant Violation: bad lower bound` —— 界矩阵自相矛盾时
            # RDKit 直接抛 `RuntimeError`,不走返回码。不接住的话整个 `--write`
            # 半路断掉,而写出去的文件是**截断**的:头一版就这么写的,全量跑
            # 只写了 7950 条就停了,而我当时在后面接了 `tail -1`,一点没看见。
            # 两种都计进 `failed` 并打出来,免得分母悄悄变小。
            try:
                ok = AllChem.EmbedMolecule(m, ps) == 0
            except RuntimeError:
                ok = False
            if not ok:
                failed += 1
                continue
            m.SetProp("_Name", f"第{lineno}条")
            w.write(m)
            n += 1
            if limit and n >= limit:
                break
    print(f"写了 {n} 条三维记录到 {out_path};外部实现嵌不出来的 {failed} 条")
    return n


def canon(smiles, drop_stereo=False):
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


def blocks(path):
    """`(第几条, 这一段原文)`,按 `$$$$` 切。"""
    buf = []
    idx = 0
    for line in open(path, encoding="utf-8"):
        if line.rstrip("\n") == "$$$$":
            yield idx, "".join(buf)
            buf = []
            idx += 1
        else:
            buf.append(line)


def compare(sdf_path, ours_path):
    print(f"外部实现:RDKit {rdkit.__version__}")
    mine = {}
    for line in open(ours_path, encoding="utf-8"):
        i, smi, _data = line.rstrip("\n").split("\t", 2)
        mine[int(i)] = smi

    same = stereo_diff = with_chiral = with_ez = skipped = 0
    failures = []
    for idx, block in blocks(sdf_path):
        got = mine.get(idx)
        if got is None:
            failures.append(f"第 {idx} 条:我方一行输出都没有")
            continue
        # 真值:外部实现自己从这份坐标读回来的东西
        ref = Chem.MolFromMolBlock(block, removeHs=False)
        if ref is None:
            skipped += 1
            continue
        Chem.AssignStereochemistryFrom3D(ref)
        if any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in ref.GetAtoms()):
            with_chiral += 1
        if any(
            b.GetStereo() not in (Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY)
            for b in ref.GetBonds()
        ):
            with_ez += 1

        want = canon(Chem.MolToSmiles(ref))
        if got.startswith("<"):
            failures.append(f"第 {idx} 条:外部实现读得出,我方 {got}")
            continue
        if canon(got) is None:
            failures.append(f"第 {idx} 条:我方写出的 `{got}` 外部实现读不了")
            continue
        if canon(got) == want:
            same += 1
            continue
        if canon(got, drop_stereo=True) == canon(Chem.MolToSmiles(ref), drop_stereo=True):
            stereo_diff += 1
            failures.append(f"第 {idx} 条:骨架对、立体不同 —— 我方 {canon(got)},外部 {want}")
            continue
        failures.append(f"第 {idx} 条:我方读成 {canon(got)},外部实现读成 {want}")

    hard = [f for f in failures if "立体不同" not in f]
    print(f"逐条一致 {same};骨架对但立体不同 {stereo_diff}(上限 {MAX_STEREO_DIFF});"
          f"读成别的分子 {len(hard)};外部实现自己读不了 {skipped}")
    print(f"  参照侧带四面体的 {with_chiral} 条(下限 {MIN_WITH_CHIRAL});"
          f"带顺反的 {with_ez} 条(下限 {MIN_WITH_EZ})")
    for f in failures[:8]:
        print(f"  ✗ {f}")
    if hard:
        print("\n读出来不是同一个分子。")
        return 1
    if stereo_diff > MAX_STEREO_DIFF:
        print(f"\n只有立体不同的涨到 {stereo_diff} 条,超过上限 {MAX_STEREO_DIFF}")
        return 1
    if same + stereo_diff < MIN_RECORDS:
        print(f"\n只比过 {same + stereo_diff} 条,低于下限 {MIN_RECORDS} —— "
              "多半是写那份三维文件时半路断了,判据在一份截断的文件上跑")
        return 1
    if with_chiral < MIN_WITH_CHIRAL:
        print(f"\n参照侧只有 {with_chiral} 条带四面体,低于下限 {MIN_WITH_CHIRAL} —— "
              "手性那一档被喂空了")
        return 1
    if with_ez < MIN_WITH_EZ:
        print(f"\n参照侧只有 {with_ez} 条带顺反,低于下限 {MIN_WITH_EZ} —— 顺反那一档被喂空了")
        return 1
    print("\n从三维坐标读回来逐条一致。")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", nargs=2, metavar=("OUT", "CORPUS"))
    ap.add_argument("--compare", nargs=2, metavar=("SDF", "OURS"))
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()
    if args.write:
        write_sdf(args.write[0], args.write[1], args.limit)
        return 0
    if args.compare:
        return compare(args.compare[0], args.compare[1])
    ap.error("要么 --write 要么 --compare")
    return 2


if __name__ == "__main__":
    sys.exit(main())
