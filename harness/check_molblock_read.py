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

# 立体化学眼下**不比**

我方的读取器还不给立体赋值(二维靠楔形、三维靠坐标,两者都要用对称等价类,
那在 L1 之上)。所以这里两侧都把立体抹掉再比,并且**把带立体的条数打出来** ——
免得"零分歧"读起来像是立体也守住了。那一档补上之后这个开关就该去掉。

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


def flat(smiles):
    """抹掉立体、再去掉显式氢之后的规范串。读不了时给 None。

    **去氢那一步不能省。** 外部实现的 SMILES 解析器会保留"承载方向键的显式氢"
    (`[H]/[O+]=c1…` 里那个 H 是一个原子),而我方眼下不读立体、写出的串没有
    方向键,同一个氢就被合并进邻居 —— 两边原子数一个 9 一个 8,比出来是"不同的
    分子",而那是**比法**造成的,不是读法。实测 8831 条里 26 条栽在这上面。
    """
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
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
        want = flat(Chem.MolToSmiles(ref))
        mine = flat(got)
        if mine is None:
            unreadable += 1
            failures.append(f"第 {lineno} 行 {smi}:我方写出的 `{got}` 外部实现读不了")
        elif mine == want:
            same += 1
        else:
            diff += 1
            failures.append(f"第 {lineno} 行 {smi}:我方读成 {mine},外部实现读成 {want}")

    print(f"读回来一致 {same};不一致 {diff};读不了/写不出 {unreadable};"
          f"外部实现自己读不了 {skipped}")
    print(f"  外部实现写成了 V3000、我方明确拒收的 {v3000} 条(上限 {MAX_V3000})")
    print(f"  其中带立体标记的 {with_stereo} 条 —— **立体这一档眼下没比**"
          f"(我方读取器还不赋值,见模块文档)")
    if failures:
        for f in failures[:8]:
            print(f"  ✗ {f}")
    if diff or unreadable:
        print("\n读出来不是同一个分子。")
        return 1
    if v3000 > MAX_V3000:
        print(f"\nV3000 那一档涨到 {v3000} 条,超过上限 {MAX_V3000} —— "
              "要么语料变了,要么该把 V3000 也读起来了")
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
