#!/usr/bin/env python3
"""用外部实现裁判:画出来的楔形,别人读出来是不是同一个构型。

    cargo run -p omgkit-depict --release --example dump_molblock -- \\
        harness/corpus/large.smi > /tmp/blocks.txt
    python3 harness/check_wedge_readback.py /tmp/blocks.txt

# 为什么判官必须是外部实现

`stereo::assign_wedges` 是"试 Up/Down,取**反读回来**对的那一个"构造出来的,而
反读用的就是 `read_chirality` —— 两者共谋,拿它们的往返去检验是**空过的**,
函数自己的文档就写着这一点。

真正要问的是另一个问题:**别人照着这张图读,读出来是不是同一个分子。** 楔形
读法在一处是**有约定分歧**的:三个画出来的邻居若全挤在中心的同一侧(最大空隙
> 180°,中心落在它们围出的三角形**外面**),"隐式氢在楔形反面"这条读法与四面体
的读法不再等价,不同实现读出**对映体**。这类错误拓扑完全正确、线条毫无毛病,
只有分子是镜像的 —— 自己验自己永远发现不了。

# 按**中心**比,不按分子比

按分子比太粗:一个分子里可能既有如实报了 `unwedged` 的中心、又有没报却画错的
中心,按分子会把后者藏在前者后面。所以逐个原子比 CIP 码。

# 四种结局要分开

- **一致**:图画对了。
- **不一致,但已报 `unwedged`**:图上本来就没画出那个中心,如实说过了。丢信息,
  不好,但下游拿到的是"未指定",不会当真。
- **外部实现读不出**:图上画了,而判官给不出 CIP 码。**这不等于画错了** ——
  RDKit 的 2D→手性只认 S(16) 与 Se(34) 这两种三配位中心
  (`Chirality.cpp`:`anum != 16 && anum != 34 && tnzDegree != 4` 就跳过),
  三配位的**磷**不在它的支持范围。这一档单独计数,不当违例,也不藏起来。
- **读成了相反的构型**:图说自己画对了,别人读出来是**对映体**。
  **这一档必须是 0。**

**把"读不出"和"读反了"混成一档是不行的**:前者是读者的覆盖面,后者是我们画错了,
危害差着量级。

双键顺反两边都先清掉:RDKit 会从 2D 坐标读出顺反,而输入 SMILES 里未指定顺反的
双键会因此"凭空多出"立体信息 —— 那是系统性差异,不是画错。
"""

import argparse
import collections
import pathlib
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


def cip_codes(mol):
    """`{原子下标: 'R'/'S'}`,没定的不进表。读不了返回 None。"""
    if mol is None:
        return None
    m = Chem.Mol(mol)
    # 只留四面体 —— 双键顺反见模块文档
    for b in m.GetBonds():
        b.SetStereo(Chem.BondStereo.STEREONONE)
        b.SetBondDir(Chem.BondDir.NONE)
    try:
        Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    except Exception:
        return None
    return {
        a.GetIdx(): a.GetPropsAsDict()["_CIPCode"]
        for a in m.GetAtoms()
        if a.HasProp("_CIPCode")
    }


def blocks(text):
    """切成 (行号, SMILES, unwedged 集合, molblock)。"""
    for chunk in text.split("$$$$\n"):
        if not chunk.strip():
            continue
        lines = chunk.splitlines()
        head = next((i for i, l in enumerate(lines) if l.startswith(">>> ")), None)
        if head is None:
            continue
        lineno, smi = lines[head][4:].split("\t", 1)
        unw = set()
        body = head + 1
        if lines[body].startswith("#unwedged"):
            rest = lines[body][len("#unwedged") :].strip()
            unw = {int(x) for x in rest.split(",") if x}
            body += 1
        yield lineno, smi, unw, "\n".join(lines[body:]) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("dump", type=pathlib.Path, help="dump_molblock 的输出")
    ap.add_argument("--show", type=int, default=20, help="最多列几例")
    args = ap.parse_args()

    tally = collections.Counter()
    bad = []
    for lineno, smi, unw, block in blocks(args.dump.read_text()):
        ref = Chem.MolFromSmiles(smi)
        got = Chem.MolFromMolBlock(block)
        if ref is None:
            tally["输入外部实现读不了(跳过)"] += 1
            continue
        if got is None:
            tally["**导出的 molblock 读不了**"] += 1
            bad.append((lineno, smi, "molblock 解析失败"))
            continue
        a, b = cip_codes(ref), cip_codes(got)
        if a is None or b is None:
            tally["指派手性失败(跳过)"] += 1
            continue
        # 原子下标两边一致:molblock 按本实现的原子序写,外部实现按序读
        for k in set(a) | set(b):
            if a.get(k) == b.get(k):
                tally["构型一致"] += 1
            elif k in unw:
                tally["不一致,但已报 unwedged"] += 1
            elif b.get(k) is None:
                # 画了,但判官给不出 CIP 码 —— 是它的覆盖面,不是我们画错了
                tally["外部实现读不出(见模块文档)"] += 1
            else:
                tally["**读成了相反的构型**"] += 1
                bad.append(
                    (lineno, smi, f"原子 {k}:本实现画成 {b.get(k)},应当是 {a.get(k)}")
                )

    print(f"{'档':<30}{'中心数':>8}")
    for k, v in tally.most_common():
        print(f"{k:<30}{v:>8}")
    if bad:
        print(f"\n=== 前 {args.show} 例 ===")
        for lineno, smi, why in bad[: args.show]:
            print(f"  行 {lineno}: {smi}\n      {why}")
        if len(bad) > args.show:
            print(f"  …… 另有 {len(bad) - args.show} 例")
    hard = tally["**读成了相反的构型**"] + tally["**导出的 molblock 读不了**"]
    print(f"\n{'全部通过' if hard == 0 else f'共 {hard} 处'}")
    return 1 if hard else 0


if __name__ == "__main__":
    sys.exit(main())
