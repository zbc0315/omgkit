#!/usr/bin/env python3
"""写二维 molblock 的 **Python 绑定**:与 Rust 侧逐字节相同。

# 这条判据守什么

`Mol.to_molblock_2d()` 要走一整条链:净化 → 感知顺反 → 排二维布局 → 指派楔形 →
**把画出来的那个分子**(可能补过一根显式 C–H)交给写出器。任何一步接错,交出去
的都是一张画错构型的图 —— 而线条本身看着一点毛病没有。

所以不另造真值,直接比:同一条 SMILES,`omgkit-depict/examples/dump_molblock`
走 Rust、`Mol.to_molblock_2d()` 走 Python,**写出来的 molblock 必须逐字节相同**。
布局全程无随机数,所以"逐字节"是可以要求的。

Rust 那一侧已经由 `check_wedge_readback.py` 与外部实现比过(它拿 RDKit 从二维 +
楔形指派手性,逐个中心比 CIP 码),这一侧与它逐字节相同,那条外部判据就继承了
过来。

# 楔形那一档要单独有个下限

"逐字节相同"在**两边都不画楔形**时同样成立 —— 那正是这条判据要抓的故障。所以
比过的块里带楔形码(键块第四列的 1 或 6)的条数配一条下限。

用法:

    cargo run -q -p omgkit-depict --release --example dump_molblock -- \\
        harness/corpus/large.smi > /tmp/blocks.txt
    python3 harness/check_python_molblock2d.py /tmp/blocks.txt
"""
import argparse
import sys

import omgkit

# 比过的块里**带楔形**的条数下限。
#
# Rust 那侧只导"画得出构型或如实报了画不出"的分子,实测 311 条,**每一条都带
# 楔形**。(495 是 `check_wedge_readback.py` 数的**中心**数,不是分子数 ——
# 两个数不是一回事,别混。)
MIN_WITH_WEDGE = 280


def blocks(path):
    """Rust 那侧的输出:`>>> 行号\\t原串` + `#unwedged …` + molblock + `$$$$`。"""
    smi = None
    lineno = None
    buf = []
    for line in open(path, encoding="utf-8"):
        if line.startswith(">>> "):
            lineno, smi = line[4:].rstrip("\n").split("\t", 1)
            buf = []
        elif line.startswith("#"):
            # Rust 侧的诊断行(`#unwedged` 等)不是 molblock 的内容。
            # 按前缀整体跳过,而不是只认 `#unwedged` —— 那边加一行诊断,
            # 这边就会把它当成 molblock 的一行,报出来的分歧指向错的东西。
            continue
        elif line.rstrip("\n") == "$$$$":
            yield lineno, smi, "".join(buf)
        else:
            buf.append(line)


def has_wedge(block):
    """键块第四列里有没有 1 或 6。"""
    lines = block.split("\n")
    if len(lines) < 4:
        return False
    try:
        na, nb = int(lines[3][0:3]), int(lines[3][3:6])
    except ValueError:
        return False
    return any(ln[9:12].strip() in ("1", "6") for ln in lines[4 + na : 4 + na + nb])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("blocks", help="dump_molblock 的输出")
    ap.add_argument("--min-checked", type=int, default=280)
    args = ap.parse_args()

    print(f"  omgkit wheel:{omgkit.__file__}")
    checked = with_wedge = 0
    failures = []
    for lineno, smi, want in blocks(args.blocks):
        try:
            got = omgkit.parse_smiles(smi).to_molblock_2d()
        except ValueError as e:
            failures.append(f"第 {lineno} 行 {smi}:Python 写不出来({e})")
            continue
        if got != want:
            # 只报第一处不同,整块贴出来没法看
            gl, wl = got.split("\n"), want.split("\n")
            at = next((i for i in range(min(len(gl), len(wl))) if gl[i] != wl[i]), None)
            where = f"第 {at} 行:Python `{gl[at]}` ≠ Rust `{wl[at]}`" if at is not None \
                else f"行数不同:Python {len(gl)},Rust {len(wl)}"
            failures.append(f"第 {lineno} 行 {smi}:{where}")
            continue
        checked += 1
        if has_wedge(want):
            with_wedge += 1

    print(f"逐字节相同 {checked} 条;不一致 {len(failures)} 条")
    print(f"  其中带楔形的 {with_wedge} 条(下限 {MIN_WITH_WEDGE})")
    for f in failures[:8]:
        print(f"  ✗ {f}")
    if failures:
        print("\nPython 绑定画出来的不是 Rust 画出来的那张图。")
        return 1
    if with_wedge < MIN_WITH_WEDGE:
        print(f"\n带楔形的只有 {with_wedge} 条,低于下限 {MIN_WITH_WEDGE} —— "
              "楔形那一档被喂空了,两边一起不画也是这个样子")
        return 1
    if checked < args.min_checked:
        print(f"\n只比过 {checked} 条,低于下限 {args.min_checked} —— 判据被喂空了")
        return 1
    print("\n与 Rust 侧逐字节相同。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
