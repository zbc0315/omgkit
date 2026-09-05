#!/usr/bin/env python3
"""**二维结构式的 Python 绑定与 Rust 侧逐字节相同。**

    cargo run -p omgkit-depict --release --example dump_depict2d -- \\
        harness/corpus/large.smi > /tmp/two.jsonl
    .venv/bin/python harness/check_python_depict2d.py /tmp/two.jsonl

绑定这一层是**转发**,不是重新实现 —— 所以判据是最强的那种:同一个分子、
同一套规范、同一组配色,两边吐出来的 SVG 一个字节都不许差。

# 为什么值得单判一条

绑定里能出的错都是静默的:忘了跑顺反感知(顺式和反式画成同一张图)、
`aromatic_fill=False` 却把底色接上了、两个颜色参数接反了、规范名映射错一套
(拿 ChemDraw 的键长画了 ACS 的图)。这些都画得出图,而且图看着都对。

# 拼错的名字与写坏的颜色都要报错

静默退回默认值的话,拼错一个字母就得到一张"看着对但不是你要的"的图。
颜色这一条**关着底色时也要报错** —— 只在用得上时才校验,那么关着底色调好
的一组颜色,开的那天才发现拼错了。
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys

import omgkit

# 分母闸:见 harness/denominator.py 的开头一段。
MIN_MOLECULES = 300
MIN_SVGS = 1800
MIN_STYLES = 2


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("jsonl", type=pathlib.Path)
    args = ap.parse_args()

    bad: list[str] = []
    n_mol = 0
    n_svg = 0
    styles: set[str] = set()

    with args.jsonl.open() as fh:
        for raw in fh:
            if not raw.strip():
                continue
            rec = json.loads(raw)
            smi = rec["smiles"]
            n_mol += 1
            mol = omgkit.parse_smiles(smi)
            centre, edge = rec["custom"]["centre"], rec["custom"]["edge"]
            for style, v in rec["styles"].items():
                styles.add(style)
                got = {
                    "plain": mol.to_svg(style),
                    "fill": mol.to_svg(style, aromatic_fill=True),
                    "custom": mol.to_svg(
                        style, aromatic_fill=True, fill_centre=centre, fill_edge=edge
                    ),
                }
                for tag, s in got.items():
                    n_svg += 1
                    if s != v[tag]:
                        bad.append(f"{smi}/{style}/{tag}: Python 与 Rust 的 SVG 不同")
                # 关着底色时给了颜色也不许生效 —— 颜色参数只在开着时说话
                if mol.to_svg(style, fill_centre=centre, fill_edge=edge) != v["plain"]:
                    bad.append(f"{smi}/{style}: 关着底色却让颜色参数生效了")
                n_svg += 1

    # 拼错的规范名、写坏的颜色都要抛 ValueError
    m = omgkit.parse_smiles("c1ccccc1")
    for what, call in [
        ("规范名", lambda: m.to_svg("ACS 1996")),
        ("中心色", lambda: m.to_svg(fill_centre="白")),
        ("外缘色", lambda: m.to_svg(fill_edge="#abc")),
        # 关着底色时颜色照样要校验
        ("关着底色时的颜色", lambda: m.to_svg(aromatic_fill=False, fill_edge="nope")),
    ]:
        try:
            call()
            bad.append(f"{what}写错了却没报错")
        except ValueError:
            pass

    print(f"分子 {n_mol}、SVG {n_svg} 张、规范 {len(styles)} 套")
    empty = []
    if n_mol < MIN_MOLECULES:
        empty.append(f"分子只有 {n_mol},下限 {MIN_MOLECULES}")
    if n_svg < MIN_SVGS:
        empty.append(f"SVG 只有 {n_svg} 张,下限 {MIN_SVGS}")
    if len(styles) < MIN_STYLES:
        empty.append(f"规范只有 {len(styles)} 套,下限 {MIN_STYLES}")
    if empty:
        print("判据没东西可判:", file=sys.stderr)
        for e in empty:
            print(f"  {e}", file=sys.stderr)
        return 1
    if bad:
        print(f"分歧 {len(bad)} 条(最多列 20):", file=sys.stderr)
        for b in bad[:20]:
            print(f"  {b}", file=sys.stderr)
        return 1
    print("全部通过")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
