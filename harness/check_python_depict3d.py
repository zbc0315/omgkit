#!/usr/bin/env python3
"""**三维图的 Python 绑定与 Rust 侧逐字节相同。**

    cargo run -p omgkit-depict --release --example dump_depict3d -- \\
        harness/corpus/large.smi 400 > /tmp/three.jsonl
    .venv/bin/python harness/check_python_depict3d.py /tmp/three.jsonl

绑定这一层是**转发**,不是重新实现 —— 所以判据是最强的那种:同一个分子、
同一套样式,两边吐出来的 SVG 一个字节都不许差。

# 为什么值得单判一条

绑定里能出的错都是静默的:样式名映射错一档(拿球棍的参数画了棍状)、
把没补氢的那份分子传下去(三维图里氢是看得见的实体,少了就是另一个分子)、
默认样式改了名。这些都画得出图,而且图看着都对。

`depiction_3d_report` 里的原子落点同样逐项比 —— 那是公开 API,下游拿它在图上
加标注,错位了标注就贴到别的原子上。
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys

import omgkit

# 分母闸:见 harness/check_depict3d.py 的同名一段。
MIN_MOLECULES = 300
MIN_STYLES = 4


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("jsonl", type=pathlib.Path)
    args = ap.parse_args()

    bad: list[str] = []
    n_mol = 0
    n_svg = 0
    styles: set[str] = set()

    for raw in args.jsonl.read_text().splitlines():
        if not raw.strip():
            continue
        rec = json.loads(raw)
        smi = rec["smiles"]
        n_mol += 1
        conf = omgkit.parse_smiles(smi).conformer()

        # 坐标先对上。对不上的话下面比 SVG 就是在比两个不同的构象,
        # 报出来的分歧全是假的。
        got = conf.coords
        want = rec["coords"]
        if len(got) != len(want) or any(
            abs(a - b) > 1e-12 for g, w in zip(got, want) for a, b in zip(g, w)
        ):
            bad.append(f"{smi}: Python 与 Rust 生成的构象就不是同一个,SVG 没法比")
            continue

        for style_name, st in rec["styles"].items():
            styles.add(style_name)
            n_svg += 1
            if conf.to_svg(style_name) != st["svg"]:
                bad.append(f"{smi}/{style_name}: Python 与 Rust 的 SVG 不同")
            rep = conf.depiction_3d_report(style_name)
            if rep["style"] != style_name:
                bad.append(f"{smi}/{style_name}: report 报的样式是 {rep['style']!r}")
            if abs(rep["width"] - st["width"]) > 1e-12 or abs(rep["height"] - st["height"]) > 1e-12:
                bad.append(f"{smi}/{style_name}: 画布尺寸对不上")
            if rep["degenerate"] != st["degenerate"]:
                bad.append(f"{smi}/{style_name}: 视角退化标志对不上")
            for i, (a, b) in enumerate(zip(rep["atoms"], st["placed"])):
                if (abs(a["x"] - b[0]) > 1e-12 or abs(a["y"] - b[1]) > 1e-12
                        or abs(a["radius"] - b[2]) > 1e-12 or abs(a["depth"] - b[3]) > 1e-12):
                    bad.append(f"{smi}/{style_name}: 原子 {i} 的落点对不上")
                    break

    # 样式名拼错必须报错,不许静默退回默认 —— 否则拼错一个字母就得到一张
    # "看着对但不是你要的那个样式"的图。
    try:
        omgkit.parse_smiles("CCO").conformer().to_svg("ball-n-stick")
        bad.append("拼错的样式名没报错")
    except ValueError:
        pass

    print(f"分子 {n_mol}、SVG {n_svg} 张、样式 {len(styles)} 档")
    empty = []
    if n_mol < MIN_MOLECULES:
        empty.append(f"分子只有 {n_mol},下限 {MIN_MOLECULES}")
    if len(styles) < MIN_STYLES:
        empty.append(f"样式只有 {len(styles)} 档,下限 {MIN_STYLES}")
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
