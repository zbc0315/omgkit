#!/usr/bin/env python3
"""**芳香环底色:多边形与渐变从 SVG 反读,外部真值是 RDKit 数出来的芳香环。**

    cargo run -p omgkit-depict --release --example dump_depict2d -- \\
        harness/corpus/large.smi 400 > /tmp/two.jsonl
    .venv/bin/python harness/check_depict2d.py /tmp/two.jsonl

判官不看 `Scene`,只看**产品真正吐出来的那段 SVG** —— 从里面把 `<polygon>`
和 `<radialGradient>` 读回来,再独立算一遍该是什么样。

# 五件事,各由不同的东西当真值

| 判什么 | 真值从哪来 |
|---|---|
| 铺了几块底色、每块几个角 | **RDKit** 的环信息与芳香标记,一行本库的代码都不经过 |
| 底色在不在最底层 | SVG 里的元素次序 —— 最后一块多边形要排在第一条线/字之前 |
| 高光落在哪 | 从多边形顶点自己算质心、自己挑最靠左上的角,取中点 |
| 开底色有没有动到别的东西 | 把多边形与 `<defs>` 抠掉之后,与不开底色那一份**逐字节**比 |
| 两个颜色接没接上 | 自定义配色那一份里,焦点与外缘各是各的颜色 |

最后那一条是这份判官里最有分量的:布局、线条、文字、画布尺寸一律不许因为
铺底色而变。一句"底色不影响别的"是声称,抠掉再逐字节比才是判据。

# 高光的期望值**不引用实现里的常数**

实现里光源方向写成 `LIGHT`(左上的单位向量)。判据这边直接写"x 与 y 都往
小走",不去 import 那个常数 —— 引用它的话,把光源翻到右下两边一起动,
这条判据永远打不红。

# 只判 ACS 那一套?不,两套都判

两套规范的键长差 2.08 倍,底色的半径、焦点全都跟着变。只判一套的话,
"半径写死成某个数"这类错在另一套上才现形。
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import re
import sys

from rdkit import Chem
from rdkit import RDLogger

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from denominator import corpus_size, verdict  # noqa: E402

RDLogger.DisableLog("rdApp.*")

# 分母闸:见 harness/denominator.py。
MIN_MOLECULES = 300
MIN_AROMATIC_MOLECULES = 150
MIN_FILLS = 300
MIN_RING_SIZES = 2
MIN_STYLES = 2

# 从 SVG 里读回来的坐标只有两位小数,而质心是六个顶点平均出来的 ——
# 误差按 0.005/√6 量级走。给到 0.02 pt,比舍入宽一档、比任何真错窄得多。
TOL_PT = 0.02

# **两个角一样靠左上时,判据分不出实现挑了哪一个。**
#
# 实现按全精度的坐标挑,判据只拿得到 SVG 里那两位小数。正六边形转到某个角度
# 上会出现两个顶点的"靠左上"程度差在舍入以内 —— 实测语料里
# `N1(C2(OCC)CCCCC(C1=O)2)C3=CC=CC=C3` 的苯环就是,两个顶点差 1e-9 量级。
#
# 那一档**这条判据判不了**,只能接受"是打平的那几个角之一"。别把它读成
# "高光可以随便落" —— 只有真打平时才放宽,而放宽了多少下面会打印出来。
TIE_PT = 0.05

POLY = re.compile(r'<polygon points="([^"]*)" fill="url\(#([^)]*)\)"/>')
GRAD = re.compile(
    r'<radialGradient id="([^"]*)" gradientUnits="userSpaceOnUse" '
    r'cx="([-\d.]+)" cy="([-\d.]+)" r="([-\d.]+)" fx="([-\d.]+)" fy="([-\d.]+)">'
    r'<stop offset="0" stop-color="([^"]*)"/>'
    r'<stop offset="1" stop-color="([^"]*)"/>'
    r"</radialGradient>"
)


def aromatic_rings(smi: str) -> list[int] | None:
    """RDKit 数出来的芳香环,返回每个环的大小。读不了就 `None`。

    # 判的是"环上每根**键**都芳香",不是"每个**原子**都芳香"

    先前写的是按原子。那是错的,而且是**实测出来的错**:

        [H]/N=c/1\\[nH]nc-2c(n1)-c3cccc4c3c2ccc4

    这个分子的 SSSR 里有一个五元环,五个原子全芳香,而环上有一根键是单键
    (SMILES 里那个 `-`)—— 那个环本身不是芳香环,它只是被三个芳香环围出来的
    一圈。按原子判会要求给它铺底色,全语料上这样的分子有一批。

    经典的例子是联苯烯(biphenylene):中间那个四元环四个角都是芳香碳,
    环本身当然不芳香。

    RDKit 的芳香性感知与 SSSR 都是它自己那一套,与本库没有一行共用代码 ——
    谓词读起来与实现那边同一句话,底下算的却是两回事。
    """
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    ri = m.GetRingInfo()
    return sorted(
        len(r)
        for r, br in zip(ri.AtomRings(), ri.BondRings())
        if all(m.GetBondWithIdx(b).GetIsAromatic() for b in br)
    )


def polygons(svg: str) -> list[tuple[list[tuple[float, float]], str]]:
    out = []
    for pts, gid in POLY.findall(svg):
        vs = []
        for tok in pts.split():
            x, y = tok.split(",")
            vs.append((float(x), float(y)))
        out.append((vs, gid))
    return out


def gradients(svg: str) -> dict[str, tuple[float, float, float, float, float, str, str]]:
    out = {}
    for gid, cx, cy, r, fx, fy, c0, c1 in GRAD.findall(svg):
        out[gid] = (float(cx), float(cy), float(r), float(fx), float(fy), c0, c1)
    return out


def stripped(svg: str) -> str:
    """抠掉多边形与 `<defs>` 块之后剩下的那份 SVG。

    **按行删,不是把匹配替换成空串** —— 替换成空串会在原地留一个空行,于是
    每个带芳香环的分子都"不同",而那是判据自己造出来的分歧(实测:先前那版
    在全语料上报了 13202 条,一条真的都没有)。
    """
    out = re.sub(r"<defs>\n.*?</defs>\n", "", svg, flags=re.S)
    return "".join(l + "\n" for l in out.splitlines() if not POLY.fullmatch(l.strip()))


def check_one(
    smi: str,
    style: str,
    rec: dict,
    want_rings: list[int],
    bad: list[str],
    ties: list[int],
) -> int:
    """判一个分子的一套规范,返回这一份里数出来的底色块数。"""
    plain, fill, custom = rec["plain"], rec["fill"], rec["custom"]
    tag = f"{smi}/{style}"

    # ① 不开底色的那一份,一个多边形、一个渐变都不许有
    if "<polygon" in plain or "<radialGradient" in plain:
        bad.append(f"{tag}: 没开底色却写出了底色")

    polys = polygons(fill)
    grads = gradients(fill)

    # ② 块数与每块的角数 —— 真值是 RDKit
    got_rings = sorted(len(vs) for vs, _ in polys)
    if got_rings != want_rings:
        bad.append(f"{tag}: 底色是 {got_rings},RDKit 数出来的芳香环是 {want_rings}")

    # ③ 图层:最后一块底色要排在第一条线/楔形/文字之前
    if polys:
        last = fill.rfind("<polygon")
        for t in ("<line", "<path", "<text"):
            first = fill.find(t)
            if first >= 0 and first < last:
                bad.append(f"{tag}: 有底色画在了 {t} 后面,会盖住它")
                break

    # ④ 逐块的几何
    for vs, gid in polys:
        if gid not in grads:
            bad.append(f"{tag}: 多边形引用了不存在的渐变 {gid}")
            continue
        cx, cy, r, fx, fy, c0, c1 = grads[gid]
        gx = sum(v[0] for v in vs) / len(vs)
        gy = sum(v[1] for v in vs) / len(vs)
        if abs(cx - gx) > TOL_PT or abs(cy - gy) > TOL_PT:
            bad.append(f"{tag}: 渐变圆心 ({cx},{cy}) 不是顶点质心 ({gx:.2f},{gy:.2f})")
        far = max(math.hypot(v[0] - gx, v[1] - gy) for v in vs)
        if abs(r - far) > TOL_PT:
            bad.append(f"{tag}: 渐变半径 {r} 不是外接圆半径 {far:.2f}")
        # 最靠左上 = x 与 y 都往小走(画布 y 向下)。**不引用实现里的 LIGHT**。
        lit = [-(v[0] - gx) - (v[1] - gy) for v in vs]
        top = max(lit)
        tied = [v for v, l in zip(vs, lit) if top - l <= TIE_PT]
        if len(tied) > 1:
            ties[0] += 1
        ok = any(
            abs(fx - (v[0] + gx) / 2) <= TOL_PT and abs(fy - (v[1] + gy) / 2) <= TOL_PT
            for v in tied
        )
        if not ok:
            best = tied[0]
            wx, wy = (best[0] + gx) / 2, (best[1] + gy) / 2
            bad.append(
                f"{tag}: 高光在 ({fx},{fy}),该在最靠左上的角与环心的中点 ({wx:.2f},{wy:.2f})"
            )
        if (c0, c1) != ("#ffffff", "#add8e6"):
            bad.append(f"{tag}: 默认配色写成了 {c0} → {c1}")

    # ⑤ 铺底色只多了底色:抠掉之后与不开的那份逐字节相同
    if stripped(fill) != plain:
        bad.append(f"{tag}: 抠掉底色之后与不开底色的那份不同 —— 铺底色动到了别的东西")
    if stripped(custom) != plain:
        bad.append(f"{tag}: 自定义配色那份抠掉底色之后与不开底色的那份不同")

    # ⑥ 颜色接没接上,而且没写反
    want_c0, want_c1 = rec["want_custom"]
    for gid, g in gradients(custom).items():
        if (g[5], g[6]) != (want_c0, want_c1):
            bad.append(f"{tag}: 自定义配色是 {g[5]} → {g[6]},要的是 {want_c0} → {want_c1}")
            break

    return len(polys)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("jsonl", type=pathlib.Path)
    ap.add_argument("--corpus", type=pathlib.Path, default=None)
    ap.add_argument("--cap", type=int, default=15)
    args = ap.parse_args()

    bad: list[str] = []
    n_mol = 0
    n_aromatic = 0
    n_fills = 0
    ring_sizes: set[int] = set()
    styles: set[str] = set()
    # 有几块底色的"最靠左上"在舍入以内打平了 —— 那几块高光只判到"是打平的
    # 那几个角之一"。打印出来,免得放宽了多少没人知道。
    ties = [0]

    for raw in args.jsonl.read_text().splitlines():
        if not raw.strip():
            continue
        rec = json.loads(raw)
        smi = rec["smiles"]
        want = aromatic_rings(smi)
        if want is None:
            # RDKit 读不了这条 SMILES,那就没有外部真值可比。跳过要计进
            # 分母闸,不能静默。
            continue
        n_mol += 1
        if want:
            n_aromatic += 1
        ring_sizes.update(want)
        custom = (rec["custom"]["centre"], rec["custom"]["edge"])
        for style, v in rec["styles"].items():
            styles.add(style)
            n_fills += check_one(smi, style, {**v, "want_custom": custom}, want, bad, ties)

    print(f"分子 {n_mol}、其中带芳香环的 {n_aromatic}、规范 {len(styles)} 套")
    print(f"  底色 {n_fills} 块,环大小 {sorted(ring_sizes)}")
    print(f"  最靠左上的角在舍入以内打平、只判到\"是其中之一\"的:{ties[0]} 块")

    empty = []
    if n_mol < MIN_MOLECULES:
        empty.append(f"分子只有 {n_mol},下限 {MIN_MOLECULES}")
    if n_aromatic < MIN_AROMATIC_MOLECULES:
        empty.append(f"带芳香环的分子只有 {n_aromatic},下限 {MIN_AROMATIC_MOLECULES}")
    if n_fills < MIN_FILLS:
        empty.append(f"底色只有 {n_fills} 块,下限 {MIN_FILLS}")
    if len(ring_sizes) < MIN_RING_SIZES:
        empty.append(f"芳香环只有 {len(ring_sizes)} 种大小,下限 {MIN_RING_SIZES}")
    if len(styles) < MIN_STYLES:
        empty.append(f"规范只有 {len(styles)} 套,下限 {MIN_STYLES}")
    if args.corpus is not None:
        why = verdict(corpus_size(args.corpus), n_mol, n_mol, args.cap)
        if why:
            empty.append(why)
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
