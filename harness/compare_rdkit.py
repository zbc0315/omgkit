#!/usr/bin/env python3
"""把 omgkit 与 RDKit 画的同一个分子并排放在一张图上,供人眼比对。

    cargo run -p omgkit-depict --features raster --example draw -- <目录>
    python3 harness/compare_rdkit.py <目录>

判据守得住"环内双键有没有画到环外""立体中心有没有画出来"这类**可判定**的性质;
守不住"这张图好不好看"。后者只能看,而单看一张图看不出好坏 —— 要有个参照。

# 比的是什么,不比什么

两边**用同一个键长**(omgkit 出图时 300 dpi、ACS 键长 14.4 pt,合 60 px),所以
并排看到的粗细、字号、疏密都是可比的。比不了的是布局本身的取向:两边各自有
规范朝向的做法,同一个分子摆的角度可能差 90°,那不是谁对谁错。

RDKit 那一侧的规范:ACS 一栏用 RDKit 自带的 ACS 1996 模式,ChemDraw 一栏用
RDKit 默认值 —— RDKit 没有 ChemDraw 规范,拿默认值当参照。
"""

import io
import sys
from pathlib import Path

from PIL import Image, ImageChops, ImageDraw, ImageFont
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

# omgkit 出图用 300 dpi(scale = 300/72)
SCALE = 300.0 / 72.0


def styles(d):
    """规范表从 `styles.tsv` 读,**不在这里硬编码**。

    键长写死在两边的话,改了 `Style` 的键长,两个面板就不再同尺度,而并排图
    看着依旧正常 —— 那种对照会把人引到错的结论上。ACS 模式那一列按规范名认。
    """
    f = d / "styles.tsv"
    if not f.exists():
        sys.exit(f"找不到 {f} —— 先跑 example draw 生成它")
    out = {}
    for line in f.read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) != 4:
            sys.exit(f"{f} 这行不是四列:{line!r}")
        tag, name, bond_pt, label_pt = parts
        out[tag] = (name, float(bond_pt) * SCALE, float(label_pt) * SCALE, name.startswith("ACS"))
    if not out:
        sys.exit(f"{f} 是空的")
    return out
PAD = 24  # 面板四周留白(px)
BAR = 64  # 顶部标题条高度(px)


# 画布给得足够大,再把白边裁掉 —— RDKit 的自适应画布(-1, -1)不认
# `fixedBondLength`,给出来的图键长对不上,并排就没法比了
CANVAS = (3200, 2400)


def font(size):
    """找一个排得了中文的真字体。

    PIL 自带的点阵字体在 300 dpi 的图上小得看不清;而 Arial 一类西文字体
    **排中文会全变成方框** —— 那种图看着像渲染坏了。
    """
    for p in (
        "/System/Library/Fonts/STHeiti Light.ttc",
        "/System/Library/Fonts/Hiragino Sans GB.ttc",
        "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    ):
        if Path(p).exists():
            try:
                f = ImageFont.truetype(p, size)
            except OSError:
                continue
            if f.getbbox("规范")[2] > 0:  # 真的排得出中文才算数
                return f
    return ImageFont.load_default()


def rdkit_png(smi, bond_px, label_px, acs):
    """RDKit 画一张,键长压到 `bond_px` —— 与 omgkit 那一侧同尺度才好比。"""
    m = Chem.MolFromSmiles(smi)
    if m is None:
        raise ValueError(f"RDKit 解析不了:{smi}")
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(m)
    d = rdMolDraw2D.MolDraw2DCairo(*CANVAS)
    o = d.drawOptions()
    if acs:
        # ACS 1996 模式按分子的平均键长定线宽;键长本身随后再压到我们要的值
        rdMolDraw2D.SetACS1996Mode(o, rdMolDraw2D.MeanBondLength(m))
    o.fixedBondLength = bond_px
    # 字号**没有**去对齐:试过设 `baseFontSize`,RDKit 的 ACS 模式不理会它。
    # 所以两边的标签大小各按各自的实现 —— 左边是 ACS 规范的 10 pt(合键长的
    # 0.69),右边是 RDKit 自己选的,明显小一号。比图的时候要知道这一点。
    _ = label_px
    rdMolDraw2D.PrepareAndDrawMolecule(d, m)
    d.FinishDrawing()
    return d.GetDrawingText()


def trim(im):
    """裁掉四周的白边。两边的画布边距各有各的算法,不裁就没法按内容对齐。"""
    im = im.convert("RGB")
    bg = Image.new("RGB", im.size, (255, 255, 255))
    box = ImageChops.difference(im, bg).getbbox()
    return im.crop(box) if box else im


def panel(im, w, h, title, sub, f_title, f_sub):
    """把一张图放进固定大小的面板,顶上写标题。"""
    out = Image.new("RGB", (w, h), (255, 255, 255))
    dr = ImageDraw.Draw(out)
    dr.text((PAD, 14), title, fill=(20, 20, 20), font=f_title)
    dr.text((PAD, 40), sub, fill=(120, 120, 120), font=f_sub)
    dr.line([(0, BAR - 1), (w, BAR - 1)], fill=(220, 220, 220), width=1)
    # 内容在面板里居中
    x = (w - im.width) // 2
    y = BAR + (h - BAR - im.height) // 2
    out.paste(im, (max(x, 0), max(y, BAR)))
    return out


def main():
    d = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
    table = d / "mols.tsv"
    if not table.exists():
        sys.exit(f"找不到 {table} —— 先跑 example draw 生成它")
    mols = []
    for line in table.read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) != 2:
            sys.exit(f"{table} 这行不是两列:{line!r}")
        mols.append(parts)
    if not mols:
        sys.exit(f"{table} 是空的")

    f_title, f_sub = font(26), font(19)
    style_table = styles(d)
    made = 0
    for name, smi in mols:
        for tag, (style_name, bond_px, label_px, acs) in style_table.items():
            mine = d / f"{name}.{tag}.png"
            if not mine.exists():
                sys.exit(f"缺 {mine} —— 先跑 example draw")
            raw_a = Image.open(mine)
            a = trim(raw_a)
            # 全白的图 `getbbox()` 返回 None,`trim` 会把整张原样还回来 ——
            # 于是空图会安安静静地拼进对照表。**这一侧才是要防的那一侧**:
            # 脚本存在的意义就是发现 omgkit 画错或没画
            if a.size == raw_a.size:
                sys.exit(f"{mine}:一个非白像素都没有,或者白边为零 —— 图不对")
            raw = Image.open(io.BytesIO(rdkit_png(smi, bond_px, label_px, acs)))
            b = trim(raw)
            # 顶到画布边就说明被裁掉了一截 —— 那种图看着基本正常,只是少了一块,
            # 极容易当成"omgkit 画得更全"
            if b.width >= raw.width - 2 or b.height >= raw.height - 2:
                sys.exit(f"{name}.{tag}:RDKit 那张顶到画布边了,把 CANVAS 调大")

            w = max(a.width, b.width) + 2 * PAD
            h = max(a.height, b.height) + BAR + 2 * PAD
            left = panel(a, w, h, "omgkit", style_name, f_title, f_sub)
            right = panel(
                b,
                w,
                h,
                "RDKit",
                "ACS 1996 模式" if acs else "默认",
                f_title,
                f_sub,
            )

            out = Image.new("RGB", (w * 2 + 1, h), (255, 255, 255))
            out.paste(left, (0, 0))
            out.paste(right, (w + 1, 0))
            ImageDraw.Draw(out).line([(w, 0), (w, h)], fill=(200, 200, 200), width=1)
            out.save(d / f"{name}.{tag}.compare.png")
            made += 1
        print(f"{name:16} 两套规范各一张")
    print(f"\n共 {made} 张对照图,写在 {d}")


if __name__ == "__main__":
    main()
