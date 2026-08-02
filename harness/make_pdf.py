#!/usr/bin/env python3
"""把一个目录里的图拼成一份多页 PDF,一页一张,方便翻着看。

    python3 harness/make_pdf.py <图目录> <输出.pdf> [--sort desc] [--diag <draw 的输出>]

# 为什么不是把目录里所有文件都塞进去

`draw` 对每个分子每套规范出三种格式(svg/png/jpg),画的是**同一张图**。三种
都放进 PDF 只是把同一页印三遍。所以:目录里有 `*.compare.png`(omgkit 与 RDKit
并排)就只用它 —— 它左半边就是 omgkit 自己那张;没有就用 `*.png`。
"""

import subprocess
import sys
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

DPI = 150
PAGE = (int(11.69 * DPI), int(8.27 * DPI))  # A4 横放
MARGIN = int(0.4 * DPI)
HEAD = int(0.34 * DPI)


def font(size):
    for p in (
        "/System/Library/Fonts/STHeiti Light.ttc",
        "/System/Library/Fonts/Hiragino Sans GB.ttc",
        "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
    ):
        if Path(p).exists():
            try:
                f = ImageFont.truetype(p, size)
            except OSError:
                continue
            if f.getbbox("规范")[2] > 0:
                return f
    return ImageFont.load_default()


def diagnostics(path):
    """跑一遍 `draw`,把每张图的诊断数字取回来,印在页眉上。

    翻图的时候最想知道的就是"这张有没有画不好的地方",而那正是判据守不住、
    只能靠人看的部分 —— 数字放在图边上,才对得起来。
    """
    out = {}
    try:
        r = subprocess.run(
            ["cargo", "run", "-q", "--release", "--example", "draw",
             "--features", "raster", "--", str(path)],
            capture_output=True, text=True, check=True, cwd=Path(__file__).parent.parent,
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        return out
    for line in r.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 4:
            # parts = [名字, 规范, 尺寸, "pt", 退化…, 未解冲突…, …]
            out[f"{parts[0]}.{parts[1]}"] = f"{parts[2]} pt   " + " ".join(parts[4:])
    return out


def page(im, caption, f):
    sheet = Image.new("RGB", PAGE, (255, 255, 255))
    dr = ImageDraw.Draw(sheet)
    if caption:
        dr.text((MARGIN, MARGIN // 2), caption, fill=(40, 40, 40), font=f)
    box = (PAGE[0] - 2 * MARGIN, PAGE[1] - 2 * MARGIN - HEAD)
    # **也放大。** 只缩不放的话,小分子在整页里只占一角,翻起来全是白纸。
    # 上限 3 倍:源图是 300 dpi 出的,放到 3 倍还有 100 dpi,线条不至于糊。
    k = min(box[0] / im.width, box[1] / im.height, 3.0)
    im = im.resize((max(int(im.width * k), 1), max(int(im.height * k), 1)), Image.LANCZOS)
    sheet.paste(im, ((PAGE[0] - im.width) // 2, HEAD + (box[1] - im.height) // 2 + MARGIN))
    return sheet


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    d, out = Path(sys.argv[1]), Path(sys.argv[2])
    desc = "--sort" in sys.argv and "desc" in sys.argv
    want_diag = "--diag" in sys.argv

    files = sorted(d.glob("*.compare.png")) or sorted(
        p for p in d.glob("*.png") if not p.name.endswith(".compare.png")
    )
    if not files:
        sys.exit(f"{d} 里一张 png 都没有")
    files.sort(key=lambda p: p.name, reverse=desc)

    diag = diagnostics(d) if want_diag else {}
    f = font(26)
    pages = []
    for p in files:
        stem = p.name.replace(".compare.png", "").replace(".png", "")
        cap = stem
        if stem in diag:
            cap = f"{stem}    {diag[stem]}"
        im = Image.open(p).convert("RGB")
        if im.width < 40 or im.height < 40:
            sys.exit(f"{p} 只有 {im.size},不像一张图")
        pages.append(page(im, cap, f))

    pages[0].save(
        out, save_all=True, append_images=pages[1:], resolution=float(DPI), format="PDF"
    )
    print(f"{len(pages)} 页,写在 {out}")


if __name__ == "__main__":
    main()
