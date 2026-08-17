#!/usr/bin/env python3
"""把并排比对图拼成一份 `gallery.pdf`,供人眼一次翻完。

    cargo run -p omgkit-depict --release --features raster --example draw -- <目录>
    python3 harness/compare_rdkit.py <目录>
    python3 harness/make_gallery.py <目录>

# 为什么这一步要有脚本

前两步一直有,**第三步没有** —— 先前那份 `gallery.pdf` 是随手拼的,全树搜不到
任何生成它的代码。于是出图流水线的最后一环不可复现:图变了没人知道该怎么重出,
而重出来的东西与上一份是不是同一个口径也说不清。这个文件补上那一环。

# 页序按 `mols.tsv`,不按文件名

`mols.tsv` 是 `draw` 例子写出来的清单,顺序就是它内置那 17 个分子的排法 ——
那个排法是按"能考出什么"挑的(芳环、稠环、桥环、糖、甾体、β-内酰胺、对映体、
顺反,每类一个代表)。按文件名排会把这个次序打散成字母序,翻起来就没有脉络了。

# 两个方向都要查,不能只查一头

第一版只查"清单里列了 → 盘上有没有",**反过来那一头没查**。实测:把 `mols.tsv`
截成 1 行(`draw` 带分子名参数跑一次就是这个效果,`docs/guide/depict.md` 正式
记着这个用法),盘上 34 张图原封不动,拼出来是**一本 3 页的册子,退出码 0**
—— 32 页无声蒸发。而第一版的注释正写着"不许静默少拼几页"。

所以两头都查:清单里有而盘上没有,报;盘上有而清单没列,也报。

# 重跑必须逐字节相同

`harness/README.md` 的准入门槛写着"重跑一遍必须逐字节相同",`docs/dev/building.md`
更把"跑两遍比一比"当作发现迭代序缺陷的唯一手段。第一版首页盖了挂钟时间,加上
Pillow 自己往 PDF 里写的 `/CreationDate`,两次跑**差 1158642 字节**(实测)。

现在:首页那行"生成"取的是**输入图里最新的那个 mtime**,不是挂钟;写 PDF 时
把两个日期钉成同一个值。两次跑逐字节相同。

顺带修掉一个更隐蔽的毛病:首页的版本号取的是**拼册时**的仓库版本,而册子里画的
是**跑 draw 那一刻**的布局。两者可以差很远 —— 带哈希的假出处比没有出处更骗人。
所以那一行明写"拼册时",而"图出自何时"由 mtime 那一行负责。
"""

import subprocess
import sys
import time
from pathlib import Path

from PIL import Image, ImageDraw

from compare_rdkit import font

#: 首页文字四周的内缩(px)。首页是这个脚本自己画的,与比对图的 `PAD`/`BAR`
#: 没有关系 —— 第一版的注释说它"与比对图同一个 300 dpi 口径",那是编的:
#: 比对图的留白是 `compare_rdkit.PAD = 24`,而这份 PDF 是按 72 dpi 写出来的
#: (MediaBox 与像素数一一对应)。
COVER_PAD = 48
#: 首页字号。比对图自己的标题是 26 px、副标题 19 px(`compare_rdkit.main`),
#: 首页比它们大一档,翻开第一眼就读得清。
COVER_PT = 44
#: 首页宽度(px)。**不跟最宽的那一页走** —— 那会得到一个 3177×256 的畸形长条,
#: 与后面 415×153 起步的页面完全不成体系。
COVER_W = 1600


def git_describe(repo: Path) -> str:
    """仓库当前的版本描述。取不到就如实说取不到,不编一个。"""
    try:
        out = subprocess.run(
            ["git", "-C", str(repo), "describe", "--always", "--dirty", "--abbrev=12"],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
    except (OSError, subprocess.SubprocessError) as e:
        return f"(取不到:{e})"
    if out.returncode != 0:
        return "(不是 git 仓库,或 git 不可用)"
    return out.stdout.strip() or "(git describe 没有输出)"


def cover(lines: list[str]) -> Image.Image:
    """首页:白底,左上角逐行写下出处。"""
    f = font(COVER_PT)
    step = round(COVER_PT * 1.6)
    img = Image.new("RGB", (COVER_W, COVER_PAD * 2 + step * len(lines)), "white")
    d = ImageDraw.Draw(img)
    for i, line in enumerate(lines):
        d.text((COVER_PAD, COVER_PAD + step * i), line, fill="black", font=f)
    return img


def read_tsv(path: Path, cols: int) -> list[list[str]]:
    """读一张制表符分隔的表,列数不对就当场报错 —— 与 `compare_rdkit` 同一套做法。"""
    if not path.exists():
        sys.exit(f"找不到 {path} —— 先跑 example draw 生成它")
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) != cols:
            sys.exit(f"{path} 这行不是 {cols} 列:{line!r}")
        rows.append(parts)
    if not rows:
        sys.exit(f"{path} 是空的")
    return rows


def main() -> None:
    d = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
    repo = Path(__file__).resolve().parent.parent

    names = [row[0] for row in read_tsv(d / "mols.tsv", 2)]
    if len(set(names)) != len(names):
        dup = sorted({n for n in names if names.count(n) > 1})
        sys.exit(f"{d / 'mols.tsv'} 里有重名,册子会重页:{dup}")
    styles = read_tsv(d / "styles.tsv", 4)

    want = [(n, s[0]) for n in names for s in styles]
    paths = [d / f"{n}.{tag}.compare.png" for n, tag in want]

    # **两个方向都查。** 单查一头就会静默少拼 —— 见模块文档。
    missing = [p for p in paths if not p.exists()]
    if missing:
        sys.exit(
            f"缺 {len(missing)} 张比对图,先把前两步跑完:"
            + "".join(f"\n  {p.name}" for p in missing[:10])
        )
    extra = sorted(set(d.glob("*.compare.png")) - set(paths))
    if extra:
        sys.exit(
            f"盘上有 {len(extra)} 张比对图不在 mols.tsv × styles.tsv 里 —— "
            f"清单陈旧还是目录里混了别的批次?册子只会收 {len(paths)} 张:"
            + "".join(f"\n  {p.name}" for p in extra[:10])
        )

    # **`.convert()` 不能省。** 它把解码提前到这里 —— 截断或损坏的 PNG 当场就炸,
    # 而不是等 `save()` 写到一半、在原地留下一个看着正常的半截 PDF。
    pages = [Image.open(p).convert("RGB") for p in paths]

    # 出处只取自输入,不取挂钟 —— 重跑才逐字节相同,见模块文档
    stamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(max(p.stat().st_mtime for p in paths)))
    head = cover(
        [
            "omgkit 出图册",
            "omgkit 与 RDKit 并排,同一个键长",
            "",
            f"图出自   {stamp}",
            f"拼册时仓库版本   {git_describe(repo)}",
            f"内容   {len(names)} 个分子 × {len(styles)} 套规范 = {len(pages)} 页",
            "",
        ]
        + [f"{tag}   {full}   键长 {bl} pt   字号 {fs} pt" for tag, full, bl, fs in styles]
    )

    # 先写临时文件再改名:中途失败不会在原地留下一个半截的 `gallery.pdf`
    pdf = d / "gallery.pdf"
    tmp = d / "gallery.pdf.part"
    head.save(
        tmp,
        # 临时文件后缀是 `.part`,PIL 靠扩展名猜不出格式,得明说
        format="PDF",
        save_all=True,
        append_images=pages,
        # 钉死,否则 Pillow 会往里写 `time.gmtime()`,两次跑就对不上
        creationDate="D:20000101000000Z",
        modDate="D:20000101000000Z",
    )
    tmp.replace(pdf)
    print(f"{pdf}:{1 + len(pages)} 页(首页 1 + 比对图 {len(pages)});图出自 {stamp}")


if __name__ == "__main__":
    main()
