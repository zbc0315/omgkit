#!/usr/bin/env python3
"""生成文档里的插图。**每一个分子都由 omgkit 自己画**,本脚本只做排版。

# 分工

`omgkit-depict` 的 `draw` 例子把每条 SMILES 画成一张独立的 SVG(零依赖、
无随机数,同一个分子每次逐字节相同)。本脚本把那些 SVG 摆进网格、或者串成
`A + B → C + D` 的反应式,加上加号与箭头 —— 那两个是排版符号,不是化学。

**不用任何别的化学库出图。** 文档里的结构式必须是这个库自己画出来的:
画错了要在文档上看得见,而不是被另一个库的图盖住。

用法:

    python3 docs/figures/make_figures.py

输出到 `docs/assets/`。要改图就改下面那张清单,重跑一遍。
"""
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs" / "assets"

# 用 ChemDraw 那套规范:键长 30 pt,标签只占键长的三分之一,图看着舒展。
# ACS 1996 是 14.4 pt,适合排进论文正文,放到网页上标签会显得很挤。
STYLE = "cd"

# ── 图里用到的全部分子 ───────────────────────────────────────────────
# 名字只在本脚本内部用来引用,不出现在图上。
MOLECULES = {
    # 画廊:每一个都考一类布局 —— 稠环、桥环、糖、甾体、β-内酰胺
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "caffeine": "CN1C=NC2=C1C(=O)N(C)C(=O)N2C",
    "penicillin": "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O",
    "cholesterol": "CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C",
    "glucose": "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
    "morphine": "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5",  # 桥环,故意的
    "nicotine": "CN1CCC[C@H]1c1cccnc1",
    "camphor": "CC1(C)[C@@H]2CC[C@@]1(C)C(=O)C2",  # 桥环,故意的
    "ascorbic": "OC[C@H](O)[C@H]1OC(=O)C(O)=C1O",
    # 立体:两组对照
    "L-alanine": "C[C@H](N)C(=O)O",
    "D-alanine": "C[C@@H](N)C(=O)O",
    "trans-butene": "C/C=C/C",
    "cis-butene": "C/C=C\\C",
    # 酯化:酸 + 醇 → 酯 + 水
    "benzoic-acid": "OC(=O)c1ccccc1",
    "ethanol": "CCO",
    "ethyl-benzoate": "CCOC(=O)c1ccccc1",
    "water": "O",
    # Boc 脱保护:形式副产物是叔丁基碳酸,不是实际拿到的二氧化碳加异丁烯
    "boc-amine": "CC(C)(C)OC(=O)NCc1ccccc1",
    "free-amine": "NCc1ccccc1",
    "tbu-carbonic": "CC(C)(C)OC(=O)O",
}


def draw_all(tmp: Path) -> None:
    """把清单里每条 SMILES 交给 omgkit 画一遍。"""
    args = [f"{name}={smi}" for name, smi in MOLECULES.items()]
    cmd = [
        "cargo", "run", "-q", "--release", "-p", "omgkit-depict",
        "--features", "raster", "--example", "draw", "--", str(tmp), *args,
    ]
    r = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"omgkit 出图失败:\n{r.stderr}")
    # **画不干净的要说出来。** `draw` 每行末尾报四个计数(退化/未解冲突/交叉/
    # 未画手性)。文档里的图应当全是 0 —— 不是 0 就该换分子或者去修布局,
    # 而不是把一张读不出构型的图摆进文档。
    #
    # 桥环那两个是例外:`degraded.svg` 专门拿它们展示"画不好会说出来"。
    # 把预期之内的那两条排掉,否则每次跑都有两行警告,人很快就不看警告了 ——
    # 那时真的出问题也照样被划过去。
    expected = {"camphor", "morphine"}
    bad = [ln for ln in r.stdout.splitlines()
           if STYLE in ln.split() and ln.split()[0] not in expected
           and not ln.rstrip().endswith("退化0 未解冲突0 交叉0 未画手性0")]
    for ln in bad:
        print(f"  ⚠ 画得不干净:{ln}")


def load(tmp: Path, name: str):
    """读一张 omgkit 画的 SVG,拆成 `(宽, 高, 内容)`。

    外层 `<svg>` 与白底那一行都去掉 —— 拼进大图之后,底色由大图统一铺,
    每张小图各铺一块白的话,图与图之间会看到接缝。
    """
    text = (tmp / f"{name}.{STYLE}.svg").read_text(encoding="utf-8")
    m = re.search(r'width="([\d.]+)" height="([\d.]+)"', text)
    w, h = float(m.group(1)), float(m.group(2))
    body = re.sub(r"^<svg[^>]*>\n", "", text)
    body = re.sub(r"^<rect[^>]*fill=\"#fff\"/>\n", "", body)
    body = body.replace("</svg>", "").strip()
    return w, h, body


def place(body: str, dx: float, dy: float) -> str:
    return f'<g transform="translate({dx:.2f},{dy:.2f})">\n{body}\n</g>\n'


def wrap(items: str, w: float, h: float) -> str:
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{w:.0f}" height="{h:.0f}" '
        f'viewBox="0 0 {w:.2f} {h:.2f}">\n'
        f'<rect width="{w:.2f}" height="{h:.2f}" fill="#fff"/>\n'
        f"{items}</svg>\n"
    )


def cap_width(text: str, size: float = 9.0) -> float:
    """图注大概多宽。汉字按一个字号、拉丁字母按 0.55 个字号估。

    估宽只用来给画布留地方,估宽一点点不要紧,估窄了最右那一格的图注会被
    切掉半个字 —— 那是实测撞见过的。
    """
    return sum(size if ord(c) > 0x2E80 else size * 0.55 for c in text)


def caption(text: str, cx: float, y: float, size: float = 9.0) -> str:
    """图下面那行小字。字体与 omgkit 画标签用的是同一族,免得一张图里两种字。"""
    return (
        f'<text x="{cx:.2f}" y="{y:.2f}" font-family="Helvetica,Arial,sans-serif" '
        f'font-size="{size:.1f}" fill="#333" text-anchor="middle">{text}</text>\n'
    )


def grid(tmp: Path, cells, cols: int, lang: int, pad: float = 16.0, cap_h: float = 16.0) -> str:
    """把若干 `(分子名, (英文图注, 中文图注))` 摆成网格,每格居中。"""
    loaded = [(load(tmp, name), cap[lang]) for name, cap in cells]
    rows = [loaded[i:i + cols] for i in range(0, len(loaded), cols)]
    col_w = max(max(w, cap_width(cap)) for (w, _, _), cap in loaded) + pad * 2
    row_hs = [max(h for (_, h, _), _ in r) + pad * 2 + cap_h for r in rows]

    items, y = "", 0.0
    for r, rh in zip(rows, row_hs):
        for i, ((w, h, body), cap) in enumerate(r):
            cx = i * col_w + col_w / 2
            items += place(body, cx - w / 2, y + (rh - cap_h - h) / 2)
            if cap:
                items += caption(cap, cx, y + rh - 5)
        y += rh
    return wrap(items, col_w * cols, y)


def reaction(tmp: Path, left, right, lang: int, pad: float = 18.0, cap_h: float = 16.0) -> str:
    """串成 `A + B → C + D`。两侧都是 `(分子名, (英文图注, 中文图注))` 的列表。"""
    ARROW, GAP, PLUS = 44.0, 14.0, 16.0

    def side(names):
        return [(load(tmp, n), cap[lang]) for n, cap in names]

    ls, rs = side(left), side(right)
    height = max(h for (_, h, _), _ in ls + rs) + pad * 2 + cap_h
    mid = (height - cap_h) / 2

    items, x = "", pad
    def emit(entries, x):
        nonlocal items
        for i, ((w, h, body), cap) in enumerate(entries):
            if i:
                items += caption("+", x + PLUS / 2, mid + 4, 14.0)
                x += PLUS + GAP
            items += place(body, x, mid - h / 2)
            if cap:
                items += caption(cap, x + w / 2, height - 5)
            x += w + GAP
        return x

    x = emit(ls, x)
    # 箭头:一根线加一个实心三角。图上唯一不是 omgkit 画的东西,是排版符号。
    y0 = mid
    items += (f'<line x1="{x:.2f}" y1="{y0:.2f}" x2="{x + ARROW - 6:.2f}" y2="{y0:.2f}" '
              f'stroke="#333" stroke-width="1.2"/>\n')
    items += (f'<path d="M{x + ARROW:.2f},{y0:.2f} L{x + ARROW - 8:.2f},{y0 - 3.5:.2f} '
              f'L{x + ARROW - 8:.2f},{y0 + 3.5:.2f} Z" fill="#333"/>\n')
    x += ARROW + GAP
    x = emit(rs, x)
    # 最后一格的图注是居中在那个分子上的,分子窄而图注宽时会伸到画布外面
    (last_w, _, _), last_cap = rs[-1]
    overhang = max(0.0, (cap_width(last_cap) - last_w) / 2)
    return wrap(items, x + pad - GAP + overhang, height)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        print("omgkit 画图中…")
        draw_all(tmp)

        # 一张图两套图注:英文页与中文页各用各的,分子本身是同一批 SVG。
        # 中文名对中国的化学家更直接,英文名是文档站的主语言 —— 两边都不将就。
        for lang, suffix in ((0, ""), (1, ".zh")):
            figures = {
                # 画廊:每一格考一类布局 —— 稠环、含氮杂环、β-内酰胺、糖、烯醇内酯
                f"gallery{suffix}.svg": grid(tmp, [
                    ("aspirin", ("Aspirin", "阿司匹林")),
                    ("caffeine", ("Caffeine", "咖啡因")),
                    ("nicotine", ("Nicotine", "尼古丁")),
                    ("penicillin", ("Penicillin G", "青霉素 G")),
                    ("glucose", ("Glucose", "葡萄糖")),
                    ("ascorbic", ("Ascorbic acid", "抗坏血酸")),
                ], cols=3, lang=lang),
                # 大骨架:甾体四环 + 侧链,八个手性中心全画出了构型
                f"cholesterol{suffix}.svg": grid(tmp, [("cholesterol", ("", ""))], cols=1, lang=lang),
                # 桥环:平面上没有好解,omgkit 如实报退化而不是假装画好了
                f"degraded{suffix}.svg": grid(tmp, [
                    ("camphor", ("Camphor · reports degraded=1", "樟脑 · 报 degraded=1")),
                    ("morphine", ("Morphine · reports degraded=1", "吗啡 · 报 degraded=1")),
                ], cols=2, lang=lang),
                # 立体:楔形与顺反各一组对照
                f"stereo{suffix}.svg": grid(tmp, [
                    ("L-alanine", ("L-alanine", "L-丙氨酸")),
                    ("D-alanine", ("D-alanine", "D-丙氨酸")),
                    ("trans-butene", ("trans-2-butene", "反-2-丁烯")),
                    ("cis-butene", ("cis-2-butene", "顺-2-丁烯")),
                ], cols=4, lang=lang),
                # 反应模板:模板只描述反应中心,分子其余部分自动跟着走
                f"esterification{suffix}.svg": reaction(tmp,
                    [("benzoic-acid", ("Benzoic acid", "苯甲酸")),
                     ("ethanol", ("Ethanol", "乙醇"))],
                    [("ethyl-benzoate", ("Ethyl benzoate", "苯甲酸乙酯")),
                     ("water", ("Water · reconstructed", "水 · 收口得到"))], lang=lang),
                # 副产物收口:模板只写主产物,丢掉的原子由原子账收口成分子
                f"byproduct{suffix}.svg": reaction(tmp,
                    [("boc-amine", ("Boc-protected benzylamine", "Boc 保护的苄胺"))],
                    [("free-amine", ("Benzylamine", "苄胺")),
                     ("tbu-carbonic", ("tert-Butyl carbonic acid · formal byproduct",
                                       "叔丁基碳酸 · 形式副产物"))], lang=lang),
            }
            for name, svg in figures.items():
                (OUT / name).write_text(svg, encoding="utf-8")
                print(f"  {name}  {len(svg) // 1024} KB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
