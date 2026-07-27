"""六面板图:分子内不可预知(论点 A)与盐带来的契约缺口(论点 B)。

上排 A:同一条模板产物是一片还是两片由底物决定;误差随环大小线性增长;
        两种模板写法没有哪一种对两类底物都给出正确的分子数。
下排 B:旁观离子被静默丢弃;盐在语料里的普遍程度;契约缺口的位置。

# 两处画图上的坑

**分子的标注不能交给 RDKit 画。** `MolsToGridImage` 的 legend 走的是 RDKit 自己的
字体栈,不含中文字体,标注会静默变成方框或空串。所以这里只让它画结构,文字一律
用 matplotlib 的坐标轴标题写。

**勾叉符号不能直接用。** 中文字体多半没有 U+2713 / U+2717 的字形,matplotlib 会
告警并画成方框。改用"正确 / 错"两个汉字加配色。

数据来自 measure_intra_salt.py 写出的 results/intra_salt.json。
"""

import argparse
import io
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from PIL import Image
from rdkit import Chem, RDLogger
from rdkit.Chem import Draw

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

OK, BAD, GREY = "#2f855a", "#c53030", "#4a5568"


def pick_font():
    from matplotlib import font_manager

    for name in ("PingFang SC", "Heiti SC", "STHeiti", "Arial Unicode MS", "Songti SC"):
        try:
            font_manager.findfont(name, fallback_to_default=False)
            return name
        except Exception:
            continue
    return None


def draw_mol(ax, smi, title, color="#1a202c", size=(460, 330)):
    """只让 RDKit 画结构,标注一律交给 matplotlib —— 理由见模块文档。"""
    m = Chem.MolFromSmiles(smi)
    img = Draw.MolToImage(m, size=size)
    if isinstance(img, bytes):
        img = Image.open(io.BytesIO(img))
    ax.imshow(np.asarray(img))
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_color("#e2e8f0")
    ax.set_title(title, fontsize=9.5, color=color, pad=4)


def header(fig, spec, title, sub=None, mono=None):
    """在面板自己的格子里画标题。

    不能用 `ax.text(..., y > 1)` 把标题顶到坐标轴外面 —— 那会越过 GridSpec 的
    格子边界,与总标题或相邻面板叠在一起。这里给标题单开一行不可见的坐标轴。
    """
    ax = fig.add_subplot(spec)
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0, 0.92, title, fontsize=11, fontweight="bold", va="top")
    # 副标题的位置要按等宽块**实际占几行**算。写死一个 y,等宽块一换行就压到
    # 副标题上 —— 而叠字在图里不报错,只是看不清。
    y = 0.58
    if mono:
        ax.text(0, y, mono, fontsize=7.2, family="monospace", color=GREY, va="top")
        y -= 0.30 * mono.count("\n") + 0.26
    else:
        y = 0.45
    if sub:
        ax.text(0, y, sub, fontsize=9, color=GREY, va="top")


def panel_a(fig, spec, d):
    sub = GridSpecFromSubplotSpec(
        3, 2, subplot_spec=spec, height_ratios=[0.5, 1, 1], hspace=0.30, wspace=0.05
    )
    # 模板串必须从数据里取,不能写死 —— 写死过一次,换了模板之后图上标的
    # 与实际量的就成了两条,而图看不出来
    tpl = d["templates"]["plain"].replace(">>", "\n  >>  ")
    header(
        fig, sub[0, :],
        "(a) 同一条断酰胺逆向模板(取自语料,出现 295 次),切成几片由底物决定",
        "模板作者在写模板时无从知道目标是不是环",
        tpl,
    )
    lac = next(r for r in d["ring_sweep"] if r["n"] == 7)
    cells = [
        ("CCC(=O)NCc1ccccc1", "底物:开链酰胺", "#1a202c"),
        ("CCC(=O)O.NCc1ccccc1", "→ 两个分子", OK),
        (lac["smiles"], "底物:七元内酰胺", "#1a202c"),
        ("NCCCCCC(=O)O", "→ 一个分子", OK),
    ]
    for k, (smi, title, col) in enumerate(cells):
        draw_mol(fig.add_subplot(sub[1 + k // 2, k % 2]), smi, title, col)


def panel_b(ax, d):
    ns = [r["n"] for r in d["ring_sweep"]]
    extra_plain = [r["rd_plain"]["heavy"] - r["want_heavy"] for r in d["ring_sweep"]]
    extra_group = [r["rd_grouped"]["heavy"] - r["want_heavy"] for r in d["ring_sweep"]]
    extra_og = [r["og_plain"]["heavy"] - r["want_heavy"] for r in d["ring_sweep"]]
    ax.plot(ns, extra_plain, "o-", color=BAD, lw=2, ms=6, label="RDKit,模板写作 >>A.B")
    ax.plot(ns, extra_group, "s--", color=GREY, lw=1.6, ms=5, label="RDKit,模板写作 >>(A.B)")
    ax.plot(ns, extra_og, "^-", color=OK, lw=2, ms=6, label="omgkit,两种写法均同")
    ax.axhline(0, color="#a0aec0", lw=0.8)
    ax.annotate(
        "多出的原子数 = n − 4\n即模板匹配的 4 个环原子之外的部分,\n被整份复制进了两片",
        xy=(15, 11), xytext=(5.4, 15.5), fontsize=9, color=BAD,
        arrowprops=dict(arrowstyle="->", color=BAD, lw=1),
    )
    ax.set_xlabel("内酰胺环大小 n")
    ax.set_ylabel("产物比底物多出的重原子数")
    ax.set_title("(b) 误差随环增大而线性增长", fontsize=11, fontweight="bold", loc="left")
    ax.legend(fontsize=8.5, loc="lower right", framealpha=0.95)
    ax.grid(alpha=0.25)
    ax.set_ylim(-2, 24)


def panel_c(ax, d):
    im = d["intermolecular"]
    lac = next(r for r in d["ring_sweep"] if r["n"] == 7)
    rows = [
        ("RDKit   >>A.B",
         (im["rd_plain"]["mols"] == 2, f"{im['rd_plain']['mols']} 片"),
         (lac["rd_plain"]["mols"] == 1, f"{lac['rd_plain']['mols']} 片,多 3 个原子")),
        ("RDKit   >>(A.B)",
         (im["rd_grouped"]["mols"] == 2, f"{im['rd_grouped']['mols']} 片"),
         (lac["rd_grouped"]["mols"] == 1, f"{lac['rd_grouped']['mols']} 片")),
        ("omgkit  两种写法",
         (im["og_plain"]["mols"] == 2, f"{im['og_plain']['mols']} 片"),
         (lac["og_plain"]["mols"] == 1, f"{lac['og_plain']['mols']} 片")),
    ]
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    ax.set_title("(c) 责任落在谁身上", fontsize=11, fontweight="bold", loc="left")
    ax.text(4.0, 9.0, "分子间底物\n应得 2 片", ha="center", fontsize=9.5, color=GREY)
    ax.text(7.8, 9.0, "分子内底物\n应得 1 片", ha="center", fontsize=9.5, color=GREY)
    for i, (name, inter, intra) in enumerate(rows):
        y = 6.9 - i * 1.75
        # 中文没有等宽字形,整体用 sans;只有模板写法那一段是 ASCII,对齐靠空格
        ax.text(0.05, y, name, fontsize=9.8, va="center")
        for x, (ok, txt) in ((4.0, inter), (7.8, intra)):
            col = OK if ok else BAD
            ax.add_patch(
                plt.Rectangle((x - 1.65, y - 0.58), 3.3, 1.16,
                              facecolor=col, alpha=0.13, edgecolor=col, lw=1.2)
            )
            ax.text(x, y, f"{'正确' if ok else '错'} · {txt}", ha="center", va="center",
                    fontsize=9.5, color=col, fontweight="bold")
    ax.text(
        0.05, 1.15,
        "没有哪一种模板写法对两类底物都给出正确的分子数。\n"
        "把分子数定义为“重写后图的连通分量数”,则两类同时正确 ——\n"
        "决定权因此从模板作者手里回到底物。",
        fontsize=9.5, color="#1a202c", linespacing=1.5,
    )


def panel_d(fig, spec, d):
    sub = GridSpecFromSubplotSpec(
        2, 2, subplot_spec=spec, height_ratios=[0.32, 1], hspace=0.12, wspace=0.05
    )
    header(
        fig, sub[0, :],
        "(d) 完全不连通的旁观组分被静默丢弃(明文约定)",
        "搬运是从已匹配原子出发的遍历,走不到孤立组分;两个实现同样如此。底物是拼的——语料里没有",
    )
    cells = [
        ("O=C1CCCCCN1.Cl", "底物:内酰胺 + 旁观 HCl(9 个重原子)", "#1a202c"),
        ("NCCCCCC(=O)O", "→ 两个引擎都给出这个,HCl 消失了", BAD),
    ]
    for k, (smi, title, col) in enumerate(cells):
        draw_mol(fig.add_subplot(sub[1, k]), smi, title, col, size=(560, 400))


def panel_e(ax, d, claim):
    """订正版:那 3625 条"含小离子/带电片段"里,旁观反离子有几条。

    原先这一格画的是"反应物侧含小离子或带电片段 = 7.2%",用来说盐很常见。
    那个数字把**真反应物**(甲醇钠、氨、鏻盐、格氏试剂)也算了进去。逐条对着
    记录的参与反应分子集合判过之后,旁观反离子是 0 条 —— 语料给不出这一档。
    """
    t = claim["tally"]
    n = t["记录"]
    items = [
        ("含小离子 /\n带电片段", t["含小离子/带电片段·合计"], "#cbd5e0"),
        ("其中:是\n参与反应的分子", t["含小离子/带电片段·是参与反应的"], "#2b6cb0"),
        ("其中:旁观反离子", t.get("含小离子/带电片段·旁观(真反离子那一档)", 0), BAD),
    ]
    xs = np.arange(len(items))
    vals = [100 * v / n for _, v, _ in items]
    bars = ax.bar(xs, vals, color=[c for _, _, c in items], width=0.58)
    for b, (_, v, _) in zip(bars, items):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.12,
                f"{v} 条\n{100 * v / n:.1f}%", ha="center", fontsize=8.8)
    ax.set_xticks(xs)
    ax.set_xticklabels([k for k, _, _ in items], fontsize=9)
    ax.set_ylabel(f"占 {n} 条记录的比例 / %")
    ax.set_ylim(0, max(vals) * 1.45)
    ax.set_title("(e) 订正:那 7.2% 全是真反应物,旁观反离子 0 条",
                 fontsize=11, fontweight="bold", loc="left")
    ax.grid(alpha=0.25, axis="y")


def panel_f(ax, d):
    ax.axis("off")
    ax.set_title("(f) 契约缺口出现在反应物侧,与分子内同源", fontsize=11,
                 fontweight="bold", loc="left")
    ax.text(
        0.0, 0.93,
        "盐在记录里是一个分子、多个组分:\n"
        "     Cl.NCc1ccccc1  是一个物种,不是两个反应物\n\n"
        "于是“N 个反应物模板 ↔ N 个输入分子”的契约,\n"
        "写不出同时碰阳离子与阴离子的模板 ——\n"
        "     季铵盐化的逆向      成盐 / 解离\n"
        "     抗衡离子交换        两性离子的质子转移\n\n"
        "而组分括号 (A.B) 不是无损的等价改写:\n"
        "它把“这两片本该是两个独立物种”这条信息吃掉,\n"
        "输出只剩连通分量,于是分不清\n"
        "     断出来的另一半   与   旁观的反离子\n\n"
        "两条都指向同一件事 ——\n"
        "片数不是模板的性质,是(模板, 底物)共同的性质。\n\n"
        "但要说清楚:以上是表达力论证,只依赖契约的形状。\n"
        "本语料量不出它 —— 抽取时多组分物种已被拆开,\n"
        "50016 条记录的输入分子全是单组分,见 (e)。",
        fontsize=9.5, va="top", linespacing=1.62, family="sans-serif",
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default=os.path.join(ROOT, "results", "intra_salt.json"))
    ap.add_argument("--claim", default=os.path.join(ROOT, "results", "salt_claim.json"),
                    help="audit_salt_claim.py 的输出,(e) 用它画订正后的拆分")
    ap.add_argument("--out", default=os.path.join(ROOT, "figures", "intra_and_salt.png"))
    args = ap.parse_args()

    font = pick_font()
    if font:
        plt.rcParams["font.sans-serif"] = [font]
        plt.rcParams["axes.unicode_minus"] = False

    d = json.load(open(args.src))
    claim = json.load(open(args.claim))
    fig = plt.figure(figsize=(17, 10.2))
    gs = GridSpec(2, 3, figure=fig, hspace=0.38, wspace=0.22,
                  left=0.035, right=0.985, top=0.905, bottom=0.055)
    panel_a(fig, gs[0, 0], d)
    panel_b(fig.add_subplot(gs[0, 1]), d)
    panel_c(fig.add_subplot(gs[0, 2]), d)
    panel_d(fig, gs[1, 0], d)
    panel_e(fig.add_subplot(gs[1, 1]), d, claim)
    panel_f(fig.add_subplot(gs[1, 2]), d)
    fig.suptitle(
        "反应模板应用:产物切成几片由底物而非模板决定 —— 分子内(上排)与旁观组分(下排)",
        fontsize=13.5, y=0.975,
    )
    fig.savefig(args.out, dpi=170)
    print(f"已写出 {args.out}")


if __name__ == "__main__":
    main()
