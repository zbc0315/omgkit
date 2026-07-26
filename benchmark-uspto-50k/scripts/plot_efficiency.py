"""效率图:四个子图 —— omgkit 正向 / omgkit 逆向 / rdkit 正向 / rdkit 逆向。

横轴运行时间,纵轴反应物的重原子数,每个散点一条反应,散点之上再叠分布等高线。

# 横轴为什么取对数

耗时跨三个数量级(中位数几十微秒、尾部几毫秒)。线性轴上 99% 的点会挤在最左边
一条竖线里,图就只剩尾部那几个点有信息。取对数之后主体与尾部同时看得见。

# 纵轴一律是**反应物**的重原子数

四个子图用同一个纵坐标才好横向比。逆向的输入其实是产物,但产物与反应物的重原子
数高度相关(差的是离去基团),而统一纵轴带来的可比性更值钱。另出一张按各自
**输入**规模作纵轴的补充图,免得这个取舍变成暗坑。

# 等高线是核密度估计

散点有五万个,重叠之后看不出密度。等高线用高斯核密度估计,在 log10(耗时) × 重原子
数这个平面上算 —— 与横轴的刻度一致,否则等高线会被拉成一条缝。
"""

import argparse
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

PANELS = [
    ("omgkit_fwd", "omgkit 正向"),
    ("omgkit_retro", "omgkit 逆向"),
    ("rdkit_fwd", "RDKit 正向"),
    ("rdkit_retro", "RDKit 逆向"),
]


def pick_font():
    """挑一个装得上中文的字体;找不到就退回英文标题。"""
    from matplotlib import font_manager

    for name in ("PingFang SC", "Heiti SC", "STHeiti", "Arial Unicode MS", "Songti SC"):
        try:
            font_manager.findfont(name, fallback_to_default=False)
            return name
        except Exception:
            continue
    return None


def load(path, ykey):
    xs = {k: [] for k, _ in PANELS}
    ys = {k: [] for k, _ in PANELS}
    hit = {k: [] for k, _ in PANELS}
    with open(path) as fh:
        for line in fh:
            r = json.loads(line)
            if "err" in r:
                continue
            for key, _ in PANELS:
                v = r.get(key)
                if not v or "t_min" not in v:
                    continue
                y = r[ykey["retro" if key.endswith("retro") else "fwd"]]
                if y <= 0:
                    continue
                xs[key].append(v["t_min"] * 1e6)  # 微秒
                ys[key].append(y)
                hit[key].append(bool(v["hit"]))
    return xs, ys, hit


def draw(xs, ys, out, title, ylabel, cn):
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.5), sharex=True, sharey=True)
    allx = np.concatenate([np.asarray(xs[k]) for k, _ in PANELS])
    xlim = (max(allx.min() * 0.7, 1.0), allx.max() * 1.4)
    ally = np.concatenate([np.asarray(ys[k]) for k, _ in PANELS])
    ylim = (0, np.percentile(ally, 99.9) * 1.05)

    for ax, (key, label) in zip(axes.ravel(), PANELS):
        x = np.asarray(xs[key], dtype=float)
        y = np.asarray(ys[key], dtype=float)
        lx = np.log10(x)
        # 重原子数是整数,直接画会叠成一条条横线;加 ±0.35 的抖动只影响观感,
        # 统计量(中位线、等高线)仍用原值算
        rng = np.random.default_rng(0)
        ax.scatter(
            x,
            y + rng.uniform(-0.35, 0.35, size=y.shape),
            s=1.4,
            alpha=0.08,
            color="#2b6cb0",
            linewidths=0,
            rasterized=True,
        )

        # 等高线:在 log10(耗时) × 重原子数 上做核密度估计
        n = min(len(lx), 20000)
        idx = np.linspace(0, len(lx) - 1, n).astype(int)
        kde = gaussian_kde(np.vstack([lx[idx], y[idx]]))
        gx = np.linspace(np.log10(xlim[0]), np.log10(xlim[1]), 160)
        gy = np.linspace(ylim[0], ylim[1], 160)
        gxx, gyy = np.meshgrid(gx, gy)
        dens = kde(np.vstack([gxx.ravel(), gyy.ravel()])).reshape(gxx.shape)
        levels = np.unique(np.percentile(dens[dens > 0], [55, 75, 88, 95, 98.5]))
        ax.contourf(10**gxx, gyy, dens, levels=[*levels, dens.max()], colors=[
            "#fff5eb", "#fee6ce", "#fdd0a2", "#fdae6b", "#f16913"
        ][: len(levels)], alpha=0.35)
        ax.contour(10**gxx, gyy, dens, levels=levels, colors="#c05621", linewidths=0.9)

        med = np.median(x)
        ax.axvline(med, color="#718096", ls="--", lw=1.0)
        ax.set_xscale("log")
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        rate = 100.0 * sum(1 for h in HITS[key] if h) / len(HITS[key])
        ax.set_title(
            f"{label}   n={len(x)}   中位 {med:.0f} µs   命中 {rate:.2f}%"
            if cn
            else f"{label}   n={len(x)}   median {med:.0f} us   hit {rate:.2f}%",
            fontsize=11,
        )
        ax.grid(alpha=0.25, which="both")

    for ax in axes[1]:
        ax.set_xlabel("单次 run_reactants 耗时 / µs(对数轴)" if cn else "run_reactants time / us (log)")
    for ax in axes[:, 0]:
        ax.set_ylabel(ylabel)
    fig.suptitle(title, fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out, dpi=170)
    print(f"已写出 {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default=os.path.join(ROOT, "results", "bench.jsonl"))
    ap.add_argument("--outdir", default=os.path.join(ROOT, "figures"))
    args = ap.parse_args()

    font = pick_font()
    if font:
        plt.rcParams["font.sans-serif"] = [font]
        plt.rcParams["axes.unicode_minus"] = False
    cn = font is not None

    global HITS
    # 主图:纵轴一律是反应物重原子数
    xs, ys, HITS = load(args.src, {"fwd": "n_heavy_r", "retro": "n_heavy_r"})
    draw(
        xs,
        ys,
        os.path.join(args.outdir, "efficiency.png"),
        "USPTO-50k:run_reactants 耗时 × 反应物规模" if cn else "USPTO-50k: run_reactants time x reactant size",
        "反应物重原子数" if cn else "reactant heavy atoms",
        cn,
    )
    # 补充图:纵轴换成各自**输入**的规模(逆向的输入是产物)
    xs2, ys2, HITS = load(args.src, {"fwd": "n_heavy_r", "retro": "n_heavy_p"})
    draw(
        xs2,
        ys2,
        os.path.join(args.outdir, "efficiency_by_input.png"),
        "补充:纵轴换成各自输入的规模(逆向的输入是产物)" if cn else "Supplement: y = input size",
        "输入分子重原子数" if cn else "input heavy atoms",
        cn,
    )


if __name__ == "__main__":
    main()
