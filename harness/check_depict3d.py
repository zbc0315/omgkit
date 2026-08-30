#!/usr/bin/env python3
"""**三维分子图的外部判官** —— 从我们吐出来的那段 SVG 反读,拿外部实现比。

    cargo run -p omgkit-depict --release --example dump_depict3d -- \\
        harness/corpus/large.smi 400 > /tmp/three.jsonl
    .venv/bin/python harness/check_depict3d.py /tmp/three.jsonl

# 为什么读 SVG,而不是读我们导出的中间量

导出的 jsonl 里同时有 `placed`(每个原子落在画布哪里)和 `svg`。判官**以 SVG
为准**:`placed` 是我们自己算的,拿它当真值等于自己跟自己比。SVG 才是交到读图
的人手里的那个东西 —— 圆画在哪、什么颜色、谁盖着谁,全在里面。`placed` 只被
判一件事:它与图里的圆对不对得上(那是个公开 API,对不上就是错的)。

# 外部真值分别是谁

| 判的东西 | 外部真值 |
|---|---|
| 主轴(视角) | RDKit `ComputeCanonicalTransform(ignoreHs=False)`,**以及** numpy `eigh` |
| 球半径 | RDKit 周期表的 `GetRvdw` |
| 元素颜色 | `harness/params/jmol_colors.tsv`(Jmol 的表,两源核对过) |
| 正交投影、深度序、并排方向 | numpy 独立复算 |

RDKit 与 numpy 两条独立的路**必须互相吻合**才用得上,判官第一件事就是核这个 ——
两者不一致时它们谁都不能当真值,当场退非零。

# 判什么

1. 视角矩阵是真旋转(`RᵀR = I`、`det = +1`)。**镜像会把每一个手性中心画反,
   而图上看不出来**,所以这一条是硬的。
2. 视角的三根轴与外部主轴逐根同向(允许整根反号)。主轴简并的分子跳过 ——
   那时轴本来就不唯一,跳过的个数会打印出来。
3. 图里的圆是坐标的**正交投影**:任意两个圆的位移 = 比例尺 × 旋转后的位移
   (y 翻向)。按位移比,所以不需要知道我们把画布平移了多少。
4. 球半径 = 样式的 vdW 比例 × RDKit 的 `GetRvdw` × 比例尺。
5. 球颜色 = Jmol 表。
6. **画序是深度序**:投影上重叠的两个球,靠前的那个必须画在后面。深度由判官
   从坐标和旋转矩阵自己算。
7. **键的两半各随自己那一端的元素上色**。
8. **多重键并排的方向垂直于键在屏幕上的投影**(线框样式上判,那一档没有球,
   起点不必截,几何是闭式的)。
9. 把原子倒序重编号,SVG **逐字节相同**。
10. 图元一个都不出画布。
11. `placed` 与图里的圆一致。

# 分母闸

每一条都带一个下限。**喂空的判据不许打印"全部通过"** —— 回归常从分母那侧
进来:导出脚本换了个字段名、语料读空了、样式表少了一档,判据会一条都没得判
而全绿。下限的数取自当前实测的三分之一左右,留出语料波动的余量。
"""

from __future__ import annotations

import argparse
import collections
import json
import pathlib
import re
import sys

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolTransforms as rdT
from rdkit.Geometry import Point3D

# ---- 分母闸 ----
MIN_MOLECULES = 300
MIN_CIRCLES = 5_000
MIN_OVERLAP_PAIRS = 500
MIN_STICK_HALVES = 5_000
MIN_MULTIPLE_BONDS = 100
MIN_ELEMENTS = 5
MIN_AXES_COMPARED = 200

# 投影/半径判等的容差(磅)。SVG 里的坐标是按 `{:.2}` 写的,所以本底误差是
# 半个百分位 = 0.005;取 0.02 留三倍余量,而真出错时偏差是整整一个原子的量级。
TOL_PT = 0.02
# 深度差小于它就当成平局,不判次序。
#
# **不引用实现里的量化精度**(那是被测的东西,一起改就永远打不红),取一个
# 独立的、化学上毫无意义的线:1e-4 Å 是一根键的十万分之一。比这还细地追究
# "谁在前面",量的是浮点噪声,不是遮挡。平局的对数会打印出来。
DEPTH_TIE = 1e-4

# 主轴同向的容差:|cos| 与 1 的差。
TOL_AXIS = 1e-6
# 主轴间隔小于这个(相对最大特征值)就认为简并,与 Rust 侧的 DEGENERATE_TOL
# **不共享常量是故意的**:判据的期望值引用被测常量的话,那个常量一改两边一起
# 动,变异永远打不红。这里取一个更宽的线,宁可多跳过几个也不拿不唯一的轴去比。
TOL_DEGEN = 1e-4

CIRCLE_RE = re.compile(
    r'<circle cx="([-\d.]+)" cy="([-\d.]+)" r="([\d.]+)" fill="url\(#s([0-9a-f]{6})\)"'
)
LINE_RE = re.compile(
    r'<line x1="([-\d.]+)" y1="([-\d.]+)" x2="([-\d.]+)" y2="([-\d.]+)" '
    r'stroke="url\(#c[^"]*_([0-9a-f]{6})\)" stroke-width="([\d.]+)"'
)


def load_palette(root: pathlib.Path) -> dict[int, tuple[int, int, int]]:
    """Jmol 的表。表外的元素(以及通配原子)是 deeppink。"""
    out: dict[int, tuple[int, int, int]] = {}
    path = root / "harness/params/jmol_colors.tsv"
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        z, _sym, hexv = line.split("\t")
        out[int(z)] = tuple(int(hexv[i:i + 2], 16) for i in (0, 2, 4))  # type: ignore[misc]
    return out


UNKNOWN_COLOR = (0xFF, 0x14, 0x93)
# omgkit 给"元素表里没有范德华半径"的原子顶的半径,见 three.rs 的 UNKNOWN_RVDW。
UNKNOWN_RVDW = 1.7

# **四套样式的参数,判官自己存一份。**
#
# 起初这里是从导出的 jsonl 里读 `ball_vdw_frac` 的 —— 那是**被测的常量**,
# 一改两边一起动:实测把球棍的球半径从 23% 改成 25%,判官一声不吭地全绿。
# 判据的期望值不许引用被测常量。
#
# 数出自 Jmol 自己文档里的 standard rendering styles(Jmol Wiki, *Rendering*):
#   spacefill 100% / wireframe 0.15 + spacefill 23% / wireframe 0.3 / wireframe 0.01
# 比例尺与并排间距不是 Jmol 的,是本库自己定的,一并钉住 —— 动它们要过这一关,
# 那正是"改样式必须是有意的"这句话的落点。
STYLES = {
    "space-filling":  {"ball": 1.00, "stick": 0.00, "spacing": 0.00, "scale": 24.0},
    "ball-and-stick": {"ball": 0.23, "stick": 0.15, "spacing": 0.35, "scale": 36.0},
    "stick":          {"ball": 0.00, "stick": 0.30, "spacing": 0.00, "scale": 36.0},
    "wireframe":      {"ball": 0.00, "stick": 0.01, "spacing": 0.25, "scale": 36.0},
}

# 键级 → 并排画几根圆柱。芳香键一根,见 three.rs 的 `cylinders`。
CYLINDERS = {"Double": 2, "Triple": 3, "Quadruple": 3}


def parse_svg(svg: str) -> tuple[list, list]:
    """按**出现次序**取出圆与线 —— 次序就是画序,这一条是要判的东西之一。"""
    body = svg.split("</defs>", 1)[-1]
    circles = [
        (float(m[1]), float(m[2]), float(m[3]),
         tuple(int(m[4][i:i + 2], 16) for i in (0, 2, 4)), m.start())
        for m in CIRCLE_RE.finditer(body)
    ]
    lines = [
        (float(m[1]), float(m[2]), float(m[3]), float(m[4]),
         tuple(int(m[5][i:i + 2], 16) for i in (0, 2, 4)), float(m[6]), m.start())
        for m in LINE_RE.finditer(body)
    ]
    return circles, lines


def external_axes(coords: np.ndarray):
    """外部主轴(降序),两条独立的路各算一遍。

    返回 `(numpy 的轴, RDKit 的轴, numpy 的特征值)`。前两项都是按特征值降序排的
    行向量;第三项用来判断哪两根轴简并 —— 简并方向上两条路本来就可以不一致,
    那不算矛盾。
    """
    ctr = coords.mean(axis=0)
    d = coords - ctr
    vals, vecs = np.linalg.eigh(d.T @ d)
    np_axes = vecs[:, ::-1].T
    np_vals = vals[::-1]

    mol = Chem.RWMol()
    for _ in range(len(coords)):
        mol.AddAtom(Chem.Atom(6))
    conf = Chem.Conformer(len(coords))
    for i, p in enumerate(coords):
        conf.SetAtomPosition(i, Point3D(*(float(x) for x in p)))
    mol.AddConformer(conf)
    rd = np.array(rdT.ComputeCanonicalTransform(mol.GetConformer(), None, False, False))[:3, :3]
    return np_axes, rd, np_vals


def check(path: pathlib.Path, root: pathlib.Path) -> int:
    palette = load_palette(root)
    pt = Chem.GetPeriodicTable()
    bad: list[str] = []
    n_mol = 0
    n_circles = 0
    n_overlap = 0
    n_halves = 0
    n_multiple = 0
    n_axes = 0
    n_degenerate = 0
    n_ambiguous_half = 0
    n_depth_tie = 0
    n_coincident = 0
    elements: set[int] = set()
    styles_seen: set[str] = set()

    def fail(msg: str) -> None:
        if len(bad) < 30:
            bad.append(msg)

    for raw in path.read_text().splitlines():
        if not raw.strip():
            continue
        rec = json.loads(raw)
        smi = rec["smiles"]
        z = rec["z"]
        coords = np.array(rec["coords"], dtype=float)
        bonds = rec["bonds"]
        elements.update(z)
        n_mol += 1

        np_axes, rd_axes, np_vals = external_axes(coords)
        # 两条外部路先互核。不一致的话谁都当不了真值。
        for k in range(3):
            if abs(abs(float(np_axes[k] @ rd_axes[k])) - 1.0) > 1e-5:
                # 简并方向上两者可以不一致,那不算矛盾
                gap = min(
                    abs(np_vals[k] - np_vals[j]) for j in range(3) if j != k
                ) / max(abs(np_vals[0]), 1e-300)
                if gap > TOL_DEGEN:
                    fail(f"{smi}: numpy 与 RDKit 的第 {k} 根主轴不一致 —— 两个外部真值打架了")

        for style_name, st in rec["styles"].items():
            styles_seen.add(style_name)
            spec = STYLES.get(style_name)
            if spec is None:
                fail(f"{smi}: 冒出来一档判官不认识的样式 {style_name!r}")
                continue
            # 导出的那几个数与判官自己那份必须一致 —— 不一致就是样式表被动过,
            # 而下面所有的期望值都建立在判官这一份上。
            for key, got in (("ball", st["ball_vdw_frac"]), ("stick", st["stick_radius_a"]),
                             ("spacing", st["spacing_a"]), ("scale", st["scale"])):
                if abs(got - spec[key]) > 1e-12:
                    fail(f"{smi}/{style_name}: 样式的 {key} 是 {got},判官那份写的是 {spec[key]}")
            st = dict(st)
            st["ball_vdw_frac"] = spec["ball"]
            st["stick_radius_a"] = spec["stick"]
            st["spacing_a"] = spec["spacing"]
            st["scale"] = spec["scale"]
            rot = np.array(st["rot"], dtype=float)
            centre = np.array(st["centre"], dtype=float)
            scale = st["scale"]
            svg = st["svg"]

            # ① 真旋转
            if not np.allclose(rot @ rot.T, np.eye(3), atol=1e-9):
                fail(f"{smi}/{style_name}: 视角矩阵不正交")
            det = float(np.linalg.det(rot))
            if abs(det - 1.0) > 1e-9:
                fail(f"{smi}/{style_name}: 行列式 {det:.6f},不是 +1 —— 分子被镜像了")

            # ② 主轴与外部一致
            if st["degenerate"]:
                n_degenerate += 1
            else:
                gaps = [
                    abs(np_vals[k] - np_vals[k + 1]) / max(abs(np_vals[0]), 1e-300)
                    for k in range(2)
                ]
                if min(gaps) > TOL_DEGEN:
                    n_axes += 1
                    for k in range(3):
                        cos = float(rot[k] @ np_axes[k])
                        if abs(abs(cos) - 1.0) > TOL_AXIS:
                            fail(
                                f"{smi}/{style_name}: 第 {k} 根轴与外部主轴差 "
                                f"|cos|={abs(cos):.9f}"
                            )
                        cos_rd = float(rot[k] @ rd_axes[k])
                        if abs(abs(cos_rd) - 1.0) > 1e-5:
                            fail(f"{smi}/{style_name}: 第 {k} 根轴与 RDKit 主轴不同向")

            proj = (coords - centre) @ rot.T            # (右, 上, 前)
            depth = proj[:, 2]

            circles, lines = parse_svg(svg)

            # ③④⑤ 圆:位置、半径、颜色
            if st["ball_vdw_frac"] > 0:
                if len(circles) != len(z):
                    fail(f"{smi}/{style_name}: 图里 {len(circles)} 个圆,分子有 {len(z)} 个原子")
                    continue
                obs = np.array([[c[0], c[1]] for c in circles])
                want_rel = np.stack([proj[:, 0] * scale, -proj[:, 1] * scale], axis=1)
                # 圆的次序是画序(按深度),不是原子序 —— 按位置认原子
                groups = match_circles(obs, want_rel)
                if groups is None:
                    fail(f"{smi}/{style_name}: 图里的圆认不回原子 —— 投影对不上")
                    continue
                n_circles += len(circles)

                def want_ball(atom: int) -> tuple[float, tuple[int, int, int]]:
                    zi = int(z[atom])
                    rvdw = pt.GetRvdw(zi) if zi > 0 else 0.0
                    if rvdw <= 0.0:
                        rvdw = UNKNOWN_RVDW
                    return (
                        st["ball_vdw_frac"] * rvdw * scale,
                        palette.get(zi, UNKNOWN_COLOR),
                    )

                # ④⑤ 半径与颜色。重合的一组比多重集,单个的就是逐项比。
                pos_of: dict[int, int] = {}
                for cidx, aidx in groups:
                    # 半径带容差比:SVG 写的是两位小数,而期望值是全精度算出来的。
                    # 颜色是逐字节比 —— 那是整数,没有容差可言。
                    got = sorted((circles[c][3], circles[c][2]) for c in cidx)
                    exp = sorted((want_ball(a)[1], want_ball(a)[0]) for a in aidx)
                    same = len(got) == len(exp) and all(
                        g[0] == e[0] and abs(g[1] - e[1]) <= 0.006
                        for g, e in zip(got, exp)
                    )
                    if not same:
                        fail(
                            f"{smi}/{style_name}: 原子 {aidx} 处画出来的(颜色,半径)"
                            f"是 {got},该是 {[(c, round(r, 3)) for c, r in exp]}"
                        )
                    if len(aidx) == 1:
                        pos_of[aidx[0]] = cidx[0]
                    else:
                        n_coincident += len(aidx)

                # ⑥ 画序 = 深度序(只在位置不重合的原子之间判)
                for i in range(len(z)):
                    for j in range(i):
                        if i not in pos_of or j not in pos_of:
                            continue
                        gap = abs(depth[i] - depth[j])
                        if gap < DEPTH_TIE:
                            n_depth_tie += 1
                            continue
                        ri, rj = want_ball(i)[0], want_ball(j)[0]
                        if float(np.hypot(*(want_rel[i] - want_rel[j]))) >= ri + rj:
                            continue
                        n_overlap += 1
                        near, far = (i, j) if depth[i] > depth[j] else (j, i)
                        if pos_of[near] < pos_of[far]:
                            fail(
                                f"{smi}/{style_name}: 原子 {near} 比 {far} 靠前 "
                                f"(深度差 {gap:.4f} Å),却画在了前面"
                            )

                # ⑪ placed 与图里的圆一致。SVG 里的坐标是 `{:.2}` 写出去的,
                # 所以这一条只能判到**半个百分位**(0.005 磅)—— 判不了更细的
                # 差别是这条路的固有上限,不是容差放松:图里就只有那么多位。
                placed = np.array(st["placed"], dtype=float)
                for atom, slot in pos_of.items():
                    cx, cy, r, _, _ = circles[slot]
                    if (abs(placed[atom][0] - cx) > 0.006
                            or abs(placed[atom][1] - cy) > 0.006
                            or abs(placed[atom][2] - r) > 0.006):
                        fail(
                            f"{smi}/{style_name}: placed[{atom}] = "
                            f"({placed[atom][0]:.3f},{placed[atom][1]:.3f},"
                            f"r={placed[atom][2]:.3f}) 与图里的圆 ({cx},{cy},r={r}) 对不上"
                        )

            # ⑦ 键的两半各随一端的颜色
            if st["stick_radius_a"] > 0 and lines:
                base = canvas_positions(proj, scale, circles, lines)
                if base is not None:
                    n_ok, n_amb = check_half_colours(
                        base, bonds, z, lines, palette, st, scale, fail, smi, style_name
                    )
                    n_halves += n_ok
                    n_ambiguous_half += n_amb

                    # ⑧ 并排方向 ⊥ 投影(只在没有球的样式上判,几何是闭式的)
                    if st["ball_vdw_frac"] == 0 and st["spacing_a"] > 0:
                        n_multiple += check_multiple(
                            base, bonds, lines, st, scale, fail, smi, style_name
                        )

            # ⑩ 不出画布
            for c in circles:
                if (c[0] - c[2] < -0.01 or c[0] + c[2] > st["width"] + 0.01
                        or c[1] - c[2] < -0.01 or c[1] + c[2] > st["height"] + 0.01):
                    fail(f"{smi}/{style_name}: 圆 ({c[0]:.1f},{c[1]:.1f})±{c[2]:.1f} 出了画布")
            for ln in lines:
                hw = ln[5] / 2
                for (x, y) in ((ln[0], ln[1]), (ln[2], ln[3])):
                    if (x - hw < -0.01 or x + hw > st["width"] + 0.01
                            or y - hw < -0.01 or y + hw > st["height"] + 0.01):
                        fail(f"{smi}/{style_name}: 线端 ({x:.1f},{y:.1f}) 出了画布")

            # ⑨ 重编号逐字节相同
            if svg != st["svg_renumbered"]:
                fail(f"{smi}/{style_name}: 把原子倒序重编号之后画出来不一样了")

    # ---- 分母闸 ----
    floors = [
        ("分子", n_mol, MIN_MOLECULES),
        ("圆(球)", n_circles, MIN_CIRCLES),
        ("投影上重叠的球对", n_overlap, MIN_OVERLAP_PAIRS),
        ("认得出归属的半根键", n_halves, MIN_STICK_HALVES),
        ("并排画的多重键", n_multiple, MIN_MULTIPLE_BONDS),
        ("见到的元素种类", len(elements), MIN_ELEMENTS),
        ("拿去与外部主轴比的视角", n_axes, MIN_AXES_COMPARED),
    ]
    empty = [f"{name} 只有 {got},下限 {want}" for name, got, want in floors if got < want]

    print(f"分子 {n_mol}、样式 {len(styles_seen)} 档:{'、'.join(sorted(styles_seen))}")
    print(f"  圆 {n_circles}、重叠球对 {n_overlap}、半根键 {n_halves}"
          f"(认不清归属而跳过 {n_ambiguous_half})、并排多重键 {n_multiple}")
    print(f"  深度差在 {DEPTH_TIE} Å 以内、不判次序的球对 {n_depth_tie};"
          f"投影上与别人重合、只比多重集的原子 {n_coincident}")
    print(f"  与外部主轴逐根比过 {n_axes} 个视角;主轴简并而跳过 {n_degenerate} 个")
    print(f"  元素 {len(elements)} 种")

    if empty:
        print("判据没东西可判:", file=sys.stderr)
        for e in empty:
            print(f"  {e}", file=sys.stderr)
        return 1
    if bad:
        print(f"分歧 {len(bad)} 条(最多列 30):", file=sys.stderr)
        for b in bad:
            print(f"  {b}", file=sys.stderr)
        return 1
    print("全部通过")
    return 0


def match_circles(obs: np.ndarray, want_rel: np.ndarray):
    """把图里的圆认回原子。返回按位置分的组 `[(圆下标们, 原子下标们), ...]`。

    画布的平移量判官不知道(它由包围盒定),所以先按中位数解出来。

    # 为什么要分组,而不是一一对应

    **两个原子在投影上完全重合是会发生的**:分子有一张镜面、而主轴恰好把那张
    镜面摆成屏幕平面时,镜面两侧的原子投影到同一点。实测 400 个分子里有 9 个
    是这样(三维坐标一个都不重合,是投影重合)。硬要一一对应就会认不回来,
    而那不是实现的毛病。

    重合的一组里,判官只比**半径与颜色的多重集**,不判组内谁先画 —— 那时
    "谁盖着谁"在图上本来就看不出来。组数与跳过的个数会打印出来。
    """
    shift = np.median(obs, axis=0) - np.median(want_rel, axis=0)
    want = want_rel + shift

    # 按**距离**聚簇,不按坐标四舍五入到固定格点。
    #
    # 格点法有个静默的坑,是撞出来的:两个氢的投影相距 0.0018 磅 —— 图里根本
    # 分不开(SVG 只写两位小数),而它们恰好跨在格点边界的两侧,于是被分进两组,
    # 一组分到两个圆、另一组一个都没有,判官报"认不回原子"。
    #
    # 门槛取 TOL_PT:比这更近的两个原子,在两位小数的 SVG 里本来就是同一个点。
    cluster = list(range(len(want)))

    def root(x: int) -> int:
        while cluster[x] != x:
            cluster[x] = cluster[cluster[x]]
            x = cluster[x]
        return x

    for a in range(len(want)):
        for b in range(a):
            if float(np.hypot(*(want[a] - want[b]))) <= TOL_PT:
                cluster[root(a)] = root(b)

    atom_groups: dict = collections.defaultdict(list)
    for i, _ in enumerate(want):
        atom_groups[root(i)].append(i)
    centres = {k: want[v].mean(axis=0) for k, v in atom_groups.items()}

    circle_groups: dict = collections.defaultdict(list)
    for j, o in enumerate(obs):
        best, bd = None, None
        for k, c in centres.items():
            d = float(np.hypot(*(c - o)))
            if bd is None or d < bd:
                best, bd = k, d
        if bd is None or bd > TOL_PT:
            return None
        circle_groups[best].append(j)

    if set(atom_groups) != set(circle_groups):
        return None
    out = []
    for k, atoms in atom_groups.items():
        circles = circle_groups[k]
        if len(circles) != len(atoms):
            return None
        out.append((circles, atoms))
    return out


def canvas_positions(proj, scale, circles, lines) -> np.ndarray | None:
    """每个原子在画布上的位置。

    画布的平移量由包围盒定,判官不知道 —— 只能从图里解出来。

    有球的样式按圆心的中位数解。没球的样式(棍状、线框)先前按"线端的最小值
    对齐原子的最小值"解,**那是错的**:多重键并排的圆柱偏到了原子中心之外,
    包围盒的边界是那几根偏出去的圆柱,不是最边上的原子。差出来的那半个间距
    让每一根多重键都报"图里没有",而实现一点毛病没有。

    改成投票:每个"线端 − 原子"的差都是一个候选平移量,取被最多原子认可的
    那个。单键的端点正落在原子中心上,所以正确的平移量总能拿到多数票,而
    偏出去的圆柱各自散开,凑不出第二个峰。
    """
    rel = np.stack([proj[:, 0] * scale, -proj[:, 1] * scale], axis=1)
    if circles:
        shift = np.median(np.array([[c[0], c[1]] for c in circles]), axis=0) - np.median(rel, axis=0)
        return rel + shift
    if not lines:
        return None
    ends = np.array([[ln[0], ln[1]] for ln in lines] + [[ln[2], ln[3]] for ln in lines])
    votes = collections.Counter()
    for e in ends:
        for r in rel:
            votes[(round(float(e[0] - r[0]), 2), round(float(e[1] - r[1]), 2))] += 1
    (dx, dy), n = votes.most_common(1)[0]
    if n < 2:
        return None
    return rel + np.array([dx, dy])


def check_half_colours(base, bonds, z, lines, palette, st, scale, fail, smi, style_name):
    """每一段圆柱的颜色 = 它所属那半根键、那一端的元素色。

    # 半根键的位置要**连并排的偏移一起**算出来

    先前这里只按"原子中心 → 键中点"那条不偏移的线去认,于是多重键那两根偏出去
    的圆柱谁都认不上,再被最近的另一根键**唯一**认走 —— 报出来的是颜色错了,
    而实现一点毛病没有。16 条假分歧全是这么来的。

    一段同时贴着好几半键(投影上重合)时**跳过**并记一笔:含糊的归属只会造出
    假分歧。跳过的段数会打印出来。
    """
    halves = []
    for a, b, order in bonds:
        pa, pb = base[a], base[b]
        d = pb - pa
        flat = float(np.hypot(*d))
        k = CYLINDERS.get(order, 1)
        if st["spacing_a"] <= 0 or flat < 1e-9:
            k = 1
        perp = np.array([d[1], -d[0]]) / flat if flat > 1e-9 else np.zeros(2)
        for idx in range(k):
            off = (idx - (k - 1) / 2) * st["spacing_a"] * scale
            qa, qb = pa + perp * off, pb + perp * off
            mid = (qa + qb) / 2
            halves.append((qa, mid, palette.get(int(z[a]), UNKNOWN_COLOR)))
            halves.append((qb, mid, palette.get(int(z[b]), UNKNOWN_COLOR)))

    ok = amb = 0
    for ln in lines:
        m = np.array([(ln[0] + ln[2]) / 2, (ln[1] + ln[3]) / 2])
        hits = [h for h in halves if point_on_segment(m, h[0], h[1])]
        if len(hits) != 1:
            amb += 1
            continue
        ok += 1
        if ln[4] != hits[0][2]:
            fail(
                f"{smi}/{style_name}: 半根键画成了 "
                f"#{ln[4][0]:02x}{ln[4][1]:02x}{ln[4][2]:02x},该随那一端的元素色 "
                f"#{hits[0][2][0]:02x}{hits[0][2][1]:02x}{hits[0][2][2]:02x}"
            )
    return ok, amb


def point_on_segment(p, a, b, tol=0.3) -> bool:
    d = b - a
    L2 = float(d @ d)
    if L2 < 1e-12:
        return bool(np.hypot(*(p - a)) < tol)
    t = float((p - a) @ d) / L2
    if t < -0.05 or t > 1.05:
        return False
    return bool(np.hypot(*(p - (a + t * d))) < tol)


def check_multiple(base, bonds, lines, st, scale, fail, smi, style_name):
    """多重键并排的两(三)根圆柱,偏移方向必须垂直于键在屏幕上的投影。

    只在没有球的样式上判:那一档不用截,起点就是原子中心,期望位置是闭式的。
    """
    n = 0
    ends = np.array([[ln[0], ln[1]] for ln in lines])
    for a, b, order in bonds:
        k = CYLINDERS.get(order)
        if k is None:
            continue
        pa, pb = base[a], base[b]
        d = pb - pa
        flat = float(np.hypot(*d))
        if flat < 1e-9:
            continue          # 键正对观察者,本来就只画一根
        perp = np.array([d[1], -d[0]]) / flat
        n += 1
        for idx in range(k):
            off = (idx - (k - 1) / 2) * st["spacing_a"] * scale
            want = pa + perp * off
            if float(np.min(np.hypot(ends[:, 0] - want[0], ends[:, 1] - want[1]))) > TOL_PT:
                fail(
                    f"{smi}/{style_name}: {order} 键 {a}-{b} 的第 {idx} 根圆柱该从 "
                    f"({want[0]:.2f},{want[1]:.2f}) 起步,图里没有"
                )
    return n


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("jsonl", type=pathlib.Path)
    args = ap.parse_args()
    return check(args.jsonl, pathlib.Path(__file__).resolve().parent.parent)


if __name__ == "__main__":
    raise SystemExit(main())
