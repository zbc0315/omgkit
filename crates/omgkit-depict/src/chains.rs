//! 链与取代基:给一个原子周围**还没放置**的邻居分配方向。
//!
//! # 规则
//!
//! 理想夹角按度数定:度数 ≤ 3 用 [`Style::chain_angle_deg`](crate::style::Style)
//! (两套内置规范都是 120°),度数 ≥ 4 用 360/度数。
//!
//! 度数 2 的原子**必须**取 120° 而不是 180° —— 否则链会拉成直线,那既不是化学
//! 惯例,也让后面的可旋转键无处可翻。
//!
//! 新邻居放进**最大的空闲扇区**并在其中均分。这条规则在各种度数下都退化得对:
//!
//! | 已占方向 | 空闲扇区 | 做法 |
//! |---|---|---|
//! | 0 个(起点) | 360° | 从一个固定角度起,按理想夹角铺开 |
//! | 1 个(链上一步) | 360° | 取 ±理想夹角,符号由锯齿决定 |
//! | ≥2 个 | < 360° | 在最大空隙内均分 |
//!
//! 中间那一档不能套用"扇区居中":360° 扇区的中心正好是已占方向的反向,那就把
//! 链拉成了直线。
//!
//! # 锯齿的符号必须可复现
//!
//! 直链每一步左右交替才成锯齿。交替本身要有个起点,而那个起点若取自原子的**存储
//! 下标**,同一个分子换种 SMILES 写法就会得到镜像的图。这里由父原子的符号翻转
//! 得到,起点则来自规范秩 —— 写法就影响不到。

use std::collections::BTreeMap;

use omgkit_core::{BondOrder, MolBuilder};

use crate::geom::{segments_cross, Point2, BOND_LEN};
use crate::style::Style;

/// 起点原子第一根键的方向。
///
/// 取 30° 是惯例:六元环按这个角度画出来是"尖朝上"的标准姿态,直链则呈水平
/// 锯齿。取 0° 会让直链变成一条水平线加上下交替的折点,观感上偏斜。
const SEED_ANGLE: f64 = std::f64::consts::FRAC_PI_6;

/// 一次放置的结果。
pub(crate) struct Placed {
    /// 原子
    pub atom: u32,
    /// 坐标
    pub at: Point2,
    /// 它自己的锯齿符号,传给它的子代取反
    pub zig: i8,
}

/// 给 `a` 周围还没放置的邻居 `todo` 分配坐标。
///
/// `todo` 必须**已按规范秩排好**;顺序决定谁分到哪个方向,拿存储下标排就会
/// 引入写法依赖。
pub(crate) fn place_neighbours(
    mol: &MolBuilder,
    a: u32,
    pos: &BTreeMap<u32, Point2>,
    todo: &[u32],
    ranks: &[u32],
    style: &Style,
    zig: i8,
) -> Vec<Placed> {
    if todo.is_empty() {
        return Vec::new();
    }
    let center = pos[&a];

    // 已占方向:所有已经放好的邻居。
    //
    // **排序要量化后再按规范秩打破平局。** 两个方向在数学上相等、浮点上差
    // 1e-16 时,直接比大小会让它们的先后取决于算到那一步的运算次序 —— 而
    // `allocate` 挑最大空隙时用的正是这个次序,于是同一个分子换个写法,某个
    // 取代基就换了个方向挂。实测:一个稠三环的甲基因此差了 120°。
    const QUANT: f64 = 1e9;
    let mut occ: Vec<((i64, u32), f64)> = mol
        .neighbors(a)
        .filter_map(|(n, _)| {
            pos.get(&n).map(|p| {
                let t = (*p - center).angle();
                #[allow(clippy::cast_possible_truncation)]
                (((t * QUANT).round() as i64, ranks[n as usize]), t)
            })
        })
        .collect();
    occ.sort_unstable_by_key(|x| x.0);
    let occupied: Vec<f64> = occ.iter().map(|x| x.1).collect();

    let ideal = ideal_angle(mol, a, style);
    let dirs = allocate(&occupied, todo.len(), ideal, zig);

    debug_assert_eq!(dirs.len(), todo.len(), "方向数必须与待放邻居数相等");

    // 已经占住的位置,以及已经画出来的键。新原子不许落在前者上、新键不许与
    // 后者交叉 —— 见 [`free_direction`]。
    let mut taken: Vec<Point2> = pos.values().copied().collect();
    let mut drawn: Vec<(Point2, Point2)> = mol
        .bonds()
        .iter()
        .filter_map(|b| Some((*pos.get(&b.begin)?, *pos.get(&b.end)?)))
        .collect();
    let mut out = Vec::with_capacity(todo.len());
    for (&atom, theta) in todo.iter().zip(dirs) {
        let theta = free_direction(center, theta, &taken, &drawn);
        let at = center + Point2::new(BOND_LEN, 0.0).rotated(theta);
        taken.push(at);
        drawn.push((center, at));
        out.push(Placed {
            atom,
            at,
            // 子代取反,直链就走出锯齿
            zig: -zig,
        });
    }
    out
}

/// 从 `ideal` 出发,找一个不会与已放好的原子重合的方向。
///
/// # 为什么宁可歪着也不重合
///
/// 两个原子叠在同一点上时,它们各自的键首尾相接 —— **图上就多出一个分子里
/// 没有的环**,而读者没有任何办法看出那个环是假的。角度偏离理想值只是难看,
/// 不会让人读错结构。
///
/// 实测:一个三萜的两个甲基落在同一个栅格点上,图上凭空出现一个三元环,三条
/// 边正好都是一个键长。全语料上这种重合占 6%,而且距离全是**正好 0**:布局
/// 走的是 30° 栅格上的单位步长,两条支路撞到同一个格点是系统性的,不是浮点抖动。
///
/// 按 30° 一档往两边试,与整张图的栅格一致;五档之内都腾不开就退回 `ideal`,
/// 交给消冲突,消不掉再如实报进 `unresolved`。
fn free_direction(center: Point2, ideal: f64, taken: &[Point2], drawn: &[(Point2, Point2)]) -> f64 {
    const STEP: f64 = std::f64::consts::FRAC_PI_6;
    /// 多近算重合。取键长的十分之一 —— 真正分得开的两个位置至少差半个键长。
    const TOL: f64 = 0.1;
    let at = |t: f64| center + Point2::new(BOND_LEN, 0.0).rotated(t);
    let clear = |t: f64| {
        let p = at(t);
        !taken.iter().any(|q| p.dist(*q) < TOL)
    };
    // 新键与已画的键交叉。共端点不算 —— 那是相邻的键,`segments_cross` 已经放过。
    let uncrossed = |t: f64| {
        let p = at(t);
        !drawn.iter().any(|(u, v)| segments_cross(center, p, *u, *v))
    };

    // 候选:理想方向,然后按 30° 一档往两边铺开
    let mut cands = Vec::with_capacity(11);
    cands.push(ideal);
    for k in 1..=5 {
        cands.push(ideal + STEP * f64::from(k));
        cands.push(ideal - STEP * f64::from(k));
    }

    // 两轮:先要"既不重合也不交叉",都腾不开就退而只求"不重合"。
    // **重合排在交叉前面** —— 重合会凭空造出一个假环,交叉只是难读。
    if let Some(t) = cands.iter().find(|t| clear(**t) && uncrossed(**t)) {
        return *t;
    }
    cands.iter().copied().find(|t| clear(*t)).unwrap_or(ideal)
}

/// 一个原子周围相邻两根键的理想夹角(弧度)。
///
/// # sp 的原子要画成直线
///
/// 只看度数是不够的:氰基的碳、炔碳、累积双键的中心碳都是 **sp 杂化,键角
/// 180°**,而它们的度数是 2 —— 按度数给 120° 的话,`R—C≡N` 会画成折的。
/// 这不是好看不好看的问题,是**画错了**:读者从图上读到的键角与分子的实际
/// 几何不符,而线条本身看着一点毛病没有。
///
/// 判据是键级不是度数:有三键,或者有两根双键(累积双键),就是 sp。
fn ideal_angle(mol: &MolBuilder, a: u32, style: &Style) -> f64 {
    let mut doubles = 0usize;
    let mut triple = false;
    for (_, bi) in mol.neighbors(a) {
        match mol.bonds()[bi as usize].order {
            BondOrder::Triple => triple = true,
            BondOrder::Double => doubles += 1,
            _ => {}
        }
    }
    if triple || doubles >= 2 {
        return std::f64::consts::PI;
    }
    let degree = mol.degree(a);
    if degree <= 3 {
        style.chain_angle_deg.to_radians()
    } else {
        std::f64::consts::TAU / degree as f64
    }
}

/// 把 `n` 个新方向分配到已占方向 `occupied`(已排序)之外的空隙里。
fn allocate(occupied: &[f64], n: usize, ideal: f64, zig: i8) -> Vec<f64> {
    let sign = if zig >= 0 { 1.0 } else { -1.0 };

    match occupied.len() {
        // 起点:从固定角度铺开
        0 => (0..n).map(|k| SEED_ANGLE + ideal * k as f64).collect(),

        // 链上一步:**不能取扇区中心**,那是 180°,会把链拉直
        1 => {
            let base = occupied[0];
            (0..n)
                .map(|k| {
                    // 第一个取 ±ideal,其余向另一侧交替铺开
                    let step = (k as f64 / 2.0).floor() + 1.0;
                    let s = if k % 2 == 0 { sign } else { -sign };
                    base + s * ideal * step
                })
                .collect()
        }

        // 已有两个以上方向:找最大空隙,在里面均分
        _ => {
            let (start, gap) = largest_gap(occupied);
            // n 个新方向把空隙切成 n+1 份
            (0..n)
                .map(|k| start + gap * (k as f64 + 1.0) / (n as f64 + 1.0))
                .collect()
        }
    }
}

/// 已排序角度序列中最大的空隙:返回(空隙起始角, 空隙大小)。
fn largest_gap(sorted: &[f64]) -> (f64, f64) {
    let n = sorted.len();
    debug_assert!(n >= 2);
    let mut best = (
        sorted[n - 1],
        sorted[0] + std::f64::consts::TAU - sorted[n - 1],
    );
    for i in 0..n - 1 {
        let g = sorted[i + 1] - sorted[i];
        if g > best.1 {
            best = (sorted[i], g);
        }
    }
    best
}

#[cfg(test)]
mod tests {
    #[test]
    fn no_two_atoms_are_drawn_on_the_same_point() {
        // **重合是最糟的一种画错。** 两个原子叠在一起时,它们各自的键首尾相接,
        // 图上就多出一个分子里没有的环 —— 读者没有任何办法看出那个环是假的。
        //
        // 实测:下面第一个三萜的两个甲基落在同一个栅格点上,图上凭空出现一个
        // 三元环,三条边正好都是一个键长。全语料上这种重合曾占 6%,而且距离
        // 全是**正好 0** —— 布局走的是 30° 栅格上的单位步长,两条支路撞到同
        // 一个格点是系统性的,不是浮点抖动。
        for smi in [
            "CC([CH]1CC[C]2(CC[C]3(C)[C]4(C)[CH](CC[CH]3[CH]12)[C]1(C)[CH](CC4)C([CH](CC1)O)(C)C)CO)=C",
            "[O-][N+](=O)C1=CC(=CC=C1Cl)S(=O)(=O)C2=CC=C(Cl)C(=C2)[N+]([O-])=O",
            "CC(C)(C)c1ccccc1",
            "CC(=O)Oc1ccccc1C(=O)O",
        ] {
            for style in &Style::ALL {
                let mut m = omgkit_io::smiles::parse(smi).unwrap();
                omgkit_chem::pipeline::sanitize(&mut m).unwrap();
                let d = crate::generate(&m, style);
                for i in 0..d.coords.len() {
                    for j in (i + 1)..d.coords.len() {
                        let dist = d.coords[i].dist(d.coords[j]);
                        assert!(
                            dist > 0.1,
                            "[{}] {smi}:原子 {i} 与 {j} 相距 {dist:.4} 个键长 —— 画在同一点上了",
                            style.name
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn an_sp_atom_is_drawn_straight() {
        // 氰基的碳、炔碳、累积双键的中心碳都是 sp 杂化,键角 **180°**。它们的
        // 度数是 2,只按度数给 120° 的话 `R—C≡N` 会画成折的 —— 那是画错了,
        // 而线条本身看不出毛病。
        for (smi, centre) in [
            ("CC#N", 1u32), // 乙腈:C1 是 sp
            ("CC#CC", 1),   // 2-丁炔
            ("CC=C=CC", 2), // 累积双键的中心碳
            ("N#CC(C)(C)C#N", 1),
        ] {
            let mut m = omgkit_io::smiles::parse(smi).unwrap();
            omgkit_chem::pipeline::sanitize(&mut m).unwrap();
            let d = crate::generate(&m, &Style::ACS_1996);
            let nbrs: Vec<u32> = m.neighbors(centre).map(|(n, _)| n).collect();
            assert!(nbrs.len() >= 2, "{smi}:中心该有两个邻居");
            let (p, q) = (d.coords[nbrs[0] as usize], d.coords[nbrs[1] as usize]);
            let c = d.coords[centre as usize];
            let (u, v) = ((p - c).normalized(), (q - c).normalized());
            let deg = u.dot(v).clamp(-1.0, 1.0).acos().to_degrees();
            assert!(
                (deg - 180.0).abs() < 1e-6,
                "{smi}:原子 {centre} 是 sp,键角却画成了 {deg:.1}°"
            );
        }
    }

    /// 判据里也要按规范秩打破平局 —— 与实现取同一个来源
    fn canonical(m: &MolBuilder) -> Vec<u32> {
        omgkit_io::canon::canonical_ranks(m)
    }

    use super::*;
    use crate::style::Style;

    const TOL: f64 = 1e-9;

    fn prep(smi: &str) -> MolBuilder {
        let mut m = omgkit_io::smiles::parse(smi).unwrap();
        omgkit_chem::pipeline::sanitize(&mut m).unwrap();
        m
    }

    /// 两个方向之间的夹角,取 [0, π]
    fn between(u: Point2, v: Point2) -> f64 {
        let c = u.normalized().dot(v.normalized()).clamp(-1.0, 1.0);
        c.acos()
    }

    #[test]
    fn a_chain_zigzags_instead_of_running_straight() {
        // 度数 2 的原子若取 180°,链就成了一条直线 —— 既不是化学惯例,
        // 也让后面的可旋转键无处可翻。这条守的正是那个 180°。
        let m = prep("CCCCC");
        let style = Style::ACS_1996;
        let mut pos: BTreeMap<u32, Point2> = BTreeMap::new();
        pos.insert(0, Point2::ORIGIN);
        let mut zig = 1i8;
        for a in 0..4u32 {
            let out = place_neighbours(&m, a, &pos, &[a + 1], &canonical(&m), &style, zig);
            pos.insert(out[0].atom, out[0].at);
            zig = out[0].zig;
        }
        assert_eq!(pos.len(), 5);
        for i in 1..4u32 {
            let ang = between(pos[&(i - 1)] - pos[&i], pos[&(i + 1)] - pos[&i]);
            assert!(
                (ang - 120f64.to_radians()).abs() < TOL,
                "第 {i} 个原子处的键角是 {:.1}°,应当是 120°",
                ang.to_degrees()
            );
        }
        // 锯齿:相邻两步的转向必须相反,否则会绕成圆弧
        let turn = |i: u32| (pos[&i] - pos[&(i - 1)]).cross(pos[&(i + 1)] - pos[&i]);
        assert!(turn(1) * turn(2) < 0.0, "第 1、2 步没有交替转向");
        assert!(turn(2) * turn(3) < 0.0, "第 2、3 步没有交替转向");
    }

    #[test]
    fn every_bond_is_one_unit_long() {
        let m = prep("CC(C)(C)C");
        let style = Style::ACS_1996;
        let mut pos: BTreeMap<u32, Point2> = BTreeMap::new();
        pos.insert(1, Point2::ORIGIN);
        let mut todo: Vec<u32> = m.neighbors(1).map(|(n, _)| n).collect();
        todo.sort_unstable();
        for p in place_neighbours(&m, 1, &pos, &todo, &canonical(&m), &style, 1) {
            pos.insert(p.atom, p.at);
        }
        for n in todo {
            let d = pos[&n].dist(pos[&1]);
            assert!((d - BOND_LEN).abs() < TOL, "键长 {d}");
        }
    }

    #[test]
    fn four_substituents_are_spread_not_stacked() {
        // 季碳:四根键必须分开。理想夹角在度数 4 时该退到 90°,若仍用 120°,
        // 四个方向只能铺满 360° 中的 360°—— 会有两个重叠。
        let m = prep("CC(C)(C)C");
        let style = Style::ACS_1996;
        let mut pos: BTreeMap<u32, Point2> = BTreeMap::new();
        pos.insert(1, Point2::ORIGIN);
        let mut todo: Vec<u32> = m.neighbors(1).map(|(n, _)| n).collect();
        todo.sort_unstable();
        assert_eq!(todo.len(), 4, "季碳应当有四个邻居");
        let out = place_neighbours(&m, 1, &pos, &todo, &canonical(&m), &style, 1);
        for i in 0..out.len() {
            for j in (i + 1)..out.len() {
                let ang = between(out[i].at, out[j].at);
                assert!(
                    ang > 45f64.to_radians(),
                    "第 {i}、{j} 个取代基只差 {:.1}°,挤在一起了",
                    ang.to_degrees()
                );
            }
        }
    }

    #[test]
    fn a_new_branch_goes_into_the_largest_free_sector() {
        // 已经占了两个方向时,新的必须落进**最大的空隙**。落进小空隙不会报错,
        // 只会让图上一边挤一边空。
        let m = prep("CC(C)C");
        let style = Style::ACS_1996;
        let mut pos: BTreeMap<u32, Point2> = BTreeMap::new();
        pos.insert(1, Point2::ORIGIN);
        // 手工把两个邻居摆在 0° 和 60°,留下一个 300° 的大空隙
        pos.insert(0, Point2::new(1.0, 0.0));
        pos.insert(2, Point2::new(0.5, 3f64.sqrt() / 2.0));
        let out = place_neighbours(&m, 1, &pos, &[3], &canonical(&m), &style, 1);
        let ang = out[0].at.angle().rem_euclid(std::f64::consts::TAU);
        // 大空隙是 60° → 360°,中点在 210°
        assert!(
            (ang - 210f64.to_radians()).abs() < 1e-6,
            "新支链落在 {:.1}°,应当落在最大空隙的中点 210°",
            ang.to_degrees()
        );
    }

    #[test]
    fn the_largest_gap_wraps_around_the_seam() {
        // 空隙搜索必须绕过 ±π 的接缝。漏掉环绕的那一段不会报错,只会在某些
        // 角度组合下把新键塞进一个其实很窄的缝里。
        //
        // **数据必须让环绕的那一段真的是最大空隙。** 第一版用的是
        // `[-3.0, -2.9, 3.0]`,那里最大的其实是 -2.9→3.0 这个**内部**空隙,
        // 于是把环绕逻辑破坏掉,这条照样绿 —— 名字说的和断言测的是两回事。
        // 三个方向挤在 0 附近,环绕那一段才是最大的。
        let sorted = vec![-0.1, 0.0, 0.1];
        let (start, gap) = largest_gap(&sorted);
        assert!(
            (start - 0.1).abs() < TOL,
            "空隙应当从最后一个方向 0.1 起,实得 {start}"
        );
        let want = -0.1 + std::f64::consts::TAU - 0.1;
        assert!((gap - want).abs() < TOL, "空隙大小应当是 {want},实得 {gap}");

        // 顺带守住"内部空隙也要能选中",免得改成只看环绕
        let inner = vec![-3.0, -2.9, 3.0];
        let (s2, g2) = largest_gap(&inner);
        assert!(
            (s2 - (-2.9)).abs() < TOL && (g2 - 5.9).abs() < TOL,
            "内部空隙没选对"
        );
    }
}
