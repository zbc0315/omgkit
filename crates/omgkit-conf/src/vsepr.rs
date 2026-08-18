//! 一个中心原子的**取代基怎么排布** —— 只管排布,不管角度大小。
//!
//! # 与 v1 的分工不同,这一条是记着教训的
//!
//! v1 拿**杂化去推角度**:二价氧被感知判成 `Sp2` → steric 3 → 三角平面 → 120°,
//! 而实测真值跨 104.0~117.4°(醇 106.7、酯氧 115.6、芳醚氧 117.0)——
//! **4925 个角全错**,占当时全部键角违例的 26.3%。
//!
//! 这一版**角度一律来自实测表**([`crate::params::angle`]),杂化只用来回答
//! 另一个问题:**几个取代基之间怎么摆** —— 共平面(sp²)还是交错(sp³)。
//! 那才是杂化真正说的事,而且这件事表里查不到(表只给一个角,不给排布)。

use omgkit_core::Hybridization;

/// 取代基绕"父键"的排布方式。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Arrangement {
    /// 直线(sp):只有一个取代基,与父键成 180°,扭转角无意义。
    Linear,
    /// 共平面(sp²/芳香):中心与它的取代基**必须共面**,
    /// 所以取代基的扭转角只能取 0° 或 180°。
    Planar,
    /// 交错(sp³):取代基绕父键均分,反式那个在 180°。
    Tetrahedral,
    /// 别的(配位数 ≥ 5 等):均分,只求不重叠。
    Spread,
}

/// 由杂化与配位数定排布方式。
///
/// **配位数 ≥ 5 的中心走 [`Arrangement::Spread`]** —— 那是"均分了事",
/// 不是真的三角双锥/八面体。语料里 65 个分子(0.74%)有这种中心,
/// 一期**明确不保证**它们摆得对,判据要单独计数(见方案 §4.5)。
#[must_use]
pub fn arrangement(hyb: Hybridization, degree: usize) -> Arrangement {
    if degree >= 5 {
        return Arrangement::Spread;
    }
    match hyb {
        Hybridization::Sp => Arrangement::Linear,
        Hybridization::Sp2 => Arrangement::Planar,
        Hybridization::Sp3 => Arrangement::Tetrahedral,
        // **sp²d 是平面四方形(90°/180°),这个构造法给不出来。**
        // 头一版把它映到四面体 —— 方向是错的:四面体给 109.47°,
        // 平面四方形一个 109.47° 都没有。语料里现在没有 sp²d 的中心,
        // 但映错了迟早咬人,所以走 `Spread`("均分了事",而且一期明确不保证它)。
        Hybridization::Sp2d => Arrangement::Spread,
        // 感知没给的:按配位数猜,2 当直线、3 当平面、其余当四面体
        _ => match degree {
            0..=2 => Arrangement::Linear,
            3 => Arrangement::Planar,
            _ => Arrangement::Tetrahedral,
        },
    }
}

/// **要摆 `n` 个取代基时,它们的扭转角各取多少**(弧度)。
///
/// # 返回的个数可能**少于** `n`,那是故意的
///
/// 一种排布能放下的取代基数是有上限的:直线 1 个、共平面 2 个、四面体 3 个
/// (都是在父键之外算的)。要多了就**只给得下的那几个**,剩下的由调用方计数。
///
/// 头一版是硬凑够 `n` 个,于是 `Planar` 要 3 个时给出 `[180°, 0°, 180°]` ——
/// 第 0 个和第 2 个拿到**同样的键长、同样的键角、同样的扭转角**,NeRF 摆出
/// **同一个坐标**。实测 `[Zn](C)(C)(C)C`(Zn 被感知成 Sp2、配位 4)两个甲基碳
/// 相距 **0.000000 Å**,而 `place()` 报的是 `complete = true`、`degenerate = 0` ——
/// 两个原子重合对 `1/r¹²` 就是无穷大,后续优化器直接起不来,
/// 而这正是本 crate 文档列的"力场救不回来的四条"之一。
///
/// 扭转角是 `祖父–父–中心–取代基` 那个二面角。约定:
///
/// - [`Arrangement::Linear`]:只摆一个,扭转角无意义(角是 180°,
///   NeRF 里 `sin(角)=0`,扭转项自动消掉),给 π 占位;
/// - [`Arrangement::Planar`]:**只能是 0 或 π** —— 别的值会让 sp² 中心离面。
///   所以至多给 2 个:`[π, 0]`(先反式,后顺式);
/// - [`Arrangement::Tetrahedral`]:从 π 起按 `2π/3` 均分 → `[π, π/3, 5π/3]`,
///   这就是标准的**交错**构象(至多 3 个);
/// - [`Arrangement::Spread`]:`n+1` 个方向均分(把父键也算一个)。
///
/// **反式(π)永远排第一** —— 链要伸展开,而伸展的那一支应当优先。
#[must_use]
pub fn child_torsions(arr: Arrangement, n: usize) -> Vec<f64> {
    use std::f64::consts::{PI, TAU};
    match arr {
        Arrangement::Linear => vec![PI; n.min(1)],
        // 共平面中心除父键外只放得下 2 个 —— 要第 3 个就得离面,那不再是平面中心
        Arrangement::Planar => (0..n.min(2))
            .map(|k| if k % 2 == 0 { PI } else { 0.0 })
            .collect(),
        // 四面体除父键外放得下 3 个;第 4 个会与第 1 个重合(π + 3·2π/3 = π)
        Arrangement::Tetrahedral => (0..n.min(3))
            .map(|k| {
                #[allow(clippy::cast_precision_loss)]
                let t = PI + (k as f64) * TAU / 3.0;
                t.rem_euclid(TAU)
            })
            .collect(),
        Arrangement::Spread => (0..n)
            .map(|k| {
                #[allow(clippy::cast_precision_loss)]
                let t = PI + (k as f64) * TAU / ((n + 1) as f64);
                t.rem_euclid(TAU)
            })
            .collect(),
    }
}

/// **按排布算:表角 `θ` 之下,兄弟角与 `θ` 应当差多少**(弧度)。
///
/// 判据要用它 —— 拿兄弟角直接去比表值是错的口径。两种排布各有各的式子:
///
/// | 排布 | 兄弟角 | 与 θ 之差 | 何时为 0 |
/// |---|---|---|---|
/// | `Planar` | `min(2θ, 2π − 2θ)` | 见下 | θ = 120° |
/// | `Tetrahedral` | `arccos(cos²θ + sin²θ·cos120°)` | 见下 | θ = 109.4712° |
/// | `Linear` | 没有兄弟 | 0 | 恒 |
///
/// **平面那一行被改过两次,两次都是判据逼的**:
///
/// 1. 头一版写的边界是 0,结果氮那个中心(表值 115.60°)实得 128.80°、
///    超出 13.2° 被判红 —— 而 `360 − 2×115.6 = 128.8`,几何对、**公式错**;
/// 2. 改成 `|2π − 3θ|` 之后仍然只在 **θ ≥ 90°** 时成立。两个取代基都在平面内、
///    各与父键成 θ、分居两侧,夹角是 `2θ`;只有 `2θ > π` 时那个夹角才折回成
///    `2π − 2θ`。θ = 60° 时真值是 60°,而 `|2π − 3θ| = 180°` ——
///    **方向是只把判据变绿**,那种错最难发现。
#[must_use]
pub fn expected_sibling_skew(arr: Arrangement, theta: f64) -> f64 {
    use std::f64::consts::TAU;
    match arr {
        // 夹角是 2θ,超过 π 才折回成 2π − 2θ
        Arrangement::Planar => ((TAU - 2.0 * theta).min(2.0 * theta) - theta).abs(),
        Arrangement::Tetrahedral => sibling_skew(theta),
        Arrangement::Linear | Arrangement::Spread => 0.0,
    }
}

/// **四面体排布下,表角 `θ` 推出来的兄弟角与 `θ` 差多少**(弧度)。
///
/// 4 个取代基有 6 个夹角、而模掉整体转动只有 5 个自由度 —— 超定。
/// 构造法让"父–子"那几个精确等于 θ,兄弟之间的是推出来的:
/// 扭转差 120° 时 `cos φ = cos²θ + sin²θ·cos120°`。
///
/// 只有 θ **恰好** 109.4712° 时 φ = θ。表里的 109.4° 差 +0.1423°(可忽略);
/// θ 一远就不行:120° 差 −22.82°、104.5° 差 +9.45°、98.2° 差 +19.80°。
///
/// **判据要用它**:拿兄弟角去比表值是错的口径,该比的是"偏差有没有超过这个解析值"。
#[must_use]
pub fn sibling_skew(theta: f64) -> f64 {
    let phi = theta
        .cos()
        .mul_add(theta.cos(), theta.sin().powi(2) * -0.5)
        .clamp(-1.0, 1.0)
        .acos();
    (phi - theta).abs()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::{angle_at, dihedral, place_nerf, Point3};
    use std::f64::consts::PI;

    /// **sp² 中心必须共面。** 这是 `Planar` 那条存在的全部理由 ——
    /// 扭转角取 0/π 之外的任何值都会把它摆离面。
    #[test]
    fn a_planar_centre_really_comes_out_flat() {
        let g = Point3::new(0.0, 1.2, 0.3);
        let p = Point3::ORIGIN;
        let c = Point3::new(1.4, 0.0, 0.0);
        let ang = f64::to_radians(120.0);
        let ts = child_torsions(Arrangement::Planar, 2);
        assert_eq!(ts.len(), 2);
        let x: Vec<Point3> = ts
            .iter()
            .map(|t| place_nerf(g, p, c, 1.4, ang, *t).expect("摆得出来"))
            .collect();
        // 离面角:p、x0、x1 三个取代基与中心 c 应当共面 → 改进二面角 ≈ 0 或 π
        let d = dihedral(p, x[0], c, x[1]).expect("改进二面角");
        assert!(
            d.abs() < 1e-9 || (d.abs() - PI).abs() < 1e-9,
            "sp² 中心离面了,改进二面角 {:.4}°",
            d.to_degrees()
        );
        // 两个取代基之间也该是 120°
        let a = angle_at(x[0], c, x[1]).expect("角");
        assert!(
            (a.to_degrees() - 120.0).abs() < 1e-9,
            "取代基之间 {:.4}°",
            a.to_degrees()
        );
    }

    /// **sp³ 中心的三个取代基该两两 109.47°**(在角也取 109.471° 时)。
    #[test]
    fn a_tetrahedral_centre_spreads_its_three_children_evenly() {
        let g = Point3::new(0.0, 1.2, 0.3);
        let p = Point3::ORIGIN;
        let c = Point3::new(1.54, 0.0, 0.0);
        let ang = f64::to_radians(109.471);
        let ts = child_torsions(Arrangement::Tetrahedral, 3);
        assert_eq!(ts.len(), 3);
        let x: Vec<Point3> = ts
            .iter()
            .map(|t| place_nerf(g, p, c, 1.54, ang, *t).expect("摆得出来"))
            .collect();
        let mut n = 0;
        for i in 0..3 {
            for j in (i + 1)..3 {
                let a = angle_at(x[i], c, x[j]).expect("角").to_degrees();
                assert!((a - 109.471).abs() < 1e-3, "第 {i}/{j} 对是 {a:.4}°");
                n += 1;
            }
        }
        assert_eq!(n, 3);
    }

    /// **反式必须排第一** —— 只摆一个取代基时它得在 180°,链才伸展。
    #[test]
    fn the_first_child_is_always_anti() {
        for arr in [
            Arrangement::Linear,
            Arrangement::Planar,
            Arrangement::Tetrahedral,
            Arrangement::Spread,
        ] {
            let ts = child_torsions(arr, 1);
            assert_eq!(ts.len(), 1, "{arr:?} 要一个就该给一个");
            assert!(
                (ts[0] - PI).abs() < 1e-12,
                "{arr:?} 的第一个该是 180°,给的是 {:.3}°",
                ts[0].to_degrees()
            );
        }
    }

    /// 平面排布**只许**给 0 或 π,别的值一律是 bug。
    #[test]
    fn planar_never_hands_out_an_off_plane_torsion() {
        for n in 0..4 {
            for t in child_torsions(Arrangement::Planar, n) {
                assert!(
                    t.abs() < 1e-12 || (t - PI).abs() < 1e-12,
                    "平面排布给了 {:.3}°",
                    t.to_degrees()
                );
            }
        }
    }

    /// **同一个中心上,两个取代基不许拿到同一个扭转角。**
    ///
    /// 拿到同一个,加上同样的键长键角,NeRF 就会把它们摆到**同一个坐标**上。
    /// 上面那条 `planar_never_hands_out_an_off_plane_torsion` 断的是"值 ∈ {0, π}",
    /// 而 `[π, 0, π]` 照过 —— **它守不住这件事**,得单独有一条。
    #[test]
    fn no_two_children_ever_get_the_same_torsion() {
        for arr in [
            Arrangement::Linear,
            Arrangement::Planar,
            Arrangement::Tetrahedral,
            Arrangement::Spread,
        ] {
            for n in 0..8 {
                let ts = child_torsions(arr, n);
                for i in 0..ts.len() {
                    for j in (i + 1)..ts.len() {
                        let d = (ts[i] - ts[j]).abs();
                        let d = d.min(std::f64::consts::TAU - d);
                        assert!(
                            d > 1e-9,
                            "{arr:?} 要 {n} 个时第 {i} 与第 {j} 个都是 {:.3}° —— 两个原子会重合",
                            ts[i].to_degrees()
                        );
                    }
                }
            }
        }
    }

    /// **放不下就少给,不许硬凑。** 各排布在父键之外能放的上限。
    #[test]
    fn an_arrangement_only_hands_out_as_many_as_it_can_hold() {
        for (arr, cap) in [
            (Arrangement::Linear, 1),
            (Arrangement::Planar, 2),
            (Arrangement::Tetrahedral, 3),
        ] {
            for n in 0..8 {
                assert_eq!(
                    child_torsions(arr, n).len(),
                    n.min(cap),
                    "{arr:?} 要 {n} 个时该给 {} 个",
                    n.min(cap)
                );
            }
        }
    }

    /// **兄弟角容差必须与实际摆出来的几何对得上,而且不许无缘无故变大。**
    ///
    /// 这条把 `expected_sibling_skew` 与 `child_torsions` + NeRF **真的摆一遍**
    /// 的结果比对 —— 公式是判据的尺子,尺子错了判据就废了,而它错过两次
    /// (平面那一支先写成 0,再写成只在 θ ≥ 90° 成立的 `|2π−3θ|`)。
    #[test]
    fn the_sibling_tolerance_matches_the_geometry_it_claims_to_describe() {
        let g = Point3::new(0.0, 1.2, 0.3);
        let p = Point3::ORIGIN;
        let c = Point3::new(1.5, 0.0, 0.0);
        let mut checked = 0;
        for arr in [Arrangement::Planar, Arrangement::Tetrahedral] {
            for deg in [60.0_f64, 90.0, 104.5, 109.4712, 115.6, 120.0, 150.0] {
                let th = deg.to_radians();
                let ts = child_torsions(arr, 2);
                let x: Vec<Point3> = ts
                    .iter()
                    .map(|t| place_nerf(g, p, c, 1.4, th, *t).expect("摆得出来"))
                    .collect();
                assert_eq!(x.len(), 2, "{arr:?} 该给两个");
                let got = angle_at(x[0], c, x[1]).expect("兄弟角");
                let want = expected_sibling_skew(arr, th);
                // **容差 1e-6 不是 1e-9。** θ = 90° 时两个取代基恰好反向,
                // 兄弟角是 180°,而 `acos` 在 ±1 附近导数发散:点积上 1e-16 的
                // 浮点误差被放大成 **1.5e-8 rad** 的角度误差。那是数值条件,不是公式错。
                assert!(
                    ((got - th).abs() - want).abs() < 1e-6,
                    "{arr:?} θ={deg}°:实摆兄弟角 {:.4}°,偏差 {:.4}°,而公式说 {:.4}°",
                    got.to_degrees(),
                    (got - th).abs().to_degrees(),
                    want.to_degrees()
                );
                checked += 1;
            }
        }
        assert!(checked >= 14, "只验了 {checked} 组");
    }

    /// 配位数 ≥ 5 一律走 `Spread`,不管杂化写的是什么。
    #[test]
    fn five_or_more_neighbours_always_spread() {
        for h in [
            Hybridization::Sp3,
            Hybridization::Sp2,
            Hybridization::Sp3d,
            Hybridization::Unspecified,
        ] {
            for d in 5..8 {
                assert_eq!(arrangement(h, d), Arrangement::Spread, "{h:?} 配位 {d}");
            }
        }
    }

    /// `Spread` 要几个就给几个 —— 它是"均分了事",没有容量上限。
    ///
    /// 别的排布**有**上限,见 `an_arrangement_only_hands_out_as_many_as_it_can_hold`。
    #[test]
    fn spread_always_hands_out_exactly_what_was_asked() {
        for n in 0..8 {
            assert_eq!(
                child_torsions(Arrangement::Spread, n).len(),
                n,
                "Spread 要 {n} 个"
            );
        }
    }
}
