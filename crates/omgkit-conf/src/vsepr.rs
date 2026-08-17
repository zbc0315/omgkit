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
        Hybridization::Sp3 | Hybridization::Sp2d => Arrangement::Tetrahedral,
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
/// 扭转角是 `祖父–父–中心–取代基` 那个二面角。约定:
///
/// - [`Arrangement::Linear`]:只摆一个,扭转角无意义(角是 180°,
///   NeRF 里 `sin(角)=0`,扭转项自动消掉),给 π 占位;
/// - [`Arrangement::Planar`]:**只能是 0 或 π** —— 别的值会让 sp² 中心离面。
///   两个取代基就是 `[π, 0]`(先反式,后顺式);
/// - [`Arrangement::Tetrahedral`]:从 π 起按 `2π/3` 均分 → `[π, π/3, 5π/3]`,
///   这就是标准的**交错**构象;
/// - [`Arrangement::Spread`]:`n+1` 个方向均分(把父键也算一个)。
///
/// **反式(π)永远排第一** —— 链要伸展开,而伸展的那一支应当优先。
#[must_use]
pub fn child_torsions(arr: Arrangement, n: usize) -> Vec<f64> {
    use std::f64::consts::{PI, TAU};
    match arr {
        Arrangement::Linear => vec![PI; n.min(1)],
        Arrangement::Planar => (0..n).map(|k| if k % 2 == 0 { PI } else { 0.0 }).collect(),
        Arrangement::Tetrahedral => (0..n)
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

    /// 要几个就给几个,不多不少。
    #[test]
    fn the_count_always_matches_what_was_asked() {
        for arr in [
            Arrangement::Planar,
            Arrangement::Tetrahedral,
            Arrangement::Spread,
        ] {
            for n in 0..5 {
                assert_eq!(child_torsions(arr, n).len(), n, "{arr:?} 要 {n} 个");
            }
        }
        // Linear 最多一个 —— 直线中心塞不下第二个取代基
        assert_eq!(child_torsions(Arrangement::Linear, 3).len(), 1);
    }
}
