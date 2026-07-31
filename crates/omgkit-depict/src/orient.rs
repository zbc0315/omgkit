//! 规范朝向:同一个分子的任何写法,摆出来**逐点相同**。
//!
//! # 为什么要有这一步
//!
//! 布局本身给出的朝向取决于起手环的原子顺序,虽然形状是确定的,摆放的角度却
//! 不是。反应式里同一个分子出现两次却歪着不同的角度,读起来很别扭。
//!
//! 更要紧的是判据:先前那条"写法无关"用的是**两两距离的多重集**当指纹,而那
//! 个指纹**镜像不变** —— 两种写法给出互为镜像的图,它一点都发现不了。规范朝向
//! 把镜像这一维也定死之后,判据可以升级成最强的形式:**按规范秩排好的坐标序列
//! 逐点相同**。
//!
//! # 只在 30° 的整数倍上转
//!
//! 主轴对齐(PCA)会把图转到一个任意角度,于是环上那些本来齐整的 120°、30°
//! 全歪了。改成在 24 个候选姿态里挑(12 个 30° 倍数的旋转 × 镜不镜像),既拿到
//! 横向的版式,又保住了齐整的角度。这也是 IUPAC 与各家工具箱的通行做法
//! (Mayfield 在 RDKit UGM 2016 上把它叫 "bond snapping to 30°")。
//!
//! # 镜像是安全的,但顺序有讲究
//!
//! 镜像会把手性画反 —— 所以这一步必须排在**楔形指派之前**。楔形是照最终坐标
//! 算出来的,先摆正再指派,构型自然是对的;反过来做就会把已经画好的楔形悬空。
//! 顺反在镜像下不变(同侧还是同侧),不受影响。

use crate::geom::Point2;

/// 候选姿态的角度步长:30°。
const STEP: f64 = std::f64::consts::FRAC_PI_6;

/// 比较坐标时的量化精度。浮点直接比大小会让"选哪个姿态"取决于最后一位,
/// 而那一位取决于运算次序 —— 同一个分子的不同写法就会挑到不同的姿态。
const QUANT: f64 = 1e6;

/// 把一张图摆成规范姿态。就地修改。
///
/// `ranks` 是规范秩;平局全靠它打破,拿存储下标打破会引入写法依赖。
pub(crate) fn canonicalise(coords: &mut [Point2], ranks: &[u32]) {
    if coords.len() < 2 {
        return;
    }
    // 先挪到质心,让旋转与镜像都绕着一个与编号无关的点做
    let n = coords.len() as f64;
    let c = coords.iter().fold(Point2::ORIGIN, |s, p| s + *p) * (1.0 / n);
    for p in coords.iter_mut() {
        *p = *p - c;
    }

    // 按规范秩排好的下标 —— 指纹要按这个顺序取,才与原子编号无关
    let mut order: Vec<usize> = (0..coords.len()).collect();
    order.sort_by_key(|i| (ranks[*i], *i));

    let mut best: Option<(Key, Vec<Point2>)> = None;
    for mirror in [false, true] {
        for k in 0..12 {
            let cand: Vec<Point2> = coords
                .iter()
                .map(|p| {
                    let q = if mirror { Point2::new(p.x, -p.y) } else { *p };
                    q.rotated(STEP * f64::from(k))
                })
                .collect();
            let key = key_of(&cand, &order);
            // 不用 `Option::is_none_or` —— 它到 1.82 才稳定,工作区 MSRV 是 1.75
            let take = match &best {
                None => true,
                Some((b, _)) => key < *b,
            };
            if take {
                best = Some((key, cand));
            }
        }
    }
    let (_, chosen) = best.expect("12×2 个候选里必有一个");
    coords.copy_from_slice(&chosen);
}

/// 一个姿态的排序键。**越小越好。**
///
/// 先要横向的版式(高比宽小),再拿量化后的坐标序列打破平局 —— 后者对高度
/// 对称的分子(苯、萘)是唯一的决定因素,那时若不定死就会随写法漂移。
type Key = (i64, Vec<(i64, i64)>);

fn key_of(pts: &[Point2], order: &[usize]) -> Key {
    let (mut x0, mut y0, mut x1, mut y1) = (f64::MAX, f64::MAX, f64::MIN, f64::MIN);
    for p in pts {
        x0 = x0.min(p.x);
        y0 = y0.min(p.y);
        x1 = x1.max(p.x);
        y1 = y1.max(p.y);
    }
    let (w, h) = (x1 - x0, y1 - y0);
    // 高宽比:越小越扁,越扁越好。量化之后同一个比值不会因末位差别分出高下。
    let flat = ((h / w.max(1e-9)) * QUANT).round() as i64;
    let seq: Vec<(i64, i64)> = order
        .iter()
        .map(|i| {
            (
                (pts[*i].x * QUANT).round() as i64,
                (pts[*i].y * QUANT).round() as i64,
            )
        })
        .collect();
    (flat, seq)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{generate, style::Style};
    use omgkit_core::MolBuilder;

    fn prep(smi: &str) -> MolBuilder {
        let mut m = omgkit_io::smiles::parse(smi).unwrap();
        omgkit_chem::pipeline::sanitize(&mut m).unwrap();
        m
    }

    /// 按规范秩排好的坐标序列 —— 与原子编号无关的完整指纹。
    fn canonical_coords(smi: &str, style: &Style) -> Vec<(i64, i64)> {
        let m = prep(smi);
        let ranks = omgkit_io::canon::canonical_ranks(&m);
        let d = generate(&m, style);
        let mut order: Vec<usize> = (0..d.coords.len()).collect();
        order.sort_by_key(|i| (ranks[*i], *i));
        order
            .iter()
            .map(|i| {
                (
                    (d.coords[*i].x * 1e4).round() as i64,
                    (d.coords[*i].y * 1e4).round() as i64,
                )
            })
            .collect()
    }

    #[test]
    fn different_writings_give_literally_the_same_coordinates() {
        // **判据升级到最强的形式。** 先前那条用两两距离的多重集当指纹,而那个
        // 指纹镜像不变 —— 两种写法给出互为镜像的图,它一点都发现不了。规范朝向
        // 把镜像这一维也定死之后,可以直接比坐标。
        let groups = [
            vec![
                "CC(=O)Oc1ccccc1C(=O)O",
                "O=C(C)Oc1ccccc1C(O)=O",
                "OC(=O)c1ccccc1OC(C)=O",
            ],
            vec!["c1ccc2ccccc2c1", "c1ccc2c(c1)cccc2", "c1cc2ccccc2cc1"],
            vec!["CC(C)(C)c1ccccc1", "c1ccccc1C(C)(C)C"],
            vec!["CCCCO", "OCCCC"],
            vec!["CN1C=NC2=C1C(=O)N(C)C(=O)N2C", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"],
        ];
        for style in &Style::ALL {
            for ws in &groups {
                let a = canonical_coords(ws[0], style);
                for w in &ws[1..] {
                    assert_eq!(
                        a,
                        canonical_coords(w, style),
                        "[{}] {w} 与 {} 坐标不同",
                        style.name,
                        ws[0]
                    );
                }
            }
        }
    }

    #[test]
    fn the_picture_comes_out_wider_than_tall() {
        // 期刊版式是横向的。竖着摆不报错,只是排版时要转一下 —— 而"要不要转"
        // 这件事本该在出图时就定好。
        for smi in ["CCCCCCCCCC", "c1ccc2ccccc2c1", "CC(=O)Oc1ccccc1C(=O)O"] {
            let d = generate(&prep(smi), &Style::ACS_1996);
            let (mut x0, mut y0, mut x1, mut y1) = (f64::MAX, f64::MAX, f64::MIN, f64::MIN);
            for p in &d.coords {
                x0 = x0.min(p.x);
                y0 = y0.min(p.y);
                x1 = x1.max(p.x);
                y1 = y1.max(p.y);
            }
            assert!(
                x1 - x0 >= y1 - y0 - 1e-9,
                "{smi} 摆成了竖的:{:.2} 宽 × {:.2} 高",
                x1 - x0,
                y1 - y0
            );
        }
    }

    #[test]
    fn rings_keep_their_tidy_angles() {
        // 主轴对齐(PCA)会把图转到任意角度,环上本来齐整的 30°/120° 就全歪了。
        // 只在 30° 的整数倍上转才保得住。这条量的正是键的角度。
        let m = prep("c1ccccc1");
        let d = generate(&m, &Style::ACS_1996);
        for b in m.bonds() {
            let v = d.coords[b.end as usize] - d.coords[b.begin as usize];
            let deg = v.angle().to_degrees().rem_euclid(30.0);
            let off = deg.min(30.0 - deg);
            assert!(
                off < 1e-6,
                "键 {}–{} 偏离 30° 的整数倍 {off:.4}°",
                b.begin,
                b.end
            );
        }
    }

    #[test]
    fn orienting_does_not_move_atoms_relative_to_each_other() {
        // 摆正只该是刚体变换(旋转 + 可能的镜像)。若不小心写成了缩放或投影,
        // 键长会跟着变 —— 而"键长全等"那条判据在别处,这里单独守一次形状。
        let m = prep("CC(=O)Oc1ccccc1C(=O)O");
        let ranks = omgkit_io::canon::canonical_ranks(&m);
        let d = generate(&m, &Style::ACS_1996);
        let mut moved = d.coords.clone();
        canonicalise(&mut moved, &ranks);
        // 再摆一次必须是不动点 —— 否则规范化本身不收敛
        for (a, b) in d.coords.iter().zip(&moved) {
            assert!(a.dist(*b) < 1e-9, "已经摆正的图再摆一次又动了");
        }
    }
}
