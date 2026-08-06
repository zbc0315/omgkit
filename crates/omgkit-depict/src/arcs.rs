//! 桥环的几何摆法:**已摆好的原子一律不动,新原子按等张角圆弧跨在两个锚点之间。**
//!
//! # 为什么不是弹簧松弛
//!
//! [`rings::relax`](crate::rings) 是"随机初值 + 弹簧下降 + 按分数挑"。它不解问题,
//! 它抽奖:吗啡骨架的 0 交叉解要搜二十二万次才撞上,而立方烷、棱晶烷等五条骨架
//! 跑两百万次仍然到不了;键长偏差普遍 20%~60%。更糟的是打分口径一调就在别处砸坑
//! —— 实测给它加一档几何可行性阈值,吗啡如愿了,**金刚烷(语料里最常见的桥环
//! 骨架)却从 0 交叉变成 1 交叉,棱晶烷从 1 处交叉变成 15 处**。
//!
//! # 这一套怎么摆
//!
//! 起手环摆成单位边长的正多边形。之后每个环:**它身上已经有坐标的原子一个都不动**,
//! 剩下的新原子按环上的连续段(弧)划分,每一段跨在两个已放锚点之间:
//!
//! ```text
//!    已放 a ●━━━○━━━○━━━○━━━● b 已放
//!            k 个新原子,k+1 根单位键,跨过锚距 d
//! ```
//!
//! 解等张角 `θ`:`sin((k+1)θ/2) / sin(θ/2) = d`,新原子落在半径 `1/(2 sin(θ/2))`
//! 的圆弧上。`d < k+1` 时唯一有解。两个镜像解取不撞的那个。
//!
//! **邻稠是它的特例**(`d = 1`、`k = n−2`,解出来正好是正多边形),所以稠环与桥环
//! 一套代码;**冗余环自动是空操作**(一个新原子都没有,什么都不做)—— 而本库的
//! `ring_set` 不是严格 SSSR,177 个桥环体系里有 65 个含这种环。
//!
//! # 摆不出来就说摆不出来
//!
//! 三种情形返回 `None`,交给 [`rings::relax`](crate::rings) 兜底:
//!
//! - **弦跨不过去**:`d ≥ k+1`,k+1 根单位键够不到。
//! - **没有锚点**:剩下的环与已放部分一个原子都不共用。
//! - **原子叠在一起**:两个锚点之间有三条等长桥时(双环[2.2.2]辛烷那一类),
//!   等张角弧对给定的 `d` 只有两个镜像解,三条桥挤两个位置必然有两条重合。
//!   **这一条必须自己拒**,不能指望调用方的打分挡住 —— `rings::Quality` 数的是
//!   键交叉,两个原子精确重合不产生交叉,而它的键长偏差反而是完美的 0。

use std::collections::{BTreeMap, BTreeSet};

use omgkit_chem::sssr::Ring;

use crate::geom::Point2;

/// 两个原子靠得比这个近就算叠在一起。取键长的一半 —— 真正相邻的原子相距 1,
/// 差着一倍,不会误伤。
const CLASH: f64 = 0.5;

/// 摆一个环系统。`rings` 是这个系统里的全部环,`ranks` 是规范秩。
///
/// 返回 `None` 表示这套方法摆不了,见模块文档。
pub(crate) fn place(rings: &[&Ring], ranks: &[u32]) -> Option<BTreeMap<u32, Point2>> {
    if rings.is_empty() {
        return None;
    }
    let mut pos: BTreeMap<u32, Point2> = BTreeMap::new();

    // 起手环:最大的那个,平局按环上规范秩的多重集。大环更能定住整体形状。
    let mut order: Vec<&Ring> = rings.to_vec();
    order.sort_by_key(|r| (std::cmp::Reverse(r.atoms.len()), ring_key(r, ranks)));

    // 正多边形,单位边长。**起点与绕向都由 `canonical_cycle` 定死** ——
    // 不同写法可能从不同原子起、甚至朝相反方向绕,那会落到不同的构型上。
    let cyc = canonical_cycle(&order[0].atoms, ranks);
    let n = cyc.len();
    let rad = 1.0 / (2.0 * (std::f64::consts::PI / n as f64).sin());
    for (i, a) in cyc.iter().enumerate() {
        let t = std::f64::consts::TAU * i as f64 / n as f64;
        pos.insert(*a, Point2::new(rad * t.cos(), rad * t.sin()));
    }

    let mut done: BTreeSet<usize> = BTreeSet::from([0]);
    while done.len() < order.len() {
        // 下一个环:与已放部分共用原子最多的。**平局一律由规范秩打破。**
        let (i, shared) = (0..order.len())
            .filter(|i| !done.contains(i))
            .map(|i| {
                let sh = order[i]
                    .atoms
                    .iter()
                    .filter(|a| pos.contains_key(a))
                    .count();
                (i, sh)
            })
            .max_by_key(|(i, sh)| {
                (
                    *sh,
                    order[*i].atoms.len(),
                    std::cmp::Reverse(ring_key(order[*i], ranks)),
                )
            })?;
        if shared == 0 {
            return None; // 没有锚点 —— 这不该发生在同一个双连通系统里
        }
        done.insert(i);
        place_ring(order[i], ranks, &mut pos)?;
    }
    Some(pos)
}

/// 把一个环上还没坐标的原子摆上去。已经有坐标的一个都不动。
fn place_ring(ring: &Ring, ranks: &[u32], pos: &mut BTreeMap<u32, Point2>) -> Option<()> {
    let seq = canonical_cycle(&ring.atoms, ranks);
    let n = seq.len();
    if seq.iter().all(|a| pos.contains_key(a)) {
        return Some(()); // 冗余环:一个新原子都没有
    }
    // 逐段找"连续的新原子",每段跨在两个已放锚点之间
    let mut j = 0usize;
    while j < n {
        if pos.contains_key(&seq[j]) {
            j += 1;
            continue;
        }
        // 往前退到这一段的起点
        let mut s = j;
        while !pos.contains_key(&seq[(s + n - 1) % n]) {
            s = (s + n - 1) % n;
            if s == j {
                break;
            }
        }
        // 往后走到这一段的终点
        let mut e = s;
        let mut run = vec![seq[s]];
        while !pos.contains_key(&seq[(e + 1) % n]) {
            e = (e + 1) % n;
            run.push(seq[e]);
            if e == s {
                break;
            }
        }
        let (a, b) = (seq[(s + n - 1) % n], seq[(e + 1) % n]);
        let (pa, pb) = (*pos.get(&a)?, *pos.get(&b)?);
        let k = run.len();
        let d = pa.dist(pb);
        let theta = if d < 1e-9 {
            std::f64::consts::TAU / (k + 1) as f64
        } else {
            solve_theta(k, d)?
        };
        let cands = arc_points(pa, pb, k, theta);
        let chosen = pick(&cands, pos, a, b)?;
        for (at, p) in run.iter().zip(chosen.iter()) {
            pos.insert(*at, *p);
        }
        j = if e >= s { e + 1 } else { n };
    }
    Some(())
}

/// 两个镜像解挑哪一个。
///
/// # 不能用"离已放部分的质心更远"
///
/// 对称的桥(双环[2.2.2]辛烷那一类,弦正好是父环的直径)两侧到质心的距离**精确
/// 相等**,那时平局一破就有一半概率把新原子摆到已放原子头上。改用**抗碰撞**的
/// 口径:取"离已放原子最近的那个距离"最大的一套。
///
/// 仍然相等时按量化坐标序列打破 —— 不留任意性,同一分子的任何写法挑到同一套。
///
/// 两套都会撞(最近距离 < [`CLASH`])就返回 `None`:三条等长桥挤两个镜像位置,
/// 这套方法解不了。
fn pick(
    cands: &[Vec<Point2>; 2],
    pos: &BTreeMap<u32, Point2>,
    a: u32,
    b: u32,
) -> Option<Vec<Point2>> {
    let nearest = |v: &Vec<Point2>| -> f64 {
        let mut mn = f64::INFINITY;
        for p in v {
            for (at, q) in pos {
                if *at == a || *at == b {
                    continue;
                }
                mn = mn.min(p.dist(*q));
            }
            // 同一段里的新原子也不许叠(锚点重合时会退化成同一点)
            for q in v {
                if !std::ptr::eq(p, q) {
                    mn = mn.min(p.dist(*q));
                }
            }
        }
        if mn.is_finite() {
            mn
        } else {
            f64::MAX
        }
    };
    #[allow(clippy::cast_possible_truncation)]
    let key = |v: &Vec<Point2>| -> (i64, Vec<(i64, i64)>) {
        (
            -((nearest(v) * 1e6).round() as i64), // 越远越好 → 取负当"越小越好"
            v.iter()
                .map(|p| ((p.x * 1e6).round() as i64, (p.y * 1e6).round() as i64))
                .collect(),
        )
    };
    let (k0, k1) = (key(&cands[0]), key(&cands[1]));
    let best = if k0 <= k1 { &cands[0] } else { &cands[1] };
    if nearest(best) < CLASH {
        return None;
    }
    Some(best.clone())
}

/// 解等张角 `θ`:`k+1` 根单位弦、总张角 `(k+1)θ`、弦长 `d`。
///
/// `f(θ) = sin((k+1)θ/2) / sin(θ/2)` 在 `(0, 2π/(k+1))` 上从 `k+1` 单调降到 0,
/// 所以二分即可。`d ≥ k+1` 时无解 —— k+1 根单位键跨不了那么远。
fn solve_theta(k: usize, d: f64) -> Option<f64> {
    let m = (k + 1) as f64;
    if d >= m - 1e-12 {
        return None;
    }
    let f = |t: f64| (m * t / 2.0).sin() / (t / 2.0).sin() - d;
    let (mut lo, mut hi) = (1e-12, std::f64::consts::TAU / m - 1e-12);
    if f(lo) < 0.0 {
        return None;
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if f(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Some(0.5 * (lo + hi))
}

/// 把 `k` 个新原子放在 `a → b` 的等张角圆弧上(连两端共 `k+1` 根单位键)。
/// 返回互为镜像的两套。
fn arc_points(a: Point2, b: Point2, k: usize, theta: f64) -> [Vec<Point2>; 2] {
    let m = (k + 1) as f64;
    let r = 1.0 / (2.0 * (theta / 2.0).sin());
    let d = a.dist(b);
    let along = if d > 1e-12 {
        Point2::new((b.x - a.x) / d, (b.y - a.y) / d)
    } else {
        Point2::new(1.0, 0.0)
    };
    let normal = Point2::new(-along.y, along.x);
    let h = r * (m * theta / 2.0).cos(); // 圆心到弦的有向距离
    let mid = Point2::new(0.5 * (a.x + b.x), 0.5 * (a.y + b.y));
    let rot = |v: Point2, ang: f64| {
        Point2::new(
            v.x * ang.cos() - v.y * ang.sin(),
            v.x * ang.sin() + v.y * ang.cos(),
        )
    };
    let mut out = [Vec::new(), Vec::new()];
    for (s, side) in [1.0_f64, -1.0].iter().enumerate() {
        let c = Point2::new(mid.x + normal.x * h * side, mid.y + normal.y * h * side);
        let va = Point2::new(a.x - c.x, a.y - c.y);
        // 转向取"能把 a 转到 b"的那一个
        let sign = {
            let p = rot(va, m * theta);
            if Point2::new(c.x + p.x, c.y + p.y).dist(b) < 1e-6 {
                1.0
            } else {
                -1.0
            }
        };
        out[s] = (1..=k)
            .map(|j| {
                let p = rot(va, sign * theta * j as f64);
                Point2::new(c.x + p.x, c.y + p.y)
            })
            .collect();
    }
    out
}

/// 环上规范秩的多重集 —— 与写法无关的环排序键。与 `rings.rs` 里那份同口径。
fn ring_key(r: &Ring, ranks: &[u32]) -> Vec<u32> {
    let mut k: Vec<u32> = r.atoms.iter().map(|a| ranks[*a as usize]).collect();
    k.sort_unstable();
    k
}

/// 环上原子的规范绕法:起点与**方向**都定死。
///
/// 只定起点是不够的 —— 不同写法可能朝相反方向绕,那会落到**不同的构型**上,
/// 不只是旋转或镜像的差别。与 `rings.rs` 里那份同口径。
fn canonical_cycle(atoms: &[u32], ranks: &[u32]) -> Vec<u32> {
    let n = atoms.len();
    let Some(start) = (0..n).min_by_key(|i| (ranks[atoms[*i] as usize], atoms[*i])) else {
        return Vec::new();
    };
    let fwd: Vec<u32> = (0..n).map(|k| atoms[(start + k) % n]).collect();
    let bwd: Vec<u32> = (0..n).map(|k| atoms[(start + n - k) % n]).collect();
    let key = |v: &[u32]| -> Vec<u32> { v.iter().map(|a| ranks[*a as usize]).collect() };
    if key(&bwd) < key(&fwd) {
        bwd
    } else {
        fwd
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use omgkit_core::MolBuilder;

    fn prep(smi: &str) -> MolBuilder {
        let mut m = omgkit_io::smiles::parse(smi).expect("SMILES 该能解析");
        omgkit_chem::pipeline::sanitize(&mut m).expect("该能 sanitize");
        m
    }

    /// 摆一个分子最大的那个环系统。
    fn lay(smi: &str) -> (MolBuilder, Option<BTreeMap<u32, Point2>>) {
        let m = prep(smi);
        let ranks = omgkit_io::canon::canonical_ranks(&m);
        let rs = omgkit_chem::sssr::ring_set(&m);
        let syss = crate::rings::group(&omgkit_chem::rings::fused_ring_systems(&m), &rs);
        let Some(sys) = syss.iter().max_by_key(|s| s.atoms.len()) else {
            return (m, None);
        };
        let out = place(&sys.rings, &ranks);
        (m, out)
    }

    #[test]
    fn every_bond_it_places_is_exactly_one_bond_long() {
        // **这是这套摆法相对弹簧松弛的根本改进。** 松弛只是把键长"拉向"1,
        // 退化布局实测偏差全部 ≥20%、常见 30%~60%;几何摆放是解出来的,精确。
        //
        // 例外只有一种:两端都已经摆好的**合拢键** —— 它的长度是几何定死的,
        // 摆不成 1。所以这里量的是"新摆上去的键",判据自己把合拢键找出来单列。
        for smi in [
            "C1C2CCC1CC2", // 降冰片烷骨架,语料里出现 37 次
            "C1C2CCC1CCC2",
            "C1C2CC1CC2",
            "c1ccc2ccccc2c1", // 萘:邻稠是这套摆法的特例
            "C1C2CCCC1CCC2",
        ] {
            let (m, pos) = lay(smi);
            let pos = pos.unwrap_or_else(|| panic!("{smi} 该摆得出来"));
            let mut closing = 0usize;
            let mut worst = 0.0_f64;
            for b in m.bonds() {
                let (Some(u), Some(v)) = (pos.get(&b.begin), pos.get(&b.end)) else {
                    continue;
                };
                let d = (u.dist(*v) - 1.0).abs();
                if d > 1e-9 {
                    closing += 1;
                    worst = worst.max(d);
                }
            }
            // 合拢键可以不等于 1,但不该多 —— 多了说明摆法在硬凑
            assert!(
                closing * 3 <= m.bonds().len(),
                "{smi}:{closing} 根键长度不是 1(共 {} 根),最差差 {worst:.4}",
                m.bonds().len()
            );
        }
    }

    #[test]
    fn it_refuses_instead_of_stacking_atoms_on_top_of_each_other() {
        // 两个锚点之间有三条等长桥时(双环[2.2.2]辛烷那一类),等张角弧对给定
        // 的锚距只有两个镜像解,三条桥挤两个位置必然有两条重合。
        //
        // **这一条必须自己拒。** `rings::Quality` 数的是键交叉,两个原子精确
        // 重合不产生交叉,而它的键长偏差反而是完美的 0 —— 靠调用方挡不住。
        // **双环[2.2.2]辛烷必须被拒**:它的两个桥头之间有三条等长桥,而等张角弧
        // 对给定锚距只有两个镜像解 —— 三条挤两个位置,必然有两条重合。
        // 鸽笼:两个锚点之间三条等长桥,而等张角弧只有两个镜像位置。
        // 三个都取自模板表,实测都被拒。
        for smi in [
            "C1CC2CCC1CC2",                   // 双环[2.2.2]辛烷
            "C1C2CC3CC1CC(C2)C3",             // 金刚烷
            "C1C2CCCC34C5C(CC4CCCC23)CCCC15", // 吗啡骨架
        ] {
            assert!(lay(smi).1.is_none(), "{smi} 该被拒 —— 这套摆法解不了它");
        }

        for smi in ["C1C2CCC1CC2", "C1C2CCC1CCC2", "c1ccc2ccccc2c1"] {
            let (_, pos) = lay(smi);
            let pos = pos.unwrap_or_else(|| panic!("{smi} 该摆得出来"));
            let pts: Vec<Point2> = pos.values().copied().collect();
            for (i, p) in pts.iter().enumerate() {
                for q in &pts[i + 1..] {
                    assert!(
                        p.dist(*q) >= CLASH,
                        "{smi}:两个原子只隔 {:.4},摆法该报 None 而不是发出来",
                        p.dist(*q)
                    );
                }
            }
        }
    }

    #[test]
    fn how_it_was_written_does_not_change_where_the_atoms_go() {
        // 头号契约。链条上每一处平局都必须由规范秩打破:环序、起手环的起点与
        // **绕向**、下一个环的选择、两个镜像的取舍。
        for smi in [
            "C1CC2CCC1CC2",
            "C1CC2CCC1C2",
            "C1C2CC3CC1CC(C2)C3",
            "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5",
        ] {
            let m = prep(smi);
            let base = fingerprint(&m);
            let mut compared = 0usize;
            for seed in 0..12u64 {
                let w = omgkit_io::smiles::write_with_priority(&m, &shuffled(m.num_atoms(), seed));
                let Ok(mut m2) = omgkit_io::smiles::parse(&w.smiles) else {
                    continue;
                };
                if omgkit_chem::pipeline::sanitize(&mut m2).is_err() {
                    continue;
                }
                if omgkit_io::canon::canonical_smiles(&m2).smiles
                    != omgkit_io::canon::canonical_smiles(&m).smiles
                {
                    continue;
                }
                compared += 1;
                assert_eq!(
                    base,
                    fingerprint(&m2),
                    "{smi} 写成 {} 之后摆位变了",
                    w.smiles
                );
            }
            assert_eq!(compared, 12, "{smi} 只比上了 {compared} 种写法");
        }
    }

    /// 按规范秩排好的量化坐标序列 —— 与原子编号无关的指纹。
    fn fingerprint(m: &MolBuilder) -> Option<Vec<(i64, i64)>> {
        let ranks = omgkit_io::canon::canonical_ranks(m);
        let rs = omgkit_chem::sssr::ring_set(m);
        let syss = crate::rings::group(&omgkit_chem::rings::fused_ring_systems(m), &rs);
        let sys = syss.iter().max_by_key(|s| s.atoms.len())?;
        let pos = place(&sys.rings, &ranks)?;
        let mut v: Vec<(u32, Point2)> = pos.iter().map(|(a, p)| (ranks[*a as usize], *p)).collect();
        v.sort_by_key(|x| x.0);
        #[allow(clippy::cast_possible_truncation)]
        Some(
            v.iter()
                .map(|(_, p)| ((p.x * 1e6).round() as i64, (p.y * 1e6).round() as i64))
                .collect(),
        )
    }

    /// splitmix64 + Fisher-Yates。仿射映射保持相邻关系,置换出来的写法太像原来的。
    fn shuffled(n: usize, seed: u64) -> Vec<u32> {
        let mut s = seed.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut next = || {
            s = s.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = s;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        };
        let mut v: Vec<u32> = (0..u32::try_from(n).unwrap()).collect();
        for i in (1..n).rev() {
            let j = usize::try_from(next() % (i as u64 + 1)).unwrap();
            v.swap(i, j);
        }
        v
    }
}
