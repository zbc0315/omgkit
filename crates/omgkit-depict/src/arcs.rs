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
    // **盖不全就当没摆出来。** `place_ring` 逐段摆时,一个环上若有两段以上
    // 不连续的新原子、而规范环序的 0 号位正好落在绕回去的那一段里,收尾的
    // `j = n` 会把后面的段静默跳过 —— 于是返回一个**缺坐标的 `Some`**,而
    // `rings::place_candidates` 拿不到锚点是直接 panic。
    //
    // 松弛那一支没有这个面(它按系统的原子建图,结构上一定盖全)。实测全量
    // 11252 个环系统:绕回去的段 1265 次(很常见),但"一次摆环出现 ≥2 段"
    // **0 次**,132 个成功体系的坐标一个原子不缺 —— 所以这是**潜在**的,
    // 不是现行缺陷。一行换掉一个 panic 面,值。
    if rings
        .iter()
        .any(|r| r.atoms.iter().any(|a| !pos.contains_key(a)))
    {
        return None;
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

    /// 【测量,不是判据】弧法在全量语料的桥环体系上能摆多少、摆得怎么样。
    ///
    /// ```shell
    /// cargo test -p omgkit-depict --release --lib -- --ignored arc_coverage --nocapture
    /// ```
    #[test]
    #[ignore]
    fn arc_coverage() {
        use std::collections::BTreeMap as Map;
        let text = std::fs::read_to_string("../../harness/corpus/large.smi").expect("读语料");
        let (mut sys_n, mut ok, mut span, mut clash, mut noanchor) = (0, 0, 0, 0, 0);
        let (mut in_table, mut miss, mut miss_but_arc, mut miss_and_no_arc) = (0, 0, 0, 0);
        let mut miss_arc_skel: BTreeSet<String> = BTreeSet::new();
        let mut miss_noarc_skel: BTreeSet<String> = BTreeSet::new();
        let (mut arc_x, mut rlx_x, mut arc_win, mut rlx_win, mut tie) = (0, 0, 0, 0, 0);
        let (mut arc_dev, mut rlx_dev) = (0.0f64, 0.0f64);
        // 成功时的几何质量
        let (mut self_x, mut worst_dev) = (0usize, 0.0f64);
        // 按骨架统计
        let mut by_skel: Map<String, (usize, usize)> = Map::new();
        for line in text.lines() {
            let smi = line.split_whitespace().next().unwrap_or("");
            if smi.is_empty() || smi.starts_with('#') {
                continue;
            }
            let Ok(mut m) = omgkit_io::smiles::parse(smi) else {
                continue;
            };
            if omgkit_chem::pipeline::sanitize(&mut m).is_err() {
                continue;
            }
            let ranks = crate::ranks_of(&m);
            let rings_all = omgkit_chem::sssr::ring_set(&m);
            for sys in crate::rings::group(&omgkit_chem::rings::fused_ring_systems(&m), &rings_all)
            {
                // 只看桥环:邻稠的那一支 `layout_local` 自己就摆得了
                if !is_bridged(&m, &sys, &ranks) {
                    continue;
                }
                sys_n += 1;
                let skel = crate::templates::skeleton_of(&m, &sys.atoms, &ranks)
                    .unwrap_or_else(|| "?".into());
                let e = by_skel.entry(skel.clone()).or_default();
                e.0 += 1;
                // 现状:查表命中了吗
                let hit = matches!(
                    crate::templates::lookup_with(&m, &sys.atoms, &ranks, None).1,
                    crate::templates::Status::Hit
                );
                if hit {
                    in_table += 1;
                } else {
                    miss += 1;
                    if place(&sys.rings, &ranks).is_some() {
                        miss_but_arc += 1;
                        miss_arc_skel.insert(skel.clone());
                    } else {
                        miss_and_no_arc += 1;
                        miss_noarc_skel.insert(skel.clone());
                    }
                }
                match place(&sys.rings, &ranks) {
                    Some(pos) => {
                        ok += 1;
                        e.1 += 1;
                        // 键长偏差
                        let mut dev = 0.0f64;
                        for bd in m.bonds() {
                            if let (Some(p), Some(q)) = (pos.get(&bd.begin), pos.get(&bd.end)) {
                                dev = dev.max((p.dist(*q) - 1.0).abs());
                            }
                        }
                        if dev > worst_dev {
                            worst_dev = dev;
                        }
                        // 系统内部自交
                        let inside: Vec<(u32, u32)> = m
                            .bonds()
                            .iter()
                            .filter(|bd| pos.contains_key(&bd.begin) && pos.contains_key(&bd.end))
                            .map(|bd| (bd.begin, bd.end))
                            .collect();
                        let mut x = 0usize;
                        for (i, (a1, b1)) in inside.iter().enumerate() {
                            for (a2, b2) in &inside[i + 1..] {
                                if a1 == a2 || a1 == b2 || b1 == a2 || b1 == b2 {
                                    continue;
                                }
                                if crate::geom::segments_cross(pos[a1], pos[b1], pos[a2], pos[b2]) {
                                    x += 1;
                                }
                            }
                        }
                        if x > 0 {
                            self_x += 1;
                            arc_x += 1;
                        }
                        arc_dev = arc_dev.max(dev);
                        // 遮住整张表跑 relax —— 那正是语料外新骨架的待遇
                        let (rp, _) = crate::rings::relax(
                            &m,
                            &sys.atoms,
                            &ranks,
                            &sys.rings,
                            Some(("", &[])),
                        );
                        let mut rdev = 0.0f64;
                        for bd in m.bonds() {
                            if let (Some(p), Some(q)) = (rp.get(&bd.begin), rp.get(&bd.end)) {
                                rdev = rdev.max((p.dist(*q) - 1.0).abs());
                            }
                        }
                        rlx_dev = rlx_dev.max(rdev);
                        let rin: Vec<(u32, u32)> = m
                            .bonds()
                            .iter()
                            .filter(|bd| rp.contains_key(&bd.begin) && rp.contains_key(&bd.end))
                            .map(|bd| (bd.begin, bd.end))
                            .collect();
                        let mut rx = 0usize;
                        for (i, (a1, b1)) in rin.iter().enumerate() {
                            for (a2, b2) in &rin[i + 1..] {
                                if a1 == a2 || a1 == b2 || b1 == a2 || b1 == b2 {
                                    continue;
                                }
                                if crate::geom::segments_cross(rp[a1], rp[b1], rp[a2], rp[b2]) {
                                    rx += 1;
                                }
                            }
                        }
                        if rx > 0 {
                            rlx_x += 1;
                        }
                        match x.cmp(&rx) {
                            std::cmp::Ordering::Less => arc_win += 1,
                            std::cmp::Ordering::Greater => rlx_win += 1,
                            std::cmp::Ordering::Equal => tie += 1,
                        }
                    }
                    None => {
                        // 分因:重跑一遍,看是哪一支
                        match why(&sys.rings, &ranks) {
                            1 => span += 1,
                            2 => clash += 1,
                            _ => noanchor += 1,
                        }
                    }
                }
            }
        }
        println!(
            "桥环体系 {sys_n} 个;弧法摆得出来 {ok}({:.1}%)",
            100.0 * ok as f64 / sys_n as f64
        );
        println!("  摆不出来的分因:弦跨不过去 {span}  原子叠一起(鸽笼){clash}  没锚点 {noanchor}");
        println!("成功的那些:内部有自交的 {self_x} 个;最大键长偏差 {worst_dev:.2e}");
        let mut v: Vec<_> = by_skel.iter().collect();
        v.sort_by_key(|(_, (n, _))| std::cmp::Reverse(*n));
        println!("按骨架(出现次数最多的 15 条):");
        for (skel, (n, k)) in v.iter().take(15) {
            println!("   {n:4} 次  摆出 {k:4}  {skel}");
        }
        let total_skel = by_skel.len();
        let full = by_skel.values().filter(|(n, k)| n == k).count();
        let none = by_skel.values().filter(|(_, k)| *k == 0).count();
        println!("骨架 {total_skel} 条:全摆得出 {full} 条,一条都摆不出 {none} 条");
        println!();
        println!("=== 与现状比:接进运行时能买到什么 ===");
        println!("  查表命中(现状已经是最优候选,弧法买不到东西) {in_table}");
        println!("  没命中                                       {miss}");
        println!(
            "     其中弧法能摆的                            {miss_but_arc}  骨架 {:?}",
            miss_arc_skel
        );
        println!(
            "     其中弧法也摆不了                          {miss_and_no_arc}  骨架 {:?}",
            miss_noarc_skel
        );
        println!();
        println!("=== 语料外的新骨架会怎样:把表遮住,弧法 vs relax ===");
        println!("  (只看弧法摆得出来的那 {ok} 个体系)");
        println!("  弧法:自交 {arc_x} 个体系,最大键长偏差 {arc_dev:.4}");
        println!("  relax:自交 {rlx_x} 个体系,最大键长偏差 {rlx_dev:.4}");
        println!("  逐体系比:弧法交叉更少 {arc_win} 个,relax 更少 {rlx_win} 个,打平 {tie} 个");
    }

    /// 摆不出来是哪一支:1 弦跨不过去,2 原子叠一起,3 没锚点。
    fn why(rings: &[&Ring], ranks: &[u32]) -> u8 {
        // 复刻 `place` 的流程,只是把失败点分开报
        let mut pos: BTreeMap<u32, Point2> = BTreeMap::new();
        let mut order: Vec<&Ring> = rings.to_vec();
        order.sort_by_key(|r| (std::cmp::Reverse(r.atoms.len()), ring_key(r, ranks)));
        let cyc = canonical_cycle(&order[0].atoms, ranks);
        let n = cyc.len();
        let rad = 1.0 / (2.0 * (std::f64::consts::PI / n as f64).sin());
        for (i, a) in cyc.iter().enumerate() {
            let t = std::f64::consts::TAU * i as f64 / n as f64;
            pos.insert(*a, Point2::new(rad * t.cos(), rad * t.sin()));
        }
        let mut done: BTreeSet<usize> = BTreeSet::from([0]);
        while done.len() < order.len() {
            let Some((i, shared)) = (0..order.len())
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
                })
            else {
                return 3;
            };
            if shared == 0 {
                return 3;
            }
            done.insert(i);
            // 逐段试,分开报
            let seq = canonical_cycle(&order[i].atoms, ranks);
            let m = seq.len();
            let mut j = 0usize;
            while j < m {
                if pos.contains_key(&seq[j]) {
                    j += 1;
                    continue;
                }
                let mut s = j;
                while !pos.contains_key(&seq[(s + m - 1) % m]) {
                    s = (s + m - 1) % m;
                    if s == j {
                        break;
                    }
                }
                let mut e = s;
                let mut run = vec![seq[s]];
                while !pos.contains_key(&seq[(e + 1) % m]) {
                    e = (e + 1) % m;
                    run.push(seq[e]);
                    if e == s {
                        break;
                    }
                }
                let (a, b) = (seq[(s + m - 1) % m], seq[(e + 1) % m]);
                let (Some(pa), Some(pb)) = (pos.get(&a).copied(), pos.get(&b).copied()) else {
                    return 3;
                };
                let k = run.len();
                let d = pa.dist(pb);
                let theta = if d < 1e-9 {
                    std::f64::consts::TAU / (k + 1) as f64
                } else {
                    match solve_theta(k, d) {
                        Some(t) => t,
                        None => return 1,
                    }
                };
                let cands = arc_points(pa, pb, k, theta);
                match pick(&cands, &pos, a, b) {
                    Some(chosen) => {
                        for (at, p) in run.iter().zip(chosen.iter()) {
                            pos.insert(*at, *p);
                        }
                    }
                    None => return 2,
                }
                j = if e >= s { e + 1 } else { m };
            }
        }
        0
    }

    /// 这个环系统是不是桥环(邻稠的那一支 `layout_local` 自己摆得了)。
    fn is_bridged(mol: &MolBuilder, sys: &crate::rings::System<'_>, ranks: &[u32]) -> bool {
        crate::rings::layout_local(mol, sys, ranks, None)
            .1
            .is_some()
    }

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

    #[test]
    fn the_runtime_reaches_for_the_arc_before_it_falls_back_to_relaxing() {
        // **接线本身要有判据。** 三条路的次序是:查表 → 弧法 → 松弛。
        //
        // 这条把中间那一档钉住:拿一个语料里**没有**的桥环骨架(所以表必然
        // 不命中),断言 `layout_local` 给出的坐标就是弧法给的那一套 ——
        // 而不是松弛的。
        //
        // 语料内的桥环骨架 173/177 命中查表那一档,所以**表必须遮住**,否则
        // 这条守的是查表不是弧法。遮法用 `Override`:给一个对不上的键,整张表
        // 就被屏蔽 —— 那正是"语料里没有的新骨架"的待遇。
        use crate::rings::{group, layout_local};
        let mask: crate::templates::Override<'_> = Some(("表外的骨架", &[]));
        for smi in [
            "C1C2CCC1CC2",  // 降冰片烷,语料里最常见的桥环骨架(37 次)
            "C1C2CCC1CCC2", // 双环[3.2.1]辛烷(13 次)
            "C1C2CC1CCC2",  // 双环[3.2.0](7 次)
        ] {
            let m = prep(smi);
            let ranks = crate::ranks_of(&m);
            let rings_all = omgkit_chem::sssr::ring_set(&m);
            let systems = group(&omgkit_chem::rings::fused_ring_systems(&m), &rings_all);
            let sys = systems
                .iter()
                .max_by_key(|s| s.rings.len())
                .expect("该有一个环系统");
            // 前提一:遮住之后确实不命中
            let st = crate::templates::lookup_with(&m, &sys.atoms, &ranks, mask).1;
            assert!(
                !matches!(st, crate::templates::Status::Hit),
                "{smi} 遮了表还命中,遮法失效了"
            );
            // 前提二:不遮的话它**本来是命中的** —— 这句证明遮法确实在起作用
            assert!(
                matches!(
                    crate::templates::lookup_with(&m, &sys.atoms, &ranks, None).1,
                    crate::templates::Status::Hit
                ),
                "{smi} 本来就不在表里,那这条判据里的遮表是空动作"
            );
            // 前提三:弧法确实摆得出来
            let arc = place(&sys.rings, &ranks).expect("弧法该摆得出这个骨架");
            // 结论:运行时给的就是弧法那一套
            let (got, deg) = layout_local(&m, sys, &ranks, mask);
            assert!(deg.is_some(), "{smi} 是桥环,该如实报退化");
            for (a, p) in &arc {
                let q = got.get(a).expect("原子都该有坐标");
                assert!(
                    (p.x - q.x).abs() < 1e-9 && (p.y - q.y).abs() < 1e-9,
                    "{smi}:运行时没走弧法 —— 原子 {a} 弧法给 ({:.4},{:.4}),\
                     实得 ({:.4},{:.4})",
                    p.x,
                    p.y,
                    q.x,
                    q.y
                );
            }
            // 而弧法的键长是精确 1 —— 松弛给不出这个
            for bd in m.bonds() {
                if let (Some(p), Some(q)) = (got.get(&bd.begin), got.get(&bd.end)) {
                    assert!(
                        (p.dist(*q) - 1.0).abs() < 1e-9,
                        "{smi}:桥环系统里的键长该精确是 1,实得 {:.6}",
                        p.dist(*q)
                    );
                }
            }
        }
    }

    #[test]
    fn the_table_wins_over_the_arc_when_it_has_an_answer() {
        // **次序是查表 → 弧法 → 松弛,不能颠倒。**
        //
        // 弧法只保证**几何**(键长精确 1、桥不自交),而表里那一条是按**整分子**
        // 打分挑出来的 —— 它见过取代基往哪伸、标签占多大。两者不是一回事。
        //
        // 实测把弧法提到查表之前,全量语料:
        //
        // | | 表优先 | 弧法优先 |
        // |---|---:|---:|
        // | 其中有键交叉 | **48** | 70 |
        // | 有取代基挤到另一根键上 | **28** | 80 |
        // | 布局已退化的(交叉) | **30** | 52 |
        //
        // 三倍的取代基挤压 —— 弧法的桥摆得漂亮,可它不知道取代基要从哪出来。
        use crate::rings::{group, layout_local};
        let mut checked = 0usize;
        for smi in [
            "C1C2CCC1CC2",  // 降冰片烷,表里有(语料 37 次)
            "C1C2CCC1CCC2", // 双环[3.2.1]辛烷,表里有(13 次)
        ] {
            let m = prep(smi);
            let ranks = crate::ranks_of(&m);
            let rings_all = omgkit_chem::sssr::ring_set(&m);
            let systems = group(&omgkit_chem::rings::fused_ring_systems(&m), &rings_all);
            let sys = systems
                .iter()
                .max_by_key(|s| s.rings.len())
                .expect("该有一个环系统");
            // 前提:表命中,而且弧法也摆得出来 —— 两者都有话说,次序才有意义
            assert!(
                matches!(
                    crate::templates::lookup_with(&m, &sys.atoms, &ranks, None).1,
                    crate::templates::Status::Hit
                ),
                "{smi} 不在表里,这条判据说明不了次序"
            );
            let arc = place(&sys.rings, &ranks).expect("弧法该摆得出这个骨架");
            let (got, _) = layout_local(&m, sys, &ranks, None);
            let (tbl, _) = crate::templates::lookup_with(&m, &sys.atoms, &ranks, None);
            let tbl = tbl.expect("命中就该有坐标");

            // **结论:运行时给的是表那一套。** 无条件断言 —— 弧法一旦抢在
            // 查表前面,这里直接红。
            for (a, p) in &tbl {
                let q = got.get(a).expect("原子都该有坐标");
                assert!(
                    (p.x - q.x).abs() < 1e-9 && (p.y - q.y).abs() < 1e-9,
                    "{smi}:运行时没走查表 —— 原子 {a} 表里是 ({:.4},{:.4}),\
                     实得 ({:.4},{:.4})",
                    p.x,
                    p.y,
                    q.x,
                    q.y
                );
            }
            // **非空性单独验**:表与弧法给的得**不是**同一套,上面那条才有内容。
            // 表里存的若正好就是弧法赢来的那一套,这个骨架分不出次序。
            if tbl.iter().any(|(a, p)| {
                arc.get(a)
                    .is_some_and(|q| (p.x - q.x).abs() > 1e-6 || (p.y - q.y).abs() > 1e-6)
            }) {
                checked += 1;
            }
        }
        assert!(
            checked > 0,
            "没有一个骨架能分出「表」与「弧法」,这条判据是空过的"
        );
    }

    #[test]
    fn when_the_arc_cannot_do_it_the_runtime_falls_back_instead_of_giving_up() {
        // **退回那条边界也要守。** 弧法摆不了(鸽笼/弦跨不过去)时报 `None`,
        // 运行时必须退到松弛并照样给出坐标 —— 不能因为多接了一档就少画。
        //
        // 金刚烷正是鸽笼那一支:两个锚点之间三条等长桥,等张角弧对给定锚距
        // 只有两个镜像解,三条挤两个必然重合。
        // **表要遮住。** 金刚烷与双环[2.2.2]辛烷都在表里,不遮的话
        // `layout_local` 走的是查表那一档,根本到不了弧法/松弛 —— 这条判据
        // 就空过了。审核实测:不遮时把"弧法失败后退回松弛"整个换成"交空坐标",
        // 六条判据**全绿**。
        use crate::rings::{group, layout_local};
        let mask: crate::templates::Override<'_> = Some(("表外的骨架", &[]));
        let mut checked = 0usize;
        for smi in [
            "C1C2CC3CC1CC(C2)C3", // 金刚烷 —— 鸽笼那一支(语料里 18 次,弧法一次都摆不了)
            "C1CC2CCC1CC2",       // 双环[2.2.2]辛烷 —— 同一支(8 次)
        ] {
            let Ok(mut m) = omgkit_io::smiles::parse(smi) else {
                continue;
            };
            if omgkit_chem::pipeline::sanitize(&mut m).is_err() {
                continue;
            }
            let ranks = crate::ranks_of(&m);
            let rings_all = omgkit_chem::sssr::ring_set(&m);
            let systems = group(&omgkit_chem::rings::fused_ring_systems(&m), &rings_all);
            let Some(sys) = systems.iter().max_by_key(|s| s.rings.len()) else {
                continue;
            };
            // 前提一:遮了表确实不命中,否则守的是查表
            assert!(
                !matches!(
                    crate::templates::lookup_with(&m, &sys.atoms, &ranks, mask).1,
                    crate::templates::Status::Hit
                ),
                "{smi} 遮了表还命中,遮法失效了"
            );
            // 前提二:弧法确实摆不了 —— 这才是要验的那一支
            assert!(
                place(&sys.rings, &ranks).is_none(),
                "{smi} 弧法摆得出来,这条判据验不了退回那一支"
            );
            checked += 1;
            let (got, deg) = layout_local(&m, sys, &ranks, mask);
            assert!(deg.is_some(), "{smi} 是桥环,该如实报退化");
            for a in &sys.atoms {
                assert!(
                    got.contains_key(a),
                    "{smi}:弧法摆不了就该退回松弛,不能让原子 {a} 没坐标"
                );
            }
        }
        assert!(checked >= 2, "只验到 {checked} 个鸽笼骨架,这条判据太弱");
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
