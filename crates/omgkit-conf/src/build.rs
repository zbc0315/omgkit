//! **树遍历摆链** —— 把 [`crate::geom`]、[`crate::params`]、[`crate::vsepr`] 接起来。
//!
//! # 这一步只做无环分子
//!
//! 有环的分子这一期整个不摆(如实计数),环系要等二期的距离几何。
//! 这是**范围**不是失败:摆不了就说摆不了,不硬摆。
//!
//! # 为什么链不需要距离几何
//!
//! 链是树,没有闭合约束,所以每个原子都能由(键长、键角、扭转角)**闭式**定出来
//! ([`crate::geom::place_nerf`]),`O(N)`、零迭代、不会失败。
//! RDKit 把链也塞进全分子的距离几何里,那是它 `O(N³)` 与长尾的来源之一。
//!
//! # 规范序是**传进来**的,不是这里算的
//!
//! 摆放次序决定坐标,所以次序必须确定。这里要求调用方给一份
//! `ranks`(`omgkit_io::canon::classed_ranks`),所有"挑下一个"的地方
//! 一律按 `(rank, 原子下标)` 排。**不排的话同一个分子换个 SMILES 写法就换一组坐标** ——
//! v1 为这件事红过两次(法向量按存储序累加、排序键在平局时退回存储序)。

use crate::geom::{place_nerf, Point3};
use crate::params::{self, Source};
use crate::vsepr::{arrangement, child_torsions};
use omgkit_core::MolBuilder;

/// 一次摆放的**分级计数**。
///
/// 每一条提前返回都要在这儿留下痕迹 —— 没有计数器的分支迟早会咬人
/// (v1 反复栽在静默 `continue` 上)。
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Stats {
    /// 原子总数。
    pub atoms: usize,
    /// 摆好的原子数。
    pub placed: usize,
    /// 因为**在环上**而没摆的(一期的范围之外,不是失败)。
    pub skipped_ring: usize,
    /// 因为参考原子共线、NeRF 定不出标架而没摆的。
    pub degenerate: usize,
    /// 连不上已摆部分的(多片段分子的第二个片段等)。
    pub disconnected: usize,
    /// 键长查表逐项命中的次数。
    pub bond_table: usize,
    /// 键长退到"不在环里"那一行的次数。
    pub bond_relaxed: usize,
    /// 键长只能用共价半径模型的次数。
    pub bond_model: usize,
    /// 键角查表逐项命中 / 放宽 / 兜底。
    pub angle_table: usize,
    /// 见 [`Stats::angle_table`]。
    pub angle_relaxed: usize,
    /// 见 [`Stats::angle_table`]。
    pub angle_model: usize,
    /// 配位数 ≥ 5 的中心个数。**一期明确不保证它们摆得对**,只计数。
    pub degree_ge5: usize,
    /// **兄弟角被推歪超过 5° 的中心个数。**
    ///
    /// 表只给中心一个角 θ,而 4 个取代基有 6 个夹角、模掉整体转动只有 5 个自由度 ——
    /// 超定。构造法让"父–子"那几个精确等于 θ,兄弟之间的是**推出来**的:
    /// 扭转差 120° 时 `cos φ = cos²θ + sin²θ·cos120°`。
    ///
    /// 只有 θ **恰好** 109.4712° 时 φ = θ。表里的 109.4° 差 +0.1423°,可忽略;
    /// 但 θ 离得远就不行:**θ = 120° 时差 −22.8°**、θ = 104.5° 时差 +9.5°。
    /// 所以"排布判成四面体、而表角离 109.47° 很远"的中心必须计数 ——
    /// 那种中心的兄弟角是错的,而键长判据看不见、父–子那几个角也照样精确。
    pub angle_strained: usize,
}

impl Stats {
    fn note_bond(&mut self, s: Source) {
        match s {
            Source::Table => self.bond_table += 1,
            Source::RingRelaxed => self.bond_relaxed += 1,
            Source::Model => self.bond_model += 1,
        }
    }
    fn note_angle(&mut self, s: Source) {
        match s {
            Source::Table => self.angle_table += 1,
            Source::RingRelaxed => self.angle_relaxed += 1,
            Source::Model => self.angle_model += 1,
        }
    }
}

/// 摆放结果。
#[derive(Debug, Clone)]
pub struct Placed {
    /// 逐原子坐标。没摆好的那些是 [`Point3::ORIGIN`],**必须看 [`Placed::placed`]**,
    /// 不能拿坐标本身判有没有摆好。
    pub coords: Vec<Point3>,
    /// 逐原子:摆好了没有。
    pub placed: Vec<bool>,
    /// 分级计数。
    pub stats: Stats,
}

impl Placed {
    /// 全部原子都摆好了。
    #[must_use]
    pub fn complete(&self) -> bool {
        self.stats.atoms == self.stats.placed
    }
}

/// 按 `(rank, 下标)` 排好的邻居 —— **所有挑选都走这里**,免得漏掉某处按存储序。
fn sorted_neighbors(mol: &MolBuilder, a: u32, ranks: &[u32]) -> Vec<u32> {
    let mut v: Vec<u32> = mol.neighbors(a).map(|(y, _)| y).collect();
    v.sort_by_key(|y| (ranks[*y as usize], *y));
    v
}

/// 两个原子之间那根键的键级与"所在最小环"。一期无环,环尺寸恒 0。
fn bond_between(mol: &MolBuilder, a: u32, b: u32) -> Option<omgkit_core::BondOrder> {
    mol.neighbors(a)
        .find(|(y, _)| *y == b)
        .map(|(_, bi)| mol.bonds()[bi as usize].order)
}

/// 四面体排布下,表角 `θ` 推出来的兄弟角与 `θ` 差多少(弧度)。见 [`Stats::angle_strained`]。
fn sibling_skew(theta: f64) -> f64 {
    let phi = theta
        .cos()
        .mul_add(theta.cos(), theta.sin().powi(2) * -0.5)
        .clamp(-1.0, 1.0)
        .acos();
    (phi - theta).abs()
}

/// 摆一个**无环**分子。
///
/// `ranks` 要是规范秩(见模块文档)。有环的分子这一期不摆,
/// 会把环上的原子记进 [`Stats::skipped_ring`] 并留在 `placed = false`。
///
/// **永不 panic**:任何摆不了的情形都落进计数器,坐标留 `ORIGIN` 且 `placed = false`。
#[must_use]
pub fn place(mol: &MolBuilder, ranks: &[u32]) -> Placed {
    let n = mol.num_atoms();
    let mut coords = vec![Point3::ORIGIN; n];
    let mut placed = vec![false; n];
    let mut parent: Vec<Option<u32>> = vec![None; n];
    let mut stats = Stats {
        atoms: n,
        ..Stats::default()
    };
    if n == 0 || ranks.len() < n {
        return Placed {
            coords,
            placed,
            stats,
        };
    }

    // 环上的原子一期不摆
    let rings = omgkit_chem::sssr::ring_set(mol);
    let mut on_ring = vec![false; n];
    for r in &rings {
        for a in &r.atoms {
            on_ring[*a as usize] = true;
        }
    }
    stats.skipped_ring = on_ring.iter().filter(|x| **x).count();

    for a in 0..n {
        #[allow(clippy::cast_possible_truncation)]
        if mol.neighbors(a as u32).count() >= 5 {
            stats.degree_ge5 += 1;
        }
    }

    // 根:非环原子里 (rank, 下标) 最小的那个
    let root = (0..n)
        .filter(|a| !on_ring[*a])
        .min_by_key(|a| (ranks[*a], *a));
    let Some(root) = root else {
        return Placed {
            coords,
            placed,
            stats,
        };
    };
    #[allow(clippy::cast_possible_truncation)]
    let root = root as u32;

    // ---- 起手:根放原点,第一个邻居沿 +x,其余绕 root–n1 轴铺开 ----
    placed[root as usize] = true;
    let nbrs: Vec<u32> = sorted_neighbors(mol, root, ranks)
        .into_iter()
        .filter(|x| !on_ring[*x as usize])
        .collect();
    let mut queue: std::collections::VecDeque<u32> = std::collections::VecDeque::new();
    if let Some(&n1) = nbrs.first() {
        let ord = bond_between(mol, root, n1).unwrap_or(omgkit_core::BondOrder::Single);
        let b = params::bond_length(
            mol.atoms()[root as usize].atomic_num,
            mol.atoms()[n1 as usize].atomic_num,
            ord,
            0,
        );
        stats.note_bond(b.source);
        coords[n1 as usize] = Point3::new(b.value, 0.0, 0.0);
        placed[n1 as usize] = true;
        parent[n1 as usize] = Some(root);
        queue.push_back(n1);

        // 根上剩下的邻居:以 n1 为角的另一条边,绕 root–n1 轴按扭转角铺开。
        // 参考点 `v` 只用来定扭转角的零点(整体转动本来就是自由的)。
        let v = Point3::new(0.0, 1.0, 0.0);
        let deg = mol.neighbors(root).count();
        let arr = arrangement(mol.atoms()[root as usize].hybridization, deg);
        let ts = child_torsions(arr, nbrs.len().saturating_sub(1));
        for (k, &x) in nbrs.iter().skip(1).enumerate() {
            let Some(&t) = ts.get(k) else {
                stats.degenerate += 1;
                continue;
            };
            let ord = bond_between(mol, root, x).unwrap_or(omgkit_core::BondOrder::Single);
            let bl = params::bond_length(
                mol.atoms()[root as usize].atomic_num,
                mol.atoms()[x as usize].atomic_num,
                ord,
                0,
            );
            stats.note_bond(bl.source);
            let ang = params::angle(
                mol.atoms()[root as usize].atomic_num,
                deg,
                mol.atoms()[root as usize]
                    .flags
                    .contains(omgkit_core::AtomFlags::AROMATIC),
                0,
                0,
            );
            stats.note_angle(ang.source);
            if k == 0
                && arr == crate::vsepr::Arrangement::Tetrahedral
                && sibling_skew(ang.value) > 5f64.to_radians()
            {
                stats.angle_strained += 1;
            }
            match place_nerf(
                v,
                coords[n1 as usize],
                coords[root as usize],
                bl.value,
                ang.value,
                t,
            ) {
                Some(p) => {
                    coords[x as usize] = p;
                    placed[x as usize] = true;
                    parent[x as usize] = Some(root);
                    queue.push_back(x);
                }
                None => stats.degenerate += 1,
            }
        }
    }

    // ---- BFS:每个已摆好的中心,把它没摆的邻居一次摆完 ----
    while let Some(c) = queue.pop_front() {
        let p = parent[c as usize];
        let kids: Vec<u32> = sorted_neighbors(mol, c, ranks)
            .into_iter()
            .filter(|x| !placed[*x as usize] && !on_ring[*x as usize])
            .collect();
        if kids.is_empty() {
            continue;
        }
        let Some(p) = p else {
            stats.degenerate += kids.len();
            continue;
        };
        // 祖父:`p` 身上另一个已摆好的邻居;没有就用一个虚点定扭转零点
        let g = sorted_neighbors(mol, p, ranks)
            .into_iter()
            .find(|y| *y != c && placed[*y as usize])
            .map_or_else(
                || coords[p as usize] + Point3::new(0.0, 1.0, 0.0),
                |y| coords[y as usize],
            );
        let deg = mol.neighbors(c).count();
        let arr = arrangement(mol.atoms()[c as usize].hybridization, deg);
        let ts = child_torsions(arr, kids.len());
        for (k, &x) in kids.iter().enumerate() {
            let Some(&t) = ts.get(k) else {
                stats.degenerate += 1;
                continue;
            };
            let ord = bond_between(mol, c, x).unwrap_or(omgkit_core::BondOrder::Single);
            let bl = params::bond_length(
                mol.atoms()[c as usize].atomic_num,
                mol.atoms()[x as usize].atomic_num,
                ord,
                0,
            );
            stats.note_bond(bl.source);
            let ang = params::angle(
                mol.atoms()[c as usize].atomic_num,
                deg,
                mol.atoms()[c as usize]
                    .flags
                    .contains(omgkit_core::AtomFlags::AROMATIC),
                0,
                0,
            );
            stats.note_angle(ang.source);
            if k == 0
                && arr == crate::vsepr::Arrangement::Tetrahedral
                && sibling_skew(ang.value) > 5f64.to_radians()
            {
                stats.angle_strained += 1;
            }
            match place_nerf(
                g,
                coords[p as usize],
                coords[c as usize],
                bl.value,
                ang.value,
                t,
            ) {
                Some(q) => {
                    coords[x as usize] = q;
                    placed[x as usize] = true;
                    parent[x as usize] = Some(c);
                    queue.push_back(x);
                }
                None => stats.degenerate += 1,
            }
        }
    }

    stats.placed = placed.iter().filter(|x| **x).count();
    // 没摆好、又不在环上、又不是退化的 —— 那就是连不上(多片段)
    stats.disconnected = stats
        .atoms
        .saturating_sub(stats.placed + stats.skipped_ring + stats.degenerate);
    Placed {
        coords,
        placed,
        stats,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::{angle_at, dihedral};

    /// 从 SMILES 造一个补好氢、感知过的分子,外加规范秩。
    fn prep(smi: &str) -> (MolBuilder, Vec<u32>) {
        let mut m = omgkit_io::smiles::parse(smi).expect("SMILES 该解析得了");
        omgkit_chem::pipeline::sanitize(&mut m).expect("该 sanitize 得了");
        let r = omgkit_io::canon::classed_ranks(&m);
        omgkit_chem::add_explicit_hs(&mut m, &r);
        let r = omgkit_io::canon::classed_ranks(&m);
        (m, r)
    }

    /// `sibling_skew` 的数必须与解析式对得上 —— [`Stats::angle_strained`] 全靠它。
    ///
    /// 这几个数是算出来的:`cos φ = cos²θ + sin²θ·cos120°`。
    #[test]
    fn sibling_skew_matches_the_analytic_values() {
        for (deg, want) in [
            (109.4712_f64, 0.0_f64), // 真四面体角:唯一自洽的点
            (109.4, 0.1423),         // 表里的值:可忽略
            (120.0, 22.8192),        // 远离四面体角:**差 22.8°**
            (104.5, 9.4516),
        ] {
            let got = sibling_skew(deg.to_radians()).to_degrees();
            assert!(
                (got - want).abs() < 1e-3,
                "θ={deg}° 推出来差 {got:.4}°,解析式说 {want}°"
            );
        }
        // 门槛 5° 必须把 120° 那种拦下、把 109.4° 那种放过
        assert!(sibling_skew(120f64.to_radians()) > 5f64.to_radians());
        assert!(sibling_skew(109.4f64.to_radians()) < 5f64.to_radians());
    }

    /// 无环分子必须**全部摆好**,而且每根键的长度就是表里那个值。
    ///
    /// 键长是构造法唯一"按定义精确"的量 —— 它要是不准,说明 NeRF 或查表接错了。
    #[test]
    fn every_acyclic_molecule_is_fully_placed_with_exact_bond_lengths() {
        let mut checked = 0;
        for smi in [
            "C",
            "CC",
            "CCC",
            "CCO",
            "CC(C)C",
            "CC(=O)O",
            "CCN(CC)CC",
            "C#N",
            "CC(C)(C)CC(C)(C)C",
            "OCC(O)CO",
            "CC(=O)OC",
            "FC(F)(F)C(Cl)Br",
        ] {
            let (m, r) = prep(smi);
            let out = place(&m, &r);
            assert!(out.complete(), "{smi} 没摆全:{:?}", out.stats);
            assert_eq!(out.stats.skipped_ring, 0, "{smi} 不该有环");
            for (bi, b) in m.bonds().iter().enumerate() {
                let want = params::bond_length(
                    m.atoms()[b.begin as usize].atomic_num,
                    m.atoms()[b.end as usize].atomic_num,
                    b.order,
                    0,
                )
                .value;
                let got = out.coords[b.begin as usize].dist(out.coords[b.end as usize]);
                assert!(
                    ((got - want) / want).abs() < 1e-12,
                    "{smi} 第 {bi} 根键 {got:.6} 该是 {want:.6}"
                );
                checked += 1;
            }
        }
        assert!(checked > 100, "只验了 {checked} 根键");
    }

    /// **有环的分子这一期不摆,但要如实报**:环上的原子计进 `skipped_ring`,
    /// 而且**永不 panic**。
    #[test]
    fn a_ring_is_reported_not_mangled_and_never_panics() {
        for smi in [
            "c1ccccc1",
            "C1CCCCC1",
            "c1ccncc1",
            "C1CC1",
            "c1ccc2ccccc2c1",
        ] {
            let (m, r) = prep(smi);
            let out = place(&m, &r);
            assert!(out.stats.skipped_ring > 0, "{smi} 该报有环");
            assert!(!out.complete(), "{smi} 一期不该摆全");
            // 摆好的那些坐标必须是有限数
            for (i, p) in out.coords.iter().enumerate() {
                if out.placed[i] {
                    assert!(p.is_finite(), "{smi} 第 {i} 个原子坐标不是有限数");
                }
            }
        }
    }

    /// **sp³ 骨架该是交错的**:丁烷主链的二面角要么 180° 要么 ±60°,
    /// 不能是 0°(全重叠,那是能量最高的构象)。
    #[test]
    fn an_sp3_backbone_comes_out_staggered_not_eclipsed() {
        let (m, r) = prep("CCCC");
        let out = place(&m, &r);
        assert!(out.complete());
        // 找四个碳
        let c: Vec<u32> = (0..m.num_atoms())
            .filter(|i| m.atoms()[*i].atomic_num == 6)
            .map(|i| u32::try_from(i).unwrap())
            .collect();
        assert_eq!(c.len(), 4);
        // 主链顺序:按连接找一条路径
        let mut chain = vec![c[0]];
        while chain.len() < 4 {
            let last = *chain.last().unwrap();
            let nxt = m
                .neighbors(last)
                .map(|(y, _)| y)
                .find(|y| c.contains(y) && !chain.contains(y));
            match nxt {
                Some(y) => chain.push(y),
                None => break,
            }
        }
        if chain.len() == 4 {
            let d = dihedral(
                out.coords[chain[0] as usize],
                out.coords[chain[1] as usize],
                out.coords[chain[2] as usize],
                out.coords[chain[3] as usize],
            )
            .expect("二面角")
            .to_degrees()
            .abs();
            assert!(d > 45.0, "主链二面角 {d:.1}° —— 太接近重叠构象");
        }
    }

    /// **sp² 中心必须共面。** 乙烯的四个氢与两个碳该在一个平面上。
    #[test]
    fn an_sp2_centre_stays_flat_in_a_real_molecule() {
        let (m, r) = prep("C=C");
        let out = place(&m, &r);
        assert!(out.complete(), "{:?}", out.stats);
        let c: Vec<u32> = (0..m.num_atoms())
            .filter(|i| m.atoms()[*i].atomic_num == 6)
            .map(|i| u32::try_from(i).unwrap())
            .collect();
        assert_eq!(c.len(), 2);
        let hs: Vec<u32> = m
            .neighbors(c[0])
            .map(|(y, _)| y)
            .filter(|y| *y != c[1])
            .collect();
        assert_eq!(hs.len(), 2, "乙烯每个碳该有两个氢");
        // H–C=C–H 的二面角必须是 0 或 180
        for &h in &hs {
            let h2 = m
                .neighbors(c[1])
                .map(|(y, _)| y)
                .find(|y| *y != c[0])
                .expect("另一端的氢");
            let d = dihedral(
                out.coords[h as usize],
                out.coords[c[0] as usize],
                out.coords[c[1] as usize],
                out.coords[h2 as usize],
            )
            .expect("二面角")
            .to_degrees()
            .abs();
            assert!(
                d < 1e-6 || (d - 180.0).abs() < 1e-6,
                "乙烯离面了,H–C=C–H = {d:.4}°"
            );
        }
    }

    /// 键角必须是表里那个值(无环、非芳香的中心)。
    #[test]
    fn bond_angles_match_the_table() {
        let (m, r) = prep("CCC");
        let out = place(&m, &r);
        assert!(out.complete());
        let mid = (0..m.num_atoms())
            .find(|i| {
                m.atoms()[*i].atomic_num == 6
                    && m.neighbors(u32::try_from(*i).unwrap()).count() == 4
            })
            .map(|i| u32::try_from(i).unwrap());
        let mid = mid.expect("丙烷中间那个碳");
        let want = params::angle(6, 4, false, 0, 0).value.to_degrees();
        let nb: Vec<u32> = m.neighbors(mid).map(|(y, _)| y).collect();
        let mut n = 0;
        let mut worst = 0.0f64;
        for i in 0..nb.len() {
            for j in (i + 1)..nb.len() {
                let got = angle_at(
                    out.coords[nb[i] as usize],
                    out.coords[mid as usize],
                    out.coords[nb[j] as usize],
                )
                .expect("角")
                .to_degrees();
                worst = worst.max((got - want).abs());
                n += 1;
            }
        }
        assert_eq!(n, 6, "四配位中心该有六对");
        // **六个角不可能都等于表里那一个值。** 4 个取代基有 6 个夹角、
        // 而模掉整体转动只有 5 个自由度 —— 超定。构造法让"父–子"那几个
        // 精确等于 θ,兄弟之间的是**推出来**的:`cos φ = cos²θ + sin²θ·cos120°`。
        // θ = 109.4° 时 φ = 109.5423°,差 **+0.1423°**(解析式与实测逐位吻合);
        // 只有 θ 恰好 109.4712° 时才为 0。
        //
        // 头一版这条判据要求六个全等于表值,**是我写错了期望**,不是代码错。
        assert!(
            worst < 0.15,
            "最大偏差 {worst:.4}° —— 该只有兄弟角那 0.1423°"
        );
        assert!(
            worst > 0.10,
            "偏差 {worst:.4}° 太小了 —— 兄弟角那 0.1423° 哪去了?"
        );
    }

    /// **同一分子换个 SMILES 写法,判决不能变。**
    ///
    /// 不要求坐标逐位相同(那是 v1 付过大代价的),但**几何质量**必须一样:
    /// 摆好几个、键长最大误差、最小非键距离 —— 这几个数一位都不许动。
    /// 存储序泄漏、排序键撞车这类 bug 改变的正是它们。
    #[test]
    fn the_verdict_does_not_depend_on_how_the_molecule_is_written() {
        for group in [
            vec!["CCO", "OCC", "C(C)O"],
            vec!["CC(C)C", "C(C)(C)C"],
            vec!["CC(=O)O", "OC(C)=O", "O=C(O)C"],
        ] {
            let mut seen: Option<(usize, u64, u64)> = None;
            for smi in &group {
                let (m, r) = prep(smi);
                let out = place(&m, &r);
                // 键长最大相对误差
                let mut worst = 0.0f64;
                for b in m.bonds() {
                    let want = params::bond_length(
                        m.atoms()[b.begin as usize].atomic_num,
                        m.atoms()[b.end as usize].atomic_num,
                        b.order,
                        0,
                    )
                    .value;
                    let got = out.coords[b.begin as usize].dist(out.coords[b.end as usize]);
                    worst = worst.max(((got - want) / want).abs());
                }
                // 最小非键距离
                let mut mind = f64::INFINITY;
                for i in 0..m.num_atoms() {
                    for j in (i + 1)..m.num_atoms() {
                        let (iu, ju) = (u32::try_from(i).unwrap(), u32::try_from(j).unwrap());
                        if m.neighbors(iu).any(|(y, _)| y == ju) {
                            continue;
                        }
                        mind = mind.min(out.coords[i].dist(out.coords[j]));
                    }
                }
                let key = (out.stats.placed, worst.to_bits(), mind.to_bits());
                match seen {
                    None => seen = Some(key),
                    Some(k) => assert_eq!(
                        key, k,
                        "{smi} 的判决与同组第一个写法不同(摆好数/键长最大误差/最小非键距离)"
                    ),
                }
            }
        }
    }
}
