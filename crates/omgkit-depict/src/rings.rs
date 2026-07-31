//! 环系统布局 —— 整个算法里最难的一段。
//!
//! # 分三档,而且**第三档如实承认自己是退化解**
//!
//! | 形状 | 做法 |
//! |---|---|
//! | 单环 | 正多边形 |
//! | 邻稠(逐个环只与已放置部分共用**一根键**) | 沿那根键把新多边形拼到外侧 |
//! | 桥环 / 笼状 | 规则给不出好解 —— 退化到弹簧松弛,并记进 [`Degradation`] |
//!
//! 第三档是所有工具箱共同的软肋(见 Mayfield, RDKit UGM 2016:桥环与拥挤小环
//! 是 11 类障碍中反复失手的两类)。**这里不假装它成功了**:退化的地方明确
//! 报出来,下游可以选择拒绝渲染或人工介入。悄悄给一张看着还行、其实构型
//! 读不出来的图,比明说"这一块我画不好"糟得多。
//!
//! # 平局一律按规范秩打破
//!
//! 从哪个环起手、共用键取哪一根 —— 这些选择只要沾上原子的**存储下标**,同一个
//! 分子换一种 SMILES 写法就会得到另一张图。全部改用规范秩,写法就影响不到结果。

use std::collections::{BTreeMap, BTreeSet};

use omgkit_chem::sssr::Ring;
use omgkit_core::MolBuilder;

use crate::geom::{regular_polygon, Point2, BOND_LEN};

/// 布局中不得不退化的地方。
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Degradation {
    /// 桥环或笼状体系:没有规则能给出平面上的好解,坐标由弹簧松弛得到。
    ///
    /// 环内键角、键长都不再保证,重叠也可能消不掉。
    BridgedRingSystem {
        /// 涉及的原子
        atoms: Vec<u32>,
    },
}

/// 一个环系统连同落在它里面的 SSSR 环。
pub(crate) struct System<'a> {
    pub atoms: Vec<u32>,
    pub rings: Vec<&'a Ring>,
}

/// 把 SSSR 环按所属的稠环系统归类。
///
/// `fused_ring_systems` 用的是双连通分解,所以**螺环与单键相连的环各成一个
/// 系统**(实测:螺[4.4]壬烷给出两个系统,共用那个螺原子;联苯给出两个系统,
/// 中间那根键不自成系统)。这对布局是好事 —— 各自摆好再接起来即可。
pub(crate) fn group<'a>(systems: &[Vec<u32>], rings: &'a [Ring]) -> Vec<System<'a>> {
    systems
        .iter()
        .map(|atoms| {
            let set: BTreeSet<u32> = atoms.iter().copied().collect();
            System {
                atoms: atoms.clone(),
                rings: rings
                    .iter()
                    .filter(|r| r.atoms.iter().all(|a| set.contains(a)))
                    .collect(),
            }
        })
        .collect()
}

/// 在**局部坐标系**里给一个环系统布局。返回逐原子坐标与(可能的)退化记录。
///
/// 调用方拿到之后再整体平移旋转到该去的位置,见 [`place_at`]。
pub(crate) fn layout_local(
    mol: &MolBuilder,
    sys: &System<'_>,
    ranks: &[u32],
) -> (BTreeMap<u32, Point2>, Option<Degradation>) {
    let mut pos: BTreeMap<u32, Point2> = BTreeMap::new();

    if sys.rings.is_empty() {
        // 环感知说这里有环、SSSR 却一个都没给出来。不猜,直接走退化。
        return (relax(mol, &sys.atoms, ranks), Some(bridged(&sys.atoms)));
    }

    // 起手环:先按大小(大环更能定住整体形状),再按规范秩 —— 不看存储下标
    let mut order: Vec<&Ring> = sys.rings.clone();
    order.sort_by_key(|r| (std::cmp::Reverse(r.atoms.len()), ring_key(r, ranks)));
    // **起手环怎么摆,必须完全由规范秩决定。**
    //
    // SSSR 给出的环原子顺序依赖存储序:同一个环,不同写法可能从不同的原子起、
    // 甚至朝相反方向绕。直接拿它去对多边形顶点,两种写法就落在**不同的构型**上
    // (不只是旋转或镜像)—— 后面的取代基方向、消冲突翻哪根键全跟着分岔。
    //
    // 实测:阿司匹林的两种写法,一种消冲突一次没翻,另一种翻了两次,最后坐标
    // 对不上。而"两两距离的多重集"那种指纹**看不出来**,因为点集确实全等。
    let first = canonical_cycle(&order[0].atoms, ranks);
    for (a, p) in first.iter().zip(regular_polygon(first.len(), 0.0)) {
        pos.insert(*a, p);
    }

    let mut placed: BTreeSet<usize> = BTreeSet::from([0]);
    let mut degraded = None;

    // 反复找"只与已放置部分共用一根键"的环,拼上去
    while placed.len() < order.len() {
        let mut best: Option<(usize, u32, u32)> = None;
        for (i, r) in order.iter().enumerate() {
            if placed.contains(&i) {
                continue;
            }
            let shared: Vec<u32> = r
                .atoms
                .iter()
                .copied()
                .filter(|a| pos.contains_key(a))
                .collect();
            if shared.len() != 2 {
                continue; // 共用 1 个是螺(不会同系统)、>2 个是桥
            }
            let (u, v) = (shared[0], shared[1]);
            if !adjacent_in_ring(r, u, v) {
                continue; // 共用两个原子却不相邻 —— 那也是桥
            }
            let key = (ring_key(r, ranks), ranks[u as usize].min(ranks[v as usize]));
            // 不用 `Option::is_none_or` —— 它到 1.82 才稳定,而工作区 MSRV 是 1.75
            let better = match best {
                None => true,
                Some((bi, bu, bv)) => {
                    key < (
                        ring_key(order[bi], ranks),
                        ranks[bu as usize].min(ranks[bv as usize]),
                    )
                }
            };
            if better {
                best = Some((i, u, v));
            }
        }

        let Some((i, u, v)) = best else {
            // 剩下的环都不是邻稠 —— 桥环。整个系统交给松弛。
            // 桥环一经识别就**丢掉部分结果**,整个系统重排 —— 见 relax 的注释
            degraded = Some(bridged(&sys.atoms));
            return (relax(mol, &sys.atoms, ranks), degraded);
        };

        fuse_on_bond(order[i], u, v, &mut pos);
        placed.insert(i);
    }

    (pos, degraded)
}

/// 把一个环的原子序列旋转/翻转到**规范起点与规范方向**。
///
/// 起点取规范秩最小的原子;方向取"沿着走一圈得到的秩序列字典序更小"的那一边。
/// 两个自由度都被定死,同一个环无论怎么写都得到同一个序列。
fn canonical_cycle(atoms: &[u32], ranks: &[u32]) -> Vec<u32> {
    let n = atoms.len();
    let start = (0..n)
        .min_by_key(|i| (ranks[atoms[*i] as usize], atoms[*i]))
        .expect("环非空");
    let fwd: Vec<u32> = (0..n).map(|k| atoms[(start + k) % n]).collect();
    let bwd: Vec<u32> = (0..n).map(|k| atoms[(start + n - k) % n]).collect();
    let key = |v: &[u32]| -> Vec<u32> { v.iter().map(|a| ranks[*a as usize]).collect() };
    if key(&bwd) < key(&fwd) {
        bwd
    } else {
        fwd
    }
}

fn bridged(atoms: &[u32]) -> Degradation {
    let mut a = atoms.to_vec();
    a.sort_unstable();
    Degradation::BridgedRingSystem { atoms: a }
}

/// 环的确定性排序键:环上规范秩的**有序**多重集。
///
/// 用规范秩而不是原子下标,同一分子的不同写法才会选出同一个起手环。
fn ring_key(r: &Ring, ranks: &[u32]) -> Vec<u32> {
    let mut k: Vec<u32> = r.atoms.iter().map(|a| ranks[*a as usize]).collect();
    k.sort_unstable();
    k
}

fn adjacent_in_ring(r: &Ring, u: u32, v: u32) -> bool {
    let n = r.atoms.len();
    (0..n).any(|i| {
        let (a, b) = (r.atoms[i], r.atoms[(i + 1) % n]);
        (a == u && b == v) || (a == v && b == u)
    })
}

/// 沿已放置的键 `u–v` 把环 `r` 拼到外侧。
fn fuse_on_bond(r: &Ring, u: u32, v: u32, pos: &mut BTreeMap<u32, Point2>) {
    let n = r.atoms.len();
    // 把环的原子序列转到以 u 开头、v 紧随其后
    let start = r.atoms.iter().position(|a| *a == u).expect("u 在环上");
    let forward = r.atoms[(start + 1) % n] == v;
    let seq: Vec<u32> = (0..n)
        .map(|k| {
            let i = if forward { start + k } else { start + n - k };
            r.atoms[i % n]
        })
        .collect();
    debug_assert_eq!(seq[0], u);
    debug_assert_eq!(seq[1], v);

    let (pu, pv) = (pos[&u], pos[&v]);
    let mid = (pu + pv) * 0.5;
    let along = (pv - pu).normalized();
    let normal = Point2::new(-along.y, along.x);
    // 边心距:边长 s 的正 n 边形,中心到边的距离是 s / (2 tan(π/n))
    let apothem = BOND_LEN / (2.0 * (std::f64::consts::PI / n as f64).tan());

    // 两个候选中心,取**远离已放置质心**的那个 —— 新环要长在外侧
    let anchor = centroid(pos.values().copied());
    let c1 = mid + normal * apothem;
    let c2 = mid - normal * apothem;
    let center = if c1.dist(anchor) > c2.dist(anchor) {
        c1
    } else {
        c2
    };

    // 转向的正负:取能把 pu 转到 pv 的那一个。这一步顺带验证了几何 ——
    // 边心距或法线写错的话,两个方向都转不到 pv,debug 下会直接断言失败。
    let step = std::f64::consts::TAU / n as f64;
    let sign = if pu.rotated_about(center, step).dist(pv) < 1e-6 {
        1.0
    } else {
        debug_assert!(
            pu.rotated_about(center, -step).dist(pv) < 1e-6,
            "拼环的几何不自洽:两个转向都到不了对面那个原子"
        );
        -1.0
    };

    for (k, a) in seq.iter().enumerate().skip(2) {
        pos.insert(*a, pu.rotated_about(center, sign * step * k as f64));
    }
}

fn centroid(pts: impl Iterator<Item = Point2>) -> Point2 {
    let mut sum = Point2::ORIGIN;
    let mut n = 0.0;
    for p in pts {
        sum = sum + p;
        n += 1.0;
    }
    if n == 0.0 {
        Point2::ORIGIN
    } else {
        sum * (1.0 / n)
    }
}

/// 桥环的兜底:弹簧松弛。
///
/// 键长拉向 [`BOND_LEN`],非键原子互斥。给不出标准键角,也不保证消得掉重叠 ——
/// 这正是调用方要把它记进 [`Degradation`] 的原因。
fn relax(mol: &MolBuilder, atoms: &[u32], ranks: &[u32]) -> BTreeMap<u32, Point2> {
    // **原子按规范秩排序,不按存储下标。** 初值、乃至浮点求和的次序都因此固定,
    // 于是同一个分子的任何写法得到同一张图。
    //
    // 这里刻意**不**接受"贪心走到一半"的部分结果当种子。那个部分结果依赖遍历
    // 顺序,拿它当初值会把写法依赖直接带进来 —— 实测:苊的两种写法就是这样
    // 给出了两个不同形状,而萘、菲、蒽因为太对称,根本触发不到,看着一切正常。
    let mut atoms: Vec<u32> = atoms.to_vec();
    atoms.sort_by_key(|a| (ranks[*a as usize], *a));
    let atoms = &atoms[..];
    let idx: BTreeMap<u32, usize> = atoms.iter().enumerate().map(|(i, a)| (*a, i)).collect();
    let n = atoms.len();

    let r = BOND_LEN * n as f64 / std::f64::consts::TAU.max(1.0);
    let mut p: Vec<Point2> = (0..n)
        .map(|i| Point2::new(r, 0.0).rotated(std::f64::consts::TAU * i as f64 / n as f64))
        .collect();

    let bonded: Vec<(usize, usize)> = mol
        .bonds()
        .iter()
        .filter_map(|b| Some((*idx.get(&b.begin)?, *idx.get(&b.end)?)))
        .collect();

    for _ in 0..400 {
        let mut force = vec![Point2::ORIGIN; n];
        for &(i, j) in &bonded {
            let d = p[j] - p[i];
            let len = d.norm().max(1e-6);
            let f = d.normalized() * ((len - BOND_LEN) * 0.35);
            force[i] = force[i] + f;
            force[j] = force[j] - f;
        }
        for i in 0..n {
            for j in (i + 1)..n {
                let d = p[j] - p[i];
                let len = d.norm().max(1e-6);
                if len < BOND_LEN * 1.2 {
                    let f = d.normalized() * ((BOND_LEN * 1.2 - len) * 0.25);
                    force[i] = force[i] - f;
                    force[j] = force[j] + f;
                }
            }
        }
        for i in 0..n {
            p[i] = p[i] + force[i];
        }
    }

    atoms.iter().copied().zip(p).collect()
}

/// 把一个局部布局整体搬到位:让 `anchor` 落在 `at`,并让整体质心朝 `dir`。
///
/// 两个自由度(平移 + 旋转)刚好被这两个条件定死,不留任意性。
pub(crate) fn place_at(
    local: &BTreeMap<u32, Point2>,
    anchor: u32,
    at: Point2,
    dir: Point2,
) -> BTreeMap<u32, Point2> {
    let a = local[&anchor];
    let c = centroid(local.values().copied());
    let from = (c - a).normalized();
    let to = dir.normalized();
    // from 是零向量只可能出现在"质心恰好落在锚点上"的对称情形,那时转多少都一样
    let theta = if from.norm() < 1e-9 {
        0.0
    } else {
        to.angle() - from.angle()
    };
    local
        .iter()
        .map(|(k, p)| (*k, (*p - a).rotated(theta) + at))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use omgkit_chem::{rings::fused_ring_systems, sssr::ring_set};

    fn prep(smi: &str) -> MolBuilder {
        let mut m = omgkit_io::smiles::parse(smi).unwrap();
        omgkit_chem::pipeline::sanitize(&mut m).unwrap();
        m
    }

    fn layout(smi: &str) -> (BTreeMap<u32, Point2>, Option<Degradation>) {
        let m = prep(smi);
        let ranks = omgkit_io::canon::canonical_ranks(&m);
        let rings = ring_set(&m);
        let sys = group(&fused_ring_systems(&m), &rings);
        layout_local(&m, &sys[0], &ranks)
    }

    fn bond_lengths(m: &MolBuilder, pos: &BTreeMap<u32, Point2>) -> Vec<f64> {
        m.bonds()
            .iter()
            .filter_map(|b| Some(pos.get(&b.begin)?.dist(*pos.get(&b.end)?)))
            .collect()
    }

    #[test]
    fn a_single_ring_is_a_regular_polygon() {
        for smi in ["C1CC1", "C1CCC1", "c1ccccc1", "C1CCCCCC1"] {
            let m = prep(smi);
            let (pos, deg) = layout(smi);
            assert_eq!(deg, None, "{smi} 不该退化");
            assert_eq!(pos.len(), m.num_atoms(), "{smi} 有原子没放上");
            for d in bond_lengths(&m, &pos) {
                assert!((d - BOND_LEN).abs() < 1e-9, "{smi} 键长 {d}");
            }
        }
    }

    #[test]
    fn ortho_fused_rings_share_exactly_one_bond_and_keep_unit_bonds() {
        // 萘、吲哚、芴 —— 逐个环只与已放置部分共用一根键的典型
        for smi in [
            "c1ccc2ccccc2c1",
            "c1ccc2[nH]ccc2c1",
            "c1ccc2c(c1)Cc1ccccc1-2",
        ] {
            let m = prep(smi);
            let (pos, deg) = layout(smi);
            assert_eq!(deg, None, "{smi} 不该退化");
            let ring_atoms: BTreeSet<u32> = pos.keys().copied().collect();
            for b in m.bonds() {
                if ring_atoms.contains(&b.begin) && ring_atoms.contains(&b.end) {
                    let d = pos[&b.begin].dist(pos[&b.end]);
                    assert!((d - BOND_LEN).abs() < 1e-9, "{smi} 环内键长 {d}");
                }
            }
            // 稠环的原子必须两两分开 —— 拼错方向会让新环叠回旧环上,而键长
            // 全都还是 1.0,只看键长发现不了
            let pts: Vec<Point2> = pos.values().copied().collect();
            for i in 0..pts.len() {
                for j in (i + 1)..pts.len() {
                    assert!(pts[i].dist(pts[j]) > 0.5, "{smi} 有两个原子挤在一起");
                }
            }
        }
    }

    #[test]
    fn a_bridged_system_says_so_instead_of_pretending() {
        // 双环[2.2.2]辛烷:三个六元环两两共用不止一根键,平面上没有好解。
        // 判据不是"画得好看",是"**如实说自己画不好**"。
        let (pos, deg) = layout("C1CC2CCC1CC2");
        assert!(
            matches!(deg, Some(Degradation::BridgedRingSystem { .. })),
            "桥环必须报退化,得到 {deg:?}"
        );
        assert_eq!(pos.len(), 8, "退化也要把每个原子都放上");
        for p in pos.values() {
            assert!(p.x.is_finite() && p.y.is_finite(), "退化解不能给出 NaN");
        }
    }

    #[test]
    fn the_layout_does_not_depend_on_how_the_ring_was_written() {
        // **这条测试挑分子要挑对。** 萘、菲、蒽换写法都给出同一形状 —— 但那
        // 不是因为算法写法无关,而是因为它们太对称,起手环已经被"按大小降序"
        // 定死,平局判据压根没被触发。拿它们当判据是走过场(实测:把 ring_key
        // 改回用存储下标,这三个仍然全绿)。
        //
        // 苊是桥式系统,会落到 relax() 兜底,而兜底以"贪心走到哪一步"为初值
        // —— 写法依赖真正藏在那里。用它才守得住。
        let shapes: Vec<Vec<i64>> = ["C1Cc2cccc3cccc1c23", "c1cc2CCc3cccc(c1)c23"]
            .iter()
            .map(|smi| shape_key(smi))
            .collect();
        assert_eq!(shapes[0], shapes[1], "同一分子的两种写法给出了不同形状");
    }

    /// 形状指纹:两两距离排序后的多重集。与原子编号、平移旋转都无关。
    fn shape_key(smi: &str) -> Vec<i64> {
        let m = prep(smi);
        let ranks = omgkit_io::canon::canonical_ranks(&m);
        let rings = ring_set(&m);
        let sys = group(&fused_ring_systems(&m), &rings);
        let s = sys.iter().max_by_key(|s| s.atoms.len()).expect("有环系统");
        let (pos, _) = layout_local(&m, s, &ranks);
        let pts: Vec<Point2> = pos.values().copied().collect();
        let mut ds: Vec<i64> = (0..pts.len())
            .flat_map(|i| ((i + 1)..pts.len()).map(move |j| (i, j)))
            .map(|(i, j)| (pts[i].dist(pts[j]) * 1e4).round() as i64)
            .collect();
        ds.sort_unstable();
        ds
    }
}
