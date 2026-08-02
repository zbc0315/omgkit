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
        return (
            relax(mol, &sys.atoms, ranks, &sys.rings),
            Some(bridged(&sys.atoms)),
        );
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
            return (relax(mol, &sys.atoms, ranks, &sys.rings), degraded);
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
fn relax(mol: &MolBuilder, atoms: &[u32], ranks: &[u32], rings: &[&Ring]) -> BTreeMap<u32, Point2> {
    // **先查表。** 松弛是局部下降,5 个初值本身就常常给出自交的解 —— 实测最常见
    // 的 8 个骨架里 5 个自交,双环[2.2.2]辛烷和金刚烷都在内。表里存的是同一个
    // `quality` 口径下搜得久得多的结果,见 [`crate::templates`]。
    if let Some(p) = crate::templates::lookup(mol, atoms, ranks) {
        return p;
    }
    // **原子按规范秩排序,不按存储下标。** 初值、乃至浮点求和的次序都因此固定,
    // 于是同一个分子的任何写法得到同一张图。
    //
    // 这里刻意**不**接受"贪心走到一半"的部分结果当种子。那个部分结果依赖遍历
    // 顺序,拿它当初值会把写法依赖直接带进来 —— 实测:苊的两种写法就是这样
    // 给出了两个不同形状,而萘、菲、蒽因为太对称,根本触发不到,看着一切正常。
    let mut sorted: Vec<u32> = atoms.to_vec();
    sorted.sort_by_key(|a| (ranks[*a as usize], *a));

    // **多起点。** 松弛是局部下降,落到哪个局部极小全看初值。单一初值下
    // 实测 177 个桥环系统里 172 个(97%)自身有键交叉 —— 那不是消冲突没做好,
    // 消冲突根本够不着:环系统是 2-连通的刚性块,翻转动不了它内部的相对位置。
    //
    // 换几个初值再挑最好的,算法一个字不用改。每个初值都由规范秩派生,
    // 挑选时的平局用量化坐标序列打破,写法无关这条不受影响。
    let mut best: Option<(Quality, BTreeMap<u32, Point2>)> = None;
    for seed in 0..SEEDS {
        let out = relax_from(mol, &sorted, seed, rings, ranks);
        let key = quality(mol, &out, ranks);
        let take = match &best {
            None => true,
            Some((b, _)) => key < *b,
        };
        if take {
            best = Some((key, out));
        }
    }
    best.expect("SEEDS 至少为 1").1
}

/// 试几个初值。**每多一个都要有个说法**,而且要拿全量语料量过。
///
/// 试过第 6 个(多边形起手 + 新原子放进最大空隙):交叉多消 6 处,写法无关
/// 却多 6 处违例。**写法无关是本库的头号契约,不拿它换**,所以没要。
const SEEDS: usize = 5;

/// 一个松弛解好不好:(系统内自交的键对数, 最大键长偏差, 量化坐标序列)。
///
/// **越小越好。**
type Quality = (usize, i64, Vec<(i64, i64)>);

/// 给一个松弛解打分,口径见 [`Quality`]。
///
/// 键长偏差排第二是因为这条路径上"键长全等"本来就不成立(实测全部 177 个
/// 桥环系统 relax 之后偏差都 ≥20%),但同样糟的两个解里该挑偏差小的。
/// 第三项是平局兜底:不留任意性,同一个分子的任何写法挑到同一个解。
fn quality(mol: &MolBuilder, pos: &BTreeMap<u32, Point2>, ranks: &[u32]) -> Quality {
    let live: Vec<(u32, Point2, Point2)> = mol
        .bonds()
        .iter()
        .enumerate()
        .filter_map(|(i, b)| {
            Some((
                u32::try_from(i).ok()?,
                *pos.get(&b.begin)?,
                *pos.get(&b.end)?,
            ))
        })
        .collect();
    let mut cross = 0usize;
    for (k, (_, u1, v1)) in live.iter().enumerate() {
        for (_, u2, v2) in &live[k + 1..] {
            if crate::geom::segments_cross(*u1, *v1, *u2, *v2) {
                cross += 1;
            }
        }
    }
    #[allow(clippy::cast_possible_truncation)]
    let dev = live
        .iter()
        .map(|(_, u, v)| ((u.dist(*v) - BOND_LEN).abs() * 1e6).round() as i64)
        .max()
        .unwrap_or(0);
    // **按规范秩排,不按原子下标。** `BTreeMap` 的迭代序是下标序,拿它当平局
    // 兜底就把写法依赖带了进来 —— 两种写法会在同分的几个解里挑到不同的那个。
    let mut by_rank: Vec<(u32, Point2)> =
        pos.iter().map(|(a, p)| (ranks[*a as usize], *p)).collect();
    by_rank.sort_by_key(|x| x.0);
    #[allow(clippy::cast_possible_truncation)]
    let seq: Vec<(i64, i64)> = by_rank
        .iter()
        .map(|(_, p)| ((p.x * 1e6).round() as i64, (p.y * 1e6).round() as i64))
        .collect();
    (cross, dev, seq)
}

/// 从第 `seed` 个初值出发做一遍松弛。`atoms` 已按规范秩排好。
fn relax_from(
    mol: &MolBuilder,
    atoms: &[u32],
    seed: usize,
    rings: &[&Ring],
    ranks: &[u32],
) -> BTreeMap<u32, Point2> {
    let idx: BTreeMap<u32, usize> = atoms.iter().enumerate().map(|(i, a)| (*a, i)).collect();
    let n = atoms.len();

    let bonded: Vec<(usize, usize)> = mol
        .bonds()
        .iter()
        .filter_map(|b| Some((*idx.get(&b.begin)?, *idx.get(&b.end)?)))
        .collect();

    // 初值 4:**最大的那个环先摆成正多边形**,其余原子沿着已放好的邻居向外
    // 铺开。前四个初值都是"所有原子摆在一个圆上",拓扑上太像,弹簧下降往往
    // 收敛到同一批坏极小;这个起手的形状不一样,实测它才是降幅的主要来源。
    if seed >= 4 {
        if let Some(p) = polygon_seed(mol, atoms, rings, ranks, &idx) {
            return settle(p, n, &bonded, atoms);
        }
    }

    // 其余初值全部由规范秩派生,不看存储下标:
    //   0 圆环、规范秩序      1 圆环、逆序
    //   2 圆环、BFS 序        3 圆环、隔一个取一个(把成键的原子在圆上分开)
    let order: Vec<usize> = match seed {
        1 => (0..n).rev().collect(),
        2 => bfs_order(n, &bonded),
        3 => (0..n).step_by(2).chain((1..n).step_by(2)).collect(),
        _ => (0..n).collect(),
    };
    let r = BOND_LEN * n as f64 / std::f64::consts::TAU.max(1.0);
    let mut p: Vec<Point2> = vec![Point2::ORIGIN; n];
    for (slot, &i) in order.iter().enumerate() {
        p[i] = Point2::new(r, 0.0).rotated(std::f64::consts::TAU * slot as f64 / n as f64);
    }

    settle(p, n, &bonded, atoms)
}

/// 弹簧松弛本体:键长拉到 1,靠得太近的推开。400 步。
fn settle(
    mut p: Vec<Point2>,
    n: usize,
    bonded: &[(usize, usize)],
    atoms: &[u32],
) -> BTreeMap<u32, Point2> {
    for _ in 0..400 {
        let mut force = vec![Point2::ORIGIN; n];
        for &(i, j) in bonded {
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

/// 初值:系统里最大的那个环摆成正多边形,其余原子沿已放好的邻居向外铺开。
///
/// 放不出来(系统里一个环都没有)时返回 `None`,退回圆环初值。
fn polygon_seed(
    mol: &MolBuilder,
    atoms: &[u32],
    rings: &[&Ring],
    ranks: &[u32],
    idx: &BTreeMap<u32, usize>,
) -> Option<Vec<Point2>> {
    // 起手环:先按大小,再按规范秩 —— 平局不许看存储下标
    let first = rings
        .iter()
        .filter(|r| r.atoms.iter().all(|a| idx.contains_key(a)))
        .min_by_key(|r| (std::cmp::Reverse(r.atoms.len()), ring_key(r, ranks)))?;
    let cyc = canonical_cycle(&first.atoms, ranks);

    let n = atoms.len();
    let mut p = vec![Point2::ORIGIN; n];
    let mut placed = vec![false; n];
    for (a, q) in cyc.iter().zip(regular_polygon(cyc.len(), 0.0)) {
        let i = *idx.get(a)?;
        p[i] = q;
        placed[i] = true;
    }

    // 其余的沿已放好的邻居向外铺:方向取"背离已放好那堆的质心"
    loop {
        let next = atoms.iter().enumerate().find(|(i, a)| {
            !placed[*i]
                && mol
                    .neighbors(**a)
                    .any(|(nb, _)| idx.get(&nb).is_some_and(|j| placed[*j]))
        });
        let Some((i, a)) = next else { break };
        // **锚点按规范秩挑,不看 `neighbors` 的存储序。** 有两个已放好的邻居
        // 可选时,拿存储序挑就把写法依赖直接带了进来 —— 实测全量语料的写法
        // 无关违例会从 129 涨到 349。
        let anchor = mol
            .neighbors(*a)
            .filter_map(|(nb, _)| Some((ranks[nb as usize], nb, *idx.get(&nb)?)))
            .filter(|(_, _, j)| placed[*j])
            .min()
            .map(|(_, _, j)| j)?;
        // 背离已放好那堆的质心
        let dir = {
            let (mut c, mut k) = (Point2::ORIGIN, 0.0_f64);
            for (j, on) in placed.iter().enumerate() {
                if *on {
                    c = c + p[j];
                    k += 1.0;
                }
            }
            let away = (p[anchor] - c * (1.0 / k.max(1.0))).normalized();
            if away.norm() < 1e-9 {
                0.0
            } else {
                away.angle()
            }
        };
        p[i] = p[anchor] + Point2::new(BOND_LEN, 0.0).rotated(dir);
        placed[i] = true;
    }
    // 还有没连上的(理论上不会 —— 环系统是连通的),摊在圆上兜底
    for (i, on) in placed.iter().enumerate() {
        if !on {
            p[i] = Point2::new(BOND_LEN * n as f64, 0.0)
                .rotated(std::f64::consts::TAU * i as f64 / n as f64);
        }
    }
    Some(p)
}

/// 从 0 号原子出发的 BFS 序。邻接表按下标升序,而下标已经是规范秩序,
/// 所以这个序也与写法无关。
fn bfs_order(n: usize, bonded: &[(usize, usize)]) -> Vec<usize> {
    let mut adj = vec![Vec::new(); n];
    for &(i, j) in bonded {
        adj[i].push(j);
        adj[j].push(i);
    }
    for a in &mut adj {
        a.sort_unstable();
    }
    let mut seen = vec![false; n];
    let mut out = Vec::with_capacity(n);
    for start in 0..n {
        if seen[start] {
            continue;
        }
        seen[start] = true;
        let mut q = std::collections::VecDeque::from([start]);
        while let Some(x) = q.pop_front() {
            out.push(x);
            for &y in &adj[x] {
                if !seen[y] {
                    seen[y] = true;
                    q.push_back(y);
                }
            }
        }
    }
    out
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
    #[test]
    fn a_bridged_system_is_relaxed_from_several_starts() {
        // 松弛是局部下降,落到哪个局部极小全看初值。单一初值下实测 177 个桥环
        // 系统里 172 个自身有键交叉 —— 而**消冲突根本够不着**:环系统是 2-连通
        // 的刚性块,翻转只动挂在外面的子树,动不了它内部的相对位置。
        //
        // 这条要求多起点确实起作用:把候选砍到一个,下面这些分子的桥环系统
        // 内部就会出现自交。
        let mut won = 0;
        for smi in [
            "CC1(C)[C@@H]2CC[C@@]1(C)C(=O)C2",                          // 樟脑
            "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5", // 吗啡
            "CN1[C@H]2CC[C@@H]1C[C@@H](C2)OC(=O)C(CO)c1ccccc1",         // 阿托品
        ] {
            let mut m = omgkit_io::smiles::parse(smi).unwrap();
            omgkit_chem::pipeline::sanitize(&mut m).unwrap();
            let ranks = omgkit_io::canon::canonical_ranks(&m);
            let systems = omgkit_chem::rings::fused_ring_systems(&m);
            let rs = omgkit_chem::sssr::ring_set(&m);
            let mut checked = 0;
            for sys in group(&systems, &rs) {
                if sys.rings.is_empty() || sys.atoms.len() < 6 {
                    continue;
                }
                let (pos, deg) = layout_local(&m, &sys, &ranks);
                // **只算真正走了松弛那条路的系统。** 邻稠系统走的是正多边形
                // 拼接,拿它去和强行松弛比,当然赢 —— 那样这条判据就是空过的。
                if deg.is_none() {
                    continue;
                }
                let single = relax_from(
                    &m,
                    &{
                        let mut a = sys.atoms.clone();
                        a.sort_by_key(|x| (ranks[*x as usize], *x));
                        a
                    },
                    0,
                    &sys.rings,
                    &ranks,
                );
                let (best, _, _) = quality(&m, &pos, &ranks);
                let (one, _, _) = quality(&m, &single, &ranks);
                assert!(
                    best <= one,
                    "{smi}:多起点挑出来的解({best} 处自交)还不如单起点({one} 处)"
                );
                // `best <= one` 是恒真的(初值 0 本来就在候选里),单靠它这条
                // 判据是空过的。真正要守的是**多起点确实赢过单起点**。
                if best < one {
                    won += 1;
                }
                checked += 1;
            }
            assert!(checked >= 1, "{smi}:一个环系统都没查到");
        }
        assert!(
            won >= 1,
            "多起点在这三个桥环分子上一次都没赢过单起点 —— 那它就是白跑的"
        );
    }

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

/// 桥环骨架坐标表的**生成器**。平时不跑。
///
/// ```shell
/// cargo test -p omgkit-depict --release --lib -- --ignored regenerate_templates --nocapture
/// ```
///
/// 把输出贴进 [`crate::templates`]。与 `harness/gen_elements.py` 生成
/// `element_data.rs` 是同一个路子:**生成脚本进版本库,产物也进版本库**,
/// 谁都能重跑一遍核对。
#[cfg(test)]
mod generator {
    use super::*;
    use crate::geom::Point2;

    const TOP: usize = 24;
    const TRIES: usize = 20_000;

    fn splitmix(state: &mut u64) -> u64 {
        *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = *state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    #[test]
    #[ignore]
    fn regenerate_templates() {
        // 一、扫语料,按出现次数排出最常见的桥环骨架
        let text = std::fs::read_to_string("../../harness/corpus/large.smi").unwrap();
        let mut freq: BTreeMap<String, usize> = BTreeMap::new();
        for line in text.lines() {
            let smi = line.split_whitespace().next().unwrap_or("");
            if smi.is_empty() {
                continue;
            }
            let Ok(mut m) = omgkit_io::smiles::parse(smi) else {
                continue;
            };
            if omgkit_chem::pipeline::sanitize(&mut m).is_err() {
                continue;
            }
            if m.num_atoms() < 2 {
                continue;
            }
            let ranks = omgkit_io::canon::canonical_ranks(&m);
            let rs = omgkit_chem::sssr::ring_set(&m);
            for sys in group(&omgkit_chem::rings::fused_ring_systems(&m), &rs) {
                let (_, deg) = layout_local(&m, &sys, &ranks);
                if deg.is_none() {
                    continue;
                }
                if let Some(k) = crate::templates::skeleton_of(&m, &sys.atoms, &ranks) {
                    *freq.entry(k).or_default() += 1;
                }
            }
        }
        let mut v: Vec<(String, usize)> = freq.into_iter().collect();
        v.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));

        // 二、对每个骨架跑带扰动的多起点,按现成的 Quality 挑最好的
        println!("// 本表由 rings.rs 的 `regenerate_templates` 生成,勿手改。");
        println!("pub(crate) const TABLE: &[(&str, &[(f64, f64)])] = &[");
        let mut kept = 0usize;
        for (skel, n) in v.iter().take(TOP) {
            let Ok(mut m) = omgkit_io::smiles::parse(skel) else {
                continue;
            };
            if omgkit_chem::pipeline::sanitize(&mut m).is_err() {
                continue;
            }
            let ranks = omgkit_io::canon::canonical_ranks(&m);
            let atoms: Vec<u32> = (0..u32::try_from(m.num_atoms()).unwrap()).collect();
            let mut sorted = atoms.clone();
            sorted.sort_by_key(|a| (ranks[*a as usize], *a));
            let cnt = sorted.len();
            let idx: BTreeMap<u32, usize> =
                sorted.iter().enumerate().map(|(i, a)| (*a, i)).collect();
            let bonded: Vec<(usize, usize)> = m
                .bonds()
                .iter()
                .filter_map(|b| Some((*idx.get(&b.begin)?, *idx.get(&b.end)?)))
                .collect();

            let rs = omgkit_chem::sssr::ring_set(&m);
            let sys = group(&omgkit_chem::rings::fused_ring_systems(&m), &rs);
            let Some(s) = sys.iter().max_by_key(|s| s.atoms.len()) else {
                continue;
            };
            let base = relax(&m, &s.atoms, &ranks, &s.rings);
            let mut best = (quality(&m, &base, &ranks), base);

            let mut st = 0x51ED_270B_D5AB_C0DEu64 ^ (cnt as u64);
            let r = BOND_LEN * cnt as f64 / std::f64::consts::TAU;
            for _ in 0..TRIES {
                let mut p = vec![Point2::ORIGIN; cnt];
                for (i, q) in p.iter_mut().enumerate() {
                    let j = (splitmix(&mut st) % 1000) as f64 / 1000.0 - 0.5;
                    let t = std::f64::consts::TAU * (i as f64 + j * 3.0) / cnt as f64;
                    let rad = r * (1.0 + ((splitmix(&mut st) % 1000) as f64 / 1000.0 - 0.5) * 0.6);
                    *q = Point2::new(rad, 0.0).rotated(t);
                }
                let out = settle(p, cnt, &bonded, &sorted);
                let q = quality(&m, &out, &ranks);
                if q < best.0 {
                    best = (q, out);
                }
            }
            // 按**骨架自己的规范秩**存坐标,查表时才对得上
            let mut by_rank: Vec<(u32, Point2)> = best
                .1
                .iter()
                .map(|(a, p)| (ranks[*a as usize], *p))
                .collect();
            by_rank.sort_by_key(|x| x.0);
            print!("    (\"{skel}\", &[");
            for (_, p) in &by_rank {
                print!("({:.6}, {:.6}), ", p.x, p.y);
            }
            println!("]),   // 出现 {n} 次,自交 {}", best.0 .0);
            kept += 1;
        }
        println!("];");
        println!("// 共 {kept} 条");
    }
}
