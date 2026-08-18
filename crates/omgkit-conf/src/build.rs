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
use crate::vsepr::{arrangement, child_torsions, expected_sibling_skew};

/// 两个连通片段之间沿 +x 留的空隙(Å)。
const FRAGMENT_GAP: f64 = 5.0;

/// 兄弟角被推歪超过这个角度就计一笔(弧度)。判据里那个 `STRAIN_WARN` 是同一个数。
const STRAIN_RADIANS: f64 = 0.087_266_462_599_716_47; // 5°
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
    /// **挂在环上、因此整片都没摆的原子。**
    ///
    /// 分量是在"去掉环原子"的子图上算的,所以苯的每个氢都会变成一个独立"片段"。
    /// 头一版把它们摆到 `原点 + shift` 并标成 `placed = true` —— 与真正键连的
    /// 环碳完全脱开,实测苯 6 个氢沿 x 轴一字排开、每个间隔 5 Å。
    /// 那是在撒谎:模块的规矩是"摆不了就说摆不了"。
    pub skipped_ring_attached: usize,
    /// **因为排布放不下这么多取代基而没摆的。**
    ///
    /// `Planar` 除父键之外只放得下 2 个、`Linear` 1 个、`Tetrahedral` 3 个。
    /// 要更多就少给,多出来的记在这里 —— 头一版是硬凑一个**重复**的扭转角,
    /// 于是两个原子被摆到同一个坐标上(`[Zn](C)(C)(C)C` 实测相距 0.000000 Å,
    /// 而 `complete()` 报 true)。
    pub skipped_arrangement: usize,
    /// **在被拒原子的下游、因此 BFS 到不了的原子。**
    ///
    /// 这个数是**走出来**的(从直接被拒的那些原子出发遍历),不是拿总数减出来的。
    pub skipped_downstream: usize,
    /// **既没摆好、也说不出原因的原子。判据把它钉在 0。**
    ///
    /// 上面每一个桶都是**真的在对应分支里累加**的,所以这一条是真残差,
    /// 不是恒等式。头一版把"剩下的"直接赋给 `skipped_hypervalent`,
    /// 于是守恒判据变成代数恒等式:往 BFS 里插一条静默 `continue` 丢掉 191 个原子,
    /// 守恒闸、覆盖率闸一个都不响。
    pub unaccounted: usize,
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
    /// **因为中心配位数 ≥ 5 而没摆的原子数。**
    ///
    /// 方案 §4.5 写明一期不保证这一档:`vsepr` 对它只有"均分了事"的 `Spread`,
    /// 不是真的三角双锥/八面体。实测钴的六配位中心用 `Spread` 摆出来的角是
    /// **56.25°**,而表里查不到、退到的是 109.47° —— **差 53.22°**。
    /// 与其硬摆一个错的,不如按本仓的规矩:**摆不了就说摆不了**。
    pub skipped_hypervalent: usize,
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
    /// **摆放树上每个原子的父亲**(根与没摆的是 `None`)。
    ///
    /// 判据需要它:一个中心上"父–子"那几个角**按构造精确等于表值**,
    /// 而兄弟之间的是推出来的、允许有解析偏差。分不清哪个是父亲,
    /// 就只能把那个偏差容差**减在所有原子对上** —— 实测那样会让某些中心的
    /// 容差涨到 180°,角判据在那儿等于整个关掉。
    pub parent: Vec<Option<u32>>,
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

/// **参考点共线时的兜底**:造一个垂直于 `p→c` 轴的点当祖父。
///
/// NeRF 要 `g`、`p`、`c` 不共线才定得出标架。而链上真有共线的:
/// 炔基(`C#C`)、腈、累积双键的中心,它们的父与祖父就在一条直线上。
/// 实测语料里 5 个原子因此摆不出来 —— 数量小,但那是个洞,不能留。
///
/// 垂直方向取哪个都行(它只定扭转角的零点,而绕轴的整体转动本来就是自由的),
/// 但**必须确定**:这里固定挑与轴最不平行的那个坐标轴去叉乘。
fn perpendicular_reference(p: Point3, c: Point3) -> Point3 {
    let u = (c - p).normalized().unwrap_or(Point3::new(1.0, 0.0, 0.0));
    // 与 u 最不平行的坐标轴
    let ax = if u.x.abs() <= u.y.abs() && u.x.abs() <= u.z.abs() {
        Point3::new(1.0, 0.0, 0.0)
    } else if u.y.abs() <= u.z.abs() {
        Point3::new(0.0, 1.0, 0.0)
    } else {
        Point3::new(0.0, 0.0, 1.0)
    };
    let perp = u
        .cross(ax)
        .normalized()
        .unwrap_or(Point3::new(0.0, 1.0, 0.0));
    p + perp
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

/// 摆一个**无环**分子。
///
/// `ranks` 要是规范秩(见模块文档)。有环的分子这一期不摆,
/// 会把环上的原子记进 [`Stats::skipped_ring`] 并留在 `placed = false`;
/// **挂在环上的那些片段也不摆**(见 [`Stats::skipped_ring_attached`])。
///
/// **永不 panic**:任何摆不了的情形都落进计数器,坐标留 `ORIGIN` 且 `placed = false`。
#[must_use]
#[allow(clippy::too_many_lines)]
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
        // **这条也要留痕**:模块文档说每一条提前返回都得留下痕迹,
        // 头一版这里一个计数器都不动,于是 `place(9 个原子, &[])` 报的是
        // "9 个原子、摆好 0 个",而每一个原因桶都是 0 —— 账对不上却没人说。
        stats.unaccounted = n;
        return Placed {
            coords,
            placed,
            parent,
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

    // **逐个连通片段摆。** 语料里带 `.` 的多片段分子(盐、共晶)不少 ——
    // 只从一个根 BFS 的话,别的片段一个原子都摆不到:实测 235 个原子因此丢掉,
    // 而它们会被记成"连不上",看着像 bug 其实是没实现。
    //
    // **连通分量要先算好,不能拿"还没摆的原子"当新片段的根。**
    // 头一版就是那么写的,结果:被拒绝的超配位中心、以及退化点,
    // 它们的取代基也是"还没摆的",于是被当成新片段**整体平移走** ——
    // 而那些原子跟已摆好的部分是**连着的**,键当场拉断。
    // 实测键长最大相对误差从 2.4e-15 跳到 **30 Å**(判据当场逮住)。
    let mut comp = vec![usize::MAX; n];
    let mut ncomp = 0;
    for a in 0..n {
        if on_ring[a] || comp[a] != usize::MAX {
            continue;
        }
        let mut stack = vec![a];
        comp[a] = ncomp;
        while let Some(x) = stack.pop() {
            #[allow(clippy::cast_possible_truncation)]
            for (y, _) in mol.neighbors(x as u32) {
                let y = y as usize;
                if !on_ring[y] && comp[y] == usize::MAX {
                    comp[y] = ncomp;
                    stack.push(y);
                }
            }
        }
        ncomp += 1;
    }

    // **挨着环的片段整个不摆。** 分量是在"去掉环原子"的子图上算的,
    // 于是苯的 6 个氢会各自变成一个"片段",被摆到 `原点 + shift` 上 ——
    // 与它真正键连的环碳完全脱开。头一版把它们标成 `placed = true`:
    // 实测苯 `placed = 6/12`,6 个氢沿 x 轴一字排开、每个间隔 5 Å,
    // 而它们各自的碳一个都没摆。甲苯则是整块甲基飘走。
    //
    // 那是在对调用方**撒谎** —— 模块的规矩是"摆不了就说摆不了,不硬摆"。
    // 二期把环摆出来之后,这些片段会以刚体装配的方式接回去(方案 §4.4)。
    let mut ring_attached = vec![false; ncomp];
    for a in 0..n {
        if on_ring[a] || comp[a] == usize::MAX {
            continue;
        }
        #[allow(clippy::cast_possible_truncation)]
        if mol.neighbors(a as u32).any(|(y, _)| on_ring[y as usize]) {
            ring_attached[comp[a]] = true;
        }
    }
    let mut skipped_attached = vec![false; n];
    for a in 0..n {
        if !on_ring[a] && comp[a] != usize::MAX && ring_attached[comp[a]] {
            skipped_attached[a] = true;
            stats.skipped_ring_attached += 1;
        }
    }

    // 直接被拒的原子(超配位/排布放不下/参考点退化)—— 它们**下游**的原子
    // BFS 也到不了,后面要单独走一遍算清楚,不能拿"剩下的都算它"糊过去。
    let mut refused_direct = vec![false; n];

    let mut shift = 0.0f64;
    for (cid, &attached) in ring_attached.iter().enumerate() {
        if attached {
            continue;
        }
        // 根:本分量里 (rank, 下标) 最小的那个
        let Some(root) = (0..n)
            .filter(|a| comp[*a] == cid)
            .min_by_key(|a| (ranks[*a], *a))
        else {
            continue;
        };
        #[allow(clippy::cast_possible_truncation)]
        let root = root as u32;
        let before: Vec<bool> = placed.clone();

        // ---- 起手:根放原点,第一个邻居沿 +x,其余绕 root–n1 轴铺开 ----
        placed[root as usize] = true;
        let nbrs: Vec<u32> = sorted_neighbors(mol, root, ranks)
            .into_iter()
            .filter(|x| !on_ring[*x as usize])
            .collect();
        let mut queue: std::collections::VecDeque<u32> = std::collections::VecDeque::new();
        // **超配位中心当根时也一个取代基都不摆。**
        //
        // 头一版把 `deg >= 5` 的判断放在摆 `n1` **之后**,于是同一件事有两个答案:
        // 那个中心当根时第一个邻居照摆(`FS(F)(F)(F)(F)F` 实测 `placed = 2`),
        // 经 BFS 到达时一个都不摆。走哪条路取决于规范秩,与化学无关。
        let root_deg = mol.neighbors(root).count();
        if root_deg >= 5 {
            stats.skipped_hypervalent += nbrs.len();
            for &x in &nbrs {
                refused_direct[x as usize] = true;
            }
        } else if let Some(&n1) = nbrs.first() {
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
            let deg = root_deg;
            let arr = arrangement(mol.atoms()[root as usize].hybridization, deg);
            let ts = child_torsions(arr, nbrs.len().saturating_sub(1));
            note_strain(&mut stats, mol, root, deg, arr, ts.len());
            for (k, &x) in nbrs.iter().skip(1).enumerate() {
                let Some(&t) = ts.get(k) else {
                    // **原因要分清楚。** 排布放不下(`Planar` 除父键外只放得下 2 个)
                    // 与"退化"是两回事 —— 头一版把它记成退化,名字是错的,
                    // 而且它头一版根本不存在:是硬凑一个**重复**的扭转角,
                    // 把两个原子摆到同一个坐标上。
                    stats.skipped_arrangement += 1;
                    refused_direct[x as usize] = true;
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
                let ang = angle_at_centre(&mut stats, mol, root, deg);
                match place_nerf(
                    v,
                    coords[n1 as usize],
                    coords[root as usize],
                    bl.value,
                    ang,
                    t,
                ) {
                    Some(p) => {
                        coords[x as usize] = p;
                        placed[x as usize] = true;
                        parent[x as usize] = Some(root);
                        queue.push_back(x);
                    }
                    None => {
                        stats.degenerate += 1;
                        refused_direct[x as usize] = true;
                    }
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
                for &x in &kids {
                    refused_direct[x as usize] = true;
                }
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
            // **配位数 ≥ 5 的中心不摆它的取代基。** `Spread` 只是"均分了事",
            // 实测钴的六配位中心给出 56.25° 而表值 109.47° —— 差 53.22°。
            // 硬摆一个错的不如如实说摆不了(方案 §4.5)。
            if deg >= 5 {
                stats.skipped_hypervalent += kids.len();
                for &x in &kids {
                    refused_direct[x as usize] = true;
                }
                continue;
            }
            let arr = arrangement(mol.atoms()[c as usize].hybridization, deg);
            let ts = child_torsions(arr, kids.len());
            note_strain(&mut stats, mol, c, deg, arr, ts.len());
            for (k, &x) in kids.iter().enumerate() {
                let Some(&t) = ts.get(k) else {
                    stats.skipped_arrangement += 1;
                    refused_direct[x as usize] = true;
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
                let ang = angle_at_centre(&mut stats, mol, c, deg);
                // 共线定不出标架时,换一个垂直的参考点再试一次(炔基/腈/累积双键)
                let got = place_nerf(g, coords[p as usize], coords[c as usize], bl.value, ang, t)
                    .or_else(|| {
                        place_nerf(
                            perpendicular_reference(coords[p as usize], coords[c as usize]),
                            coords[p as usize],
                            coords[c as usize],
                            bl.value,
                            ang,
                            t,
                        )
                    });
                match got {
                    Some(q) => {
                        coords[x as usize] = q;
                        placed[x as usize] = true;
                        parent[x as usize] = Some(c);
                        queue.push_back(x);
                    }
                    None => {
                        stats.degenerate += 1;
                        refused_direct[x as usize] = true;
                    }
                }
            }
        }

        // **把这一片挪到上一片右边。** 平移量要减掉本片自己的 min x ——
        // 头一版写的是 `shift = 上一片的 max x + 5`,没减新片的 min x,
        // 而片段从自己的根长出去时**也向 −x 伸**(实测单个片段的 min x 能到 −18.7 Å),
        // 于是两片会叠上:语料里最狠的一对相距 0.5 Å。
        let (mut minx, mut maxx) = (f64::INFINITY, f64::NEG_INFINITY);
        for i in 0..n {
            if placed[i] && !before[i] {
                minx = minx.min(coords[i].x);
                maxx = maxx.max(coords[i].x);
            }
        }
        if minx.is_finite() {
            let dx = shift - minx;
            for i in 0..n {
                if placed[i] && !before[i] {
                    coords[i] = coords[i] + Point3::new(dx, 0.0, 0.0);
                }
            }
            shift = maxx + dx + FRAGMENT_GAP;
        }
    }

    stats.placed = placed.iter().filter(|x| **x).count();

    // **下游原子要走出来,不能拿"剩下的都算它"顶。**
    //
    // 头一版这里是 `skipped_hypervalent = atoms − (placed + skipped_ring + degenerate)`,
    // 于是判据里那条守恒式变成**代数恒等式** —— 往 BFS 里插一条静默 `continue`
    // 丢掉 191 个原子,守恒闸、覆盖率闸、"连不上"闸**一个都没响**,
    // 唯一挡住的是那个总额上限常数,而它的报错还指向了一条毫不相干的路。
    //
    // 现在从直接被拒的那些原子出发走一遍,真正算出下游有多少;
    // 剩下的进 `unaccounted`,判据把它钉在 0 —— 那才是一道真闸。
    let mut blocked = refused_direct.clone();
    let mut stack: Vec<usize> = (0..n).filter(|i| refused_direct[*i]).collect();
    while let Some(x) = stack.pop() {
        #[allow(clippy::cast_possible_truncation)]
        for (y, _) in mol.neighbors(x as u32) {
            let y = y as usize;
            if !placed[y] && !on_ring[y] && !skipped_attached[y] && !blocked[y] {
                blocked[y] = true;
                stack.push(y);
            }
        }
    }
    stats.skipped_downstream = (0..n)
        .filter(|i| blocked[*i] && !refused_direct[*i])
        .count();
    stats.unaccounted = stats.atoms.saturating_sub(
        stats.placed
            + stats.skipped_ring
            + stats.skipped_ring_attached
            + stats.degenerate
            + stats.skipped_hypervalent
            + stats.skipped_arrangement
            + stats.skipped_downstream,
    );
    Placed {
        coords,
        placed,
        parent,
        stats,
    }
}

/// 查这个中心的键角,顺带记下它走的哪一级表。
fn angle_at_centre(stats: &mut Stats, mol: &MolBuilder, c: u32, deg: usize) -> f64 {
    let a = params::angle(
        mol.atoms()[c as usize].atomic_num,
        deg,
        mol.atoms()[c as usize]
            .flags
            .contains(omgkit_core::AtomFlags::AROMATIC),
        0,
        0,
    );
    stats.note_angle(a.source);
    a.value
}

/// 记一笔"这个中心的兄弟角被推歪了多少"。**只有真有兄弟(≥2 个子)时才算。**
///
/// 头一版只数 `Tetrahedral` 那一档,漏掉的正是偏得最狠的:平面中心的表角一旦
/// 不是 120°,兄弟角被强制偏 `|2π − 3θ|`,θ = 109.47° 时是 **31.6°** ——
/// 比四面体那一支能产生的任何值都大。实测漏报了大约一半(400 对真值约 817)。
///
/// 另外头一版在 `k == 0` 时无条件计数,哪怕那个中心只有一个子、根本没有兄弟。
fn note_strain(
    stats: &mut Stats,
    mol: &MolBuilder,
    c: u32,
    deg: usize,
    arr: crate::vsepr::Arrangement,
    n_children: usize,
) {
    if n_children < 2 {
        return;
    }
    let theta = params::angle(
        mol.atoms()[c as usize].atomic_num,
        deg,
        mol.atoms()[c as usize]
            .flags
            .contains(omgkit_core::AtomFlags::AROMATIC),
        0,
        0,
    )
    .value;
    if expected_sibling_skew(arr, theta) > STRAIN_RADIANS {
        stats.angle_strained += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::{angle_at, dihedral};
    use crate::vsepr::sibling_skew;

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

    /// **超配位中心一个取代基都不摆 —— 不管它是不是根。**
    ///
    /// 头一版把 `deg >= 5` 的判断放在摆第一个邻居**之后**,于是同一件事有两个答案:
    /// 那个中心当根时第一个邻居照摆,经 BFS 到达时一个都不摆。
    /// 走哪条路取决于规范秩 —— **与化学无关**。
    ///
    /// # 根那条路要**造**出来才测得到
    ///
    /// 根取的是分量里 `(rank, 下标)` 最小的原子,而 `FS(F)(F)(F)(F)F` 里
    /// 最小的是某个氟(配位 1)—— 硫从来当不上根,那条路测不到。
    /// `ranks` 本来就是参数,所以这里直接递一份把硫排在最前的秩。
    ///
    /// 判的是**契约**而不是"摆好几个":**没有任何原子的父亲是超配位中心**。
    #[test]
    fn a_hypervalent_centre_refuses_its_substituents_whether_or_not_it_is_the_root() {
        for smi in ["FS(F)(F)(F)(F)F", "CS(F)(F)(F)(F)F", "F[Co](F)(F)(F)(F)F"] {
            let (m, natural) = prep(smi);
            let n = m.num_atoms();
            // 找一个配位 ≥5 的中心
            let hv = (0..n)
                .map(|i| u32::try_from(i).expect("下标"))
                .find(|i| m.neighbors(*i).count() >= 5)
                .expect("该有一个超配位中心");
            // 两份秩:自然的那份,以及**把超配位中心排到最前**的那份(逼它当根)
            let mut forced = natural.clone();
            for r in &mut forced {
                *r += 1;
            }
            forced[hv as usize] = 0;
            for (which, ranks) in [("自然秩", &natural), ("逼它当根", &forced)] {
                let out = place(&m, ranks);
                for x in 0..n {
                    assert_ne!(
                        out.parent[x],
                        Some(hv),
                        "{smi}({which}):原子 {x} 的父亲是超配位中心 {hv} —— 它一个子都不该有"
                    );
                }
                assert!(
                    out.stats.skipped_hypervalent > 0,
                    "{smi}({which}) 该报超配位"
                );
                assert_eq!(
                    out.stats.unaccounted, 0,
                    "{smi}({which}) 的账对不上:{:?}",
                    out.stats
                );
            }
        }
    }

    /// **两个不相连的片段不许叠在一起。**
    ///
    /// 片段是沿 +x 一片一片挪开的,而平移量必须减掉**新片自己的 min x** ——
    /// 片段从根长出去时也向 −x 伸(实测单片的 min x 能到 −18.7 Å)。
    /// 头一版只写了"上一片的 max x + 5",于是两片会叠上。
    ///
    /// 这条测试是变异验证逼出来的:把那个减法去掉,当时**十二个变异里唯独它
    /// 一个闸都没响**。语料判据现在也有一条(`MIN_FRAGMENT_GAP`),
    /// 这里是逐分子的那道。
    #[test]
    fn two_disconnected_fragments_never_overlap() {
        let mut checked = 0;
        for smi in [
            // 取自 harness/corpus/large.smi
            "CCCCOP(O)(O)=O.CCCCOP(O)(=O)OCCCC",
            "NC(N)=O.OC(=O)C(O)=O",
            "CCCCCCCCCCCCCCCCCC(=O)O.CCCCCCCCCCCCCCCCCCO",
            "CC(=O)O.CC(=O)O.CC(=O)O",
        ] {
            let (m, r) = prep(smi);
            let out = place(&m, &r);
            assert!(out.complete(), "{smi} 该摆全:{:?}", out.stats);
            // 分量:测试自己算,不用 place 里那份
            let n = m.num_atoms();
            let mut comp = vec![usize::MAX; n];
            let mut nc = 0usize;
            for a in 0..n {
                if comp[a] != usize::MAX {
                    continue;
                }
                let mut st = vec![a];
                comp[a] = nc;
                while let Some(x) = st.pop() {
                    for (y, _) in m.neighbors(u32::try_from(x).expect("下标")) {
                        if comp[y as usize] == usize::MAX {
                            comp[y as usize] = nc;
                            st.push(y as usize);
                        }
                    }
                }
                nc += 1;
            }
            assert!(nc >= 2, "{smi} 该是多片段的");
            let mut worst = f64::INFINITY;
            for i in 0..n {
                for j in (i + 1)..n {
                    if comp[i] == comp[j] {
                        continue;
                    }
                    worst = worst.min(out.coords[i].dist(out.coords[j]));
                    checked += 1;
                }
            }
            assert!(
                worst >= 4.0,
                "{smi}:两个片段只隔 {worst:.3} Å —— 平移量减掉新片的 min x 了吗"
            );
        }
        assert!(checked > 500, "只验了 {checked} 对");
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
            // 摆好的那些坐标必须是有限数
            for (i, p) in out.coords.iter().enumerate() {
                if out.placed[i] {
                    assert!(p.is_finite(), "{smi} 第 {i} 个原子坐标不是有限数");
                }
            }
            // **两端都说摆好了的键,键长必须是表里那个值。**
            //
            // 头一版这条测试只查"报了有环 / 没摆全 / 坐标有限",**一根键长都不查** ——
            // 而键长正是当时唯一坏掉的东西:分量是在"去掉环原子"的子图上算的,
            // 于是苯的 6 个氢各成一个"片段"、被摆到 5 Å 间隔的一排上并标成
            // `placed = true`,与各自的碳完全脱开。
            for b in m.bonds() {
                let (x, y) = (b.begin as usize, b.end as usize);
                if !out.placed[x] || !out.placed[y] {
                    continue;
                }
                let wantb = params::bond_length(
                    m.atoms()[x].atomic_num,
                    m.atoms()[y].atomic_num,
                    b.order,
                    0,
                )
                .value;
                let got = out.coords[x].dist(out.coords[y]);
                assert!(
                    ((got - wantb) / wantb).abs() < 1e-12,
                    "{smi}:两端都摆好的键 {x}–{y} 长 {got:.6},该是 {wantb:.6}"
                );
            }
            // **不许摆一个连不回去的原子。** 摆好的原子,它的每个非环邻居也得摆好 ——
            // 不然那个坐标接在哪儿都不知道。二期把环摆出来之后这条自然还成立。
            for i in 0..m.num_atoms() {
                if !out.placed[i] {
                    continue;
                }
                let iu = u32::try_from(i).expect("下标");
                for (y, _) in m.neighbors(iu) {
                    assert!(
                        out.placed[y as usize],
                        "{smi}:原子 {i} 说摆好了,可它的邻居 {y} 没摆 —— 这个坐标接不回去"
                    );
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
        // 找不出主链就该红,**不能静默跳过** —— 那会让这条测试当场变恒真
        assert_eq!(chain.len(), 4, "没走出四个碳的主链");
        {
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
        // **断契约,不断"当前的偏差落在哪个区间"。**
        //
        // 头一版这里是 `0.10 < worst < 0.15` —— 下界断的是"**当前构造法的缺陷
        // 必须存在**"。哪天把扭转角改成让 6 个夹角最小二乘均摊(完全合法、
        // 而且更好),worst 会掉到 0.04,这条测试当场红,而代码是变好了。
        // 那种断言是能力长进的绊脚石,不是判据。
        //
        // 该断的是两条契约:
        //   父–子:精确等于表值(机器精度);
        //   兄弟 :精确等于解析式给的那个值。
        // 4 个取代基有 6 个夹角、模掉整体转动只有 5 个自由度 —— 超定,
        // 所以兄弟角只能是推出来的:`cos φ = cos²θ + sin²θ·cos120°`。
        let par = out.parent[mid as usize];
        let mut n_pc = 0;
        let mut n_sib = 0;
        for i in 0..nb.len() {
            for j in (i + 1)..nb.len() {
                let got = angle_at(
                    out.coords[nb[i] as usize],
                    out.coords[mid as usize],
                    out.coords[nb[j] as usize],
                )
                .expect("角")
                .to_degrees();
                if par == Some(nb[i]) || par == Some(nb[j]) {
                    assert!(
                        (got - want).abs() < 1e-9,
                        "父–子角 {got:.9}° 该精确等于表值 {want:.9}°"
                    );
                    n_pc += 1;
                } else {
                    let sib = want + sibling_skew(want.to_radians()).to_degrees();
                    assert!(
                        (got - sib).abs() < 1e-6,
                        "兄弟角 {got:.6}° 该等于解析式给的 {sib:.6}°"
                    );
                    n_sib += 1;
                }
            }
        }
        assert_eq!(n_pc + n_sib, 6);
        assert!(n_pc >= 1 && n_sib >= 1, "父–子 {n_pc} 对、兄弟 {n_sib} 对");
        assert!(worst < 0.15, "最大偏差 {worst:.4}°");
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
