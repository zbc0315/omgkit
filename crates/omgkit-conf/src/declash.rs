//! **消撞** —— 把叠在一起的原子推开,而且**只用扭转角**。
//!
//! # 为什么只转扭转角
//!
//! [`crate::build::place`] 摆出来的键长精确到 `3e-15`、键角精确落在解析边界上,
//! 这两项是这一期唯一"按定义正确"的东西。任何**平移单个原子**的消撞都会把它们弄脏 ——
//! 而绕一根键转动**整个子树**是刚体运动,键长键角**逐位不变**。
//! 判据 `declash_never_touches_a_bond_length_or_angle` 钉的就是这条。
//!
//! 代价是消撞能做的事有上界:固定键角撑不开的拥挤,转扭转角也撑不开
//! (实测真分子会把 C–C–C 开到 123.5°,而表里写的是 109.4°)。所以这一步
//! **不承诺**消干净,只承诺**只会变好、不会变坏**。
//!
//! # 三条让它跑得动的结构性质
//!
//! 1. **树上剪掉一条边,横跨切口的边有且仅有那一条。**
//!    于是要排除的 1-2/1-3 对可以直接写出来:`(b,c)`、`(x,c) x∈N(b)`、`(b,y) y∈N(c)` ——
//!    **不需要 `O(N²)` 的拓扑距离矩阵**。(证明:1-3 对 `i–k–j` 跨切口 ⟹ 两条边之一是 `(b,c)` ⟹ `k∈{b,c}`。)
//! 2. **绕轴转动的距离下界是闭式的。** 动点 `i` 与定点 `j` 在轴坐标系里是 `(ρ,z)`,
//!    `d² = Δz² + ρᵢ² + ρⱼ² − 2ρᵢρⱼcos(·)`,取遍所有转角的最小值是 `√(Δz² + (ρᵢ−ρⱼ)²)`。
//!    这个下界 ≥ 目标间距的对**永远撞不上**,一次性筛掉,不用每个候选角重算。
//! 3. **只看跨切口的对就够了。** 不跨切口的对这次转动一个都不动,所以拿"跨切口那部分"
//!    的 `(最小间距比, 罚和)` 做取舍,与拿全局的做取舍**给出同一个单调性保证**(见下)。
//!
//! # 单调性:只会变好,不会变坏
//!
//! 记跨切口的对为 `P`、其余为 `Q`(这次转动中 `Q` 一个都不变)。接受一次转动的条件是
//! `min_P r` **严格变大**,或者 `min_P r` 不变而 `Σ_P 罚` 变小。于是:
//!
//! - `min_P r ↑` ⟹ `min(min_P, min_Q)` 不减;
//! - `min_P r` 不变 ⟹ 全局最小不变,而全局罚和 `Σ_P + Σ_Q` 严格变小。
//!
//! 所以**全局的 `(最小间距比, −罚和) 按字典序单调不减**,而且两项都有界(比值封顶 1、罚和 ≥ 0),
//! 加上扫描轮数上限就必然停机。这里的比较**不带 epsilon 松弛** —— 松一点点就不是单调了,
//! 而"因为浮点噪声漏掉一次本可以接受的转动"只是少优化一点,不是错。

use crate::build::Placed;
use crate::geom::Point3;
use crate::params::vdw_radius;
use omgkit_core::{BondOrder, MolBuilder};

/// 目标间距 = 这个系数 × 两个 vdW 半径之和。
///
/// **与判据 `conf_audit` 的 `MIN_VDW_FRAC` 是同一个数** —— 消撞照着判据的尺子推,
/// 两边不许各写各的。
pub const TARGET_FRAC: f64 = 0.75;

/// 每根键试几个扭转偏移(均分 360°,第 0 个是"不动")。
///
/// # 这个数是量出来的,而且量到了一个"别再调"的结论
///
/// 全语料 1278 个无环分子,非键最小距离的**地板**(其余分辨率下的中位/代价见下):
///
/// | 候选 | 地板 Å | 中位 Å | 低于 1.6 Å 的分子 | μs/分子 | 收敛轮数 |
/// |---|---|---|---|---|---|
/// | 12 | **1.222** | 2.089 | 6 | 63.3 | 6 |
/// | 18 | 1.463 | 2.081 | 3 | 77.7 | 18 |
/// | **24** | **1.489** | 2.056 | 5 | **70.5** | **9** |
/// | 36 | 1.497 | 2.012 | 6 | 75.1 | 7 |
/// | 48 | 1.514 | 2.002 | 6 | 74.5 | 15 |
/// | 60 | 1.407 | 1.978 | 7 | 80.1 | 10 |
/// | 72 | 1.568 | 1.972 | 4 | 94.4 | 13 |
///
/// 读出来两件事:
///
/// 1. **12 太粗是真的** —— 它的地板 1.222 是全场最低,而 ≥18 一律落在 1.40~1.57;
/// 2. **≥18 之间的差别是噪声,不是趋势** —— 地板 1.463/1.489/1.497/1.514/1.407/1.568
///    不单调,"低于 1.6 Å 的分子数" 3/5/6/6/7/4 也不单调。挑其中最好的那个
///    等于在拟合贪心搜索的噪声。**唯一单调的两个量(中位、代价)都是越细越差。**
///
/// 所以取等价那一档里最便宜、收敛最快的 24(15°,仍含全部交错/共平面位:
/// 60/15 = 4、180/15 = 12)。**不要再往细调** —— 剩下的地板差距是键角撑不开,
/// 不是搜索不够细:候选加到 72 仍有 4 个分子低于 1.6 Å。
const N_CANDIDATES: usize = 24;

/// 最多扫几轮。**这是保险,不该顶到** —— 顶到了 [`DeclashStats::capped`] 会说。
///
/// 实测 24 候选下全语料最多用 9 轮就自己收敛了(把上限放到 200 也还是 9)。
/// 留到 16 是余量。注意轮数上限**不是**可调的旋钮:12 候选下放到 20 与放到 6
/// 给出逐位相同的结果 —— 算法自己会停。
const MAX_SWEEPS: usize = 16;

/// 一次消撞的分级计数。
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct DeclashStats {
    /// 可转的键有几根。
    pub rotatable: usize,
    /// 因为剪不断(在环上)而没算的键。
    pub in_ring: usize,
    /// 实际扫了几轮。
    pub sweeps: usize,
    /// **顶到轮数上限了** —— 还有能接受的转动没做完。
    ///
    /// 判据要看这个,不许自己去写"`sweeps >= 6`"那种硬编码:
    /// 上限一改,那种写法当场变成谎话。
    pub capped: bool,
    /// 接受了几次转动。
    pub moves: usize,
    /// 一共算了多少次原子对距离 —— 代价就看这个数。
    pub pair_evals: u64,
    /// **接受过的转动里,跨切口最小间距比最多掉了多少。按接受准则必须恒 ≤ 0。**
    ///
    /// 这个数存在的理由是变异验证逼出来的:把接受准则改成"只看罚和、不管最小值"
    /// (即去掉那道守卫),**全语料的单调性判据一个都没逮住** ——
    /// 因为逐键的最小值变差常常被分子别处更差的接触**遮住**,全局最小根本没动。
    /// 判据没写错,是它**天然测不到**逐键的那条准则。
    ///
    /// 所以把准则自己算的那两个数直接报出来,判据就能钉住它了。
    pub worst_move_regress: f64,
}

/// 全分子非键距离的普查。
///
/// `O(N²)`,**消撞内部不用它**(那边只看跨切口的对),这里是给判据和测试用的。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Contact {
    /// 最小非键距离(Å)。一个非键对都没有时是 `f64::INFINITY`。
    pub min_dist: f64,
    /// 最小的 `d / (TARGET_FRAC × (vdW_i + vdW_j))`。同样可能是 `INFINITY`。
    pub min_ratio: f64,
    /// `Σ max(0, 1 − d/d₀)²` —— 与 [`declash`] 内部用的罚函数是同一个式子。
    ///
    /// 有它,测试才能钉住[`declash`]真正承诺的那条**字典序**单调
    /// (`min_ratio` 不减;它一位不差时罚和不增),而不是"撞的对数不增" ——
    /// 后者是**假的**:把最狠的一处撞让开,常常换来两处更轻的
    /// (实测 `CC(C)(C)C(C)(C)C(C)(C)C` 从 8 对涨到 9 对,而最小间距比是升的)。
    /// 这正是方案里"目标要看最狠的一处,不能看撞的个数"那一条。
    pub penalty: f64,
    /// 低于目标间距的原子对数。**只是个报数,不是判据**(见 [`Contact::penalty`])。
    pub below: usize,
    /// 参与统计的原子对数(拓扑距离 ≥ 3 且两端都摆好了)。
    pub pairs: usize,
}

impl Default for Contact {
    fn default() -> Self {
        Self {
            min_dist: f64::INFINITY,
            min_ratio: f64::INFINITY,
            penalty: 0.0,
            below: 0,
            pairs: 0,
        }
    }
}

/// 普查非键距离:**拓扑距离 ≥ 3** 才算(1-2 与 1-3 是键长键角管的,转扭转角也动不了)。
#[must_use]
pub fn survey(mol: &MolBuilder, out: &Placed) -> Contact {
    let n = mol.num_atoms().min(out.coords.len()).min(out.placed.len());
    // 拓扑距离 ≤ 2 的邻接(含自己)
    let mut near = vec![vec![false; n]; n];
    for (i, row) in near.iter_mut().enumerate() {
        row[i] = true;
        let Ok(iu) = u32::try_from(i) else { continue };
        for (y, _) in mol.neighbors(iu) {
            if (y as usize) < n {
                row[y as usize] = true;
            }
            for (z, _) in mol.neighbors(y) {
                if (z as usize) < n {
                    row[z as usize] = true;
                }
            }
        }
    }
    let mut c = Contact::default();
    for (i, row) in near.iter().enumerate() {
        for (j, &close) in row.iter().enumerate().skip(i + 1) {
            if close || !out.placed[i] || !out.placed[j] {
                continue;
            }
            let d = out.coords[i].dist(out.coords[j]);
            let d0 = target(mol, i, j);
            let r = d / d0;
            c.pairs += 1;
            c.min_dist = c.min_dist.min(d);
            c.min_ratio = c.min_ratio.min(r);
            if r < 1.0 {
                c.penalty += (1.0 - r) * (1.0 - r);
                c.below += 1;
            }
        }
    }
    c
}

/// 两个原子该离多远。
fn target(mol: &MolBuilder, i: usize, j: usize) -> f64 {
    TARGET_FRAC * (vdw_radius(mol.atoms()[i].atomic_num) + vdw_radius(mol.atoms()[j].atomic_num))
}

/// 绕过 `origin`、方向 `u`(须已归一化)的轴转 `θ`(给的是 `cos`/`sin`)。罗德里格斯公式。
fn rotate_about(p: Point3, origin: Point3, u: Point3, cos: f64, sin: f64) -> Point3 {
    let v = p - origin;
    origin + v * cos + u.cross(v) * sin + u * (u.dot(v) * (1.0 - cos))
}

/// 跨切口那部分的 `(最小间距比, 罚和)`。**间距比封顶 1** —— 分得够开就不再往开里推,
/// 不然分子会被越拉越长而毫无意义。
fn score(pairs: &[(usize, usize, f64)], pos: &[Point3]) -> (f64, f64) {
    let mut min_r = 1.0f64;
    let mut sum = 0.0f64;
    for &(i, j, d0) in pairs {
        let r = pos[i].dist(pos[j]) / d0;
        min_r = min_r.min(r);
        if r < 1.0 {
            sum += (1.0 - r) * (1.0 - r);
        }
    }
    (min_r, sum)
}

/// **消撞**:反复绕可转的单键转子树,直到没有能让情况变好的转动。
///
/// `ranks` 要与 [`crate::build::place`] 用的是同一份规范秩 —— 键的处理次序按它排,
/// 不然同一个分子换个 SMILES 写法会走出不同的结果。
///
/// 就地改 `out.coords`。**键长与键角一位都不会变**(刚体转动)。
#[must_use]
#[allow(clippy::too_many_lines)]
pub fn declash(mol: &MolBuilder, ranks: &[u32], out: &mut Placed) -> DeclashStats {
    let mut stats = DeclashStats::default();
    let n = mol.num_atoms();
    if n < 4 || out.coords.len() < n || out.placed.len() < n || ranks.len() < n {
        return stats;
    }

    // 可转的键:单键、两端都摆好了。次序按规范秩,**不按存储序**。
    let mut cand: Vec<(u32, u32)> = Vec::new();
    for b in mol.bonds() {
        let (x, y) = (b.begin as usize, b.end as usize);
        if b.order != BondOrder::Single || x >= n || y >= n || !out.placed[x] || !out.placed[y] {
            continue;
        }
        // `b` 端取规范秩小的那个,`c` 端取大的 —— 转的是 `c` 这一侧的子树
        let (lo, hi) = if (ranks[x], x) <= (ranks[y], y) {
            (x, y)
        } else {
            (y, x)
        };
        cand.push((
            u32::try_from(lo).unwrap_or(0),
            u32::try_from(hi).unwrap_or(0),
        ));
    }
    cand.sort_by_key(|&(b, c)| (ranks[b as usize], b, ranks[c as usize], c));

    // 每根键**只预存拓扑**:动的那一侧、定的那一侧、要排除的 1-2/1-3 对。
    //
    // **转轴与原子对筛选一律现算,不许预存。** 头一版把 `origin`/`axis` 在开头
    // 按初始坐标算好存进来,结果转过别的键之后 `b`、`c` 自己也动了 ——
    // 拿旧轴转子树,相对整体仍是刚体,**相对 `b` 就不是**,键当场拉断:
    // 实测键长从 1.523 Å 变成 6.116 Å(判据 `declash_never_touches_a_bond_length_or_angle` 逮的)。
    // 同理那个"转一圈够不着"的筛选也依赖当前几何,预存下来就是错的筛。
    struct Job {
        b: usize,
        c: usize,
        moving: Vec<usize>,
        fixed: Vec<usize>,
        /// 至多 `deg(b)+deg(c)−1` 个,已按 `(小,大)` 归一化。
        excl: Vec<(usize, usize)>,
    }
    let mut jobs: Vec<Job> = Vec::new();

    for &(b, c) in &cand {
        let (bu, cu) = (b as usize, c as usize);
        // `c` 那一侧:从 `c` 出发、不许经过 `b`,只走摆好了的原子
        let mut side = vec![false; n];
        side[cu] = true;
        let mut stack = vec![cu];
        let mut reached_b = false;
        while let Some(x) = stack.pop() {
            let Ok(xu) = u32::try_from(x) else { continue };
            for (y, _) in mol.neighbors(xu) {
                let y = y as usize;
                if y >= n || !out.placed[y] || y == bu || side[y] {
                    continue;
                }
                side[y] = true;
                stack.push(y);
            }
            if mol.neighbors(xu).any(|(y, _)| y as usize == bu) && x != cu {
                reached_b = true;
            }
        }
        if reached_b {
            // 剪不断 —— 这根键在环上,转它会把环拉断
            stats.in_ring += 1;
            continue;
        }
        let s: Vec<usize> = (0..n).filter(|i| out.placed[*i] && side[*i]).collect();
        let t: Vec<usize> = (0..n).filter(|i| out.placed[*i] && !side[*i]).collect();
        // 两侧都得有 ≥2 个原子:只剩轴上那一个的话转动是恒等变换
        if s.len() < 2 || t.len() < 2 {
            continue;
        }
        // 1-2 / 1-3:树上剪一条边,跨切口的边只有 `(b,c)` 这一条,所以要排除的对能直接写出来 ——
        // 1-2 就是 `(b,c)`;1-3 是 `(x,c) x∈N(b)` 与 `(b,y) y∈N(c)`,**两端必有一个是 b 或 c**。
        let norm2 = |x: usize, y: usize| (x.min(y), x.max(y));
        let mut excl: Vec<(usize, usize)> = vec![norm2(bu, cu)];
        for (base, other) in [(bu, cu), (cu, bu)] {
            let Ok(x) = u32::try_from(base) else { continue };
            for (y, _) in mol.neighbors(x) {
                let y = y as usize;
                if y < n && y != other {
                    excl.push(norm2(y, other));
                }
            }
        }
        excl.sort_unstable();
        excl.dedup();

        // 动的取小的那一侧(相对几何一样,而候选偏移是绕整圈均分的,能达到的构型集合相同)
        let (moving, fixed) = if s.len() <= t.len() { (s, t) } else { (t, s) };
        stats.rotatable += 1;
        jobs.push(Job {
            b: bu,
            c: cu,
            moving,
            fixed,
            excl,
        });
    }

    if jobs.is_empty() {
        return stats;
    }

    #[allow(clippy::cast_precision_loss)]
    let step = std::f64::consts::TAU / N_CANDIDATES as f64;
    let mut scratch: Vec<Point3> = out.coords.clone();
    let mut pairs: Vec<(usize, usize, f64)> = Vec::new();
    for sweep in 0..MAX_SWEEPS {
        stats.sweeps = sweep + 1;
        let mut moved = false;
        for job in &jobs {
            // **轴现算** —— 上一根键转过之后 b、c 自己可能已经挪了位置
            let Some(axis) = (out.coords[job.c] - out.coords[job.b]).normalized() else {
                continue;
            };
            let origin = out.coords[job.c];
            // 轴坐标:每个原子的 (到轴的距离 ρ, 沿轴的坐标 z)
            let cyl = |p: Point3| -> (f64, f64) {
                let v = p - origin;
                let z = v.dot(axis);
                ((v - axis * z).norm(), z)
            };
            // 要看的原子对 **也现算**:那个"转一圈够不着"的下界依赖当前几何
            pairs.clear();
            for &i in &job.moving {
                let (ri, zi) = cyl(out.coords[i]);
                for &j in &job.fixed {
                    if (i == job.b || i == job.c || j == job.b || j == job.c)
                        && job.excl.binary_search(&(i.min(j), i.max(j))).is_ok()
                    {
                        continue;
                    }
                    let d0 = target(mol, i, j);
                    let (rj, zj) = cyl(out.coords[j]);
                    // 转遍一圈能达到的最小距离 = √(Δz² + (ρᵢ−ρⱼ)²) —— 够不着就永远撞不上
                    if (zi - zj).hypot(ri - rj) >= d0 {
                        continue;
                    }
                    pairs.push((i, j, d0));
                }
            }
            #[allow(clippy::cast_possible_truncation)]
            {
                stats.pair_evals += (job.moving.len() * job.fixed.len()) as u64;
            }
            if pairs.is_empty() {
                continue;
            }
            // 基准 = 现在的样子,**直接拿原坐标算**(转 0° 再算会差最后几位,
            // 那点差会把"不带 epsilon 的严格单调"毁掉)
            let (base_r, base_sum) = score(&pairs, &out.coords);
            stats.pair_evals += pairs.len() as u64;
            if base_r >= 1.0 {
                continue; // 这根键上已经没有撞的了
            }
            // 现任最优:一开始就是"不动"。候选按下标从小到大试,只有**严格**更好才换 ——
            // 于是平手一律归"不动",结果与遍历次序无关。
            let (mut best_k, mut best_r, mut best_sum) = (0usize, base_r, base_sum);
            for k in 1..N_CANDIDATES {
                #[allow(clippy::cast_precision_loss)]
                let (sin, cos) = (step * k as f64).sin_cos();
                scratch.copy_from_slice(&out.coords);
                for &i in &job.moving {
                    scratch[i] = rotate_about(out.coords[i], origin, axis, cos, sin);
                }
                let (r, sum) = score(&pairs, &scratch);
                stats.pair_evals += pairs.len() as u64;
                // **不带 epsilon**:要么最小间距比严格变大,要么它一位不差而罚和严格变小
                if r > best_r || (r >= best_r && sum < best_sum) {
                    best_k = k;
                    best_r = r;
                    best_sum = sum;
                }
            }
            if best_k > 0 {
                #[allow(clippy::cast_precision_loss)]
                let (sin, cos) = (step * best_k as f64).sin_cos();
                for &i in &job.moving {
                    out.coords[i] = rotate_about(out.coords[i], origin, axis, cos, sin);
                }
                stats.moves += 1;
                stats.worst_move_regress = stats.worst_move_regress.max(base_r - best_r);
                moved = true;
            }
        }
        if !moved {
            break;
        }
        stats.capped = sweep + 1 == MAX_SWEEPS;
    }
    stats
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::build::place;
    use crate::geom::{angle_at, dihedral};
    use crate::params;

    fn prep(smi: &str) -> (MolBuilder, Vec<u32>) {
        let mut m = omgkit_io::smiles::parse(smi).expect("SMILES 该解析得了");
        omgkit_chem::pipeline::sanitize(&mut m).expect("该 sanitize 得了");
        let r = omgkit_io::canon::classed_ranks(&m);
        omgkit_chem::add_explicit_hs(&mut m, &r);
        let r = omgkit_io::canon::classed_ranks(&m);
        (m, r)
    }

    /// **这是消撞存在的全部前提:键长与键角一位都不许变。**
    ///
    /// 刚体转动理应如此,但"理应"不算数 —— 一旦哪天有人把消撞改成推单个原子,
    /// 这条会当场红。语料判据看的是全语料的上限,这条看的是逐根键。
    #[test]
    fn declash_never_touches_a_bond_length_or_angle() {
        let mut checked = 0;
        let mut total_moves = 0;
        for smi in [
            "CC(C)CC(CC(C)C)OB(OC(CC(C)C)CC(C)C)OC(CC(C)C)CC(C)C", // 判据报的最狠那个
            "CC(C)(C)CC(C)(C)C",
            "CCCCCCCC",
            "OCC(O)CO",
            "CC(=O)OC",
            "FC(F)(F)C(Cl)Br",
        ] {
            let (m, r) = prep(smi);
            let mut out = place(&m, &r);
            assert!(out.complete(), "{smi} 该摆全");
            let before = out.coords.clone();
            let st = declash(&m, &r, &mut out);
            total_moves += st.moves;

            for b in m.bonds() {
                let (x, y) = (b.begin as usize, b.end as usize);
                let was = before[x].dist(before[y]);
                let now = out.coords[x].dist(out.coords[y]);
                assert!(
                    ((now - was) / was).abs() < 1e-12,
                    "{smi}:键长从 {was:.9} 变成 {now:.9}"
                );
                checked += 1;
            }
            for k in 0..m.num_atoms() {
                let ku = u32::try_from(k).expect("下标");
                let nb: Vec<usize> = m.neighbors(ku).map(|(y, _)| y as usize).collect();
                for i in 0..nb.len() {
                    for j in (i + 1)..nb.len() {
                        let (Some(was), Some(now)) = (
                            angle_at(before[nb[i]], before[k], before[nb[j]]),
                            angle_at(out.coords[nb[i]], out.coords[k], out.coords[nb[j]]),
                        ) else {
                            continue;
                        };
                        assert!(
                            (now - was).abs() < 1e-9,
                            "{smi}:{k} 上的键角从 {:.6}° 变成 {:.6}°",
                            was.to_degrees(),
                            now.to_degrees()
                        );
                        checked += 1;
                    }
                }
            }
        }
        assert!(checked > 300, "只验了 {checked} 项 —— 这条快变恒真了");
        // **按整组算,不按单个分子算。** 头一版逐个分子要求"至少动一次",
        // 结果全反式的正辛烷当场红 —— 它本来就一处不撞,消撞**不动它才是对的**。
        // 这里要防的是"消撞根本没跑,于是键长键角当然没变"那种恒真。
        assert!(
            total_moves > 10,
            "整组只转了 {total_moves} 次 —— 消撞是不是没跑"
        );
    }

    /// **只会变好,不会变坏。** 承诺的是那条**字典序**:`min_ratio` 不减,
    /// 它一位不差时罚和不增。
    ///
    /// **不是**"撞的对数不增" —— 头一版这条断言就是那么写的,当场红:
    /// `CC(C)(C)C(C)(C)C(C)(C)C` 从 8 对涨到 9 对,而最小间距比是升的。
    /// 把最狠的一处让开、换来两处更轻的,**正是方案要的**(目标看最狠的一处,
    /// 不看个数)。错的是断言,不是代码。
    #[test]
    fn declash_never_makes_the_worst_contact_worse() {
        let mut improved = 0;
        for smi in [
            "CC(C)CC(CC(C)C)OB(OC(CC(C)C)CC(C)C)OC(CC(C)C)CC(C)C",
            "CC(C)(C)CC(C)(C)C",
            "CCCCCCCC",
            "CC(C)(C)C(C)(C)C(C)(C)C",
            "OCC(O)CO",
            "CCOCCOCCOCC",
            "N",
            "C",
            "CC",
        ] {
            let (m, r) = prep(smi);
            let mut out = place(&m, &r);
            let was = survey(&m, &out);
            let st = declash(&m, &r, &mut out);
            let now = survey(&m, &out);
            assert_eq!(now.pairs, was.pairs, "{smi}:统计的原子对数不该变");
            // 逐次转动的那条守卫 —— 全局判据遮得住它,这里直接看准则自己算的数
            assert!(
                st.worst_move_regress <= 0.0,
                "{smi}:有一次转动让跨切口的最小间距比掉了 {:.9}",
                st.worst_move_regress
            );
            // **不带 epsilon** —— 承诺是逐位的:Q 那部分的坐标一位没动,
            // P 那部分与取舍时算的是同一串运算
            assert!(
                now.min_ratio >= was.min_ratio,
                "{smi}:最小间距比从 {:.17} 掉到 {:.17}",
                was.min_ratio,
                now.min_ratio
            );
            if now.min_ratio == was.min_ratio {
                assert!(
                    now.penalty <= was.penalty + 1e-12,
                    "{smi}:最小间距比没动,罚和却从 {:.9} 涨到 {:.9}",
                    was.penalty,
                    now.penalty
                );
            }
            if now.min_ratio > was.min_ratio + 1e-9 {
                improved += 1;
            }
        }
        assert!(
            improved >= 3,
            "只有 {improved} 个分子变好了 —— 消撞是不是没起作用"
        );
    }

    /// **不许转双键、三键、芳香键** —— 转了就是另一个立体异构体。
    ///
    /// # 分子必须真的挤,不然这条是空转的
    ///
    /// 头一版用的是丁烯、顺丁烯二酸这种,变异验证当场揭穿:**把单键那道过滤去掉,
    /// 这条测试照样绿** —— 那几个分子的双键上根本没有撞,所以"允许转"和"不许转"
    /// 给出同一个结果。下面这几个是从语料里挑的**四取代、两端都带大基团**的烯烃,
    /// 双键两侧顺位的取代基一定撞,一旦放行就会被转走。
    #[test]
    fn a_double_bond_is_never_rotated() {
        let mut checked = 0;
        for smi in [
            // 全都取自 harness/corpus/large.smi
            "C(#N)S/C(=C(/[N+](=O)[O-])\\Cl)/C(=C(Cl)Cl)Cl",
            "FC(F)=C(F)C(F)(F)C(F)(Cl)C(F)(F)C(F)(Cl)C(F)(F)Cl",
            "C(=C(/C(F)(F)F)\\N)(\\C(=[N-])C(F)(F)F)/F",
            "CN(C)C=C(C=[N+](C)C)C(F)(F)F",
            "CCN(CC)C(=O)/C(=C/COP(=O)(OC)OC)Cl",
        ] {
            let (m, r) = prep(smi);
            let mut out = place(&m, &r);
            assert!(out.complete(), "{smi} 该摆全");
            // 先确认这个分子**真的挤** —— 不挤的话下面验的就是个恒真命题
            let c0 = survey(&m, &out);
            assert!(
                c0.below > 0,
                "{smi} 一处都不撞 —— 拿它验'不许转双键'是空转的"
            );
            let before = out.coords.clone();
            let st = declash(&m, &r, &mut out);
            assert!(st.moves > 0, "{smi} 消撞一次都没转 —— 同样是空转");
            for b in m.bonds() {
                if b.order == BondOrder::Single {
                    continue;
                }
                let (x, y) = (b.begin as usize, b.end as usize);
                let (xu, yu) = (b.begin, b.end);
                let Some(g) = m.neighbors(xu).map(|(z, _)| z).find(|z| *z != yu) else {
                    continue;
                };
                let Some(h) = m.neighbors(yu).map(|(z, _)| z).find(|z| *z != xu) else {
                    continue;
                };
                let (Some(was), Some(now)) = (
                    dihedral(before[g as usize], before[x], before[y], before[h as usize]),
                    dihedral(
                        out.coords[g as usize],
                        out.coords[x],
                        out.coords[y],
                        out.coords[h as usize],
                    ),
                ) else {
                    continue;
                };
                // **按绕回比** —— 180° 与 −180° 是同一个角,直接相减会诬告
                let raw = (now - was).abs();
                let diff = raw.min(std::f64::consts::TAU - raw);
                assert!(
                    diff < 1e-9,
                    "{smi}:双/三键 {x}={y} 上的二面角从 {:.4}° 转到了 {:.4}°",
                    was.to_degrees(),
                    now.to_degrees()
                );
                checked += 1;
            }
        }
        // 只有两端各自还有别的邻居的双/三键才验得了二面角(末端的 =O、#N 验不了)
        assert!(checked >= 7, "只验了 {checked} 根双/三键 —— 这条快变恒真了");
    }

    /// **环上的键剪不断,不许转** —— 转了会把环拉开。
    ///
    /// # 这条测试为什么要手工摆坐标
    ///
    /// 一期**根本不摆环上的原子**,所以那道守卫在真实流程里一次都走不到 ——
    /// 变异验证把它整个删掉,全语料一个数都不变。它是给二期铺的路,
    /// 而没被测试盯住的防御代码等于没有。
    ///
    /// 于是这里手工造一个"环上的原子也都摆好了"的局面:环己烷的六个碳放成平面六边形,
    /// 氢挂在各自碳的外侧。几何不像真分子无所谓 —— 要验的是**拓扑判断**:
    /// 六根环键必须全被认成剪不断,而且一根键都不许被拉长。
    #[test]
    fn a_ring_bond_is_refused_because_the_cut_does_not_disconnect() {
        let (m, r) = prep("C1CCCCC1");
        let n = m.num_atoms();
        // 环上的六个碳,按连接顺序走一圈
        let carbons: Vec<u32> = (0..n)
            .filter(|i| m.atoms()[*i].atomic_num == 6)
            .map(|i| u32::try_from(i).expect("下标"))
            .collect();
        assert_eq!(carbons.len(), 6);
        let mut ring = vec![carbons[0]];
        while ring.len() < 6 {
            let last = *ring.last().expect("非空");
            let nxt = m
                .neighbors(last)
                .map(|(y, _)| y)
                .find(|y| carbons.contains(y) && !ring.contains(y));
            match nxt {
                Some(y) => ring.push(y),
                None => break,
            }
        }
        assert_eq!(ring.len(), 6, "没走出一圈六元环");

        // 正六边形:边长 = 外接圆半径
        let bond = 1.523_f64;
        let mut coords = vec![Point3::ORIGIN; n];
        let mut placed = vec![false; n];
        for (k, &c) in ring.iter().enumerate() {
            #[allow(clippy::cast_precision_loss)]
            let a = std::f64::consts::TAU * (k as f64) / 6.0;
            coords[c as usize] = Point3::new(bond * a.cos(), bond * a.sin(), 0.0);
            placed[c as usize] = true;
        }
        // 氢:沿径向朝外、上下各一个
        for &c in &ring {
            let p = coords[c as usize];
            let out_dir = p.normalized().expect("环心不在碳上");
            for (k, (h, _)) in m
                .neighbors(c)
                .filter(|(y, _)| !ring.contains(y))
                .enumerate()
            {
                let z = if k % 2 == 0 { 1.0 } else { -1.0 };
                let d = Point3::new(out_dir.x * 0.5, out_dir.y * 0.5, z * 0.866);
                let d = d.normalized().expect("方向非零");
                coords[h as usize] = p + d * 1.09;
                placed[h as usize] = true;
            }
        }
        assert!(placed.iter().all(|x| *x), "手工摆的坐标没盖住全部原子");

        let mut out = Placed {
            coords,
            placed,
            stats: crate::build::Stats::default(),
        };
        let before = out.coords.clone();
        let st = declash(&m, &r, &mut out);
        assert_eq!(
            st.in_ring, 6,
            "六根环键该全被认成剪不断,实得 {}",
            st.in_ring
        );
        for b in m.bonds() {
            let (x, y) = (b.begin as usize, b.end as usize);
            let was = before[x].dist(before[y]);
            let now = out.coords[x].dist(out.coords[y]);
            assert!(
                ((now - was) / was).abs() < 1e-12,
                "键长从 {was:.9} 变成 {now:.9} —— 环被拉开了"
            );
        }
    }

    /// **同一分子换个写法,消撞的结果判决不能变。**
    ///
    /// 键的处理次序按规范秩排的理由就在这里 —— 按存储序排的话,
    /// 同一个分子换个原子顺序就会走出不同的转动序列。
    ///
    /// # 写法是**生成**的,不是手挑的
    ///
    /// 头一版是手写几组"同分子不同 SMILES",变异验证当场揭穿:**把排序键换成存储序,
    /// 这条测试照样绿** —— 手挑的那几组恰好给出同一个结果(还有一组我写成了
    /// 两个一模一样的串,纯空转)。现在改成用 [`omgkit_io::smiles::write_with_priority`]
    /// 按几种不同的优先级把同一个分子重写一遍,存储序是**真的**不同的,
    /// 而且下面会先断言它确实不同,不然这条测试还是空转。
    #[test]
    fn declash_does_not_depend_on_how_the_molecule_is_written() {
        let mut forms_checked = 0;
        for smi in [
            // 取自 harness/corpus/large.smi,都是挤到消撞非动不可的
            "CC(C)CC(CC(C)C)OB(OC(CC(C)C)CC(C)C)OC(CC(C)C)CC(C)C",
            "C(#N)S/C(=C(/[N+](=O)[O-])\\Cl)/C(=C(Cl)Cl)Cl",
            "CCCCCCCCCCCC(=O)OC(C)C(=O)NC(C)(C)CC(C)(C)C",
            "CCN(CC)C(=O)/C(=C/COP(=O)(OC)OC)Cl",
        ] {
            let mut base = omgkit_io::smiles::parse(smi).expect("SMILES 该解析得了");
            omgkit_chem::pipeline::sanitize(&mut base).expect("该 sanitize 得了");
            let n = base.num_atoms();
            // 几种不同的优先级 → 几种不同的写法
            let mut forms = vec![smi.to_string()];
            for style in 0..3u32 {
                #[allow(clippy::cast_possible_truncation)]
                let pri: Vec<u32> = (0..n)
                    .map(|i| match style {
                        0 => (n - 1 - i) as u32,
                        1 => ((i * 7 + 3) % n) as u32,
                        _ => ((i * 13 + 5) % n) as u32,
                    })
                    .collect();
                forms.push(omgkit_io::smiles::write_with_priority(&base, &pri).smiles);
            }

            let mut seen: Option<(u64, usize, usize, u64)> = None;
            let mut orders: Vec<Vec<u8>> = Vec::new();
            for form in &forms {
                let (m, r) = prep(form);
                orders.push(m.atoms().iter().map(|a| a.atomic_num).collect());
                let mut out = place(&m, &r);
                let st = declash(&m, &r, &mut out);
                let c = survey(&m, &out);
                let key = (
                    c.min_ratio.to_bits(),
                    c.below,
                    st.moves,
                    c.penalty.to_bits(),
                );
                match seen {
                    None => seen = Some(key),
                    Some(k) => assert_eq!(key, k, "{smi} 的写法 {form} 给出了不同的消撞判决"),
                }
                forms_checked += 1;
            }
            // **这些写法的存储序必须真的不同** —— 都一样的话上面验的是恒真命题
            assert!(
                orders.iter().any(|o| *o != orders[0]),
                "{smi} 的几种写法给出了同一个存储序 —— 这条测试是空转的"
            );
        }
        assert!(forms_checked >= 16, "只验了 {forms_checked} 个写法");
    }

    /// **停机。** 扫到没有可接受的转动就该自己停,不该顶到轮数上限。
    #[test]
    fn declash_converges_before_the_sweep_cap() {
        for smi in [
            "CCCCCC",
            "OCC(O)CO",
            "CC(C)(C)CC(C)(C)C",
            "CC(C)CC(CC(C)C)OB(OC(CC(C)C)CC(C)C)OC(CC(C)C)CC(C)C",
        ] {
            let (m, r) = prep(smi);
            let mut out = place(&m, &r);
            let st = declash(&m, &r, &mut out);
            // 看 `capped`,**不看 `sweeps < MAX_SWEEPS`** —— 后者把上限的值抄了一份,
            // 上限一改这条就变成谎话
            assert!(!st.capped, "{smi} 扫满了 {} 轮还没收敛", st.sweeps);
        }
    }

    /// 小分子、单原子、空的:不许 panic,也不该乱动。
    #[test]
    fn tiny_and_degenerate_inputs_are_left_alone() {
        for smi in ["C", "N", "O", "CC", "[H][H]", "C#N"] {
            let (m, r) = prep(smi);
            let mut out = place(&m, &r);
            let before = out.coords.clone();
            let st = declash(&m, &r, &mut out);
            assert_eq!(st.moves, 0, "{smi} 上没有可转的键,却动了");
            for (i, p) in out.coords.iter().enumerate() {
                assert_eq!(*p, before[i], "{smi} 第 {i} 个原子被挪了");
            }
        }
        // ranks 给短了也不许 panic
        let (m, _) = prep("CCCC");
        let mut out = place(&m, &omgkit_io::canon::classed_ranks(&m));
        let st = declash(&m, &[], &mut out);
        assert_eq!(st.moves, 0);
    }

    /// `survey` 的口径必须与判据一致:**拓扑距离 ≥ 3**。
    ///
    /// 丙烷的 1-3(两端的碳)不算,而两端碳上的氢与另一端的氢(1-4 及更远)要算。
    #[test]
    fn survey_counts_only_pairs_three_bonds_apart() {
        let (m, r) = prep("CCC");
        let out = place(&m, &r);
        let c = survey(&m, &out);
        // 丙烷 C3H8 共 11 个原子 = 55 对;减去键(10 对)与 1-3 对
        let n = m.num_atoms();
        let mut close = 0;
        for i in 0..n {
            for j in (i + 1)..n {
                let (iu, ju) = (
                    u32::try_from(i).expect("下标"),
                    u32::try_from(j).expect("下标"),
                );
                let bonded = m.neighbors(iu).any(|(y, _)| y == ju);
                let geminal = m
                    .neighbors(iu)
                    .any(|(y, _)| m.neighbors(y).any(|(z, _)| z == ju));
                if bonded || geminal {
                    close += 1;
                }
            }
        }
        assert_eq!(c.pairs, n * (n - 1) / 2 - close, "统计的对数与手数的对不上");
        assert!(c.pairs > 0);
        // 1-3 的碳(2.5 Å 上下)绝不能混进来 —— 混进来的话最小距离会掉到 2.6 以下
        assert!(
            c.min_dist > 1.5,
            "丙烷的最小非键距离 {:.3} —— 1-3 对漏进来了?",
            c.min_dist
        );
        // 目标间距用的是同一个系数。**容差 1e-6 不是 1e-9** —— 元素表里 `rvdw` 是 `f32`,
        // 1.2f32 转成 f64 是 1.2000000476837158,乘出来 1.8000000715255737。
        let d0 = TARGET_FRAC * (params::vdw_radius(1) + params::vdw_radius(1));
        assert!(
            (d0 - 1.8).abs() < 1e-6,
            "H–H 的目标间距该是 1.8 Å,算出来 {d0}"
        );
    }
}
