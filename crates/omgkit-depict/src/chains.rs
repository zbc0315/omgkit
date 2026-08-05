//! 链与取代基:给一个原子周围**还没放置**的邻居分配方向。
//!
//! # 规则
//!
//! 理想夹角按度数定:度数 ≤ 3 用 [`Style::chain_angle_deg`](crate::style::Style)
//! (两套内置规范都是 120°),度数 ≥ 4 用 360/度数。
//!
//! 度数 2 的原子**必须**取 120° 而不是 180° —— 否则链会拉成直线,那既不是化学
//! 惯例,也让后面的可旋转键无处可翻。
//!
//! 新邻居放进**最大的空闲扇区**并在其中均分。这条规则在各种度数下都退化得对:
//!
//! | 已占方向 | 空闲扇区 | 做法 |
//! |---|---|---|
//! | 0 个(起点) | 360° | 从一个固定角度起,按理想夹角铺开 |
//! | 1 个(链上一步) | 360° | 取 ±理想夹角,符号由锯齿决定 |
//! | ≥2 个 | < 360° | 在最大空隙内均分 |
//!
//! 中间那一档不能套用"扇区居中":360° 扇区的中心正好是已占方向的反向,那就把
//! 链拉成了直线。
//!
//! # 锯齿的符号必须可复现
//!
//! 直链每一步左右交替才成锯齿。交替本身要有个起点,而那个起点若取自原子的**存储
//! 下标**,同一个分子换种 SMILES 写法就会得到镜像的图。这里由父原子的符号翻转
//! 得到,起点则来自规范秩 —— 写法就影响不到。

use std::collections::BTreeMap;

use omgkit_core::{BondOrder, MolBuilder};

use crate::geom::{segments_cross, Point2, BOND_LEN};
use crate::style::Style;

/// 起点原子第一根键的方向。
///
/// 取 30° 是惯例:六元环按这个角度画出来是"尖朝上"的标准姿态,直链则呈水平
/// 锯齿。取 0° 会让直链变成一条水平线加上下交替的折点,观感上偏斜。
const SEED_ANGLE: f64 = std::f64::consts::FRAC_PI_6;

/// 一次放置的结果。
pub(crate) struct Placed {
    /// 原子
    pub atom: u32,
    /// 坐标
    pub at: Point2,
    /// 它自己的锯齿符号,传给它的子代取反
    pub zig: i8,
}

/// 给 `a` 周围还没放置的邻居 `todo` 分配坐标。
///
/// `todo` 必须**已按规范秩排好**;顺序决定谁分到哪个方向,拿存储下标排就会
/// 引入写法依赖。
pub(crate) fn place_neighbours(
    mol: &MolBuilder,
    a: u32,
    pos: &BTreeMap<u32, Point2>,
    todo: &[u32],
    ranks: &[u32],
    style: &Style,
    zig: i8,
) -> Vec<Placed> {
    if todo.is_empty() {
        return Vec::new();
    }
    let center = pos[&a];

    // 已占方向:所有已经放好的邻居。
    //
    // **排序要量化后再按规范秩打破平局。** 两个方向在数学上相等、浮点上差
    // 1e-16 时,直接比大小会让它们的先后取决于算到那一步的运算次序 —— 而
    // `allocate` 挑最大空隙时用的正是这个次序,于是同一个分子换个写法,某个
    // 取代基就换了个方向挂。实测:一个稠三环的甲基因此差了 120°。
    //
    // **而且要先化到 `[0, 2π)`。** `angle()` 走的是 `atan2`,值域是 `(-π, π]`
    // —— **−180° 与 +180° 是同一个方向,却排在序列的两头**。末位差 4.4e-16 就
    // 足以决定它落在哪一端,于是同一组方向在两种写法下排出两个不同的序列,
    // `largest_gap` 看到的空隙序列跟着不同,取代基差 120°。
    //
    // 实测(`C[C]1(CCC[C]2(C)[CH]1CCC3=C2C=C(O)C=C3)C(O)=O`):
    //
    // ```text
    // 写法 A: occ = [-3.14159265358979312, -1.04719755119659808, 1.04719755119659763]
    // 写法 B: occ = [-1.04719755119659808,  1.04719755119659763, 3.14159265358979267]
    // ```
    //
    // 化到 `[0, 2π)` 之后两边都是 `[1.047, 3.142, 5.236]`,一模一样。
    const QUANT: f64 = 1e9;
    let mut occ: Vec<((i64, u32), f64)> = mol
        .neighbors(a)
        .filter_map(|(n, _)| {
            pos.get(&n).map(|p| {
                let t = (*p - center).angle().rem_euclid(std::f64::consts::TAU);
                // **化到 `[0, 2π)` 还不够:断点只是从 ±π 挪到了 0/2π。**
                // 一个 −1.8e-16 的角化出来是 `6.28318530717958534`(将近 2π),
                // 而 0 与 2π 同样是一个方向 —— 两种写法照样排出两个序列。
                // 实测踩到过两遍,第二遍就是这个:
                //
                // ```text
                // 写法 A: occ = [2.0944, 4.1888, 6.28318530717958534]
                // 写法 B: occ = [0,      2.0944, 4.1888]
                // ```
                //
                // 贴着 2π 的一律掐回 0。容差取 1e-9:真正不同的两个方向至少
                // 差 30°,而浮点噪声在 1e-16 量级。
                let t = if std::f64::consts::TAU - t < 1e-9 {
                    0.0
                } else {
                    t
                };
                #[allow(clippy::cast_possible_truncation)]
                (((t * QUANT).round() as i64, ranks[n as usize]), t)
            })
        })
        .collect();
    occ.sort_unstable_by_key(|x| x.0);
    let occupied: Vec<f64> = occ.iter().map(|x| x.1).collect();

    let ideal = ideal_angle(mol, a, style);
    let mut dirs = allocate(&occupied, todo.len(), ideal, zig);

    // **只有一个已占方向时,±理想角两侧都是合法的,而 `allocate` 按锯齿的符号
    // 盲选。** 挑错一侧,挂在环上的臂就朝环卷回去 —— 臂上的取代基撞到环,只能
    // 按 30° 一档挪,挪出来就是 120∓30 = 90° 或 150°。
    //
    // 实测:阿司匹林乙酰基那个 sp² 碳,三个角是 90/120/150,而它三个都该是
    // 120°;整条乙酰基臂折回来贴着苯环。
    //
    // 所以两侧都算一遍**拥挤度**,挑空的那边。直链两侧一样空,分不出高下时
    // 保持锯齿的选择 —— 锯齿因此不受影响。
    // 已经占住的位置,以及已经画出来的键。新原子不许落在前者上、新键不许与
    // 后者交叉 —— 见 [`free_direction`]。
    let mut taken: Vec<Point2> = pos.values().copied().collect();
    let mut drawn: Vec<(Point2, Point2)> = mol
        .bonds()
        .iter()
        .filter_map(|b| Some((*pos.get(&b.begin)?, *pos.get(&b.end)?)))
        .collect();

    if occupied.len() == 1 && !todo.is_empty() {
        // 两侧都算一遍**拥挤度**,挑空的那边;分不出高下时保持锯齿的选择 ——
        // 直链两侧一样空,锯齿因此不受影响。
        //
        // 试过一个更保守的版本:"只在这一侧确实会被迫歪角时才换边"(数一数
        // `free_direction` 会挪走几个)。它**救不了阿司匹林的 ACS 那张** ——
        // 乙酰基那个 sp² 碳仍是 90/120/150,因为羰基氧的理想位置在放它的那一刻
        // 还没被占,是后面的原子挤过来的。拥挤度看的是整体,才拦得住。
        //
        // 拥挤度:新位置到每个已放好的原子的平方反比之和(与 RDKit 的 density
        // 同一个口径)。**量化之后再比** —— 直接比浮点会让"分不出高下"取决于
        // 末位,而那一位取决于运算次序,写法一换就可能翻边。
        #[allow(clippy::cast_possible_truncation)]
        let crowd = |ds: &[f64]| -> i64 {
            let mut sum = 0.0_f64;
            for t in ds {
                let p = center + Point2::new(BOND_LEN, 0.0).rotated(*t);
                for q in pos.values() {
                    sum += 1.0 / (p.dist(*q).powi(2) + 1e-6);
                }
            }
            (sum * 1e6).round() as i64
        };
        let mirror: Vec<f64> = dirs.iter().map(|t| 2.0 * occupied[0] - t).collect();
        if crowd(&mirror) < crowd(&dirs) {
            dirs = mirror;
        }
    }

    debug_assert_eq!(dirs.len(), todo.len(), "方向数必须与待放邻居数相等");

    let mut out = Vec::with_capacity(todo.len());
    for (&atom, theta) in todo.iter().zip(dirs) {
        let theta = free_direction(center, theta, &occupied, &taken, &drawn);
        let at = center + Point2::new(BOND_LEN, 0.0).rotated(theta);
        taken.push(at);
        drawn.push((center, at));
        out.push(Placed {
            atom,
            at,
            // 子代取反,直链就走出锯齿
            zig: -zig,
        });
    }
    out
}

/// 从 `ideal` 出发,找一个不会与已放好的原子重合的方向。
///
/// # 为什么宁可歪着也不重合
///
/// 两个原子叠在同一点上时,它们各自的键首尾相接 —— **图上就多出一个分子里
/// 没有的环**,而读者没有任何办法看出那个环是假的。角度偏离理想值只是难看,
/// 不会让人读错结构。
///
/// 实测:一个三萜的两个甲基落在同一个栅格点上,图上凭空出现一个三元环,三条
/// 边正好都是一个键长。全语料上这种重合占 6%,而且距离全是**正好 0**:布局
/// 走的是 30° 栅格上的单位步长,两条支路撞到同一个格点是系统性的,不是浮点抖动。
///
/// 按 30° 一档往两边试,与整张图的栅格一致;五档之内都腾不开就退回 `ideal`,
/// 交给消冲突,消不掉再如实报进 `unresolved`。
fn free_direction(
    center: Point2,
    ideal: f64,
    occupied: &[f64],
    taken: &[Point2],
    drawn: &[(Point2, Point2)],
) -> f64 {
    const STEP: f64 = std::f64::consts::FRAC_PI_6;
    /// 多近算重合。取键长的十分之一 —— 真正分得开的两个位置至少差半个键长。
    const TOL: f64 = 0.1;
    let at = |t: f64| center + Point2::new(BOND_LEN, 0.0).rotated(t);
    let clear = |t: f64| {
        let p = at(t);
        !taken.iter().any(|q| p.dist(*q) < TOL)
    };
    // 挪出来的方向落在哪一侧,决定键角是变宽还是变窄。**60° 不只是难看** ——
    // 链上出现一个 60° 的拐角,看着像旁边有个三元环,那是让人读错结构。
    //
    // 但**不能因此把窄的那一侧拒掉**:被拒的那个方向往往是唯一不撞的,拒了
    // 就换来一处碰撞。实测硬拒的代价是未解冲突 +494、干净率 −2.8 个百分点,
    // 只换来 107 处窄角 —— 亏的。
    //
    // 所以只调顺序不拒绝:同样偏离一档时,先试角度更宽的那一侧。
    //
    // # 试过"越接近 120° 越好",亏了
    //
    // "越宽越好"有个副作用:180° 正是最宽的,于是新键与已有的键连成一条直线,
    // 那个二度原子在图上**根本看不见**(顶点处没有拐角)。实测全量语料 148 张
    // 图(0.8%)因此出现"骨架原子被摆成 180°"。
    //
    // 于是试过把排序键换成 `|最窄夹角 − 理想值|`,让 60° 与 180° 同等地躲。
    // **全量语料上是亏的**:
    //
    // | | 越宽越好 | 越接近 120° |
    // |---|---:|---:|
    // | 骨架原子 180° | 148 | **140**(−8) |
    // | 键角不过窄 | **288** | 356(+68) |
    // | 未解冲突 | **1161** | 1199(+38) |
    // | 写法无关 | **257** | 260(+3) |
    // | 干净 | **91.5%** | 91.3% |
    //
    // 拿 8 处共线换 68 处窄角加 38 处冲突,不划算 —— "更宽"同时也意味着"离
    // 别的原子更远",这才是它压住碰撞的原因。那 148 处共线由渲染那边补符号
    // 兜底(见 `render::is_collinear`),并如实记在审计的质量分档里。
    let narrowest = |t: f64| {
        occupied
            .iter()
            .map(|o| {
                let d = (t - o).rem_euclid(std::f64::consts::TAU);
                d.min(std::f64::consts::TAU - d)
            })
            .fold(std::f64::consts::TAU, f64::min)
    };
    // 新键与已画的键交叉。共端点不算 —— 那是相邻的键,`segments_cross` 已经放过。
    let uncrossed = |t: f64| {
        let p = at(t);
        !drawn.iter().any(|(u, v)| segments_cross(center, p, *u, *v))
    };

    // 候选:理想方向,然后按 30° 一档往两边铺开。同一档里角度宽的排前面。
    //
    // # 试过在这里加"对侧那个同样理想的位置",没用
    //
    // 想法是:理想位置被占时,先试它关于已占方向的镜像(仍是精确的理想角),
    // 再考虑偏离。**变异验证说它不吃劲** —— 去掉之后角度判据照样绿,而全量
    // 语料上去掉它反而更好(窄角 209 → 180、交叉 83 → 78、干净 +14)。
    //
    // 原因是"对侧那个理想位置"通常正被兄弟取代基占着。真正管用的是**上游**
    // 那步:在 `place_neighbours` 里比较两侧的拥挤度、整条臂换边。
    let mut ranked: Vec<(u32, i64, f64)> = vec![(0, 0, ideal)];
    for k in 1..=5u32 {
        for sign in [1.0, -1.0] {
            let t = ideal + STEP * f64::from(k) * sign;
            #[allow(clippy::cast_possible_truncation)]
            let wide = -(narrowest(t) * 1e6).round() as i64; // 取负 → 宽的排前
            ranked.push((k, wide, t));
        }
    }
    ranked.sort_by_key(|c| (c.0, c.1));
    let cands: Vec<f64> = ranked.into_iter().map(|c| c.2).collect();

    // 两轮:先要"既不重合也不交叉",都腾不开就退而只求"不重合"。
    // **重合排在交叉前面** —— 重合会凭空造出一个假环,交叉只是难读。
    if let Some(t) = cands.iter().find(|t| clear(**t) && uncrossed(**t)) {
        return *t;
    }
    cands.iter().copied().find(|t| clear(*t)).unwrap_or(ideal)
}

/// 一个原子周围相邻两根键的理想夹角(弧度)。
///
/// # sp 的原子要画成直线
///
/// 只看度数是不够的:氰基的碳、炔碳、累积双键的中心碳都是 **sp 杂化,键角
/// 180°**,而它们的度数是 2 —— 按度数给 120° 的话,`R—C≡N` 会画成折的。
/// 这不是好看不好看的问题,是**画错了**:读者从图上读到的键角与分子的实际
/// 几何不符,而线条本身看着一点毛病没有。
///
/// 判据是键级不是度数:有三键,或者有两根双键(累积双键),就是 sp。
fn ideal_angle(mol: &MolBuilder, a: u32, style: &Style) -> f64 {
    let mut doubles = 0usize;
    let mut triple = false;
    for (_, bi) in mol.neighbors(a) {
        match mol.bonds()[bi as usize].order {
            BondOrder::Triple => triple = true,
            BondOrder::Double => doubles += 1,
            _ => {}
        }
    }
    if triple || doubles >= 2 {
        return std::f64::consts::PI;
    }
    let degree = mol.degree(a);
    if degree <= 3 {
        style.chain_angle_deg.to_radians()
    } else {
        std::f64::consts::TAU / degree as f64
    }
}

/// 把 `n` 个新方向分配到已占方向 `occupied`(已排序)之外的空隙里。
fn allocate(occupied: &[f64], n: usize, ideal: f64, zig: i8) -> Vec<f64> {
    let sign = if zig >= 0 { 1.0 } else { -1.0 };

    match occupied.len() {
        // 起点:从固定角度铺开
        0 => (0..n).map(|k| SEED_ANGLE + ideal * k as f64).collect(),

        // 链上一步:**不能取扇区中心**,那是 180°,会把链拉直
        1 => {
            let base = occupied[0];
            (0..n)
                .map(|k| {
                    // 第一个取 ±ideal,其余向另一侧交替铺开
                    let step = (k as f64 / 2.0).floor() + 1.0;
                    let s = if k % 2 == 0 { sign } else { -sign };
                    base + s * ideal * step
                })
                .collect()
        }

        // 已有两个以上方向:找最大空隙,在里面均分
        _ => {
            let (start, gap) = largest_gap(occupied);
            // n 个新方向把空隙切成 n+1 份
            (0..n)
                .map(|k| start + gap * (k as f64 + 1.0) / (n as f64 + 1.0))
                .collect()
        }
    }
}

/// 已排序角度序列中最大的空隙:返回(空隙起始角, 空隙大小)。
///
/// # 比大小必须先量化,平局必须有绝对判据
///
/// 三个已占方向恰好各差 120° 时(稠环上挂一个取代基就是这样),**三个空隙在
/// 数学上精确相等**。先前直接拿 `>` 比浮点,于是"哪个最大"由末位决定 —— 而
/// 末位取决于环坐标是按什么次序算出来的,同一个分子换种写法就换一个答案。
///
/// 实测:`C1CCN2[C@@H](C1)C=CC3=C2CCCC3=O` 的稠合碳上,三个已占方向算出来是
///
/// ```text
/// 写法 A: -2.09439510239319571, -0.00000000000000067, 2.09439510239319526
/// 写法 B: -2.09439510239319615, -0.00000000000000067, 2.09439510239319526
/// ```
///
/// 只差 4.4e-16,而补出来的那个氢因此挂到了 **120° 外的另一个扇区**。这一处
/// 是"写法无关"违例里相当大的一块 —— 它不是布局挑错了,是根本没在挑。
///
/// 所以:空隙量化到 1e-9 再比;仍然并列时取**起始角最小**的那个扇区 —— 那是
/// 与写法无关的绝对判据,与本文件里 `occ` 的排序、`mitre_end` 的量化同一个路子。
fn largest_gap(sorted: &[f64]) -> (f64, f64) {
    /// 量化的刻度。真正不等的两个空隙至少差 30°(栅格步长),而浮点噪声在
    /// 1e-15 量级 —— 中间空得很,取哪个数量级都一样。
    const QUANT: f64 = 1e9;
    #[allow(clippy::cast_possible_truncation)]
    let q = |x: f64| (x * QUANT).round() as i64;

    let n = sorted.len();
    debug_assert!(n >= 2);
    let mut cands: Vec<(i64, i64, f64, f64)> = Vec::with_capacity(n);
    let wrap_start = sorted[n - 1];
    let wrap = sorted[0] + std::f64::consts::TAU - wrap_start;
    cands.push((q(wrap), q(wrap_start), wrap_start, wrap));
    for i in 0..n - 1 {
        let g = sorted[i + 1] - sorted[i];
        cands.push((q(g), q(sorted[i]), sorted[i], g));
    }
    // 空隙大的在前;并列取起始角最小的
    cands.sort_by_key(|c| (std::cmp::Reverse(c.0), c.1));
    let c = cands[0];
    (c.2, c.3)
}

#[cfg(test)]
mod tests {
    use super::largest_gap;

    /// splitmix64 + Fisher–Yates。仿射式的"置换"搅不动东西 —— 审计里记过那个坑。
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
            let j = (next() % (i as u64 + 1)) as usize;
            v.swap(i, j);
        }
        v
    }

    #[test]
    fn the_same_direction_always_gets_the_same_angle() {
        // `angle()` 走 `atan2`,值域 `(-π, π]` —— **−180° 与 +180° 是同一个
        // 方向,却排在序列的两头**。末位差 4.4e-16 就足以决定它落在哪一端,
        // 于是同一组已占方向在两种写法下排出两个不同的序列,`largest_gap` 看到
        // 的空隙序列跟着不同,取代基差 120°。
        //
        // 实测这个分子:
        //
        // ```text
        // 写法 A: occ = [-3.14159265358979312, -1.04719755119659808, 1.04719755119659763]
        // 写法 B: occ = [-1.04719755119659808,  1.04719755119659763, 3.14159265358979267]
        // ```
        //
        // 化到 `[0, 2π)` 之后两边都是 `[1.047, 3.142, 5.236]`。全量语料上这一处
        // 让写法无关违例从 **77 降到 23**。
        // 前两个踩的是 ±π 那个断点,第三个踩的是 **0/2π** 那个 —— 化到
        // `[0, 2π)` 只把断点挪了个地方,贴着 2π 的角要再掐回 0。
        for smi in [
            "C[C]1(CCC[C]2(C)[CH]1CCC3=C2C=C(O)C=C3)C(O)=O",
            "C[C]1(CC[CH]2C(=C1)CC[CH]3[C]2(C)CCC[C]3(C)C(O)=O)C=C",
            "CC(C)C1=CC[CH]2C(=C1)CC[CH]3[C]2(C)CCC[C]3(C)C(O)=O",
        ] {
            let mut m = crate::tests_prep(smi);
            omgkit_io::stereo::perceive_bond_stereo(&mut m);
            let ranks = omgkit_io::canon::canonical_ranks(&m);
            let fp = |x: &MolBuilder, r: &[u32]| {
                let c = crate::generate(x, &crate::style::Style::ACS_1996).coords;
                let mut v: Vec<(u32, i64, i64)> = (0..c.len())
                    .map(|i| {
                        (
                            r[i],
                            (c[i].x * 1e4).round() as i64,
                            (c[i].y * 1e4).round() as i64,
                        )
                    })
                    .collect();
                v.sort_unstable();
                v
            };
            let want = fp(&m, &ranks);
            let mut compared = 0usize;
            for seed in 0..16u64 {
                let w = omgkit_io::smiles::write_with_priority(&m, &shuffled(m.num_atoms(), seed));
                let Ok(mut m2) = omgkit_io::smiles::parse(&w.smiles) else {
                    continue;
                };
                if omgkit_chem::pipeline::sanitize(&mut m2).is_err() {
                    continue;
                }
                omgkit_io::stereo::perceive_bond_stereo(&mut m2);
                if omgkit_io::canon::canonical_smiles(&m2).smiles
                    != omgkit_io::canon::canonical_smiles(&m).smiles
                {
                    continue;
                }
                let r2 = omgkit_io::canon::canonical_ranks(&m2);
                assert_eq!(
                    fp(&m2, &r2),
                    want,
                    "{smi}:换成 {} 之后摆得不一样了",
                    w.smiles
                );
                compared += 1;
            }
            assert!(compared > 0, "{smi}:一次都没比成 —— 判据空过了");
        }
    }

    #[test]
    fn a_three_way_tie_of_gaps_is_not_broken_by_the_last_bit() {
        // 稠环上的取代基:三个已占方向恰好各差 120°,**三个空隙精确相等**。
        // 先前拿 `>` 直接比浮点,谁"最大"由末位决定 —— 而末位取决于环坐标是按
        // 什么次序算出来的,同一个分子换种写法就换一个扇区,取代基差 120°。
        //
        // 实测那两组数只差 4.4e-16(见 `largest_gap` 的文档)。这里把那个量级
        // 的扰动加在每一个位置上,结果必须一个样。
        let base = [
            -2.094_395_102_393_195_7_f64,
            -0.000_000_000_000_000_67,
            2.094_395_102_393_195_3,
        ];
        let want = largest_gap(&base);
        for i in 0..3 {
            for eps in [-4.4e-16, 4.4e-16, -1e-15, 1e-15] {
                let mut v = base;
                v[i] += eps;
                v.sort_by(|a, b| a.partial_cmp(b).expect("非 NaN"));
                let got = largest_gap(&v);
                assert!(
                    (got.0 - want.0).abs() < 1e-9,
                    "第 {i} 个方向抖动 {eps:e} 之后挑了另一个扇区:{:.6} → {:.6}",
                    want.0,
                    got.0
                );
            }
        }
    }

    #[test]
    fn a_genuinely_larger_gap_still_wins() {
        // 上一条只说"平局要稳",不能顺手把"真的更大"也压掉 —— 那样就成了
        // "永远取第一个扇区"。
        let v = [0.0_f64, 1.0, 1.2];
        let (start, gap) = largest_gap(&v);
        assert!(
            (start - 1.2).abs() < 1e-9,
            "该取 1.2 起那个最大的空隙,实得起点 {start:.4}"
        );
        assert!((gap - (std::f64::consts::TAU - 1.2)).abs() < 1e-9);
    }

    #[test]
    fn an_arm_hanging_off_a_ring_keeps_its_ideal_angles() {
        // `allocate` 在"只有一个已占方向"时按锯齿的符号取 ±理想角,**那个符号
        // 不看旁边有没有东西**。挑错一侧,挂在环上的臂就朝环卷回去,臂上的
        // 取代基撞到环,只能按 30° 一档挪 —— 挪出来就是 120∓30 = 90° 或 150°。
        //
        // 实测:阿司匹林(ChemDraw 规范)乙酰基那个 sp² 碳,三个角是
        // **90 / 120 / 150**,三个都该是 120°;整条乙酰基臂折回来贴着苯环。
        //
        // 修法是在偏离理想角**之前**先试"对侧"那个同样理想的位置(把理想方向
        // 关于已占方向镜像,镜像保角)。这条判据守的就是"能不歪就不歪"。
        for smi in [
            "CC(=O)Oc1ccccc1C(=O)O", // 阿司匹林
            "CC(=O)Nc1ccc(O)cc1",    // 扑热息痛
            "CC(=O)Oc1ccccc1",       // 乙酸苯酯
            "COc1ccccc1OC(C)=O",     // 两个取代基挤在邻位
        ] {
            for style in &Style::ALL {
                let mut m = omgkit_io::smiles::parse(smi).unwrap();
                omgkit_chem::pipeline::sanitize(&mut m).unwrap();
                omgkit_io::stereo::perceive_bond_stereo(&mut m);
                let d = crate::generate(&m, style);
                for a in 0..u32::try_from(m.num_atoms()).unwrap() {
                    let n: Vec<u32> = m.neighbors(a).map(|(x, _)| x).collect();
                    if n.len() < 2 {
                        continue;
                    }
                    let c = d.coords[a as usize];
                    for i in 0..n.len() {
                        for j in (i + 1)..n.len() {
                            let u = (d.coords[n[i] as usize] - c).normalized();
                            let v = (d.coords[n[j] as usize] - c).normalized();
                            let deg = u.dot(v).clamp(-1.0, 1.0).acos().to_degrees();
                            // 允许的角是**这个原子自己的理想角的整数倍**:度 3
                            // 只许 120,度 4 许 90 与 180,sp 只许 180。
                            // 拿一张白名单(120/180/90)是不行的 —— 度 3 的
                            // 原子上 90° 会被放过,而那正是要抓的毛病。
                            let ideal = ideal_angle(&m, a, style).to_degrees();
                            let ok = (1..=6)
                                .map(|k| ideal * f64::from(k))
                                .take_while(|t| *t <= 180.5)
                                .any(|t| (deg - t).abs() < 1.0);
                            assert!(
                                ok,
                                "[{}] {smi}:{}-{a}-{} 的夹角是 {deg:.1}°,不是标准角 —— \
                                 理想位置被占时该先试对侧,而不是按 30° 一档歪",
                                style.name, n[i], n[j]
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn avoiding_a_taken_spot_does_not_pinch_the_angle_to_sixty_degrees() {
        // 位置被占了要挪,而挪的方向落在哪一侧决定键角是变宽还是变窄。
        // **60° 不只是难看** —— 链上出现一个 60° 的拐角,看着像旁边有个三元环,
        // 那是让人读错结构。同样偏离一档时先试宽的那一侧就能躲开。
        //
        // 实测:氮芥的两条 `N—CH₂—CH₂—Cl` 臂上量到过 60.1°。
        for smi in [
            "CC(CCCN(CCCl)CCCl)NC1=C2C=CC(=CC2=NC=C1)Cl",
            "ClCCN(CCCl)CCCl",
            "CC(C)(C)CC(C)(C)C",
        ] {
            for style in &Style::ALL {
                let mut m = omgkit_io::smiles::parse(smi).unwrap();
                omgkit_chem::pipeline::sanitize(&mut m).unwrap();
                let d = crate::generate(&m, style);
                for a in 0..u32::try_from(m.num_atoms()).unwrap() {
                    let nbrs: Vec<u32> = m.neighbors(a).map(|(n, _)| n).collect();
                    if nbrs.len() != 2 {
                        continue;
                    }
                    let c = d.coords[a as usize];
                    let u = (d.coords[nbrs[0] as usize] - c).normalized();
                    let v = (d.coords[nbrs[1] as usize] - c).normalized();
                    let deg = u.dot(v).clamp(-1.0, 1.0).acos().to_degrees();
                    assert!(
                        deg > 89.0,
                        "[{}] {smi}:{}–{a}–{} 的夹角只有 {deg:.1}°",
                        style.name,
                        nbrs[0],
                        nbrs[1]
                    );
                }
            }
        }
    }

    #[test]
    fn no_two_atoms_are_drawn_on_the_same_point() {
        // **重合是最糟的一种画错。** 两个原子叠在一起时,它们各自的键首尾相接,
        // 图上就多出一个分子里没有的环 —— 读者没有任何办法看出那个环是假的。
        //
        // 实测:下面第一个三萜的两个甲基落在同一个栅格点上,图上凭空出现一个
        // 三元环,三条边正好都是一个键长。全语料上这种重合曾占 6%,而且距离
        // 全是**正好 0** —— 布局走的是 30° 栅格上的单位步长,两条支路撞到同
        // 一个格点是系统性的,不是浮点抖动。
        for smi in [
            "CC([CH]1CC[C]2(CC[C]3(C)[C]4(C)[CH](CC[CH]3[CH]12)[C]1(C)[CH](CC4)C([CH](CC1)O)(C)C)CO)=C",
            "[O-][N+](=O)C1=CC(=CC=C1Cl)S(=O)(=O)C2=CC=C(Cl)C(=C2)[N+]([O-])=O",
            "CC(C)(C)c1ccccc1",
            "CC(=O)Oc1ccccc1C(=O)O",
        ] {
            for style in &Style::ALL {
                let mut m = omgkit_io::smiles::parse(smi).unwrap();
                omgkit_chem::pipeline::sanitize(&mut m).unwrap();
                let d = crate::generate(&m, style);
                for i in 0..d.coords.len() {
                    for j in (i + 1)..d.coords.len() {
                        let dist = d.coords[i].dist(d.coords[j]);
                        assert!(
                            dist > 0.1,
                            "[{}] {smi}:原子 {i} 与 {j} 相距 {dist:.4} 个键长 —— 画在同一点上了",
                            style.name
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn an_sp_atom_is_drawn_straight() {
        // 氰基的碳、炔碳、累积双键的中心碳都是 sp 杂化,键角 **180°**。它们的
        // 度数是 2,只按度数给 120° 的话 `R—C≡N` 会画成折的 —— 那是画错了,
        // 而线条本身看不出毛病。
        for (smi, centre) in [
            ("CC#N", 1u32), // 乙腈:C1 是 sp
            ("CC#CC", 1),   // 2-丁炔
            ("CC=C=CC", 2), // 累积双键的中心碳
            ("N#CC(C)(C)C#N", 1),
        ] {
            let mut m = omgkit_io::smiles::parse(smi).unwrap();
            omgkit_chem::pipeline::sanitize(&mut m).unwrap();
            let d = crate::generate(&m, &Style::ACS_1996);
            let nbrs: Vec<u32> = m.neighbors(centre).map(|(n, _)| n).collect();
            assert!(nbrs.len() >= 2, "{smi}:中心该有两个邻居");
            let (p, q) = (d.coords[nbrs[0] as usize], d.coords[nbrs[1] as usize]);
            let c = d.coords[centre as usize];
            let (u, v) = ((p - c).normalized(), (q - c).normalized());
            let deg = u.dot(v).clamp(-1.0, 1.0).acos().to_degrees();
            assert!(
                (deg - 180.0).abs() < 1e-6,
                "{smi}:原子 {centre} 是 sp,键角却画成了 {deg:.1}°"
            );
        }
    }

    /// 判据里也要按规范秩打破平局 —— 与实现取同一个来源
    fn canonical(m: &MolBuilder) -> Vec<u32> {
        omgkit_io::canon::canonical_ranks(m)
    }

    use super::*;
    use crate::style::Style;

    const TOL: f64 = 1e-9;

    fn prep(smi: &str) -> MolBuilder {
        let mut m = omgkit_io::smiles::parse(smi).unwrap();
        omgkit_chem::pipeline::sanitize(&mut m).unwrap();
        m
    }

    /// 两个方向之间的夹角,取 [0, π]
    fn between(u: Point2, v: Point2) -> f64 {
        let c = u.normalized().dot(v.normalized()).clamp(-1.0, 1.0);
        c.acos()
    }

    #[test]
    fn a_chain_zigzags_instead_of_running_straight() {
        // 度数 2 的原子若取 180°,链就成了一条直线 —— 既不是化学惯例,
        // 也让后面的可旋转键无处可翻。这条守的正是那个 180°。
        let m = prep("CCCCC");
        let style = Style::ACS_1996;
        let mut pos: BTreeMap<u32, Point2> = BTreeMap::new();
        pos.insert(0, Point2::ORIGIN);
        let mut zig = 1i8;
        for a in 0..4u32 {
            let out = place_neighbours(&m, a, &pos, &[a + 1], &canonical(&m), &style, zig);
            pos.insert(out[0].atom, out[0].at);
            zig = out[0].zig;
        }
        assert_eq!(pos.len(), 5);
        for i in 1..4u32 {
            let ang = between(pos[&(i - 1)] - pos[&i], pos[&(i + 1)] - pos[&i]);
            assert!(
                (ang - 120f64.to_radians()).abs() < TOL,
                "第 {i} 个原子处的键角是 {:.1}°,应当是 120°",
                ang.to_degrees()
            );
        }
        // 锯齿:相邻两步的转向必须相反,否则会绕成圆弧
        let turn = |i: u32| (pos[&i] - pos[&(i - 1)]).cross(pos[&(i + 1)] - pos[&i]);
        assert!(turn(1) * turn(2) < 0.0, "第 1、2 步没有交替转向");
        assert!(turn(2) * turn(3) < 0.0, "第 2、3 步没有交替转向");
    }

    #[test]
    fn every_bond_is_one_unit_long() {
        let m = prep("CC(C)(C)C");
        let style = Style::ACS_1996;
        let mut pos: BTreeMap<u32, Point2> = BTreeMap::new();
        pos.insert(1, Point2::ORIGIN);
        let mut todo: Vec<u32> = m.neighbors(1).map(|(n, _)| n).collect();
        todo.sort_unstable();
        for p in place_neighbours(&m, 1, &pos, &todo, &canonical(&m), &style, 1) {
            pos.insert(p.atom, p.at);
        }
        for n in todo {
            let d = pos[&n].dist(pos[&1]);
            assert!((d - BOND_LEN).abs() < TOL, "键长 {d}");
        }
    }

    #[test]
    fn four_substituents_are_spread_not_stacked() {
        // 季碳:四根键必须分开。理想夹角在度数 4 时该退到 90°,若仍用 120°,
        // 四个方向只能铺满 360° 中的 360°—— 会有两个重叠。
        let m = prep("CC(C)(C)C");
        let style = Style::ACS_1996;
        let mut pos: BTreeMap<u32, Point2> = BTreeMap::new();
        pos.insert(1, Point2::ORIGIN);
        let mut todo: Vec<u32> = m.neighbors(1).map(|(n, _)| n).collect();
        todo.sort_unstable();
        assert_eq!(todo.len(), 4, "季碳应当有四个邻居");
        let out = place_neighbours(&m, 1, &pos, &todo, &canonical(&m), &style, 1);
        for i in 0..out.len() {
            for j in (i + 1)..out.len() {
                let ang = between(out[i].at, out[j].at);
                assert!(
                    ang > 45f64.to_radians(),
                    "第 {i}、{j} 个取代基只差 {:.1}°,挤在一起了",
                    ang.to_degrees()
                );
            }
        }
    }

    #[test]
    fn a_new_branch_goes_into_the_largest_free_sector() {
        // 已经占了两个方向时,新的必须落进**最大的空隙**。落进小空隙不会报错,
        // 只会让图上一边挤一边空。
        let m = prep("CC(C)C");
        let style = Style::ACS_1996;
        let mut pos: BTreeMap<u32, Point2> = BTreeMap::new();
        pos.insert(1, Point2::ORIGIN);
        // 手工把两个邻居摆在 0° 和 60°,留下一个 300° 的大空隙
        pos.insert(0, Point2::new(1.0, 0.0));
        pos.insert(2, Point2::new(0.5, 3f64.sqrt() / 2.0));
        let out = place_neighbours(&m, 1, &pos, &[3], &canonical(&m), &style, 1);
        let ang = out[0].at.angle().rem_euclid(std::f64::consts::TAU);
        // 大空隙是 60° → 360°,中点在 210°
        assert!(
            (ang - 210f64.to_radians()).abs() < 1e-6,
            "新支链落在 {:.1}°,应当落在最大空隙的中点 210°",
            ang.to_degrees()
        );
    }

    #[test]
    fn the_largest_gap_wraps_around_the_seam() {
        // 空隙搜索必须绕过 ±π 的接缝。漏掉环绕的那一段不会报错,只会在某些
        // 角度组合下把新键塞进一个其实很窄的缝里。
        //
        // **数据必须让环绕的那一段真的是最大空隙。** 第一版用的是
        // `[-3.0, -2.9, 3.0]`,那里最大的其实是 -2.9→3.0 这个**内部**空隙,
        // 于是把环绕逻辑破坏掉,这条照样绿 —— 名字说的和断言测的是两回事。
        // 三个方向挤在 0 附近,环绕那一段才是最大的。
        let sorted = vec![-0.1, 0.0, 0.1];
        let (start, gap) = largest_gap(&sorted);
        assert!(
            (start - 0.1).abs() < TOL,
            "空隙应当从最后一个方向 0.1 起,实得 {start}"
        );
        let want = -0.1 + std::f64::consts::TAU - 0.1;
        assert!((gap - want).abs() < TOL, "空隙大小应当是 {want},实得 {gap}");

        // 顺带守住"内部空隙也要能选中",免得改成只看环绕
        let inner = vec![-3.0, -2.9, 3.0];
        let (s2, g2) = largest_gap(&inner);
        assert!(
            (s2 - (-2.9)).abs() < TOL && (g2 - 5.9).abs() < TOL,
            "内部空隙没选对"
        );
    }
}
