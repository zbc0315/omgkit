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

/// 挑方向时一个候选有多坏:`(逐位重合的对数, 会不会压出窄角, 标签碰撞深度)`。
///
/// `(逐位重合的对数, 会不会压出窄角, 会不会与已画的键交叉, 标签碰撞深度)`。
/// 越小越好,按字典序比 —— 四位的次序就是这个项目立的严重度,见
/// [`free_direction`]。
type Cost = (usize, u8, u8, i64);

/// 一整块环系统的坐标:原子编号 → 位置。
pub(crate) type Block = BTreeMap<u32, Point2>;

/// 一次放置的结果。
pub(crate) struct Placed {
    /// 原子
    pub atom: u32,
    /// 坐标
    pub at: Point2,
    /// 它自己的锯齿符号,传给它的子代取反
    pub zig: i8,
    /// 这个邻居身后那**一整块环系统**的坐标 —— 挑方向时就已经算好的那一份。
    ///
    /// **调用方必须原样用,不许重算。** 重算会重新在两个镜像里挑,而挑的时候
    /// `pos` 已经变了(兄弟摆上去了),挑出来的可能是另一个镜像 —— 那么前瞻
    /// 累积进"已占"的就是一块**根本没被画出来的**坐标,后面每一步都建在假的
    /// 占位上。
    pub block: Option<Block>,
}

/// 布局那边传进来的上下文。打包成一个结构,是为了让 [`place_neighbours`] 的
/// 参数不超过 clippy 的上限。
///
/// 前四项是**整个分量共用**的;`blocks` 不是 —— 它每摆一个枢纽就重建一次,
/// 内容取决于那一刻还有哪些系统没摆。
pub(crate) struct Env<'a> {
    pub mol: &'a MolBuilder,
    pub ranks: &'a [u32],
    pub style: &'a Style,
    /// 碰撞半径,与 [`crate::refine::radii`] 同一套口径
    pub radii: &'a [f64],
    /// 成键的原子对 —— 它们本来就靠在一起,不算撞。与 `layout` 那边 `Around`
    /// 用的是**同一个集合**,免得"口径与 `refine` 一致"这句话变成两套实现。
    pub bonded: &'a std::collections::BTreeSet<(u32, u32)>,
    /// **坐标不是从理想 30° 栅格算出来的那些原子。**
    ///
    /// 指 [`crate::rings::layout_local`] **记了退化**的那些环系统里的原子 ——
    /// 查模板表命中、走弧法、跑松弛,三档都算。它们的坐标落在任意角度上,
    /// 于是那里的"最窄键角"量出来的不是"偏离理想几何多少",而是那一步碰巧
    /// 给了什么。见 `place_neighbours` 里 `pinch_floor` 的注释。
    ///
    /// **只装环系统自己的原子,不沿链传下去。** 传下去是错的:枢纽自己不在
    /// 那个系统里的话,它的**所有**已占方向都是 `free_direction`/`place_clear`
    /// 按 30° 档位挑出来的,相对角规整,地板照样有意义。实测传下去的代价:
    /// 非环枢纽上冒出 **21 张**图带着 30.00°/60.00° 的拐角(角度精确落在栅格上,
    /// 正说明那里的几何是规整的),而硬判据 `键角不过窄` 看不见它们 ——
    /// 那条判据只查干净的图,而离网必然意味着退化。
    pub off_grid: &'a std::collections::BTreeSet<u32>,
    /// 邻居 → 它身后那个**还没摆过的环系统**的局部坐标。
    ///
    /// 由调用方**预先算好**([`crate::rings::layout_local`] 很贵:要查模板表、
    /// 走弧法、必要时跑松弛),这里只做刚体变换。算好的那一份最后也要被真正
    /// 用上,见 [`Placed::block`]。
    pub blocks: &'a BTreeMap<u32, Block>,
}

/// 给 `a` 周围还没放置的邻居 `todo` 分配坐标。
///
/// `todo` 必须**已按规范秩排好**;顺序决定谁分到哪个方向,拿存储下标排就会
/// 引入写法依赖。
pub(crate) fn place_neighbours(
    env: &Env<'_>,
    a: u32,
    pos: &BTreeMap<u32, Point2>,
    todo: &[u32],
    zig: i8,
) -> Vec<Placed> {
    let (mol, ranks, style) = (env.mol, env.ranks, env.style);
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
    // **带上原子编号**:整块前瞻要按编号去查碰撞半径,也要跳过成键的那些对。
    let mut taken: Vec<(u32, Point2)> = pos.iter().map(|(k, v)| (*k, *v)).collect();
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
        // 拥挤度:新位置到每个已放好的原子的**平方反比**之和。RDKit 的 density
        // 是同一个想法但**不是同一个口径** —— 它累加的是 `1.0 / d`
        // (`EmbeddedFrag.cpp:855`,一次方)。这里用平方反比让近处的惩罚更陡。
        // **量化之后再比** —— 直接比浮点会让"分不出高下"取决于
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
    // **兄弟一摆上,就要算进"已占方向"。**
    //
    // `occupied` 是循环**开始前**的快照,而同一批兄弟摆下去只进 `taken`/`drawn`
    // —— 于是 `free_direction` 里的 `narrowest(t)`(既管"同一档先试更宽那侧"
    // 的排序,也是判"这个方向会不会把键角压窄"的唯一依据)**看不见刚摆下的
    // 兄弟**。而真正撞上的那些标签,八成正是同一个原子上的兄弟。
    let mut occ_live = occupied.clone();
    // **只在度数 ≤ 3 的枢纽上守 89°。** 与审计 `键角不过窄` 大体同一个适用范围
    // —— 度 4 的理想角本来就是 90°、度 6 是 60°,拿 89° 去拒它们等于拒掉正解。
    //
    // "大体"是因为还差两小条,都查过没有实际影响,记在这儿:
    // ① 审计还会跳过三元环内角,这里没跳。三元环整块作为环系统进来,`todo`
    //    里不会出现枢纽的同环邻居,所以走不到。
    // ② 这里的 `mol` 是**摘细后**的副本(η 配位只留一根代表键),审计量的是
    //    `drawn`。η5 的铁在这儿是度 2、在那儿是度 10。无害:摘细一定记了
    //    `HaptoCoordination` 退化,审计那一档只查干净的图。
    // **坐标不是从理想栅格算出来的枢纽,不守这一条。**
    //
    // 89° 这道地板是给理想几何设的:方向按 30° 一档铺开,挪两档就是 60°,
    // 而 60° 的拐角看着像旁边有个三元环。桥环那些系统的坐标(模板命中、弧法、
    // 松弛都算)不在那个栅格上,那里的"最窄角"量出来的是那一步碰巧给了什么
    // —— 拿它去拒方向,买不到"别让人读错",却要付真实的交叉。
    //
    // 这与上面那条"度 4 的理想角本来就是 90°,拿 89° 去拒它们等于拒掉正解"
    // 是同一个道理:**地板只在它标定的那套几何里有意义**。
    let pinch_floor =
        (mol.degree(a) <= 3 && !env.off_grid.contains(&a)).then(|| 89f64.to_radians());
    for (&atom, theta) in todo.iter().zip(dirs) {
        let look = Lookahead {
            env,
            atom,
            local: env.blocks.get(&atom),
        };
        let (theta, block) =
            free_direction(center, theta, &occ_live, &taken, &drawn, &look, pinch_floor);
        occ_live.push(theta);
        let at = center + Point2::new(BOND_LEN, 0.0).rotated(theta);
        taken.push((atom, at));
        drawn.push((center, at));
        // **兄弟摆过的那一整块也要记进"已占"。**
        //
        // `place_neighbours` 是一次算完全部方向才返回的 —— 调用方要等它返回
        // 之后才往 `pos` 里插。所以不在这儿累积的话,同一个枢纽上"这个环撞
        // 那个环"从头到尾**看不见**。实测代价是键交叉 72 → 110(+53%)。
        if let Some(b) = &block {
            for (k, p) in b {
                if *k != atom {
                    taken.push((*k, *p));
                }
            }
            for bd in mol.bonds() {
                if let (Some(u), Some(v)) = (b.get(&bd.begin), b.get(&bd.end)) {
                    drawn.push((*u, *v));
                }
            }
        }
        out.push(Placed {
            atom,
            at,
            // 子代取反,直链就走出锯齿
            zig: -zig,
            block,
        });
    }
    out
}

/// 挑方向时要用的「这个邻居身后是不是挂着一整块」。
struct Lookahead<'a> {
    env: &'a Env<'a>,
    /// 正在挑方向的那个邻居 —— 它同时是那一块的锚点
    atom: u32,
    /// 那一块的局部坐标;没有块就是 `None`,一切退化成从前的行为
    local: Option<&'a Block>,
}

impl Lookahead<'_> {
    /// 这个邻居落在 `at` 时,它身后那一块会摆成什么样、有多坏。
    ///
    /// 返回 `(代价, 那一块的坐标)`。代价是 `(逐位重合的对数, 量化的碰撞深度)`,
    /// **重合对数排在深度前面** —— 理由见 [`free_direction`]。
    ///
    /// **没有块时也要打分** —— 那时这个原子自己就是只有一个点的块,见函数体。
    /// (先前这里返回零代价,于是"挑最小"退化成"取第一个";那正是这一轮改掉
    /// 的东西。)
    fn cost(
        &self,
        center: Point2,
        at: Point2,
        taken: &[(u32, Point2)],
    ) -> ((usize, i64), Option<Block>) {
        let Some(local) = self.local else {
            // **没有块时,这个原子自己就是一块**(只有一个点)。
            //
            // 先前这里直接返回零代价,于是挑档位那一步只查了"坐标重不重合"
            // (`clear`,阈值 0.1 个键长)与"键交不交叉" —— **标签会不会压到
            // 别人身上,一个字都没问**。
            let one: Block = [(self.atom, at)].into_iter().collect();
            return (block_cost(self.env, &one, taken), None);
        };
        let dir = (at - center).normalized();
        let mut best: Option<((usize, i64), Block)> = None;
        for cand in crate::rings::place_candidates(self.env.mol, local, self.atom, at, dir) {
            let c = block_cost(self.env, &cand, taken);
            // **不许裸比浮点。** 深度已经量化成 i64,平局留前一个;而
            // `place_candidates` 的返回序是定死的,所以这是个规范的选择。
            // 单环的两个镜像本来就是同一个点集,深度只差最后一位。
            let better = match &best {
                None => true,
                Some((old, _)) => c < *old,
            };
            if better {
                best = Some((c, cand));
            }
        }
        let (c, cand) = best.expect("`place_candidates` 恒返回两个候选");
        (c, Some(cand))
    }
}

/// 一块待放置的坐标压在已占部分上有多重:`(逐位重合的对数, 量化的碰撞深度)`。
///
/// 口径与 [`crate::refine`] 一致:半径来自标签,**成键的一对不算**(相邻原子
/// 本来就靠在一起,锚点与枢纽正是这样一对)。
fn block_cost(env: &Env<'_>, cand: &Block, taken: &[(u32, Point2)]) -> (usize, i64) {
    /// 多近算"画在同一点上" —— 与硬判据 `原子不重合` 同一个阈值。
    const SAME: f64 = 0.05;
    let mut same = 0usize;
    let mut parts: Vec<f64> = Vec::new();
    for (i, p) in cand {
        for (j, q) in taken {
            if i == j || env.bonded.contains(&((*i).min(*j), (*i).max(*j))) {
                continue;
            }
            let d = p.dist(*q);
            if d < SAME {
                same += 1;
            }
            let want = env.radii[*i as usize] + env.radii[*j as usize];
            if d < want {
                parts.push((want - d).powi(2));
            }
        }
    }
    // **先排序再求和。** 两种写法给的是同一个几何、同一个多重集,但迭代序是
    // 存储序 —— 不排的话和会差最后一位,而下面的平局判定就靠这一位。
    parts.sort_by(f64::total_cmp);
    let depth: f64 = parts.iter().sum();
    #[allow(clippy::cast_possible_truncation)]
    let q = (depth * 1e9).round() as i64;
    (same, q)
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
/// # 只看一个原子是不够的:它身后可能挂着一整个环系统
///
/// `clear` 问的是"**这个原子**落在这儿撞不撞",可这个邻居往往是一整块环系统的
/// 接口 —— 环摆下去占的是十来个格点,而这一步对此一无所知。后果是硬判据
/// `原子不重合` 上 57 处里的一大族:六配位金属挂吡啶,配体按 60° 分开,配体氮
/// 到金属一个键长、环心到氮又一个键长,于是**两个相邻环心相距正好 2、外接圆
/// 半径各 1 —— 精确相切**,而切点正落在栅格上,两个环各有一个邻位碳逐位重合。
///
/// 所以候选要连**整块**一起打分,见 [`Lookahead::cost`]。块的坐标是真算出来的
/// (`layout_local` 的结果做刚体变换),不是估计。
///
/// # 打分不能写成布尔的「整块不撞」
///
/// 试过。一票否决会拿"蹭一下"去换"精确重合":某个取代基为了躲开一处深度
/// 0.0538 的轻微重叠跳到别的档,**给后面的兄弟挖了坑**;轮到兄弟时一个候选都
/// 满足不了,于是掉进兜底那一档退回 `ideal`,而 `ideal` 恰恰是唯一逐位重合的
/// 那个。全量语料实测:重合 57 → 20 的同时**新坏 12 张**,签名全是"度 4 枢纽、
/// 60°、环心距 2.000",与被修好的那族是同一个几何。
///
/// 改成"在每一轮里挑**最不坏**的"就没有这条兜底路了:重合 57 → 8,新坏 0。
///
/// 排序键里**重合对数必须排在深度前面**,而且这一位是吃劲的:只按深度跑一遍
/// 全量语料,`原子不重合` **8 → 10**(坏的是语料第 3469、3472 行,两个吡啶
/// 配合物)。道理是一次精确相切的深度只有 0.25 —— 比"挤到三个原子上"那种和
/// 还小,光看深度会**主动选中相切**。
fn free_direction(
    center: Point2,
    ideal: f64,
    occupied: &[f64],
    taken: &[(u32, Point2)],
    drawn: &[(Point2, Point2)],
    look: &Lookahead<'_>,
    // 低于这个夹角就算"把键角压窄了";`None` = 这个枢纽不守这一条
    pinch_floor: Option<f64>,
) -> (f64, Option<Block>) {
    const STEP: f64 = std::f64::consts::FRAC_PI_6;
    /// 多近算重合。取键长的十分之一 —— 真正分得开的两个位置至少差半个键长。
    const TOL: f64 = 0.1;
    let at = |t: f64| center + Point2::new(BOND_LEN, 0.0).rotated(t);
    let clear = |t: f64| {
        let p = at(t);
        !taken.iter().any(|(_, q)| p.dist(*q) < TOL)
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
    // 初值取 π:`d.min(TAU − d)` 的值域就是 `[0, π]`,`occupied` 为空(枢纽还
    // 没有任何已放邻居)时压根不存在键角。**别用 TAU** —— 那是个超出值域的
    // 哨兵,今天两个消费方都只问"是不是够小"所以无害,但将来若有人加一条
    // "别太宽/别 180°",TAU 会静悄悄地通过。
    let narrowest = |t: f64| {
        occupied
            .iter()
            .map(|o| {
                let d = (t - o).rem_euclid(std::f64::consts::TAU);
                d.min(std::f64::consts::TAU - d)
            })
            .fold(std::f64::consts::PI, f64::min)
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

    // **一轮。** 只有"落点不与已放原子重合"是硬门槛 —— 那会凭空造出一个假环,
    // 读者没有任何办法看出它是假的。其余全进排序键,按严重度排。
    //
    // **先前这里是两轮**:第一轮要"既不重合也不交叉",腾不开才退到只要"不重合"。
    // 那等于把"不交叉"当成硬门槛、排在"不压窄键角"**前面** —— 与本库的严重度
    // 次序正好相反。实测:第一轮只要有解第二轮就不跑,于是第 4134 行画出了
    // 60° 的拐角,而第二轮里明明有 `(same=0, pinch=0)` 的候选。
    //
    // 不是"取第一个",是**取最不坏的那个** —— 挂着环系统时算整块,
    // 没挂时算这个原子自己的标签(见 `Lookahead::cost`)。
    let pick = |pred: &dyn Fn(f64) -> bool| -> Option<(f64, Option<Block>)> {
        let mut best: Option<(Cost, f64, Option<Block>)> = None;
        for &t in cands.iter().filter(|t| pred(**t)) {
            let ((same, depth), block) = look.cost(center, at(t), taken);
            // **四位,次序就是这个项目自己立的严重度**:
            // ① 逐位重合 —— 图上凭空多一个环,读者没办法看出它是假的
            // ② 把键角压到 89° 以下 —— 60° 的拐角看着像旁边有个三元环,读错结构
            // ③ 与已画的键交叉 —— 看得见的丑,但信息没错
            // ④ 标签压上去 —— 难看,但读者知道那里有两个原子
            //
            // **③ 那一位是唯一没量化的。** 它是布尔的,由 `segments_cross` 用
            // 1e-9 的带宽判side;坐标在 30° 栅格上走单位步长,噪声在 1e-16 量级,
            // 而刚体变换保号 —— 结构上安全。真要哪天不安全了,断层会在这里。
            // **量化之后再比。** 这一位排在 `depth` 前面,翻它翻的不是末位,
            // 是整条支路重画;而裸比浮点会让"够不够 89°"取决于算到这一步的
            // 运算次序 —— 本文件里 `occ` 排序、`largest_gap`、`mitre_end` 都在
            // 这上面栽过。实测语料里 `|narrowest − 89°|` 落在 1e-3 弧度内的
            // 候选是 0 个,但松弛出来的桥环坐标不在 30° 栅格上,边界不是
            // 结构性安全的。
            #[allow(clippy::cast_possible_truncation)]
            let q = |x: f64| (x * 1e6).round() as i64;
            let pinch = u8::from(pinch_floor.is_some_and(|f| q(narrowest(t)) < q(f)));
            let key = (same, pinch, u8::from(!uncrossed(t)), depth);
            let better = match &best {
                None => true,
                Some((old, _, _)) => key < *old,
            };
            if better {
                let done = key == (0, 0, 0, 0);
                best = Some((key, t, block));
                if done {
                    break; // 一点不坏,后面不可能更好
                }
            }
        }
        best.map(|(_, t, block)| (t, block))
    };
    if let Some(hit) = pick(&|t| clear(t)) {
        return hit;
    }
    // 全都腾不开:退回理想方向,块也照这个方向摆
    let (_, block) = look.cost(center, at(ideal), taken);
    (ideal, block)
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
mod tests_ordering {
    use crate::style::Style;

    fn prep(smi: &str) -> omgkit_core::MolBuilder {
        let mut m = omgkit_io::smiles::parse(smi).unwrap();
        omgkit_chem::pipeline::sanitize(&mut m).unwrap();
        omgkit_io::stereo::perceive_bond_stereo(&mut m);
        m
    }

    /// 一张图上最窄的键角(度)。只看度数 2..=3 的原子、跳过三元环内角 ——
    /// 与硬判据 `no_angle_is_pinched` 同一个口径。
    fn narrowest(smi: &str, style: &Style) -> (f64, usize) {
        let m = prep(smi);
        let d = crate::generate(&m, style);
        let g = d.drawn(&m);
        let mut worst = 180.0f64;
        for a in 0..u32::try_from(g.num_atoms()).expect("原子数超出 u32") {
            let nb: Vec<u32> = g.neighbors(a).map(|(n, _)| n).collect();
            if !(2..=3).contains(&nb.len()) {
                continue;
            }
            for i in 0..nb.len() {
                for j in (i + 1)..nb.len() {
                    if g.neighbors(nb[i]).any(|(n, _)| n == nb[j]) {
                        continue;
                    }
                    let u = (d.coords[nb[i] as usize] - d.coords[a as usize]).normalized();
                    let v = (d.coords[nb[j] as usize] - d.coords[a as usize]).normalized();
                    worst = worst.min(u.dot(v).clamp(-1.0, 1.0).acos().to_degrees());
                }
            }
        }
        (worst, d.crossings.len())
    }

    /// **压窄键角比画出交叉更严重** —— 挑方向那一步的排序键必须这么排。
    ///
    /// 语料第 4134 行(一个季碳挂两条仲丁基加羧基)。先前「不交叉」是硬门槛、
    /// 排在「不压窄键角」前面,于是这张图两套规范都在原子 10 处画出 **60.0°**
    /// 的拐角 —— 那看着像旁边有个三元环,是**读错结构**。把交叉降到排序键的
    /// 第三位之后是 90.0°,而且这张图**一处交叉都没多付**。
    ///
    /// 变异:把 `crossed` 挪到 `pinch` 前面(或挪到 `same` 前面),当场红。
    #[test]
    fn a_pinched_angle_outranks_a_bond_crossing() {
        for style in &Style::ALL {
            let (worst, cross) = narrowest("CCC(CC)C(O)(C(CC)CC)C(O)=O", style);
            assert!(
                worst >= 89.0,
                "[{}] 键角被压到 {worst:.1}° —— 60° 的拐角看着像个三元环",
                style.name
            );
            assert_eq!(cross, 0, "[{}] 这张图本不该有交叉", style.name);
        }
    }

    /// 但交叉**仍然要算** —— 它只是排在窄角后面,不是不管了。
    ///
    /// 语料第 2284 行。把 `crossed` 从排序键里整个去掉,这张图会多出一处交叉
    /// (实测全语料 `有键交叉` 113 → **208**)。
    #[test]
    fn a_bond_crossing_still_counts_it_is_only_ranked_lower() {
        for style in &Style::ALL {
            let (_, cross) = narrowest(
                "c1cc(oc1)C2=[N+]([C@@H]3CCCC[C@@H]3[N+](=C2)[O-])[O-]",
                style,
            );
            assert_eq!(
                cross, 0,
                "[{}] 交叉那一位从排序键里掉了 —— 这张图多出了交叉",
                style.name
            );
        }
    }

    /// **松弛出来的坐标上不守 89° 那道地板。**
    ///
    /// 语料第 780 行,一个桥环缩醛内酯。它的环系统是松弛出来的,坐标不在 30°
    /// 栅格上 —— 那里的"最窄角"量出来的是松弛器碰巧给了什么,拿它去拒方向
    /// 买不到"别让人读错",却要付真实的交叉。
    ///
    /// 变异:把 `pinch_floor` 上那个 `!env.off_grid.contains(&a)` 去掉,
    /// 这张图当场从 0 处交叉变成 1 处。全语料一起看:`有键交叉` 113 → **128**,
    /// 而 `只是标着退化,没有别的毛病` 290 → **275** —— 15 张原本无瑕的桥环图
    /// 有了看得见的交叉,那正是这道门挡住的东西。
    #[test]
    fn a_relaxed_ring_system_does_not_get_the_ninety_degree_floor() {
        for style in &Style::ALL {
            let smi = "C[C@@H]1[C@H]2[C@H]3C[C@@H](O1)O[C@@H]2OC=C3C(=O)OC";
            let m = prep(smi);
            let d = crate::generate(&m, style);
            // 前提:这张图的布局**确实**退化了,否则这条判据验的不是它想验的东西
            assert!(
                !d.degraded.is_empty(),
                "[{}] 这个分子的布局没退化 —— 选错例子了",
                style.name
            );
            assert_eq!(
                d.crossings.len(),
                0,
                "[{}] 松弛出来的坐标上守了 89° 地板,换来一处交叉",
                style.name
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::largest_gap;

    /// 这个分子里,有没有一个原子从**外面**接着好几个独立的环系统。
    ///
    /// 判据的前提要自己成立:接不着好几个环系统的分子,根本走不到"挑方向时
    /// 要不要把整块算进去"那一步,判据就空过了。
    fn rings_hanging_off_one_atom(m: &omgkit_core::MolBuilder) -> usize {
        let systems = omgkit_chem::rings::fused_ring_systems(m);
        (0..u32::try_from(m.num_atoms()).expect("原子数超出 u32"))
            .filter(|a| !systems.iter().any(|s| s.contains(a)))
            .map(|a| {
                systems
                    .iter()
                    .filter(|s| m.neighbors(a).any(|(n, _)| s.contains(&n)))
                    .count()
            })
            .max()
            .unwrap_or(0)
    }

    fn drawn_without_overlap(smi: &str, want_rings: usize) {
        let mut m = omgkit_io::smiles::parse(smi).expect("SMILES 该能解析");
        omgkit_chem::pipeline::sanitize(&mut m).expect("该能 sanitize");
        let n = rings_hanging_off_one_atom(&m);
        assert!(
            n >= want_rings,
            "{smi} 只有一个原子外接 {n} 个环系,少于 {want_rings} —— 这条判据空过了"
        );
        for style in &crate::style::Style::ALL {
            let d = crate::generate(&m, style);
            for i in 0..d.coords.len() {
                for j in (i + 1)..d.coords.len() {
                    let dist = d.coords[i].dist(d.coords[j]);
                    assert!(
                        dist >= 0.05,
                        "{}:原子 {i} 与 {j} 相距 {dist:.4} 个键长",
                        style.name
                    );
                }
            }
        }
    }

    /// 挂在同一个原子上的两个标签,不许压在一起。
    ///
    /// 分子取自真语料(第 4000 行):对甲苯磺酰二氯胺,磺酰硫上挂着两个 `O`
    /// 和一个 `N`。
    ///
    /// **消冲突够不着这一类。** 全量语料上"两个标签盒真的叠上"的 520 对里,
    /// **425 对(81.7%)拓扑距离是 2** —— 两个原子挂在同一个原子上,夹角由
    /// 这个文件定死,而 `refine` 的算子全是等距变换(绕某根键镜像),**动不了
    /// 它们的相对位置**。所以这一步不管,就没人管了。
    ///
    /// 先前挑方向只查两件事:坐标重不重合(`clear`,阈值 0.1 个键长)、新键
    /// 交不交叉。**标签会不会压到别人身上,一个字都没问。**
    ///
    /// 变异:把 [`Lookahead::cost`] 里 `local` 为 `None` 那一支换回 `((0, 0), None)`
    /// (单原子不打分)→ 这条当场红,报的是**第 6226 行那个双甲磺酸酯**的
    /// `原子 2 与 3 挂在 1 上,两个标签盒叠着`。
    ///
    /// **三个分子里只有第三个吃劲**,前两个在这个变异下也不叠 —— 它们是回归
    /// 锚(磺酰、磺酰胺这一族别的摆法),不是主角。如实说,免得后来的人以为
    /// 换掉第一个也照样红。
    #[test]
    fn two_labels_on_the_same_atom_do_not_get_stacked_on_each_other() {
        for smi in [
            "Cc1ccc(cc1)S(=O)(=O)N(Cl)Cl",
            "CCOC(=O)CNS(=O)(=O)C1=CC=CC=C1",
            "CS(=O)(=O)OCCCCCCCCCOS(C)(=O)=O",
        ] {
            labels_do_not_overlap(smi);
        }
    }

    /// 躲标签的时候,不许把键角压到 89° 以下。
    ///
    /// 分子取自真语料(第 367 行):苯环上两个**间位**硝基,外加一个环己硫醚。
    ///
    /// 变异打红的那一处是 `原子 17 处 18–17–19` —— 17 是第二个硝基的氮、
    /// 18/19 是它自己那两个氧。**被压窄的是一个硝基内部的 O–N–O**,正是
    /// "同一个原子上两个标签"那一类,与两个硝基之间无关。
    ///
    /// 上一条那个"单原子也打分"单独上会**把好图画坏**:挑档位那一步在**所有
    /// 档位**里取代价最小,理想角那一档不再优先,于是为了让标签分开而偏两档
    /// —— 120° 就成了 60°。全量语料实测 168 处窄角,其中 **156 处正好 60.0°**,
    /// 而这个文件自己写着「60° 的拐角看着像旁边有个三元环,那是让人读错结构」。
    ///
    /// 所以代价要四位:`(逐位重合, 会不会压出窄角, 会不会交叉, 标签碰撞深度)` —— 正是这个
    /// 项目立的严重度次序。补上 `pinch` 那一位之后窄角 168 → **8**;后来把交叉
    /// 也收进这条键(排在 `pinch` 之后),窄角进一步归 **0**。
    ///
    /// 两个变异各红一次:
    /// - 把 `pick` 的排序键去掉 `pinch` 那一位 → 这条红。
    /// - 把 `occ_live` 换回循环开始前的 `occupied` 快照 → 也红。后者是**先前
    ///   就有的缺陷**:兄弟摆上去只进 `taken`,从不进 `occupied`,而 `narrowest`
    ///   只看 `occupied` —— 于是"这个方向会不会压窄"看不见刚摆下的兄弟,而真正
    ///   撞上的正是兄弟。
    #[test]
    fn dodging_a_label_must_not_pinch_the_bond_angle() {
        let smi = "[O-][N+](=O)C1=CC(=C(S[CH]2CCCC[CH]2Cl)C=C1)[N+]([O-])=O";
        let mut m = omgkit_io::smiles::parse(smi).expect("SMILES 该能解析");
        omgkit_chem::pipeline::sanitize(&mut m).expect("该能 sanitize");
        for style in &crate::style::Style::ALL {
            let d = crate::generate(&m, style);
            // **拿被画的那个分子。** 为画构型补的显式氢也在里面,而 `coords`
            // 的下标是相对它的 —— 拿原分子的度数会与适用范围对不上号。
            let grown = d.drawn(&m);
            let mut looked = 0usize;
            for a in 0..u32::try_from(grown.num_atoms()).expect("原子数超出 u32") {
                let nbrs: Vec<u32> = grown.neighbors(a).map(|(n, _)| n).collect();
                // 与审计 `键角不过窄` 同一个适用范围:度数 ≤ 3;三元环内角不算
                if nbrs.len() < 2 || nbrs.len() > 3 {
                    continue;
                }
                let c = d.coords[a as usize];
                for i in 0..nbrs.len() {
                    for j in (i + 1)..nbrs.len() {
                        if grown.neighbors(nbrs[i]).any(|(n, _)| n == nbrs[j]) {
                            continue;
                        }
                        looked += 1;
                        let u = (d.coords[nbrs[i] as usize] - c).normalized();
                        let v = (d.coords[nbrs[j] as usize] - c).normalized();
                        let deg = u.dot(v).clamp(-1.0, 1.0).acos().to_degrees();
                        assert!(
                            deg >= 89.0,
                            "{}:原子 {a} 处 {}–{a}–{} 的夹角只有 {deg:.1}°",
                            style.name,
                            nbrs[i],
                            nbrs[j]
                        );
                    }
                }
            }
            // **前提要自己成立**:度数过滤与三元环跳过若把所有对都排除掉,
            // 上面那圈断言一次都不跑,判据就空过了。
            assert!(looked > 0, "{}:一对键角都没查到,判据空过了", style.name);
        }
    }

    /// 两端都有标签的原子对,盒不许重叠 —— 与审计那一档同一个口径。
    fn labels_do_not_overlap(smi: &str) {
        let mut m = omgkit_io::smiles::parse(smi).expect("SMILES 该能解析");
        omgkit_chem::pipeline::sanitize(&mut m).expect("该能 sanitize");
        for style in &crate::style::Style::ALL {
            let d = crate::generate(&m, style);
            let grown = d.drawn(&m);
            // **前提要自己成立**:这个分子得真有两个带标签的原子挂在同一个
            // 原子上,否则这条判据查的不是它想查的东西。
            let mut geminal = 0usize;
            for a in 0..u32::try_from(grown.num_atoms()).expect("原子数超出 u32") {
                let labelled: Vec<u32> = grown
                    .neighbors(a)
                    .map(|(n, _)| n)
                    .filter(|n| crate::render::label_at(&grown, *n, style, &d.coords).is_some())
                    .collect();
                geminal += labelled.len() * labelled.len().saturating_sub(1) / 2;
                for i in 0..labelled.len() {
                    for j in (i + 1)..labelled.len() {
                        let (x, y) = (labelled[i], labelled[j]);
                        let la =
                            crate::render::label_at(&grown, x, style, &d.coords).expect("刚才筛过");
                        let lb =
                            crate::render::label_at(&grown, y, style, &d.coords).expect("刚才筛过");
                        let dv = d.coords[x as usize] - d.coords[y as usize];
                        assert!(
                            dv.x.abs() >= la.half_w + lb.half_w
                                || dv.y.abs() >= la.half_h + lb.half_h,
                            "{}:{smi} 的原子 {x} 与 {y} 挂在 {a} 上,两个标签盒叠着",
                            style.name
                        );
                    }
                }
            }
            assert!(
                geminal >= 2,
                "{smi} 只有 {geminal} 对同枢纽的带标签原子,判据空过"
            );
        }
    }

    /// 从一根键挂出去的环,不许摞在兄弟环身上。
    ///
    /// 分子取自真语料(第 6398 行):六配位钴挂四个吡啶 + 两个硫氰酸根。
    /// 配体按 60° 分开、配体氮到金属一个键长、环心到氮又一个键长,于是**两个
    /// 相邻环心相距正好 2、外接圆半径各 1 —— 精确相切**,切点又正落在 30° 栅格
    /// 上,两个环各有一个邻位碳**逐位重合**。
    ///
    /// 挑方向那一步先前只看"配体氮这一个原子撞不撞",看不见它身后那一整个环。
    /// 全量语料上这一族占 `原子不重合` 57 处里的一大半;补上整块前瞻之后 57 → 8。
    ///
    /// 变异:把 `layout_component` 里建 `blocks` 的那一段换成空表 —— 前瞻拿不到
    /// 块,`Lookahead::cost` 退回"这个原子自己就是一块"那一支。实测这条当场红:
    /// `原子 8 与 20 相距 0.0000`,而下面那条仍是绿的。
    ///
    /// (这段先前写的是那一支"恒返回 `((0,0), None)`"。那是**更早**的行为:
    /// 它现在返回单点块的真实代价,而改掉它正是同一轮做的事 —— 见
    /// [`super::Lookahead::cost`] 里那段注释。变异的结论不受影响,前提写错了。)
    #[test]
    fn a_ring_hanging_off_a_bond_is_not_dropped_onto_its_neighbour() {
        drawn_without_overlap(
            "N#CS[Co](SC#N)([N+]1=CC=CC=C1)([N+]2=CC=CC=C2)([N+]3=CC=CC=C3)[N+]4=CC=CC=C4",
            4,
        );
    }

    /// 「逐位重合的对数」必须排在「碰撞深度」前面。
    ///
    /// 分子取自真语料(第 3469 行):钴上挂四个甲基吡啶 + 两个硫氰酸根。
    ///
    /// 这条守的是 [`super::block_cost`] 返回值里的**第一位**。光比深度是不够的
    /// —— 一次**精确相切**只有一对原子重合,深度 `(0.5 − 0)² = 0.25`;而"整块
    /// 蹭到三四个原子上"虽然一处都没重合,深度和却更大。于是纯按深度会**主动
    /// 选中相切**,也就是主动选中"两个原子画在同一点上"。
    ///
    /// 变异:把 `block_cost` 的返回值从 `(same, q)` 改成 `(0, q)`(只按深度)。
    /// 实测这条当场红(`原子 12 与 20 相距 0.0000`),另外两条仍绿;全量语料上
    /// 那个变异让 `原子不重合` 8 → 10,坏的正是这个分子与第 3472 行的锌盐。
    #[test]
    fn how_many_atoms_land_on_top_of_each_other_outranks_how_deep_they_press() {
        drawn_without_overlap(
            "CC1=[N+](C=CC=C1)[Co](SC#N)(SC#N)([N+]2=C(C)C=CC=C2)([N+]3=C(C)C=CC=C3)\
             [N+]4=C(C)C=CC=C4",
            4,
        );
    }

    /// 一个候选都不干净时,要挑**最不坏**的,不许退回理想方向。
    ///
    /// 分子取自真语料(第 7794 行):三苯基甲基苯基酮,枢纽是个度 4 的季碳,
    /// 四个方向上挂着四个苯环。
    ///
    /// 这条守的是 [`super::free_direction`] 里那个**不能写成布尔**的判断。写成
    /// 布尔("整块撞了就否决")的话:前一个取代基为了躲开一处深度 0.0538 的
    /// 轻微重叠跳到别的档,**给后面的兄弟挖了坑**;轮到兄弟时一个候选都满足
    /// 不了,于是掉进兜底那一档退回 `ideal` —— 而 `ideal` 恰恰是唯一逐位重合的
    /// 那个。全量语料实测那一版**新坏 12 张**,签名全是"度 4 枢纽、60°、环心距
    /// 2.000"。
    ///
    /// 变异:把 `free_direction` 里的 `pick` 换成布尔式 —— 只收 `c == (0, 0)` 的
    /// 候选,一个都没有就让这一轮空手而归(于是掉进最后的 `ideal`)。实测
    /// 这条当场红:`原子 14 与 26 相距 0.0000`,而上面那条仍是绿的 ——
    /// 两条判据守的**不是**同一件事。
    #[test]
    fn when_nothing_is_clean_the_least_bad_direction_wins_instead_of_the_ideal_one() {
        drawn_without_overlap("O=C(C1=CC=CC=C1)C(C2=CC=CC=C2)(C3=CC=CC=C3)C4=CC=CC=C4", 3);
    }

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
        let mut checked = 0usize;
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
                            checked += 1;
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
        assert!(checked > 0, "一个键角都没查到,判据空过了");
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

    /// 老判据里没有任何环系统要摆,`blocks` 给空表 —— [`super::Lookahead::cost`]
    /// 于是走"这个原子自己就是一块"那一支,只按单点算代价。
    fn env<'a>(
        m: &'a omgkit_core::MolBuilder,
        ranks: &'a [u32],
        style: &'a Style,
        radii: &'a [f64],
        bonded: &'a std::collections::BTreeSet<(u32, u32)>,
        blocks: &'a BTreeMap<u32, super::Block>,
    ) -> super::Env<'a> {
        // 判据里的分子都是从理想栅格算出来的(没有桥环松弛),`off_grid` 给空表
        // —— 于是 `pinch_floor` 与从前一样只由度数决定。
        static NONE: std::sync::OnceLock<std::collections::BTreeSet<u32>> =
            std::sync::OnceLock::new();
        super::Env {
            mol: m,
            ranks,
            style,
            radii,
            bonded,
            off_grid: NONE.get_or_init(Default::default),
            blocks,
        }
    }

    #[test]
    fn a_chain_zigzags_instead_of_running_straight() {
        // 度数 2 的原子若取 180°,链就成了一条直线 —— 既不是化学惯例,
        // 也让后面的可旋转键无处可翻。这条守的正是那个 180°。
        let m = prep("CCCCC");
        let style = Style::ACS_1996;
        let radii = crate::refine::radii(&m, &style);
        let bonded: std::collections::BTreeSet<(u32, u32)> = m
            .bonds()
            .iter()
            .map(|b| (b.begin.min(b.end), b.begin.max(b.end)))
            .collect();
        let blocks: BTreeMap<u32, super::Block> = BTreeMap::new();
        let mut pos: BTreeMap<u32, Point2> = BTreeMap::new();
        pos.insert(0, Point2::ORIGIN);
        let mut zig = 1i8;
        for a in 0..4u32 {
            let out = place_neighbours(
                &env(&m, &canonical(&m), &style, &radii, &bonded, &blocks),
                a,
                &pos,
                &[a + 1],
                zig,
            );
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
        let radii = crate::refine::radii(&m, &style);
        let bonded: std::collections::BTreeSet<(u32, u32)> = m
            .bonds()
            .iter()
            .map(|b| (b.begin.min(b.end), b.begin.max(b.end)))
            .collect();
        let blocks: BTreeMap<u32, super::Block> = BTreeMap::new();
        let mut pos: BTreeMap<u32, Point2> = BTreeMap::new();
        pos.insert(1, Point2::ORIGIN);
        let mut todo: Vec<u32> = m.neighbors(1).map(|(n, _)| n).collect();
        todo.sort_unstable();
        for p in place_neighbours(
            &env(&m, &canonical(&m), &style, &radii, &bonded, &blocks),
            1,
            &pos,
            &todo,
            1,
        ) {
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
        let radii = crate::refine::radii(&m, &style);
        let bonded: std::collections::BTreeSet<(u32, u32)> = m
            .bonds()
            .iter()
            .map(|b| (b.begin.min(b.end), b.begin.max(b.end)))
            .collect();
        let blocks: BTreeMap<u32, super::Block> = BTreeMap::new();
        let mut pos: BTreeMap<u32, Point2> = BTreeMap::new();
        pos.insert(1, Point2::ORIGIN);
        let mut todo: Vec<u32> = m.neighbors(1).map(|(n, _)| n).collect();
        todo.sort_unstable();
        assert_eq!(todo.len(), 4, "季碳应当有四个邻居");
        let out = place_neighbours(
            &env(&m, &canonical(&m), &style, &radii, &bonded, &blocks),
            1,
            &pos,
            &todo,
            1,
        );
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
        let radii = crate::refine::radii(&m, &style);
        let bonded: std::collections::BTreeSet<(u32, u32)> = m
            .bonds()
            .iter()
            .map(|b| (b.begin.min(b.end), b.begin.max(b.end)))
            .collect();
        let blocks: BTreeMap<u32, super::Block> = BTreeMap::new();
        let mut pos: BTreeMap<u32, Point2> = BTreeMap::new();
        pos.insert(1, Point2::ORIGIN);
        // 手工把两个邻居摆在 0° 和 60°,留下一个 300° 的大空隙
        pos.insert(0, Point2::new(1.0, 0.0));
        pos.insert(2, Point2::new(0.5, 3f64.sqrt() / 2.0));
        let out = place_neighbours(
            &env(&m, &canonical(&m), &style, &radii, &bonded, &blocks),
            1,
            &pos,
            &[3],
            1,
        );
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
