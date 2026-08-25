//! 子结构匹配的回溯搜索。
//!
//! # 查询原子的处理顺序决定一切
//!
//! 回溯搜索的代价几乎完全由"先试哪个查询原子"决定。两条规则:
//!
//! 1. **每个新原子都要与已映射的原子相连**。这样每加一个原子就立刻能用键去
//!    筛候选,候选集从"整个分子"缩到"某个已映射原子的邻居"。做不到这点的话,
//!    每一层都要遍历全部目标原子,层数一多就是指数。
//! 2. **稀有的先试**。候选少的查询原子先钉住,能在浅层就把大部分分支砍掉。
//!    候选数由元素约束估出来(见 `allowed_elements`)—— 元素是最强也最便宜的
//!    那一维筛子。
//!
//! 第 1 条是正确性无关但性能致命的;第 2 条是纯启发式,但差别很大:
//! `CCCCCCCCBr` 配到 500 个碳 + 1 个溴的链上,只看度数要 1.96 ms
//! (链上度数几乎都一样,挑中的是碳,几百个碳各试一遍),看候选数是 31 µs。
//!
//! # 递归 SMARTS 的求值会重复
//!
//! `$(...)` 的求值本身就是一次完整的匹配。同一个 (子模式, 目标原子) 组合在
//! 一次搜索里会被问到很多次,所以要缓存 —— 否则一个 `[$(...)]` 就能把复杂度
//! 乘上一整轮搜索。

use std::collections::HashMap;

use omgkit_core::{BondDirection, BondOrder, BondStereo, MolBuilder};
use omgkit_io::smarts::{
    allowed_elements, atom_matches, bond_matches, required_chirality, BondExpr, BondPrim,
    BondProps, QueryMol,
};
use omgkit_io::stereo;

use crate::props::MolProps;

/// 扩展时给候选查询原子打的分,越大越先取:
/// (与已放置原子的连边数, 候选数的逆序, 度数)
type OrderKey = (usize, std::cmp::Reverse<usize>, usize);

/// 匹配的选项。
#[derive(Debug, Clone, Copy)]
pub struct MatchOptions {
    /// 最多返回多少个匹配。0 表示不限。
    ///
    /// 只问"有没有"时设成 1 —— 高度对称的分子上匹配数可以爆炸,
    /// 而判断"含不含某个子结构"根本不需要把它们全枚举出来。
    pub max_matches: usize,
    /// 是否按**目标原子集合**去重。
    ///
    /// 对称等价的映射会给出同一组原子的多个排列。关心"命中了哪些原子"时
    /// 应当开启;关心"有多少种映射方式"时不能开。
    pub uniquify: bool,
    /// 查询里写的手性与顺反算不算数。
    ///
    /// # 为什么摆成显式选项
    ///
    /// 这个开关的影响面很大:2000 个分子 × 776 条模式里,开与不开结果有差别的
    /// 组合占 **23%**。这种规模的分歧不能由一个默认值悄悄定下,所以做成必须
    /// 填的字段,逼每个调用点自己表态。
    ///
    /// 默认 **`true`(判)**:作者写了 `[C@]` 或 `/C=C/` 却被忽略时不会报错,
    /// 只会悄悄多出一批匹配 —— 那种错要等到下游才暴露,而且很难回溯。
    ///
    /// 例外是 [`run_reactants`](crate::run_reactants),它显式关掉:反应模板是
    /// 跨工具流通的东西,读得更严会让现成的模板不再出产物,而"少了产物"
    /// 比"多了产物"更难发现。
    pub use_chirality: bool,
}

impl Default for MatchOptions {
    fn default() -> Self {
        Self {
            max_matches: 0,
            uniquify: true,
            use_chirality: true,
        }
    }
}

/// 一次匹配:`mapping[i]` 是查询原子 `i` 对应的目标原子。
pub type Mapping = Vec<u32>;

/// 找出查询在分子中的全部匹配。
///
/// `props` 必须由**同一个** `mol` 算出来 —— 两者不一致时匹配结果没有意义,
/// 而且不会报错。
/// 一次搜索花了多少工夫。
///
/// # 为什么要把这个数交出来
///
/// 回溯搜索天然是指数的,全靠剪枝压住;剪枝失效时**结果照样全对,只是慢**。
/// 守这件事只能量"工夫",而先前量的是**墙钟**——
/// `crates/omgkit-match/tests/scaling.rs` 的两条增长曲线判据在 8/16/32 元模式上
/// 测出 213 µs / 440 µs / 1374 µs,微秒级的数放在共享 CI 机器上,
/// 2026-08-25 直接把一个**改的是别的 crate**的提交打红了(涨幅 1.61 > 阈值 1.6)。
///
/// `candidate_tests` 是**整数、确定、与机器无关**的:同一份输入永远同一个数。
/// 墙钟仍然留着,但降成一道很粗的闸(常数因子崩了才拦),细的那道交给这个数。
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SearchStats {
    /// 回溯时**考虑过**的候选原子次数(在 `used` / 可行性筛之前就计)。
    ///
    /// 剪枝失效的典型样子是候选集本身变大 —— 定序断开时某一层会退化成
    /// 全分子扫描,这个数立刻按分子大小翻倍。
    pub candidate_tests: u64,
}

/// 找出全部匹配,**并交出这次搜索花了多少工夫**。
///
/// 语义与 [`substructure_matches`] 完全一致,只是多返回一个 [`SearchStats`]。
#[must_use]
pub fn substructure_matches_counted(
    query: &QueryMol,
    mol: &MolBuilder,
    props: &MolProps,
    opts: MatchOptions,
) -> (Vec<Mapping>, SearchStats) {
    if query.num_atoms() == 0 || query.num_atoms() > mol.num_atoms() {
        return (Vec::new(), SearchStats::default());
    }
    let mut ctx = Ctx {
        mol,
        props,
        candidate_tests: 0,
        recursive_cache: HashMap::new(),
    };
    let counts = candidate_counts(query, props);
    let order = search_order(query, &counts);
    let mut out = Vec::new();
    let mut seen: std::collections::HashSet<Vec<u32>> = std::collections::HashSet::new();

    let mut mapping = vec![u32::MAX; query.num_atoms()];
    let mut used = vec![false; mol.num_atoms()];
    extend(
        query,
        &order,
        0,
        &mut mapping,
        &mut used,
        &mut ctx,
        &opts,
        &mut seen,
        &mut out,
    );
    (
        out,
        SearchStats {
            candidate_tests: ctx.candidate_tests,
        },
    )
}

/// 找出查询在分子中的全部匹配。
///
/// `props` 必须由**同一个** `mol` 算出来 —— 两者不一致时匹配结果没有意义,
/// 而且不会报错。
#[must_use]
pub fn substructure_matches(
    query: &QueryMol,
    mol: &MolBuilder,
    props: &MolProps,
    opts: MatchOptions,
) -> Vec<Mapping> {
    substructure_matches_counted(query, mol, props, opts).0
}

/// 查询在分子中是否至少有一个匹配,且查询原子 0 落在 `root` 上。
///
/// 递归 SMARTS `$(...)` 的语义就是这个:该原子要能作为子模式的**首原子**
/// 匹配上。
fn matches_rooted(query: &QueryMol, root: u32, ctx: &mut Ctx) -> bool {
    if query.num_atoms() == 0 || query.num_atoms() > ctx.mol.num_atoms() {
        return false;
    }
    let counts = candidate_counts(query, ctx.props);
    let order = search_order(query, &counts);
    // 把查询原子 0 排到最前,它才能被钉在 root 上
    let order = if order.first() == Some(&0) {
        order
    } else {
        let mut o = vec![0u32];
        o.extend(order.into_iter().filter(|&x| x != 0));
        o
    };

    let mut mapping = vec![u32::MAX; query.num_atoms()];
    let mut used = vec![false; ctx.mol.num_atoms()];
    let opts = MatchOptions {
        max_matches: 1,
        uniquify: false,
        // 递归 `$(...)` 里的立体照判 —— 子模式与外层是同一套语义,
        // 外层判而里层不判会让 `[$([C@](...))]` 静默失效
        use_chirality: true,
    };
    let mut seen = std::collections::HashSet::new();
    let mut out = Vec::new();

    // 手工放置第一个原子,其余交给通用的扩展
    if !atom_feasible(query, order[0], root, &mapping, ctx) {
        return false;
    }
    mapping[order[0] as usize] = root;
    used[root as usize] = true;
    extend(
        query,
        &order,
        1,
        &mut mapping,
        &mut used,
        ctx,
        &opts,
        &mut seen,
        &mut out,
    );
    !out.is_empty()
}

struct Ctx<'a> {
    mol: &'a MolBuilder,
    props: &'a MolProps,
    /// 搜索工作量,见 [`SearchStats`]
    candidate_tests: u64,
    /// (子模式指针, 目标原子) → 是否匹配。
    ///
    /// 用指针当键是安全的:子模式活在查询树里,整个搜索期间不会移动。
    recursive_cache: HashMap<(usize, u32), bool>,
}

/// 每个查询原子的**候选数上界**:目标里有多少个原子的元素可能对得上。
///
/// 元素是最强也最便宜的一维筛子。推不出元素约束时返回原子总数(即"不筛")——
/// 估计只能偏大,偏小会让起点挑错,但不会导致漏解:候选数只用于排序。
fn candidate_counts(query: &QueryMol, props: &MolProps) -> Vec<usize> {
    let mut by_element = [0usize; 256];
    for a in &props.atoms {
        by_element[a.atomic_num as usize] += 1;
    }
    let total = props.atoms.len();
    (0..query.num_atoms())
        .map(|i| match allowed_elements(&query.atoms[i]) {
            Some(set) => set.iter().map(|&z| by_element[z as usize]).sum(),
            None => total,
        })
        .collect()
}

/// 决定查询原子的处理顺序。见模块文档。
fn search_order(query: &QueryMol, counts: &[usize]) -> Vec<u32> {
    let n = query.num_atoms();
    let topo = &query.topology;
    let mut order = Vec::with_capacity(n);
    let mut placed = vec![false; n];

    // 起点:**候选最少**的原子,并列时取度数大的。
    //
    // 按度数挑是不够的:`CCCCCCCCBr` 这样的链上度数几乎都一样,挑中的会是
    // 某个碳,于是几百个碳都要各试一遍。按候选数挑会直接钉在溴上 —— 只有
    // 一个候选,整条链一次就走完。
    while order.len() < n {
        let start = (0..n as u32)
            .filter(|&a| !placed[a as usize])
            .min_by_key(|&a| (counts[a as usize], std::cmp::Reverse(topo.degree(a))))
            .expect("还有未放置的原子");
        order.push(start);
        placed[start as usize] = true;

        // 之后每次都挑"与已放置原子相连最多"的那个 —— 连接越多,
        // 加进来时能立刻验证的键就越多,剪枝越早
        loop {
            let mut best: Option<(u32, OrderKey)> = None;
            for a in 0..n as u32 {
                if placed[a as usize] {
                    continue;
                }
                let links = topo
                    .neighbors(a)
                    .filter(|&(o, _)| placed[o as usize])
                    .count();
                if links == 0 {
                    continue;
                }
                // 已连边多的优先(能立刻验的键多),其次候选少的
                let key: OrderKey = (links, std::cmp::Reverse(counts[a as usize]), topo.degree(a));
                if best.map_or(true, |(_, k)| key > k) {
                    best = Some((a, key));
                }
            }
            match best {
                Some((a, _)) => {
                    order.push(a);
                    placed[a as usize] = true;
                }
                // 该连通片段走完了,回到外层挑下一个片段的起点
                None => break,
            }
        }
    }
    order
}

/// 目标原子 `t` 能否承载查询原子 `q`(只看原子本身的条件)。
fn atom_feasible(query: &QueryMol, q: u32, t: u32, _mapping: &[u32], ctx: &mut Ctx) -> bool {
    let props = ctx.props.atoms[t as usize];
    // 先分离出递归子模式的求值,避免在闭包里再借 ctx
    let mut resolve = |sub: &QueryMol| {
        let key = (sub as *const QueryMol as usize, t);
        if let Some(&v) = ctx.recursive_cache.get(&key) {
            return v;
        }
        // 先放一个占位,防止病态模式自引用时无限递归
        ctx.recursive_cache.insert(key, false);
        let v = matches_rooted(sub, t, ctx);
        ctx.recursive_cache.insert(key, v);
        v
    };
    atom_matches(&query.atoms[q as usize], &props, &mut resolve)
}

/// 查询里 `q` 与已映射原子之间的每条键,在目标里都要存在且满足键查询。
fn bonds_feasible(query: &QueryMol, q: u32, t: u32, mapping: &[u32], ctx: &Ctx) -> bool {
    for (other_q, qbond) in query.topology.neighbors(q) {
        let mapped = mapping[other_q as usize];
        if mapped == u32::MAX {
            continue; // 对端还没映射,这条键留到那时再验
        }
        let Some(tbond) = ctx.mol.bond_between(t, mapped) else {
            return false; // 查询要求相连,目标里却没有这条键
        };
        if !bond_ok(query, qbond, tbond, mapping, ctx) {
            return false;
        }
    }
    true
}

fn bond_ok(query: &QueryMol, qbond: u32, tbond: u32, mapping: &[u32], ctx: &Ctx) -> bool {
    let qb = query.topology.bonds()[qbond as usize];
    let tb = ctx.mol.bonds()[tbond as usize];
    let mut props: BondProps = ctx.props.bonds[tbond as usize];

    // 配位键的方向有语义。查询表达式是相对**查询里书写的朝向**解析的,
    // 所以要问:目标键的给体,是不是对应查询键的 begin 端?
    if tb.order == BondOrder::Dative {
        let want_donor = mapping[qb.begin as usize];
        props.dative_forward = want_donor != u32::MAX && tb.begin == want_donor;
    }
    bond_matches(&query.bonds[qbond as usize], &props)
}

/// 查询侧的顺反,从**表达式树**里读。
///
/// 不能复用 [`stereo::raw_cis_trans`]:查询的拓扑里键级是占位的,方向也存在
/// `BondExpr` 里而不是 `direction` 字段。真分子那一侧才走那个函数。
fn query_cis_trans(query: &QueryMol, bond: u32) -> Option<(BondStereo, [u32; 2])> {
    if !expr_has(&query.bonds[bond as usize], BondPrim::Double) {
        return None;
    }
    let b = query.topology.bonds()[bond as usize];
    let (ra, da) = query_outward(query, b.begin, b.end)?;
    let (rb, dbi) = query_outward(query, b.end, b.begin)?;
    Some((
        if da == dbi {
            BondStereo::Cis
        } else {
            BondStereo::Trans
        },
        [ra, rb],
    ))
}

/// `end` 那一侧写了方向的邻居,方向换算到"从 `end` 向外"。
fn query_outward(query: &QueryMol, end: u32, other: u32) -> Option<(u32, BondDirection)> {
    query
        .topology
        .neighbors(end)
        .filter(|&(o, _)| o != other)
        .find_map(|(o, bi)| {
            let e = &query.bonds[bi as usize];
            let raw = if expr_has(e, BondPrim::UpRight) {
                BondDirection::UpRight
            } else if expr_has(e, BondPrim::DownRight) {
                BondDirection::DownRight
            } else {
                return None;
            };
            let tb = query.topology.bonds()[bi as usize];
            Some((o, if tb.begin == end { raw } else { raw.flipped() }))
        })
}

/// 表达式里(在合取意义下)是否出现了某个基元。
///
/// `,` 与 `!` 下面的方向要求彼此矛盾,取任一个都不对 —— 这类写法罕见,
/// 放行比猜一个好。
fn expr_has(e: &BondExpr, want: BondPrim) -> bool {
    match e {
        BondExpr::Prim(p) => *p == want,
        BondExpr::And(parts) => parts.iter().any(|x| expr_has(x, want)),
        BondExpr::Or(_) | BondExpr::Not(_) => false,
    }
}

/// 校验映射是否满足查询里写的双键顺反要求。
///
/// # 为什么不能在 `bond_matches` 里判
///
/// 与手性同源:`/` 相对键自己的 `begin → end` 朝向,查询与底物的朝向不同,
/// 直接比就成了"比写法" —— `F/C=C/F` 与 `C(\F)=C/F` 是同一个分子,直接比
/// 会一个匹配、一个不匹配。
///
/// 而且单独一条 `/` 什么都不表示:顺反是**双键加两条参照键**的性质。
///
/// # 参照原子不同时要翻
///
/// 两边挑中的参照原子未必对应。查询的参照映过去若不是底物记的那个,说明落在
/// 双键的另一侧,顺反要翻一次;两端都不同就翻两次,等于没翻。
fn cis_trans_ok(query: &QueryMol, mapping: &Mapping, mol: &MolBuilder) -> bool {
    for qb in 0..query.topology.num_bonds() as u32 {
        let Some((want, qrefs)) = query_cis_trans(query, qb) else {
            continue;
        };
        let b = query.topology.bonds()[qb as usize];
        let (tb0, tb1) = (mapping[b.begin as usize], mapping[b.end as usize]);
        let Some((_, ti)) = mol.neighbors(tb0).find(|&(o, _)| o == tb1) else {
            return false;
        };
        let Some((got, trefs)) = stereo::raw_cis_trans(mol, ti) else {
            // 底物这根双键没有成对的方向 —— 未指定,配不上写死构型的查询
            return false;
        };
        // 底物那根键的 begin 未必对应查询的 begin,参照对要跟着摆正
        let tb = mol.bonds()[ti as usize];
        let trefs = if tb.begin == tb0 {
            trefs
        } else {
            [trefs[1], trefs[0]]
        };
        let mut flips = 0;
        for (i, &qr) in qrefs.iter().enumerate() {
            if mapping[qr as usize] != trefs[i] {
                flips += 1;
            }
        }
        let effective = if flips % 2 == 1 { flipped(got) } else { got };
        if effective != want {
            return false;
        }
    }
    true
}

fn flipped(s: BondStereo) -> BondStereo {
    match s {
        BondStereo::Cis => BondStereo::Trans,
        BondStereo::Trans => BondStereo::Cis,
        other => other,
    }
}

/// 校验映射是否满足查询里写的手性要求。
///
/// # 为什么不能在 `atom_matches` 里判
///
/// 手性标记相对**各自分子的邻居存储顺序**。查询与底物的存储顺序不同,直接比
/// 原始标记就是拿两个参照系里的值去比,结论可以正好相反。要算宇称,必须知道
/// 查询的每个邻居映到了底物的哪个原子,而那要等映射齐了。
///
/// # 查询原子度 < 3 时不判构型
///
/// 四面体要三个显式邻居(加原子自身)才定得下构型。度更小的查询是**欠定**的,
/// 此时只要求底物"有手性",不判是哪一个。
fn chirality_ok(query: &QueryMol, mapping: &Mapping, mol: &MolBuilder) -> bool {
    for (qi, expr) in query.atoms.iter().enumerate() {
        let Some(want) = required_chirality(expr) else {
            continue;
        };
        if !want.is_tetrahedral() {
            continue;
        }
        let t = mapping[qi];
        let got = mol.atoms()[t as usize].chiral_tag;
        if !got.is_tetrahedral() {
            return false;
        }
        let qn: Vec<u32> = query
            .topology
            .neighbors(qi as u32)
            .map(|(o, _)| o)
            .collect();
        // 欠定:只要求"有手性"。
        //
        // 不能把"一个邻居都没写"(`[C@H]`)特判成永不匹配 —— 这类模式落在**真**
        // 手性中心上本就该命中。**非真**中心上的标记该不该留是解析与净化的事
        // (本库选择保留,见 [`omgkit_io::stereo`]),不由匹配逻辑代为决定。
        if qn.len() < 3 {
            continue;
        }
        // 查询邻居映到底物,再与底物的存储序比宇称。
        //
        // 底物可能比查询多一个邻居(度 3 的查询配度 4 的底物)。**那一个的位置
        // 会改变宇称**,只拿映射到的三个去比会得到相反的构型。所以补成完整的
        // 置换:未映射的那个按它在底物里的存储位置留在原处,查询侧则把它放在
        // 末尾(查询没写它,只能排在写出来的三个之后)。
        let imaged: Vec<u32> = qn.iter().map(|&o| mapping[o as usize]).collect();
        let stored: Vec<u32> = mol.neighbors(t).map(|(o, _)| o).collect();
        let extra: Vec<u32> = stored
            .iter()
            .copied()
            .filter(|o| !imaged.contains(o))
            .collect();
        if stored.len() != imaged.len() + extra.len() {
            continue;
        }
        let mut query_side = imaged.clone();
        query_side.extend(extra);
        let Some(odd) = permutation_is_odd(&query_side, &stored) else {
            continue;
        };
        let effective = if odd { got.inverted() } else { got };
        if effective != want {
            return false;
        }
    }
    true
}

/// `from` → `to` 置换的宇称。两者不是同一多重集时返回 `None`。
fn permutation_is_odd(from: &[u32], to: &[u32]) -> Option<bool> {
    if from.len() != to.len() {
        return None;
    }
    let mut cur = from.to_vec();
    let mut swaps = 0usize;
    for i in 0..to.len() {
        if cur[i] == to[i] {
            continue;
        }
        let j = (i + 1..cur.len()).find(|&j| cur[j] == to[i])?;
        cur.swap(i, j);
        swaps += 1;
    }
    Some(swaps % 2 == 1)
}

#[allow(clippy::too_many_arguments)]
fn extend(
    query: &QueryMol,
    order: &[u32],
    depth: usize,
    mapping: &mut Mapping,
    used: &mut [bool],
    ctx: &mut Ctx,
    opts: &MatchOptions,
    seen: &mut std::collections::HashSet<Vec<u32>>,
    out: &mut Vec<Mapping>,
) {
    if opts.max_matches != 0 && out.len() >= opts.max_matches {
        return;
    }
    if depth == order.len() {
        // 手性与顺反都要等映射齐了才判 —— 见各自的函数文档
        if opts.use_chirality
            && (!chirality_ok(query, mapping, ctx.mol) || !cis_trans_ok(query, mapping, ctx.mol))
        {
            return;
        }
        if opts.uniquify {
            let mut key = mapping.clone();
            key.sort_unstable();
            if !seen.insert(key) {
                return;
            }
        }
        out.push(mapping.clone());
        return;
    }

    let q = order[depth];
    // 候选来自"某个已映射邻居的邻居";没有已映射邻居时才退化成全分子扫描
    let anchor = query
        .topology
        .neighbors(q)
        .map(|(o, _)| o)
        .find(|&o| mapping[o as usize] != u32::MAX);

    let candidates: Vec<u32> = match anchor {
        Some(o) => ctx
            .mol
            .neighbors(mapping[o as usize])
            .map(|(other, _)| other)
            .collect(),
        None => (0..ctx.mol.num_atoms() as u32).collect(),
    };

    for t in candidates {
        // 工作量的计量单位,见 [`SearchStats`]。放在 `used` 判断**之前** ——
        // 剪枝失效的典型样子正是"候选集本身变大",那些候选多半立刻被 `used`
        // 或 `atom_feasible` 挡掉,放在后面就数不到。
        ctx.candidate_tests += 1;
        if used[t as usize] {
            continue;
        }
        if !atom_feasible(query, q, t, mapping, ctx) {
            continue;
        }
        mapping[q as usize] = t;
        if bonds_feasible(query, q, t, mapping, ctx) {
            used[t as usize] = true;
            extend(query, order, depth + 1, mapping, used, ctx, opts, seen, out);
            used[t as usize] = false;
        }
        mapping[q as usize] = u32::MAX;
        if opts.max_matches != 0 && out.len() >= opts.max_matches {
            return;
        }
    }
}
