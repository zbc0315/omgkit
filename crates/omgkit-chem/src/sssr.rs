//! 环集感知 —— 给出分子的一组具体的环,供芳香性感知使用。
//!
//! # 这个集合到底是什么
//!
//! 算出来的是:**一组最小环基,外加所有"不能表示为更短环之和"的 Horton 候选环**。
//!
//! 判据(按环长升序处理,同长度的一起判):
//!
//! > 环 C 入选 ⟺ C 不在"严格更短的环"张成的 GF(2) 空间中。
//!
//! 这个集合夹在两个熟知的对象之间:
//!
//! ```text
//! 最小环基  ⊆  本模块的环集  ⊆  相关环(relevant cycles,即全部最小环基之并)
//! ```
//!
//! 三者在绝大多数分子上重合,但在**高对称的笼状/带状/螺环**体系上会分开:
//!
//! - 相关环可以有**指数多个**(n 个螺环丁烷串成的环有 2ⁿ 个六元环),
//!   所以"返回全部相关环"这个 API 形状在一般图上根本不可行
//! - 本模块只走**单棵 BFS 最短路树**的 Horton 候选集,因而会漏掉一部分相关环
//!
//! **不要把这个集合当作相关环使用。** 需要完备的相关环时要走 Vismara 的
//! 原型/环族表示,那是另一件事。
//!
//! # 与差分基准的已知分歧
//!
//! 基准走的是"经典 SSSR + 同尺寸可替换环补回",那条补回路径在高对称体系上
//! 并不完全。本模块在笼状体系上给出的环**更多**:
//!
//! | 分子 | 基准 | 本模块 |
//! |---|---|---|
//! | `C12C3C4C1C3C24` | 7 | 9 |
//! | `C1C23CC14CC(C2)(C3)C4` | 3 | 9 |
//! | `C1C2CC3CC4CC1C1CC2CC3CC4C1` | 4 | 6 |
//!
//! **这些分歧不会传导到下游**:上表每一条的芳香原子数与芳香键数两边完全一致
//! (笼状体系本就不芳香;`c12c3c4c1c3c24` 这条可芳香的也一致)。差分测试把它们
//! 登记在案并持续盯着,见 `differential_l2_ringset`。
//!
//! # 集合是规范的
//!
//! 与候选的枚举顺序无关:每一级判完再整组并入基,而"并入后的张成"就是全部长度
//! ≤ 当前值的环张成,与并入次序无关。实测 7553 条含环分子、每条随机重排原子
//! 编号 3 次,环集无一改变。
//!
//! # 算法与三处复杂度控制
//!
//! **Horton 候选集 + 更短者张成过滤**:对每个顶点做一次 BFS,对每条边取
//! `SP(v→x) + (x,y) + SP(y→v)`,要求两条最短路除 v 外不相交。
//!
//! 三处控制缺一不可,少任何一处都会在某种形状上退化:
//!
//! | 控制 | 针对的形状 | 不做的话 |
//! |---|---|---|
//! | 逐双连通分量 | 多环系分子 | O(原子数 × 键数) |
//! | 按环长迭代加深 | 单个大稠合体系 | 长链并苯 O(n³) |
//! | 平衡 + 首步剪枝 | 大环 | 800 元大环 1.1 秒 |
//!
//! **迭代加深**:一旦基的秩达到该分量的圈秩,更长的环必然落在更短者的张成中,
//! 不可能入选,于是停止加深。BFS 同时限深到 `⌈上限/2⌉`。
//!
//! **平衡剪枝**:只取两条最短路长度差 ≤ 1 的 (顶点, 边) 对。Horton 定理保证
//! 每个最小环基成员都能由某个顶点"对面"的边平衡地生成。
//!
//! ——但**这一条现在跑不到**,而且是结构上跑不到,见
//! `horton_candidates` 里那处 `debug_assert`(私有项,搜文件名即可)。
//!
//! **首步剪枝**:两条最短路若从 v 迈出的第一步相同,必然相交。这个判断是 O(1),
//! 而回溯出路径再判是 O(路径长)。

use omgkit_core::{BondOrder, MolBuilder};

/// 一个环。原子与键都按**环上顺序**排列。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Ring {
    /// 环上的原子,首尾相接
    pub atoms: Vec<u32>,
    /// 环上的键;`bonds[i]` 连接 `atoms[i]` 与 `atoms[(i+1) % len]`
    pub bonds: Vec<u32>,
}

impl Ring {
    /// 环的大小(原子数 = 键数)
    #[must_use]
    pub fn len(&self) -> usize {
        self.atoms.len()
    }

    /// 恒为假 —— 环至少有 3 个原子。提供它只是为了满足 clippy 的成对约定。
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }
}

/// 求分子的环集 —— 具体含义见[模块文档](self)。
///
/// 配位键不参与成环,与 [`crate::rings`] 一致。
///
/// 返回顺序:按(环大小, 最小原子号)升序,保证确定性。
#[must_use]
pub fn ring_set(mol: &MolBuilder) -> Vec<Ring> {
    ring_set_counted(mol).0
}

/// 环搜索**做了多少事**。整数、确定,debug 与 release 逐位相同。
///
/// # 为什么要数,而不是量耗时
///
/// 复杂度是"做了多少事"的性质,墙钟只是它乘上一个会抖的常数。而"每原子耗时
/// 涨幅"这个形状还漏掉一整类退化:按比例整体变慢时涨幅纹丝不动。
/// (同一条教训在 `omgkit-match` 的 `scaling.rs` 上栽过一次,那里已经换成
/// 数工作量并钉死绝对值。)
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct SearchStats {
    /// BFS 出队次数 —— 全部分量、全部加深轮次、全部起点合计
    pub bfs_visits: u64,
    /// 考察过的(起点, 边)对数。计在 `edge_done` 判重**之后**、
    /// 各条剪枝**之前** —— 计在剪枝之后的话,被剪掉的那些数不到,
    /// 而剪枝失效正是要守的东西。
    pub edge_tests: u64,
    /// 回溯两条最短路时走过的顶点数合计。
    ///
    /// **这一项专门盯剪枝。** `bfs_visits` 与 `edge_tests` 都记在剪枝**之前**,
    /// 所以平衡剪枝(两条路长度差 ≤ 1)与首步剪枝(第一步相同必相交)失效时
    /// 那两个数纹丝不动 —— 而模块文档说得很清楚,少了它们复杂度会从 O(n²)
    /// 退化到 O(n³)。回溯是它们保护的那个 O(路径长) 操作,数它才碰得到。
    pub path_steps: u64,
}

/// 与 [`ring_set`] 同一条路,外加工作量计数。
#[must_use]
pub fn ring_set_counted(mol: &MolBuilder) -> (Vec<Ring>, SearchStats) {
    let active: Vec<bool> = mol
        .bonds()
        .iter()
        .map(|b| b.order != BondOrder::Dative)
        .collect();
    let adj = crate::rings::Adjacency::build(mol, &active);

    let mut stats = SearchStats::default();
    let mut out: Vec<Ring> = Vec::new();
    for comp_bonds in crate::rings::biconnected_bond_components(&adj) {
        out.extend(component_ring_set(mol, &comp_bonds, &mut stats));
    }

    out.sort_by_key(|r| {
        (
            r.atoms.len(),
            r.atoms.iter().copied().min().unwrap_or(u32::MAX),
            r.atoms.iter().copied().max().unwrap_or(u32::MAX),
        )
    });
    (out, stats)
}

// ---------------------------------------------------------------------------

/// 分量内的局部图:原子与键都重编号到 0..n / 0..m,让位集短、缓存友好。
struct Component {
    /// 局部原子号 → 全局原子号
    atoms: Vec<u32>,
    /// 局部键号 → 全局键号
    bonds: Vec<u32>,
    /// CSR:每个局部原子的 (局部邻居, 局部键号)
    offset: Vec<u32>,
    nbr: Vec<(u32, u32)>,
    /// 每条局部边的两个端点。**必须预存** —— 现查是 O(原子数 × 度数),
    /// 而它在最内层循环里被调用,会让候选生成整体变成平方。
    edge_ends: Vec<(u32, u32)>,
}

impl Component {
    fn build(mol: &MolBuilder, comp_bonds: &[u32]) -> Self {
        let mut atoms: Vec<u32> = comp_bonds
            .iter()
            .flat_map(|&bi| {
                let b = mol.bonds()[bi as usize];
                [b.begin, b.end]
            })
            .collect();
        atoms.sort_unstable();
        atoms.dedup();

        // 全局原子号 → 局部号。分量通常很小,线性查找反而比建哈希快,
        // 但分量也可能有几十个原子,故用一个按全局号排序的表做二分。
        let local_of = |g: u32| atoms.binary_search(&g).expect("端点必在分量内") as u32;

        let n = atoms.len();
        let mut degree = vec![0u32; n];
        for &bi in comp_bonds {
            let b = mol.bonds()[bi as usize];
            degree[local_of(b.begin) as usize] += 1;
            degree[local_of(b.end) as usize] += 1;
        }
        let mut offset = vec![0u32; n + 1];
        for i in 0..n {
            offset[i + 1] = offset[i] + degree[i];
        }
        let mut cursor = offset[..n].to_vec();
        let mut nbr = vec![(0u32, 0u32); offset[n] as usize];
        let mut edge_ends = Vec::with_capacity(comp_bonds.len());
        for (lb, &bi) in comp_bonds.iter().enumerate() {
            let b = mol.bonds()[bi as usize];
            let (u, v) = (local_of(b.begin), local_of(b.end));
            edge_ends.push((u, v));
            let lb = lb as u32;
            nbr[cursor[u as usize] as usize] = (v, lb);
            cursor[u as usize] += 1;
            nbr[cursor[v as usize] as usize] = (u, lb);
            cursor[v as usize] += 1;
        }

        Self {
            atoms,
            bonds: comp_bonds.to_vec(),
            offset,
            nbr,
            edge_ends,
        }
    }

    fn n_atoms(&self) -> usize {
        self.atoms.len()
    }

    fn n_bonds(&self) -> usize {
        self.bonds.len()
    }

    fn neighbors(&self, a: u32) -> &[(u32, u32)] {
        let s = self.offset[a as usize] as usize;
        let e = self.offset[a as usize + 1] as usize;
        &self.nbr[s..e]
    }
}

/// 定长位集,按 64 位分块。
#[derive(Clone, PartialEq, Eq)]
struct BitSet(Vec<u64>);

impl BitSet {
    fn new(n_bits: usize) -> Self {
        Self(vec![0; n_bits.div_ceil(64)])
    }

    fn set(&mut self, i: usize) {
        self.0[i / 64] |= 1u64 << (i % 64);
    }

    fn xor_with(&mut self, other: &Self) {
        for (a, b) in self.0.iter_mut().zip(&other.0) {
            *a ^= b;
        }
    }

    fn is_zero(&self) -> bool {
        self.0.iter().all(|&w| w == 0)
    }

    /// 最高位的下标;全零时返回 `None`
    fn leading(&self) -> Option<usize> {
        self.0
            .iter()
            .enumerate()
            .rev()
            .find(|(_, &w)| w != 0)
            .map(|(i, &w)| i * 64 + (63 - w.leading_zeros() as usize))
    }
}

/// 增量维护的 GF(2) 基,按主元位下标索引。
struct Gf2Basis {
    /// 主元位 → 该主元对应的向量
    rows: Vec<Option<BitSet>>,
}

impl Gf2Basis {
    fn new(n_bits: usize) -> Self {
        Self {
            rows: vec![None; n_bits],
        }
    }

    /// 把 `v` 对当前基做约简,返回余量(全零表示 `v` 在张成中)。
    fn reduce(&self, v: &BitSet) -> BitSet {
        let mut cur = v.clone();
        while let Some(lead) = cur.leading() {
            match &self.rows[lead] {
                Some(row) => cur.xor_with(row),
                None => break,
            }
        }
        cur
    }

    /// 把 `v` 并入基;返回是否真的让秩增加了(即 `v` 原先不在张成中)。
    fn insert(&mut self, v: &BitSet) -> bool {
        let reduced = self.reduce(v);
        match reduced.leading() {
            Some(lead) => {
                self.rows[lead] = Some(reduced);
                true
            }
            None => false,
        }
    }
}

/// 一个候选环
struct Candidate {
    /// 按环上顺序的局部原子号
    atoms: Vec<u32>,
    /// 按环上顺序的局部键号
    bonds: Vec<u32>,
    /// 键集合的位表示,用于 GF(2) 运算
    bits: BitSet,
    /// 去重与排序用的规范键:排序后的局部原子号
    key: Vec<u32>,
}

fn component_ring_set(mol: &MolBuilder, comp_bonds: &[u32], stats: &mut SearchStats) -> Vec<Ring> {
    let c = Component::build(mol, comp_bonds);

    // 快路径:分量内所有顶点度数都是 2 ⇒ 这个分量本身就是一条简单环。
    //
    // 一般路径要对**每个**顶点做一次覆盖整个分量的 BFS,在大环上就是 O(n²):
    // 400 元大环要 8 毫秒。而大环恰恰是最常见的这种形状 —— 大环肽、大环内酯、
    // 冠醚的取代基都是桥边,双连通分解之后剩下的正是一条纯环。
    if let Some(ring) = single_cycle_component(&c) {
        return vec![ring];
    }

    // 分量是连通的,圈秩 = 边数 - 点数 + 1
    let rank_target = c.n_bonds() + 1 - c.n_atoms();

    // 迭代加深 —— 见模块文档"只找到必要的长度为止"
    let mut limit = INITIAL_LENGTH_LIMIT;
    loop {
        let capped = limit >= c.n_atoms();
        let (rings, rank) = cycles_up_to(&c, limit.min(c.n_atoms()), stats);
        if rank >= rank_target || capped {
            debug_assert!(
                rank >= rank_target,
                "搜遍全部长度仍未张满圈空间:秩 {rank} < 圈秩 {rank_target}"
            );
            return rings;
        }
        limit *= 2;
    }
}

/// 起步的环长上限。绝大多数分子的环都在这个范围内,一轮就收敛。
const INITIAL_LENGTH_LIMIT: usize = 8;

/// 分量若是一条简单环(所有顶点度数为 2),直接把它作为唯一的环返回。
///
/// 双连通分量的圈秩 = 边数 − 点数 + 1;全部度数为 2 时边数 = 点数,圈秩为 1,
/// 而这条环就是分量本身,不需要任何搜索。
fn single_cycle_component(c: &Component) -> Option<Ring> {
    let n = c.n_atoms();
    if n < 3 || c.n_bonds() != n {
        return None;
    }
    if (0..n as u32).any(|a| c.neighbors(a).len() != 2) {
        return None;
    }

    // 沿环走一圈,顺带确认它确实是一条环而不是多个分离的圈
    let mut atoms = Vec::with_capacity(n);
    let mut bonds = Vec::with_capacity(n);
    let (mut prev, mut cur) = (u32::MAX, 0u32);
    for _ in 0..n {
        atoms.push(c.atoms[cur as usize]);
        let (next, lb) = *c
            .neighbors(cur)
            .iter()
            .find(|&&(y, _)| y != prev)
            .expect("度数为 2 的顶点必有一个非来路的邻居");
        bonds.push(c.bonds[lb as usize]);
        prev = cur;
        cur = next;
    }
    // 走 n 步必须回到起点,否则分量不是单一环
    if cur != 0 {
        return None;
    }
    Some(Ring { atoms, bonds })
}

/// 求"长度 ≤ `limit`"的全部相关环,并返回它们张成的秩。
///
/// 秩达到圈秩即可停止加深:更长的环必然落在更短者的张成中,按定义不 relevant。
fn cycles_up_to(c: &Component, limit: usize, stats: &mut SearchStats) -> (Vec<Ring>, usize) {
    let mut cands = horton_candidates(c, limit, stats);
    // 同长度必须整组判定 —— 见模块文档
    cands.sort_by(|a, b| a.atoms.len().cmp(&b.atoms.len()).then(a.key.cmp(&b.key)));

    let mut basis = Gf2Basis::new(c.n_bonds());
    let mut rank = 0usize;
    let mut out: Vec<Ring> = Vec::new();
    let mut i = 0usize;
    while i < cands.len() {
        let len = cands[i].atoms.len();
        let mut j = i;
        while j < cands.len() && cands[j].atoms.len() == len {
            j += 1;
        }

        // 先全判:只对照**严格更短**的环张成
        for cand in &cands[i..j] {
            if !basis.reduce(&cand.bits).is_zero() {
                out.push(Ring {
                    atoms: cand.atoms.iter().map(|&a| c.atoms[a as usize]).collect(),
                    bonds: cand.bonds.iter().map(|&b| c.bonds[b as usize]).collect(),
                });
            }
        }
        // 再全加
        for cand in &cands[i..j] {
            if basis.insert(&cand.bits) {
                rank += 1;
            }
        }
        i = j;
    }
    (out, rank)
}

/// Horton 候选环,只生成长度 ≤ `limit` 的。
///
/// 每个顶点做一次**限深** BFS。深度上限取 `limit.div_ceil(2)`:Horton 定理里
/// 每个最小环基成员 C 都存在顶点 v ∈ C,使得 C 由 v 出发的两条**近乎等长**的
/// 最短路加一条边构成,故两条路各不超过 `⌈|C|/2⌉`。限深因此不会漏掉任何
/// 长度 ≤ limit 的相关环。
///
/// 这个限深是把代价从 O(分量点数 × 分量边数 × 路径长)压回来的关键。
fn horton_candidates(c: &Component, limit: usize, stats: &mut SearchStats) -> Vec<Candidate> {
    let n = c.n_atoms();
    let m = c.n_bonds();
    let depth_limit = limit.div_ceil(2) as u32;

    let mut seen: std::collections::BTreeSet<Vec<u32>> = std::collections::BTreeSet::new();
    let mut out: Vec<Candidate> = Vec::new();

    let mut dist = vec![u32::MAX; n];
    let mut parent = vec![u32::MAX; n];
    let mut parent_bond = vec![u32::MAX; n];
    // 每个点的最短路从 v 出发迈出的**第一步**是哪个邻居。
    // 两条路若第一步就相同,它们必然相交,不必回溯整条路径去发现这件事。
    let mut branch = vec![u32::MAX; n];
    let mut queue: Vec<u32> = Vec::with_capacity(n);
    // 标记某原子是否在当前考察的一条路径上,用于 O(路径长) 的相交判定
    let mut on_path = vec![false; n];
    // 本轮 BFS 内已考察过的边,避免同一条边从两个端点各来一次
    let mut edge_done = vec![false; m];

    for v in 0..n as u32 {
        // 复位只动上一轮碰过的下标 —— 整体复位是 O(点数),乘上外层循环就是平方
        for &x in &queue {
            dist[x as usize] = u32::MAX;
        }
        queue.clear();
        dist[v as usize] = 0;
        parent[v as usize] = u32::MAX;
        branch[v as usize] = u32::MAX;
        queue.push(v);
        let mut head = 0;
        while head < queue.len() {
            let x = queue[head];
            head += 1;
            stats.bfs_visits += 1;
            if dist[x as usize] >= depth_limit {
                continue; // 不再向外扩
            }
            for &(y, bi) in c.neighbors(x) {
                if dist[y as usize] == u32::MAX {
                    dist[y as usize] = dist[x as usize] + 1;
                    parent[y as usize] = x;
                    parent_bond[y as usize] = bi;
                    branch[y as usize] = if x == v { y } else { branch[x as usize] };
                    queue.push(y);
                }
            }
        }

        // 只看两端都在本轮 BFS 覆盖范围内的边。
        // 直接遍历 BFS 队列本身 —— 每个顶点克隆一次队列是 O(点数²) 的分配噪声。
        for &x in queue.iter() {
            for &(_, lb) in c.neighbors(x) {
                if edge_done[lb as usize] {
                    continue;
                }
                edge_done[lb as usize] = true;
                stats.edge_tests += 1;

                let (a, b) = c.edge_ends[lb as usize];
                if dist[a as usize] == u32::MAX || dist[b as usize] == u32::MAX {
                    continue;
                }
                let (da, db) = (dist[a as usize], dist[b as usize]);
                let len = (da + db + 1) as usize;
                if len < 3 || len > limit {
                    continue;
                }
                // 只取**平衡**的 (顶点, 边) 对:两条最短路长度差 ≤ 1。
                //
                // Horton 定理里,每个最小环基成员 C 都存在顶点 v ∈ C,使得 C 由
                // v 出发的两条近乎等长的最短路加"对面"那条边构成 —— 对面的边
                // 必然平衡。不平衡的组合要么给出同一个环(由别的 v 平衡地生成),
                // 要么根本不是简单环。
                //
                // # 但这一支**跑不到**,而且是结构上跑不到
                //
                // 无权 BFS 里,任意一条边的两个端点距离差 ≤ 1 —— 这是 BFS 的
                // 基本性质。限深只会让越界的点 `dist` 保持 `MAX`,而那种点上面
                // 那条 `dist == u32::MAX` 已经滤掉了。所以 `|da − db| ≥ 2`
                // 不可能出现。
                //
                // 实测:全语料 8830 个分子 + 五种合成形状(theta、并苯、大环、
                // 苯稠大环、共边双大环),这一支**一次都没进过**;把它整个关掉,
                // 三个工作量计数器在每一档上**逐位相同**。
                //
                // 这段注释先前写着"这个剪枝是把大环从立方压回来的关键 ……
                // 不加的话 200 原子的大环要 26 ms" —— 那个数是别的版本上量的,
                // 现在纯大环走的是开头那条快路径,压根到不了这里。
                //
                // 留着守卫 + 一条 `debug_assert`:哪天最短路换成**带权**的
                // (做"抑制二度顶点"那个修法时就要换成 Dijkstra),这个不变量
                // 就没了,而那时它必须重新起作用。
                debug_assert!(
                    da.abs_diff(db) <= 1,
                    "无权 BFS 里一条边两端的距离差应当 ≤ 1,实得 {da} 与 {db}"
                );
                if da.abs_diff(db) > 1 {
                    continue;
                }
                // 两条最短路的**第一步**若相同,它们必然相交,这个候选无效。
                //
                // 这一条是 O(1),而回溯出路径再判相交是 O(路径长)。少了它,
                // 纯大环里每个顶点都要为**同侧**的每条边白白回溯两次 ——
                // 整体从 O(n²) 退化成 O(n³):800 元大环要 1.1 秒。
                if branch[a as usize] == branch[b as usize] {
                    continue;
                }

                let px = trace(&parent, &parent_bond, v, a);
                let py = trace(&parent, &parent_bond, v, b);
                stats.path_steps += (px.len() + py.len()) as u64;
                // 两条路径除 v 外不得相交
                for &(t, _) in &px {
                    on_path[t as usize] = true;
                }
                let disjoint = py[..py.len() - 1]
                    .iter()
                    .all(|&(t, _)| !on_path[t as usize]);
                for &(t, _) in &px {
                    on_path[t as usize] = false;
                }
                if !disjoint {
                    continue;
                }

                // 环上顺序:a → … → v → … → b,再由边 (b,a) 闭合
                let mut atoms: Vec<u32> = px.iter().map(|&(t, _)| t).collect();
                atoms.extend(py[..py.len() - 1].iter().rev().map(|&(t, _)| t));
                let bonds = ring_bonds(&px, &py, lb);
                debug_assert_eq!(atoms.len(), len);
                debug_assert_eq!(bonds.len(), len);

                let mut key = atoms.clone();
                key.sort_unstable();
                if !seen.insert(key.clone()) {
                    continue;
                }

                let mut bits = BitSet::new(m);
                for &bd in &bonds {
                    bits.set(bd as usize);
                }
                out.push(Candidate {
                    atoms,
                    bonds,
                    bits,
                    key,
                });
            }
        }
        for &x in queue.iter() {
            for &(_, lb) in c.neighbors(x) {
                edge_done[lb as usize] = false;
            }
        }
    }
    out
}

/// 从 `t` 沿父指针回溯到 `v`,产出 `[(t, t到父的键), …, (v, u32::MAX)]`
fn trace(parent: &[u32], parent_bond: &[u32], v: u32, t: u32) -> Vec<(u32, u32)> {
    let mut out = Vec::new();
    let mut cur = t;
    while cur != v {
        out.push((cur, parent_bond[cur as usize]));
        cur = parent[cur as usize];
    }
    out.push((v, u32::MAX));
    out
}

/// 把两条路径与闭合边拼成与 `atoms`(a→…→v→…→b)对齐的键序列。
///
/// `bonds[i]` 必须连接 `atoms[i]` 与 `atoms[i+1]`。
fn ring_bonds(px: &[(u32, u32)], py: &[(u32, u32)], closing: u32) -> Vec<u32> {
    let mut bonds: Vec<u32> = px[..px.len() - 1].iter().map(|&(_, b)| b).collect();
    // py[..len-1] 逆序即 c₁…c_k(c₁ 是 v 的子、c_k 是另一端点);
    // 连接 c_i 与 c_{i+1} 的键正是 c_{i+1} 自己的 parent_bond,
    // 连接 v 与 c₁ 的是 c₁ 的 parent_bond —— 所以逐个压进来正好对齐。
    for &(_, b) in py[..py.len() - 1].iter().rev() {
        bonds.push(b);
    }
    bonds.push(closing);
    bonds
}

#[cfg(test)]
mod tests {
    use omgkit_io::smiles;
    use std::collections::BTreeSet;

    use super::*;

    /// 环集,表示成"原子集合的集合" —— 环内顺序是遍历产物,不参与比较
    fn ring_atom_sets(smi: &str) -> BTreeSet<BTreeSet<u32>> {
        let m = smiles::parse(smi).unwrap_or_else(|e| panic!("{}", e.render()));
        ring_set(&m)
            .into_iter()
            .map(|r| r.atoms.into_iter().collect())
            .collect()
    }

    fn ring_sizes(smi: &str) -> Vec<usize> {
        let mut v: Vec<usize> = ring_atom_sets(smi).iter().map(BTreeSet::len).collect();
        v.sort_unstable();
        v
    }

    #[test]
    fn acyclic_has_no_rings() {
        assert!(ring_sizes("CCO").is_empty());
        assert!(ring_sizes("CC(C)(C)CC").is_empty());
    }

    #[test]
    fn simple_rings() {
        assert_eq!(ring_sizes("C1CC1"), vec![3]);
        assert_eq!(ring_sizes("c1ccccc1"), vec![6]);
        assert_eq!(ring_sizes("C1CCCCCCC1"), vec![8]);
    }

    /// 稠环:萘 2 个、菲 3 个、芘 4 个
    #[test]
    fn fused_aromatics() {
        assert_eq!(ring_sizes("c1ccc2ccccc2c1"), vec![6, 6]);
        assert_eq!(ring_sizes("c1ccc2c(c1)ccc1ccccc12"), vec![6, 6, 6]);
        assert_eq!(ring_sizes("c1cc2ccc3cccc4ccc(c1)c2c34"), vec![6, 6, 6, 6]);
    }

    /// 环数**超过圈秩**的体系。
    ///
    /// 立方烷 8 原子 12 键 ⇒ 圈秩 5,但有 6 个四元面,彼此地位对等,
    /// 丢掉任何一个都没道理 —— 相关环的定义正是为此。
    #[test]
    fn cages_have_more_rings_than_the_cyclomatic_number() {
        for (smi, name, expect) in [
            ("C12C3C4C1C5C4C3C25", "立方烷", vec![4; 6]),
            ("C12C3C1C1C2C31", "棱晶烷", vec![3, 3, 4, 4, 4]),
            ("C1CC2CCC1CC2", "双环[2.2.2]辛烷", vec![6, 6, 6]),
            ("C1C2CC3CC1CC(C2)C3", "金刚烷", vec![6, 6, 6, 6]),
            ("C1CC2CCC1C2", "降冰片烷", vec![5, 5]),
        ] {
            assert_eq!(ring_sizes(smi), expect, "{name}");
        }
    }

    /// 环集必须与**原子编号**无关 —— 这是它能当规范量使用的前提。
    ///
    /// 做法是真的把原子重排一遍再算,而不是挑几个"看起来等价"的 SMILES:
    /// `C1CC2CCC1CC2`(双环[2.2.2])与 `C1CC2CCC(C1)CC2`(双环[3.2.1])
    /// 并不是同一个分子,环集本来就该不同。
    #[test]
    fn ring_set_is_independent_of_atom_numbering() {
        /// 按 `perm` 重排原子:新分子的第 i 个原子是旧分子的 `perm[i]`
        fn permute(m: &omgkit_core::MolBuilder, perm: &[u32]) -> omgkit_core::MolBuilder {
            let mut inv = vec![0u32; m.num_atoms()];
            for (new, &old) in perm.iter().enumerate() {
                inv[old as usize] = new as u32;
            }
            let mut out = omgkit_core::MolBuilder::with_capacity(m.num_atoms(), m.num_bonds());
            for &old in perm {
                out.add_atom_data(m.atoms()[old as usize]);
            }
            for b in m.bonds() {
                let mut nb = *b;
                nb.begin = inv[b.begin as usize];
                nb.end = inv[b.end as usize];
                out.add_bond_data(nb).expect("重排后端点仍合法");
            }
            out
        }

        for smi in [
            "c1ccccc1",
            "c1ccc2ccccc2c1",
            "C12C3C4C1C5C4C3C25",         // 立方烷:6 个面互相对等
            "C12C3C1C1C2C31",             // 棱晶烷
            "C1C2CC3CC1CC(C2)C3",         // 金刚烷
            "C1CC2CCC1CC2",               // 双环[2.2.2]辛烷
            "c1cc2ccc3cccc4ccc(c1)c2c34", // 芘
            "C1CCC2(CC1)CCCCC2",          // 螺环
            "c1ccccc1Cc1ccccc1",          // 两个独立环系
        ] {
            let m = smiles::parse(smi).unwrap_or_else(|e| panic!("{}", e.render()));
            let n = m.num_atoms() as u32;
            let base: BTreeSet<BTreeSet<u32>> = ring_set(&m)
                .into_iter()
                .map(|r| r.atoms.into_iter().collect())
                .collect();

            // 几种确定性的重排:反转、循环移位、奇偶交错
            let perms: Vec<Vec<u32>> = vec![
                (0..n).rev().collect(),
                (0..n).map(|i| (i + n / 2) % n).collect(),
                (0..n).step_by(2).chain((1..n).step_by(2)).collect(),
            ];
            for perm in perms {
                let pm = permute(&m, &perm);
                // 换回原始编号再比对
                let got: BTreeSet<BTreeSet<u32>> = ring_set(&pm)
                    .into_iter()
                    .map(|r| r.atoms.into_iter().map(|a| perm[a as usize]).collect())
                    .collect();
                assert_eq!(got, base, "{smi}:重排 {perm:?} 后环集改变");
            }
        }
    }

    /// 环上的原子与键必须真的首尾相接 —— 芳香性感知要靠这个顺序取键。
    #[test]
    fn ring_atoms_and_bonds_are_in_cycle_order() {
        for smi in [
            "c1ccccc1",
            "c1ccc2ccccc2c1",
            "C1C2CC3CC1CC(C2)C3",
            "C12C3C4C1C5C4C3C25",
            "C1CCC2(CC1)CCCCC2",
        ] {
            let m = smiles::parse(smi).unwrap();
            for r in ring_set(&m) {
                assert_eq!(r.atoms.len(), r.bonds.len(), "{smi}: 原子数应等于键数");
                for i in 0..r.len() {
                    let a = r.atoms[i];
                    let b = r.atoms[(i + 1) % r.len()];
                    let bond = m.bonds()[r.bonds[i] as usize];
                    assert!(
                        (bond.begin == a && bond.end == b) || (bond.begin == b && bond.end == a),
                        "{smi}: 键 {} 未连接 {a} 与 {b}",
                        r.bonds[i]
                    );
                }
                // 原子不能重复 —— 必须是简单环
                let uniq: BTreeSet<u32> = r.atoms.iter().copied().collect();
                assert_eq!(uniq.len(), r.atoms.len(), "{smi}: 环上原子重复");
            }
        }
    }

    /// 配位键不参与成环。
    #[test]
    fn dative_bonds_do_not_close_rings() {
        let mut m = smiles::parse("C1CCCCC1").unwrap();
        assert_eq!(ring_set(&m).len(), 1);
        m.bond_mut(0).unwrap().set_order(BondOrder::Dative);
        assert!(
            ring_set(&m).is_empty(),
            "把一条环键改成配位键后,该环不应再成立"
        );
    }

    /// 多个互不相连的环系各自独立计数。
    #[test]
    fn separate_ring_systems_are_all_found() {
        assert_eq!(ring_sizes("c1ccccc1.c1ccccc1"), vec![6, 6]);
        assert_eq!(ring_sizes("c1ccccc1Cc1ccccc1"), vec![6, 6]);
        assert_eq!(ring_sizes("C1CC1CCC1CCC1"), vec![3, 4]);
    }
}
