//! 规范化排序:给原子一个**与输入编号无关**的全序。
//!
//! 有了它,[`smiles::write_with_priority`](crate::smiles::write_with_priority)
//! 就能写出规范 SMILES —— 同一个分子无论原子怎么编号,都得到同一个字符串。
//!
//! # 判据可以自证
//!
//! 规范化是整条管线里唯一**不需要外部参照**就能验证的部分:把原子随机重排,
//! 规范秩必须不变。这条性质抓得住绝大多数排序错误,而且跑得起量。
//!
//! # 三件事:细化、立体、打破对称
//!
//! **一、颜色细化(1-WL)**:先按"与编号无关的原子属性"分格,再反复用
//! "邻居分布在哪些格里"细分,直到不再分裂。这一步是确定性的,得到的稳定
//! 划分与输入编号无关。
//!
//! **二、把立体信息也喂进细化**。纯图细化看不见手性,内消旋型的分子因此留下
//! 一个致命模糊:两个手性中心在图上完全等价,但互换它们的那个自同构**反转
//! 手性**。办法是把标记表达成"相对邻居等价类顺序"的宇称 —— 等价类与编号
//! 无关,这个宇称也就无关,可以当成新的原子属性再喂回去,做到不动点。
//!
//! **三、打破对称**:细化停下时若还有格含多个原子,说明它们在 1-WL 意义下
//! 不可区分。把第一个多元格的成员逐个试作起点,取字典序最小的串;试的过程中
//! 顺带发现自同构,同轨道的成员直接跳过。
//!
//! # 为什么不能朴素地一轮轮细化
//!
//! "每轮重算所有原子的新颜色"实现起来只要十几行,但轮数是**图的直径量级**。
//! 一条 11200 个原子的长链(语料里的"很多个苯环串起来"正是这个形状)要跑
//! 五千多轮,整体退化成平方。
//!
//! 本模块用的是分裂器工作表:每次只拿**一个格**当分裂器,只碰它的邻居。
//! 配合"分裂后把除最大块以外的都入表"这条规则(Hopcroft 的技巧),每个原子
//! 一生中最多进入 O(log n) 个分裂器,总代价 O((n + m) log n)。
//!
//! # 打破对称为什么要枚举而不是任取
//!
//! 一格之内的原子在 1-WL 意义下不可区分,但**未必真的等价** —— 稳定划分可以
//! 粗于自同构轨道。任取一个的话,取谁就成了输入编号的函数,同一个分子换个
//! 编号能得到两个不同的规范串。
//!
//! 这在真实分子上确实发生:8839 条语料里有 **7 条**取不同起点会写出不同的串
//! (用 [`tie_break_matters`] 可以重新量)。取最小值把这个自由度消掉,代价是
//! |第一格| 倍 —— 通常只有几个原子;高度对称的分子格大,但那时各分支写出的串
//! 本来就相同,多花的是重复功而不是错。
//!
//! 这个数**随写出能力增长**:规范串里每多一类信息(例如双键方向键),
//! 能区分起点的分子就多一批。写出得越细,这一步越必要。
//!
//! 更深层次的并列仍是任取。重排不变测试正是冲着这一点去的:一旦失败,
//! 就是这里要补。

use std::collections::BTreeSet;

use omgkit_core::{AtomFlags, MolBuilder};

/// 一个原子在某次分裂里的签名:连向分裂器的键型多重集(已排序)。
type Signature = Vec<u8>;
/// 待分裂的原子及其签名
type Touched = Vec<(u32, Signature)>;

/// 计算规范秩。返回 `rank[atom]`,取值是 `0..num_atoms` 的一个排列。
///
/// 秩越小越靠前。要写规范 SMILES 请直接用 [`canonical_smiles`] ——
/// 它还会处理立体标记,只拿秩去写会漏掉那一步。
#[must_use]
pub fn canonical_ranks(mol: &MolBuilder) -> Vec<u32> {
    if mol.num_atoms() == 0 {
        return Vec::new();
    }
    let mut p = Partition::new(mol);
    p.refine_with_stereo(mol);
    p.break_all_ties(mol);
    p.ranks()
}

/// 对称等价类。同类的原子在颜色细化(1-WL)意义下互不可区分。
///
/// 类编号本身没有意义,只有"是否同类"有意义。返回值对输入编号不敏感。
///
/// # 别拿它单独判"是不是真手性中心"
///
/// "两个邻居同类 ⇒ 标记没有内容"这条推断**不成立**。1,4-二取代环己烷的两个
/// 中心各自看邻居都不可区分,但两者合起来区分顺式与反式 —— 按那条推断会把
/// 顺反两个分子判成同一个。相互依赖的立体中心要用迭代判准,见
/// [`canonical_smiles`] 的说明。
///
/// 这里给出的是**必要条件而非充分条件**:同类即"单独看不可区分",不等于
/// "整体上无内容"。
#[must_use]
pub fn symmetry_classes(mol: &MolBuilder) -> Vec<u32> {
    if mol.num_atoms() == 0 {
        return Vec::new();
    }
    let mut p = Partition::new(mol);
    p.refine_with_stereo(mol);
    p.class_ids()
}

/// 写出规范 SMILES:同一个分子无论原子怎么编号,都得到同一个字符串。
///
/// # 无内容的立体标记会被抹掉
///
/// 四个甲基上的 `@` 表达不了任何东西 —— 换两个甲基得到的是同一个分子,标记
/// 却翻转了。这类标记由 [`stereo::genuine_tetrahedral`](crate::stereo::genuine_tetrahedral)
/// 判出来并抹掉。
///
/// 判准本身很难写对,而写错的代价是**不对称**的:
///
/// - 判松了,输出多一个没有内容的标记,分子仍然是对的
/// - 判严了,会抹掉**真的**立体信息,把两个不同的分子塌成同一个串
///
/// 所以那条判准刻意偏保守。单看"两个邻居同类"不足以判非真:1,4-二取代环己烷
/// 的两个中心各自看邻居都不可区分,合起来却区分顺式与反式,只按邻居同类判会
/// 把顺反写成同一个串。判准因此额外要求"等价支路里没有别的手性中心",
/// 正是为了放行这一对。
///
/// 规范性本身**不依赖**这一步:打破对称时取最小值已经保证了唯一性
/// (实测:抹与不抹,重排不变测试都通过)。抹掉只是让输出不带无意义的标记。
#[must_use]
pub fn canonical_smiles(mol: &MolBuilder) -> crate::smiles::Written {
    if mol.num_atoms() == 0 {
        return crate::smiles::write(mol);
    }
    // 抹掉没有内容的立体标记。判准偏保守,见本函数文档。
    let cleaned = drop_uninformative_stereo(mol);
    let mol = &cleaned;

    let mut base = Partition::new(mol);
    base.refine_with_stereo(mol);

    let Some(cell) = base.first_non_singleton() else {
        // 细化已经把所有原子分开,没有可挑的余地
        return crate::smiles::write_with_priority_styled(
            mol,
            &base.ranks(),
            crate::smiles::WriteStyle::Canonical,
        );
    };

    // 第一个多元格的成员逐个试作起点,取字典序最小的串。
    //
    // 这一格里的原子在 1-WL 意义下不可区分,但**未必真的等价** —— 稳定划分
    // 可以粗于自同构轨道。笼状多环就做得到:挑不同的原子起头,写出的串不同,
    // 于是同一个分子换个编号就有两个规范形式。取最小值把这个自由度消掉。
    //
    // # 靠发现自同构来剪枝
    //
    // 光是"每个成员都试一遍"会在高度对称的分子上退化:一个 n 元大环的第一格
    // 就是全部 n 个原子,于是要跑 n 遍,整体平方。而那 n 遍算的是同一个答案。
    //
    // 剪枝的依据来自枚举自身:两个起点若写出**同一个串**,把两次标号复合起来
    // 就得到一个自同构 —— 它把第一个起点映到第二个。该自同构轨道里的原子
    // 全都不必再试,因为从它们出发必然得到同一个串。大环因此只要试两次。
    let members: Vec<u32> = base.order[base.start[cell]..base.end[cell]].to_vec();
    let n = mol.num_atoms();
    let mut orbit = UnionFind::new(n);
    let mut tried: Vec<u32> = Vec::new();
    let mut best: Option<(crate::smiles::Written, Vec<u32>)> = None;

    for &a in &members {
        // 与某个试过的起点同轨道 —— 结果必然一样,跳过
        if tried.iter().any(|&t| orbit.same(t as usize, a as usize)) {
            continue;
        }
        let mut p = base.clone();
        p.split_off_atom(a, cell);
        p.refine_with_stereo(mol);
        p.break_all_ties(mol);
        let written = crate::smiles::write_with_priority_styled(
            mol,
            &p.ranks(),
            crate::smiles::WriteStyle::Canonical,
        );

        if let Some((prev, prev_order)) = &best {
            if prev.smiles == written.smiles {
                // 同串 ⇒ 两次标号复合出一个自同构,合并它的所有轮换
                for (x, y) in prev_order.iter().zip(p.order.iter()) {
                    orbit.union(*x as usize, *y as usize);
                }
            }
        }
        let replace = best
            .as_ref()
            .map_or(true, |(b, _)| written.smiles < b.smiles);
        if replace {
            best = Some((written, p.order.clone()));
        }
        tried.push(a);
    }
    best.expect("格非空").0
}

/// 抹掉不携带信息的四面体标记,返回处理后的分子。
///
/// [`canonical_smiles`] 与 [`tie_break_matters`] 共用 —— 两者必须看到**同一个**
/// 分子,否则后者量的就不是前者实际会遇到的情形。
fn drop_uninformative_stereo(mol: &MolBuilder) -> MolBuilder {
    let genuine = crate::stereo::genuine_tetrahedral(mol);
    let mut out = mol.clone();
    for (i, &g) in genuine.iter().enumerate() {
        if g {
            continue;
        }
        if let Some(a) = out.atom_mut(i as u32) {
            if a.chiral_tag.is_tetrahedral() {
                a.chiral_tag = omgkit_core::ChiralTag::Unspecified;
            }
        }
    }
    out
}

/// 诊断:打破对称这一步对这个分子**是否真的影响结果**。
///
/// 第一格里的原子在 1-WL 意义下不可区分,但未必真的等价。全都等价时,
/// 取哪个起点写出的串都一样,枚举取最小只是重复功;不全等价时,任取一个就会
/// 让规范串成为输入编号的函数 —— 那才是必须枚举的理由。
///
/// 这个函数回答的正是后一种情形出现了没有,[模块文档](self) 里那个"7 条"
/// 就是用它数出来的 —— 该说法可以随时重新量,而不是只能相信。
#[must_use]
pub fn tie_break_matters(mol: &MolBuilder) -> bool {
    if mol.num_atoms() == 0 {
        return false;
    }
    // 与 canonical_smiles 走同一条预处理,否则量的不是同一件事
    let cleaned = drop_uninformative_stereo(mol);
    let mol = &cleaned;
    let mut base = Partition::new(mol);
    base.refine_with_stereo(mol);
    let Some(cell) = base.first_non_singleton() else {
        return false;
    };
    let members: Vec<u32> = base.order[base.start[cell]..base.end[cell]].to_vec();
    let mut first: Option<String> = None;
    for &a in &members {
        let mut p = base.clone();
        p.split_off_atom(a, cell);
        p.refine_with_stereo(mol);
        p.break_all_ties(mol);
        let s = crate::smiles::write_with_priority_styled(
            mol,
            &p.ranks(),
            crate::smiles::WriteStyle::Canonical,
        )
        .smiles;
        match &first {
            None => first = Some(s),
            Some(f) if *f != s => return true,
            Some(_) => {}
        }
    }
    false
}

/// 并查集,用来把发现的自同构闭成轨道。
struct UnionFind(Vec<usize>);

impl UnionFind {
    fn new(n: usize) -> Self {
        Self((0..n).collect())
    }

    fn find(&mut self, mut x: usize) -> usize {
        while self.0[x] != x {
            self.0[x] = self.0[self.0[x]]; // 路径压缩
            x = self.0[x];
        }
        x
    }

    fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra != rb {
            self.0[ra] = rb;
        }
    }

    fn same(&mut self, a: usize, b: usize) -> bool {
        self.find(a) == self.find(b)
    }
}

/// 把四面体标记换算到"相对邻居等价类顺序"的参照系,得到一个与输入编号
/// 无关的取值。非四面体、或取代基不可区分时返回 0。
///
/// 这是让立体信息参与细化的关键:标记本身相对**存储序**,而存储序随建键
/// 顺序而变,直接拿来当不变量会把编号信息偷渡进规范化。
fn stereo_descriptor(mol: &MolBuilder, a: u32, classes: &[u32]) -> u8 {
    let at = mol.atoms()[a as usize];
    if !at.chiral_tag.is_tetrahedral() {
        return 0;
    }
    let nbrs: Vec<(u32, u32)> = mol
        .neighbors(a)
        .map(|(other, bond)| (classes[other as usize], bond))
        .collect();
    // 有两个邻居同类,标记就没有内容 —— 换这两个取代基得到的是同一个分子
    let mut cs: Vec<u32> = nbrs.iter().map(|&(c, _)| c).collect();
    cs.sort_unstable();
    if cs.windows(2).any(|w| w[0] == w[1]) {
        return 0;
    }

    let storage: Vec<u32> = nbrs.iter().map(|&(_, b)| b).collect();
    let mut by_class = nbrs;
    by_class.sort_unstable_by_key(|&(c, _)| c);
    let class_order: Vec<u32> = by_class.iter().map(|&(_, b)| b).collect();

    let odd = crate::smiles::permutation_is_odd(&storage, &class_order).unwrap_or(false);
    let tag = if odd {
        at.chiral_tag.inverted()
    } else {
        at.chiral_tag
    };
    tag as u8
}

/// 与输入编号无关的原子属性,用作细化的起点。
///
/// 只能放**分子决定的**量。放进任何随编号而变的东西(比如原子下标),
/// 整个规范化就失去意义,而且失效方式很隐蔽 —— 结果照样是个全序,
/// 只是换个编号就变了。
fn initial_invariant(mol: &MolBuilder, a: u32) -> (u8, u8, u32, i8, u8, u16, bool) {
    let at = mol.atoms()[a as usize];
    // 键级和乘 2 取整:芳香键的 1.5 才不会被截断
    let bond_sum2: u32 = mol
        .neighbors(a)
        .map(|(_, b)| (mol.bonds()[b as usize].order.as_double() * 2.0) as u32)
        .sum();
    (
        at.atomic_num,
        mol.degree(a) as u8,
        bond_sum2,
        at.formal_charge,
        at.num_explicit_hs.saturating_add(at.num_implicit_hs),
        at.isotope,
        at.flags.contains(AtomFlags::AROMATIC),
    )
}

/// 有序划分。格是 [`Partition::order`] 上的连续区间,格与格之间的先后
/// 就是最终的秩序。
#[derive(Clone)]
struct Partition {
    /// 全部原子,按格分段排列
    order: Vec<u32>,
    /// 原子 → 它在 `order` 里的下标
    pos: Vec<usize>,
    /// 原子 → 所属格
    cell_of: Vec<usize>,
    /// 格 → 区间起点(同时也定义了格的先后)
    start: Vec<usize>,
    /// 格 → 区间终点(不含)
    end: Vec<usize>,
    /// 待用作分裂器的格。按**区间起点**取,保证处理顺序也与编号无关。
    pending: BTreeSet<(usize, usize)>,
}

impl Partition {
    fn new(mol: &MolBuilder) -> Self {
        let n = mol.num_atoms();
        let mut order: Vec<u32> = (0..n as u32).collect();
        order.sort_by_key(|&a| initial_invariant(mol, a));

        let mut pos = vec![0usize; n];
        let mut cell_of = vec![0usize; n];
        let (mut start, mut end) = (Vec::new(), Vec::new());

        let mut i = 0;
        while i < n {
            let key = initial_invariant(mol, order[i]);
            let cell = start.len();
            let lo = i;
            while i < n && initial_invariant(mol, order[i]) == key {
                pos[order[i] as usize] = i;
                cell_of[order[i] as usize] = cell;
                i += 1;
            }
            start.push(lo);
            end.push(i);
        }

        let pending = start.iter().enumerate().map(|(c, &s)| (s, c)).collect();
        Self {
            order,
            pos,
            cell_of,
            start,
            end,
            pending,
        }
    }

    fn size(&self, c: usize) -> usize {
        self.end[c] - self.start[c]
    }

    /// 反复取分裂器细分,直到划分稳定。
    fn refine(&mut self, mol: &MolBuilder) {
        // 每个原子对每条邻边最多贡献一次签名,签名缓冲复用以免逐格重分配
        let mut sig: Touched = Vec::new();
        let mut mark: Vec<usize> = vec![usize::MAX; mol.num_atoms()];

        while let Some(&(_, splitter)) = self.pending.iter().next() {
            self.pending.remove(&(self.start[splitter], splitter));

            // 收集"与分裂器相邻"的原子,以及它们连过去的键型多重集
            sig.clear();
            for i in self.start[splitter]..self.end[splitter] {
                let x = self.order[i];
                for (nbr, bond) in mol.neighbors(x) {
                    let code = mol.bonds()[bond as usize].order as u8;
                    let slot = mark[nbr as usize];
                    if slot == usize::MAX {
                        mark[nbr as usize] = sig.len();
                        sig.push((nbr, vec![code]));
                    } else {
                        sig[slot].1.push(code);
                    }
                }
            }
            for (a, codes) in &mut sig {
                codes.sort_unstable();
                mark[*a as usize] = usize::MAX;
            }

            // 按所属格归拢,再逐格分裂
            let mut by_cell: Vec<(usize, Touched)> = Vec::new();
            let mut cell_slot: Vec<usize> = Vec::new();
            for (a, codes) in sig.drain(..) {
                let c = self.cell_of[a as usize];
                if cell_slot.len() <= c {
                    cell_slot.resize(c + 1, usize::MAX);
                }
                if cell_slot[c] == usize::MAX {
                    cell_slot[c] = by_cell.len();
                    by_cell.push((c, Vec::new()));
                }
                by_cell[cell_slot[c]].1.push((a, codes));
            }

            for (c, touched) in by_cell {
                self.split_cell(c, touched);
            }
        }
    }

    /// 把格 `c` 按签名分裂。未被触及的原子签名视作空,排在最前。
    ///
    /// 只搬动被触及的原子,代价 O(|touched| log |touched|) —— 若连未触及的
    /// 也要扫一遍,大格会把整体拖成平方。
    fn split_cell(&mut self, c: usize, mut touched: Touched) {
        if touched.len() == self.size(c) && touched.iter().all(|(_, s)| *s == touched[0].1) {
            return; // 整格同签名,不分裂
        }
        touched.sort_by(|a, b| a.1.cmp(&b.1).then(a.0.cmp(&b.0)));

        // 把被触及的原子搬到区间尾部。未处理的原子必然还在 boundary 之前,
        // 因为尾部区间里装的正好是已处理过的那些。
        let mut boundary = self.end[c];
        for &(a, _) in &touched {
            boundary -= 1;
            let pa = self.pos[a as usize];
            let moved = self.order[boundary];
            self.order[pa] = moved;
            self.pos[moved as usize] = pa;
            self.order[boundary] = a;
            self.pos[a as usize] = boundary;
        }
        // 尾部此刻正好是被触及的那些原子,按签名顺序重写一遍
        for (k, (a, _)) in touched.iter().enumerate() {
            self.order[boundary + k] = *a;
            self.pos[*a as usize] = boundary + k;
        }

        // 未触及的残留仍是格 c(可能为空);尾部按签名切成若干新格
        let old_end = self.end[c];
        let was_pending = self.pending.remove(&(self.start[c], c));
        let mut pieces: Vec<usize> = Vec::new();
        if boundary > self.start[c] {
            self.end[c] = boundary;
            pieces.push(c);
        }

        let mut k = 0;
        while k < touched.len() {
            let mut j = k + 1;
            while j < touched.len() && touched[j].1 == touched[k].1 {
                j += 1;
            }
            let lo = boundary + k;
            let hi = boundary + j;
            let cell = if pieces.is_empty() && boundary == self.start[c] {
                // 整格都被触及:第一块沿用原编号,免得留下空格
                self.end[c] = hi;
                c
            } else {
                self.start.push(lo);
                self.end.push(hi);
                self.start.len() - 1
            };
            for i in lo..hi {
                self.cell_of[self.order[i] as usize] = cell;
            }
            pieces.push(cell);
            k = j;
        }
        debug_assert_eq!(self.end[*pieces.last().expect("至少一块")], old_end);

        // Hopcroft:除最大块外全部入表。原格本就待处理时,所有块都要入表 ——
        // 它此前作为分裂器的效力还没兑现,不能被"最大块"这条规则吞掉。
        let largest = pieces
            .iter()
            .copied()
            .max_by_key(|&p| self.size(p))
            .expect("至少一块");
        for p in pieces {
            if was_pending || p != largest {
                self.pending.insert((self.start[p], p));
            }
        }
    }

    /// 细化到连立体信息也用尽。
    ///
    /// 纯图细化看不见手性,于是内消旋型的分子会留下一个致命的模糊:两个手性
    /// 中心在图上完全等价,但把它们互换的那个自同构**反转手性**。打破对称先
    /// 挑中哪一个,就决定了最后写出的是 `@` 还是 `@@` —— 同一个分子换个编号
    /// 就得到两个规范串。
    ///
    /// 出路是把标记表达成"相对**邻居等价类**顺序"的宇称。等价类与输入编号
    /// 无关,这个宇称因而也无关,可以当成一个新的原子属性再喂回细化。
    /// 一轮细化可能让等价类变细,于是要反复做到不动点。
    fn refine_with_stereo(&mut self, mol: &MolBuilder) {
        loop {
            self.refine(mol);
            let classes = self.class_ids();
            let cells: Vec<usize> = (0..self.start.len())
                .filter(|&c| self.size(c) > 1)
                .collect();
            let before = self.start.len();
            for c in cells {
                if self.size(c) <= 1 {
                    continue;
                }
                let touched: Touched = self.order[self.start[c]..self.end[c]]
                    .iter()
                    .map(|&a| (a, vec![stereo_descriptor(mol, a, &classes)]))
                    .collect();
                self.split_cell(c, touched);
            }
            if self.start.len() == before {
                return; // 没有格被立体信息分开,到不动点了
            }
        }
    }

    /// 逐个打破对称,直到每格只剩一个原子。
    fn break_all_ties(&mut self, mol: &MolBuilder) {
        while let Some(cell) = self.first_non_singleton() {
            self.split_off_first(cell);
            self.refine_with_stereo(mol);
        }
    }

    /// 每原子的最终秩(格全为单元素时才有意义)。
    fn ranks(&self) -> Vec<u32> {
        let mut r = vec![0u32; self.order.len()];
        for (i, &a) in self.order.iter().enumerate() {
            r[a as usize] = i as u32;
        }
        r
    }

    /// 每原子的等价类编号。用格的区间起点当编号 —— 与输入编号无关,
    /// 且同格的原子必然取到同一个值。
    fn class_ids(&self) -> Vec<u32> {
        self.cell_of
            .iter()
            .map(|&cell| self.start[cell] as u32)
            .collect()
    }

    /// 最靠前的多原子格。全是单元素格时返回 `None`。
    fn first_non_singleton(&self) -> Option<usize> {
        (0..self.start.len())
            .filter(|&c| self.size(c) > 1)
            .min_by_key(|&c| self.start[c])
    }

    /// 把格 `c` 的第一个原子单独提出来,排在该格之前。
    fn split_off_first(&mut self, c: usize) {
        self.split_off_atom(self.order[self.start[c]], c);
    }

    /// 把 `a` 从格 `c` 里单独提出来,排在该格之前。
    fn split_off_atom(&mut self, a: u32, c: usize) {
        debug_assert_eq!(self.cell_of[a as usize], c, "原子不在该格里");
        // 先把 a 换到格首,再切掉格首
        let lo = self.start[c];
        let pa = self.pos[a as usize];
        let head = self.order[lo];
        self.order[lo] = a;
        self.pos[a as usize] = lo;
        self.order[pa] = head;
        self.pos[head as usize] = pa;

        let cell = self.start.len();
        self.start.push(lo);
        self.end.push(lo + 1);
        self.cell_of[a as usize] = cell;
        self.start[c] = lo + 1;

        // 两块都要重新参与细化:提出来的那个成了新的分裂器,
        // 残留的那格也得重新算一遍与它的关系
        self.pending.insert((self.start[cell], cell));
        self.pending.insert((self.start[c], c));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles;

    fn ranks_of(smi: &str) -> Vec<u32> {
        let m = smiles::parse(smi).unwrap_or_else(|e| panic!("{smi}: {}", e.render()));
        canonical_ranks(&m)
    }

    /// 秩必须是 0..n 的一个排列 —— 有重复或有空缺都会让写出乱套。
    #[test]
    fn ranks_are_a_permutation() {
        for smi in [
            "C",
            "CCO",
            "c1ccccc1",
            "OC(=O)c1ccccc1N",
            "CCO.CCN",
            "C1CC2CCC1CC2",
            "CC(C)(C)C",
        ] {
            let r = ranks_of(smi);
            let mut sorted = r.clone();
            sorted.sort_unstable();
            let expect: Vec<u32> = (0..r.len() as u32).collect();
            assert_eq!(sorted, expect, "{smi} 的秩不是一个排列:{r:?}");
        }
    }

    /// 对称等价的原子会被细化归到同一格,靠打破对称才分开 —— 这条确认
    /// 打破对称确实在跑,而不是初始不变量就已经把所有原子分开了。
    #[test]
    fn symmetric_molecules_need_tie_breaking() {
        // 苯的六个碳在 1-WL 下完全不可区分
        let m = smiles::parse("c1ccccc1").unwrap();
        let mut p = Partition::new(&m);
        p.refine(&m);
        assert!(
            p.first_non_singleton().is_some(),
            "苯细化之后应当仍有多原子的格"
        );
        // 但最终的秩仍是完整的排列
        assert_eq!(canonical_ranks(&m).len(), 6);
    }

    /// 细化本身要有分辨力:甲苯的环上原子按到甲基的距离分开。
    #[test]
    fn refinement_separates_inequivalent_atoms() {
        let m = smiles::parse("Cc1ccccc1").unwrap();
        let mut p = Partition::new(&m);
        p.refine(&m);
        // 甲基碳、连接碳、邻、间、对 —— 五类
        let cells: std::collections::BTreeSet<usize> =
            (0..m.num_atoms()).map(|a| p.cell_of[a]).collect();
        assert_eq!(cells.len(), 5, "甲苯应细分成 5 类,实际 {}", cells.len());
    }
}
