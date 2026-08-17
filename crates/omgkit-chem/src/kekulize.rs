//! Kekulize —— 把芳香键还原为交替的单双键(净化第 5 步)。
//!
//! # 问题形式
//!
//! 这是芳香子图上的**完美匹配**:先算出哪些原子"还差一根双键"(候选原子),
//! 再给每个候选原子恰好配一根双键,两端都必须是候选原子且键上带芳香标志。
//! 配不上就回溯。
//!
//! # 解不唯一,但化学量唯一
//!
//! 一个芳香体系通常有多个同样合法的 Kekulé 式(萘有两个)。选中哪一个取决于
//! 遍历顺序,没有化学含义。本模块保证的是:
//!
//! - 芳香标志全部清除,不残留芳香键
//! - 每个原子的显式价与隐式氢正确(这些跨不同 Kekulé 式不变)
//! - 每个原子的双键数正确(完美匹配覆盖同一批顶点)
//! - kekulize 前后总价不变
//! - **确定性**:同一输入恒给同一结果,否则规范化输出无从稳定
//!
//! # 无解与搜索未尽是两回事
//!
//! 回溯有一个随体系规模伸缩的上限(`backtrack_budget`)。触到上限时返回
//! [`KekulizeError::SearchBudgetExhausted`],表示**结果未知**;只有搜索真正
//! 穷尽后才返回 [`KekulizeError::CannotKekulize`]。把这两者混为一谈,会把
//! "没找到"说成"不存在"。
//!
//! # 尚未支持:芳香体系中的通配原子
//!
//! 通配原子身份未知,可接可不接双键,在匹配里是"可选顶点"。正确的表述是
//! "求一个饱和全部必需顶点的匹配",一般图上需要 Blossom 算法。在实现之前
//! 显式拒绝([`KekulizeError::UnsupportedAromaticDummy`]),不给出可能错误的
//! 结果。
//!
//! 判据是"体系里**存在**通配原子",而不是"该通配原子自身芳香":
//! `c1cc[*]cc1` 里的 `*` 不带芳香标志,却同样参与配对,当成普通原子处理会把
//! 有解报成无解。
//!
//! # 芳香标志与键级是两个独立的量
//!
//! `mark_dbond_cands` 把芳香键的**键级**改成单键,但**保留芳香标志** ——
//! 搜索阶段正是靠这个标志判断"哪些键可以变双键"。标志要到最后
//! `mark_atoms_bonds` 才统一清除。两者混淆会让匹配无从下手。

use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::collections::VecDeque;

use omgkit_core::{element, AtomFlags, BondFlags, BondOrder, MolBuilder};

use crate::valence::{explicit_valence_nonstrict, implicit_hs_nonstrict, total_valence_nonstrict};

/// 回溯次数的安全上限,随芳香体系规模伸缩。
///
/// 分子的芳香体系通常只有几十个原子,穷尽搜索完全跑得动,所以这个上限只是
/// 防病态输入的安全网,取得足够大以致实际不会触发。真触到时返回
/// [`KekulizeError::SearchBudgetExhausted`],绝不冒充"无解"。
fn backtrack_budget(system_size: usize) -> usize {
    10_000 + 1_000 * system_size
}

/// Kekulize 失败。
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum KekulizeError {
    /// 搜索已**穷尽**,该芳香体系不存在合法的 Kekulé 结构
    CannotKekulize {
        /// 仍需要双键却配不上的原子
        atoms: Vec<u32>,
    },
    /// 触到安全上限,搜索未穷尽 —— 结果**未知**,既非"有解"也非"无解"。
    SearchBudgetExhausted {
        /// 该芳香体系的原子数
        system_size: usize,
        /// 已用的回溯次数
        backtracks: usize,
    },
    /// 不在任何环中的原子被标记为芳香
    NonRingAromaticAtom {
        /// 该原子
        atom: u32,
    },
    /// Kekulize 前后总价发生了变化 —— 内部一致性被破坏
    ValenceChanged {
        /// 该原子
        atom: u32,
        /// 之前的总价
        before: i32,
        /// 之后的总价
        after: i32,
    },
    /// 芳香体系中出现通配原子;该路径尚未实现,拒绝给出可能错误的结果。
    ///
    /// 通配原子身份未知,**可接可不接**双键,因此它是匹配里的"可选顶点"。
    /// 正确的表述是:求一个**饱和全部必需顶点**的匹配,可选顶点覆盖与否不限;
    /// 一般图上需要 Blossom 算法。在实现之前显式拒绝,不给出可能错误的结果 ——
    /// 当成普通原子硬算会把"有解"报成"无解"。
    UnsupportedAromaticDummy {
        /// 该芳香体系中的某个通配原子
        atom: u32,
    },
}

impl core::fmt::Display for KekulizeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::CannotKekulize { atoms } => {
                write!(f, "该芳香体系不存在合法 Kekulé 结构,未配平的原子:{atoms:?}")
            }
            Self::SearchBudgetExhausted {
                system_size,
                backtracks,
            } => write!(
                f,
                "kekulize 搜索触到安全上限({backtracks} 次回溯,体系 {system_size} 原子),\
                 结果未知 —— 这不等于无解"
            ),
            Self::NonRingAromaticAtom { atom } => {
                write!(f, "原子 #{atom} 不在环中却被标记为芳香")
            }
            Self::ValenceChanged {
                atom,
                before,
                after,
            } => write!(f, "kekulize 改变了原子 #{atom} 的总价:{before} → {after}"),
            Self::UnsupportedAromaticDummy { atom } => {
                write!(f, "原子 #{atom} 是芳香体系中的通配原子;该路径尚未实现")
            }
        }
    }
}

impl std::error::Error for KekulizeError {}

/// 把分子中的芳香键还原为交替单双键。
///
/// 调用前必须已完成第 3 步(价键)与第 4 步(环感知)——
/// 后者提供 [`AtomFlags::IN_RING`],用于"非环原子不得为芳香"的校验。
///
/// # Errors
/// 见 [`KekulizeError`]。
pub fn kekulize(mol: &mut MolBuilder) -> Result<(), KekulizeError> {
    // -- 是否有活可干 --
    let has_aromatic = mol
        .bonds()
        .iter()
        .any(|b| b.flags.contains(BondFlags::AROMATIC) || b.order == BondOrder::Aromatic)
        || (0..mol.num_atoms() as u32).any(|i| is_aromatic_atom(mol, i));
    if !has_aromatic {
        return Ok(());
    }

    // -- 总价快照,收尾时比对 --
    let before: Vec<i32> = (0..mol.num_atoms() as u32)
        .map(|i| total_valence_nonstrict(mol, i))
        .collect();

    // 一次分配,所有稠环体系共用 —— 见 `Scratch`
    let mut scratch = Scratch::new(mol.num_atoms(), mol.num_bonds());

    for system in crate::rings::fused_ring_systems(mol) {
        // 芳香体系里只要**存在**通配原子就拒绝,不论它自身是否芳香。
        // `c1cc[*]cc1` 里的 `*` 不芳香,但它同样参与配对(身份未知 →
        // 可接可不接双键);当成普通原子处理会把有解判成无解。
        if let Some(&d) = system
            .iter()
            .find(|&&a| mol.atoms()[a as usize].atomic_num == 0)
        {
            if system.iter().any(|&a| is_aromatic_atom(mol, a)) {
                return Err(KekulizeError::UnsupportedAromaticDummy { atom: d });
            }
        }
        kekulize_fused(mol, &system, &mut scratch)?;
    }
    scratch.debug_assert_clean();

    mark_atoms_bonds(mol)?;

    for i in 0..mol.num_atoms() as u32 {
        let after = total_valence_nonstrict(mol, i);
        if after != before[i as usize] {
            return Err(KekulizeError::ValenceChanged {
                atom: i,
                before: before[i as usize],
                after,
            });
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------

/// 原子自身带芳香标志,或任一关联键是芳香键。
fn is_aromatic_atom(mol: &MolBuilder, idx: u32) -> bool {
    if mol.atoms()[idx as usize]
        .flags
        .contains(AtomFlags::AROMATIC)
    {
        return true;
    }
    mol.neighbors(idx).any(|(_, bi)| {
        let b = mol.bonds()[bi as usize];
        b.flags.contains(BondFlags::AROMATIC) || b.order == BondOrder::Aromatic
    })
}

/// 该原子关联的键下标,按**键下标升序**产出。
///
/// 升序是邻接索引本身的性质:它按键的插入序串链。
fn bonds_of(mol: &MolBuilder, idx: u32) -> impl Iterator<Item = u32> + '_ {
    mol.neighbors(idx).map(|(_, bi)| bi)
}

/// 搜索结果三态 —— "没搜完"与"无解"必须分开。
enum SearchOutcome {
    Solved,
    NoSolution,
    BudgetExhausted { backtracks: usize },
}

/// 跨稠环体系复用的暂存空间。
///
/// 一个分子可以有多个互不相连的稠环体系,多环化合物、聚合物尤其如此。
/// 若每个体系都新分配一套按原子数/键数定长的数组,总开销就是
/// O(体系数 × 分子规模) —— 一个隐式的平方项。
///
/// 所以数组只分配一次,每个体系用完后**只清理自己碰过的下标**
/// ([`Scratch::clear`]),清理量正比于该体系的规模而非分子的规模。
struct Scratch {
    /// 可接双键的原子(按原子下标)
    cands: Vec<bool>,
    /// 原子是否属于当前体系
    in_all: Vec<bool>,
    /// 本体系已加过双键的键
    dbnd_adds: Vec<bool>,
    /// 本轮已加过双键的键 —— 判定无解时据此整体回滚
    local_added: Vec<bool>,
    /// 每原子在"已处理"序列中的出现次数
    done_count: Vec<u32>,
    /// 每原子在待处理队列中的出现次数
    astack_count: Vec<u32>,
}

impl Scratch {
    fn new(n_atoms: usize, n_bonds: usize) -> Self {
        Self {
            cands: vec![false; n_atoms],
            in_all: vec![false; n_atoms],
            dbnd_adds: vec![false; n_bonds],
            local_added: vec![false; n_bonds],
            done_count: vec![0; n_atoms],
            astack_count: vec![0; n_atoms],
        }
    }

    /// 把本体系碰过的下标恢复成初始值,交给下一个体系。
    ///
    /// 覆盖面的论证:按原子寻址的四个数组只会写 `all_atms` 里的原子
    /// (邻居入队前先过 `in_all` 过滤);按键寻址的两个数组只会写体系内
    /// 两原子之间的键,那些键必定关联 `all_atms` 中的原子。所以扫一遍
    /// `all_atms` 及其关联键就够了。
    ///
    /// 漏掉任何一项都会造成体系之间的**静默串扰**:前一个体系留下的候选
    /// 标记会让后一个体系少配或多配一根双键。
    fn clear(&mut self, mol: &MolBuilder, all_atms: &[u32]) {
        for &a in all_atms {
            let i = a as usize;
            self.cands[i] = false;
            self.in_all[i] = false;
            self.done_count[i] = 0;
            self.astack_count[i] = 0;
            for (_, bi) in mol.neighbors(a) {
                self.dbnd_adds[bi as usize] = false;
                self.local_added[bi as usize] = false;
            }
        }
    }

    /// 调试期自检:全部体系处理完后,暂存空间必须回到初始状态。
    ///
    /// 刻意放在末尾而不是每个体系入口 —— 后者是 O(体系数 × 分子规模),
    /// 会让 debug 构建自己变成平方,而测试正跑在 debug 上。漏清理的项会一直
    /// 留到最后(`clear` 只会写回默认值,不会再置位),所以末尾查同样抓得到。
    fn debug_assert_clean(&self) {
        debug_assert!(self.cands.iter().all(|&x| !x), "cands 未清理干净");
        debug_assert!(self.in_all.iter().all(|&x| !x), "in_all 未清理干净");
        debug_assert!(self.dbnd_adds.iter().all(|&x| !x), "dbnd_adds 未清理干净");
        debug_assert!(
            self.local_added.iter().all(|&x| !x),
            "local_added 未清理干净"
        );
        debug_assert!(
            self.done_count.iter().all(|&c| c == 0),
            "done_count 未清理干净"
        );
        debug_assert!(
            self.astack_count.iter().all(|&c| c == 0),
            "astack_count 未清理干净"
        );
    }
}

fn kekulize_fused(
    mol: &mut MolBuilder,
    all_atms: &[u32],
    scratch: &mut Scratch,
) -> Result<(), KekulizeError> {
    let done = mark_dbond_cands(mol, all_atms, scratch);
    let outcome = kekulize_worker(mol, all_atms, done, scratch);

    // 错误载荷要在清理之前取出来
    let result = match outcome {
        SearchOutcome::Solved => Ok(()),
        SearchOutcome::NoSolution => {
            let mut atoms: Vec<u32> = all_atms
                .iter()
                .copied()
                .filter(|&i| scratch.cands[i as usize])
                .collect();
            atoms.sort_unstable(); // 载荷必须确定,与体系内的遍历顺序无关
            Err(KekulizeError::CannotKekulize { atoms })
        }
        SearchOutcome::BudgetExhausted { backtracks } => {
            Err(KekulizeError::SearchBudgetExhausted {
                system_size: all_atms.len(),
                backtracks,
            })
        }
    };

    scratch.clear(mol, all_atms);
    result
}

/// 标出"还差一根双键"的候选原子(仅非通配原子路径)。
///
/// 返回已判定无需处理的原子;候选原子写进 `scratch.cands`。同时把芳香键的
/// **键级**改成单键,芳香**标志**保留 —— 搜索阶段要靠标志判断哪些键可变双键。
fn mark_dbond_cands(mol: &mut MolBuilder, all_atms: &[u32], scratch: &mut Scratch) -> Vec<u32> {
    let mut done: Vec<u32> = Vec::new();

    let has_aromatic = all_atms
        .iter()
        .any(|&a| mol.atoms()[a as usize].atomic_num == 0 || is_aromatic_atom(mol, a));
    if !has_aromatic {
        return done;
    }

    let mut make_single: Vec<u32> = Vec::new();

    for &a in all_atms {
        let atom = mol.atoms()[a as usize];
        if atom.atomic_num != 0 && !is_aromatic_atom(mol, a) {
            done.push(a);
            continue;
        }

        let mut sbo: i32 = 0;
        let mut n_to_ignore: i32 = 0;
        for bi in bonds_of(mol, a) {
            let b = mol.bonds()[bi as usize];
            let aromatic_flagged = b.flags.contains(BondFlags::AROMATIC)
                && matches!(
                    b.order,
                    BondOrder::Single | BondOrder::Double | BondOrder::Aromatic
                );
            if aromatic_flagged {
                sbo += 1;
                make_single.push(bi);
            } else {
                let contrib = b.valence_contribution_to(a).round() as i32;
                sbo += contrib;
                if contrib == 0 {
                    n_to_ignore += 1;
                }
            }
        }

        // -- 非通配原子:判断能否接双键 --
        let ev = explicit_valence_nonstrict(mol, a);
        let implicit = i32::from(implicit_hs_nonstrict(mol, a, ev));
        let total_hs = i32::from(atom.num_explicit_hs) + implicit;
        sbo += total_hs;

        let z = atom.atomic_num;
        let valens = element::by_atomic_num(z).map_or(&[-1i8][..], |e| e.valences);
        let mut dv = valens.first().map_or(-1, |&v| i32::from(v));

        let mut chrg = i32::from(atom.formal_charge);
        if element::is_early_atom(z) {
            chrg = -chrg; // 周期表中位于碳左侧的元素,电荷对价的影响方向相反
        }
        if z == 6 && chrg > 0 {
            chrg = -chrg; // 带正电的碳同理
        }
        dv += chrg;

        let tbo = ev + implicit;
        // 第 6 步 FINDRADICALS 排在 kekulize 之后,故管线内此值必为 0。
        // 仍读字段而非写死 —— 见 `valence` 模块文档中同样的理由。
        let n_radicals = i32::from(atom.num_radical_electrons);
        let degree = mol.degree(a) as i32;
        let total_degree = degree + implicit - n_to_ignore;

        let mut vi = 1usize;
        while tbo > dv && vi < valens.len() && valens[vi] > 0 {
            dv = i32::from(valens[vi]) + chrg;
            vi += 1;
        }

        // 芳香 N-氧化物(如 O=n1ccccc1)—— 只有在跳过第 1 步时才会走到这里
        if tbo == 5
            && sbo == 4
            && dv == 3
            && total_degree == 3
            && n_radicals == 0
            && chrg == 0
            && total_hs == 0
            && matches!(z, 7 | 15 | 33)
        {
            dv = 5;
        }

        if total_degree + n_radicals >= dv {
            continue;
        }
        // 第一项:当前键序 + 1 恰好补齐价态。
        // 第二项:目前没有自由基,但若允许一个,该原子同样能接双键 ——
        // 仅对显式指定了氢数的原子成立。
        let can_take_double = dv == sbo + 1 + n_radicals
            || (n_radicals == 0 && atom.flags.contains(AtomFlags::NO_IMPLICIT) && dv == sbo + 2);
        if can_take_double {
            scratch.cands[a as usize] = true;
        }
    }

    // 注意:只改键级,保留芳香标志 —— 匹配阶段要靠它判断哪些键可以变双键
    for bi in make_single {
        if let Some(mut b) = mol.bond_mut(bi) {
            b.set_order(BondOrder::Single);
        }
    }
    done
}

/// 已处理原子的序列。
///
/// 顺序必须保留 —— 回溯要按"首次/末次出现位置"截断。同时维护一份出现
/// **次数**,让"是否已处理"从 O(n) 线性扫描降到 O(1) 查表。
///
/// 为什么是次数而不是布尔:同一个原子会多次入列 —— 回溯把它们推回待处理
/// 队列后又会被重新处理一遍,`back_track` 里用 `rposition` 正是因为这个。
struct DoneList<'a> {
    order: Vec<u32>,
    /// 借自 [`Scratch`];进入时必须全零
    count: &'a mut [u32],
}

impl<'a> DoneList<'a> {
    fn new(count: &'a mut [u32], initial: Vec<u32>) -> Self {
        for &a in &initial {
            count[a as usize] += 1;
        }
        Self {
            order: initial,
            count,
        }
    }

    fn push(&mut self, a: u32) {
        self.order.push(a);
        self.count[a as usize] += 1;
    }

    fn contains(&self, a: u32) -> bool {
        self.count[a as usize] > 0
    }

    fn len(&self) -> usize {
        self.order.len()
    }

    /// 只保留前 `keep` 个
    fn truncate(&mut self, keep: usize) {
        for &a in &self.order[keep..] {
            self.count[a as usize] -= 1;
        }
        self.order.truncate(keep);
    }
}

/// 待处理原子队列,同样带出现次数 —— 理由见 [`DoneList`]。
struct AtomStack<'a> {
    queue: VecDeque<u32>,
    /// 借自 [`Scratch`];进入时必须全零
    count: &'a mut [u32],
}

impl<'a> AtomStack<'a> {
    fn new(count: &'a mut [u32]) -> Self {
        Self {
            queue: VecDeque::new(),
            count,
        }
    }

    fn pop_front(&mut self) -> Option<u32> {
        let a = self.queue.pop_front()?;
        self.count[a as usize] -= 1;
        Some(a)
    }

    fn push_front(&mut self, a: u32) {
        self.queue.push_front(a);
        self.count[a as usize] += 1;
    }

    fn push_back(&mut self, a: u32) {
        self.queue.push_back(a);
        self.count[a as usize] += 1;
    }

    fn contains(&self, a: u32) -> bool {
        self.count[a as usize] > 0
    }

    fn is_empty(&self) -> bool {
        self.queue.is_empty()
    }
}

/// 匹配搜索主循环。返回三态结果:配平成功 / 确定无解 / 触到搜索上限。
fn kekulize_worker(
    mol: &mut MolBuilder,
    all_atms: &[u32],
    done: Vec<u32>,
    scratch: &mut Scratch,
) -> SearchOutcome {
    // 拆开借用:这几个数组彼此独立,可以同时持有可变引用
    let Scratch {
        cands,
        in_all,
        dbnd_adds,
        local_added,
        done_count,
        astack_count,
    } = scratch;

    let mut done = DoneList::new(done_count, done);
    let mut astack = AtomStack::new(astack_count);
    let mut options: HashMap<u32, VecDeque<u32>> = HashMap::new();
    let mut btmoves: Vec<u32> = Vec::new();
    let mut last_opt: Option<u32> = None;
    let mut num_bt = 0usize;
    let budget = backtrack_budget(all_atms.len());

    for &a in all_atms {
        in_all[a as usize] = true;
    }
    let mut sorted_atms = all_atms.to_vec();
    sorted_atms.sort_unstable();

    while done.len() < sorted_atms.len() || !astack.is_empty() {
        let curr = match astack.pop_front() {
            Some(c) => c,
            None => match sorted_atms.iter().copied().find(|&a| !done.contains(a)) {
                Some(c) => c,
                None => break,
            },
        };
        done.push(curr);

        let c_cand = cands[curr as usize];
        let mut opts: VecDeque<u32> = match options.get(&curr) {
            Some(o) => o.clone(),
            None => {
                let mut nbrs: Vec<u32> = mol
                    .neighbors(curr)
                    .map(|(x, _)| x)
                    .filter(|&x| in_all[x as usize] && !done.contains(x))
                    .collect();
                nbrs.sort_unstable();
                nbrs.dedup();

                let mut lstack: Vec<u32> = Vec::new();
                let mut o: VecDeque<u32> = VecDeque::new();
                for nbr in nbrs {
                    if !astack.contains(nbr) {
                        lstack.push(nbr);
                    }
                    if c_cand && cands[nbr as usize] {
                        let bi = mol.bond_between(curr, nbr).expect("邻居必有键");
                        // 只有带芳香标志的键才可以变成双键
                        if mol.bonds()[bi as usize].flags.contains(BondFlags::AROMATIC) {
                            o.push_back(nbr);
                        }
                    }
                }
                for a in lstack {
                    astack.push_back(a);
                }
                o
            }
        };

        if !c_cand {
            continue;
        }

        if let Some(ncnd) = opts.pop_front() {
            let bi = mol.bond_between(curr, ncnd).expect("选项必有键");
            if let Some(mut b) = mol.bond_mut(bi) {
                b.set_order(BondOrder::Double);
                b.set_direction(omgkit_core::BondDirection::None);
            }
            cands[curr as usize] = false;
            cands[ncnd as usize] = false;
            dbnd_adds[bi as usize] = true;
            local_added[bi as usize] = true;

            match options.entry(curr) {
                Entry::Occupied(mut e) => {
                    if opts.is_empty() {
                        e.remove();
                        btmoves.pop();
                        last_opt = btmoves.last().copied();
                    } else {
                        e.insert(opts);
                    }
                }
                Entry::Vacant(e) => {
                    // 还有别的选项没试 —— 记下回溯点
                    if !opts.is_empty() {
                        last_opt = Some(curr);
                        btmoves.push(curr);
                        e.insert(opts);
                    }
                }
            }
        } else if mol.atoms()[curr as usize].atomic_num != 0 {
            // 判定失败时必须撤销本次留下的所有改动,否则分子被改坏
            let undo = |mol: &mut MolBuilder| {
                for (bi, &added) in local_added.iter().enumerate() {
                    if added {
                        if let Some(mut b) = mol.bond_mut(bi as u32) {
                            b.set_order(BondOrder::Single);
                        }
                    }
                }
            };
            match last_opt {
                // 还有未试过的分支 —— 继续搜
                Some(lo) if num_bt < budget => {
                    back_track(mol, lo, &mut done, &mut astack, cands, dbnd_adds);
                    num_bt += 1;
                }
                // 没有任何未试分支了:搜索已穷尽,确实无解
                None => {
                    undo(mol);
                    return SearchOutcome::NoSolution;
                }
                // 触到安全上限:结果未知,不能冒充"无解"
                Some(_) => {
                    undo(mol);
                    return SearchOutcome::BudgetExhausted { backtracks: num_bt };
                }
            }
        }
    }
    SearchOutcome::Solved
}

/// 退回到 `last_opt` 那次选择之前的状态。
fn back_track(
    mol: &mut MolBuilder,
    last_opt: u32,
    done: &mut DoneList,
    astack: &mut AtomStack,
    cands: &mut [bool],
    dbnd_adds: &mut [bool],
) {
    // tdone 用**首次**出现位置截断;回栈范围用**末次**出现位置 —— 与 C++ 一致
    let first = done.order.iter().position(|&x| x == last_opt).unwrap_or(0);
    let last = done
        .order
        .iter()
        .rposition(|&x| x == last_opt)
        .unwrap_or(done.len().saturating_sub(1));

    // 必须先把 [last..] 推回待处理队列,再截断 —— 顺序反了就读不到了
    for &a in done.order[last..].iter().rev() {
        astack.push_front(a);
    }
    done.truncate(first);

    // 若某端已在 lastOpt 之前处理完,这根双键不必撤销。
    // 截断之后 `done` 就是原来的 tdone,查表即可。
    let to_undo: Vec<usize> = dbnd_adds
        .iter()
        .enumerate()
        .filter(|&(_, &added)| added)
        .map(|(bi, _)| bi)
        .filter(|&bi| {
            let b = mol.bonds()[bi];
            !done.contains(b.begin) && !done.contains(b.end)
        })
        .collect();
    for bi in to_undo {
        let b = mol.bonds()[bi];
        dbnd_adds[bi] = false;
        if let Some(mut bm) = mol.bond_mut(bi as u32) {
            bm.set_order(BondOrder::Single);
        }
        cands[b.begin as usize] = true;
        cands[b.end as usize] = true;
    }
}

/// 收尾:清除全部芳香标志,并把吡咯型氮的显式氢转回隐式。
fn mark_atoms_bonds(mol: &mut MolBuilder) -> Result<(), KekulizeError> {
    for bi in 0..mol.num_bonds() as u32 {
        if let Some(mut b) = mol.bond_mut(bi) {
            b.flags_mut().remove(BondFlags::AROMATIC);
        }
    }

    for i in 0..mol.num_atoms() as u32 {
        let atom = mol.atoms()[i as usize];
        if !atom.flags.contains(AtomFlags::AROMATIC) {
            continue;
        }
        if !atom.flags.contains(AtomFlags::IN_RING) {
            return Err(KekulizeError::NonRingAromaticAtom { atom: i });
        }
        let fix_pyrrole = matches!(atom.atomic_num, 7 | 15)
            && atom.formal_charge == 0
            && atom.num_explicit_hs == 1;
        if let Some(a) = mol.atom_mut(i) {
            a.flags.remove(AtomFlags::AROMATIC);
            if fix_pyrrole {
                // 吡咯型 [nH] 的显式氢要转回隐式,否则它会一直挂在那里,
                // 影响后续的价键计算与输出
                a.flags.remove(AtomFlags::NO_IMPLICIT);
                a.num_explicit_hs = 0;
            }
        }
        if fix_pyrrole {
            let ev = explicit_valence_nonstrict(mol, i);
            let ih = implicit_hs_nonstrict(mol, i, ev);
            if let Some(a) = mol.atom_mut(i) {
                a.num_implicit_hs = ih;
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use omgkit_core::BondFlags;
    use omgkit_io::smiles;

    use super::*;
    use crate::{clean_up, perceive_rings, update_property_cache};

    /// 跑到第 5 步为止的管线。
    fn pipeline(smi: &str) -> Result<MolBuilder, KekulizeError> {
        let mut m = smiles::parse(smi).unwrap_or_else(|e| panic!("{}", e.render()));
        clean_up(&mut m);
        update_property_cache(&mut m).expect("价键校验应通过");
        let _ = perceive_rings(&mut m);
        kekulize(&mut m)?;
        Ok(m)
    }

    fn orders(smi: &str) -> Vec<BondOrder> {
        pipeline(smi)
            .unwrap()
            .bonds()
            .iter()
            .map(|b| b.order)
            .collect()
    }

    // -- 结构不变量 --

    /// 多个稠环体系共用一份 [`Scratch`],所以"每个体系用完只清理自己碰过的
    /// 下标"必须**清干净**。漏一项就会造成体系之间的静默串扰 —— 前一个环
    /// 留下的候选标记会让后一个环少配或多配一根双键。
    ///
    /// 判据:每个苯环恰好 3 根双键,N 个环就是 3N 根,且每个芳香碳恰好
    /// 摊到 1 根。串扰会立刻打破这个计数。
    #[test]
    fn independent_ring_systems_do_not_contaminate_each_other() {
        for n in [1usize, 2, 3, 5, 8, 13, 30] {
            let smi = vec!["c1ccccc1"; n].join("C");
            let m = pipeline(&smi).unwrap_or_else(|e| panic!("{n} 个环: {e}"));

            let doubles = m
                .bonds()
                .iter()
                .filter(|b| b.order == BondOrder::Double)
                .count();
            assert_eq!(doubles, 3 * n, "{n} 个苯环应共有 {} 根双键", 3 * n);

            // 每个环碳恰好摊到一根双键 —— 连接用的 CH2 一根都不该有
            for a in 0..m.num_atoms() as u32 {
                let d = m
                    .neighbors(a)
                    .filter(|&(_, bi)| m.bonds()[bi as usize].order == BondOrder::Double)
                    .count();
                let expect = usize::from(m.degree(a) >= 2 && m.atoms()[a as usize].atomic_num == 6);
                // 连接碳的度也是 2,用"是否在环中"区分
                let in_ring = m.atoms()[a as usize]
                    .flags
                    .contains(omgkit_core::AtomFlags::IN_RING);
                let expect = if in_ring { expect } else { 0 };
                assert_eq!(d, expect, "{n} 个环:原子 #{a} 的双键数不对");
            }
        }
    }

    /// 用**独立分子**的结果作对照:把 N 个苯环放进同一个分子(片段形式)后,
    /// 每个环得到的结构必须与单独 kekulize 一模一样。
    ///
    /// 注意不能按键下标切块 —— 解析器把**所有片段**的环闭合键统一追加到
    /// 末尾(按环标号排序),键表是交错的。要按端点所属片段归组。
    #[test]
    fn ring_system_result_is_independent_of_siblings() {
        const RING: usize = 6;
        let alone = orders("c1ccccc1");
        assert_eq!(alone.len(), RING);

        for n in [2usize, 4, 9] {
            let smi = vec!["c1ccccc1"; n].join(".");
            let m = pipeline(&smi).unwrap();
            for k in 0..n {
                let lo = (k * RING) as u32;
                let hi = lo + RING as u32;
                let ours: Vec<BondOrder> = m
                    .bonds()
                    .iter()
                    .filter(|b| (lo..hi).contains(&b.begin))
                    .map(|b| b.order)
                    .collect();
                assert_eq!(
                    ours, alone,
                    "{n} 个片段:第 {k} 个环的结果与单独 kekulize 不同"
                );
            }
        }
    }

    /// **证明无解**才是回溯法的最坏情形 —— 找一个解通常随便撞就中,
    /// 而要断言不存在解就得把搜索空间走遍。
    ///
    /// 这里喂的是一族逐渐变大的**无解**稠合芳香体系:线性并苯骨架,总原子数
    /// 为奇,故必无完美匹配。
    ///
    /// 判据有二:
    ///
    /// 1. 必须返回 [`KekulizeError::CannotKekulize`](搜索已穷尽,确实无解),
    ///    **不能**是 `SearchBudgetExhausted`(没搜完,结果未知)—— 后者等于把
    ///    一个能判定的问题判成了"不知道"
    /// 2. 规模翻倍,代价不该爆炸
    ///
    /// 奇偶约束沿并苯带局部传播,很快就能证否,代价随规模线性增长。
    /// 这条测试盯着这个性质别悄悄失效。
    #[test]
    fn insoluble_systems_are_proven_insoluble_not_abandoned() {
        // 线性并苯骨架,4n+3 个碳 —— 原子数为奇,必无完美匹配
        fn odd_fused(n: usize) -> String {
            let rn = |k: usize| {
                if k < 10 {
                    k.to_string()
                } else {
                    format!("%{k}")
                }
            };
            let mut s = String::from("c1ccc2");
            for k in 2..n {
                s.push_str(&format!("cc{}", rn(k + 1)));
            }
            s.push_str(&format!("ccccc{}", rn(n)));
            for k in (2..n).rev() {
                s.push_str(&format!("cc{}", rn(k)));
            }
            s.push_str("cc1");
            s
        }

        for n in [2usize, 3, 4, 8, 16, 32] {
            let smi = odd_fused(n);
            let err = pipeline(&smi).expect_err(&format!("n={n}:该体系原子数为奇,不该能 kekulize"));
            assert!(
                matches!(err, KekulizeError::CannotKekulize { .. }),
                "n={n}({smi}):应判定为确定无解,实际 {err}"
            );
        }
    }

    /// 芳香标志必须**全部**清除,键级里也不能残留 `Aromatic`。
    #[test]
    fn no_aromatic_residue() {
        for smi in [
            "c1ccccc1",
            "c1ccncc1",
            "c1cc[nH]c1",
            "c1ccc2ccccc2c1",
            "CN1C=NC2=C1C(=O)N(C)C(=O)N2C",
            "CC(=O)Oc1ccccc1C(=O)O",
            "c1ccc2ccccc2c1.c1ccccc1",
        ] {
            let m = pipeline(smi).unwrap();
            assert!(
                m.bonds()
                    .iter()
                    .all(|b| b.order != BondOrder::Aromatic
                        && !b.flags.contains(BondFlags::AROMATIC)),
                "{smi}: 仍有芳香键"
            );
            assert!(
                m.atoms()
                    .iter()
                    .all(|a| !a.flags.contains(AtomFlags::AROMATIC)),
                "{smi}: 仍有芳香原子"
            );
        }
    }

    /// **确定性**:同一输入必须给出同一结果。
    /// 没有这一条,L3 的规范化输出不可能稳定。
    #[test]
    fn is_deterministic() {
        for smi in [
            "c1ccccc1",
            "c1ccc2ccccc2c1",
            "c1ccc2c(c1)ccc1ccccc12",
            "CN1C=NC2=C1C(=O)N(C)C(=O)N2C",
        ] {
            let first = orders(smi);
            for _ in 0..5 {
                assert_eq!(orders(smi), first, "{smi}: 结果不确定");
            }
        }
    }

    /// kekulize 不得改变任何原子的总价 —— 这是"结构合法"的充要条件之一。
    /// 实现内部已有此校验,这里从外部再证一次。
    #[test]
    fn total_valence_is_preserved() {
        for smi in ["c1ccccc1", "c1cc[nH]c1", "c1ccncc1", "c1ccc2ccccc2c1"] {
            let mut m = smiles::parse(smi).unwrap();
            clean_up(&mut m);
            update_property_cache(&mut m).unwrap();
            let _ = perceive_rings(&mut m);
            let before: Vec<i32> = (0..m.num_atoms() as u32)
                .map(|i| total_valence_nonstrict(&m, i))
                .collect();
            kekulize(&mut m).unwrap();
            let after: Vec<i32> = (0..m.num_atoms() as u32)
                .map(|i| total_valence_nonstrict(&m, i))
                .collect();
            assert_eq!(before, after, "{smi}: 总价被改变");
        }
    }

    /// 苯环:6 个碳应交替出现 3 根双键、3 根单键。
    /// **不断言是哪 3 根** —— 两个 Kekulé 式都合法。
    #[test]
    fn benzene_gets_three_double_bonds() {
        let o = orders("c1ccccc1");
        assert_eq!(o.iter().filter(|&&x| x == BondOrder::Double).count(), 3);
        assert_eq!(o.iter().filter(|&&x| x == BondOrder::Single).count(), 3);
    }

    #[test]
    fn naphthalene_gets_five_double_bonds() {
        let o = orders("c1ccc2ccccc2c1");
        assert_eq!(o.iter().filter(|&&x| x == BondOrder::Double).count(), 5);
    }

    /// 吡咯 `[nH]`:kekulize 后显式氢应转回隐式,
    /// 否则那个氢会一直挂着,影响后续价键与输出。
    #[test]
    fn pyrrole_explicit_h_becomes_implicit() {
        let m = pipeline("c1cc[nH]c1").unwrap();
        let n = m
            .atoms()
            .iter()
            .find(|a| a.atomic_num == 7)
            .expect("应有氮");
        assert_eq!(n.num_explicit_hs, 0, "显式氢应清零");
        assert!(
            !n.flags.contains(AtomFlags::NO_IMPLICIT),
            "应恢复推断隐式氢"
        );
        assert_eq!(n.num_implicit_hs, 1, "氢应转为隐式");
    }

    #[test]
    fn non_aromatic_molecule_is_untouched() {
        let before = smiles::parse("CC(=O)O").unwrap();
        let after = pipeline("CC(=O)O").unwrap();
        assert_eq!(
            before.bonds().iter().map(|b| b.order).collect::<Vec<_>>(),
            after.bonds().iter().map(|b| b.order).collect::<Vec<_>>()
        );
    }

    /// 芳香体系里的通配原子需要"可选顶点"匹配,
    /// omgkit 未实现 —— 必须**显式报错**,而不是静默给出可能错误的结果。
    ///
    /// 注意:**SMILES 写不出芳香通配原子**(`*` 没有小写形式,`[*]` 与芳香
    /// 邻居之间的键也只会是单键),所以这条路径从 SMILES 进不来 ——
    /// 语料里 0 条触发正是因为如此。守卫是给将来的 MOL/SDF 等入口留的,
    /// 这里直接构造该状态来验证它确实会开火。
    #[test]
    fn aromatic_dummy_is_rejected_loudly() {
        let mut m = smiles::parse("c1ccccc1").unwrap();
        clean_up(&mut m);
        update_property_cache(&mut m).unwrap();
        let _ = perceive_rings(&mut m);
        // 把环上一个碳换成通配原子,保留其芳香标志
        m.atom_mut(0).unwrap().atomic_num = 0;

        assert!(
            matches!(
                kekulize(&mut m),
                Err(KekulizeError::UnsupportedAromaticDummy { atom: 0 })
            ),
            "芳香通配原子必须被显式拒绝"
        );
    }

    // -- 芳香原子的氢数:SMILES 通常不写,靠推断 --

    /// **芳香声称的合法性只有 kekulize 能验证。**
    ///
    /// `c1cncc1` 是五元环 + 一个不带 H 的氮。第 3 步会顺利通过 ——
    /// 氮按 2 根芳香键(1.5+1.5=3)算,恰好等于默认价 3,推得 0 个隐式氢。
    /// 但这样只有 5 个 π 电子,不成芳香;直到 kekulize 找不到完美匹配才暴露。
    ///
    /// 这条是"kekulize 不能完全惰性化"的直接证据:它同时承担**验证**职责。
    #[test]
    fn invalid_aromatic_claim_is_caught_only_by_kekulize() {
        let mut m = smiles::parse("c1cncc1").unwrap();
        clean_up(&mut m);
        // 第 3 步顺利通过,还给氮推出了 0 个隐式氢
        let v = update_property_cache(&mut m).expect("第 3 步不该失败");
        let n_idx = m.atoms().iter().position(|a| a.atomic_num == 7).unwrap();
        assert_eq!(v.implicit_hs[n_idx], 0, "无 H 的芳香氮应推得 0 个隐式氢");

        let _ = perceive_rings(&mut m);
        assert!(
            matches!(kekulize(&mut m), Err(KekulizeError::CannotKekulize { .. })),
            "非法芳香声称必须在 kekulize 处被拒"
        );
    }

    /// 芳香杂原子的氢数由推断得出,写法上并不对称:
    /// 吡咯**必须**写 `[nH]`,而呋喃/噻吩不需要 —— O/S 有两对孤对,
    /// 拿一对去成环仍余一对,不必补氢。
    #[test]
    fn aromatic_heteroatom_hydrogens_are_inferred() {
        for (smi, sym, want_h) in [
            ("c1cc[nH]c1", 7u8, 1u8), // 吡咯:显式写的 H
            ("c1ccoc1", 8, 0),        // 呋喃:不需要 H
            ("c1ccsc1", 16, 0),       // 噻吩:不需要 H
            ("c1ccncc1", 7, 0),       // 吡啶:六元环氮不需要 H
        ] {
            let mut m = smiles::parse(smi).unwrap();
            clean_up(&mut m);
            let v = update_property_cache(&mut m).unwrap_or_else(|e| panic!("{smi}: {e}"));
            let i = m.atoms().iter().position(|a| a.atomic_num == sym).unwrap();
            let total = v.implicit_hs[i] + m.atoms()[i].num_explicit_hs;
            assert_eq!(total, want_h, "{smi}: 杂原子总氢数不对");
            let _ = perceive_rings(&mut m);
            assert!(kekulize(&mut m).is_ok(), "{smi}: 应能 kekulize");
        }
    }

    /// 融合碳有 3 根芳香键(1.5×3 = 4.5),超过碳的默认价 4。
    /// `calculateExplicitValence` 的芳香分支会把它"吸附"回 4,从而推出 0 个氢。
    /// 没有那段逻辑,融合碳会被误判成需要补氢。
    #[test]
    fn fusion_carbons_get_no_hydrogen() {
        let mut m = smiles::parse("c1ccc2ccccc2c1").unwrap();
        clean_up(&mut m);
        let v = update_property_cache(&mut m).unwrap();
        let with_h = v.implicit_hs.iter().filter(|&&h| h == 1).count();
        let without_h = v.implicit_hs.iter().filter(|&&h| h == 0).count();
        assert_eq!((with_h, without_h), (8, 2), "萘应是 8 个 CH + 2 个融合碳");
    }

    /// 通配原子在芳香环里、但自身不芳香 —— 同样必须走"未实现"分支。
    ///
    /// `c1cc[*]cc1`:存在合法解(把 `*` 也算作可接双键的候选),
    /// 当成普通原子硬算会得出"无解"这个**错误结论**。
    #[test]
    fn non_aromatic_dummy_in_aromatic_ring_is_also_rejected() {
        let mut m = smiles::parse("c1cc[*]cc1").unwrap();
        clean_up(&mut m);
        update_property_cache(&mut m).unwrap();
        let _ = perceive_rings(&mut m);
        assert!(
            matches!(
                kekulize(&mut m),
                Err(KekulizeError::UnsupportedAromaticDummy { .. })
            ),
            "应报'未实现',而不是错误的'无解'"
        );
    }
}
