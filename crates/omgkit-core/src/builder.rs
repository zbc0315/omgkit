//! 可变的单分子构造器。
//!
//! 架构上刻意与 [`MolBatch`](crate::MolBatch) 分开:编辑操作(解析建图、
//! 反应产物构建)需要随意增删原子和键,而列式布局要为此付出高昂代价。
//! 于是分成两段 —— `MolBuilder` 负责建,`freeze` 之后进入不可变的列式批,
//! 所有算法都跑在后者上。
//!
//! # 邻接索引
//!
//! `MolBuilder` 自带一份**始终有效**的邻接索引,[`neighbors`](MolBuilder::neighbors)
//! 与 [`degree`](MolBuilder::degree) 都是 O(度数)。没有它的话,"取某原子的
//! 邻居"只能扫全部键,而化学算法几乎每一步都在做这件事,整体就退化成
//! O(原子数 × 键数)。
//!
//! ## 为什么是半边链表,不是 CSR
//!
//! [`MolBatch`](crate::MolBatch) 用 CSR,因为它不可变、只被扫描。`MolBuilder` 是**增量构建**
//! 的:CSR 每加一条键都要重排整个邻接数组,解析一个 E 条键的分子就变成
//! O(E·(V+E))。半边链表加边是 O(1),代价是遍历时跳指针 —— 在编辑期的
//! 规模下无所谓,真正吃吞吐的批量算法跑在 `MolBatch` 上。
//!
//! 这个分工是架构层面的:**`MolBuilder` 为编辑优化,`MolBatch` 为扫描优化。**
//!
//! ## 索引不会失效
//!
//! 索引在 [`add_bond_data`](MolBuilder::add_bond_data) 里同步维护,不是缓存,
//! 没有脏标记,也不需要谁记得去刷新。为此 [`bond_mut`](MolBuilder::bond_mut)
//! 返回的 [`BondMut`] **不暴露端点** —— 改端点就是改拓扑,只能走建边接口。
//!
//! ```
//! use omgkit_core::{MolBuilder, BondOrder};
//!
//! // 乙醇 CCO
//! let mut b = MolBuilder::new();
//! let c0 = b.add_atom(6);
//! let c1 = b.add_atom(6);
//! let o  = b.add_atom(8);
//! b.add_bond(c0, c1, BondOrder::Single).unwrap();
//! b.add_bond(c1, o,  BondOrder::Single).unwrap();
//! assert_eq!(b.num_atoms(), 3);
//! assert_eq!(b.num_bonds(), 2);
//! ```

use crate::error::{Error, Result};
use crate::types::{
    AtomFlags, BondDirection, BondFlags, BondOrder, BondStereo, ChiralTag, Hybridization,
};

/// 单个原子的可变数据。
///
/// 字段与 [`MolBatch`](crate::MolBatch) 的列一一对应 —— 增加字段时两处必须同步。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AtomData {
    /// 原子序数。0 表示 SMILES 通配原子 `*`。
    pub atomic_num: u8,
    /// 形式电荷
    pub formal_charge: i8,
    /// 同位素质量数。0 表示未指定(即天然丰度)。
    pub isotope: u16,
    /// 方括号中显式书写的氢数。仅当 [`AtomFlags::NO_IMPLICIT`] 置位时有意义。
    pub num_explicit_hs: u8,
    /// 隐式氢数。由 L2 的价键计算填充,解析阶段恒为 0。
    pub num_implicit_hs: u8,
    /// 自由基电子数。由 L2 第 6 步 `FINDRADICALS` 填充,在那之前恒为 0。
    ///
    /// 隐式氢推断会读它,所以两者的先后顺序有实际后果 ——
    /// 净化管线里第 6 步排在第 3 步**之后**,故第 3 步看到的必然是 0。
    pub num_radical_electrons: u8,
    /// 反应原子映射号(SMILES 中的 `:n`)。0 表示无映射。
    pub atom_map: u16,
    /// 立体标记的几何类别
    pub chiral_tag: ChiralTag,
    /// 立体标记的类内排列序号,与 [`chiral_tag`](Self::chiral_tag) 配套。
    ///
    /// 0 表示未指定。四面体的两种排列由 `chiral_tag` 自身表达,此处恒为 0 ——
    /// 两处都记就有了两个可以互相矛盾的真相来源。
    ///
    /// # 相对**邻居的存储顺序**
    ///
    /// 序号的含义依赖"配体按什么顺序排列",而解析会重排邻居(环闭合键统一
    /// 追加到末尾)—— 不说清楚"相对什么顺序",这个字段就没有意义。
    ///
    /// | 类别 | 多面体有几个顶点 | 本字段相对什么顺序 |
    /// |---|---|---|
    /// | [`ChiralTag::SquarePlanar`] | 4 | 邻居的**存储顺序**(解析时已归一) |
    /// | [`ChiralTag::TrigonalBipyramidal`] | 5 | 同上 |
    /// | [`ChiralTag::Octahedral`] | 6 | 同上 |
    /// | 以上三类但**顶点比邻居多两个及以上** | — | 书写时的字面值;写出时整个丢掉 |
    /// | [`ChiralTag::Allene`] (`@AL`) | 4(来自累积双键**两端**) | 四个配体的**存储顺序**(解析时已归一) |
    ///
    /// 丙二烯那一行的"配体"不是这个原子自己的邻居 —— 中心只有两个邻居,
    /// 四个配体是两端端原子上的取代基。取不到四个(标记没落在丙二烯中心上、
    /// 某一端两个配体相同)时序号归 0,写出侧整个丢掉。
    ///
    /// 顶点比邻居**多一个**时照样归一:方括号里的氢、或者一个空的配位位置,
    /// 也占一个顶点而不在邻居序列里,补一个占位的进去即可 —— 它排在存储序的
    /// **最前**,书写序里则落在"自身位置"。补法在 `omgkit-io` 那一侧
    /// (`smiles::coordination_ligands`),解析与写出共用同一份。
    ///
    /// **多两个及以上**就不归一了:那几个顶点在 SMILES 里全落在同一处、彼此
    /// 分不开,换算出来的序号不唯一。那时写出侧**整个丢掉**这个标记:丢掉是
    /// 老实的,瞎写一个序号是撒谎。
    ///
    /// 换算见 [`crate::polyhedron::renumber`];那三张表是从 RDKit 2025.09.2
    /// 穷举量出来的(72 / 2400 / 21600 种写法,零反例),不是照着规范条文写的。
    ///
    pub stereo_perm: u8,
    /// 杂化状态。由净化第 9 步填充,在那之前为 `Unspecified`。
    pub hybridization: Hybridization,
    /// 标志位
    pub flags: AtomFlags,
}

impl AtomData {
    /// 构造一个只指定了元素的原子,其余字段取默认值。
    #[must_use]
    pub fn new(atomic_num: u8) -> Self {
        Self {
            atomic_num,
            formal_charge: 0,
            isotope: 0,
            num_explicit_hs: 0,
            num_implicit_hs: 0,
            num_radical_electrons: 0,
            atom_map: 0,
            chiral_tag: ChiralTag::Unspecified,
            stereo_perm: 0,
            hybridization: Hybridization::Unspecified,
            flags: AtomFlags::NONE,
        }
    }
}

impl Default for AtomData {
    fn default() -> Self {
        Self::new(0)
    }
}

/// 单条键的可变数据。端点为**分子内局部**原子下标。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BondData {
    /// 起点原子(局部下标)
    pub begin: u32,
    /// 终点原子(局部下标)
    pub end: u32,
    /// 键级
    pub order: BondOrder,
    /// SMILES 方向键 `/` `\`
    pub direction: BondDirection,
    /// 双键立体。由 omgkit-io 的 stereo::perceive_bond_stereo 从 direction 感知,
    /// 与 [`stereo_atoms`](Self::stereo_atoms) 配套。
    pub stereo: BondStereo,
    /// 顺/反的**参照原子**:`[begin 侧一个, end 侧一个]`。
    ///
    /// 只有 `stereo` 不是 [`BondStereo::None`] 时有意义,那时两个下标都合法。
    /// 无参照时是 [`BondData::NO_STEREO_ATOM`]。
    ///
    /// # 为什么必须存
    ///
    /// "顺"与"反"离开参照就没有意义 —— 一根双键两端各有两个取代基,说
    /// "同侧"总得回答"谁和谁同侧"。
    ///
    /// 这也是它与 [`BondData::direction`] 的分工:方向是**写法**,依附于
    /// 某根单键,那根键被删掉信息就没了;顺反是双键**自己**的属性,只要
    /// 两个参照原子还在就一直成立。图编辑之后要重新写出方向键,靠的是这一对。
    pub stereo_atoms: [u32; 2],
    /// 标志位
    pub flags: BondFlags,
}

impl BondData {
    /// [`BondData::stereo_atoms`] 里表示"没有参照原子"的值。
    pub const NO_STEREO_ATOM: u32 = u32::MAX;

    /// 构造一条指定端点与键级的键。
    #[must_use]
    pub fn new(begin: u32, end: u32, order: BondOrder) -> Self {
        Self {
            begin,
            end,
            order,
            direction: BondDirection::None,
            stereo: BondStereo::None,
            stereo_atoms: [Self::NO_STEREO_ATOM; 2],
            flags: BondFlags::NONE,
        }
    }

    /// 本键对端点 `atom` 的**价贡献**。
    ///
    /// 配位键的贡献是**不对称**的:对起点(给体)算 0,对终点(受体)算 1。
    /// 这与直觉相反,也与 [`BondOrder::as_double`] 不同 —— 后者对配位键
    /// 一律算 1。搞混会让隐式氢推断在有机金属分子上系统性出错。
    ///
    /// `atom` 不是本键端点时返回 0。
    #[must_use]
    pub fn valence_contribution_to(&self, atom: u32) -> f32 {
        if atom != self.begin && atom != self.end {
            return 0.0;
        }
        if self.order == BondOrder::Dative && atom != self.end {
            return 0.0; // 给体不计
        }
        self.order.as_double()
    }

    /// 给定一端,返回另一端。若 `from` 不是本键端点则返回 `None`。
    #[must_use]
    pub fn other_end(&self, from: u32) -> Option<u32> {
        if from == self.begin {
            Some(self.end)
        } else if from == self.end {
            Some(self.begin)
        } else {
            None
        }
    }
}

/// 半边编号的空标记。
const NO_HALF: u32 = u32::MAX;

/// 可变的单分子构造器。
///
/// 自带邻接索引;见[模块文档](self)。
#[derive(Debug, Clone, Default)]
pub struct MolBuilder {
    atoms: Vec<AtomData>,
    bonds: Vec<BondData>,
    name: Option<String>,

    // -- 邻接索引:每原子一条半边链表 --
    //
    // 键 `bi` 拆成两条半边:`2*bi` 挂在 `begin` 上,`2*bi+1` 挂在 `end` 上。
    // 由半边编号 `h` 可反推:键号 `h >> 1`,`h & 1 == 0` 表示自己是 begin 侧,
    // 于是邻居就是另一端。这样不必为邻居再存一份原子号。
    /// 每原子的链首半边,`NO_HALF` 表示孤立原子
    first_half: Vec<u32>,
    /// 每原子的链尾半边,用于 O(1) 尾插 —— 尾插保证遍历顺序 = 键的插入顺序
    last_half: Vec<u32>,
    /// 每条半边的后继,`NO_HALF` 表示链尾。下标即半边编号
    next_half: Vec<u32>,
    /// 每原子的度(不含隐式氢)
    degree: Vec<u32>,
}

impl MolBuilder {
    /// 空分子。
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// 预分配容量的空分子。解析器已知大致规模时用这个避免反复扩容。
    #[must_use]
    pub fn with_capacity(n_atoms: usize, n_bonds: usize) -> Self {
        Self {
            atoms: Vec::with_capacity(n_atoms),
            bonds: Vec::with_capacity(n_bonds),
            name: None,
            first_half: Vec::with_capacity(n_atoms),
            last_half: Vec::with_capacity(n_atoms),
            next_half: Vec::with_capacity(n_bonds * 2),
            degree: Vec::with_capacity(n_atoms),
        }
    }

    /// 追加一个只指定元素的原子,返回其局部下标。
    pub fn add_atom(&mut self, atomic_num: u8) -> u32 {
        self.add_atom_data(AtomData::new(atomic_num))
    }

    /// 追加一个完整指定的原子,返回其局部下标。
    pub fn add_atom_data(&mut self, atom: AtomData) -> u32 {
        let idx = self.atoms.len() as u32;
        self.atoms.push(atom);
        self.first_half.push(NO_HALF);
        self.last_half.push(NO_HALF);
        self.degree.push(0);
        idx
    }

    /// 追加一条键,返回其局部下标。
    ///
    /// # Errors
    /// 端点越界、或两端相同(自环)时返回错误。分子图中不存在自环,
    /// 让它静默通过会在 CSR 构建时产生难以追查的度数异常。
    pub fn add_bond(&mut self, begin: u32, end: u32, order: BondOrder) -> Result<u32> {
        self.add_bond_data(BondData::new(begin, end, order))
    }

    /// 追加一条完整指定的键,返回其局部下标。
    ///
    /// 同时把两条半边挂进端点的邻接链表,O(1)。
    ///
    /// # Errors
    /// 同 [`add_bond`](Self::add_bond)。
    ///
    /// # Panics
    /// 键数超过 `u32::MAX / 2` 时 —— 半边编号会与哨兵冲突。这需要 32 GB
    /// 仅存键表,实际不可达,但宁可响亮地停下也不要静默算错。
    pub fn add_bond_data(&mut self, bond: BondData) -> Result<u32> {
        let n = self.atoms.len() as u32;
        if bond.begin >= n || bond.end >= n {
            return Err(Error::AtomIndexOutOfRange {
                index: bond.begin.max(bond.end),
                num_atoms: n,
            });
        }
        if bond.begin == bond.end {
            return Err(Error::SelfLoop { atom: bond.begin });
        }
        assert!(
            self.bonds.len() < (u32::MAX / 2) as usize,
            "键数超出半边编号能表示的范围"
        );
        let idx = self.bonds.len() as u32;
        self.bonds.push(bond);
        self.link_half(bond.begin, idx * 2);
        self.link_half(bond.end, idx * 2 + 1);
        Ok(idx)
    }

    /// 交换一条键的两个端点。
    ///
    /// 端点顺序对**配位键**有语义(`begin` 是给电子的一端),这个接口就是为它
    /// 准备的。其余键型的端点顺序只是书写痕迹,交换它没有意义。
    ///
    /// # 为什么是整体重建索引
    ///
    /// 半边链表是单向的:把两条半边在两个原子的链表之间对调,要各自找到前驱
    /// 再改指针,还要顾及链首链尾 —— 几十行且容易出错。而整体重建只有几行,
    /// 显然正确,并且保持"邻居顺序 = 键的插入顺序"这条不变量(重建按键号
    /// 递增走,每个原子的相对顺序不变)。
    ///
    /// 代价是 O(原子数 + 键数)。这个接口的调用场景是净化第 2 步,
    /// 实测 8839 条语料里只触发 2 次 —— 为它做精细的指针手术不划算。
    ///
    /// # Errors
    /// 键下标越界时返回 [`Error::BondIndexOutOfRange`]。
    pub fn swap_bond_ends(&mut self, bond: u32) -> Result<()> {
        let num_bonds = self.bonds.len() as u32;
        let b = self
            .bonds
            .get_mut(bond as usize)
            .ok_or(Error::BondIndexOutOfRange {
                index: bond,
                num_bonds,
            })?;
        std::mem::swap(&mut b.begin, &mut b.end);
        self.rebuild_index();
        Ok(())
    }

    /// 按当前的键表重建邻接索引。
    fn rebuild_index(&mut self) {
        let n = self.atoms.len();
        self.first_half.clear();
        self.first_half.resize(n, NO_HALF);
        self.last_half.clear();
        self.last_half.resize(n, NO_HALF);
        self.degree.clear();
        self.degree.resize(n, 0);
        self.next_half.clear();

        for i in 0..self.bonds.len() {
            let (begin, end) = (self.bonds[i].begin, self.bonds[i].end);
            self.link_half(begin, (i * 2) as u32);
            self.link_half(end, (i * 2 + 1) as u32);
        }
    }

    /// 把半边 `half` 尾插到 `atom` 的链表上。
    ///
    /// 必须按半边编号递增调用 —— `next_half` 的下标就是半边编号。
    fn link_half(&mut self, atom: u32, half: u32) {
        debug_assert_eq!(
            self.next_half.len() as u32,
            half,
            "半边必须按编号顺序追加,否则 next_half 的下标语义就断了"
        );
        self.next_half.push(NO_HALF);

        let a = atom as usize;
        let tail = self.last_half[a];
        if tail == NO_HALF {
            self.first_half[a] = half;
        } else {
            self.next_half[tail as usize] = half;
        }
        self.last_half[a] = half;
        self.degree[a] += 1;
    }

    /// 遍历某原子的邻居,产出 `(邻居原子下标, 键下标)`。O(度数)。
    ///
    /// 顺序即**键的插入顺序**,与 [`MolView::neighbors`](crate::MolView::neighbors)
    /// 一致 —— 手性语义依赖于这个顺序,两处必须给出同样的序列。
    ///
    /// 原子下标越界时产出空序列。
    #[must_use]
    pub fn neighbors(&self, atom: u32) -> Neighbors<'_> {
        let head = self
            .first_half
            .get(atom as usize)
            .copied()
            .unwrap_or(NO_HALF);
        Neighbors {
            mol: self,
            half: head,
        }
    }

    /// 原子的度(不含隐式氢)。O(1);越界返回 0。
    #[must_use]
    pub fn degree(&self, atom: u32) -> usize {
        self.degree.get(atom as usize).copied().unwrap_or(0) as usize
    }

    /// 连接 `a` 与 `b` 的键下标。O(min 度数);不相邻时返回 `None`。
    #[must_use]
    pub fn bond_between(&self, a: u32, b: u32) -> Option<u32> {
        // 从度数小的一端出发
        let from = if self.degree(a) <= self.degree(b) {
            a
        } else {
            b
        };
        let to = if from == a { b } else { a };
        self.neighbors(from)
            .find(|&(nbr, _)| nbr == to)
            .map(|(_, bi)| bi)
    }

    /// 原子数
    #[must_use]
    pub fn num_atoms(&self) -> usize {
        self.atoms.len()
    }

    /// 键数
    #[must_use]
    pub fn num_bonds(&self) -> usize {
        self.bonds.len()
    }

    /// 只读访问全部原子
    #[must_use]
    pub fn atoms(&self) -> &[AtomData] {
        &self.atoms
    }

    /// 只读访问全部键
    #[must_use]
    pub fn bonds(&self) -> &[BondData] {
        &self.bonds
    }

    /// 可变访问单个原子
    pub fn atom_mut(&mut self, idx: u32) -> Option<&mut AtomData> {
        self.atoms.get_mut(idx as usize)
    }

    /// 可变访问单条键的**属性**。
    ///
    /// 返回的 [`BondMut`] 刻意不暴露端点 —— 端点即拓扑,而邻接索引是随建边
    /// 增量维护的。若允许就地改端点,索引会在无人察觉的情况下失效,而这类
    /// bug 只会在很久以后以"某个原子少了个邻居"的形式冒出来。
    pub fn bond_mut(&mut self, idx: u32) -> Option<BondMut<'_>> {
        self.bonds
            .get_mut(idx as usize)
            .map(|bond| BondMut { bond })
    }

    /// 分子名(SDF 的标题行、SMILES 后跟的名字)
    #[must_use]
    pub fn name(&self) -> Option<&str> {
        self.name.as_deref()
    }

    /// 设置分子名
    pub fn set_name(&mut self, name: impl Into<String>) {
        self.name = Some(name.into());
    }

    /// 重新计算邻接索引并与当前索引比对。仅供测试使用。
    ///
    /// 索引是增量维护的,任何维护逻辑的错误都会静默地表现为"图连错了"。
    /// 这个自检把它变成显式失败。
    #[doc(hidden)]
    #[must_use]
    pub fn adjacency_index_is_consistent(&self) -> bool {
        for a in 0..self.atoms.len() as u32 {
            let expected: Vec<(u32, u32)> = self
                .bonds
                .iter()
                .enumerate()
                .filter_map(|(bi, b)| b.other_end(a).map(|o| (o, bi as u32)))
                .collect();
            let actual: Vec<(u32, u32)> = self.neighbors(a).collect();
            if expected != actual || self.degree(a) != expected.len() {
                return false;
            }
        }
        true
    }
}

/// [`MolBuilder::neighbors`] 的迭代器,沿半边链表前进。
#[derive(Debug, Clone)]
pub struct Neighbors<'a> {
    mol: &'a MolBuilder,
    half: u32,
}

impl Iterator for Neighbors<'_> {
    /// `(邻居原子下标, 键下标)`
    type Item = (u32, u32);

    fn next(&mut self) -> Option<Self::Item> {
        if self.half == NO_HALF {
            return None;
        }
        let h = self.half;
        self.half = self.mol.next_half[h as usize];

        let bi = h >> 1;
        let bond = self.mol.bonds[bi as usize];
        // h 为偶 ⇒ 挂在 begin 上 ⇒ 邻居是 end
        let nbr = if h & 1 == 0 { bond.end } else { bond.begin };
        Some((nbr, bi))
    }
}

/// 键属性的可变句柄,由 [`MolBuilder::bond_mut`] 取得。
///
/// **不提供修改端点的途径**,理由见 [`MolBuilder::bond_mut`]。
#[derive(Debug)]
pub struct BondMut<'a> {
    bond: &'a mut BondData,
}

impl BondMut<'_> {
    /// 读回当前键的完整数据
    #[must_use]
    pub fn get(&self) -> BondData {
        *self.bond
    }

    /// 设置键级
    pub fn set_order(&mut self, order: BondOrder) {
        self.bond.order = order;
    }

    /// 设置方向键标记(`/` `\`)
    pub fn set_direction(&mut self, direction: BondDirection) {
        self.bond.direction = direction;
    }

    /// 设置双键立体
    pub fn set_stereo(&mut self, stereo: BondStereo) {
        self.bond.stereo = stereo;
    }

    /// 设置顺反的参照原子。
    ///
    /// 与 [`set_stereo`](Self::set_stereo) 配套 —— 顺反离开参照没有意义,
    /// 两者要一起写。
    pub fn set_stereo_atoms(&mut self, atoms: [u32; 2]) {
        self.bond.stereo_atoms = atoms;
    }

    /// 可变访问标志位
    pub fn flags_mut(&mut self) -> &mut BondFlags {
        &mut self.bond.flags
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn build_ethanol() {
        let mut b = MolBuilder::new();
        let c0 = b.add_atom(6);
        let c1 = b.add_atom(6);
        let o = b.add_atom(8);
        b.add_bond(c0, c1, BondOrder::Single).unwrap();
        b.add_bond(c1, o, BondOrder::Single).unwrap();

        assert_eq!(b.num_atoms(), 3);
        assert_eq!(b.num_bonds(), 2);
        assert_eq!(b.atoms()[2].atomic_num, 8);
    }

    #[test]
    fn rejects_out_of_range_endpoint() {
        let mut b = MolBuilder::new();
        b.add_atom(6);
        let err = b.add_bond(0, 5, BondOrder::Single).unwrap_err();
        assert!(matches!(
            err,
            Error::AtomIndexOutOfRange {
                index: 5,
                num_atoms: 1
            }
        ));
    }

    #[test]
    fn rejects_self_loop() {
        let mut b = MolBuilder::new();
        b.add_atom(6);
        let err = b.add_bond(0, 0, BondOrder::Single).unwrap_err();
        assert!(matches!(err, Error::SelfLoop { atom: 0 }));
    }

    #[test]
    fn bond_other_end() {
        let bond = BondData::new(3, 7, BondOrder::Double);
        assert_eq!(bond.other_end(3), Some(7));
        assert_eq!(bond.other_end(7), Some(3));
        assert_eq!(bond.other_end(5), None);
    }
}

#[cfg(test)]
mod adjacency_tests {
    use super::*;

    /// 异丁烷 CC(C)C:中心碳连着三个甲基
    fn isobutane() -> MolBuilder {
        let mut m = MolBuilder::new();
        for _ in 0..4 {
            m.add_atom(6);
        }
        m.add_bond(0, 1, BondOrder::Single).unwrap();
        m.add_bond(1, 2, BondOrder::Single).unwrap();
        m.add_bond(1, 3, BondOrder::Single).unwrap();
        m
    }

    #[test]
    fn neighbors_and_degree() {
        let m = isobutane();
        assert_eq!(
            m.neighbors(1).collect::<Vec<_>>(),
            vec![(0, 0), (2, 1), (3, 2)]
        );
        assert_eq!(m.neighbors(0).collect::<Vec<_>>(), vec![(1, 0)]);
        assert_eq!(m.degree(1), 3);
        assert_eq!(m.degree(0), 1);
    }

    #[test]
    fn isolated_atom_has_no_neighbors() {
        let mut m = MolBuilder::new();
        m.add_atom(10); // Ne
        assert_eq!(m.neighbors(0).count(), 0);
        assert_eq!(m.degree(0), 0);
    }

    #[test]
    fn out_of_range_atom_is_empty_not_panic() {
        let m = isobutane();
        assert_eq!(m.neighbors(99).count(), 0);
        assert_eq!(m.degree(99), 0);
        assert_eq!(m.bond_between(99, 0), None);
    }

    /// 邻居顺序必须是**键的插入顺序** —— 手性判定依赖于此。
    /// 这条与 `MolView::neighbors` 的同名不变量是一对,两边不能各说各话。
    #[test]
    fn neighbor_order_is_bond_insertion_order() {
        let mut m = MolBuilder::new();
        for _ in 0..4 {
            m.add_atom(6);
        }
        // 故意让中心原子在键里时而作 begin、时而作 end
        m.add_bond(3, 0, BondOrder::Single).unwrap();
        m.add_bond(0, 1, BondOrder::Single).unwrap();
        m.add_bond(2, 0, BondOrder::Single).unwrap();

        assert_eq!(
            m.neighbors(0).map(|(a, _)| a).collect::<Vec<_>>(),
            vec![3, 1, 2],
            "无论中心原子在哪一端,顺序都应是键的插入顺序"
        );
        assert_eq!(
            m.neighbors(0).map(|(_, b)| b).collect::<Vec<_>>(),
            vec![0, 1, 2]
        );
    }

    #[test]
    fn bond_between_finds_edges_from_either_side() {
        let m = isobutane();
        assert_eq!(m.bond_between(1, 0), Some(0));
        assert_eq!(m.bond_between(0, 1), Some(0));
        assert_eq!(m.bond_between(1, 3), Some(2));
        assert_eq!(m.bond_between(0, 2), None, "0 与 2 不相邻");
    }

    /// 索引是增量维护的,拒边之后必须原样不动 —— 否则一次失败的
    /// `add_bond` 会悄悄污染整张图。
    #[test]
    fn rejected_bond_leaves_index_untouched() {
        let mut m = isobutane();
        assert!(m.add_bond(1, 1, BondOrder::Single).is_err());
        assert!(m.add_bond(0, 99, BondOrder::Single).is_err());
        assert_eq!(m.degree(1), 3);
        assert_eq!(m.num_bonds(), 3);
        assert!(m.adjacency_index_is_consistent());
    }

    /// 建图过程中的**每一步**索引都必须自洽,不只是最后一步。
    #[test]
    fn index_stays_consistent_through_incremental_build() {
        let mut m = MolBuilder::new();
        assert!(m.adjacency_index_is_consistent());
        for i in 0..12u32 {
            m.add_atom(6);
            assert!(m.adjacency_index_is_consistent(), "加原子 {i} 后失配");
            if i > 0 {
                m.add_bond(i - 1, i, BondOrder::Single).unwrap();
                assert!(m.adjacency_index_is_consistent(), "加键 {i} 后失配");
            }
        }
        // 再补几条成环的键,制造非线性拓扑
        m.add_bond(0, 11, BondOrder::Single).unwrap();
        m.add_bond(3, 8, BondOrder::Single).unwrap();
        assert!(m.adjacency_index_is_consistent());
    }

    /// 克隆出来的分子必须带着一份同样有效的索引。
    #[test]
    fn clone_carries_a_valid_index() {
        let m = isobutane().clone();
        assert!(m.adjacency_index_is_consistent());
        assert_eq!(m.degree(1), 3);
    }

    /// 改键级、改标志都不该动到拓扑。
    #[test]
    fn property_edits_do_not_disturb_topology() {
        let mut m = isobutane();
        let mut b = m.bond_mut(1).unwrap();
        b.set_order(BondOrder::Double);
        b.flags_mut().insert(BondFlags::AROMATIC);
        b.set_direction(BondDirection::UpRight);
        assert!(m.adjacency_index_is_consistent());
        assert_eq!(m.bonds()[1].order, BondOrder::Double);
        assert_eq!(m.bonds()[1].direction, BondDirection::UpRight);
    }
}

#[cfg(test)]
mod valence_contrib_tests {
    use super::*;

    /// 配位键的价贡献不对称:起点(给体)算 0,终点(受体)算 1。
    #[test]
    fn dative_contribution_is_asymmetric() {
        let d = BondData::new(3, 7, BondOrder::Dative);
        assert_eq!(d.valence_contribution_to(3), 0.0, "给体不计价");
        assert_eq!(d.valence_contribution_to(7), 1.0, "受体计 1");
        assert_eq!(d.valence_contribution_to(9), 0.0, "非端点");
    }

    #[test]
    fn normal_bonds_are_symmetric() {
        for (order, v) in [
            (BondOrder::Single, 1.0),
            (BondOrder::Double, 2.0),
            (BondOrder::Triple, 3.0),
            (BondOrder::Aromatic, 1.5),
        ] {
            let b = BondData::new(1, 2, order);
            assert_eq!(b.valence_contribution_to(1), v);
            assert_eq!(b.valence_contribution_to(2), v);
        }
    }
}
