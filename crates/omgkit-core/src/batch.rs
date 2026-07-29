//! 列式分子批 —— omgkit 的核心数据结构。
//!
//! # 为什么是列式
//!
//! 常见的分子图实现用邻接表 + 堆分配的原子对象
//! (`Code/GraphMol/ROMol.h:55`),顶点属性是**堆上的 `Atom` 指针** ——
//! 遍历一个 40 原子分子的原子属性就是 40 次随机访存。
//!
//! [`MolBatch`] 把同一属性的所有原子放进一个连续数组(SoA),邻接用 CSR。
//! 这一个决策同时买到:
//!
//! 1. **顺序访存** —— 属性遍历是一次连续扫描,不是逐原子随机跳转
//! 2. **零拷贝** —— 列可直接暴露为 numpy / Arrow buffer
//! 3. **廉价并行** —— 按分子偏移切分,无共享状态
//!
//! 这三条都是**布局的直接后果**,不需要额外测量就成立。至于这套布局能带来
//! 多少端到端加速,那是另一回事,要拿基准说话,不写在这里。
//!
//! # 索引约定
//!
//! 批内部一律使用**全局**下标(跨分子连续)。这是 CSR 的自然形式。
//! 面向用户的**局部**下标(分子内 0 起)由 [`MolView`] 换算,
//! 它持有 `&MolBatch` 与分子号,不复制任何数据。

use crate::builder::{BondData, MolBuilder};
use crate::error::{Error, Result};
use crate::types::{
    AtomFlags, BondDirection, BondFlags, BondOrder, BondStereo, ChiralTag, Hybridization,
};
use crate::view::MolView;

/// 一批分子的不可变列式存储。
///
/// 用 [`MolBatchBuilder`] 构造。
#[derive(Debug, Clone, Default)]
pub struct MolBatch {
    // ---- 分子边界(长度均为 n_mols + 1)----
    pub(crate) mol_atom_offset: Vec<u32>,
    pub(crate) mol_bond_offset: Vec<u32>,
    pub(crate) names: Vec<Option<String>>,

    // ---- 每原子列(长度均为 n_atoms)----
    pub(crate) atomic_num: Vec<u8>,
    pub(crate) formal_charge: Vec<i8>,
    pub(crate) isotope: Vec<u16>,
    pub(crate) num_explicit_hs: Vec<u8>,
    pub(crate) num_implicit_hs: Vec<u8>,
    pub(crate) num_radical_electrons: Vec<u8>,
    pub(crate) atom_map: Vec<u16>,
    pub(crate) chiral_tag: Vec<ChiralTag>,
    pub(crate) stereo_perm: Vec<u8>,
    pub(crate) hybridization: Vec<Hybridization>,
    pub(crate) atom_flags: Vec<AtomFlags>,

    // ---- CSR 邻接 ----
    /// 长度 n_atoms + 1
    pub(crate) nbr_offset: Vec<u32>,
    /// 长度 2 * n_bonds,存全局原子下标
    pub(crate) nbr_atom: Vec<u32>,
    /// 长度 2 * n_bonds,存全局键下标
    pub(crate) nbr_bond: Vec<u32>,

    // ---- 每键列(长度均为 n_bonds;端点为全局原子下标)----
    pub(crate) bond_begin: Vec<u32>,
    pub(crate) bond_end: Vec<u32>,
    pub(crate) bond_order: Vec<BondOrder>,
    pub(crate) bond_direction: Vec<BondDirection>,
    pub(crate) bond_stereo: Vec<BondStereo>,
    /// 顺反的参照原子,与端点一样是**全局**下标;
    /// [`BondData::NO_STEREO_ATOM`] 表示无参照,那个值不做基址平移。
    pub(crate) bond_stereo_atoms: Vec<[u32; 2]>,
    pub(crate) bond_flags: Vec<BondFlags>,
}

impl MolBatch {
    /// 分子数
    #[must_use]
    pub fn num_mols(&self) -> usize {
        self.mol_atom_offset.len().saturating_sub(1)
    }

    /// 批中是否没有分子
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.num_mols() == 0
    }

    /// 全批原子总数
    #[must_use]
    pub fn num_atoms(&self) -> usize {
        self.atomic_num.len()
    }

    /// 全批键总数
    #[must_use]
    pub fn num_bonds(&self) -> usize {
        self.bond_order.len()
    }

    /// 取第 `idx` 个分子的零拷贝视图。
    #[must_use]
    pub fn mol(&self, idx: u32) -> Option<MolView<'_>> {
        if (idx as usize) < self.num_mols() {
            Some(MolView::new(self, idx))
        } else {
            None
        }
    }

    /// 取第 `idx` 个分子的视图。
    ///
    /// # Errors
    /// 下标越界时返回 [`Error::MolIndexOutOfRange`]。
    pub fn try_mol(&self, idx: u32) -> Result<MolView<'_>> {
        self.mol(idx).ok_or(Error::MolIndexOutOfRange {
            index: idx,
            num_mols: self.num_mols() as u32,
        })
    }

    /// 遍历全部分子。
    pub fn iter(&self) -> impl Iterator<Item = MolView<'_>> + '_ {
        (0..self.num_mols() as u32).map(move |i| MolView::new(self, i))
    }

    // ---- 列的原始访问 ----
    // 这些是将来零拷贝暴露给 numpy / Arrow 的接口。切片形式保证了
    // 无论上层怎么包装,底下始终是一块连续内存。

    /// 原子序数列(全局)
    #[must_use]
    pub fn atomic_nums(&self) -> &[u8] {
        &self.atomic_num
    }

    /// 形式电荷列(全局)
    #[must_use]
    pub fn formal_charges(&self) -> &[i8] {
        &self.formal_charge
    }

    /// 键级列(全局)
    #[must_use]
    pub fn bond_orders(&self) -> &[BondOrder] {
        &self.bond_order
    }

    /// CSR 邻接偏移(长度 `num_atoms() + 1`)
    #[must_use]
    pub fn nbr_offsets(&self) -> &[u32] {
        &self.nbr_offset
    }

    /// 每分子的原子起始偏移(长度 `num_mols() + 1`)
    #[must_use]
    pub fn mol_atom_offsets(&self) -> &[u32] {
        &self.mol_atom_offset
    }
}

/// [`MolBatch`] 的流式构造器。
///
/// 解析器边解析边 [`push`](Self::push),最后一次性 [`finish`](Self::finish)
/// 构建 CSR —— CSR 只做一遍计数排序,比每分子单独建图快得多。
#[derive(Debug, Clone, Default)]
pub struct MolBatchBuilder {
    batch: MolBatch,
}

impl MolBatchBuilder {
    /// 空构造器。
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// 追加一个分子,返回其在批中的下标。
    ///
    /// # Errors
    /// 累计原子/键/分子数超出 `u32` 索引上限时返回 [`Error::BatchTooLarge`]。
    pub fn push(&mut self, mol: &MolBuilder) -> Result<u32> {
        let b = &mut self.batch;

        // 首次 push 时补上起始偏移 0
        if b.mol_atom_offset.is_empty() {
            b.mol_atom_offset.push(0);
            b.mol_bond_offset.push(0);
        }

        let atom_base = b.atomic_num.len();
        let n_atoms_after = atom_base + mol.num_atoms();
        let n_bonds_after = b.bond_order.len() + mol.num_bonds();
        if n_atoms_after > u32::MAX as usize {
            return Err(Error::BatchTooLarge {
                what: "atoms",
                count: n_atoms_after,
            });
        }
        if n_bonds_after > u32::MAX as usize {
            return Err(Error::BatchTooLarge {
                what: "bonds",
                count: n_bonds_after,
            });
        }
        let mol_idx = b.names.len();
        if mol_idx >= u32::MAX as usize {
            return Err(Error::BatchTooLarge {
                what: "molecules",
                count: mol_idx + 1,
            });
        }

        for a in mol.atoms() {
            b.atomic_num.push(a.atomic_num);
            b.formal_charge.push(a.formal_charge);
            b.isotope.push(a.isotope);
            b.num_explicit_hs.push(a.num_explicit_hs);
            b.num_implicit_hs.push(a.num_implicit_hs);
            b.num_radical_electrons.push(a.num_radical_electrons);
            b.atom_map.push(a.atom_map);
            b.chiral_tag.push(a.chiral_tag);
            b.stereo_perm.push(a.stereo_perm);
            b.hybridization.push(a.hybridization);
            b.atom_flags.push(a.flags);
        }

        // 端点由局部下标平移为全局下标
        let base = atom_base as u32;
        for bd in mol.bonds() {
            b.bond_begin.push(base + bd.begin);
            b.bond_end.push(base + bd.end);
            b.bond_order.push(bd.order);
            b.bond_direction.push(bd.direction);
            b.bond_stereo.push(bd.stereo);
            // 无参照的哨兵值不能平移,平移过的哨兵就不再是哨兵了
            b.bond_stereo_atoms.push(bd.stereo_atoms.map(|a| {
                if a == BondData::NO_STEREO_ATOM {
                    a
                } else {
                    base + a
                }
            }));
            b.bond_flags.push(bd.flags);
        }

        b.names.push(mol.name().map(str::to_owned));
        b.mol_atom_offset.push(n_atoms_after as u32);
        b.mol_bond_offset.push(n_bonds_after as u32);

        Ok(mol_idx as u32)
    }

    /// 构建 CSR 并冻结为 [`MolBatch`]。
    #[must_use]
    pub fn finish(mut self) -> MolBatch {
        if self.batch.mol_atom_offset.is_empty() {
            self.batch.mol_atom_offset.push(0);
            self.batch.mol_bond_offset.push(0);
        }
        build_csr(&mut self.batch);
        self.batch
    }
}

/// 由 `bond_begin` / `bond_end` 构建 CSR 邻接。
///
/// **邻居顺序保持键的插入顺序**,不按邻居下标排序。这不是实现细节而是
/// 语义要求:SMILES 的 `@` / `@@` 手性含义依赖于邻居在字符串中出现的先后,
/// 一旦重排,手性解释就全错了。
fn build_csr(b: &mut MolBatch) {
    let n_atoms = b.atomic_num.len();
    let n_bonds = b.bond_begin.len();

    let mut offset = vec![0u32; n_atoms + 1];
    for i in 0..n_bonds {
        offset[b.bond_begin[i] as usize + 1] += 1;
        offset[b.bond_end[i] as usize + 1] += 1;
    }
    for i in 1..=n_atoms {
        offset[i] += offset[i - 1];
    }

    let total = 2 * n_bonds;
    let mut nbr_atom = vec![0u32; total];
    let mut nbr_bond = vec![0u32; total];
    let mut cursor = offset[..n_atoms].to_vec();

    // 按键下标递增填充 —— 这正是"插入顺序"
    for bi in 0..n_bonds {
        let (x, y) = (b.bond_begin[bi], b.bond_end[bi]);

        let p = cursor[x as usize] as usize;
        nbr_atom[p] = y;
        nbr_bond[p] = bi as u32;
        cursor[x as usize] += 1;

        let p = cursor[y as usize] as usize;
        nbr_atom[p] = x;
        nbr_bond[p] = bi as u32;
        cursor[y as usize] += 1;
    }

    b.nbr_offset = offset;
    b.nbr_atom = nbr_atom;
    b.nbr_bond = nbr_bond;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::BondOrder;

    /// 乙醇 CCO
    fn ethanol() -> MolBuilder {
        let mut m = MolBuilder::new();
        let c0 = m.add_atom(6);
        let c1 = m.add_atom(6);
        let o = m.add_atom(8);
        m.add_bond(c0, c1, BondOrder::Single).unwrap();
        m.add_bond(c1, o, BondOrder::Single).unwrap();
        m.set_name("ethanol");
        m
    }

    /// 苯 c1ccccc1(此处只建拓扑,芳香性属于 L2)
    fn benzene() -> MolBuilder {
        let mut m = MolBuilder::new();
        for _ in 0..6 {
            m.add_atom(6);
        }
        for i in 0..6u32 {
            m.add_bond(i, (i + 1) % 6, BondOrder::Aromatic).unwrap();
        }
        m
    }

    #[test]
    fn empty_batch() {
        let b = MolBatchBuilder::new().finish();
        assert_eq!(b.num_mols(), 0);
        assert_eq!(b.num_atoms(), 0);
        assert_eq!(b.num_bonds(), 0);
        assert!(b.is_empty());
        assert!(b.mol(0).is_none());
    }

    #[test]
    fn single_molecule_shape() {
        let mut bb = MolBatchBuilder::new();
        bb.push(&ethanol()).unwrap();
        let b = bb.finish();

        assert_eq!(b.num_mols(), 1);
        assert_eq!(b.num_atoms(), 3);
        assert_eq!(b.num_bonds(), 2);
        assert_eq!(b.atomic_nums(), &[6, 6, 8]);
        assert_eq!(b.mol_atom_offsets(), &[0, 3]);
    }

    #[test]
    fn multiple_molecules_get_global_indices() {
        let mut bb = MolBatchBuilder::new();
        let i0 = bb.push(&ethanol()).unwrap();
        let i1 = bb.push(&benzene()).unwrap();
        let i2 = bb.push(&ethanol()).unwrap();
        let b = bb.finish();

        assert_eq!((i0, i1, i2), (0, 1, 2));
        assert_eq!(b.num_mols(), 3);
        assert_eq!(b.num_atoms(), 3 + 6 + 3);
        assert_eq!(b.num_bonds(), 2 + 6 + 2);
        assert_eq!(b.mol_atom_offsets(), &[0, 3, 9, 12]);

        // 第二个分子的键端点必须已平移到全局
        assert_eq!(b.bond_begin[2], 3, "苯的第一条键应从全局原子 3 起");
    }

    #[test]
    fn csr_is_consistent() {
        let mut bb = MolBatchBuilder::new();
        bb.push(&ethanol()).unwrap();
        bb.push(&benzene()).unwrap();
        let b = bb.finish();

        // 偏移单调、末项等于 2*键数
        assert_eq!(b.nbr_offset.len(), b.num_atoms() + 1);
        assert!(b.nbr_offset.windows(2).all(|w| w[0] <= w[1]));
        assert_eq!(*b.nbr_offset.last().unwrap() as usize, 2 * b.num_bonds());

        // 度数之和 = 2 * 键数(握手定理)
        let deg_sum: u32 = (0..b.num_atoms())
            .map(|i| b.nbr_offset[i + 1] - b.nbr_offset[i])
            .sum();
        assert_eq!(deg_sum as usize, 2 * b.num_bonds());

        // 每条邻接项都能在 bond 列里对上
        for a in 0..b.num_atoms() as u32 {
            for k in b.nbr_offset[a as usize]..b.nbr_offset[a as usize + 1] {
                let nbr = b.nbr_atom[k as usize];
                let bond = b.nbr_bond[k as usize] as usize;
                let (x, y) = (b.bond_begin[bond], b.bond_end[bond]);
                assert!(
                    (x == a && y == nbr) || (y == a && x == nbr),
                    "原子 {a} 的邻接项 {nbr}(键 {bond})与键端点 ({x},{y}) 不符"
                );
            }
        }
    }

    /// 邻居顺序必须是键插入顺序 —— 手性解释依赖它。
    #[test]
    fn csr_preserves_bond_insertion_order() {
        let mut m = MolBuilder::new();
        for _ in 0..5 {
            m.add_atom(6);
        }
        // 刻意让中心原子 0 的邻居按 4,2,3,1 的顺序连上
        m.add_bond(0, 4, BondOrder::Single).unwrap();
        m.add_bond(0, 2, BondOrder::Single).unwrap();
        m.add_bond(0, 3, BondOrder::Single).unwrap();
        m.add_bond(0, 1, BondOrder::Single).unwrap();

        let mut bb = MolBatchBuilder::new();
        bb.push(&m).unwrap();
        let b = bb.finish();

        let nbrs: Vec<u32> = (b.nbr_offset[0]..b.nbr_offset[1])
            .map(|k| b.nbr_atom[k as usize])
            .collect();
        assert_eq!(nbrs, vec![4, 2, 3, 1], "邻居顺序被重排会破坏手性语义");
    }

    #[test]
    fn names_are_preserved() {
        let mut bb = MolBatchBuilder::new();
        bb.push(&ethanol()).unwrap();
        bb.push(&benzene()).unwrap();
        let b = bb.finish();
        assert_eq!(b.names[0].as_deref(), Some("ethanol"));
        assert_eq!(b.names[1], None);
    }

    #[test]
    fn molecule_with_no_bonds() {
        // 单原子分子(如 [Na+])不能把 CSR 建坏
        let mut m = MolBuilder::new();
        m.add_atom(11);
        let mut bb = MolBatchBuilder::new();
        bb.push(&m).unwrap();
        let b = bb.finish();

        assert_eq!(b.num_atoms(), 1);
        assert_eq!(b.num_bonds(), 0);
        assert_eq!(b.nbr_offset, vec![0, 0]);
    }
}
