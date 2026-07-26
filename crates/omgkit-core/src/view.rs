//! 单分子的零拷贝视图。
//!
//! [`MolBatch`] 内部一律用全局下标,而算法和用户想要的是分子内局部下标。
//! [`MolView`] 只持有 `&MolBatch` 和分子号,把两者之间的换算收在一处 ——
//! 不复制任何数据,构造代价是几个字段的赋值。

use crate::batch::MolBatch;
use crate::builder::{AtomData, BondData, MolBuilder};

/// 批中单个分子的零拷贝视图。
///
/// 所有下标参数与返回值都是**分子内局部**下标(0 起)。
#[derive(Debug, Clone, Copy)]
pub struct MolView<'a> {
    batch: &'a MolBatch,
    idx: u32,
    atom_base: u32,
    bond_base: u32,
    n_atoms: u32,
    n_bonds: u32,
}

impl<'a> MolView<'a> {
    pub(crate) fn new(batch: &'a MolBatch, idx: u32) -> Self {
        let i = idx as usize;
        let atom_base = batch.mol_atom_offset[i];
        let bond_base = batch.mol_bond_offset[i];
        Self {
            batch,
            idx,
            atom_base,
            bond_base,
            n_atoms: batch.mol_atom_offset[i + 1] - atom_base,
            n_bonds: batch.mol_bond_offset[i + 1] - bond_base,
        }
    }

    /// 该分子在批中的下标
    #[must_use]
    pub fn index(&self) -> u32 {
        self.idx
    }

    /// 所属的批
    #[must_use]
    pub fn batch(&self) -> &'a MolBatch {
        self.batch
    }

    /// 原子数
    #[must_use]
    pub fn num_atoms(&self) -> usize {
        self.n_atoms as usize
    }

    /// 键数
    #[must_use]
    pub fn num_bonds(&self) -> usize {
        self.n_bonds as usize
    }

    /// 是否为空分子
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.n_atoms == 0
    }

    /// 分子名
    #[must_use]
    pub fn name(&self) -> Option<&'a str> {
        self.batch.names[self.idx as usize].as_deref()
    }

    /// 局部原子下标 → 全局下标。不做越界检查。
    #[must_use]
    pub fn global_atom(&self, local: u32) -> u32 {
        self.atom_base + local
    }

    /// 全局原子下标 → 局部下标。不属于本分子时返回 `None`。
    #[must_use]
    pub fn local_atom(&self, global: u32) -> Option<u32> {
        global
            .checked_sub(self.atom_base)
            .filter(|&l| l < self.n_atoms)
    }

    /// 取原子。越界返回 `None`。
    #[must_use]
    pub fn atom(&self, local: u32) -> Option<AtomData> {
        if local >= self.n_atoms {
            return None;
        }
        let g = (self.atom_base + local) as usize;
        let b = self.batch;
        Some(AtomData {
            atomic_num: b.atomic_num[g],
            formal_charge: b.formal_charge[g],
            isotope: b.isotope[g],
            num_explicit_hs: b.num_explicit_hs[g],
            num_implicit_hs: b.num_implicit_hs[g],
            num_radical_electrons: b.num_radical_electrons[g],
            atom_map: b.atom_map[g],
            chiral_tag: b.chiral_tag[g],
            stereo_perm: b.stereo_perm[g],
            hybridization: b.hybridization[g],
            flags: b.atom_flags[g],
        })
    }

    /// 取键,端点为**局部**下标。越界返回 `None`。
    #[must_use]
    pub fn bond(&self, local: u32) -> Option<BondData> {
        if local >= self.n_bonds {
            return None;
        }
        let g = (self.bond_base + local) as usize;
        let b = self.batch;
        Some(BondData {
            begin: b.bond_begin[g] - self.atom_base,
            end: b.bond_end[g] - self.atom_base,
            order: b.bond_order[g],
            direction: b.bond_direction[g],
            stereo: b.bond_stereo[g],
            stereo_atoms: b.bond_stereo_atoms[g].map(|a| {
                if a == BondData::NO_STEREO_ATOM {
                    a
                } else {
                    a - self.atom_base
                }
            }),
            flags: b.bond_flags[g],
        })
    }

    /// 原子的度(不含隐式氢)
    #[must_use]
    pub fn degree(&self, local: u32) -> usize {
        if local >= self.n_atoms {
            return 0;
        }
        let g = (self.atom_base + local) as usize;
        (self.batch.nbr_offset[g + 1] - self.batch.nbr_offset[g]) as usize
    }

    /// 遍历某原子的邻居,产出 `(邻居局部下标, 键局部下标)`。
    ///
    /// 顺序即键的插入顺序 —— 手性语义依赖于此,详见
    /// [`batch`](crate::batch) 模块中 CSR 构建的说明。
    pub fn neighbors(&self, local: u32) -> impl Iterator<Item = (u32, u32)> + 'a {
        let (start, end) = if local < self.n_atoms {
            let g = (self.atom_base + local) as usize;
            (self.batch.nbr_offset[g], self.batch.nbr_offset[g + 1])
        } else {
            (0, 0)
        };
        let b = self.batch;
        let (abase, bbase) = (self.atom_base, self.bond_base);
        (start..end).map(move |k| {
            let k = k as usize;
            (b.nbr_atom[k] - abase, b.nbr_bond[k] - bbase)
        })
    }

    /// 遍历全部原子,产出 `(局部下标, 原子)`。
    pub fn atoms(&self) -> impl Iterator<Item = (u32, AtomData)> + '_ {
        (0..self.n_atoms).map(move |i| (i, self.atom(i).expect("下标由 n_atoms 生成")))
    }

    /// 遍历全部键,产出 `(局部下标, 键)`。
    pub fn bonds(&self) -> impl Iterator<Item = (u32, BondData)> + '_ {
        (0..self.n_bonds).map(move |i| (i, self.bond(i).expect("下标由 n_bonds 生成")))
    }

    /// 拷贝回可变的 [`MolBuilder`],用于编辑或反应产物构建。
    ///
    /// 这是 `MolBatchBuilder::push` 的逆操作;二者的往返恒等性由测试保证。
    #[must_use]
    pub fn to_builder(&self) -> MolBuilder {
        let mut m = MolBuilder::with_capacity(self.num_atoms(), self.num_bonds());
        for (_, a) in self.atoms() {
            m.add_atom_data(a);
        }
        for (_, bd) in self.bonds() {
            m.add_bond_data(bd).expect("视图中的键端点必然合法");
        }
        if let Some(n) = self.name() {
            m.set_name(n);
        }
        m
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::batch::MolBatchBuilder;
    use crate::types::{
        AtomFlags, BondDirection, BondFlags, BondOrder, BondStereo, ChiralTag, Hybridization,
    };

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

    fn acetic_acid() -> MolBuilder {
        // CC(=O)O
        let mut m = MolBuilder::new();
        let c0 = m.add_atom(6);
        let c1 = m.add_atom(6);
        let o1 = m.add_atom(8);
        let o2 = m.add_atom(8);
        m.add_bond(c0, c1, BondOrder::Single).unwrap();
        m.add_bond(c1, o1, BondOrder::Double).unwrap();
        m.add_bond(c1, o2, BondOrder::Single).unwrap();
        m
    }

    fn batch_of(mols: &[MolBuilder]) -> MolBatch {
        let mut bb = MolBatchBuilder::new();
        for m in mols {
            bb.push(m).unwrap();
        }
        bb.finish()
    }

    #[test]
    fn local_indices_are_zero_based_per_molecule() {
        let b = batch_of(&[ethanol(), acetic_acid()]);
        let m1 = b.mol(1).unwrap();

        assert_eq!(m1.num_atoms(), 4);
        assert_eq!(m1.num_bonds(), 3);
        // 第二个分子的原子在全局是 3..7,局部必须是 0..4
        assert_eq!(m1.global_atom(0), 3);
        assert_eq!(m1.local_atom(3), Some(0));
        assert_eq!(m1.local_atom(2), None, "全局原子 2 属于第一个分子");
        assert_eq!(m1.local_atom(7), None, "全局原子 7 已越过本分子");

        let bond = m1.bond(1).unwrap();
        assert_eq!((bond.begin, bond.end), (1, 2), "键端点应为局部下标");
        assert_eq!(bond.order, BondOrder::Double);
    }

    #[test]
    fn neighbors_are_local() {
        let b = batch_of(&[ethanol(), acetic_acid()]);
        let m1 = b.mol(1).unwrap();

        // 乙酸的 C1 连着 C0、O1、O2
        let mut nbrs: Vec<u32> = m1.neighbors(1).map(|(a, _)| a).collect();
        nbrs.sort_unstable();
        assert_eq!(nbrs, vec![0, 2, 3]);
        assert_eq!(m1.degree(1), 3);
        assert_eq!(m1.degree(0), 1);
    }

    #[test]
    fn neighbor_bond_indices_are_local() {
        let b = batch_of(&[ethanol(), acetic_acid()]);
        let m1 = b.mol(1).unwrap();
        for (_, bond_local) in m1.neighbors(1) {
            assert!(
                (bond_local as usize) < m1.num_bonds(),
                "键下标 {bond_local} 未换算为局部"
            );
        }
    }

    #[test]
    fn out_of_range_access_is_none() {
        let b = batch_of(&[ethanol()]);
        let m = b.mol(0).unwrap();
        assert!(m.atom(3).is_none());
        assert!(m.bond(2).is_none());
        assert_eq!(m.degree(99), 0);
        assert_eq!(m.neighbors(99).count(), 0);
    }

    /// 造一个**每个字段都取非默认值**的分子。
    ///
    /// 这里刻意用**完整的结构体字面量**(不写 `..Default::default()`):
    /// 以后给 `AtomData` 或 `BondData` 加字段时,这里会直接编译失败,逼着人把
    /// 新字段也纳入往返测试。
    ///
    /// 之所以要这道保险:漏同步列式存储不会报错,只会让某个字段在
    /// builder → batch → builder 的往返中悄悄变回默认值,而这类 bug 要到很久
    /// 以后才以"化学算错了"的形式冒出来。
    fn every_field_set() -> MolBuilder {
        let mut m = MolBuilder::new();
        m.add_atom_data(AtomData {
            atomic_num: 7,
            formal_charge: -1,
            isotope: 15,
            num_explicit_hs: 2,
            num_implicit_hs: 3,
            num_radical_electrons: 1,
            atom_map: 7,
            chiral_tag: ChiralTag::Cw,
            stereo_perm: 0,
            hybridization: Hybridization::Sp3d2,
            flags: AtomFlags::AROMATIC | AtomFlags::NO_IMPLICIT | AtomFlags::IN_RING,
        });
        m.add_atom_data(AtomData {
            atomic_num: 16,
            formal_charge: 2,
            isotope: 34,
            num_explicit_hs: 1,
            num_implicit_hs: 0,
            num_radical_electrons: 2,
            atom_map: 3,
            chiral_tag: ChiralTag::Octahedral,
            stereo_perm: 25,
            hybridization: Hybridization::Sp2d,
            flags: AtomFlags::CONJUGATED,
        });
        m.add_bond_data(BondData {
            begin: 0,
            end: 1,
            order: BondOrder::Dative,
            direction: BondDirection::DownRight,
            stereo: BondStereo::Trans,
            // 故意不是 [0, 1]:参照原子要跟着基址平移,顺序反过来才能看出
            // 平移是否作用在了正确的位置上
            stereo_atoms: [1, 0],
            flags: BondFlags::AROMATIC | BondFlags::IN_RING | BondFlags::CONJUGATED,
        })
        .expect("端点合法");
        m.set_name("每字段非默认");
        m
    }

    /// builder → batch → view → builder 必须恒等。
    /// 这是 L0 最重要的性质:列式存储不能悄悄丢字段。
    #[test]
    fn roundtrip_through_batch_is_identity() {
        let m = every_field_set();

        // 夹在两个分子中间,确保偏移换算被真正考验
        let b = batch_of(&[ethanol(), m.clone(), ethanol()]);
        let back = b.mol(1).unwrap().to_builder();

        assert_eq!(back.atoms(), m.atoms(), "原子列往返不一致");
        assert_eq!(back.bonds(), m.bonds(), "键列往返不一致");
        assert_eq!(back.name(), m.name());
    }

    /// 确认上面那个分子**真的**每个字段都偏离了默认值 —— 否则往返测试会
    /// 在一个全是默认值的分子上"通过",什么也没验证。
    #[test]
    fn every_field_is_actually_non_default() {
        let m = every_field_set();
        let default_atom = AtomData::default();
        for (i, a) in m.atoms().iter().enumerate() {
            assert_ne!(a.atomic_num, default_atom.atomic_num, "原子{i}.元素");
            assert_ne!(a.formal_charge, default_atom.formal_charge, "原子{i}.电荷");
            assert_ne!(a.isotope, default_atom.isotope, "原子{i}.同位素");
            assert_ne!(a.atom_map, default_atom.atom_map, "原子{i}.映射号");
            assert_ne!(a.chiral_tag, default_atom.chiral_tag, "原子{i}.手性");
            assert_ne!(a.hybridization, default_atom.hybridization, "原子{i}.杂化");
            assert_ne!(a.flags, AtomFlags::NONE, "原子{i}.标志");
        }
        // 两个原子合起来覆盖显式氢/隐式氢/自由基的非零取值
        assert!(m.atoms().iter().any(|a| a.num_explicit_hs != 0));
        assert!(m.atoms().iter().any(|a| a.num_implicit_hs != 0));
        assert!(m.atoms().iter().any(|a| a.num_radical_electrons != 0));
        // 排列序号对四面体恒为 0,所以只能要求"至少一个原子非零",并且
        // 两类立体标记(四面体 / 配位几何)都要出现
        assert!(m.atoms().iter().any(|a| a.stereo_perm != 0));
        assert!(m.atoms().iter().any(|a| a.chiral_tag.is_tetrahedral()));
        assert!(m.atoms().iter().any(|a| !a.chiral_tag.is_tetrahedral()));

        for (i, b) in m.bonds().iter().enumerate() {
            assert_ne!(b.order, BondOrder::Unspecified, "键{i}.键级");
            assert_ne!(b.direction, BondDirection::None, "键{i}.方向");
            assert_ne!(b.stereo, BondStereo::None, "键{i}.立体");
            // 参照原子存的是全局下标,取回来必须换算成局部 —— 少了换算,
            // 批里第二个分子往后的参照就全指向别人家的原子
            assert_eq!(b.stereo_atoms, [1, 0], "键{i}.立体参照原子");
            assert_ne!(b.flags, BondFlags::NONE, "键{i}.标志");
        }
    }

    #[test]
    fn iter_visits_every_molecule() {
        let b = batch_of(&[ethanol(), acetic_acid(), ethanol()]);
        let sizes: Vec<usize> = b.iter().map(|m| m.num_atoms()).collect();
        assert_eq!(sizes, vec![3, 4, 3]);
        assert_eq!(b.iter().count(), 3);
    }

    #[test]
    fn try_mol_reports_out_of_range() {
        let b = batch_of(&[ethanol()]);
        assert!(b.try_mol(0).is_ok());
        let e = b.try_mol(5).unwrap_err();
        assert!(matches!(
            e,
            crate::error::Error::MolIndexOutOfRange {
                index: 5,
                num_mols: 1
            }
        ));
    }
}
