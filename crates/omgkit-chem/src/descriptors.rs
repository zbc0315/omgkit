//! 图神经网络特征化要读的那一组原子/键描述符。
//!
//! # 这个模块为什么存在
//!
//! 描述符本身早就散在各处:元素表里有原子量,`AtomData` 上有形式电荷与杂化,
//! 标志位里有芳香与环成员,[`gasteiger`](crate::gasteiger) 里有部分电荷。
//! 散着的问题不是难找,是**每个调用方都要自己拼一遍**——而"总连接度含不含隐式
//! 氢""芳香键的键级填 1.5 还是另开一档"这类口径,拼十遍就会有十个说法,
//! 且分歧是静默的:模型照训,只是特征列和别人的对不上。
//!
//! 所以口径只写一处,由本模块交付。
//!
//! # 口径按外部实现钉死
//!
//! 每一项都与 RDKit 的同名取值逐位相同,由 `harness/check_descriptors.py` 拿
//! 全语料逐原子/逐键比。几处容易各拼各的:
//!
//! - **总连接度**是"显式邻居数 + 总氢数",不是度数(RDKit `GetTotalDegree`)。
//! - **总氢数**只算显式声明数与隐式推断数,**不含**图里那些独立的 `[H]` 原子
//!   (RDKit `GetTotalNumHs(includeNeighbors=false)`)。
//! - **原子量**在标了同位素时给的是**那个核素的精确质量**,没标才是标准原子量
//!   (RDKit `GetMass`)。两者在氘上差了一倍。
//!
//! # 前置:两步,而且第二步不在本 crate 里
//!
//! 分子要先 `sanitize`(填芳香、环、杂化、共轭、隐式氢),**再**跑
//! `omgkit_io::stereo::perceive_bond_stereo`(把 `/` `\` 的书写方向折算成双键
//! 自己的顺反)。第二步在 `omgkit-io` 里,而 `omgkit-chem` 是它的同级,调不到
//! ——所以本模块只**读** [`BondDescriptors::stereo`],永远不会自己去填它。
//!
//! 漏了第二步不会报错,只会让每根双键的顺反恒为"未指定"。Python 侧的
//! `Mol.sanitize()` 已经把两步并在一起,走那条路的调用方不必操心。
//!
//! # 顺反给的是 cis/trans,不是 Z/E
//!
//! `Z`/`E` 按 **CIP 优先级**定义,而 CIP 排序本仓库没有实现。这里交出的是
//! "相对记录下来的那两个参照原子"的顺反,与 RDKit
//! `SetBondStereoFromDirections` 给的 `STEREOCIS`/`STEREOTRANS` 同一口径,
//! 不是 `AssignStereochemistry` 给的 `STEREOZ`/`STEREOE`。
//!
//! 两者承载的**几何信息相同**——只要参照原子一并交出去,顺反与 Z/E 可以互相
//! 换算;差的只是标签的定法。所以 [`BondDescriptors::stereo_atoms`] 不是附赠,
//! 是这一项能用的前提。要 Z/E 的调用方得自己排 CIP。
//!
//! # 缺失是一种取值
//!
//! [`AtomDescriptors::electronegativity`] 与
//! [`AtomDescriptors::gasteiger_charge`] 都可能"算不出":前者是该元素没有公认
//! 的 Pauling 值,后者是该元素不在 Gasteiger 参数表里。两者都**如实交出缺失**
//! (`None` / 非有限值),不拿默认数顶——把"不知道"和"恰好是这个数"混成一格,
//! 下游就再也没法决定该不该屏蔽这一维。

use omgkit_core::{
    element, AtomFlags, BondData, BondFlags, BondOrder, BondStereo, ChiralTag, Hybridization,
    MolBuilder,
};

use crate::gasteiger::{gasteiger_charges, DEFAULT_ITERATIONS};

/// 单个原子的描述符。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AtomDescriptors {
    /// 元素种类(原子序数)。0 为通配原子 `*`。
    pub atomic_num: u8,
    /// 总连接度:显式邻居数 + 总氢数
    pub total_degree: u32,
    /// 形式电荷
    pub formal_charge: i8,
    /// 手性标记(几何类别,不是 CIP 的 R/S)
    pub chiral_tag: ChiralTag,
    /// 总氢数:显式声明 + 隐式推断,不含图里独立的 `[H]` 原子
    pub total_num_hs: u32,
    /// 杂化状态
    pub hybridization: Hybridization,
    /// 是否芳香
    pub is_aromatic: bool,
    /// 是否在环上
    pub is_in_ring: bool,
    /// 原子量。**标了同位素就用那个核素的精确质量**,没标才用标准原子量。
    ///
    /// 氘代化合物上两者差了一倍(1.008 对 2.0141),混用会在整列上悄悄偏移。
    /// 同位素标了但表里查不到那个质量数时,退回质量数本身(整数)——
    /// 与外部实现同一口径。
    pub mass: f64,
    /// Pauling 电负性;该元素没有公认值时为 `None`
    pub electronegativity: Option<f64>,
    /// Gasteiger 部分电荷。表外元素上为非有限值,见
    /// [`gasteiger_is_valid`](Self::gasteiger_is_valid)。
    pub gasteiger_charge: f64,
}

impl AtomDescriptors {
    /// 上一项电荷算不算得出来。
    ///
    /// 参数表覆盖不到的元素会给出 `NaN` / `inf` 并沿图扩散,所以这不是个别
    /// 金属原子的事——同一个分子里的碳也可能因此失效。特征化必须按这一位
    /// 决定要不要屏蔽电荷那一维。
    #[must_use]
    pub fn gasteiger_is_valid(&self) -> bool {
        self.gasteiger_charge.is_finite()
    }
}

/// 单根键的描述符。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BondDescriptors {
    /// 起点原子下标
    pub begin: u32,
    /// 终点原子下标
    pub end: u32,
    /// 键级
    pub order: BondOrder,
    /// 是否共轭
    pub is_conjugated: bool,
    /// 是否在环上
    pub is_in_ring: bool,
    /// 双键顺反
    pub stereo: BondStereo,
    /// 顺反是**相对这两个原子**说的;没有顺反时为 `None`。
    ///
    /// 少了它,上一项就没有意义:四取代双键上参照挑得不同,同一个几何会得出
    /// 相反的顺反值。所以两项必须一起交。
    pub stereo_atoms: Option<[u32; 2]>,
}

/// 原子量:标了同位素用精确质量,否则用标准原子量。
///
/// 三条分支都照外部实现:
/// 1. 没标同位素 → 标准原子量;
/// 2. 标了、表里有 → 该核素的精确质量;
/// 3. 标了、表里没有 → **质量数本身**(通配原子除外,它给 0)。
fn atom_mass(atomic_num: u8, isotope: u16, el: Option<&'static omgkit_core::Element>) -> f64 {
    if isotope == 0 {
        return el.map_or(0.0, |e| e.mass);
    }
    element::isotope_mass(atomic_num, isotope).unwrap_or(if atomic_num == 0 {
        0.0
    } else {
        f64::from(isotope)
    })
}

/// 逐原子算出全部原子描述符,与 `mol.atoms()` 等长、同序。
///
/// # 前置
///
/// 分子必须**净化过**。芳香、环成员、杂化、共轭、隐式氢数全部由净化填写,
/// 没净化的分子进来不报错,只会让每一项都是解析时的占位值。
#[must_use]
pub fn atom_descriptors(mol: &MolBuilder) -> Vec<AtomDescriptors> {
    let charges = gasteiger_charges(mol, DEFAULT_ITERATIONS);
    mol.atoms()
        .iter()
        .enumerate()
        .map(|(i, a)| {
            let total_num_hs = u32::from(a.num_explicit_hs) + u32::from(a.num_implicit_hs);
            let degree = u32::try_from(mol.degree(u32::try_from(i).unwrap_or(0))).unwrap_or(0);
            let el = element::by_atomic_num(a.atomic_num);
            AtomDescriptors {
                atomic_num: a.atomic_num,
                total_degree: degree + total_num_hs,
                formal_charge: a.formal_charge,
                chiral_tag: a.chiral_tag,
                total_num_hs,
                hybridization: a.hybridization,
                is_aromatic: a.flags.contains(AtomFlags::AROMATIC),
                is_in_ring: a.flags.contains(AtomFlags::IN_RING),
                mass: atom_mass(a.atomic_num, a.isotope, el),
                electronegativity: el.and_then(|e| e.electronegativity),
                gasteiger_charge: charges[i],
            }
        })
        .collect()
}

/// 逐键算出全部键描述符,与 `mol.bonds()` 等长、同序。
///
/// 前置同 [`atom_descriptors`]。
#[must_use]
pub fn bond_descriptors(mol: &MolBuilder) -> Vec<BondDescriptors> {
    mol.bonds()
        .iter()
        .map(|b| BondDescriptors {
            begin: b.begin,
            end: b.end,
            order: b.order,
            is_conjugated: b.flags.contains(BondFlags::CONJUGATED),
            is_in_ring: b.flags.contains(BondFlags::IN_RING),
            stereo: b.stereo,
            stereo_atoms: if b.stereo == BondStereo::None
                || b.stereo_atoms[0] == BondData::NO_STEREO_ATOM
                || b.stereo_atoms[1] == BondData::NO_STEREO_ATOM
            {
                None
            } else {
                Some(b.stereo_atoms)
            },
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use omgkit_io::smiles;

    /// 照抄产品入口的前置序列(Python 侧 `Mol.sanitize()` 就是这两步)。
    /// 只跑净化的话,顺反那一档恒为"未指定",而用例照样能全绿。
    fn sanitized(smi: &str) -> MolBuilder {
        let mut m = smiles::parse(smi).unwrap_or_else(|e| panic!("{smi}: {}", e.render()));
        crate::pipeline::sanitize(&mut m).unwrap_or_else(|e| panic!("{smi}: {e}"));
        omgkit_io::stereo::perceive_bond_stereo(&mut m);
        m
    }

    /// 总连接度含隐式氢,总氢数**不**含图里独立的 `[H]`。
    ///
    /// 两条口径写在一个用例里,因为它们最容易被拼成同一个数:甲烷写成 `C`
    /// 时四个氢是隐式的,写成 `[H]C([H])([H])[H]` 时是四个独立原子——
    /// **总连接度两种写法都该是 4**,而总氢数一个是 4、一个是 0。
    #[test]
    fn total_degree_counts_hydrogen_either_way_but_total_hs_does_not() {
        let implicit = atom_descriptors(&sanitized("C"));
        assert_eq!(implicit[0].total_degree, 4);
        assert_eq!(implicit[0].total_num_hs, 4);

        let explicit = atom_descriptors(&sanitized("[H]C([H])([H])[H]"));
        let carbon = explicit
            .iter()
            .find(|d| d.atomic_num == 6)
            .expect("该有一个碳");
        assert_eq!(carbon.total_degree, 4, "四个显式氢邻居也算进总连接度");
        assert_eq!(carbon.total_num_hs, 0, "独立的 [H] 原子不计入总氢数");
    }

    /// 电负性缺失与电荷失效各是一种取值,不许被顶成 0。
    #[test]
    fn missing_is_reported_as_missing() {
        let he = atom_descriptors(&sanitized("[He]"));
        assert_eq!(he[0].electronegativity, None, "氦没有公认的 Pauling 值");

        let na = atom_descriptors(&sanitized("[Na][Na]"));
        assert!(
            na.iter().any(|d| !d.gasteiger_is_valid()),
            "表外元素的电荷该报失效,实得 {:?}",
            na.iter().map(|d| d.gasteiger_charge).collect::<Vec<_>>()
        );
    }

    /// 标了同位素就得给那个核素的质量,不是标准原子量。
    ///
    /// 这一档不设用例的话,`mass` 永远读元素表也能全绿 —— 而氘代分子在语料里
    /// 是真实存在的(实测大语料里有十处),整列会悄悄差一倍。
    #[test]
    fn an_isotope_label_changes_the_mass() {
        let d = atom_descriptors(&sanitized("[2H]C"));
        let h = d.iter().find(|x| x.atomic_num == 1).expect("该有一个氢");
        assert!(
            (h.mass - 2.014_101_778).abs() < 1e-9,
            "氘的质量应是 2.0141,实得 {}",
            h.mass
        );
        let plain = atom_descriptors(&sanitized("[H]C"));
        let h0 = plain
            .iter()
            .find(|x| x.atomic_num == 1)
            .expect("该有一个氢");
        assert!(
            (h0.mass - 1.008).abs() < f64::EPSILON,
            "没标同位素该用标准原子量,实得 {}",
            h0.mass
        );
    }

    /// 苯:六根键全是芳香、全共轭、全在环上,且没有顺反。
    #[test]
    fn benzene_bonds_are_aromatic_conjugated_and_cyclic() {
        let d = bond_descriptors(&sanitized("c1ccccc1"));
        assert_eq!(d.len(), 6);
        for b in &d {
            assert_eq!(b.order, BondOrder::Aromatic);
            assert!(b.is_conjugated);
            assert!(b.is_in_ring);
            assert_eq!(b.stereo, BondStereo::None);
        }
    }

    /// 双键顺反读得回来——这一档不设用例的话,`stereo` 恒为 `None` 也能全绿。
    #[test]
    fn a_configured_double_bond_reports_its_geometry() {
        let d = bond_descriptors(&sanitized("C/C=C/C"));
        let ez: Vec<_> = d.iter().filter(|b| b.stereo != BondStereo::None).collect();
        assert_eq!(ez.len(), 1, "该有且只有一根键带顺反:{d:?}");
        assert_eq!(ez[0].order, BondOrder::Double);
        // 顺反离开参照原子没有意义,两项必须一起有
        assert!(ez[0].stereo_atoms.is_some(), "带顺反却没有参照原子");
        for b in &d {
            assert_eq!(
                b.stereo == BondStereo::None,
                b.stereo_atoms.is_none(),
                "顺反与参照原子必须同时有或同时没有:{b:?}"
            );
        }
    }
}
