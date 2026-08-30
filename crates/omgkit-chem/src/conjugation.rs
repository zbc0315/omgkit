//! 共轭标记与杂化状态(净化第 8、9 步)。
//!
//! 两步放在一个模块里,因为第 9 步要读第 8 步的结果:一个成键数为 4 的原子,
//! 带共轭键时判为 sp²,否则判为 sp³。
//!
//! # 共轭的判定
//!
//! 共轭不是全局的电子离域计算,而是一条**局部**判据:
//!
//! > 若原子 A 够格共轭,且它的一条重键(键级 ≥ 1.5)通向同样够格的 B,
//! > 同时 A 还有另一条键通向够格的 C,那么这两条键都标为共轭。
//!
//! 也就是"重键 — 原子 — 另一条键"这个片段。芳香键在起点就直接算共轭。
//!
//! # 杂化的判定
//!
//! 基于**成键数 + 孤对数**,不是几何构型:
//!
//! | 轨道数 | 杂化 |
//! |---|---|
//! | 0、1 | s |
//! | 2 | sp |
//! | 3 | sp² |
//! | 4 | sp³,但有共轭键且总度数 ≤ 3 时降为 sp² |
//! | 5 | sp³d |
//! | 6 | sp³d² |
//!
//! 立体标记若与配位数吻合,直接据此定杂化 —— 立体构型比电子计数更可靠。
//! **不吻合就不能采信**:标记只是作者的声称,配位数对不上时它是错的,这时
//! 退回电子计数。四面体要求配位数恰为 4,`@SP`/`@TB`/`@OH` 分别要求
//! 2–4 / 2–5 / 2–6。少了这道护栏,`[Pt@SP1](Cl)(Cl)(N)(N)Cl` 这种五配位的
//! 平面四方声称会被照单全收。
//!
//! # 三个容易写错的地方
//!
//! **1. 共轭候选判据对第三周期及以后的元素另有限制。** 外层电子数为 5 或 6
//! 的重元素(P、S 等)不算共轭候选,除非它是外层 6 且总度数 < 2。少了这条,
//! `Pc1ccccc1` 里的 C—P 键会被误标共轭,进而把磷判成 sp²。
//!
//! **2. 轨道数为 4 时降级为 sp² 还要看总度数。** 只有总度数 ≤ 3 才降 ——
//! 否则像 `CP1(C)=CC=CN=C1C` 里的磷会被错判。
//!
//! **3. 锕系及更重的元素(原子序数 ≥ 89)直接用总度数**,不做孤对推断 ——
//! 那些元素的价电子结构没有可靠的简单模型。八配位的钍配合物会因此落到
//! "轨道数 8"这一档,判为未确定,而不是被硬塞进 sp³d²。

use omgkit_core::{element, AtomFlags, BondFlags, Hybridization, MolBuilder};

use crate::aromaticity::count_pi_electrons;
use crate::valence::total_valence_nonstrict;

/// 标记共轭键(第 8 步)。
///
/// 先把所有键的共轭标志重置为"与芳香标志相同",再逐原子向外扩散。
pub fn set_conjugation(mol: &mut MolBuilder) {
    for bi in 0..mol.num_bonds() as u32 {
        let aromatic = mol.bonds()[bi as usize].flags.contains(BondFlags::AROMATIC);
        if let Some(mut b) = mol.bond_mut(bi) {
            b.flags_mut().set(BondFlags::CONJUGATED, aromatic);
        }
    }

    let candidate: Vec<bool> = (0..mol.num_atoms() as u32)
        .map(|i| is_conjugation_candidate(mol, i))
        .collect();

    for i in 0..mol.num_atoms() as u32 {
        mark_conjugated_around(mol, i, &candidate);
    }
}

/// 判定杂化状态(第 9 步)。必须排在 [`set_conjugation`] 之后。
pub fn set_hybridization(mol: &mut MolBuilder) {
    for i in 0..mol.num_atoms() as u32 {
        let h = hybridization_of(mol, i);
        if let Some(a) = mol.atom_mut(i) {
            a.hybridization = h;
        }
    }
}

// ---------------------------------------------------------------------------

/// 总度数 = 重原子邻居 + 全部氢
fn total_degree(mol: &MolBuilder, idx: u32) -> i32 {
    let a = mol.atoms()[idx as usize];
    mol.degree(idx) as i32 + i32::from(a.num_explicit_hs) + i32::from(a.num_implicit_hs)
}

/// 该原子是否够格参与共轭。
fn is_conjugation_candidate(mol: &MolBuilder, idx: u32) -> bool {
    let atom = mol.atoms()[idx as usize];
    let z = atom.atomic_num;
    let valences = element::by_atomic_num(z).map_or(&[-1i8][..], |e| e.valences);
    let min_valence = valences.first().map_or(-1, |&v| i32::from(v));

    // 中性且超过最低价 —— 超价就不共轭
    if atom.formal_charge == 0
        && min_valence >= 0
        && total_valence_nonstrict(mol, idx) > min_valence
    {
        return false;
    }

    // 第三周期及以后、外层电子 5 或 6 的元素(P、S 等)另有限制 —— 见模块文档
    let n_outer = i32::from(element::by_atomic_num(z).map_or(0, |e| e.outer_electrons));
    let row_ok =
        z <= 10 || (n_outer != 5 && n_outer != 6) || (n_outer == 6 && total_degree(mol, idx) < 2);

    row_ok && count_pi_electrons(mol, idx).is_some_and(|n| n > 0)
}

/// 以 `idx` 为中心,标记"重键 — 中心 — 另一条键"这个片段的两条键。
fn mark_conjugated_around(mol: &mut MolBuilder, idx: u32, candidate: &[bool]) {
    if !candidate[idx as usize] {
        return;
    }
    // 中心原子必须有 2 或 3 个取代基
    let sbo = mol.degree(idx) as i32
        + i32::from(mol.atoms()[idx as usize].num_explicit_hs)
        + i32::from(mol.atoms()[idx as usize].num_implicit_hs);
    if !(2..=3).contains(&sbo) {
        return;
    }

    let nbrs: Vec<(u32, u32)> = mol.neighbors(idx).collect();
    let mut to_mark: Vec<u32> = Vec::new();

    for &(other1, b1) in &nbrs {
        // 第一条键必须是重键,且对端够格
        if mol.bonds()[b1 as usize].valence_contribution_to(idx) < 1.5
            || !candidate[other1 as usize]
        {
            continue;
        }
        for &(other2, b2) in &nbrs {
            if b1 == b2 {
                continue;
            }
            if total_degree(mol, other2) > 3 || !candidate[other2 as usize] {
                continue;
            }
            to_mark.push(b1);
            to_mark.push(b2);
        }
    }

    for b in to_mark {
        if let Some(mut bond) = mol.bond_mut(b) {
            bond.flags_mut().insert(BondFlags::CONJUGATED);
        }
    }
}

/// 该原子是否关联任何共轭键
fn has_conjugated_bond(mol: &MolBuilder, idx: u32) -> bool {
    mol.neighbors(idx).any(|(_, bi)| {
        mol.bonds()[bi as usize]
            .flags
            .contains(BondFlags::CONJUGATED)
    })
}

/// 成键数 + 孤对数,即杂化所需的轨道数。
fn orbital_count(mol: &MolBuilder, idx: u32) -> i32 {
    let atom = mol.atoms()[idx as usize];
    let mut deg = total_degree(mol, idx);
    // 零价键与配位键的给体端不占轨道
    for (_, bi) in mol.neighbors(idx) {
        let b = mol.bonds()[bi as usize];
        if b.order == omgkit_core::BondOrder::Unspecified
            || (b.order == omgkit_core::BondOrder::Dative && b.end != idx)
        {
            deg -= 1;
        }
    }
    if atom.atomic_num <= 1 {
        return deg;
    }

    let n_outer =
        i32::from(element::by_atomic_num(atom.atomic_num).map_or(0, |e| e.outer_electrons));
    let total_valence = total_valence_nonstrict(mol, idx);
    let charge = i32::from(atom.formal_charge);
    let free_electrons = n_outer - (total_valence + charge);

    if total_valence + n_outer - charge < 8 {
        // 未满八隅体,要把自由基电子单独计入
        let n_radicals = i32::from(atom.num_radical_electrons);
        deg + (free_electrons - n_radicals) / 2 + n_radicals
    } else {
        deg + free_electrons / 2
    }
}

fn hybridization_of(mol: &MolBuilder, idx: u32) -> Hybridization {
    let atom = mol.atoms()[idx as usize];
    if atom.atomic_num == 0 {
        return Hybridization::Unspecified;
    }

    // 立体标记直接给出了配位几何,比电子计数更可靠 —— 作者写下 `@OH` 就是在
    // 断言这是个八面体中心,而八面体的电子计数(过渡金属的 d 电子、反馈键)
    // 恰恰是 `orbital_count` 最不靠谱的地方。
    //
    // **但必须先验配位数**。标记只是作者的声称,配位数对不上时它就是错的,
    // 这时宁可退回电子计数。下界一律是 2:一两个配体的中心谈不上什么几何。
    let deg = total_degree(mol, idx);
    match atom.chiral_tag {
        omgkit_core::ChiralTag::SquarePlanar if (2..=4).contains(&deg) => {
            return Hybridization::Sp2d
        }
        omgkit_core::ChiralTag::TrigonalBipyramidal if (2..=5).contains(&deg) => {
            return Hybridization::Sp3d
        }
        omgkit_core::ChiralTag::Octahedral if (2..=6).contains(&deg) => {
            return Hybridization::Sp3d2
        }
        // 四面体只在配位数确实是 4 时采信;`@` 也可能写在三配位的
        // 亚砜、膦上,那里的几何要靠孤对推断
        t if t.is_tetrahedral() && deg == 4 => return Hybridization::Sp3,
        _ => {}
    }

    // 锕系及更重的元素没有可靠的孤对推断,直接用总度数
    let norbs = if atom.atomic_num < 89 {
        orbital_count(mol, idx)
    } else {
        total_degree(mol, idx)
    };

    match norbs {
        0 | 1 => Hybridization::S,
        2 => Hybridization::Sp,
        3 => Hybridization::Sp2,
        4 => {
            // 有共轭键时降为 sp²,但总度数超过 3 的不降 —— 见模块文档
            if total_degree(mol, idx) > 3 || !has_conjugated_bond(mol, idx) {
                Hybridization::Sp3
            } else {
                Hybridization::Sp2
            }
        }
        5 => Hybridization::Sp3d,
        6 => Hybridization::Sp3d2,
        _ => Hybridization::Unspecified,
    }
}

/// 原子是否位于共轭体系中 —— 供上层查询。
#[must_use]
pub fn atom_is_conjugated(mol: &MolBuilder, idx: u32) -> bool {
    has_conjugated_bond(mol, idx)
}

/// 标记共轭原子:关联任何共轭键的原子都置 [`AtomFlags::CONJUGATED`]。
pub fn mark_conjugated_atoms(mol: &mut MolBuilder) {
    for i in 0..mol.num_atoms() as u32 {
        let c = has_conjugated_bond(mol, i);
        if let Some(a) = mol.atom_mut(i) {
            a.flags.set(AtomFlags::CONJUGATED, c);
        }
    }
}

#[cfg(test)]
mod tests {
    use omgkit_io::smiles;

    use super::*;
    use crate::valence::update_property_cache;
    use crate::{assign_radicals, clean_up, kekulize, perceive_rings, set_aromaticity};

    /// 跑完整的第 1/3/4/5/6/7/8/9 步前缀
    fn pipeline(smi: &str) -> MolBuilder {
        let mut m = smiles::parse(smi).unwrap_or_else(|e| panic!("{}", e.render()));
        clean_up(&mut m);
        update_property_cache(&mut m).expect("价键校验应通过");
        let _ = perceive_rings(&mut m);
        kekulize(&mut m).expect("应能 kekulize");
        assign_radicals(&mut m);
        set_aromaticity(&mut m);
        set_conjugation(&mut m);
        set_hybridization(&mut m);
        m
    }

    fn hybrids(smi: &str) -> Vec<Hybridization> {
        pipeline(smi)
            .atoms()
            .iter()
            .map(|a| a.hybridization)
            .collect()
    }

    fn conjugated(smi: &str) -> Vec<bool> {
        pipeline(smi)
            .bonds()
            .iter()
            .map(|b| b.flags.contains(BondFlags::CONJUGATED))
            .collect()
    }

    #[test]
    fn simple_hybridizations() {
        assert_eq!(hybrids("CC")[0], Hybridization::Sp3, "乙烷");
        assert_eq!(hybrids("C=C")[0], Hybridization::Sp2, "乙烯");
        assert_eq!(hybrids("C#C")[0], Hybridization::Sp, "乙炔");
        assert_eq!(hybrids("CO")[1], Hybridization::Sp3, "甲醇的氧");
        assert_eq!(hybrids("C#N")[1], Hybridization::Sp, "氰基的氮");
    }

    #[test]
    fn aromatic_atoms_are_sp2() {
        assert!(hybrids("c1ccccc1").iter().all(|&h| h == Hybridization::Sp2));
        assert!(hybrids("c1ccncc1").iter().all(|&h| h == Hybridization::Sp2));
    }

    /// 芳香键必定共轭。
    #[test]
    fn aromatic_bonds_are_conjugated() {
        assert!(conjugated("c1ccccc1").iter().all(|&c| c));
        assert!(conjugated("c1ccc2ccccc2c1").iter().all(|&c| c));
    }

    /// 1,3-丁二烯:三条键全部共轭。
    #[test]
    fn butadiene_is_fully_conjugated() {
        assert_eq!(conjugated("C=CC=C"), vec![true, true, true]);
    }

    /// 孤立双键不共轭 —— 共轭要求"重键 — 原子 — 另一条键"的片段两端都够格。
    #[test]
    fn isolated_double_bonds_are_not_conjugated() {
        assert_eq!(conjugated("CC=CC"), vec![false, false, false], "2-丁烯");
        assert_eq!(
            conjugated("C=CCC=C"),
            vec![false, false, false, false],
            "1,4-戊二烯:中间的 sp3 碳打断共轭"
        );
        assert_eq!(conjugated("CCO"), vec![false, false]);
    }

    /// 三键之间同样共轭。
    #[test]
    fn triple_bonds_conjugate() {
        assert_eq!(conjugated("C#CC#C"), vec![true, true, true], "丁二炔");
        assert_eq!(conjugated("C=CC#N"), vec![true, true, true], "丙烯腈");
    }

    /// 第三周期元素的限制:磷不是共轭候选,所以 C—P 键不共轭,磷保持 sp³。
    ///
    /// 少了这条限制,苯基膦的 C—P 键会被误标共轭,进而把磷判成 sp²。
    #[test]
    fn heavy_elements_do_not_extend_conjugation() {
        let c = conjugated("Pc1ccccc1");
        assert!(!c[0], "C—P 键不应共轭,实际 {c:?}");
        assert!(c[1..].iter().all(|&x| x), "苯环内部仍应全部共轭");
        assert_eq!(hybrids("Pc1ccccc1")[0], Hybridization::Sp3, "磷应为 sp³");
    }

    /// 轨道数为 4 时降级为 sp² 还要看总度数:总度数超过 3 的不降。
    #[test]
    fn sp3_downgrade_requires_low_total_degree() {
        let h = hybrids("CP1(C)=CC=CN=C1C");
        assert_eq!(h[1], Hybridization::Sp3, "磷总度数为 4,不应降为 sp²");
    }

    /// 羰基与相邻的杂原子共轭 —— 羧酸里的两条 C—O 键都共轭。
    #[test]
    fn carboxyl_is_conjugated() {
        let c = conjugated("CC(=O)O");
        assert!(c[1] && c[2], "C=O 与 C—O 都应共轭,实际 {c:?}");
    }

    /// 幂等
    #[test]
    fn is_idempotent() {
        for smi in ["c1ccccc1", "C=CC=C", "CC(=O)O", "CCO"] {
            let mut m = pipeline(smi);
            let once: Vec<_> = m.atoms().iter().map(|a| a.hybridization).collect();
            let once_b: Vec<_> = m.bonds().iter().map(|b| b.flags).collect();
            set_conjugation(&mut m);
            set_hybridization(&mut m);
            let twice: Vec<_> = m.atoms().iter().map(|a| a.hybridization).collect();
            let twice_b: Vec<_> = m.bonds().iter().map(|b| b.flags).collect();
            assert_eq!(once, twice, "{smi}: 杂化不幂等");
            assert_eq!(once_b, twice_b, "{smi}: 共轭不幂等");
        }
    }

    /// 通配原子不判杂化。
    #[test]
    fn dummy_atoms_have_no_hybridization() {
        assert_eq!(hybrids("C[*]")[1], Hybridization::Unspecified);
    }
}
