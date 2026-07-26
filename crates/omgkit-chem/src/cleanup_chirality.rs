//! 净化第 11 步:剔除**几何上不成立**的立体标记。
//!
//! # 它管的是"几何对不对",不是"是不是真手性中心"
//!
//! 这个区分很容易搞混,而且搞混之后会去错误的方向找 bug。本步骤只做一件事:
//! 检查立体标记与原子的**几何**是否自洽 —— 四面体标记要求 sp³,配位几何要求
//! 配位数落在对应的区间里。不自洽就把标记清掉。
//!
//! 它**不**判断"这个中心的四个取代基是否两两可区分"。`[C@@H]1CCCC1` 的碳是
//! sp³、四个配位,几何完全成立,本步骤原样保留它的标记 —— 尽管环上两条通路
//! 等价、它根本不是真手性中心。剔除非真手性中心要跑对称性判定,不属于净化 ——
//! 那件事由 `omgkit_io::stereo::genuine_tetrahedral` 做,写规范 SMILES 时调用。
//!
//! # 排在杂化之后
//!
//! 判据里用到杂化,所以必须排在第 9 步之后。放前面的话杂化还是 `Unspecified`,
//! 所有四面体标记都会被误判成"不是 sp³"而清光。
//!
//! # 排列序号也要落在范围内
//!
//! 每种几何的排列数是固定的(平面四方 3、三角双锥 20、八面体 30)。超出范围的
//! 序号归零 —— 解析阶段已经拦掉了写出界的输入,但分子也可以由程序构造,
//! 那条路上没有解析器把关。

use omgkit_core::{ChiralTag, Hybridization, MolBuilder};

/// 剔除几何上不成立的立体标记(第 11 步)。
///
/// 必须排在杂化判定之后 —— 见模块文档。
///
/// 返回被改动的原子数。触发面很窄,调用方若要断言"它确实开了火",拿这个数比零。
pub fn cleanup_chirality(mol: &mut MolBuilder) -> usize {
    let mut changed = 0;
    for i in 0..mol.num_atoms() as u32 {
        let a = mol.atoms()[i as usize];
        let total_degree =
            mol.degree(i) + usize::from(a.num_explicit_hs) + usize::from(a.num_implicit_hs);

        let verdict = match a.chiral_tag {
            ChiralTag::Unspecified => continue,
            // 四面体要求 sp³。标记与杂化不符时,标记是错的 —— 杂化是从
            // 成键与孤对算出来的,比作者手写的标记可靠。
            ChiralTag::Cw | ChiralTag::Ccw => {
                if a.hybridization == Hybridization::Sp3 {
                    Verdict::Keep
                } else {
                    Verdict::Clear
                }
            }
            ChiralTag::SquarePlanar => geometry_verdict(total_degree, 4, a.stereo_perm, 3),
            ChiralTag::TrigonalBipyramidal => geometry_verdict(total_degree, 5, a.stereo_perm, 20),
            ChiralTag::Octahedral => geometry_verdict(total_degree, 6, a.stereo_perm, 30),
            // `@AL` 这类轴手性的几何判据不在本步骤的范围内
            ChiralTag::Other => continue,
        };

        match verdict {
            Verdict::Keep => {}
            Verdict::Clear => {
                if let Some(m) = mol.atom_mut(i) {
                    m.chiral_tag = ChiralTag::Unspecified;
                    m.stereo_perm = 0;
                }
                changed += 1;
            }
            Verdict::ResetPerm => {
                if let Some(m) = mol.atom_mut(i) {
                    m.stereo_perm = 0;
                }
                changed += 1;
            }
        }
    }
    changed
}

enum Verdict {
    Keep,
    /// 几何不成立 —— 整个标记作废
    Clear,
    /// 几何成立但排列序号出界 —— 只把序号归零,几何类别保留
    ResetPerm,
}

/// 配位几何的判据:配位数要落在 `2..=max_degree`,序号要落在 `0..=max_perm`。
///
/// 下界一律是 2 —— 一两个配体的中心谈不上什么几何。
fn geometry_verdict(degree: usize, max_degree: usize, perm: u8, max_perm: u8) -> Verdict {
    if !(2..=max_degree).contains(&degree) {
        Verdict::Clear
    } else if perm > max_perm {
        Verdict::ResetPerm
    } else {
        Verdict::Keep
    }
}

#[cfg(test)]
mod tests {
    use omgkit_io::smiles;

    use super::*;
    use crate::{
        assign_radicals, clean_up, conjugation::set_hybridization, kekulize, perceive_rings,
        set_aromaticity, set_conjugation, valence::update_property_cache,
    };

    /// 跑到第 9 步(含杂化),即本步骤的前置条件
    fn upto_step9(smi: &str) -> MolBuilder {
        let mut m = smiles::parse(smi).unwrap_or_else(|e| panic!("{smi}: {}", e.render()));
        clean_up(&mut m);
        update_property_cache(&mut m).expect("价键");
        let _ = perceive_rings(&mut m);
        kekulize(&mut m).expect("kekulize");
        assign_radicals(&mut m);
        set_aromaticity(&mut m);
        update_property_cache(&mut m).expect("收尾价键");
        set_conjugation(&mut m);
        set_hybridization(&mut m);
        m
    }

    fn tag_after(smi: &str, atom: usize) -> ChiralTag {
        let mut m = upto_step9(smi);
        cleanup_chirality(&mut m);
        m.atoms()[atom].chiral_tag
    }

    /// sp³ 的四面体标记要保留 —— **哪怕它不是真手性中心**。
    ///
    /// `[C@@H]1CCCC1` 环上两条通路等价,这个碳其实无手性可言。但那是
    /// `stereo::genuine_tetrahedral` 该判的事,本步骤只看几何。
    #[test]
    fn sp3_tetrahedral_is_kept_even_when_not_a_real_stereocentre() {
        // 断言的是"**没被动过**",不是某个具体取值 —— 标记的具体值取决于
        // 解析时的宇称补偿,与本步骤无关
        for (smi, atom) in [("N[C@@H](C)C(=O)O", 1), ("[C@@H]1CCCC1", 0)] {
            let before = upto_step9(smi).atoms()[atom].chiral_tag;
            assert!(before.is_tetrahedral(), "{smi}:前提是它带四面体标记");
            assert_eq!(
                tag_after(smi, atom),
                before,
                "{smi}:sp³ 的四面体标记应当原样保留 —— \
                 即使它其实不是真手性中心,那也是 genuine_tetrahedral 该判的事"
            );
        }
    }

    /// 杂化不是 sp³ 的四面体标记要清掉。
    #[test]
    fn non_sp3_tetrahedral_is_cleared() {
        // 芳香碳是 sp²,写在它上面的 `@` 站不住
        let mut m = upto_step9("c1ccccc1");
        if let Some(a) = m.atom_mut(0) {
            a.chiral_tag = ChiralTag::Cw;
        }
        assert_eq!(cleanup_chirality(&mut m), 1);
        assert_eq!(m.atoms()[0].chiral_tag, ChiralTag::Unspecified);
    }

    /// 配位几何:配位数出界就清掉。
    #[test]
    fn coordination_geometry_checks_the_degree() {
        // 平面四方允许 2..=4
        assert_eq!(
            tag_after("[Pt@SP1](Cl)(Cl)(N)N", 0),
            ChiralTag::SquarePlanar
        );
        assert_eq!(
            tag_after("[Pt@SP1](Cl)(Cl)(N)(N)Cl", 0),
            ChiralTag::Unspecified,
            "五配位的平面四方声称站不住"
        );
        assert_eq!(
            tag_after("[Pt@SP1]Cl", 0),
            ChiralTag::Unspecified,
            "一配位谈不上几何"
        );
        // 八面体允许 2..=6
        assert_eq!(
            tag_after("C[Co@OH25](N)(O)(S)(P)Cl", 1),
            ChiralTag::Octahedral
        );
        assert_eq!(
            tag_after("C[Co@OH25](N)(O)(S)(P)(Cl)Br", 1),
            ChiralTag::Unspecified,
            "七配位的八面体声称站不住"
        );
    }

    /// 排列序号出界时**只归零序号**,几何类别保留。
    ///
    /// 解析器已经拦掉了写出界的输入,但程序构造的分子没有解析器把关。
    #[test]
    fn out_of_range_permutation_is_reset_but_geometry_kept() {
        let mut m = upto_step9("[Pt@SP1](Cl)(Cl)(N)N");
        if let Some(a) = m.atom_mut(0) {
            a.stereo_perm = 99; // 平面四方最大 3
        }
        assert_eq!(cleanup_chirality(&mut m), 1);
        assert_eq!(m.atoms()[0].chiral_tag, ChiralTag::SquarePlanar, "类别保留");
        assert_eq!(m.atoms()[0].stereo_perm, 0, "序号归零");
    }

    /// 幂等
    #[test]
    fn is_idempotent() {
        for smi in [
            "N[C@@H](C)C(=O)O",
            "[Pt@SP1](Cl)(Cl)(N)N",
            "[Pt@SP1]Cl",
            "CCO",
        ] {
            let mut m = upto_step9(smi);
            cleanup_chirality(&mut m);
            let once: Vec<_> = m.atoms().to_vec();
            assert_eq!(cleanup_chirality(&mut m), 0, "{smi}:第二次不该再有改动");
            assert_eq!(m.atoms(), &once[..], "{smi}:不幂等");
        }
    }
}
