//! 净化第 3 步:把隐式氢算出来**写回分子**。
//!
//! **规则本身住在 [`omgkit_core::valence`]** —— `omgkit-io` 的写出侧也要用同一条
//! 规则(去掉方括号之后,氢数由读者按它反推),两处各写一遍必然静默分岔。
//! 这里只剩"跑一遍、写回去"这件事,以及它的产出类型。
//!
//! 严格模式:超价直接判为失败,整条净化终止。

use omgkit_core::valence::{explicit_valence_of, implicit_hs_of};
use omgkit_core::MolBuilder;

pub use omgkit_core::valence::{
    explicit_valence_nonstrict, implicit_hs_nonstrict, is_aromatic_atom, total_valence_nonstrict,
    valence_shift, ValenceError, ValenceErrorKind,
};

/// 第 3 步的产出。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ValenceResult {
    /// 每原子的显式价
    pub explicit_valence: Vec<i32>,
    /// 每原子的隐式氢数
    pub implicit_hs: Vec<u8>,
}

/// 计算全部原子的显式价与隐式氢数,并就地写回 `num_implicit_hs`。
///
/// # Errors
/// 任一原子超价时返回 [`ValenceError`],此时整条净化应当判为失败。
pub fn update_property_cache(mol: &mut MolBuilder) -> Result<ValenceResult, ValenceError> {
    let n = mol.num_atoms();
    let mut explicit_valence = vec![0i32; n];
    let mut implicit_hs = vec![0u8; n];

    for i in 0..n as u32 {
        let ev = explicit_valence_of(mol, i, true)?;
        explicit_valence[i as usize] = ev;
        let ih = implicit_hs_of(mol, i, ev, true)?;
        implicit_hs[i as usize] = ih;
    }

    for (i, &h) in implicit_hs.iter().enumerate() {
        if let Some(a) = mol.atom_mut(i as u32) {
            a.num_implicit_hs = h;
        }
    }

    Ok(ValenceResult {
        explicit_valence,
        implicit_hs,
    })
}

#[cfg(test)]
mod tests {
    use omgkit_io::smiles;

    use super::*;

    fn calc(smi: &str) -> (Vec<i32>, Vec<u8>) {
        let mut m = smiles::parse(smi).unwrap_or_else(|e| panic!("{}", e.render()));
        let r = update_property_cache(&mut m).unwrap_or_else(|e| panic!("{smi}: {e}"));
        (r.explicit_valence, r.implicit_hs)
    }

    #[test]
    fn simple_organics() {
        // (显式价, 隐式氢)
        assert_eq!(calc("CCO"), (vec![1, 2, 1], vec![3, 2, 1]));
        assert_eq!(calc("C"), (vec![0], vec![4]));
        assert_eq!(calc("CC(=O)O"), (vec![1, 4, 2, 1], vec![3, 0, 0, 1]));
    }

    #[test]
    fn bracket_atoms_get_no_implicit_hs() {
        // [CH4] 的氢是作者显式给的,不再推断
        let (ev, ih) = calc("[CH4]");
        assert_eq!(ev, vec![4]);
        assert_eq!(ih, vec![0], "NO_IMPLICIT 置位时隐式氢恒为 0");
    }

    #[test]
    fn aromatic_ring_atoms() {
        // 苯:芳香碳的显式价 3(1.5+1.5),补 1 个氢
        let (ev, ih) = calc("c1ccccc1");
        assert!(ev.iter().all(|&v| v == 3), "实际 {ev:?}");
        assert!(ih.iter().all(|&h| h == 1), "实际 {ih:?}");
    }

    #[test]
    fn kekulized_form_gives_same_counts() {
        // 凯库勒式与芳香式的隐式氢数应一致
        let (_, ih_arom) = calc("c1ccccc1");
        let (_, ih_kek) = calc("C1=CC=CC=C1");
        assert_eq!(ih_arom, ih_kek);
    }

    #[test]
    fn charged_atoms_use_effective_atomic_number() {
        // [NH4+] 的有效原子序数是 7-1=6(碳),故允许 4 价
        let (ev, ih) = calc("[NH4+]");
        assert_eq!(ev, vec![4]);
        assert_eq!(ih, vec![0]);
        // 铵的非方括号写法不存在,用 [O-] 验证负电
        let (ev, _) = calc("CC(=O)[O-]");
        assert_eq!(ev[3], 1);
    }

    #[test]
    fn hypervalent_sulfur() {
        // 硫酸:S 是 6 价
        let (ev, ih) = calc("O=S(=O)(O)O");
        assert_eq!(ev[1], 6, "S 应为 6 价");
        assert_eq!(ih[1], 0);
    }

    #[test]
    fn wildcard_gets_no_implicit_hs() {
        let (_, ih) = calc("*CC");
        assert_eq!(ih[0], 0, "通配原子不补氢");
    }

    #[test]
    fn overvalent_carbon_is_rejected() {
        // 五键碳:严格模式下必须判为失败
        let mut m = smiles::parse("C(C)(C)(C)(C)C").unwrap();
        let err = update_property_cache(&mut m).expect_err("五键碳应当被拒绝");
        assert_eq!(err.atom, 0);
        assert_eq!(err.kind, ValenceErrorKind::ExplicitValenceTooHigh);
        assert_eq!(err.valence, 5);
    }

    #[test]
    fn implicit_hs_written_back_to_atoms() {
        let mut m = smiles::parse("CCO").unwrap();
        let r = update_property_cache(&mut m).unwrap();
        for (i, a) in m.atoms().iter().enumerate() {
            assert_eq!(a.num_implicit_hs, r.implicit_hs[i]);
        }
    }

    #[test]
    fn dative_bond_asymmetry_flows_through() {
        // 配位键:给体不计价,受体计 1
        let mut m = smiles::parse("CC").unwrap();
        m.bond_mut(0)
            .unwrap()
            .set_order(omgkit_core::BondOrder::Dative);
        let r = update_property_cache(&mut m).unwrap();
        assert_eq!(r.explicit_valence[0], 0, "给体(起点)显式价为 0");
        assert_eq!(r.explicit_valence[1], 1, "受体(终点)显式价为 1");
        assert_eq!(r.implicit_hs, vec![4, 3]);
    }
}
