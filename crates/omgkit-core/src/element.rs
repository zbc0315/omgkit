//! 元素周期表。
//!
//! 数据由 `harness/gen_elements.py` 生成到
//! `element_data`。**默认价表直接决定 L2 的隐式氢推断**,任何一处
//! 数字不一致都会在 ChEMBL 差分测试里表现为成千条难以定位的分歧,
//! 所以它是生成的而非手写的。

use crate::element_data::{symbol_to_atomic_num, ELEMENTS, IS_EARLY_ATOM};

/// 单个元素的静态数据。
#[derive(Debug, Clone, Copy)]
pub struct Element {
    /// 原子序数。0 为 SMILES 通配原子 `*`。
    pub atomic_num: u8,
    /// 元素符号(首字母大写形式)
    pub symbol: &'static str,
    /// 周期
    pub period: u8,
    /// 共价半径 (Å)
    pub rcov: f32,
    /// 范德华半径 (Å)
    pub rvdw: f32,
    /// 标准原子量
    pub mass: f32,
    /// 外层电子数
    pub outer_electrons: u8,
    /// 最常见同位素的质量数
    pub common_isotope: u16,
    /// 最常见同位素的精确质量
    pub common_isotope_mass: f64,
    /// 默认价列表。`-1` 表示无价约束(该元素不参与隐式氢推断)。
    pub valences: &'static [i8],
}

impl Element {
    /// 该元素是否有价约束。无约束的元素(多数金属)不推断隐式氢。
    #[must_use]
    pub fn has_valence_constraint(&self) -> bool {
        !self.valences.is_empty() && self.valences[0] != -1
    }

    /// 给定已用价数,返回应当采用的默认价。
    ///
    /// 规则:取第一个 ≥ 已用价的默认价;若全部小于已用价,
    /// 说明是超价(hypervalent),返回 `None`,由调用方决定是报错还是放行。
    #[must_use]
    pub fn default_valence_for(&self, used: i8) -> Option<i8> {
        if !self.has_valence_constraint() {
            return None;
        }
        self.valences.iter().copied().find(|&v| v >= used)
    }
}

/// 按原子序数取元素。越界返回 `None`。
#[must_use]
pub fn by_atomic_num(atomic_num: u8) -> Option<&'static Element> {
    ELEMENTS.get(atomic_num as usize)
}

/// 按元素符号取元素(大小写敏感,如 `"Cl"`)。
#[must_use]
pub fn by_symbol(symbol: &str) -> Option<&'static Element> {
    symbol_to_atomic_num(symbol).and_then(by_atomic_num)
}

/// 元素符号 → 原子序数。
#[must_use]
pub fn atomic_num_of(symbol: &str) -> Option<u8> {
    symbol_to_atomic_num(symbol)
}

/// 表中元素总数(含 0 号通配原子)。
#[must_use]
pub fn count() -> usize {
    ELEMENTS.len()
}

/// 元素是否位于周期表中碳的左侧。
///
/// 用于 kekulize 的 `markDbondCands`:早期元素的形式电荷参与默认价计算时
/// 要取反。表由 `harness/gen_elements.py` 从
/// `Code/GraphMol/Atom.cpp` 抽取。
#[must_use]
pub fn is_early_atom(atomic_num: u8) -> bool {
    IS_EARLY_ATOM
        .get(atomic_num as usize)
        .copied()
        .unwrap_or(false)
}

/// SMILES **有机子集** —— 这些元素可以不加方括号直接书写,
/// 其隐式氢由默认价推断。其余元素必须写在 `[...]` 中。
///
/// 见 OpenSMILES 规范 §3.1.5。注意芳香形式 `b c n o p s` 也属于此集合。
#[must_use]
pub fn is_organic_subset(atomic_num: u8) -> bool {
    matches!(
        atomic_num,
        5 | 6 | 7 | 8 | 15 | 16 | 9 | 17 | 35 | 53 // B C N O P S F Cl Br I
    )
}

/// 可以以芳香小写形式出现在 SMILES 中的元素。
///
/// 这是**语法**层面的白名单(解析器用),不是芳香性感知的结果 ——
/// 后者属于 L2。
#[must_use]
pub fn can_be_aromatic_lowercase(atomic_num: u8) -> bool {
    matches!(atomic_num, 5 | 6 | 7 | 8 | 15 | 16 | 33 | 34 | 52) // b c n o p s as se te
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_is_indexed_by_atomic_number() {
        for (i, e) in ELEMENTS.iter().enumerate() {
            assert_eq!(e.atomic_num as usize, i, "第 {i} 项的原子序数错位");
        }
    }

    #[test]
    fn covers_through_oganesson() {
        // 118 号 + 0 号通配原子
        assert_eq!(count(), 119, "元素表应覆盖 0..=118");
        assert_eq!(by_atomic_num(118).unwrap().symbol, "Og");
        assert_eq!(by_atomic_num(0).unwrap().symbol, "*");
    }

    #[test]
    fn symbol_roundtrip() {
        for e in ELEMENTS.iter() {
            assert_eq!(
                atomic_num_of(e.symbol),
                Some(e.atomic_num),
                "符号 {} 无法反查",
                e.symbol
            );
        }
    }

    /// 这些默认价是隐式氢推断的基石。
    /// 此测试是防止生成器回归的护栏。
    #[test]
    fn organic_subset_default_valences() {
        assert_eq!(by_symbol("C").unwrap().valences, &[4]);
        assert_eq!(by_symbol("N").unwrap().valences, &[3]);
        assert_eq!(by_symbol("O").unwrap().valences, &[2]);
        assert_eq!(by_symbol("F").unwrap().valences, &[1]);
        assert_eq!(by_symbol("B").unwrap().valences, &[3]);
        assert_eq!(by_symbol("H").unwrap().valences, &[1]);
    }

    #[test]
    fn default_valence_selection() {
        let s = by_symbol("S").unwrap();
        // S 是多价的;应选第一个 ≥ 已用价的
        assert!(
            s.valences.len() > 1,
            "S 应有多个默认价,实际 {:?}",
            s.valences
        );
        assert_eq!(s.default_valence_for(2), Some(2));
        assert_eq!(s.default_valence_for(3), Some(s.valences[1]));

        let c = by_symbol("C").unwrap();
        assert_eq!(c.default_valence_for(4), Some(4));
        // 超价:无可用默认价
        assert_eq!(c.default_valence_for(5), None);
    }

    #[test]
    fn unknown_symbol_is_none() {
        assert!(by_symbol("Xx").is_none());
        assert!(
            by_symbol("c").is_none(),
            "小写芳香符号应由调用方归一化后再查表"
        );
    }

    #[test]
    fn organic_subset_membership() {
        for sym in ["B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"] {
            let a = atomic_num_of(sym).unwrap();
            assert!(is_organic_subset(a), "{sym} 应属于有机子集");
        }
        for sym in ["Si", "Se", "Na", "Fe"] {
            let a = atomic_num_of(sym).unwrap();
            assert!(!is_organic_subset(a), "{sym} 不应属于有机子集");
        }
    }
}
