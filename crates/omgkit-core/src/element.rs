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

/// 这个原子在**三配位**状态下,第四个"取代基"是一对孤对电子 ——
/// 于是它照样可以是四面体立体中心。
///
/// 亚砜、亚磺酰胺、亚砜亚胺的 S,膦、膦氧化物的 P,以及同族的 As / Se / Te。
///
/// # 带正电的**目前**不算 —— 这是个已知的保守缺口,不是化学结论
///
/// 头一版这里写着"`[S+]` 三配位是平面的",那是**假话**:锍盐 R₃S⁺ 恰恰是
/// 三配体 + 一对孤对,构型稳定、可以拆分。真正没有孤对、四配体全在的是季铵 R₄N⁺。
///
/// 现在仍然排除带正电的,理由是**没有验证依据**:语料里带手性标记的三配位
/// 阳离子中心 0 个,RDKit 2022.09.5 也把 `C[S@+](C)CC` 的标记清成
/// `CHI_UNSPECIFIED` —— 外部判据看不见这一档,放开就是无据可依的改动。
/// 这条约定继承自 `omgkit-depict`(那边经 `check_wedge_readback.py` 验过)。
/// 要放开的话得先有能判它的判据。
///
/// # 为什么这条要放在 core
///
/// 它不是算法,是一张化学事实表,而**两个 crate 都要问同一个问题**:
/// `omgkit-depict` 从楔形反读构型时要知道"三根键也能定构型",
/// `omgkit-conf` 抽手性中心时要知道"三个邻居也算数"。
/// 各写一份的话迟早分岔,而分岔的表现是一半的中心画对了、另一半摆错了。
///
/// (实测:`omgkit-conf` 先前根本不认这一档 —— `<[u32; 4]>::try_from` 凑不够
/// 四个邻居就整个 `continue`,于是语料里 13 个分子、16 个中心的构型是掷硬币。)
#[must_use]
pub fn has_stereogenic_lone_pair(atomic_num: u8, formal_charge: i8) -> bool {
    formal_charge <= 0 && matches!(atomic_num, 15 | 16 | 33 | 34 | 52) // P S As Se Te
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn 孤对立体中心只认那五个元素() {
        for (z, chg, want) in [
            (16u8, 0i8, true), // S:亚砜
            (15, 0, true),     // P:膦
            (33, 0, true),     // As
            (34, 0, true),     // Se
            (52, 0, true),     // Te
            // 下面两条锁的是**当前的保守取舍**,不是化学结论:
            // 锍盐 R₃S⁺ 其实有孤对、构型稳定,只是外部判据看不见这一档。
            (16, 1, false), // [S+]:已知缺口,见函数文档
            (15, 1, false), // [P+]:四配位的鏻盐不走这一支
            (7, 0, false),  // N:孤对翻转太快,不当立体中心
            (6, 0, false),  // C:三配位是 sp²
            (8, 0, false),  // O:三配位是 [O+]
        ] {
            assert_eq!(
                has_stereogenic_lone_pair(z, chg),
                want,
                "元素 {z} 电荷 {chg}"
            );
        }
    }

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
