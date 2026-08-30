//! 元素周期表。
//!
//! 数据由 `harness/gen_elements.py` 生成到
//! `element_data`。**默认价表直接决定 L2 的隐式氢推断**,任何一处
//! 数字不一致都会在 ChEMBL 差分测试里表现为成千条难以定位的分歧,
//! 所以它是生成的而非手写的。

use crate::element_data::{symbol_to_atomic_num, ELEMENTS, ISOTOPE_MASSES, IS_EARLY_ATOM};

/// 单个元素的静态数据。
#[derive(Debug, Clone, Copy)]
pub struct Element {
    /// 原子序数。0 为 SMILES 通配原子 `*`。
    pub atomic_num: u8,
    /// 元素符号(首字母大写形式)
    pub symbol: &'static str,
    /// 周期
    pub period: u8,
    /// 共价半径 (Å)。
    ///
    /// **`1.9` 是上游的"未知"哨兵,不是一个测量值** —— `atomic_data.cpp:34`
    /// 原话是 `rCov (…). 1.9 if unknown.`。转录时那句声明丢了,于是表里 19 处
    /// `1.9` 里,Z=97..112 那连续 16 个是"没有数据",而 Y / Tm / Np 恰好真的
    /// 就是 1.90 —— 下游分不出这两件事。
    ///
    /// 同一个 struct 里 [`electronegativity`](Self::electronegativity) 花五行
    /// 文档反对的正是这种做法(用一个魔数当"不知道")。这里没有改成 `Option`,
    /// 是因为下游(`omgkit-conf` 的界矩阵)拿它当兜底模型用,给 `None` 就得
    /// 在那里再编一个数 —— 换汤不换药。**要改就得连同"没有共价半径的元素
    /// 该怎么建界"一起改**,那是另一件事。
    pub rcov: f32,
    /// 范德华半径 (Å)
    pub rvdw: f32,
    /// 标准原子量。
    ///
    /// **f64 而不是 f32**:上游表里写的是 `12.011`,存成 f32 再读回来是
    /// 12.0109996…,差 3e-7。这点差在化学上无关紧要,可它让"与参照逐位相同"
    /// 变成"与参照差在末几位",判据只好配一个容差 —— 容差一旦有了,就再也
    /// 分不出"存储精度"和"抄错了一位"。电负性同理。
    pub mass: f64,
    /// 外层电子数
    pub outer_electrons: u8,
    /// 最常见同位素的质量数
    pub common_isotope: u16,
    /// 最常见同位素的精确质量
    pub common_isotope_mass: f64,
    /// Pauling 电负性。
    ///
    /// `None` 表示**该元素没有公认的 Pauling 值**(稀有气体 He/Ne/Ar/Rn、
    /// Pm、Eu、Tb、Yb、Fr 等),不是"还没填"。两者必须能被下游分开:
    /// 拿一个默认数顶上去,调用方就再也看不出这一格是量出来的还是补的。
    pub electronegativity: Option<f64>,
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

/// 某个同位素的精确质量。该元素/质量数不在表里时返回 `None`。
///
/// 与"标准原子量"([`Element::mass`])是两件事:后者是天然丰度加权的平均值,
/// 与写不写同位素标注无关;这里给的是**指定那一个核素**的质量。
/// 氘(`[2H]`)的标准原子量仍是 1.008,精确质量是 2.0141 —— 差了一倍,
/// 混用会在氘代化合物上悄悄给出错的数。
///
/// 表由 `harness/gen_elements.py` 从 RDKit `atomic_data.cpp` 的
/// `isotopesAtomData` 抽取(3111 条,覆盖 113 种元素;105–109 号在源表里
/// 本来就没有数据)。
#[must_use]
pub fn isotope_mass(atomic_num: u8, mass_number: u16) -> Option<f64> {
    ISOTOPE_MASSES
        .get(atomic_num as usize)?
        .iter()
        .find(|&&(m, _)| m == mass_number)
        .map(|&(_, mass)| mass)
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
/// # 六族的 +1 也算:锍盐、硒盐、碲盐
///
/// R₃S⁺ 是三配体 + 一对孤对(S⁺ 有 5 个价电子,三根键用掉 3 个,剩 2 个正好一对),
/// 构型稳定、可以拆分。真正没有孤对、四配体全在的是季铵 R₄N⁺。
///
/// 这一档先前被整个排除,理由写的是"外部判据看不见":当时举的例子是
/// `C[S@+](C)CC` —— 而那个分子**两个甲基一模一样,本来就不是手性中心**,
/// 任何实现都会把标记清掉。换成真正的锍盐 `C[S@+](CC)CCC`,钉住的 RDKit 2025.09.2
/// 给的是 `CHI_TETRAHEDRAL_CCW`。**拿一个非手性的例子论证"判据看不见",
/// 论证的是别的事。**
///
/// 代价是实打实的:排除期间 `C[S@+](CC)CCC` 走一趟二维往返构型整个丢掉
/// (`C[S@@+](CCC)CC` → `C[S+](CCC)CC`),不报错。
///
/// **五族(P / As)的 +1 仍然不算**,而且不是保守取舍:P⁺ 只有 4 个价电子,
/// 三根键用掉 3 个,剩下的是**一个单电子**而不是一对 —— 那是膦自由基阳离子,
/// 不是稳定的立体中心。四配位的鏻盐 R₄P⁺ 本来就走四邻居那条路。
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
    match atomic_num {
        // 六族:中性(亚砜等)与 +1(锍盐、硒盐、碲盐)都是三配体 + 一对孤对
        16 | 34 | 52 => formal_charge <= 1,
        // 五族:中性的膦、胂有孤对;+1 之后剩的是单电子,不是一对
        15 | 33 => formal_charge <= 0,
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn 孤对立体中心认哪几个元素与电荷() {
        for (z, chg, want) in [
            (16u8, 0i8, true), // S:亚砜
            (15, 0, true),     // P:膦
            (33, 0, true),     // As
            (34, 0, true),     // Se
            (52, 0, true),     // Te
            // 六族的 +1 也算:锍盐 R₃S⁺ 是三配体 + 一对孤对,钉住的外部实现
            // 对**真正**的锍盐(三个取代基不全同)给的是四面体标记。
            (16, 1, true), // [S+]:锍盐
            (34, 1, true), // [Se+]:硒盐
            (52, 1, true), // [Te+]:碲盐
            // 五族的 +1 不算,而且不是保守取舍:P⁺ 三配位剩的是单电子不是一对。
            (15, 1, false), // [P+]
            (33, 1, false), // [As+]
            // 再往上就是四配位了,不走这一支
            (16, 2, false),
            (7, 0, false), // N:孤对翻转太快,不当立体中心
            (6, 0, false), // C:三配位是 sp²
            (8, 0, false), // O:三配位是 [O+]
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

    /// 电负性表的护栏。
    ///
    /// 表是从 RDKit 源码抽的,而 CI 里没有那份源码 —— 重新生成一遍再逐字节比
    /// 的闸进不了 CI,进不了 CI 就不是闸。所以这里钉两样在 CI 里跑得动的:
    ///
    /// 1. **几个课本上的值**(氟最大、铯最小、碳氧氢),抓整表错位与量纲写错;
    /// 2. **哪些元素没有值**,以及有值的元素**个数**。
    ///
    /// 第 2 条是关键:没有公认 Pauling 值的元素(稀有气体、Pm/Eu/Tb/Yb/Fr)
    /// 必须是 `None`。谁要是给它们补个"默认 2.0",这里立刻红 —— 那种补法会让
    /// 下游再也分不出"这个元素没有值"和"这个元素的值恰好是 2.0"。
    #[test]
    fn pauling_electronegativity_table() {
        let en = |sym: &str| by_symbol(sym).unwrap().electronegativity;
        // 课本值:氟最大 3.98,铯最小 0.79
        assert_eq!(en("F"), Some(3.98));
        assert_eq!(en("Cs"), Some(0.79));
        assert_eq!(en("H"), Some(2.2));
        assert_eq!(en("C"), Some(2.55));
        assert_eq!(en("N"), Some(3.04));
        assert_eq!(en("O"), Some(3.44));
        // 全表的最大值就该是氟
        let max = ELEMENTS
            .iter()
            .filter_map(|e| e.electronegativity)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(
            (max - 3.98).abs() < f64::EPSILON,
            "全表最大电负性应为氟的 3.98,实得 {max}"
        );

        // 没有公认值的:如实为 None
        for sym in ["He", "Ne", "Ar", "Rn", "Pm", "Eu", "Tb", "Yb", "Fr"] {
            assert_eq!(en(sym), None, "{sym} 没有公认的 Pauling 值,应为 None");
        }
        // 通配原子当然也没有
        assert_eq!(by_atomic_num(0).unwrap().electronegativity, None);

        let n = ELEMENTS
            .iter()
            .filter(|e| e.electronegativity.is_some())
            .count();
        assert_eq!(
            n, 93,
            "有电负性的元素个数变了 —— 表被改过,重跑 gen_elements.py"
        );
    }

    /// 同位素质量表的护栏。CI 里没有 RDKit 源码,所以钉几个核素 + 一条
    /// "标准原子量与精确质量不是一回事"的对照。
    #[test]
    fn isotope_mass_table() {
        // 氕与氘:标准原子量都是 1.008,精确质量差了一倍
        assert_eq!(by_symbol("H").unwrap().mass, 1.008);
        assert_eq!(isotope_mass(1, 1), Some(1.007_825_032));
        assert_eq!(isotope_mass(1, 2), Some(2.014_101_778));
        assert_eq!(isotope_mass(6, 12), Some(12.0));
        assert_eq!(isotope_mass(6, 13), Some(13.00335484));
        assert_eq!(isotope_mass(7, 15), Some(15.0001089));
        // 表里没有的质量数、以及源表本来就没有数据的元素(105–109)
        assert_eq!(isotope_mass(6, 200), None);
        assert_eq!(isotope_mass(105, 268), None);
        // 每个有数据的元素都得含它最常见的那个同位素 —— 分块解析漏首行时,
        // 漏掉的恰恰是排在最前的质量数(实测漏过 H-1)
        for e in ELEMENTS.iter() {
            let rows = ISOTOPE_MASSES[e.atomic_num as usize];
            if rows.is_empty() || e.common_isotope == 0 {
                continue;
            }
            assert!(
                isotope_mass(e.atomic_num, e.common_isotope).is_some(),
                "{} 的最常见同位素 {} 不在表里",
                e.symbol,
                e.common_isotope
            );
        }
        let n: usize = ISOTOPE_MASSES.iter().map(|r| r.len()).sum();
        assert_eq!(n, 3111, "同位素条目数变了 —— 表被改过,重跑 gen_elements.py");
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

    /// 多价元素的默认价表:**期望值写死,不许引用被测的那张表**。
    ///
    /// 先前这里写的是 `assert_eq!(s.default_valence_for(3), Some(s.valences[1]))`
    /// —— 表一改两边一起动。把 S 的 `[2,4,6]` 变异成 `[2,5,6]` 它照样绿,而那个
    /// 变异会让 `[SH](=O)C` 的隐式氢从 1 变 2。
    ///
    /// 而且先前只钉了六个**单值**列表,P `[3,5]` / S `[2,4,6]` / I `[1,3,5]`
    /// 一个没钉 —— 多值那档才是抄错后果最重的。
    #[test]
    fn default_valence_selection() {
        for (sym, want) in [
            ("S", &[2, 4, 6][..]),
            ("P", &[3, 5][..]),
            ("I", &[1, 3, 5][..]),
            ("Cl", &[1][..]),
            ("N", &[3][..]),
            ("C", &[4][..]),
        ] {
            let e = by_symbol(sym).unwrap();
            assert_eq!(e.valences, want, "{sym} 的默认价表");
        }

        let s = by_symbol("S").unwrap();
        // 选的是第一个 ≥ 已用价的那一档
        assert_eq!(s.default_valence_for(2), Some(2));
        assert_eq!(s.default_valence_for(3), Some(4));
        assert_eq!(s.default_valence_for(5), Some(6));
        assert_eq!(s.default_valence_for(7), None, "超过最大价就没有可用的了");

        let c = by_symbol("C").unwrap();
        assert_eq!(c.default_valence_for(4), Some(4));
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
