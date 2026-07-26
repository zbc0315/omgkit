//! SMARTS 的原子/键表达式:查询树与它的解析。
//!
//! # 三档优先级,不能想当然
//!
//! 从高到低是 `&`(可省略)> `,` > `;`,`!` 是一元且绑得最紧:
//!
//! | 写法 | 含义 |
//! |---|---|
//! | `[C,N;H1]` | (C 或 N) **且** H1 |
//! | `[C,N&H1]` | C **或** (N 且 H1) |
//! | `[C,NH1]` | 同上 —— **并置就是 `&`**,不是"顺序连接" |
//!
//! 前两行语义完全不同,而且都是常见写法。第三行最容易错:`NH1` 看起来像一个
//! 整体,实际是 `N & H1`,于是整条表达式变成 `C , (N & H1)`。
//!
//! 用优先级爬升(Pratt)处理最干净 —— 每一档一个函数,递归下降下去,
//! 没有需要维护的运算符栈。
//!
//! # 一元 `!` 只作用于紧跟其后的那一项
//!
//! `[!C;!N]` 是"既非脂肪碳、又非脂肪氮",不是 `!(C;!N)`。

use omgkit_core::ChiralTag;

/// 原子的一个查询基元。
///
/// 每个基元只回答一个问题。复合条件由 [`AtomExpr`] 的逻辑运算组合出来,
/// 基元本身不带否定 —— 否定统一走 [`AtomExpr::Not`],这样"!" 的作用范围
/// 只有一处定义,不会出现两套否定语义互相打架。
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AtomPrim {
    /// `*` —— 任意原子
    Any,
    /// `a` —— 芳香原子
    Aromatic,
    /// `A` —— 脂肪原子
    Aliphatic,
    /// 元素。`aromatic` 为 `None` 时不限芳香性(`#6` 的写法),
    /// 否则区分 `C`(脂肪碳)与 `c`(芳香碳)。
    Element {
        /// 原子序数
        z: u8,
        /// `Some(true)` = 必须芳香,`Some(false)` = 必须脂肪
        aromatic: Option<bool>,
    },
    /// `D<n>` —— 显式连接数(不含隐式氢)
    Degree(u32),
    /// `X<n>` —— 总连接数(含氢)
    TotalDegree(u32),
    /// `H<n>` —— 总氢数
    TotalHs(u32),
    /// `h<n>` —— 隐式氢数
    ImplicitHs(u32),
    /// `R<n>` —— 所属环的个数;`R` 不带数字表示"在任意环中"
    RingCount(Option<u32>),
    /// `r<n>` —— 最小环的大小;`r` 不带数字等价于 `R`
    RingSize(Option<u32>),
    /// `x<n>` —— 环键数;`x` 不带数字表示"至少有一条环键"
    RingBondCount(Option<u32>),
    /// `v<n>` —— 总价
    Valence(u32),
    /// `+` / `-` / `+2` —— 形式电荷
    Charge(i32),
    /// 前置数字 —— 同位素
    Isotope(u16),
    /// `:n` —— 反应原子映射号
    AtomMap(u16),
    /// `@` / `@@` —— 手性
    Chirality(ChiralTag),
    /// `$(...)` —— 递归 SMARTS:该原子要能作为括号内模式的**首原子**匹配上
    Recursive(Box<super::QueryMol>),
}

/// 原子的查询表达式树。
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AtomExpr {
    /// 单个基元
    Prim(AtomPrim),
    /// `!`
    Not(Box<AtomExpr>),
    /// `&` 或并置,也用于 `;`
    And(Vec<AtomExpr>),
    /// `,`
    Or(Vec<AtomExpr>),
}

/// 键的查询基元。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BondPrim {
    /// `~` —— 任意键
    Any,
    /// `-`
    Single,
    /// `=`
    Double,
    /// `#`
    Triple,
    /// `$`
    Quadruple,
    /// `:`
    Aromatic,
    /// `@` —— 环键
    InRing,
    /// `/`
    UpRight,
    /// `\`
    DownRight,
    /// `->` —— 配位键,电子对由**左**端提供
    Dative,
    /// `<-` —— 配位键,电子对由**右**端提供。
    ///
    /// 与 [`Dative`](Self::Dative) 是两个基元而不是一个带方向的标志:
    /// 键表达式匹配时不知道自己将被摆在哪个朝向上,方向必须留在基元里。
    DativeReversed,
}

/// 键的查询表达式树。
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BondExpr {
    /// 单个基元
    Prim(BondPrim),
    /// `!`
    Not(Box<BondExpr>),
    /// `&` 或并置,也用于 `;`
    And(Vec<BondExpr>),
    /// `,`
    Or(Vec<BondExpr>),
}

impl BondExpr {
    /// 未写键符号时的默认:单键或芳香键。
    ///
    /// 这是 SMARTS 与 SMILES 的一处关键差异。SMILES 里"没写键符号"要看两端
    /// 原子是否芳香来定;SMARTS 里两端可能是通配,当场定不下来,所以默认值
    /// 本身就是一个**析取式**,留到匹配时再判。
    #[must_use]
    pub fn default_bond() -> Self {
        Self::Or(vec![
            Self::Prim(BondPrim::Single),
            Self::Prim(BondPrim::Aromatic),
        ])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 表达式树是可比较的 —— 解析测试全靠这个,先确认它没退化成
    /// "只要形状像就相等"。
    #[test]
    fn expressions_compare_structurally() {
        let c = AtomExpr::Prim(AtomPrim::Element {
            z: 6,
            aromatic: Some(false),
        });
        let n = AtomExpr::Prim(AtomPrim::Element {
            z: 7,
            aromatic: Some(false),
        });
        assert_eq!(c.clone(), c.clone());
        assert_ne!(c.clone(), n.clone());
        assert_ne!(
            AtomExpr::And(vec![c.clone(), n.clone()]),
            AtomExpr::Or(vec![c.clone(), n.clone()]),
        );
        assert_ne!(
            AtomExpr::And(vec![c.clone(), n.clone()]),
            AtomExpr::And(vec![n, c]),
            "顺序不同即不同 —— 归一化是另一回事,不能在这里悄悄发生"
        );
    }

    #[test]
    fn default_bond_is_single_or_aromatic() {
        assert_eq!(
            BondExpr::default_bond(),
            BondExpr::Or(vec![
                BondExpr::Prim(BondPrim::Single),
                BondExpr::Prim(BondPrim::Aromatic),
            ])
        );
    }
}
