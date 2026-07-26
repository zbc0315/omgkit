//! SMARTS —— 子结构查询模式的解析。
//!
//! SMARTS 在语法上是 SMILES 的超集:原子和键的位置上可以写**表达式**而不只是
//! 具体的元素与键级。所以拓扑部分(分支、环闭合、片段)与 SMILES 完全同构,
//! 差别集中在方括号里。
//!
//! # 查询是外挂的,不塞进 `MolBuilder`
//!
//! [`QueryMol`] = 一个只承载**拓扑**的 [`MolBuilder`] + 逐原子逐键的查询树。
//! 分子的数据结构因此不必为查询这件事付出任何代价 —— 列式存储里不会多出
//! 一列指针,也不会有"这个字段在查询分子里含义不同"的暗坑。
//!
//! 拓扑侧的原子一律置 [`AtomFlags::HAS_QUERY`](omgkit_core::AtomFlags::HAS_QUERY),键置
//! [`BondFlags::HAS_QUERY`](omgkit_core::BondFlags::HAS_QUERY),表示"这里的值没有意义,要看查询树"。

mod bond;
mod eval;
mod expr;
pub(super) mod mol;
mod parse;
mod reaction;
mod write;

pub use expr::{AtomExpr, AtomPrim, BondExpr, BondPrim};

pub use bond::{parse_bond_expr, starts_bond_expr};
pub use eval::{allowed_elements, atom_matches, bond_matches, AtomProps, BondProps};
pub use mol::parse;
pub use parse::parse_atom_expr;
pub use reaction::{parse_reaction, Reaction};
pub use write::{atom_expr_string, bond_expr_string, write, write_reaction};

use omgkit_core::MolBuilder;

/// 一个查询模式:拓扑 + 逐原子逐键的查询树。
///
/// 相等性比的是**表示**:同样的原子序列、同样的键序列、同样的查询树。
/// 它不回答"两个模式是否等价"——那要判逻辑等价加图同构,是另一回事。
/// 这个区分刻意留在类型之外:`MolBuilder` 本身不实现 `PartialEq`,免得
/// `a == b` 被读成"同一个分子"。
#[derive(Debug, Clone)]
pub struct QueryMol {
    /// 只承载拓扑。原子的元素、键的键级等字段**没有意义**,
    /// 语义全在查询树里。
    pub topology: MolBuilder,
    /// 与 `topology.atoms()` 一一对应
    pub atoms: Vec<AtomExpr>,
    /// 与 `topology.bonds()` 一一对应
    pub bonds: Vec<BondExpr>,
}

/// 取一个原子查询要求的手性。没写手性时返回 `None`。
///
/// # 为什么要单独取出来
///
/// 手性**不能逐原子判定** —— 标记相对各自分子的邻居存储顺序,要知道查询的
/// 邻居分别映到底物的哪些原子才算得出宇称。所以求值时一律放行,由匹配器在
/// 映射完成后拿这个函数取回要求再校验。
#[must_use]
pub fn required_chirality(expr: &AtomExpr) -> Option<omgkit_core::ChiralTag> {
    fn walk(e: &AtomExpr) -> Option<omgkit_core::ChiralTag> {
        match e {
            AtomExpr::Prim(AtomPrim::Chirality(t)) => Some(*t),
            AtomExpr::And(parts) => parts.iter().find_map(walk),
            // `,` 下面的手性要求彼此矛盾,取任一个都不对;`!` 同理。
            // 这类写法罕见,放行比猜一个好。
            AtomExpr::Or(_) | AtomExpr::Not(_) => None,
            AtomExpr::Prim(_) => None,
        }
    }
    walk(expr)
}

/// 取一个原子查询里写死的映射号。没有 `:n` 时返回 `None`。
///
/// 映射号是反应里连接反应物与产物的**唯一**纽带,所以要能从查询树里取回来。
/// 它一定出现在合取的顶层(`[C:1]` 解析成 `C & :1`),或者整个表达式就是它。
#[must_use]
pub fn map_number(expr: &AtomExpr) -> Option<u16> {
    match expr {
        AtomExpr::Prim(AtomPrim::AtomMap(n)) => Some(*n),
        AtomExpr::And(parts) => parts.iter().find_map(map_number),
        // 析取式里的映射号没有确定含义,不取
        _ => None,
    }
}

impl PartialEq for QueryMol {
    fn eq(&self, other: &Self) -> bool {
        self.atoms == other.atoms
            && self.bonds == other.bonds
            && self.topology.atoms() == other.topology.atoms()
            && self.topology.bonds() == other.topology.bonds()
    }
}

impl Eq for QueryMol {}

impl QueryMol {
    /// 原子数
    #[must_use]
    pub fn num_atoms(&self) -> usize {
        self.atoms.len()
    }

    /// 键数
    #[must_use]
    pub fn num_bonds(&self) -> usize {
        self.bonds.len()
    }

    /// 查询树与拓扑是否一一对应。
    ///
    /// 两者是分开维护的两个数组,长度一旦对不上,后续的匹配就会按错位的
    /// 查询去判断原子 —— 结果看起来"只是匹配不太对",极难定位到根因。
    #[must_use]
    pub fn is_consistent(&self) -> bool {
        self.atoms.len() == self.topology.num_atoms()
            && self.bonds.len() == self.topology.num_bonds()
    }
}
