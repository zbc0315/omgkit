//! # omgkit-match
//!
//! omgkit 的 L5 层:子结构匹配。
//!
//! 给一个 SMARTS 查询和一个净化过的分子,找出所有能把查询原子映射到分子原子、
//! 且保持连接关系与全部查询条件的映射。
//!
//! ```no_run
//! use omgkit_match::{substructure_matches, MolProps};
//!
//! let mut mol = omgkit_io::smiles::parse("c1ccccc1O").unwrap();
//! // ... 净化 ...
//! let props = MolProps::compute(&mol);
//! let query = omgkit_io::smarts::parse("c[OH]").unwrap();
//! let hits = substructure_matches(&query, &mol, &props, Default::default());
//! ```
//!
//! ## 为什么必须自己写
//!
//! 通用图同构库解不了这里的问题:查询原子的条件是一棵**表达式树**,不是一个
//! 可比较的标签。而且化学特化(按元素稀有度排序、按度剪枝、环成员位图预筛)
//! 恰恰是这里最大的性能来源 —— 通用库拿不到这些信息。

pub mod matcher;
pub mod props;
pub mod react;

pub use matcher::{substructure_matches, MatchOptions};
pub use props::MolProps;
pub use react::{run_on_substrate, run_reactants, Outcome, ProductSet};
