//! # omgkit-core
//!
//! omgkit 的 L0 层:列式分子批、CSR 图与元素周期表。
//!
//! 本 crate 刻意保持**零依赖**,并且不含任何化学语义 —— 没有净化、没有
//! 芳香性感知、没有价键推断。它只回答"分子的拓扑与原子属性如何存放"这一个问题。
//! 化学语义在 `omgkit-chem`(L2),解析在 `omgkit-io`(L1)。
//!
//! ## 两种表示
//!
//! | 类型 | 可变性 | 布局 | 用途 |
//! |---|---|---|---|
//! | [`MolBuilder`] | 可变 | AoS,逐分子 | 建图:解析、反应产物构建 |
//! | [`MolBatch`] | 不可变 | SoA + CSR,跨分子 | 算法、并行、零拷贝导出 |
//!
//! 二者经 [`MolBatchBuilder::push`] / [`MolView::to_builder`] 互转,
//! 往返恒等由测试保证。
//!
//! ```
//! use omgkit_core::{BondOrder, MolBatchBuilder, MolBuilder};
//!
//! let mut ethanol = MolBuilder::new();
//! let c0 = ethanol.add_atom(6);
//! let c1 = ethanol.add_atom(6);
//! let o = ethanol.add_atom(8);
//! ethanol.add_bond(c0, c1, BondOrder::Single).unwrap();
//! ethanol.add_bond(c1, o, BondOrder::Single).unwrap();
//!
//! let mut bb = MolBatchBuilder::new();
//! bb.push(&ethanol).unwrap();
//! bb.push(&ethanol).unwrap();
//! let batch = bb.finish();
//!
//! assert_eq!(batch.num_mols(), 2);
//! assert_eq!(batch.num_atoms(), 6);
//!
//! // 列是连续内存,可零拷贝暴露给 numpy / Arrow
//! assert_eq!(batch.atomic_nums(), &[6, 6, 8, 6, 6, 8]);
//!
//! // 单分子是零拷贝视图,下标自动换算为局部
//! let m = batch.mol(1).unwrap();
//! assert_eq!(m.num_atoms(), 3);
//! assert_eq!(m.degree(1), 2);
//! ```

pub mod batch;
pub mod builder;
pub mod element;
pub mod error;
pub mod types;
pub mod valence;
pub mod view;

mod element_data;

pub use batch::{MolBatch, MolBatchBuilder};
pub use builder::{AtomData, BondData, BondMut, MolBuilder, Neighbors};
pub use element::Element;
pub use error::{Error, Result};
pub use types::{
    square_planar_renumber, AtomFlags, BondDirection, BondFlags, BondOrder, BondStereo, ChiralTag,
    Hybridization, SQUARE_PLANAR_TRANS,
};
pub use view::MolView;
