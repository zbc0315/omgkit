//! # omgkit-chem
//!
//! omgkit 的 L2 层:化学语义 —— 净化管线、芳香性、价键与隐式氢。
//!
//! ## 净化管线的 12 个步骤
//!
//! 步骤之间有依赖,**不可乱序**。"实现在"一列为空的是尚未实现的:
//!
//! | # | 步骤 | 实现在 |
//! |---|---|---|
//! | 1 | 修正硝基 / 叠氮 / 磷酰 / 高卤酸的非标准画法 | [`cleanup`] |
//! | 2 | 有机金属键改配位键 | [`organometallics`] |
//! | 3 | 价键计算 + 隐式氢推断 | [`valence`] |
//! | 4 | 环感知 | [`rings`] |
//! | 5 | 芳香环 → 交替单双键 | [`mod@kekulize`] |
//! | 6 | 自由基电子数 | [`radicals`] |
//! | 7 | 芳香性感知 | [`aromaticity`] |
//! | 8 | 共轭标记 | [`conjugation`] |
//! | 9 | 杂化标记 | [`conjugation`] |
//! | 10 | 阻转异构 | 未实现 |
//! | 11 | 剔除几何上不成立的立体标记 | [`mod@cleanup_chirality`] |
//! | 12 | 显式/隐式氢调整 | [`mod@adjust_hs`] |
//!
//! ## `removeHs` 不属于净化
//!
//! 删除显式 `[H]` 原子会改变原子数,是一个独立的编辑操作,不在这 12 步里。
//! 实现在 [`mod@remove_hs`],由调用方按需显式调用。
//!
//! ## 净化基本不改动连通性
//!
//! 只有第 2 步会把某些键改成配位键并**交换端点**;除此之外净化只改属性,
//! 不增删原子和键。这是"在解析结果的图上直接做环感知"成立的前提。

pub mod adjust_hs;
pub mod aromaticity;
pub mod cleanup;
pub mod cleanup_chirality;
pub mod conjugation;
pub mod kekulize;
pub mod organometallics;
pub mod pipeline;
pub mod radicals;
pub mod remove_hs;
pub mod rings;
pub mod sssr;
pub mod valence;

pub use adjust_hs::adjust_hs;
pub use aromaticity::set_aromaticity;
pub use cleanup::clean_up;
pub use cleanup_chirality::cleanup_chirality;
pub use conjugation::{mark_conjugated_atoms, set_conjugation, set_hybridization};
pub use kekulize::{kekulize, KekulizeError};
pub use organometallics::{cleanup_organometallics, is_metal};
pub use pipeline::{sanitize, SanitizeError};
pub use radicals::assign_radicals;
pub use remove_hs::{is_removable, remove_hs};
pub use rings::{fused_ring_systems, perceive_rings, RingPerception};
pub use sssr::{ring_set, Ring};
pub use valence::{update_property_cache, ValenceError, ValenceErrorKind, ValenceResult};
