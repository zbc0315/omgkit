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
//! ## 净化之外:两个**读**它结果的模块
//!
//! 本 crate 里还有两个模块**不是净化的步骤**,而是净化输出的消费者。放在这里
//! 是因为它们要读第 3、4、7、8、9 步填的东西(隐式氢、环、芳香、共轭、杂化),
//! 拆出去就成了一个只依赖 `omgkit-chem` 的空壳 crate。
//!
//! | 模块 | 是什么 |
//! |---|---|
//! | [`gasteiger`] | Gasteiger–Marsili 部分电荷(PEOE 迭代) |
//! | [`descriptors`] | 图神经网络特征化要读的那 16 个原子/键描述符,汇到一处交付 |
//!
//! 它们**不会**被 [`pipeline::sanitize`] 调用 —— 净化不该顺手算一批
//! 没人要的数。反过来,它们要求分子已经净化过:没净化不报错,只会让每一项都
//! 是解析时的占位值。
//!
//! [`descriptors`] 还有一处前置在**本 crate 之外**:双键顺反由
//! `omgkit_io::stereo::perceive_bond_stereo` 填,而 `omgkit-io` 是本 crate 的
//! 同级,调不到。详见该模块的文档。
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
pub mod descriptors;
pub mod explicit_hs;
pub mod gasteiger;
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
pub use descriptors::{atom_descriptors, bond_descriptors, AtomDescriptors, BondDescriptors};
pub use explicit_hs::add_explicit_hs;
pub use gasteiger::gasteiger_charges;
pub use kekulize::{kekulize, KekulizeError};
pub use organometallics::{cleanup_organometallics, is_metal};
pub use pipeline::{sanitize, SanitizeError};
pub use radicals::assign_radicals;
pub use remove_hs::{is_removable, remove_hs};
pub use rings::{fused_ring_systems, perceive_rings, RingPerception};
pub use sssr::{ring_set, Ring};
pub use valence::{
    update_property_cache, valence_shift, ValenceError, ValenceErrorKind, ValenceResult,
};
