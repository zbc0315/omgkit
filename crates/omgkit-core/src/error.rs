//! 核心层错误类型。
//!
//! L0 只关心结构性错误(下标越界、自环、规模溢出)。化学语义错误
//! (价键异常、芳香性感知失败)属于 L2,由 `omgkit-chem` 自行定义。

use core::fmt;

/// 核心层错误。
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum Error {
    /// 原子下标越界
    AtomIndexOutOfRange {
        /// 越界的下标
        index: u32,
        /// 实际原子数
        num_atoms: u32,
    },
    /// 试图添加自环。分子图中不存在自环。
    SelfLoop {
        /// 两端相同的那个原子
        atom: u32,
    },
    /// 分子下标越界
    MolIndexOutOfRange {
        /// 越界的下标
        index: u32,
        /// 实际分子数
        num_mols: u32,
    },
    /// 键下标越界
    BondIndexOutOfRange {
        /// 越界的下标
        index: u32,
        /// 实际键数
        num_bonds: u32,
    },
    /// 批规模超出 `u32` 索引上限
    BatchTooLarge {
        /// 溢出的是哪一类实体("atoms" / "bonds" / "molecules")
        what: &'static str,
        /// 实际数量
        count: usize,
    },
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::AtomIndexOutOfRange { index, num_atoms } => {
                write!(f, "原子下标 {index} 越界(共 {num_atoms} 个原子)")
            }
            Self::SelfLoop { atom } => {
                write!(f, "原子 {atom} 上出现自环;分子图中不允许自环")
            }
            Self::BondIndexOutOfRange { index, num_bonds } => {
                write!(f, "键下标 {index} 越界(共 {num_bonds} 条键)")
            }
            Self::MolIndexOutOfRange { index, num_mols } => {
                write!(f, "分子下标 {index} 越界(共 {num_mols} 个分子)")
            }
            Self::BatchTooLarge { what, count } => {
                write!(f, "批中 {what} 数量 {count} 超出 u32 索引上限")
            }
        }
    }
}

impl std::error::Error for Error {}

/// 核心层的 `Result` 别名。
pub type Result<T> = core::result::Result<T, Error>;
