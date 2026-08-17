//! # omgkit-io
//!
//! omgkit 的 L1 层:分子表示格式的解析与输出。
//!
//! 当前是 SMILES 的解析与写出。本层**只做语法**,不含化学语义 —— 芳香性
//! 感知、价键计算、隐式氢推断、环感知全部属于 L2(`omgkit-chem`)。
//!
//! ```
//! use omgkit_io::smiles;
//!
//! let mol = smiles::parse("CC(=O)O").unwrap();
//! assert_eq!(mol.num_atoms(), 4);
//! assert_eq!(mol.num_bonds(), 3);
//!
//! // 写回去
//! assert_eq!(smiles::write(&mol).smiles, "CC(=O)O");
//!
//! // 错误带精确位置
//! let err = smiles::parse("C1CC").unwrap_err();
//! assert_eq!(err.pos, 1);
//! println!("{}", err.render());
//! // C1CC
//! //  ^ 环闭合标号 1 未配对
//! ```

pub mod canon;
pub mod error;
pub mod smarts;
pub mod smiles;
pub mod stereo;

pub use error::{ParseError, ParseErrorKind, Result};
