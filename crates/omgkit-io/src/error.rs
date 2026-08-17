//! 解析错误 —— 带精确位置。
//!
//! 一条只说"语法错误"的报错,在几十个字符的 SMILES 上就已经很难定位了。
//! 精确到列的位置是手写解析器相对生成式解析器的主要收益之一,所以错误类型
//! 一开始就带位置和可打印的插字号视图。

use core::fmt;

/// 解析失败的原因。
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum ParseErrorKind {
    /// 空输入
    Empty,
    /// 出现了不该出现的字符
    UnexpectedChar(char),
    /// 输入提前结束(如方括号未闭合)
    UnexpectedEnd,
    /// 元素符号无法识别
    UnknownElement(String),
    /// 方括号原子语法错误
    BadBracketAtom(&'static str),
    /// 分支括号不匹配
    UnbalancedParen,
    /// 分支为空,如 `C()`
    EmptyBranch,
    /// 环闭合标号始终未配对
    UnclosedRingBond(u32),
    /// 环闭合两端是同一个原子,如 `C11`
    RingBondToSelf(u32),
    /// 同一对原子间重复的环闭合键
    DuplicateRingBond(u32),
    /// 环闭合两端指定了互相冲突的键级
    ConflictingRingBondOrder(u32),
    /// 键符号后面没有原子
    DanglingBond,
    /// 数字超出可表示范围
    NumberOverflow,
    /// 立体标记的排列序号超出该几何的取值范围,如 `@TB21`
    StereoPermOutOfRange {
        /// 几何的书写形式(`TH` / `AL` / `SP` / `TB` / `OH`)
        geometry: &'static str,
        /// 写出来的序号
        got: u32,
        /// 该几何允许的最大序号
        max: u32,
    },
}

impl fmt::Display for ParseErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "空的 SMILES"),
            Self::UnexpectedChar(c) => write!(f, "意外的字符 {c:?}"),
            Self::UnexpectedEnd => write!(f, "输入意外结束"),
            Self::UnknownElement(s) => write!(f, "无法识别的元素符号 {s:?}"),
            Self::BadBracketAtom(why) => write!(f, "方括号原子语法错误:{why}"),
            Self::UnbalancedParen => write!(f, "括号不匹配"),
            Self::EmptyBranch => write!(f, "空的分支"),
            Self::UnclosedRingBond(n) => write!(f, "环闭合标号 {n} 未配对"),
            Self::RingBondToSelf(n) => write!(f, "环闭合标号 {n} 的两端是同一个原子"),
            Self::DuplicateRingBond(n) => write!(f, "环闭合标号 {n} 在同一对原子间重复"),
            Self::ConflictingRingBondOrder(n) => {
                write!(f, "环闭合标号 {n} 的两端指定了冲突的键级")
            }
            Self::DanglingBond => write!(f, "键符号后面缺少原子"),
            Self::NumberOverflow => write!(f, "数值超出范围"),
            Self::StereoPermOutOfRange { geometry, got, max } => write!(
                f,
                "立体标记 @{geometry}{got} 的序号超出范围,@{geometry} 最大为 {max}"
            ),
        }
    }
}

/// 一次解析失败,携带出错位置。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseError {
    /// 失败原因
    pub kind: ParseErrorKind,
    /// 出错处的字节偏移(0 起)
    pub pos: usize,
    /// 出错时的完整输入,用于渲染上下文
    pub input: String,
}

impl ParseError {
    pub(crate) fn new(kind: ParseErrorKind, pos: usize, input: &[u8]) -> Self {
        Self {
            kind,
            pos,
            input: String::from_utf8_lossy(input).into_owned(),
        }
    }

    /// 渲染成带插字号的两行视图:
    ///
    /// ```text
    /// CC(C)C1CC
    ///         ^ 环闭合标号 1 未配对
    /// ```
    #[must_use]
    pub fn render(&self) -> String {
        // 插字号按字符数而非字节数缩进,含中日韩字符时才不会错位
        let cols = self.input[..self.pos.min(self.input.len())].chars().count();
        format!("{}\n{}^ {}", self.input, " ".repeat(cols), self.kind)
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "第 {} 字符处:{}", self.pos, self.kind)
    }
}

impl std::error::Error for ParseError {}

/// 解析结果。
pub type Result<T> = core::result::Result<T, ParseError>;
