//! 反应 SMARTS 的解析:`反应物 > 试剂 > 产物`。
//!
//! 每一段用 `.` 分成若干**模板**,每个模板是一条 [`QueryMol`]。
//!
//! # 三个切分陷阱
//!
//! **1. `>` 不总是分隔符。** 配位键写作 `->`,`[N:1]->[Cu:2]>>...` 里的第一个
//! `>` 属于键,不属于分段。判据是"前一个字符是不是 `-`"。
//!
//! **2. 括号可能是分组也可能是分支。** `(A.B)>>C` 的括号把两个片段**捆成一个
//! 模板**(1 个反应物而不是 2 个);而 `C(C)C` 的括号是分支。区分方法:分组括号
//! 从模板的第一个字符开始,且与它配对的 `)` 正好是模板的最后一个字符。
//!
//! **3. `.` 和 `>` 都可能出现在方括号里。** 递归 SMARTS `[$(C.C)]` 里的 `.`
//! 不是模板分隔符。所以切分要跟踪方括号与圆括号的深度。
//!
//! # 映射号是反应物与产物之间**唯一**的联系
//!
//! 反应物模板里带 `:n` 的原子,与产物模板里同号的原子是同一个原子。没有映射号
//! 的原子:出现在反应物侧表示"匹配到但**删掉**",出现在产物侧表示"**新建**"。
//!
//! 这条语义决定了产物生成的全部行为,所以映射号的一致性要在解析时就查
//! —— 见 [`Reaction::map_numbers`]。

use std::collections::{BTreeMap, BTreeSet};

use super::QueryMol;
use crate::error::{ParseError, ParseErrorKind as K, Result};

/// 一条反应模式。
#[derive(Debug, Clone)]
pub struct Reaction {
    /// 反应物模板
    pub reactants: Vec<QueryMol>,
    /// 试剂模板(中间那一段)。不参与产物构建,只用于筛选。
    pub agents: Vec<QueryMol>,
    /// 产物模板
    pub products: Vec<QueryMol>,
}

impl Reaction {
    /// 反应物侧与产物侧各自出现的映射号。
    ///
    /// 两侧的差集就是"会被删掉的"和"会被新建的"原子 —— 那是语义,不是错误,
    /// 所以这里只报告,不判断。
    #[must_use]
    pub fn map_numbers(&self) -> (BTreeSet<u16>, BTreeSet<u16>) {
        let collect = |ts: &[QueryMol]| -> BTreeSet<u16> {
            ts.iter()
                .flat_map(|t| t.atoms.iter())
                .filter_map(super::map_number)
                .collect()
        };
        (collect(&self.reactants), collect(&self.products))
    }

    /// 同一侧内重复的映射号。
    ///
    /// 映射号在一侧内必须唯一 —— 重复的话"产物里的 `:1` 对应哪个反应物原子"
    /// 就没有答案了,而产物构建会静默地取其中一个。
    #[must_use]
    pub fn duplicate_map_numbers(&self) -> Vec<u16> {
        let dups = |ts: &[QueryMol]| -> Vec<u16> {
            let mut count: BTreeMap<u16, usize> = BTreeMap::new();
            for t in ts {
                for a in &t.atoms {
                    if let Some(n) = super::map_number(a) {
                        *count.entry(n).or_default() += 1;
                    }
                }
            }
            count
                .into_iter()
                .filter(|&(_, c)| c > 1)
                .map(|(n, _)| n)
                .collect()
        };
        let mut out = dups(&self.reactants);
        out.extend(dups(&self.products));
        out.sort_unstable();
        out.dedup();
        out
    }
}

/// 解析一条反应 SMARTS。
///
/// # Errors
/// 语法错误、或分段数不是 2 段(`>>`)/3 段(`>`…`>`)时返回 [`ParseError`]。
pub fn parse_reaction(input: &str) -> Result<Reaction> {
    let src = input.as_bytes();
    let sections = split_sections(src)?;
    let (reactants, agents, products) = match sections.len() {
        2 => (sections[0].clone(), Vec::new(), sections[1].clone()),
        3 => (
            sections[0].clone(),
            sections[1].clone(),
            sections[2].clone(),
        ),
        n => {
            return Err(ParseError::new(
                K::BadBracketAtom(if n < 2 {
                    "反应缺少 `>>`"
                } else {
                    "反应的 `>` 分段超过三段"
                }),
                0,
                src,
            ))
        }
    };

    let build = |ranges: &[(usize, usize)]| -> Result<Vec<QueryMol>> {
        ranges
            .iter()
            .map(|&(lo, hi)| {
                let text = std::str::from_utf8(&src[lo..hi])
                    .map_err(|_| ParseError::new(K::UnexpectedEnd, lo, src))?;
                super::mol::parse(text).map_err(|e| ParseError::new(e.kind, lo + e.pos, src))
            })
            .collect()
    };

    Ok(Reaction {
        reactants: build(&reactants)?,
        agents: build(&agents)?,
        products: build(&products)?,
    })
}

/// 把输入切成若干段,每段再切成若干模板的字节区间。
fn split_sections(src: &[u8]) -> Result<Vec<Vec<(usize, usize)>>> {
    let mut sections: Vec<Vec<(usize, usize)>> = Vec::new();
    let mut current: Vec<(usize, usize)> = Vec::new();
    let (mut brackets, mut parens) = (0i32, 0i32);
    let mut start = 0usize;

    for i in 0..=src.len() {
        let b = src.get(i).copied();
        match b {
            Some(b'[') => brackets += 1,
            Some(b']') => brackets -= 1,
            Some(b'(') if brackets == 0 => parens += 1,
            Some(b')') if brackets == 0 => parens -= 1,
            _ => {}
        }
        if brackets != 0 || parens != 0 {
            continue;
        }
        // `->` 里的 `>` 是配位键,不是分段符
        let is_sep = matches!(b, Some(b'>')) && (i == 0 || src[i - 1] != b'-');
        let is_dot = matches!(b, Some(b'.'));
        let is_end = b.is_none();

        if is_sep || is_dot || is_end {
            if start < i {
                current.push(strip_group_parens(src, start, i));
            }
            start = i + 1;
            if is_sep || is_end {
                sections.push(std::mem::take(&mut current));
            }
        }
    }
    if brackets != 0 {
        return Err(ParseError::new(
            K::BadBracketAtom("方括号未闭合"),
            src.len(),
            src,
        ));
    }
    if parens != 0 {
        return Err(ParseError::new(K::UnbalancedParen, src.len(), src));
    }
    Ok(sections)
}

/// 整个模板被一对括号裹住时剥掉它 —— 那是**分组**,不是分支。
///
/// `(A.B)` 里的两个片段属于同一个模板;`C(C)C` 的括号则是分支,不能剥。
/// 判据是"第一个字符是 `(` 且与它配对的 `)` 正好是最后一个字符"。
fn strip_group_parens(src: &[u8], lo: usize, hi: usize) -> (usize, usize) {
    if src.get(lo) != Some(&b'(') || src.get(hi - 1) != Some(&b')') {
        return (lo, hi);
    }
    let (mut depth, mut brackets) = (0i32, 0i32);
    for (offset, &c) in src[lo..hi].iter().enumerate() {
        match c {
            b'[' => brackets += 1,
            b']' => brackets -= 1,
            b'(' if brackets == 0 => depth += 1,
            b')' if brackets == 0 => {
                depth -= 1;
                // 提前归零说明这对括号在中途就闭合了,是分支不是分组
                if depth == 0 && lo + offset != hi - 1 {
                    return (lo, hi);
                }
            }
            _ => {}
        }
    }
    (lo + 1, hi - 1)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn r(s: &str) -> Reaction {
        parse_reaction(s).unwrap_or_else(|e| panic!("{s}:\n{}", e.render()))
    }

    fn shape(s: &str) -> (usize, usize, usize) {
        let x = r(s);
        (x.reactants.len(), x.agents.len(), x.products.len())
    }

    /// 分段与模板切分。
    #[test]
    fn sections_and_templates() {
        assert_eq!(shape("[C:1][OH:2]>>[C:1][Cl:2]"), (1, 0, 1));
        assert_eq!(shape("[C:1].[N:2]>>[C:1][N:2]"), (2, 0, 1), "`.` 分模板");
        assert_eq!(shape("[C:1]>[Pd]>[C:1]Cl"), (1, 1, 1), "中间是试剂");
        assert_eq!(shape("[N:1]->[Cu:2]>>[N:1].[Cu:2]"), (1, 0, 2));
    }

    /// 括号分组:`(A.B)` 是**一个**模板里的两个片段。
    #[test]
    fn grouping_parens_make_one_template() {
        let x = r("([C:1].[N:2])>>[C:1][N:2]");
        assert_eq!(x.reactants.len(), 1, "括号把两个片段捆成一个模板");
        assert_eq!(x.reactants[0].num_atoms(), 2);

        // 分支括号不能被当成分组
        let x = r("[C:1](=[O:2])[OH]>>[C:1](=[O:2])Cl");
        assert_eq!(x.reactants.len(), 1);
        assert_eq!(x.reactants[0].num_atoms(), 3);
    }

    /// `->` 里的 `>` 属于配位键,不是分段符。
    ///
    /// 切错的话 `[N:1]->[Cu:2]>>...` 会被读成四段,而且报的错完全指不到根因。
    #[test]
    fn dative_arrow_is_not_a_separator() {
        let x = r("[N:1]->[Cu:2]>>[N:1].[Cu:2]");
        assert_eq!(x.reactants.len(), 1);
        assert_eq!(x.reactants[0].num_atoms(), 2, "N 和 Cu 在同一个模板里");
        assert_eq!(x.reactants[0].num_bonds(), 1);
        assert_eq!(x.products.len(), 2, "产物侧的 `.` 才是分隔符");
    }

    /// 方括号里的 `.` 与 `>` 不参与切分。
    #[test]
    fn separators_inside_brackets_are_ignored() {
        let x = r("[$(C.C):1]>>[C:1]");
        assert_eq!(x.reactants.len(), 1, "递归 SMARTS 里的 `.` 不分模板");
        assert_eq!(x.reactants[0].num_atoms(), 1);
    }

    /// 映射号:两侧的集合,以及一侧内的重复。
    #[test]
    fn map_numbers() {
        let x = r("[C:1][OH:2]>>[C:1][Cl:2]");
        let (lhs, rhs) = x.map_numbers();
        assert_eq!(lhs, [1, 2].into_iter().collect());
        assert_eq!(rhs, [1, 2].into_iter().collect());
        assert!(x.duplicate_map_numbers().is_empty());

        // 反应物侧的 OH 没映射 → 会被删掉;产物侧的 Cl 没映射 → 新建
        let x = r("[C:1][OH]>>[C:1]Cl");
        let (lhs, rhs) = x.map_numbers();
        assert_eq!(lhs, [1].into_iter().collect());
        assert_eq!(rhs, [1].into_iter().collect());

        // 同一侧内重复的映射号要报得出来
        let x = r("[C:1][C:1]>>[C:1]");
        assert_eq!(x.duplicate_map_numbers(), vec![1]);
    }

    /// 分段数不对要报错。
    #[test]
    fn wrong_section_count_is_an_error() {
        assert!(parse_reaction("CC").is_err(), "没有 `>>`");
        assert!(parse_reaction("A>B>C>D").is_err(), "四段");
    }
}
