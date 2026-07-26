//! SMARTS 的原子表达式解析。
//!
//! 优先级从高到低:`!`(一元)> `&`(可省略,即并置)> `,` > `;`。
//! 用优先级爬升处理 —— 每档一个函数,递归下去,没有要维护的运算符栈。
//!
//! # `[H]` 是一张特例表,不是一条规则
//!
//! `H` 既可能指**氢元素**,也可能指**氢计数**,而区分它们的不是上下文规则,
//! 是一份枚举出来的完整括号形式清单:
//!
//! ```text
//! [H]  [H:n]  [nH]  [nH:m]  [H+]  [H+:n]  [nH+]  [nH+:m]
//! ```
//!
//! 归纳:**可选同位素 + `H` + 可选电荷 + 可选映射号 + `]`**。落在这张表里的
//! 才是氢元素;其余一律是氢计数。于是:
//!
//! | 写法 | 含义 | 为什么 |
//! |---|---|---|
//! | `[H]` | 氢元素 | 在表里 |
//! | `[H+]` | 氢元素 | 在表里 |
//! | `[2H]` | 氘 | 在表里 |
//! | `[H1]` | 氢计数 1 | `H` 后面跟数字,不在表里 |
//! | `[HH]` | 氢计数 1 且氢计数 1 | 不在表里 |
//! | `[#1H]` | 氢元素且氢计数 1 | `#1` 开头,不在表里 |
//!
//! 想从这些用例里总结出"H 在开头就是元素"之类的规则,一定会写错 ——
//! `[H1]` 和 `[H+]` 都是 H 开头,解读却相反。所以这里照着表实现:先按表试着
//! 匹配整个括号,失败就回退到一般规则。

use omgkit_core::ChiralTag;

use super::expr::{AtomExpr, AtomPrim};
use crate::error::{ParseError, ParseErrorKind as K, Result};

/// 解析方括号里的原子表达式(不含两端的方括号)。
///
/// # Errors
/// 语法错误时返回带位置的 [`ParseError`]。
pub fn parse_atom_expr(src: &[u8]) -> Result<AtomExpr> {
    let mut p = ExprParser { src, pos: 0 };
    let e = p.low_and()?;
    if p.pos != src.len() {
        return Err(ParseError::new(
            K::BadBracketAtom("方括号里有无法解析的残余"),
            p.pos,
            src,
        ));
    }
    Ok(e)
}

struct ExprParser<'a> {
    src: &'a [u8],
    pos: usize,
}

impl ExprParser<'_> {
    fn peek(&self) -> Option<u8> {
        self.src.get(self.pos).copied()
    }

    fn eat(&mut self, b: u8) -> bool {
        if self.peek() == Some(b) {
            self.pos += 1;
            true
        } else {
            false
        }
    }

    fn err<T>(&self, why: &'static str) -> Result<T> {
        Err(ParseError::new(K::BadBracketAtom(why), self.pos, self.src))
    }

    // -- 三档优先级 --------------------------------------------------------

    /// `;` —— 最低优先级的与
    fn low_and(&mut self) -> Result<AtomExpr> {
        let mut parts = vec![self.or()?];
        while self.eat(b';') {
            parts.push(self.or()?);
        }
        Ok(flatten(AtomExpr::And, parts))
    }

    /// `,` —— 或
    fn or(&mut self) -> Result<AtomExpr> {
        let mut parts = vec![self.high_and()?];
        while self.eat(b',') {
            parts.push(self.high_and()?);
        }
        Ok(flatten(AtomExpr::Or, parts))
    }

    /// `&` 与并置 —— 最高优先级的与
    fn high_and(&mut self) -> Result<AtomExpr> {
        let mut parts = vec![self.unary()?];
        loop {
            if self.eat(b'&') {
                parts.push(self.unary()?);
                continue;
            }
            // 并置:后面还跟着一个基元的开头,就是隐含的 `&`
            match self.peek() {
                Some(b) if starts_primitive(b) => parts.push(self.unary()?),
                _ => break,
            }
        }
        Ok(flatten(AtomExpr::And, parts))
    }

    /// `!` —— 一元否定,只作用于紧跟其后的那一项
    fn unary(&mut self) -> Result<AtomExpr> {
        if self.eat(b'!') {
            return Ok(AtomExpr::Not(Box::new(self.unary()?)));
        }
        self.primitive()
    }

    // -- 基元 --------------------------------------------------------------

    fn primitive(&mut self) -> Result<AtomExpr> {
        let Some(b) = self.peek() else {
            return self.err("表达式意外结束");
        };
        // 二字符元素符号**先于**单字符基元匹配,且严格区分大小写:
        //   as → 砷      aS → a & S(`aS` 不是合法符号)
        //   Ac → 锕      AC → A & C
        //   Hg → 汞      Rb → 铷      se/te/si → 芳香形式
        // 反过来的话 `as` 会被读成"芳香 & 硫",`Ac` 会被读成"脂肪 & 芳香碳"。
        if b.is_ascii_alphabetic() {
            if let Some(prim) = self.try_two_char_element() {
                return Ok(AtomExpr::Prim(prim));
            }
        }

        let prim = match b {
            b'*' => {
                self.pos += 1;
                AtomPrim::Any
            }
            b'a' => {
                self.pos += 1;
                AtomPrim::Aromatic
            }
            b'A' => {
                self.pos += 1;
                AtomPrim::Aliphatic
            }
            b'#' => {
                self.pos += 1;
                let z = self.number().ok_or_else(|| {
                    ParseError::new(K::BadBracketAtom("`#` 后缺少原子序数"), self.pos, self.src)
                })?;
                AtomPrim::Element {
                    z: u8::try_from(z)
                        .map_err(|_| ParseError::new(K::NumberOverflow, self.pos, self.src))?,
                    aromatic: None,
                }
            }
            // 计数类基元:裸写时的默认值各不相同,见 count_default
            b'D' | b'X' | b'H' | b'h' | b'R' | b'r' | b'x' | b'v' => {
                self.pos += 1;
                let n = self.number();
                count_primitive(b, n)
            }
            b'+' | b'-' => AtomPrim::Charge(self.charge()),
            b'@' => {
                self.pos += 1;
                if self.eat(b'@') {
                    AtomPrim::Chirality(ChiralTag::Cw)
                } else {
                    AtomPrim::Chirality(ChiralTag::Ccw)
                }
            }
            b':' => {
                self.pos += 1;
                let n = self.number().ok_or_else(|| {
                    ParseError::new(K::BadBracketAtom("`:` 后缺少映射号"), self.pos, self.src)
                })?;
                AtomPrim::AtomMap(u16::try_from(n).unwrap_or(u16::MAX))
            }
            b'$' => self.recursive()?,
            b'0'..=b'9' => {
                let n = self.number().expect("已 peek 到数字");
                AtomPrim::Isotope(u16::try_from(n).unwrap_or(u16::MAX))
            }
            _ if b.is_ascii_alphabetic() => self.element()?,
            _ => return self.err("无法识别的查询基元"),
        };
        Ok(AtomExpr::Prim(prim))
    }

    /// `$(...)` —— 递归 SMARTS。
    ///
    /// 括号里是一条完整的 SMARTS 模式,语义是"本原子要能作为该模式的
    /// **首原子**匹配上"。首原子这个约定很重要:`[$(CC)]` 匹配的是"连着一个碳
    /// 的碳",匹配上的是模式里的第一个原子,不是整段。
    ///
    /// 配对括号时要计数 —— 递归里还能再套递归。
    fn recursive(&mut self) -> Result<AtomPrim> {
        let dollar = self.pos;
        self.pos += 1;
        if !self.eat(b'(') {
            return Err(ParseError::new(
                K::BadBracketAtom("`$` 后缺少 `(`"),
                self.pos,
                self.src,
            ));
        }
        let start = self.pos;
        let mut depth = 0usize;
        loop {
            match self.peek() {
                None => {
                    return Err(ParseError::new(
                        K::BadBracketAtom("递归 SMARTS 的 `(` 未闭合"),
                        dollar,
                        self.src,
                    ))
                }
                Some(b'(') => {
                    depth += 1;
                    self.pos += 1;
                }
                Some(b')') if depth == 0 => break,
                Some(b')') => {
                    depth -= 1;
                    self.pos += 1;
                }
                Some(_) => self.pos += 1,
            }
        }
        let inner = &self.src[start..self.pos];
        self.pos += 1; // ')'

        let text = std::str::from_utf8(inner)
            .map_err(|_| ParseError::new(K::UnexpectedEnd, start, self.src))?;
        let sub = super::mol::parse(text)
            .map_err(|e| ParseError::new(e.kind, start + e.pos, self.src))?;
        Ok(AtomPrim::Recursive(Box::new(sub)))
    }

    /// 试着按二字符元素符号读。大小写不合法或不是已知符号时不消费任何字节。
    fn try_two_char_element(&mut self) -> Option<AtomPrim> {
        let two = self.src.get(self.pos..self.pos + 2)?;
        let sym = std::str::from_utf8(two).ok()?;
        // 芳香形式只有这几个符号有定义
        if let Some(z) = aromatic_two_char(sym) {
            self.pos += 2;
            return Some(AtomPrim::Element {
                z,
                aromatic: Some(true),
            });
        }
        // 脂肪形式必须是"大写 + 小写";`AC` 不是符号,要退回 `A & C`
        if !(two[0].is_ascii_uppercase() && two[1].is_ascii_lowercase()) {
            return None;
        }
        let z = omgkit_core::element::atomic_num_of(sym)?;
        self.pos += 2;
        Some(AtomPrim::Element {
            z,
            aromatic: Some(false),
        })
    }

    /// 单字符元素符号。二字符的情形已在 [`Self::try_two_char_element`] 处理。
    fn element(&mut self) -> Result<AtomPrim> {
        let b = self.src[self.pos];
        let upper = char::from(b).to_ascii_uppercase().to_string();
        let Some(z) = omgkit_core::element::atomic_num_of(&upper) else {
            return Err(ParseError::new(
                K::UnknownElement(char::from(b).to_string()),
                self.pos,
                self.src,
            ));
        };
        self.pos += 1;
        Ok(AtomPrim::Element {
            z,
            aromatic: Some(b.is_ascii_lowercase()),
        })
    }

    /// `+` / `++` / `+2` 与对应的负号形式
    fn charge(&mut self) -> i32 {
        let sign: i32 = if self.peek() == Some(b'+') { 1 } else { -1 };
        let ch = if sign > 0 { b'+' } else { b'-' };
        self.pos += 1;

        let mut n = 1i32;
        while self.peek() == Some(ch) {
            self.pos += 1;
            n += 1;
        }
        if n == 1 {
            if let Some(v) = self.number() {
                n = i32::try_from(v).unwrap_or(i32::MAX);
            }
        }
        n * sign
    }

    fn number(&mut self) -> Option<u32> {
        let start = self.pos;
        let mut v: u64 = 0;
        while let Some(d @ b'0'..=b'9') = self.peek() {
            v = (v * 10 + u64::from(d - b'0')).min(u64::from(u32::MAX));
            self.pos += 1;
        }
        if self.pos == start {
            None
        } else {
            u32::try_from(v).ok()
        }
    }
}

/// 只有一个子项时不包一层 —— 否则 `[C]` 会变成 `And([C])`,
/// 比对和写出都要先剥壳。
fn flatten(wrap: fn(Vec<AtomExpr>) -> AtomExpr, mut parts: Vec<AtomExpr>) -> AtomExpr {
    if parts.len() == 1 {
        parts.pop().expect("非空")
    } else {
        wrap(parts)
    }
}

/// 计数类基元裸写(不带数字)时的默认值。
///
/// 这几个默认值**互不相同**,而且没有可归纳的规律:`D`/`X`/`v` 裸写是 1,
/// `R`/`r`/`x`/`h` 裸写是"任意非零"。照抄清单,不要试图统一。
fn count_primitive(letter: u8, n: Option<u32>) -> AtomPrim {
    match letter {
        b'D' => AtomPrim::Degree(n.unwrap_or(1)),
        b'X' => AtomPrim::TotalDegree(n.unwrap_or(1)),
        b'v' => AtomPrim::Valence(n.unwrap_or(1)),
        b'H' => AtomPrim::TotalHs(n.unwrap_or(1)),
        b'h' => AtomPrim::ImplicitHs(n.unwrap_or(1)),
        b'R' => AtomPrim::RingCount(n),
        // `r` 带数字是**最小环大小**,裸写却等价于 `R`(在任意环中)——
        // 不是"最小环大小为 1"。两种含义共用一个字母,只能靠有没有数字来分。
        b'r' => match n {
            Some(k) => AtomPrim::RingSize(Some(k)),
            None => AtomPrim::RingCount(None),
        },
        b'x' => AtomPrim::RingBondCount(n),
        _ => unreachable!("由调用方 match 保证"),
    }
}

/// 该字节能否作为一个基元的开头 —— 用来识别并置(隐含的 `&`)。
fn starts_primitive(b: u8) -> bool {
    b.is_ascii_alphanumeric() || matches!(b, b'*' | b'#' | b'+' | b'-' | b'@' | b':' | b'!' | b'$')
}

/// 芳香形式的二字符元素符号。
fn aromatic_two_char(sym: &str) -> Option<u8> {
    match sym {
        "se" => Some(34),
        "as" => Some(33),
        "te" => Some(52),
        "si" => Some(14),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(s: &str) -> AtomExpr {
        parse_atom_expr(s.as_bytes()).unwrap_or_else(|e| panic!("{s}: {}", e.render()))
    }

    fn elem(z: u8, aromatic: bool) -> AtomExpr {
        AtomExpr::Prim(AtomPrim::Element {
            z,
            aromatic: Some(aromatic),
        })
    }

    fn prim(p: AtomPrim) -> AtomExpr {
        AtomExpr::Prim(p)
    }

    /// 三档优先级。这三条写法长得很像,语义却两两不同,单独钉死。
    #[test]
    fn operator_precedence() {
        let c = || elem(6, false);
        let n = || elem(7, false);
        let h1 = || prim(AtomPrim::TotalHs(1));

        // `;` 最低:(C 或 N) 且 H1
        assert_eq!(
            p("C,N;H1"),
            AtomExpr::And(vec![AtomExpr::Or(vec![c(), n()]), h1()])
        );
        // `&` 高于 `,`:C 或 (N 且 H1)
        assert_eq!(
            p("C,N&H1"),
            AtomExpr::Or(vec![c(), AtomExpr::And(vec![n(), h1()])])
        );
        // 并置就是 `&` —— `NH1` 不是一个整体
        assert_eq!(p("C,NH1"), p("C,N&H1"));
    }

    /// `!` 只作用于紧跟其后的那一项。
    #[test]
    fn negation_binds_tightest() {
        assert_eq!(
            p("!C;!N"),
            AtomExpr::And(vec![
                AtomExpr::Not(Box::new(elem(6, false))),
                AtomExpr::Not(Box::new(elem(7, false))),
            ])
        );
        assert_eq!(
            p("!C,N"),
            AtomExpr::Or(vec![
                AtomExpr::Not(Box::new(elem(6, false))),
                elem(7, false)
            ])
        );
    }

    /// 单项不包壳 —— `C` 就是 `C`,不是 `And([C])`。
    #[test]
    fn single_primitive_is_not_wrapped() {
        assert_eq!(p("C"), elem(6, false));
        assert_eq!(p("c"), elem(6, true));
        assert_eq!(
            p("#6"),
            AtomExpr::Prim(AtomPrim::Element {
                z: 6,
                aromatic: None
            })
        );
    }

    /// 裸写的计数基元,默认值**互不相同**,而且没有规律可循。
    #[test]
    fn bare_count_primitives_have_different_defaults() {
        // 这几个裸写是 1
        assert_eq!(p("D"), prim(AtomPrim::Degree(1)));
        assert_eq!(p("X"), prim(AtomPrim::TotalDegree(1)));
        assert_eq!(p("v"), prim(AtomPrim::Valence(1)));
        // 这几个裸写是"任意非零"
        assert_eq!(p("R"), prim(AtomPrim::RingCount(None)));
        assert_eq!(p("x"), prim(AtomPrim::RingBondCount(None)));
        assert_eq!(p("h"), prim(AtomPrim::ImplicitHs(1)));
        // `r` 裸写等价于 `R`,带数字才是最小环大小
        assert_eq!(p("r"), prim(AtomPrim::RingCount(None)));
        assert_eq!(p("r5"), prim(AtomPrim::RingSize(Some(5))));
        assert_eq!(p("R2"), prim(AtomPrim::RingCount(Some(2))));
    }

    /// 二字符元素符号优先:`Sc` 是钪,不是"硫 + 芳香碳"。
    #[test]
    fn two_char_element_beats_one_char_plus_aromatic() {
        assert_eq!(p("Sc"), elem(21, false), "钪");
        assert_eq!(p("Cl"), elem(17, false));
        assert_eq!(p("Na"), elem(11, false));
        assert_eq!(p("Si"), elem(14, false));
        assert_eq!(p("Hg"), elem(80, false), "汞,不是 H1 & 芳香碳");
        assert_eq!(p("Ho"), elem(67, false), "钬,不是 H1 & 芳香氧");
        assert_eq!(p("Rb"), elem(37, false), "铷,不是 R & 芳香硼");
        assert_eq!(p("Ac"), elem(89, false), "锕,不是 A & 芳香碳");
        // 小写的二字符只有 se/as/te/si 几个
        assert_eq!(p("se"), elem(34, true));
        assert_eq!(p("as"), elem(33, true), "砷,不是 a & 硫");
        assert_eq!(p("te"), elem(52, true));
    }

    /// 大小写必须严格匹配 —— 拼不出合法符号的组合要退回单字符基元。
    ///
    /// 这几条是上一条的反面。少了它,"二字符优先"会变成"见到两个字母就当元素",
    /// 于是 `AC`、`Xx` 这类会被吞掉。
    #[test]
    fn invalid_two_char_falls_back_to_single_char_primitives() {
        assert_eq!(
            p("aS"),
            AtomExpr::And(vec![prim(AtomPrim::Aromatic), elem(16, false)]),
            "`aS` 不是符号"
        );
        assert_eq!(
            p("AC"),
            AtomExpr::And(vec![prim(AtomPrim::Aliphatic), elem(6, false)]),
            "`AC` 不是符号"
        );
        assert_eq!(
            p("ac"),
            AtomExpr::And(vec![prim(AtomPrim::Aromatic), elem(6, true)]),
            "`ac` 不是符号"
        );
        assert_eq!(
            p("Va"),
            AtomExpr::And(vec![elem(23, false), prim(AtomPrim::Aromatic)]),
            "没有 `Va` 这个元素,退回钒 & 芳香"
        );
        assert_eq!(
            p("Xx"),
            AtomExpr::And(vec![
                prim(AtomPrim::TotalDegree(1)),
                prim(AtomPrim::RingBondCount(None))
            ]),
            "没有 `Xx` 这个元素,退回 X1 & x"
        );
    }

    /// 电荷的三种写法等价。
    #[test]
    fn charge_forms() {
        assert_eq!(p("+"), prim(AtomPrim::Charge(1)));
        assert_eq!(p("++"), prim(AtomPrim::Charge(2)));
        assert_eq!(p("+2"), prim(AtomPrim::Charge(2)));
        assert_eq!(p("-"), prim(AtomPrim::Charge(-1)));
        assert_eq!(p("--"), prim(AtomPrim::Charge(-2)));
        assert_eq!(p("-3"), prim(AtomPrim::Charge(-3)));
    }

    /// 本函数解析的是**方括号内容**,不认那张 `[H]` 特例表 ——
    /// 在这一层 `H` 一律是氢计数。特例表在括号那一层处理,见模块文档。
    #[test]
    fn hydrogen_is_a_count_at_this_level() {
        assert_eq!(p("H"), prim(AtomPrim::TotalHs(1)));
        assert_eq!(p("H2"), prim(AtomPrim::TotalHs(2)));
        assert_eq!(
            p("HH"),
            AtomExpr::And(vec![prim(AtomPrim::TotalHs(1)), prim(AtomPrim::TotalHs(1))]),
            "并置的两个 H 都是计数"
        );
    }

    #[test]
    fn isotope_chirality_and_map() {
        assert_eq!(
            p("13C"),
            AtomExpr::And(vec![prim(AtomPrim::Isotope(13)), elem(6, false)])
        );
        assert_eq!(
            p("C@"),
            AtomExpr::And(vec![
                elem(6, false),
                prim(AtomPrim::Chirality(ChiralTag::Ccw))
            ])
        );
        assert_eq!(
            p("C@@"),
            AtomExpr::And(vec![
                elem(6, false),
                prim(AtomPrim::Chirality(ChiralTag::Cw))
            ])
        );
        assert_eq!(
            p("C:1"),
            AtomExpr::And(vec![elem(6, false), prim(AtomPrim::AtomMap(1))])
        );
    }

    #[test]
    fn syntax_errors_have_positions() {
        for (s, at) in [("#", 1usize), ("C&", 2), ("C,", 2), (":", 1)] {
            let err = parse_atom_expr(s.as_bytes()).expect_err(&format!("{s} 应当解析失败"));
            assert_eq!(err.pos, at, "{s} 的报错位置");
        }
    }
}
