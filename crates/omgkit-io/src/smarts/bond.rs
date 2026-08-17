//! SMARTS 的键表达式解析。
//!
//! 与原子表达式同一套三档优先级(`&` > `,` > `;`,`!` 一元),但基元少得多,
//! 而且**没有并置**:两个键符号挨着写(`=#`)不是"且",是语法错误 ——
//! 键符号本身就是单字符,并置无从区分"一个复合表达式"与"两条键"。

use super::expr::{BondExpr, BondPrim};
use crate::error::{ParseError, ParseErrorKind as K, Result};

/// 解析一段键表达式。
///
/// # Errors
/// 语法错误时返回带位置的 [`ParseError`]。
pub fn parse_bond_expr(src: &[u8]) -> Result<BondExpr> {
    let mut p = BondParser { src, pos: 0 };
    let e = p.low_and()?;
    if p.pos != src.len() {
        return Err(ParseError::new(
            K::BadBracketAtom("键表达式里有无法解析的残余"),
            p.pos,
            src,
        ));
    }
    Ok(e)
}

struct BondParser<'a> {
    src: &'a [u8],
    pos: usize,
}

impl BondParser<'_> {
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

    fn low_and(&mut self) -> Result<BondExpr> {
        let mut parts = vec![self.or()?];
        while self.eat(b';') {
            parts.push(self.or()?);
        }
        Ok(flatten(BondExpr::And, parts))
    }

    fn or(&mut self) -> Result<BondExpr> {
        let mut parts = vec![self.high_and()?];
        while self.eat(b',') {
            parts.push(self.high_and()?);
        }
        Ok(flatten(BondExpr::Or, parts))
    }

    fn high_and(&mut self) -> Result<BondExpr> {
        let mut parts = vec![self.unary()?];
        loop {
            if self.eat(b'&') {
                parts.push(self.unary()?);
                continue;
            }
            // 并置也是与,和原子表达式一样:`=!@` 就是 `= & !@`。
            //
            // 这里不存在歧义 —— 键表达式夹在两个原子之间,后面没有放第二条键
            // 的位置。所以连 `=#`(既是双键又是三键)这种自相矛盾的写法都是
            // **语法合法**的,只是永远匹配不上。
            match self.peek() {
                Some(b) if starts_bond_expr(b) => parts.push(self.unary()?),
                _ => break,
            }
        }
        Ok(flatten(BondExpr::And, parts))
    }

    fn unary(&mut self) -> Result<BondExpr> {
        if self.eat(b'!') {
            return Ok(BondExpr::Not(Box::new(self.unary()?)));
        }
        Ok(BondExpr::Prim(self.primitive()?))
    }

    fn primitive(&mut self) -> Result<BondPrim> {
        let Some(b) = self.peek() else {
            return Err(ParseError::new(
                K::BadBracketAtom("键表达式意外结束"),
                self.pos,
                self.src,
            ));
        };
        self.pos += 1;
        Ok(match b {
            b'~' => BondPrim::Any,
            // `-` 后面跟 `>` 才是配位键
            b'-' => {
                if self.eat(b'>') {
                    BondPrim::Dative
                } else {
                    BondPrim::Single
                }
            }
            b'<' => {
                if self.eat(b'-') {
                    BondPrim::DativeReversed
                } else {
                    return Err(ParseError::new(
                        K::UnexpectedChar('<'),
                        self.pos - 1,
                        self.src,
                    ));
                }
            }
            b'=' => BondPrim::Double,
            b'#' => BondPrim::Triple,
            b'$' => BondPrim::Quadruple,
            b':' => BondPrim::Aromatic,
            b'@' => BondPrim::InRing,
            b'/' => BondPrim::UpRight,
            b'\\' => BondPrim::DownRight,
            _ => {
                return Err(ParseError::new(
                    K::UnexpectedChar(char::from(b)),
                    self.pos - 1,
                    self.src,
                ))
            }
        })
    }
}

fn flatten(wrap: fn(Vec<BondExpr>) -> BondExpr, mut parts: Vec<BondExpr>) -> BondExpr {
    if parts.len() == 1 {
        parts.pop().expect("非空")
    } else {
        wrap(parts)
    }
}

/// 该字节能否作为键表达式的开头。主循环靠它区分"这里是一个键"与
/// "这里是别的东西"。
#[must_use]
pub fn starts_bond_expr(b: u8) -> bool {
    matches!(
        b,
        b'~' | b'-' | b'=' | b'#' | b'$' | b':' | b'@' | b'/' | b'\\' | b'!' | b'<'
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(s: &str) -> BondExpr {
        parse_bond_expr(s.as_bytes()).unwrap_or_else(|e| panic!("{s}: {}", e.render()))
    }

    fn prim(x: BondPrim) -> BondExpr {
        BondExpr::Prim(x)
    }

    #[test]
    fn single_primitives() {
        assert_eq!(p("-"), prim(BondPrim::Single));
        assert_eq!(p("="), prim(BondPrim::Double));
        assert_eq!(p("#"), prim(BondPrim::Triple));
        assert_eq!(p("$"), prim(BondPrim::Quadruple));
        assert_eq!(p(":"), prim(BondPrim::Aromatic));
        assert_eq!(p("~"), prim(BondPrim::Any));
        assert_eq!(p("@"), prim(BondPrim::InRing));
        assert_eq!(p("/"), prim(BondPrim::UpRight));
        assert_eq!(p("\\"), prim(BondPrim::DownRight));
    }

    /// `->` 与 `<-` 是两个字符的整体,不能被 `-` 抢走。
    #[test]
    fn dative_bonds() {
        assert_eq!(p("->"), prim(BondPrim::Dative));
        assert_eq!(p("<-"), prim(BondPrim::DativeReversed));
        // 孤立的 `<` 不合法
        assert!(parse_bond_expr(b"<").is_err());
    }

    /// 与原子表达式同一套优先级。
    #[test]
    fn operator_precedence() {
        assert_eq!(
            p("-,=;@"),
            BondExpr::And(vec![
                BondExpr::Or(vec![prim(BondPrim::Single), prim(BondPrim::Double)]),
                prim(BondPrim::InRing),
            ])
        );
        assert_eq!(
            p("-,=&@"),
            BondExpr::Or(vec![
                prim(BondPrim::Single),
                BondExpr::And(vec![prim(BondPrim::Double), prim(BondPrim::InRing)]),
            ])
        );
        assert_eq!(
            p("!@"),
            BondExpr::Not(Box::new(prim(BondPrim::InRing))),
            "非环键"
        );
    }

    /// 键表达式**有并置**,和原子表达式一样 —— 并置就是与。
    ///
    /// 键符号都是单字符,乍看并置会分不清"一个复合表达式"还是"两条键",
    /// 其实不歧义:键表达式夹在两个原子之间,后面没有放第二条键的位置。
    /// `=!@`(双键且非环键)在真实语料里是常见写法。
    ///
    /// 代价是 `=#` 这种自相矛盾的写法也语法合法 —— 它只是永远匹配不上。
    #[test]
    fn juxtaposition_is_conjunction() {
        assert_eq!(
            p("=!@"),
            BondExpr::And(vec![
                prim(BondPrim::Double),
                BondExpr::Not(Box::new(prim(BondPrim::InRing)))
            ]),
            "双键且非环键"
        );
        assert_eq!(p("-@"), p("-&@"));
        assert_eq!(
            p("=#"),
            BondExpr::And(vec![prim(BondPrim::Double), prim(BondPrim::Triple)]),
            "自相矛盾但语法合法"
        );
    }
}
