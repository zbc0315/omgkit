//! 原子标签:显示什么文本,以及它占多大地方。
//!
//! # 这是绘制反馈进布局的那个接口
//!
//! 布局本身只知道原子是个点。可原子标签要占地方,而**标签占多大是规范相关的**
//! —— ACS 1996 的 10 pt 标签落在 14.4 pt 的键上占 69% 个键长,ChemDraw 默认的
//! 同样 10 pt 标签落在 30 pt 的键上只占 33%。消冲突要用的正是这个尺寸。
//!
//! 所以本模块产出的包围盒会喂回布局,不是只给渲染用。
//!
//! # 字宽的出处
//!
//! 用 Helvetica 的 AFM 标准字宽(千分之一 em 为单位),不做估算。Arial 在设计上
//! 与 Helvetica **度量兼容**,同一张表两者通用 —— 这也是两套内置规范都取无衬线
//! 族的原因。
//!
//! 估算字宽(比如"一律按 0.6 em")的后果不是报错,是标签之间的间距系统性地
//! 偏大或偏小,而且只在某些元素上偏 —— `I` 是 0.278 em,`W` 是 0.944 em,
//! 差三倍多。

use omgkit_core::{element, AtomFlags, MolBuilder};

use crate::style::Style;

/// 上标/下标相对正文的字号比例。排版通例,与期刊规范无关,故不进 [`Style`]。
pub(crate) const SUB_SUP_SCALE: f64 = 0.6;
/// 上标基线相对正文基线的抬升,单位 em。
pub(crate) const SUP_RISE: f64 = 0.36;
/// 下标基线相对正文基线的下沉,单位 em。
pub(crate) const SUB_DROP: f64 = 0.14;
/// Helvetica 的大写字高,单位 em(AFM `CapHeight` 718)。
const CAP_HEIGHT: f64 = 0.718;

/// Helvetica 字宽,千分之一 em。
///
/// 只收原子标签真正会用到的字符:元素符号(大小写字母)、数字、正负号。
/// 表里没有的字符按 `FALLBACK_WIDTH` 计 —— 那只可能来自将来新增的记号,
/// 宁可宽一点也不要窄。
fn glyph_width(c: char) -> f64 {
    const FALLBACK_WIDTH: u32 = 600;
    let w: u32 = match c {
        'A' => 667,
        'B' => 667,
        'C' => 722,
        'D' => 722,
        'E' => 667,
        'F' => 611,
        'G' => 778,
        'H' => 722,
        'I' => 278,
        'J' => 500,
        'K' => 667,
        'L' => 556,
        'M' => 833,
        'N' => 722,
        'O' => 778,
        'P' => 667,
        'Q' => 778,
        'R' => 722,
        'S' => 667,
        'T' => 611,
        'U' => 722,
        'V' => 667,
        'W' => 944,
        'X' => 667,
        'Y' => 667,
        'Z' => 611,
        'a' => 556,
        'b' => 556,
        'c' => 500,
        'd' => 556,
        'e' => 556,
        'f' => 278,
        'g' => 556,
        'h' => 556,
        'i' => 222,
        'j' => 222,
        'k' => 500,
        'l' => 222,
        'm' => 833,
        'n' => 556,
        'o' => 556,
        'p' => 556,
        'q' => 556,
        'r' => 333,
        's' => 500,
        't' => 278,
        'u' => 556,
        'v' => 500,
        'w' => 722,
        'x' => 500,
        'y' => 500,
        'z' => 500,
        '0'..='9' => 556,
        '+' => 584,
        '-' => 333,
        '*' => 389,
        _ => FALLBACK_WIDTH,
    };
    f64::from(w) / 1000.0
}

/// 标签里的一段文本。
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Run {
    /// 正文
    Normal(String),
    /// 下标(氢原子数)
    Sub(String),
    /// 上标(形式电荷、同位素)
    Sup(String),
}

impl Run {
    fn text(&self) -> &str {
        match self {
            Run::Normal(s) | Run::Sub(s) | Run::Sup(s) => s,
        }
    }

    fn scale(&self) -> f64 {
        match self {
            Run::Normal(_) => 1.0,
            Run::Sub(_) | Run::Sup(_) => SUB_SUP_SCALE,
        }
    }

    /// 宽度,单位 em。
    fn width_em(&self) -> f64 {
        self.text().chars().map(glyph_width).sum::<f64>() * self.scale()
    }
}

/// 氢挂在元素符号的哪一侧。
///
/// `H2N–` 与 `–NH2` 是同一个原子的两种写法,差别只在键从哪边过来。挂错侧不会
/// 报错,只会让氢和键叠在一起。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HSide {
    /// 氢写在符号右边,如 `NH2`
    Right,
    /// 氢写在符号左边,如 `H2N`
    Left,
}

/// 一个原子的标签。
#[derive(Debug, Clone, PartialEq)]
pub struct Label {
    /// 从左到右的文本段
    pub runs: Vec<Run>,
    /// 相对原子中心的**半宽**,单位是键长
    pub half_w: f64,
    /// 相对原子中心的**半高**,单位是键长
    pub half_h: f64,
}

impl Label {
    /// 纯文本形式,便于测试与调试。下标上标不加标记 —— 要区分请看 `runs`。
    #[must_use]
    pub fn plain(&self) -> String {
        self.runs.iter().map(Run::text).collect()
    }
}

/// 一个原子该显示什么标签。返回 `None` 表示**不画**(骨架碳)。
///
/// # 什么时候不画
///
/// 中性、无同位素、无自由基、且**至少连着一个重原子**的碳不画 —— 那是骨架的
/// 惯例。孤立的碳(甲烷)要画成 `CH4`,否则图上什么都没有。
#[must_use]
pub fn label_for(mol: &MolBuilder, atom: u32, style: &Style, h_side: HSide) -> Option<Label> {
    let a = mol.atoms()[atom as usize];
    let hs = total_hs(mol, atom);

    let plain_carbon = a.atomic_num == 6
        && a.formal_charge == 0
        && a.isotope == 0
        && a.num_radical_electrons == 0
        && mol.degree(atom) > 0;
    if plain_carbon {
        return None;
    }

    let symbol = element::by_atomic_num(a.atomic_num).map_or("*", |e| e.symbol);

    let mut runs: Vec<Run> = Vec::new();
    let h_runs = |runs: &mut Vec<Run>| {
        if hs > 0 {
            runs.push(Run::Normal("H".into()));
            if hs > 1 {
                runs.push(Run::Sub(hs.to_string()));
            }
        }
    };

    match h_side {
        HSide::Left => {
            // `H2N–`:氢在前,同位素仍贴着符号
            h_runs(&mut runs);
            if a.isotope != 0 {
                runs.push(Run::Sup(a.isotope.to_string()));
            }
            runs.push(Run::Normal(symbol.into()));
        }
        HSide::Right => {
            if a.isotope != 0 {
                runs.push(Run::Sup(a.isotope.to_string()));
            }
            runs.push(Run::Normal(symbol.into()));
            h_runs(&mut runs);
        }
    }

    if a.formal_charge != 0 {
        runs.push(Run::Sup(charge_text(a.formal_charge)));
    }
    if a.num_radical_electrons > 0 {
        // 自由基点。用 `•` 会引入一个不在字宽表里的字符,这里用 `*` ——
        // 表里有它,宽度是真的。
        runs.push(Run::Sup(
            "*".repeat(a.num_radical_electrons.min(3) as usize),
        ));
    }

    let em = style.label_size(); // 一个 em 等于多少个键长
    let width_em: f64 = runs.iter().map(Run::width_em).sum();

    // 高度:正文占一个大写字高,上标向上、下标向下各多占一点
    let has_sup = runs.iter().any(|r| matches!(r, Run::Sup(_)));
    let has_sub = runs.iter().any(|r| matches!(r, Run::Sub(_)));
    let top = CAP_HEIGHT / 2.0 + if has_sup { SUP_RISE } else { 0.0 };
    let bottom = CAP_HEIGHT / 2.0 + if has_sub { SUB_DROP } else { 0.0 };

    Some(Label {
        runs,
        half_w: width_em * em / 2.0,
        half_h: top.max(bottom) * em,
    })
}

/// 形式电荷的上标文本:`+` / `-` / `2+` / `3-`。
fn charge_text(q: i8) -> String {
    let sign = if q > 0 { '+' } else { '-' };
    let n = q.unsigned_abs();
    if n == 1 {
        sign.to_string()
    } else {
        format!("{n}{sign}")
    }
}

/// 总氢数。
///
/// **两个字段互斥**(置了 `NO_IMPLICIT` 的那类隐式氢恒为 0),相加即总数 ——
/// 这是全仓一致的约定,`omgkit-chem` 的 `aromaticity`、`adjust_hs` 与
/// `omgkit-io` 的 `smiles::write` 三处都是这么算的。
fn total_hs(mol: &MolBuilder, atom: u32) -> u8 {
    let a = mol.atoms()[atom as usize];
    debug_assert!(
        !a.flags.contains(AtomFlags::NO_IMPLICIT) || a.num_implicit_hs == 0,
        "置了 NO_IMPLICIT 却还有隐式氢,两个字段不再互斥,相加就会重复计数"
    );
    a.num_explicit_hs.saturating_add(a.num_implicit_hs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::style::Style;

    fn prep(smi: &str) -> MolBuilder {
        let mut m = omgkit_io::smiles::parse(smi).unwrap();
        omgkit_chem::pipeline::sanitize(&mut m).unwrap();
        m
    }

    fn plain(smi: &str, atom: u32, side: HSide) -> Option<String> {
        let m = prep(smi);
        label_for(&m, atom, &Style::ACS_1996, side).map(|l| l.plain())
    }

    #[test]
    fn skeleton_carbons_are_not_drawn_but_lone_ones_are() {
        // 骨架碳不画是惯例;可孤立的碳要是也不画,图上就什么都没有了
        assert_eq!(plain("CCO", 0, HSide::Right), None, "链上的碳不该有标签");
        assert_eq!(plain("CCO", 1, HSide::Right), None);
        assert_eq!(plain("CCO", 2, HSide::Right).as_deref(), Some("OH"));
        assert_eq!(
            plain("C", 0, HSide::Right).as_deref(),
            Some("CH4"),
            "甲烷必须画出来"
        );
    }

    #[test]
    fn hydrogens_sit_on_the_side_the_bond_does_not_come_from() {
        // `H2N–` 与 `–NH2` 是同一个原子。挂错侧不报错,只是氢和键叠在一起
        assert_eq!(plain("NCC", 0, HSide::Right).as_deref(), Some("NH2"));
        assert_eq!(plain("NCC", 0, HSide::Left).as_deref(), Some("H2N"));
    }

    #[test]
    fn charge_isotope_and_radical_show_up() {
        assert_eq!(plain("[NH4+]", 0, HSide::Right).as_deref(), Some("NH4+"));
        assert_eq!(plain("[O-]C", 0, HSide::Right).as_deref(), Some("O-"));
        assert_eq!(plain("[13CH4]", 0, HSide::Right).as_deref(), Some("13CH4"));
        assert_eq!(plain("[Fe+2]", 0, HSide::Right).as_deref(), Some("Fe2+"));
    }

    #[test]
    fn the_box_is_wider_when_there_is_more_to_draw() {
        // 包围盒必须随内容变。写成常数不会报错,只会让所有标签按同一个尺寸
        // 避让 —— `I` 和 `W` 的实际宽度差三倍多
        let m = prep("NCC");
        let n = label_for(&m, 0, &Style::ACS_1996, HSide::Right).unwrap();
        let m2 = prep("[NH4+]");
        let nh4 = label_for(&m2, 0, &Style::ACS_1996, HSide::Right).unwrap();
        assert!(
            nh4.half_w > n.half_w,
            "NH4+ 应当比 NH2 宽:{} vs {}",
            nh4.half_w,
            n.half_w
        );

        // 上标要把盒子撑高
        assert!(nh4.half_h > n.half_h, "带上标的标签应当更高");
    }

    #[test]
    fn the_same_label_takes_twice_the_room_under_acs() {
        // **这是整个架构决定的落点。** 同一个标签、同样 10 pt 字号,在 ACS 上
        // 占的键长比例是 ChemDraw 默认的 2.08 倍(14.4 pt 键 vs 30 pt 键)。
        // 两边若一样大,那把 Style 传进标签就是白传,消冲突也就无从区分规范。
        let m = prep("NCC");
        let acs = label_for(&m, 0, &Style::ACS_1996, HSide::Right).unwrap();
        let cd = label_for(&m, 0, &Style::CHEMDRAW_DEFAULT, HSide::Right).unwrap();
        assert_eq!(acs.plain(), cd.plain(), "文本本身与规范无关");
        let ratio = acs.half_w / cd.half_w;
        assert!(
            (ratio - 30.0 / 14.4).abs() < 1e-9,
            "宽度比应当正好是键长之比 30/14.4 = 2.083,实得 {ratio}"
        );
    }

    #[test]
    fn glyph_widths_are_the_real_helvetica_ones() {
        // 抄错字宽不会让任何东西报错,只会让间距系统性偏一点。这里钉住几个
        // 差别最大的:I 最窄、W 最宽,相差三倍多。
        assert!((glyph_width('I') - 0.278).abs() < 1e-12);
        assert!((glyph_width('W') - 0.944).abs() < 1e-12);
        assert!((glyph_width('C') - 0.722).abs() < 1e-12);
        assert!((glyph_width('O') - 0.778).abs() < 1e-12);
        assert!((glyph_width('5') - 0.556).abs() < 1e-12);
        assert!(glyph_width('W') > glyph_width('I') * 3.0);
    }
}
