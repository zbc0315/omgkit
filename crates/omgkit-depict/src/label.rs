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

use crate::geom::Point2;
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
    pub(crate) fn text(&self) -> &str {
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
    pub(crate) fn width_em(&self) -> f64 {
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

/// 标签往哪个方向伸。
///
/// 氢(以及电荷、同位素)相对元素符号摆在哪 —— 四个方向。[`HSide`] 是它退化到
/// 左右两向的样子,`East` 对应 [`HSide::Right`],`West` 对应 [`HSide::Left`]。
///
/// # 为什么要有上下两向
///
/// 对乙酰氨基酚的酰胺氮:两根键一根指左上、一根指右上,横向分量几乎抵消。只有
/// 左右两向可选时,氢只能挤到某一根键的下面,而**正下方整片是空的**。
///
/// 方向由 [`crate::render::label_dir`] 按键的走向定。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LabelDir {
    /// 氢在符号右边,`NH2`
    East,
    /// 氢在符号左边,`H2N`
    West,
    /// 氢在符号**上方**,另起一行
    North,
    /// 氢在符号**下方**,另起一行
    South,
}

impl LabelDir {
    /// 是不是竖着排。
    #[must_use]
    pub fn is_vertical(self) -> bool {
        matches!(self, LabelDir::North | LabelDir::South)
    }
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
    /// 整串的中心该放在原子**横向偏开多少**的地方,单位是键长。
    ///
    /// # 为什么不能把整串居中
    ///
    /// 键连的是**元素符号**那个原子,不是整串。把 `OH` 整串居中,O 就落在原子
    /// 位置左边 —— 实测 ACS 下 `OH` 全宽 1.042 个键长、`O` 单独 0.540,
    /// **O 的中心离原子位置 0.251 个键长,整整四分之一根键**。
    ///
    /// `h_side` 把氢甩到键来向的反面,横着的键因此看着还对;**竖直的键就露馅**
    /// —— 它落在盒的上/下边缘、横向居中,正好落在 O 和 H 中间。
    ///
    /// 所以把整串按这个量挪开,让符号回到原子上。盒还是那个盒,只是**盒心不在
    /// 原子上**了 —— 裁键、判碰撞都要跟着算,见 `render::trim` 与 `refine::radii`。
    pub dx: f64,
    /// 同 [`Label::dx`],纵向的那一半,**布局坐标(y 向上)**。
    ///
    /// 横排时恒为 0。竖排时符号自己占一行、氢占另一行,盒心因此上下偏开。
    ///
    /// # 坐标系
    ///
    /// `Point2` 在本库同时用于布局系与画布系。**别直接读这个字段做几何** ——
    /// 用 [`Label::offset`](Label::offset)(布局系)或
    /// [`Label::offset_canvas`](Label::offset_canvas)(画布系),哪个系在调用点
    /// 就看得见了。
    pub dy: f64,
    /// 竖排时第二行(氢那行)相对符号行的纵向偏移,**布局坐标**。横排为 0。
    ///
    /// 负数表示氢在符号下面。
    pub gap: f64,
    /// 竖排时,氢那几段从第几段开始;`None` 表示横排。
    ///
    /// # 只有"纯符号 + 氢"才竖排
    ///
    /// 带电荷、同位素、自由基的标签一律横排。竖排要求**符号那一行只有符号**,
    /// 否则整行居中会把符号推偏 —— 那正是 [`Label::dx`] 当初要解决的问题,
    /// 只是搬进了行内。全量语料实测:3202 个竖排候选里 3118 个(**97.4%**)是
    /// 纯符号 + 氢,带电荷的只有 84 个,同位素与自由基一个都没有。为这 2.6%
    /// 引入"每行各自的 dx"不划算。
    pub stacked: Option<usize>,
}

impl Label {
    /// 纯文本形式,便于测试与调试。下标上标不加标记 —— 要区分请看 `runs`。
    #[must_use]
    pub fn plain(&self) -> String {
        self.runs.iter().map(Run::text).collect()
    }

    /// 盒心相对原子位置的偏移,**布局坐标(y 向上)**。
    #[must_use]
    pub fn offset(&self) -> Point2 {
        Point2::new(self.dx, self.dy)
    }

    /// 盒心相对原子位置的偏移,**画布坐标(y 向下)**。
    #[must_use]
    pub fn offset_canvas(&self) -> Point2 {
        Point2::new(self.dx, -self.dy)
    }

    /// 要落笔的每一行:`(相对原子位置的偏移(布局坐标), 这一行的字段)`。
    ///
    /// 横排只有一行。竖排两行 —— 符号一行、氢一行,**两行都横向居中在原子上**。
    ///
    /// # 为什么由这里给行,而不是让 `Primitive::Text` 支持多行
    ///
    /// 一行一个 `Text` 图元,`Primitive`、`svg`、`raster` 三处一个字都不用改。
    /// SVG 里换行还得显式给第二行 `x`(只给 `dy` 会沿上一行的笔位继续往右排,
    /// 这个坑在 `svg.rs` 里已经踩过一次),发两个图元连这个都绕开了。
    #[must_use]
    pub fn lines(&self) -> Vec<(Point2, &[Run])> {
        match self.stacked {
            None => vec![(Point2::new(self.dx, 0.0), &self.runs[..])],
            Some(at) => {
                // 符号那一行落在原子上(偏移 0),氢那一行在 `gap` 处。
                vec![
                    (Point2::new(0.0, 0.0), &self.runs[..at]),
                    (Point2::new(0.0, self.gap), &self.runs[at..]),
                ]
            }
        }
    }
}

/// 标签怎么摆:横着(氢在符号左/右)还是竖着(氢在符号上/下)。
///
/// 由 [`crate::render::label_at`] 定 —— 竖不竖排看几何([`crate::render::label_dir`]),
/// 而**左右怎么选仍归 [`crate::render::h_side`]**。两件事分开:竖排摆不成时
/// (带电荷等,见 [`Label::stacked`])要回落到横排,那时左右的选择不能丢。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LabelPlace {
    /// 横排,氢在符号的哪一侧
    Horizontal(HSide),
    /// 竖排。`below` 为真时氢在符号**下面**。
    Stacked { below: bool },
}

/// 这个原子的标签能不能竖排。
///
/// 要求:有氢可摆到第二行,且**符号那一行只有符号** —— 带电荷、同位素、
/// 自由基一律不竖排,理由见 [`Label::stacked`]。
#[must_use]
pub fn can_stack(mol: &MolBuilder, atom: u32) -> bool {
    let a = mol.atoms()[atom as usize];
    total_hs(mol, atom) > 0
        && a.formal_charge == 0
        && a.isotope == 0
        && a.num_radical_electrons == 0
}

/// 一个原子该显示什么标签。返回 `None` 表示**不画**(骨架碳)。
///
/// # 什么时候不画
///
/// 中性、无同位素、无自由基、且**至少连着一个重原子**的碳不画 —— 那是骨架的
/// 惯例。孤立的碳(甲烷)要画成 `CH4`,否则图上什么都没有。
#[must_use]
pub fn label_for(mol: &MolBuilder, atom: u32, style: &Style, place: LabelPlace) -> Option<Label> {
    let a = mol.atoms()[atom as usize];
    let plain_carbon = a.atomic_num == 6
        && a.formal_charge == 0
        && a.isotope == 0
        && a.num_radical_electrons == 0
        && mol.degree(atom) > 0;
    if plain_carbon {
        return None;
    }
    Some(build(mol, atom, style, place))
}

/// 同 [`label_for`],但**骨架碳也画出来**。
///
/// 给的是"这个原子必须看得见"的情形。目前只有一处:两根键几乎共线的二度原子
/// —— 相邻两根键连成一条直线,顶点没有拐角,图上根本看不出那里有个原子。
/// 丙二烯 `CH₃CH=C=CHCH₃` 的中心碳是 sp、键角 180°,不画符号的话整个分子读起来
/// 像顺式二烯。RDKit 也是这么办的(`DrawMol::getAtomSymbol`,注释原话是
/// "allenes need a C")。
#[must_use]
pub fn label_forced(mol: &MolBuilder, atom: u32, style: &Style, place: LabelPlace) -> Label {
    build(mol, atom, style, place)
}

fn build(mol: &MolBuilder, atom: u32, style: &Style, place: LabelPlace) -> Label {
    let a = mol.atoms()[atom as usize];
    let hs = total_hs(mol, atom);
    let symbol = element::by_atomic_num(a.atomic_num).map_or("*", |e| e.symbol);
    let em = style.label_size(); // 一个 em 等于多少个键长

    let h_runs = |runs: &mut Vec<Run>| {
        if hs > 0 {
            runs.push(Run::Normal("H".into()));
            if hs > 1 {
                runs.push(Run::Sub(hs.to_string()));
            }
        }
    };

    // ---- 竖排 ----
    // **摆不成竖排就回落横排。** 先前这里是 `debug_assert`:release 下竖排那支
    // 只 push 符号与氢,电荷、同位素、自由基一个都不发 —— `[NH4+]` 会画成
    // `NH4`,那个 `+` 静静没了;0 个氢时第二行是空切片,还会发一条空的
    // `Primitive::Text`。而 `LabelPlace` 与 `label_for` 都是公开的,调用方给什么
    // 由不得我们。回落到 `Right` 与打平时的默认一致。
    let place = match place {
        LabelPlace::Stacked { .. } if !can_stack(mol, atom) => LabelPlace::Horizontal(HSide::Right),
        other => other,
    };
    if let LabelPlace::Stacked { below } = place {
        let mut runs: Vec<Run> = vec![Run::Normal(symbol.into())];
        let at = runs.len();
        h_runs(&mut runs);

        let sym_w = Run::Normal(symbol.into()).width_em();
        let h_w: f64 = runs[at..].iter().map(Run::width_em).sum();
        let h_has_sub = runs[at..].iter().any(|r| matches!(r, Run::Sub(_)));

        // 两行中心的距离。留白照 RDKit 的 1.1 倍字高;**只有氢在上面时**才要
        // 再让开下标下沉的那一点 —— 那时下标正对着符号。氢在下面时下标朝外,
        // 让了纯属白撑高:实测 `NH₂` 两行字盒之间的白会从 0.050 涨到 0.147
        // 个键长,比 `NH` 宽出三倍,看着就是两行没对齐。
        let gap_em = CAP_HEIGHT * 1.1 + if h_has_sub && !below { SUB_DROP } else { 0.0 };
        let gap = if below { -gap_em * em } else { gap_em * em };

        // 盒:符号行占 ±CAP/2,氢行以 `gap` 为中心、占 ±CAP/2(下标再多一点)
        let half_cap = CAP_HEIGHT / 2.0 * em;
        let h_low = gap - half_cap - if h_has_sub { SUB_DROP * em } else { 0.0 };
        let h_high = gap + half_cap;
        let top = half_cap.max(h_high);
        let bottom = (-half_cap).min(h_low);

        return Label {
            runs,
            half_w: sym_w.max(h_w) * em / 2.0,
            half_h: (top - bottom) / 2.0,
            // 两行都横向居中在原子上,所以横向不用挪
            dx: 0.0,
            dy: (top + bottom) / 2.0,
            gap,
            stacked: Some(at),
        };
    }

    // ---- 横排 ----
    let LabelPlace::Horizontal(h_side) = place else {
        unreachable!("竖排在上面已经返回")
    };
    let mut runs: Vec<Run> = Vec::new();

    // 记下元素符号在整串里的位置 —— 键要接到**符号**上,不是接到整串的中心
    let mut before_sym_em = 0.0_f64;
    match h_side {
        HSide::Left => {
            // `H2N–`:氢在前,同位素仍贴着符号
            h_runs(&mut runs);
            if a.isotope != 0 {
                runs.push(Run::Sup(a.isotope.to_string()));
            }
            before_sym_em = runs.iter().map(Run::width_em).sum();
            runs.push(Run::Normal(symbol.into()));
        }
        HSide::Right => {
            if a.isotope != 0 {
                runs.push(Run::Sup(a.isotope.to_string()));
                before_sym_em = runs.iter().map(Run::width_em).sum();
            }
            runs.push(Run::Normal(symbol.into()));
            h_runs(&mut runs);
        }
    }
    let sym_w_em = Run::Normal(symbol.into()).width_em();

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

    let width_em: f64 = runs.iter().map(Run::width_em).sum();

    // 高度:正文占一个大写字高,上标向上、下标向下各多占一点
    let has_sup = runs.iter().any(|r| matches!(r, Run::Sup(_)));
    let has_sub = runs.iter().any(|r| matches!(r, Run::Sub(_)));
    let top = CAP_HEIGHT / 2.0 + if has_sup { SUP_RISE } else { 0.0 };
    let bottom = CAP_HEIGHT / 2.0 + if has_sub { SUB_DROP } else { 0.0 };

    // 整串要往哪边挪,才能让**元素符号**落在原子位置上。
    //
    // 符号中心相对整串左端在 `before_sym + sym_w/2`,整串中心在 `width/2`,
    // 所以要把整串朝相反方向挪这个差。
    let dx = (width_em / 2.0 - before_sym_em - sym_w_em / 2.0) * em;

    Label {
        runs,
        half_w: width_em * em / 2.0,
        half_h: top.max(bottom) * em,
        dx,
        dy: 0.0,
        gap: 0.0,
        stacked: None,
    }
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
        label_for(&m, atom, &Style::ACS_1996, LabelPlace::Horizontal(side)).map(|l| l.plain())
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
    fn the_element_symbol_sits_on_the_atom_not_on_the_string_centre() {
        // **键连的是元素符号那个原子,不是整串。** 把 `OH` 整串居中,O 就落在
        // 原子位置左边:实测 ACS 下 `OH` 全宽 1.042 个键长、`O` 单独 0.540,
        // O 的中心离原子 **0.251 个键长,四分之一根键**。
        //
        // `h_side` 把氢甩到键来向的反面,横着的键因此看着还对;竖直的键就露馅
        // —— 它落在盒的上/下边缘、横向居中,正好落在 O 和 H 中间。
        //
        // 判据从 `runs` **独立重算**符号该在哪,不复用 `build` 里那个中间量。
        for (smi, atom) in [
            ("CCO", 2u32),
            ("NCC", 0),
            ("CC(=O)Nc1ccc(O)cc1", 7),
            ("[NH4+]", 0),
            ("[13CH4]", 0),
        ] {
            let m = prep(smi);
            let sym = element::by_atomic_num(m.atoms()[atom as usize].atomic_num)
                .expect("元素表里有这个原子")
                .symbol;
            for style in &Style::ALL {
                for side in [HSide::Right, HSide::Left] {
                    let Some(l) = label_for(&m, atom, style, LabelPlace::Horizontal(side)) else {
                        continue;
                    };
                    let em = style.label_size();
                    let mut x = 0.0_f64;
                    let mut centre = None;
                    for r in &l.runs {
                        let w = r.width_em();
                        if matches!(r, Run::Normal(t) if t == sym) {
                            centre = Some(x + w / 2.0);
                        }
                        x += w;
                    }
                    let c = centre.expect("标签里该有元素符号");
                    // 整串左端在 `原子 + dx − half_w`
                    let off = l.dx - l.half_w + c * em;
                    assert!(
                        off.abs() < 1e-9,
                        "[{}] {smi} 原子 {atom} 标签 {}:元素符号 {sym} 的中心离原子 {off:.4} 个键长,该是 0",
                        style.name,
                        l.plain()
                    );
                }
            }
        }
    }

    #[test]
    fn the_box_is_wider_when_there_is_more_to_draw() {
        // 包围盒必须随内容变。写成常数不会报错,只会让所有标签按同一个尺寸
        // 避让 —— `I` 和 `W` 的实际宽度差三倍多
        let m = prep("NCC");
        let n = label_for(
            &m,
            0,
            &Style::ACS_1996,
            LabelPlace::Horizontal(HSide::Right),
        )
        .unwrap();
        let m2 = prep("[NH4+]");
        let nh4 = label_for(
            &m2,
            0,
            &Style::ACS_1996,
            LabelPlace::Horizontal(HSide::Right),
        )
        .unwrap();
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
        let acs = label_for(
            &m,
            0,
            &Style::ACS_1996,
            LabelPlace::Horizontal(HSide::Right),
        )
        .unwrap();
        let cd = label_for(
            &m,
            0,
            &Style::CHEMDRAW_DEFAULT,
            LabelPlace::Horizontal(HSide::Right),
        )
        .unwrap();
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
