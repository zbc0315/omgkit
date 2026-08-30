//! SVG 后端:把 [`Scene`] 序列化成字符串。**零依赖。**
//!
//! 这里只做序列化,一点几何都不算 —— 几何在 [`render`](crate::render) 里算过
//! 一遍了。加别的后端不必重写几何,而几何错了所有后端一起错、被同一批判据抓住。
//!
//! # 用 `<text>` 而不是把字转成路径
//!
//! 转成路径的 SVG 到哪都长一样,`<text>` 则依赖查看器装了什么字体。这里选
//! `<text>`,因为转路径要内嵌一套字体轮廓数据 —— 那是几十上百 KB 的表,而且
//! 它自己就需要一套判据来保证抄对了。
//!
//! 代价是**实打实的**:装不到 Arial/Helvetica 的机器上,标签宽度会与
//! [`label`](crate::label) 按 AFM 字宽算出来的不一致,于是布局时按标签留出的
//! 空隙对不上。要求到哪都像素级一致的场合,应当先转成位图再分发。

use crate::label::Run;
use crate::render::{Primitive, Scene};
use crate::style::Style;

/// 把一张 [`Scene`] 写成 SVG。
#[must_use]
pub fn to_svg(scene: &Scene, style: &Style) -> String {
    let mut s = String::with_capacity(1024 + scene.items.len() * 96);
    s.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{:.2}\" height=\"{:.2}\" \
         viewBox=\"0 0 {:.2} {:.2}\">\n",
        scene.width, scene.height, scene.width, scene.height
    ));
    // 白底:没有它,深色主题的查看器下黑线配黑底什么都看不见
    s.push_str(&format!(
        "<rect width=\"{:.2}\" height=\"{:.2}\" fill=\"#fff\"/>\n",
        scene.width, scene.height
    ));

    for item in &scene.items {
        match item {
            Primitive::Line { from, to, width } => {
                s.push_str(&format!(
                    "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" \
                     stroke=\"#000\" stroke-width=\"{:.2}\" stroke-linecap=\"round\"/>\n",
                    from.x, from.y, to.x, to.y, width
                ));
            }
            Primitive::Wedge { from, to, wide } => {
                // 实心三角:窄端在立体中心(一个点),宽端张开
                let d = (*to - *from).normalized();
                let n = crate::geom::Point2::new(-d.y, d.x) * (wide / 2.0);
                s.push_str(&format!(
                    "<path d=\"M{:.2},{:.2} L{:.2},{:.2} L{:.2},{:.2} Z\" fill=\"#000\"/>\n",
                    from.x,
                    from.y,
                    to.x + n.x,
                    to.y + n.y,
                    to.x - n.x,
                    to.y - n.y
                ));
            }
            Primitive::Hash {
                from,
                to,
                wide,
                spacing,
                width,
            } => {
                // 一叠垂直于键的短横线,从窄到宽。间距由规范的 hash_spacing 定。
                let len = from.dist(*to);
                // **`spacing` 为 0 时 `len / 0` 是 `inf`,`inf as i32` 饱和成
                // `i32::MAX`** —— 底下那个循环会画二十亿条线,不是报错,是把内存
                // 吃光。规范里给 0 是调用方的编程错误,但错在这里的表现太难查了,
                // 所以钳一道:非正的间距退回"两条线",与 `.max(2)` 同一条兜底。
                let n_lines = if *spacing > 0.0 {
                    ((len / spacing).floor() as i32).max(2)
                } else {
                    2
                };
                let d = (*to - *from).normalized();
                let perp = crate::geom::Point2::new(-d.y, d.x);
                for k in 1..=n_lines {
                    let t = f64::from(k) / f64::from(n_lines);
                    let c = *from + d * (len * t);
                    let h = perp * (wide / 2.0 * t);
                    s.push_str(&format!(
                        "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" \
                         stroke=\"#000\" stroke-width=\"{:.2}\" stroke-linecap=\"round\"/>\n",
                        c.x - h.x,
                        c.y - h.y,
                        c.x + h.x,
                        c.y + h.y,
                        width
                    ));
                }
            }
            Primitive::Text { at, runs, size } => {
                s.push_str(&format!(
                    "<text x=\"{:.2}\" y=\"{:.2}\" font-family=\"{}\" font-size=\"{:.2}\" \
                     text-anchor=\"middle\" dominant-baseline=\"central\" fill=\"#000\">",
                    at.x,
                    at.y,
                    // **字体名要转义。** 它是 `&'static str`,通常来自本仓的规范表,
                    // 但 `Style` 是公开可构造的 —— 一个带 `"` 或 `&` 的字体名直插
                    // 属性里会把整个 `<text>` 元素写坏,而写出来的 SVG 看着像是
                    // 渲染的毛病。文本内容那一侧早就在转义了,属性这一侧漏了。
                    escape(style.font_family),
                    size
                ));
                // 每一段自带 `dy`,当前基线偏移显式记着。
                //
                // **不能用空的 `<tspan dy=…>` 复位。** SVG 的 `dy` 是加在这一段
                // 的**字符**上的,段里没有字符就没有东西可加,基线复不回去 ——
                // 于是下标后面的正文一直吊在下标的高度上。实测:氨基在右侧写成
                // `H₂N` 时,`N` 比 `H` 低了一截,而这两个字母本该在同一条基线上。
                //
                // 抬升/下沉与字号缩放用的是 [`label`](crate::label) 里的同一组
                // 常数 —— 那边按它们算标签占多大,这边按它们画,两边必须一致。
                let mut cur = 0.0_f64;
                for r in runs {
                    let (text, want, fs) = match r {
                        Run::Normal(t) => (t, 0.0, *size),
                        Run::Sub(t) => (
                            t,
                            size * crate::label::SUB_DROP,
                            size * crate::label::SUB_SUP_SCALE,
                        ),
                        Run::Sup(t) => (
                            t,
                            -size * crate::label::SUP_RISE,
                            size * crate::label::SUB_SUP_SCALE,
                        ),
                    };
                    s.push_str(&format!(
                        "<tspan font-size=\"{fs:.2}\" dy=\"{:.2}\">{}</tspan>",
                        want - cur,
                        escape(text)
                    ));
                    cur = want;
                }
                s.push_str("</text>\n");
            }
        }
    }
    s.push_str("</svg>\n");
    s
}

/// XML 转义。
///
/// 元素符号与数字里没有这些字符,但**电荷用的 `+`/`-` 之外将来可能出现别的记号**,
/// 而漏转义一次就会产出坏 XML —— 那种文件有的查看器打得开、有的打不开,极难定位。
fn escape(t: &str) -> String {
    let mut out = String::with_capacity(t.len());
    for c in t.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            _ => out.push(c),
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{generate, render::scene, style::Style};
    use omgkit_core::MolBuilder;

    fn prep(smi: &str) -> MolBuilder {
        let mut m = omgkit_io::smiles::parse(smi).unwrap();
        omgkit_chem::pipeline::sanitize(&mut m).unwrap();
        m
    }

    fn svg(smi: &str, style: &Style) -> String {
        let m = prep(smi);
        to_svg(&scene(&m, &generate(&m, style), style), style)
    }

    #[test]
    fn the_output_is_well_formed_xml() {
        // 坏 XML 有的查看器打得开、有的打不开,极难定位。这里做最基本的结构检查:
        // 标签成对、属性引号配对。
        for smi in [
            "c1ccccc1",
            "CC(=O)Oc1ccccc1C(=O)O",
            "[NH4+]",
            "[13CH4]",
            "CC#N",
        ] {
            let s = svg(smi, &Style::ACS_1996);
            assert!(s.starts_with("<svg "), "{smi} 开头不对");
            assert!(s.trim_end().ends_with("</svg>"), "{smi} 结尾不对");
            assert_eq!(
                s.matches("<text").count(),
                s.matches("</text>").count(),
                "{smi} 的 <text> 没有成对"
            );
            assert_eq!(
                s.matches("<tspan").count(),
                s.matches("</tspan>").count(),
                "{smi} 的 <tspan> 没有成对"
            );
            assert_eq!(s.matches('"').count() % 2, 0, "{smi} 的属性引号数是奇数");
            assert!(!s.contains("NaN"), "{smi} 里出现了 NaN 坐标");
        }
    }

    #[test]
    fn labels_appear_and_skeleton_carbons_do_not() {
        let s = svg("CCO", &Style::ACS_1996);
        assert!(s.contains(">OH<") || s.contains(">O<"), "羟基没画出来:{s}");
        // 骨架碳不该有文字
        assert_eq!(s.matches("<text").count(), 1, "只该有羟基一个标签");
    }

    #[test]
    fn charges_and_isotopes_are_escaped_and_marked_up() {
        let s = svg("[NH4+]", &Style::ACS_1996);
        assert!(s.contains("<tspan"), "电荷应当是上标");
        assert!(s.contains('+'), "电荷符号没画出来");

        let iso = svg("[13CH4]", &Style::ACS_1996);
        assert!(iso.contains("13"), "同位素没画出来");
    }

    #[test]
    fn the_canvas_scales_with_the_style() {
        // 同一个分子,ChemDraw 默认规范的画布应当比 ACS 大近一倍
        let a = svg("c1ccc2ccccc2c1", &Style::ACS_1996);
        let c = svg("c1ccc2ccccc2c1", &Style::CHEMDRAW_DEFAULT);
        let w = |s: &str| {
            s.split("width=\"")
                .nth(1)
                .unwrap()
                .split('"')
                .next()
                .unwrap()
                .parse::<f64>()
                .unwrap()
        };
        assert!(w(&c) > w(&a) * 1.8, "画布宽度比只有 {:.2}", w(&c) / w(&a));
    }

    #[test]
    fn xml_special_characters_never_leak_through() {
        // 元素符号里现在没有这些字符,但漏转义一次就产出坏 XML。这条守的是
        // 那个函数本身,不是当下的元素表。
        assert_eq!(escape("a<b>c&d\"e'f"), "a&lt;b&gt;c&amp;d&quot;e&apos;f");
    }

    /// **属性里的字体名也要转义,间距为 0 不许把内存吃光。**
    ///
    /// 两条都是"公开可构造的 `Style` 里一个古怪的值",而两条的表现都不是报错:
    /// 前者产出坏 XML(有的查看器打得开、有的打不开),后者 `len / 0.0` 是 `inf`,
    /// `inf as i32` 饱和成 `i32::MAX`,循环画二十亿条线。
    #[test]
    fn a_hostile_style_does_not_produce_broken_xml_or_eat_the_memory() {
        let mut style = Style::ACS_1996;
        style.font_family = r#"Ari"al & <script>"#;
        style.hash_spacing_pt = 0.0;
        let svg = svg("C[C@H](N)C(=O)O", &style);
        assert!(
            !svg.contains(r#"font-family="Ari"al"#),
            "字体名没转义,`\"` 把属性截断了"
        );
        assert!(svg.contains("&amp;") && svg.contains("&lt;script&gt;"));
        assert!(svg.len() < 1_000_000, "间距为 0 画出了 {} 字节", svg.len());
    }

    #[test]
    fn a_stereocentre_is_actually_drawn_with_a_wedge() {
        // 楔形只存在 `Depiction` 里而没画出来的话,图上一个手性中心都看不出来 ——
        // 而所有"数得出来"的判据(原子数、键长、重合)照样全绿。
        let s = svg("N[C@@H](C)O", &Style::ACS_1996);
        assert!(
            s.contains("<path") || s.matches("<line").count() > 3,
            "既没有实楔形的三角也没有虚楔形的横线堆:{s}"
        );
    }

    #[test]
    fn the_two_enantiomers_do_not_produce_the_same_svg() {
        // **区分力**:对映体的坐标完全全等,唯一的差别在楔形。两张 SVG 一样,
        // 就说明楔形要么没画、要么画反了 —— 图上分不出左右手。
        let a = svg("N[C@@H](C)O", &Style::ACS_1996);
        let b = svg("N[C@H](C)O", &Style::ACS_1996);
        assert_ne!(a, b, "两个对映体画出了完全相同的图");
    }

    #[test]
    fn normal_text_stays_on_the_main_baseline() {
        // **下标之后的正文必须回到主基线上。**
        //
        // 先前是用一个空的 `<tspan dy=…>` 复位的,而 SVG 的 `dy` 是加在这一段
        // 的**字符**上的 —— 段里没有字符就没有东西可加,基线复不回去,于是
        // 下标后面的正文一直吊在下标那个高度。实测:氨基在右侧写成 `H₂N` 时,
        // `N` 比 `H` 低了一截,而这两个字母本该在同一条基线上。
        //
        // 这条把每一段的 `dy` 累加起来,要求正文段落落在偏移 0 上。
        for smi in [
            "NCc1ccccc1", // 苄胺:键从右边来,氨基写成 H₂N,下标夹在中间
            "CC(=O)Nc1ccc(O)cc1",
            "[NH4+]",  // 上标在末尾
            "[13CH4]", // 同位素在开头
            "OS(=O)(=O)O",
        ] {
            let out = svg(smi, &Style::ACS_1996);
            assert!(
                !out.contains("></tspan>"),
                "{smi}:出现了空的 <tspan> —— 它的 dy 不会生效"
            );
            for t in out.split("<text ").skip(1) {
                let body = t.split('>').skip(1).collect::<Vec<_>>().join(">");
                let body = body.split("</text>").next().expect("有结束标签");
                let mut cur = 0.0_f64;
                for seg in body.split("<tspan ").skip(1) {
                    let dy: f64 = seg
                        .split("dy=\"")
                        .nth(1)
                        .and_then(|x| x.split('"').next())
                        .and_then(|x| x.parse().ok())
                        .expect("每段都要有 dy");
                    cur += dy;
                    let fs: f64 = seg
                        .split("font-size=\"")
                        .nth(1)
                        .and_then(|x| x.split('"').next())
                        .and_then(|x| x.parse().ok())
                        .expect("每段都要有 font-size");
                    let text = seg
                        .split('>')
                        .nth(1)
                        .and_then(|x| x.split('<').next())
                        .unwrap_or("");
                    if text.is_empty() {
                        continue;
                    }
                    // 正文段(没缩小字号的)必须在主基线上
                    if (fs - Style::ACS_1996.atom_label_pt).abs() < 1e-6 {
                        assert!(
                            cur.abs() < 1e-6,
                            "{smi}:正文 {text:?} 画在了偏移 {cur:.2} 上,不在主基线"
                        );
                    } else {
                        assert!(cur.abs() > 1e-6, "{smi}:上下标 {text:?} 却没有偏移");
                    }
                }
            }
        }
    }

    #[test]
    fn there_is_a_white_background() {
        // 没有底色的话,深色主题查看器下黑线配黑底,什么都看不见
        let s = svg("CCO", &Style::ACS_1996);
        assert!(s.contains("fill=\"#fff\""), "缺少白底");
    }
}
