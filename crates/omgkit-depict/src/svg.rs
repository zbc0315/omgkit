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

use std::collections::BTreeSet;

use crate::geom::Point2;
use crate::label::Run;
use crate::render::{Primitive, Scene, LIGHT};
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

    // 三维图的球要一个径向渐变才立体。渐变只能定义在 `<defs>` 里,所以先扫一遍
    // 用到了哪些颜色 —— 一种颜色一个渐变,而不是一个球一个:一张空间填充图有
    // 几百个球,逐球定义会让文件里塞满几乎相同的 `<radialGradient>`。
    //
    // 用 `BTreeSet` 而不是 `HashSet`:`<defs>` 里的次序要与分子无关且可复现,
    // 哈希序做不到这一点。
    let mut balls: BTreeSet<[u8; 3]> = BTreeSet::new();
    let mut sticks: BTreeSet<Cyl> = BTreeSet::new();
    let mut fills: BTreeSet<Disc> = BTreeSet::new();
    for it in &scene.items {
        match it {
            Primitive::Ball { color, .. } => {
                balls.insert(*color);
            }
            Primitive::Stick {
                from,
                to,
                width,
                color,
            } => {
                sticks.insert(cyl_of(*from, *to, *width, *color));
            }
            Primitive::AromaticFill {
                centre,
                focus,
                radius,
                inner,
                outer,
                ..
            } => {
                fills.insert(disc_of(*centre, *focus, *radius, *inner, *outer));
            }
            _ => {}
        }
    }
    if !balls.is_empty() || !sticks.is_empty() || !fills.is_empty() {
        s.push_str("<defs>\n");
        for c in &balls {
            s.push_str(&sphere_gradient(*c));
        }
        for c in &sticks {
            s.push_str(&cylinder_gradient(*c));
        }
        for d in &fills {
            s.push_str(&disc_gradient(*d));
        }
        s.push_str("</defs>\n");
    }

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
            Primitive::Ball { at, r, color } => {
                // 描一圈暗边。没有它,两个同色的球贴在一起时看不出是两个 ——
                // 而空间填充图里同色球贴在一起是常态。
                s.push_str(&format!(
                    "<circle cx=\"{:.2}\" cy=\"{:.2}\" r=\"{:.2}\" fill=\"url(#{})\" \
                     stroke=\"{}\" stroke-width=\"{:.2}\"/>\n",
                    at.x,
                    at.y,
                    r,
                    gradient_id(*color),
                    hex(shade(*color, -RIM)),
                    (r * RIM_WIDTH).max(0.25)
                ));
            }
            Primitive::Stick {
                from,
                to,
                width,
                color,
            } => {
                // 圆头端帽不是装饰:圆柱的两端在投影上就是半个圆,而**棍状样式
                // 全靠这个圆头当原子球** —— 换成平头端帽,每个接头都会出现缺口。
                s.push_str(&format!(
                    "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" \
                     stroke=\"url(#{})\" stroke-width=\"{:.2}\" stroke-linecap=\"round\"/>\n",
                    from.x,
                    from.y,
                    to.x,
                    to.y,
                    cyl_id(cyl_of(*from, *to, *width, *color)),
                    width
                ));
            }
            Primitive::AromaticFill {
                poly,
                centre,
                focus,
                radius,
                inner,
                outer,
            } => {
                // 多边形按环的顶点画,**不描边**:描了边就会在键线旁边多出一圈
                // 与它错开半个线宽的色边。
                let pts: Vec<String> = poly
                    .iter()
                    .map(|p| format!("{:.2},{:.2}", p.x, p.y))
                    .collect();
                s.push_str(&format!(
                    "<polygon points=\"{}\" fill=\"url(#{})\"/>\n",
                    pts.join(" "),
                    disc_id(disc_of(*centre, *focus, *radius, *inner, *outer))
                ));
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

/// 球心那一点提亮多少(向白色插值的比例)。
const HIGHLIGHT: f64 = 0.55;
/// 球边缘压暗多少(向黑色插值的比例)。
const RIM: f64 = 0.45;
/// 描边宽度占球半径的比例。
const RIM_WIDTH: f64 = 0.04;

/// 把颜色向白(`t > 0`)或向黑(`t < 0`)插值 `|t|`。
fn shade(c: [u8; 3], t: f64) -> [u8; 3] {
    let mut out = [0u8; 3];
    for k in 0..3 {
        let v = f64::from(c[k]);
        let x = if t >= 0.0 {
            v + (255.0 - v) * t
        } else {
            v * (1.0 + t)
        };
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        {
            out[k] = x.round().clamp(0.0, 255.0) as u8;
        }
    }
    out
}

fn hex(c: [u8; 3]) -> String {
    format!("#{:02x}{:02x}{:02x}", c[0], c[1], c[2])
}

/// `<defs>` 里那个渐变的 id。取颜色本身,所以同一种颜色永远同一个 id ——
/// 与分子、与图元次序都无关。
fn gradient_id(c: [u8; 3]) -> String {
    format!("s{:02x}{:02x}{:02x}", c[0], c[1], c[2])
}

/// 一个球的径向渐变:高光偏在左上,边缘压暗。
///
/// 高光的位置(`fx`/`fy`)与光源方向是一回事 —— 全图共用一个方向,否则每个球
/// 看起来都被不同的灯照着。左上是各家分子软件的默认光位。
fn sphere_gradient(c: [u8; 3]) -> String {
    format!(
        "<radialGradient id=\"{}\" cx=\"0.5\" cy=\"0.5\" r=\"0.55\" \
         fx=\"0.32\" fy=\"0.30\">\
         <stop offset=\"0\" stop-color=\"{}\"/>\
         <stop offset=\"0.55\" stop-color=\"{}\"/>\
         <stop offset=\"1\" stop-color=\"{}\"/>\
         </radialGradient>\n",
        gradient_id(c),
        hex(shade(c, HIGHLIGHT)),
        hex(c),
        hex(shade(c, -RIM))
    )
}

/// 一根圆柱的横向渐变,用**量化过的整数**表示,好让共线的切片落进同一个渐变。
///
/// 三个字段:法向(垂直于圆柱轴、指向光源那一侧)、这条轴到原点的有向距离、
/// 半宽,外加颜色。**不含圆柱的位置沿轴的分量** —— 一根键被切成六片之后,
/// 六片的中点各不相同,而横向渐变只沿法向变化,六片该拿同一个渐变。按中点存
/// 的话会发出六个几乎相同的 `<defs>` 条目,而且片与片的接缝上会差最后一位。
type Cyl = ([i64; 2], i64, i64, [u8; 3]);

/// 渐变量化的精度(每磅多少格)。
const GRAD_QUANT: f64 = 1000.0;

fn qg(x: f64) -> i64 {
    let v = (x * GRAD_QUANT).round();
    if v > i64::MAX as f64 {
        i64::MAX
    } else if v < i64::MIN as f64 {
        i64::MIN
    } else {
        v as i64
    }
}

fn cyl_of(from: Point2, to: Point2, width: f64, color: [u8; 3]) -> Cyl {
    let d = to - from;
    let len = (d.x * d.x + d.y * d.y).sqrt();
    // 长度为 0 的一段(键的两端重合)没有轴向,法向随便取一个定值 —— 画出来是
    // 一个圆点,渐变往哪边都一样,但**必须是确定的**,否则同一分子两次跑出来
    // 的 `<defs>` 不一样。
    let mut n = if len > f64::EPSILON {
        [-d.y / len, d.x / len]
    } else {
        [1.0, 0.0]
    };
    // 法向指向光源那一侧,于是高光永远在同一边。正好垂直于光时取"x 分量为正"
    // 这个确定的分支。
    let lit = n[0] * LIGHT[0] + n[1] * LIGHT[1];
    if lit < 0.0 || (lit == 0.0 && n[0] < 0.0) {
        n = [-n[0], -n[1]];
    }
    let offset = from.x * n[0] + from.y * n[1];
    ([qg(n[0]), qg(n[1])], qg(offset), qg(width / 2.0), color)
}

fn cyl_id(c: Cyl) -> String {
    let ([nx, ny], off, w, col) = c;
    format!(
        "c{nx}_{ny}_{off}_{w}_{:02x}{:02x}{:02x}",
        col[0], col[1], col[2]
    )
    .replace('-', "m")
}

/// 圆柱的横向渐变:亮边在光源那一侧,另一侧压暗。
///
/// 这不只是好看。**纯白的氢在白底上没有渐变就是隐形的** —— 第一版画出来,
/// 每根 C–H 键都像是只画了一半(靠碳那半是灰的,靠氢那半是白的,和背景一个色)。
/// 球有一圈暗描边所以看得见,棍没有。给圆柱上渐变把这一档补上,同时让球棍图
/// 真的像是三维的。
///
/// 描边解决不了这件事:一根键被切成若干片(见
/// [`DEPTH_SLICE`](crate::three::DEPTH_SLICE)),而切片是按深度分开画的,
/// 给每一片描边会在片
/// 与片的接缝上留下一道暗痕。渐变只沿法向变化,共线的片拿到同一个渐变,
/// 接缝上一个像素都不差。
fn cylinder_gradient(c: Cyl) -> String {
    let ([nx, ny], off, w, col) = c;
    let (nx, ny, off, w) = (
        nx as f64 / GRAD_QUANT,
        ny as f64 / GRAD_QUANT,
        off as f64 / GRAD_QUANT,
        w as f64 / GRAD_QUANT,
    );
    // 渐变向量:沿法向,从暗的一侧到亮的一侧,长度是圆柱的直径。
    // 起点取"轴线往负法向让开半宽"的那一点 —— 它只由 (法向, 有向距离) 决定,
    // 与这一片在轴上的位置无关,所以共线的片给出**逐字节相同**的定义。
    let (x1, y1) = (nx * (off - w), ny * (off - w));
    let (x2, y2) = (nx * (off + w), ny * (off + w));
    format!(
        "<linearGradient id=\"{}\" gradientUnits=\"userSpaceOnUse\" \
         x1=\"{x1:.3}\" y1=\"{y1:.3}\" x2=\"{x2:.3}\" y2=\"{y2:.3}\">\
         <stop offset=\"0\" stop-color=\"{}\"/>\
         <stop offset=\"0.4\" stop-color=\"{}\"/>\
         <stop offset=\"0.75\" stop-color=\"{}\"/>\
         <stop offset=\"1\" stop-color=\"{}\"/>\
         </linearGradient>\n",
        cyl_id(c),
        hex(shade(col, -CYL_RIM)),
        hex(col),
        hex(shade(col, CYL_HIGHLIGHT)),
        hex(shade(col, -CYL_RIM))
    )
}

/// 圆柱背光那一侧压暗多少。
const CYL_RIM: f64 = 0.42;
/// 圆柱迎光那一条带提亮多少。
const CYL_HIGHLIGHT: f64 = 0.42;

/// 一块芳香环底色的径向渐变,用**量化过的整数**表示,好让同样几何的环
/// (一张图里的几个苯环)共用一个 `<defs>` 条目。
///
/// 五个字段:圆心、焦点、半径、中心色、外缘色。**位置进了 id** —— 与球那边
/// 不同:球的渐变用的是 `objectBoundingBox`(每个球自己的包围盒),同色的球
/// 无论画在哪都能共用一个;这里用的是 `userSpaceOnUse`(画布坐标),因为
/// 多边形不是圆,包围盒的比例随环的朝向变,按包围盒定的焦点会跟着歪。
type Disc = ([i64; 2], [i64; 2], i64, [u8; 3], [u8; 3]);

fn disc_of(centre: Point2, focus: Point2, radius: f64, inner: [u8; 3], outer: [u8; 3]) -> Disc {
    (
        [qg(centre.x), qg(centre.y)],
        [qg(focus.x), qg(focus.y)],
        qg(radius),
        inner,
        outer,
    )
}

fn disc_id(d: Disc) -> String {
    let ([cx, cy], [fx, fy], r, i, o) = d;
    format!(
        "r{cx}_{cy}_{fx}_{fy}_{r}_{:02x}{:02x}{:02x}_{:02x}{:02x}{:02x}",
        i[0], i[1], i[2], o[0], o[1], o[2]
    )
    .replace('-', "m")
}

/// 芳香环底色:焦点处是中心色,到外接圆上变成外缘色。
///
/// `r` 取环的外接圆半径,所以**顶点上正好是纯外缘色** —— 环的六个角颜色一致,
/// 只有高光那一侧亮起来。半径给小了角上会出现一圈突兀的截断,给大了整块都是
/// 中心色、看不出渐变。
fn disc_gradient(d: Disc) -> String {
    let ([cx, cy], [fx, fy], r, i, o) = d;
    let g = |v: i64| v as f64 / GRAD_QUANT;
    format!(
        "<radialGradient id=\"{}\" gradientUnits=\"userSpaceOnUse\" \
         cx=\"{:.3}\" cy=\"{:.3}\" r=\"{:.3}\" fx=\"{:.3}\" fy=\"{:.3}\">\
         <stop offset=\"0\" stop-color=\"{}\"/>\
         <stop offset=\"1\" stop-color=\"{}\"/>\
         </radialGradient>\n",
        disc_id(d),
        g(cx),
        g(cy),
        g(r),
        g(fx),
        g(fy),
        hex(i),
        hex(o)
    )
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

    /// **纯白的原子在白底上必须看得见。**
    ///
    /// 氢是 CPK 里的纯白 `#ffffff`,而画布也是纯白 —— 球没有暗描边、棍没有
    /// 渐变的话,那个原子和它的键在图上**根本不存在**,而 SVG 本身一点毛病没有。
    /// 这不是审美问题:读图的人会以为那个位置是空的。
    ///
    /// 实测(变异标定):把球的描边从暗色改成球本身的颜色,全语料判官
    /// `harness/check_depict3d.py` 一声不吭地全绿 —— 那条判官读的是圆心、
    /// 半径、颜色、画序,看不见"这个圆在背景上分不分得出来"。这一条补的就是
    /// 那个洞,球与圆柱两边各钉一下。
    #[test]
    fn 纯白的原子与键在白底上分得出来() {
        use crate::three::{depict, Style3D};
        let mut m = omgkit_io::smiles::parse("CO").unwrap();
        let c = omgkit_conf::pipeline::conformer_for(&mut m).unwrap();

        // 球:白球必须描一圈与白底分得开的边
        let d = depict(&m, &c.coords, &Style3D::SPACE_FILLING).unwrap();
        let svg = to_svg(&d.scene, &Style::ACS_1996);
        let white_balls: Vec<&str> = svg
            .lines()
            .filter(|l| l.contains("<circle") && l.contains("url(#sffffff)"))
            .collect();
        assert!(!white_balls.is_empty(), "甲醇该有白色的氢球,判据没东西可判");
        for line in white_balls {
            let stroke = line
                .split("stroke=\"")
                .nth(1)
                .and_then(|t| t.split('"').next())
                .expect("球该有描边");
            assert!(
                stroke != "#ffffff" && stroke != "#fff",
                "白球的描边也是白的,在白底上看不见了:{line}"
            );
        }

        // 棍:白色圆柱的渐变必须有一档不是白的
        let d = depict(&m, &c.coords, &Style3D::STICK).unwrap();
        let svg = to_svg(&d.scene, &Style::ACS_1996);
        let grads: Vec<&str> = svg
            .lines()
            .filter(|l| l.contains("<linearGradient") && l.contains("#ffffff"))
            .collect();
        assert!(!grads.is_empty(), "甲醇该有白色的 C–H 圆柱,判据没东西可判");
        for g in grads {
            let stops: Vec<&str> = g
                .split("stop-color=\"")
                .skip(1)
                .filter_map(|t| t.split('"').next())
                .collect();
            assert!(stops.len() >= 2, "白圆柱的渐变只有 {} 档:{g}", stops.len());
            assert!(
                stops.iter().any(|c| *c != "#ffffff"),
                "白圆柱的渐变从头到尾都是白的,在白底上看不见了:{g}"
            );
        }
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

    /// 一套开着底色的规范。
    fn filled(fill: crate::style::AromaticFill) -> Style {
        Style {
            aromatic_fill: Some(fill),
            ..Style::ACS_1996
        }
    }

    /// **底色在 SVG 里必须排在所有线条、楔形、文字之前。**
    ///
    /// `Scene` 那一侧已经有一条判据钉了图元次序,这一条钉的是序列化没把它
    /// 打乱 —— 两侧各钉一下:图元排对了而写出来的次序反了,图上就是一块蓝色
    /// 盖住整个环。
    #[test]
    fn 底色写在所有线条与文字之前() {
        let st = filled(crate::style::AromaticFill::DEFAULT);
        for smi in [
            "c1ccc2ccccc2c1",
            "CC(=O)Oc1ccccc1C(=O)O",
            "c1ccc2[nH]ccc2c1",
        ] {
            let out = svg(smi, &st);
            let last_fill = out.rfind("<polygon").expect("该有底色");
            for tag in ["<line", "<path", "<text"] {
                if let Some(first) = out.find(tag) {
                    assert!(
                        last_fill < first,
                        "{smi}:最后一块底色在 {last_fill},而第一个 {tag} 在 {first}"
                    );
                }
            }
            assert!(out.contains("<line"), "{smi} 一根线都没有,这一档在空过");
        }
    }

    /// **关着就一个字节都不多。** 默认规范写出来的 SVG 里既没有多边形,
    /// 也没有多出来的渐变定义。
    #[test]
    fn 不开底色就什么都不多() {
        for st in &Style::ALL {
            let out = svg("c1ccc2ccccc2c1", st);
            assert!(!out.contains("<polygon"), "{} 写出了多边形", st.name);
            assert!(!out.contains("<defs"), "{} 写出了空的 defs", st.name);
        }
    }

    /// **两个颜色真的换得动,而且换的是对的那一头。**
    ///
    /// 只查"自定义色出现在文件里"是不够的:两个 stop 写反了照样绿。这里
    /// 按 `offset` 分别查 —— 焦点(offset 0)是中心色,外缘(offset 1)是外缘色。
    #[test]
    fn 底色的两个颜色各就各位() {
        let mine = crate::style::AromaticFill {
            centre: [0xff, 0xfb, 0xe6],
            edge: [0xf5, 0xc2, 0x6b],
        };
        let out = svg("c1ccccc1", &filled(mine));
        let g = out
            .lines()
            .find(|l| l.contains("<radialGradient"))
            .expect("该有一个径向渐变");
        let stop = |off: &str| {
            g.split(&format!("<stop offset=\"{off}\" stop-color=\""))
                .nth(1)
                .and_then(|t| t.split('"').next())
                .unwrap_or_else(|| panic!("渐变里没有 offset={off} 这一档:{g}"))
                .to_string()
        };
        assert_eq!(stop("0"), "#fffbe6", "焦点那一档不是中心色");
        assert_eq!(stop("1"), "#f5c26b", "外缘那一档不是外缘色");
        // 默认色一个都不该出现 —— 参数没接上时最容易出现的表现就是"还是默认色"
        let dflt = svg("c1ccccc1", &filled(crate::style::AromaticFill::DEFAULT));
        assert!(dflt.contains("#add8e6"), "默认外缘色该是 CSS 的 lightblue");
        assert!(!out.contains("#add8e6"), "自定义配色里还留着默认的浅蓝");
    }

    /// **一张图里几个一模一样的环共用一个渐变定义。**
    ///
    /// 逐块定义的话,联苯那种分子会在 `<defs>` 里塞进两条几乎相同的渐变;
    /// 而"几乎相同"是因为位置不同 —— 位置进了 id,所以**位置不同就不该共用**,
    /// 两个方向都要查。
    #[test]
    fn 同样的环共用一个渐变位置不同的不共用() {
        let st = filled(crate::style::AromaticFill::DEFAULT);
        let n = |smi: &str| svg(smi, &st).matches("<radialGradient").count();
        assert_eq!(n("c1ccccc1"), 1, "一个环该只有一条渐变");
        // 联苯的两个苯环形状相同、位置不同 —— 各一条
        assert_eq!(
            n("c1ccc(-c2ccccc2)cc1"),
            2,
            "两个位置不同的环该各有一条渐变"
        );
        assert_eq!(
            svg("c1ccc(-c2ccccc2)cc1", &st).matches("<polygon").count(),
            2,
            "两个环该各有一块多边形"
        );
    }

    #[test]
    fn there_is_a_white_background() {
        // 没有底色的话,深色主题查看器下黑线配黑底,什么都看不见
        let s = svg("CCO", &Style::ACS_1996);
        assert!(s.contains("fill=\"#fff\""), "缺少白底");
    }
}
