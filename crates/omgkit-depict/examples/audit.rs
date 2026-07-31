//! 拿整份语料过一遍硬性质,把违例逐条报出来。
//!
//! ```shell
//! cargo run -p omgkit-depict --release --example audit -- harness/corpus/large.smi
//! ```
//!
//! # 为什么单元判据不够
//!
//! 单元判据里的分子是**手挑的**,而手挑的时候心里已经有一个模型 —— 挑出来的
//! 正好是那个模型覆盖得到的。先前"环内双键不许画到环外"那条列了八个分子,
//! 全是单环或邻稠的,于是桥环共用键那一整类漏了出去;"楔形窄端"那条列的三个
//! 分子都不含 P(V) 中心,于是"楔形落到双键上被静默丢弃"漏了出去。
//!
//! 全量过一遍不需要预先知道会错在哪。
//!
//! # 查哪几条
//!
//! | | 性质 | 前提 |
//! |---|---|---|
//! | 写法无关 | 换一种 SMILES 写法,画出来的图元一模一样 | 无 |
//! | 环内双键 | 环上双键的两条线都落在某个含它的环里 | 布局没退化 |
//! | 楔形落地 | `Depiction` 里记了几个楔形,画布上就有几个 | 无 |
//! | 楔形可读 | 画出楔形的中心,反读回来就是它该有的构型 | 无 |
//! | 键长全等 | 所有键画出来一样长 | 布局没退化 |
//! | 不出画布 | 每个图元都在画布里 | 无 |
//!
//! 带"布局没退化"前提的两条,是因为退化的坐标本身就不成形状 —— 那时要求图
//! 画得对没有意义,而**退化这件事已经报在 [`Depiction::degraded`] 里了**。

use std::collections::BTreeMap;

use omgkit_core::MolBuilder;
use omgkit_depict::{
    generate,
    geom::{point_in_polygon, Point2},
    render::{scene, Primitive, Scene},
    style::Style,
};

/// 换几种写法比对。取 3 —— 再多也只是重复同一类扰动。
const WRITINGS: usize = 3;

fn main() {
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "harness/corpus/large.smi".into());
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("读不了 {path}:{e}"));

    let mut fails: BTreeMap<&'static str, Vec<String>> = BTreeMap::new();
    let mut n_ok = 0usize;
    let mut n_skip = 0usize;
    let mut checked: BTreeMap<&'static str, usize> = BTreeMap::new();
    let mut quality: BTreeMap<&'static str, usize> = BTreeMap::new();

    for (lineno, line) in text.lines().enumerate() {
        let smi = line.split_whitespace().next().unwrap_or("");
        if smi.is_empty() {
            continue;
        }
        let Some(m) = prep(smi) else {
            n_skip += 1;
            continue;
        };
        if m.num_atoms() < 2 {
            n_skip += 1;
            continue;
        }
        n_ok += 1;
        for style in &Style::ALL {
            let d = generate(&m, style);
            let s = scene(&m, &d, style);
            let tag = format!("{}:{lineno}:{smi}", style.name);
            let clean = d.degraded.is_empty() && d.unresolved.is_empty();
            // 画得干不干净也要报 —— 判据只说"没画错",不说"画得好"
            *quality
                .entry(if !d.degraded.is_empty() {
                    "退化(桥环等)"
                } else if !d.unresolved.is_empty() {
                    "有未解冲突"
                } else if !d.crossings.is_empty() {
                    "有键交叉"
                } else {
                    "干净"
                })
                .or_default() += 1usize;

            for (name, hit, bad) in checks(&m, &d, &s, style, clean) {
                *checked.entry(name).or_default() += usize::from(hit);
                if let Some(why) = bad {
                    fails
                        .entry(name)
                        .or_default()
                        .push(format!("{tag} —— {why}"));
                }
            }
        }
    }

    println!("语料 {path}:解析成功 {n_ok},跳过 {n_skip}\n");
    let tot: usize = quality.values().sum();
    println!("出图质量({tot} 个分子×规范):");
    for (k, v) in &quality {
        println!("  {k:<16} {v:>6}  {:>5.1}%", 100.0 * *v as f64 / tot as f64);
    }
    println!();
    println!("{:<14} {:>10} {:>8}", "性质", "查到", "违例");
    let mut total = 0usize;
    for (name, n) in &checked {
        let bad = fails.get(name).map_or(0, Vec::len);
        total += bad;
        println!("{name:<14} {n:>10} {bad:>8}");
    }
    for (name, list) in &fails {
        println!("\n=== {name} 的前几例 ===");
        for x in list.iter().take(300) {
            println!("  {x}");
        }
        if list.len() > 300 {
            println!("  …… 另有 {} 例", list.len() - 300);
        }
    }
    println!(
        "\n{}",
        if total == 0 {
            "全部通过".to_string()
        } else {
            format!("共 {total} 处违例")
        }
    );
    if total > 0 {
        std::process::exit(1);
    }
}

fn prep(smi: &str) -> Option<MolBuilder> {
    let mut m = omgkit_io::smiles::parse(smi).ok()?;
    omgkit_chem::pipeline::sanitize(&mut m).ok()?;
    omgkit_io::stereo::perceive_bond_stereo(&mut m);
    Some(m)
}

/// 形状指纹:两两距离排序后的多重集。与原子编号、平移、旋转、镜像都无关。
fn shape(c: &[Point2]) -> Vec<i64> {
    let mut v: Vec<i64> = (0..c.len())
        .flat_map(|i| ((i + 1)..c.len()).map(move |j| (i, j)))
        .map(|(i, j)| (c[i].dist(c[j]) * 1e4).round() as i64)
        .collect();
    v.sort_unstable();
    v
}

/// 量化后的坐标多重集 —— 与原子编号无关。
fn quantised(c: &[Point2]) -> Vec<(i64, i64)> {
    let mut v: Vec<(i64, i64)> = c
        .iter()
        .map(|p| ((p.x * 1e4).round() as i64, (p.y * 1e4).round() as i64))
        .collect();
    v.sort_unstable();
    v
}

/// 图元的多重集。线段不分方向 —— `L A B` 与 `L B A` 是同一条线,而谁是起点
/// 取决于键的 begin/end,本来就随写法变。楔形分方向:窄端宽端不是一回事。
fn fingerprint(s: &Scene) -> Vec<String> {
    let q = |p: Point2| format!("{:.3},{:.3}", p.x, p.y);
    let mut v: Vec<String> = s
        .items
        .iter()
        .map(|it| match it {
            Primitive::Line { from, to, .. } => {
                let (x, y) = (q(*from), q(*to));
                if x <= y {
                    format!("L {x} {y}")
                } else {
                    format!("L {y} {x}")
                }
            }
            Primitive::Wedge { from, to, .. } => format!("W {} {}", q(*from), q(*to)),
            Primitive::Hash { from, to, .. } => format!("H {} {}", q(*from), q(*to)),
            Primitive::Text { at, runs, .. } => format!("T {} {runs:?}", q(*at)),
        })
        .collect();
    v.sort();
    v
}

type Check = (&'static str, bool, Option<String>);

fn checks(
    m: &MolBuilder,
    d: &omgkit_depict::Depiction,
    s: &Scene,
    style: &Style,
    clean: bool,
) -> Vec<Check> {
    vec![
        writing_independent(m, d, s, style, clean),
        ring_double_bonds(m, d, s, style, clean),
        wedges_reach_canvas(d, s),
        wedges_read_back(m, d),
        bond_lengths_equal(m, d, clean),
        inside_canvas(s),
    ]
}

/// 换写法画出来必须一模一样。
fn writing_independent(
    m: &MolBuilder,
    d0: &omgkit_depict::Depiction,
    s: &Scene,
    style: &Style,
    clean: bool,
) -> Check {
    let want = fingerprint(s);
    let n = m.num_atoms();
    // 真正比过几次、因为"换出来不是同一个分子"跳过了几次 —— 一次都没比过的话
    // 这一条是空过的,必须看得见
    let (mut compared, mut skipped) = (0usize, 0usize);
    for k in 1..=WRITINGS {
        // 用一个与规范秩无关的优先序把写法搅乱。乘法哈希,足够打散又完全确定。
        let priority: Vec<u32> = (0..n)
            .map(|i| u32::try_from((i * 2_654_435_761 * k) % (n * 7 + 13)).unwrap_or(0))
            .collect();
        let w = omgkit_io::smiles::write_with_priority(m, &priority);
        let Some(m2) = prep(&w.smiles) else { continue };
        // 换出来的必须还是同一个分子,否则比的是两个东西。
        //
        // **要用规范式比。** `write` 按存储序写,而 m2 的存储序正是被打乱过的
        // 那一套 —— 拿它去比,几乎每个分子都会判成"不是同一个"而跳过,这一条
        // 就静悄悄地什么都没查。
        if omgkit_io::canon::canonical_smiles(m).smiles
            != omgkit_io::canon::canonical_smiles(&m2).smiles
        {
            skipped += 1;
            continue;
        }
        compared += 1;
        let d2 = generate(&m2, style);
        let got = fingerprint(&scene(&m2, &d2, style));
        if got != want {
            let diff = want.iter().zip(&got).filter(|(a, b)| a != b).count();
            // **差在哪一层要分开说。** 坐标就变了是布局的问题;坐标一样只是
            // 楔形或键级不同,那是指派的问题 —— 两者的修法完全不同。
            let coords_same = quantised(&d0.coords) == quantised(&d2.coords);
            let wedge_kinds = |d: &omgkit_depict::Depiction| {
                let mut v: Vec<String> = d
                    .wedges
                    .iter()
                    .filter_map(|w| {
                        w.narrow()
                            .map(|_| format!("{w:?}").split('{').next().unwrap_or("").to_string())
                    })
                    .collect();
                v.sort();
                v
            };
            let layer = if !coords_same {
                // 形状指纹:两两距离的多重集,与平移旋转镜像都无关。形状一样
                // 只是摆位不同,病在规范朝向;形状本身变了,病在布局。
                if shape(&d0.coords) == shape(&d2.coords) {
                    "形状相同,摆位不同"
                } else {
                    "形状就变了"
                }
            } else if wedge_kinds(d0) != wedge_kinds(&d2) {
                "坐标相同,楔形不同"
            } else {
                "坐标相同,键级或落点不同"
            };
            return (
                "写法无关",
                true,
                Some(format!(
                    "[{layer}|{}] 写成 {} 之后有 {diff}/{} 处图元不同",
                    if clean {
                        "布局干净"
                    } else {
                        "布局已退化"
                    },
                    w.smiles,
                    want.len()
                )),
            );
        }
    }
    if compared == 0 {
        return (
            "写法无关",
            false,
            Some(format!("{skipped} 种改写全都不是同一个分子,这条没查成")),
        );
    }
    ("写法无关", true, None)
}

/// 环上双键的两条线都要落在某个含它的环里。
fn ring_double_bonds(
    m: &MolBuilder,
    d: &omgkit_depict::Depiction,
    s: &Scene,
    style: &Style,
    clean: bool,
) -> Check {
    if !clean {
        return ("环内双键", false, None);
    }
    let pts = canvas_pts(m, d, style);
    let rings = omgkit_chem::sssr::ring_set(m);
    // **要按图上画成什么来筛,不是按分子里记的是什么** —— 芳香键在图上是交替
    // 单双,而单键根本没有第二条线可查
    let orders = omgkit_depict::render::drawn_orders(m);
    let mut hit = false;
    for (bi, b) in m.bonds().iter().enumerate() {
        let bond_no = u32::try_from(bi).expect("键数超出 u32");
        let mine: Vec<_> = rings
            .iter()
            .filter(|r| r.bonds.contains(&bond_no))
            .collect();
        if mine.is_empty() || orders[bi] != omgkit_core::BondOrder::Double {
            continue;
        }
        let (pa, pb) = (pts[b.begin as usize], pts[b.end as usize]);
        let len = pa.dist(pb);
        if len < 1e-9 {
            continue;
        }
        let mid = (pa + pb) * 0.5;
        let axis = (pb - pa) * (1.0 / len);
        let normal = Point2::new(-axis.y, axis.x);
        // 这根键画出来的线:两端都贴着键轴、落在跨度内
        let lines: Vec<(Point2, Point2)> = s
            .items
            .iter()
            .filter_map(|it| match it {
                Primitive::Line { from, to, .. } => Some((*from, *to)),
                _ => None,
            })
            .filter(|(f, t)| {
                [f, t].iter().all(|p| {
                    let v = **p - pa;
                    v.dot(normal).abs() < 0.30 * len
                        && v.dot(axis) > -0.10 * len
                        && v.dot(axis) < 1.10 * len
                })
            })
            .collect();
        if lines.len() != 2 {
            continue; // 拥挤处会多选到别的线,那时不下判断
        }
        hit = true;
        for (f, t) in lines {
            let lm = (f + t) * 0.5;
            if (lm - mid).dot(normal).abs() < 0.02 * len {
                continue; // 骑在轴上那条
            }
            let ok = mine.iter().any(|r| {
                let poly: Vec<Point2> = r.atoms.iter().map(|a| pts[*a as usize]).collect();
                point_in_polygon(lm, &poly)
            });
            if !ok {
                return (
                    "环内双键",
                    true,
                    Some(format!("键 {bi}({}–{})有一条线画在环外", b.begin, b.end)),
                );
            }
        }
    }
    ("环内双键", hit, None)
}

/// 记了几个楔形,画布上就该有几个。
fn wedges_reach_canvas(d: &omgkit_depict::Depiction, s: &Scene) -> Check {
    let recorded = d.wedges.iter().filter(|w| w.narrow().is_some()).count();
    let drawn = s
        .items
        .iter()
        .filter(|it| matches!(it, Primitive::Wedge { .. } | Primitive::Hash { .. }))
        .count();
    if recorded != drawn {
        return (
            "楔形落地",
            true,
            Some(format!("记了 {recorded} 个,画出来 {drawn} 个")),
        );
    }
    ("楔形落地", recorded > 0, None)
}

/// 画出楔形的中心,反读回来必须还是它该有的构型。
fn wedges_read_back(m: &MolBuilder, d: &omgkit_depict::Depiction) -> Check {
    let genuine = omgkit_io::stereo::genuine_tetrahedral(m);
    let mut hit = false;
    for (i, a) in m.atoms().iter().enumerate() {
        let at = u32::try_from(i).expect("原子数超出 u32");
        if !genuine[i]
            || !matches!(
                a.chiral_tag,
                omgkit_core::ChiralTag::Cw | omgkit_core::ChiralTag::Ccw
            )
            || d.unwedged.contains(&at)
        {
            continue;
        }
        hit = true;
        let got = omgkit_depict::stereo::read_chirality(m, &d.coords, &d.wedges, at);
        if got != Some(a.chiral_tag) {
            return (
                "楔形可读",
                true,
                Some(format!(
                    "中心 {at} 画出来了,反读是 {got:?},该是 {:?}",
                    a.chiral_tag
                )),
            );
        }
    }
    ("楔形可读", hit, None)
}

/// 所有键画出来一样长。
fn bond_lengths_equal(m: &MolBuilder, d: &omgkit_depict::Depiction, clean: bool) -> Check {
    if !clean || m.num_bonds() == 0 {
        return ("键长全等", false, None);
    }
    let first = d.coords[m.bonds()[0].begin as usize].dist(d.coords[m.bonds()[0].end as usize]);
    for (bi, b) in m.bonds().iter().enumerate() {
        let l = d.coords[b.begin as usize].dist(d.coords[b.end as usize]);
        if (l - first).abs() > 1e-6 {
            return (
                "键长全等",
                true,
                Some(format!("键 {bi} 长 {l:.4},第一根长 {first:.4}")),
            );
        }
    }
    ("键长全等", true, None)
}

/// 图元不许伸到画布外。
fn inside_canvas(s: &Scene) -> Check {
    for it in &s.items {
        let pts: Vec<(Point2, f64)> = match it {
            Primitive::Line { from, to, width } => vec![(*from, *width / 2.0), (*to, *width / 2.0)],
            Primitive::Wedge { from, to, wide } | Primitive::Hash { from, to, wide, .. } => {
                vec![(*from, *wide / 2.0), (*to, *wide / 2.0)]
            }
            Primitive::Text { at, size, .. } => vec![(*at, *size / 2.0)],
        };
        for (p, r) in pts {
            if p.x - r < -0.01 || p.x + r > s.width + 0.01 || p.y - r < -0.01 || p.y + r > s.height
            {
                return (
                    "不出画布",
                    true,
                    Some(format!(
                        "图元在 ({:.2},{:.2})±{r:.2},画布 {:.2}×{:.2}",
                        p.x, p.y, s.width, s.height
                    )),
                );
            }
        }
    }
    ("不出画布", true, None)
}

fn canvas_pts(m: &MolBuilder, d: &omgkit_depict::Depiction, style: &Style) -> Vec<Point2> {
    // `bounds` 是 render 里的私有函数,这里按同样的规则重算一遍:
    // 含标签的包围盒 → 平移缩放 → y 翻转。
    let scale = style.bond_length_pt;
    let (mut x0, mut y1) = (f64::MAX, f64::MIN);
    for (i, p) in d.coords.iter().enumerate() {
        let a = u32::try_from(i).expect("原子数超出 u32");
        let (hw, hh) =
            omgkit_depict::label::label_for(m, a, style, omgkit_depict::label::HSide::Right)
                .map_or((0.0, 0.0), |l| (l.half_w, l.half_h));
        x0 = x0.min(p.x - hw);
        y1 = y1.max(p.y + hh);
    }
    d.coords
        .iter()
        .map(|p| Point2::new((p.x - x0) * scale + 8.0, (y1 - p.y) * scale + 8.0))
        .collect()
}
