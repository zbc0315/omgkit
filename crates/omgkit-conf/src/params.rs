//! 键长与键角的**查表**,以及查不到时的兜底。
//!
//! # 数从哪儿来
//!
//! `data/mmff.bonds.tsv` 与 `data/mmff.angles.tsv`,是从 **8532 个分子的
//! MMFF94 收敛几何**里量出来的(口径:ETKDGv3 种子 `0xf00d` +
//! MMFF94 `maxIters=500`,只取收敛的;生成器 `harness/measure_params.py`)。
//!
//! - 键长键:(较小元素符号, 较大元素符号, 键级, **所在最小环尺寸**),不在环里记 0;
//! - 键角键:(中心元素, 配位数, 是否芳香, **中心原子自己的最小环**, 三原子共处的最小环)。
//!
//! 每行给 计数 / 中位 / p05 / p95 / 均值。
//!
//! # 为什么 p05 / p95 也要
//!
//! 二期环系的距离几何用的是**窄区间**而不是等式(见 crate 文档):
//! 角的上下界直接取 p05 / p95。等式化会让约束系统在规格上必然自相矛盾 ——
//! 实测三元环内角"须恰好和为 180°"94 个全不成立、四配位中心只有 0.29% 精确可实现,
//! 而且 1.51% 的分子在进嵌入前界矩阵就崩。
//!
//! # 查不到怎么办
//!
//! **逐级放宽,而且每一级都说得出自己走的是哪一级**([`Source`])——
//! 判据那边要按级别计数,不能让"查不到"悄悄混进"查到了"里。

use omgkit_core::{element, BondOrder};
use std::collections::HashMap;
use std::sync::OnceLock;

/// 这个数是**怎么来的** —— 判据要按它分级计数。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Source {
    /// 表里逐项命中(含环尺寸那一维)。
    Table,
    /// 表里没有这个环尺寸,退到"不在环里"那一行。
    RingRelaxed,
    /// 表里根本没有 → 共价半径模型(键长)/ VSEPR 理想角(键角)。
    Model,
}

/// 一个查出来的参数:值 + 它是第几级来的。
#[derive(Debug, Clone, Copy)]
pub struct Param {
    /// 值(键长 Å / 键角弧度)。
    pub value: f64,
    /// 上下界(键角用;键长给 ±3%)。二期的窄区间用它。
    pub lo: f64,
    /// 见 [`Param::lo`]。
    pub hi: f64,
    /// 来源级别。
    pub source: Source,
}

/// 键级在表里的记号。表里只有 `1/2/3/ar` —— 配位键、四重键都没有,走兜底。
fn order_tag(o: BondOrder) -> Option<&'static str> {
    match o {
        BondOrder::Single => Some("1"),
        BondOrder::Double => Some("2"),
        BondOrder::Triple => Some("3"),
        BondOrder::Aromatic => Some("ar"),
        _ => None,
    }
}

/// 键级 → 共价半径模型的乘子。兜底用,数值取自 v1 的 `params.rs`。
fn order_factor(o: BondOrder) -> f64 {
    match o {
        BondOrder::Aromatic => 0.920_6,
        BondOrder::Double => 0.869_9,
        BondOrder::Triple | BondOrder::Quadruple => 0.783_3,
        _ => 1.0,
    }
}

/// 共价半径(Å)。
///
/// **f32 → f64 只在这里转一次。** core 的元素表是 `f32`,而本 crate 全程 `f64`;
/// 转换点集中在一处,免得同一个数在不同路径上转出不同的位。
#[must_use]
pub fn covalent_radius(z: u8) -> f64 {
    element::by_atomic_num(z).map_or(0.76, |e| f64::from(e.rcov))
}

/// 范德华半径(Å)。消撞判据用。
#[must_use]
pub fn vdw_radius(z: u8) -> f64 {
    element::by_atomic_num(z).map_or(1.70, |e| f64::from(e.rvdw))
}

type BondKey = (String, String, String, usize);
type AngleKey = (String, usize, u8, usize, usize);
/// 表里每行取出来的四个数:中位、p05、p95、均值。
type Row = (f64, f64, f64, f64);

fn bonds() -> &'static HashMap<BondKey, Row> {
    static T: OnceLock<HashMap<BondKey, Row>> = OnceLock::new();
    T.get_or_init(|| {
        let mut m = HashMap::new();
        for line in include_str!("../data/mmff.bonds.tsv").lines() {
            if line.starts_with('#') {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 9 {
                continue;
            }
            let (Ok(ring), Ok(med), Ok(p05), Ok(p95), Ok(mean)) = (
                f[3].parse::<usize>(),
                f[5].parse::<f64>(),
                f[6].parse::<f64>(),
                f[7].parse::<f64>(),
                f[8].parse::<f64>(),
            ) else {
                continue;
            };
            m.insert(
                (f[0].to_string(), f[1].to_string(), f[2].to_string(), ring),
                (med, p05, p95, mean),
            );
        }
        m
    })
}

fn angles() -> &'static HashMap<AngleKey, Row> {
    static T: OnceLock<HashMap<AngleKey, Row>> = OnceLock::new();
    T.get_or_init(|| {
        let mut m = HashMap::new();
        for line in include_str!("../data/mmff.angles.tsv").lines() {
            if line.starts_with('#') {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 10 {
                continue;
            }
            let (Ok(deg), Ok(ar), Ok(rs), Ok(rg), Ok(med), Ok(p05), Ok(p95), Ok(mean)) = (
                f[1].parse::<usize>(),
                f[2].parse::<u8>(),
                f[3].parse::<usize>(),
                f[4].parse::<usize>(),
                f[6].parse::<f64>(),
                f[7].parse::<f64>(),
                f[8].parse::<f64>(),
                f[9].parse::<f64>(),
            ) else {
                continue;
            };
            m.insert((f[0].to_string(), deg, ar, rs, rg), (med, p05, p95, mean));
        }
        m
    })
}

/// **一根键的目标长度**(Å),外加 ±3% 的窄区间。
///
/// `min_ring` 是这根键所在的最小环尺寸,不在环里传 0。
///
/// 逐级放宽:表里逐项命中 → 退到"不在环里"那一行 → 共价半径模型。
/// **永不返回 0 或非有限数。**
#[must_use]
pub fn bond_length(a: u8, b: u8, order: BondOrder, min_ring: usize) -> Param {
    // **查不到元素就直接走模型,不许拿碳的符号去查表。**
    //
    // 头一版回退成 `"C"`,于是一个不认识的元素会**查到碳那一行并报
    // `Source::Table`** —— 分级来源那套东西的全部意义就是"说得出自己走的哪一级",
    // 报错了级比查不到更坏。现在元素表 0~118 都有项,所以这条走不到,
    // 但回退的**方向**是反的,留着迟早咬人。
    let (Some(ea), Some(eb)) = (element::by_atomic_num(a), element::by_atomic_num(b)) else {
        return covalent_model(a, b, order);
    };
    let (sa, sb) = (ea.symbol, eb.symbol);
    // 表里的元素对按**符号字典序**排(见 measure_params.py)
    let (lo_s, hi_s) = if sa <= sb { (sa, sb) } else { (sb, sa) };
    if let Some(tag) = order_tag(order) {
        let t = bonds();
        for (ring, src) in [(min_ring, Source::Table), (0, Source::RingRelaxed)] {
            if ring == 0 && src == Source::RingRelaxed && min_ring == 0 {
                break; // 第一级就是 ring=0,别把它记成"放宽"
            }
            if let Some(&(med, _, _, _)) =
                t.get(&(lo_s.to_string(), hi_s.to_string(), tag.to_string(), ring))
            {
                return Param {
                    value: med,
                    lo: med * 0.97,
                    hi: med * 1.03,
                    source: src,
                };
            }
        }
    }
    covalent_model(a, b, order)
}

/// 查不到表时的键长兜底:共价半径之和 × 键级系数。**只有它会报 `Source::Model`。**
fn covalent_model(a: u8, b: u8, order: BondOrder) -> Param {
    let v = (covalent_radius(a) + covalent_radius(b)) * order_factor(order);
    Param {
        value: v,
        lo: v * 0.97,
        hi: v * 1.03,
        source: Source::Model,
    }
}

/// **一个键角的目标值**(弧度),外加 `[p05, p95]` 的窄区间。
///
/// 查不到就按配位数回退到 VSEPR 的理想角(2 → 180°、3 → 120°、
/// 4 → 109.47°、其余 → 109.47°),区间给 ±15°。
#[must_use]
pub fn angle(
    center: u8,
    degree: usize,
    aromatic: bool,
    ring_self: usize,
    ring_shared: usize,
) -> Param {
    // 同 `bond_length`:查不到元素就走兜底,不许拿碳的符号去查表
    let Some(ec) = element::by_atomic_num(center) else {
        return angle_model(center, degree);
    };
    let sym = ec.symbol;
    let ar = u8::from(aromatic);
    let t = angles();
    let mut tries: Vec<(AngleKey, Source)> = vec![(
        (sym.to_string(), degree, ar, ring_self, ring_shared),
        Source::Table,
    )];
    // 环那两维查不到就放宽(先放共处环,再放中心环)
    if ring_shared != 0 {
        tries.push((
            (sym.to_string(), degree, ar, ring_self, 0),
            Source::RingRelaxed,
        ));
    }
    if ring_self != 0 {
        tries.push(((sym.to_string(), degree, ar, 0, 0), Source::RingRelaxed));
    }
    for (k, src) in tries {
        if let Some(&(med, p05, p95, _)) = t.get(&k) {
            return Param {
                value: med.to_radians(),
                lo: p05.to_radians(),
                hi: p95.to_radians(),
                source: src,
            };
        }
    }
    angle_model(center, degree)
}

/// 查不到表时的键角兜底:按配位数回退到 VSEPR 的理想角。
fn angle_model(_center: u8, degree: usize) -> Param {
    let ideal: f64 = match degree {
        0..=2 => 180.0,
        3 => 120.0,
        _ => 109.471,
    };
    Param {
        value: ideal.to_radians(),
        lo: (ideal - 15.0).to_radians(),
        hi: (ideal + 15.0).to_radians(),
        source: Source::Model,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **crate 里那两份表必须与 `harness/params/` 里的逐字节相同。**
    ///
    /// 表在两个地方各存一份(crate 要 `include_str!`,harness 要给判据用),
    /// 所以必须有人盯着它们别漂。漂了的话生产用的是一份、判据用的是另一份,
    /// 而两边都"自洽",谁都不会红。
    #[test]
    fn the_embedded_tables_match_the_ones_in_harness() {
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(std::path::Path::parent)
            .expect("workspace 根");
        for name in ["mmff.bonds.tsv", "mmff.angles.tsv"] {
            let mine = std::fs::read_to_string(
                std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                    .join("data")
                    .join(name),
            )
            .expect("crate 里那份");
            let theirs = std::fs::read_to_string(root.join("harness/params").join(name))
                .expect("harness 里那份");
            assert_eq!(mine, theirs, "{name} 两份不一致 —— 生产与判据会用不同的数");
        }
    }

    /// 表真的读进来了,而且几条已知的行对得上。
    ///
    /// 这条防的是"解析写错了但每次都悄悄走兜底" —— 那种情况下所有几何都还能出,
    /// 只是全用的共价半径模型,而**判据看不出来**。
    #[test]
    fn the_tables_are_actually_parsed_and_hit() {
        assert!(bonds().len() > 100, "键长表只读到 {} 行", bonds().len());
        assert!(angles().len() > 150, "键角表只读到 {} 行", angles().len());
        // C–H 单键、不在环里:表里是 1.094
        let p = bond_length(6, 1, BondOrder::Single, 0);
        assert_eq!(p.source, Source::Table, "C–H 该查得到");
        assert!((p.value - 1.094).abs() < 1e-9, "C–H 得到 {}", p.value);
        // 芳香 C–C 在六元环里:1.397
        let p = bond_length(6, 6, BondOrder::Aromatic, 6);
        assert_eq!(p.source, Source::Table);
        assert!((p.value - 1.397).abs() < 1e-9, "芳香 C–C 得到 {}", p.value);
        // 四配位碳、不在环:109.4°
        let a = angle(6, 4, false, 0, 0);
        assert_eq!(a.source, Source::Table);
        assert!(
            (a.value.to_degrees() - 109.4).abs() < 1e-6,
            "sp³ 碳得到 {}",
            a.value.to_degrees()
        );
        assert!(a.lo < a.value && a.value < a.hi, "p05 < 中位 < p95");
    }

    /// 查不到的必须**说自己查不到**,并且给一个有限的正数。
    #[test]
    fn a_miss_says_so_and_still_returns_something_usable() {
        // 配位键在表里没有记号
        let p = bond_length(6, 8, BondOrder::Dative, 0);
        assert_eq!(p.source, Source::Model);
        assert!(p.value.is_finite() && p.value > 0.5, "{}", p.value);
        // 钨没在语料里出现过
        let p = bond_length(74, 74, BondOrder::Single, 0);
        assert_eq!(p.source, Source::Model);
        assert!(p.value.is_finite() && p.value > 0.5);
        // 配位数 7 的中心角表里没有
        let a = angle(6, 7, false, 0, 0);
        assert_eq!(a.source, Source::Model);
        assert!(a.value.is_finite() && a.value > 0.0);
    }

    /// 元素对的顺序不能影响结果 —— 表里是按符号字典序存的。
    #[test]
    fn bond_lookup_does_not_care_which_atom_comes_first() {
        for (a, b) in [(6u8, 1u8), (7, 6), (8, 6), (16, 6), (17, 6)] {
            let x = bond_length(a, b, BondOrder::Single, 0);
            let y = bond_length(b, a, BondOrder::Single, 0);
            assert!(
                (x.value - y.value).abs() < 1e-15 && x.source == y.source,
                "{a}-{b} 换个次序就变了:{} vs {}",
                x.value,
                y.value
            );
        }
    }
}
