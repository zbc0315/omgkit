//! 2D 坐标生成与分子结构绘图。
//!
//! ```
//! use omgkit_depict::{generate, style::Style};
//!
//! let mut m = omgkit_io::smiles::parse("OC(=O)c1ccccc1").unwrap();
//! omgkit_chem::pipeline::sanitize(&mut m).unwrap();
//!
//! let d = generate(&m, &Style::ACS_1996);
//! assert_eq!(d.coords.len(), m.num_atoms());
//! ```
//!
//! # 布局与绘制**不是**两件独立的事
//!
//! 几何骨架确实与规范无关:苯环是正六边形,跟用什么字号没关系。但**"算不算挤
//! 在一起"是规范相关的** —— 原子标签要占地方,而标签尺寸与键长的比例随规范变:
//!
//! | | ACS 1996 | ChemDraw 默认 |
//! |---|---|---|
//! | 键长 | 14.4 pt | 30 pt |
//! | 原子标签 | 10 pt | 10 pt |
//! | **标签占一个键长的** | **69%** | **33%** |
//!
//! 所以 [`Style`] 同时喂给布局和绘制,而 [`Depiction`] 记下产生它的规范指纹 ——
//! 拿另一套规范去渲染时可以被查出来,不会静默地挤成一团。
//!
//! # 一张图由它自己决定,不由它被写成什么样决定
//!
//! 同一个分子的任何 SMILES 写法给出**全等**的图。布局的每一处平局都按
//! [`canonical_ranks`](omgkit_io::canon::canonical_ranks) 打破,不看原子的存储
//! 下标。这一条有判据守着。
//!
//! # 画不好的地方会说出来
//!
//! 桥环、笼状体系在平面上没有好解,拥挤到一定程度的取代基也排不开。这些如实记在
//! [`Depiction::degraded`] 与 [`Depiction::unresolved`] 里,不假装成功。

#![allow(missing_docs)]

pub mod chains;
pub mod geom;
pub mod label;
pub mod layout;
pub mod orient;
// 位图输出(PNG / JPEG)。模块自己的 `//!` 已经写清楚了 —— 这里再挂一层 `///`
// 的话,两段文档会合并,而合并后整段的链接是按**外层**(crate 根)的作用域解析的,
// 于是 `[`to_png`]` 这类同模块内的链接全部解析不了,`cargo doc` 直接报错。
#[cfg(feature = "raster")]
pub mod raster;
pub mod refine;
pub mod render;
pub mod rings;
pub mod stereo;
pub mod style;
pub mod svg;

use std::collections::BTreeMap;

use omgkit_core::MolBuilder;

use geom::Point2;
use rings::Degradation;
use style::Style;

/// 一张 2D 图。
#[derive(Debug, Clone, PartialEq)]
pub struct Depiction {
    /// 逐原子坐标,下标与 [`MolBuilder`] 的原子下标一致。
    ///
    /// 单位是**键长**,不是埃 —— 2D 结构图不是比例模型。换算成 pt/px 由
    /// [`Style::bond_length_pt`] 负责。
    pub coords: Vec<Point2>,
    /// 布局中不得不退化的地方(桥环等)。
    pub degraded: Vec<Degradation>,
    /// 消冲突之后**仍然挤着**的原子对。
    pub unresolved: Vec<(u32, u32)>,
    /// 仍然交叉的键对。
    pub crossings: Vec<(u32, u32)>,
    /// 逐键的楔形指派,下标与 [`MolBuilder`] 的键下标一致。
    pub wedges: Vec<stereo::Wedge>,
    /// **没能画出构型的立体中心**。如实报出来,不假装画好了。
    pub unwedged: Vec<u32>,
    /// 产生这张图的规范名。
    pub style_name: &'static str,
    /// 规范中**影响布局**那部分的指纹。
    ///
    /// 拿另一套规范渲染时,用它可以查出错配 —— 见 [`Style::layout_fingerprint`]。
    pub style_fingerprint: u64,
}

impl Depiction {
    /// 这张图是不是按 `style` 排的。
    ///
    /// 只比布局相关的那部分:换个线宽、换个字体不会让已有坐标失效。
    #[must_use]
    pub fn matches(&self, style: &Style) -> bool {
        self.style_fingerprint == style.layout_fingerprint()
    }

    /// 有没有任何一处没画好(退化、仍在碰撞、仍有交叉)。
    #[must_use]
    pub fn is_clean(&self) -> bool {
        self.degraded.is_empty()
            && self.unresolved.is_empty()
            && self.crossings.is_empty()
            && self.unwedged.is_empty()
    }
}

/// 给分子生成 2D 坐标。
///
/// 分子**应当先净化** —— 环感知的结果决定环系统怎么划分。没净化过也能跑,
/// 但环会被当成链画出来。
#[must_use]
pub fn generate(mol: &MolBuilder, style: &Style) -> Depiction {
    let ranks = omgkit_io::canon::canonical_ranks(mol);

    // **配位键在几何上就是一根线。** 环感知按化学口径把配位键排除在环外
    // (`omgkit_chem::sssr` 的约定),而布局是几何,照那个口径走的话
    // `N->1CCCCC1` 就被当成一条链,闭环的那根键被拉到 4 个键长长。
    //
    // 所以布局用一份把配位键当单键的副本。拓扑、原子编号都不变,坐标直接对得上;
    // 画的时候仍按原分子的键级走。
    let laid = as_plain_bonds(mol);
    let mol = laid.as_ref().unwrap_or(mol);

    let mut pieces = layout::layout_all(mol, &ranks, style);
    // **分量从左到右的次序也要与写法无关。** 分量本身是按连通性收集的,次序跟着
    // 原子的存储下标走 —— 同一个盐换个写法,两个离子就左右对调,于是整张图的
    // 每一个图元都挪了位。实测:语料里 4 个盐/配合物正是这么差出来的。
    pieces.sort_by_key(|p| {
        p.pos
            .keys()
            .map(|a| ranks[*a as usize])
            .min()
            .unwrap_or(u32::MAX)
    });
    let degraded: Vec<Degradation> = pieces.iter().flat_map(|p| p.degraded.clone()).collect();

    // 分量并排摆开,再一起消冲突 —— 分量之间也可能撞上。
    //
    // **宽度要按包围盒算,不能只按原子中心。** 原子标签向两侧伸出去,只按中心
    // 排会让相邻分量的标签叠在一起,而消冲突动不了分量(它只翻可旋转键,单原子
    // 分量连键都没有)。实测:NaCl 的两个离子中心间距正好 1.0,而两个标签半径
    // 各约 0.5 —— 正好贴上。
    let radii = refine::radii(mol, style);
    let mut pos: BTreeMap<u32, Point2> = BTreeMap::new();
    let mut shift = 0.0f64;
    for p in &pieces {
        let (lo, hi) = extent(p.pos.iter().map(|(a, q)| (*q, radii[*a as usize])));
        for (a, q) in &p.pos {
            pos.insert(*a, Point2::new(q.x - lo + shift, q.y));
        }
        shift += hi - lo + PIECE_GAP;
    }

    // **顺反先摆对,再消冲突。** 反过来的话,消冲突翻的那几根键会把顺反弄反 ——
    // 见 `stereo::fix_cis_trans` 的文档。
    let mut flat = vec![Point2::ORIGIN; mol.num_atoms()];
    for (a, q) in &pos {
        flat[*a as usize] = *q;
    }
    stereo::fix_cis_trans(mol, &mut flat);
    for (a, q) in pos.iter_mut() {
        *q = flat[*a as usize];
    }

    let report = refine::relieve(mol, &mut pos, &ranks, style);

    let mut coords = vec![Point2::ORIGIN; mol.num_atoms()];
    for (a, q) in pos {
        coords[a as usize] = q;
    }

    // **摆正要排在楔形指派之前。** 规范朝向里可能含一次镜像,而镜像会把手性
    // 画反;楔形是照最终坐标算的,先摆正再指派,构型自然是对的。反过来做会把
    // 已经画好的楔形悬空 —— 而且线条本身看不出毛病。
    orient::canonicalise(&mut coords, &ranks);

    let w = stereo::assign_wedges(mol, &coords, &ranks);

    Depiction {
        coords,
        wedges: w.bonds,
        unwedged: w.unwedged,
        degraded,
        unresolved: report.unresolved,
        crossings: report.crossings,
        style_name: style.name,
        style_fingerprint: style.layout_fingerprint(),
    }
}

/// 把配位键换成单键的副本;没有配位键就返回 `None`(不必复制)。
///
/// 只动键级,拓扑与原子编号一概不变 —— 坐标因此可以直接用在原分子上。
fn as_plain_bonds(mol: &MolBuilder) -> Option<MolBuilder> {
    if !mol
        .bonds()
        .iter()
        .any(|b| b.order == omgkit_core::BondOrder::Dative)
    {
        return None;
    }
    let mut copy = MolBuilder::with_capacity(mol.num_atoms(), mol.num_bonds());
    for a in mol.atoms() {
        copy.add_atom_data(*a);
    }
    for b in mol.bonds() {
        let mut bd = *b;
        if bd.order == omgkit_core::BondOrder::Dative {
            bd.order = omgkit_core::BondOrder::Single;
        }
        copy.add_bond_data(bd).ok()?;
    }
    // 换了键级就要重新做环感知,否则新成的环没人知道
    omgkit_chem::pipeline::sanitize(&mut copy).ok()?;
    Some(copy)
}

/// 相邻两个分量之间留的空当,单位是键长。半个键长足以看出是两块东西。
const PIECE_GAP: f64 = 0.5;

/// 一组带半径的点在 x 方向上的范围(含半径)。空集给 (0, 0)。
fn extent(pts: impl Iterator<Item = (Point2, f64)>) -> (f64, f64) {
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for (p, r) in pts {
        lo = lo.min(p.x - r);
        hi = hi.max(p.x + r);
    }
    if lo.is_finite() {
        (lo, hi)
    } else {
        (0.0, 0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn prep(smi: &str) -> MolBuilder {
        let mut m = omgkit_io::smiles::parse(smi).unwrap();
        omgkit_chem::pipeline::sanitize(&mut m).unwrap();
        m
    }

    /// 形状指纹:两两距离排序后的多重集。与原子编号、平移、旋转、镜像都无关。
    fn shape_key(smi: &str, style: &Style) -> Vec<i64> {
        let d = generate(&prep(smi), style);
        let mut ds: Vec<i64> = (0..d.coords.len())
            .flat_map(|i| ((i + 1)..d.coords.len()).map(move |j| (i, j)))
            .map(|(i, j)| (d.coords[i].dist(d.coords[j]) * 1e4).round() as i64)
            .collect();
        ds.sort_unstable();
        ds
    }

    #[test]
    fn the_same_molecule_written_differently_gets_the_same_picture() {
        // 整个库的核心不变量。**它属于最终输出,不属于中间阶段** —— 布局本身
        // 会因 SSSR 给出的环原子顺序而有差异,消冲突之后才收敛。所以这条判据
        // 放在 `generate` 这一层;放在 `layout_all` 上是测错了对象。
        let groups = [
            vec![
                "CC(=O)Oc1ccccc1C(=O)O",
                "O=C(C)Oc1ccccc1C(O)=O",
                "OC(=O)c1ccccc1OC(C)=O",
            ],
            vec!["CC(C)(C)c1ccccc1", "c1ccccc1C(C)(C)C", "CC(c1ccccc1)(C)C"],
            vec!["C1CC2(CC1)CCCC2", "C1CCC2(C1)CCCC2"],
            vec!["c1ccc2ccccc2c1", "c1ccc2c(c1)cccc2"],
            vec!["CCCCO", "OCCCC"],
        ];
        for style in &Style::ALL {
            for ws in &groups {
                let keys: Vec<Vec<i64>> = ws.iter().map(|s| shape_key(s, style)).collect();
                for (w, k) in ws.iter().zip(&keys).skip(1) {
                    assert_eq!(&keys[0], k, "[{}] {w} 与 {} 形状不同", style.name, ws[0]);
                }
            }
        }
    }

    #[test]
    fn no_two_atoms_land_on_the_same_spot() {
        for smi in [
            "CC(=O)Oc1ccccc1C(=O)O",
            "OC(=O)c1ccccc1OC(C)=O",
            "c1ccc2ccccc2c1",
            "CC(C)(C)c1ccccc1",
            "CCCCCCCC",
            "C1CC2(CC1)CCCC2",
            "[Na+].[Cl-]",
            "CN1C=NC2=C1C(=O)N(C)C(=O)N2C",
        ] {
            let d = generate(&prep(smi), &Style::ACS_1996);
            for i in 0..d.coords.len() {
                for j in (i + 1)..d.coords.len() {
                    let dist = d.coords[i].dist(d.coords[j]);
                    assert!(dist > 0.3, "{smi}:原子 {i} 与 {j} 距离只有 {dist:.4}");
                }
            }
        }
    }

    #[test]
    fn every_bond_keeps_its_unit_length() {
        for smi in [
            "CC(=O)Oc1ccccc1C(=O)O",
            "c1ccc2ccccc2c1",
            "CCCCCCCC",
            "CN1C=NC2=C1C(=O)N(C)C(=O)N2C",
        ] {
            let m = prep(smi);
            let d = generate(&m, &Style::ACS_1996);
            for b in m.bonds() {
                let len = d.coords[b.begin as usize].dist(d.coords[b.end as usize]);
                assert!(
                    (len - 1.0).abs() < 1e-9,
                    "{smi} 键 {}–{} 长 {len}",
                    b.begin,
                    b.end
                );
            }
        }
    }

    #[test]
    fn disconnected_components_do_not_sit_on_top_of_each_other() {
        let m = prep("[Na+].[Cl-]");
        let d = generate(&m, &Style::ACS_1996);
        assert!(d.coords[0].dist(d.coords[1]) > 1.0, "两个离子挨得太近");
    }

    #[test]
    fn a_depiction_knows_which_style_made_it() {
        // 按 A 规范排版、拿 B 规范渲染会挤成一团,而且不报错。指纹让它可查。
        let d = generate(&prep("CCO"), &Style::ACS_1996);
        assert!(d.matches(&Style::ACS_1996));
        assert!(!d.matches(&Style::CHEMDRAW_DEFAULT), "换了规范却认为匹配");
        assert_eq!(d.style_name, "ACS Document 1996");

        // 只改渲染项不该让坐标失效
        let mut only_render = Style::ACS_1996;
        only_render.line_width_pt = 3.0;
        assert!(d.matches(&only_render), "改线宽不该让已有坐标失效");
    }

    #[test]
    fn trouble_is_reported_rather_than_hidden() {
        // 桥环没有平面好解 —— 必须报出来
        let d = generate(&prep("C1CC2CCC1CC2"), &Style::ACS_1996);
        assert!(!d.degraded.is_empty(), "桥环应当记进 degraded");
        assert!(!d.is_clean());

        // 一般分子应当是干净的
        let ok = generate(&prep("CCO"), &Style::ACS_1996);
        assert!(ok.is_clean(), "乙醇不该有任何问题:{ok:?}");
    }
}
