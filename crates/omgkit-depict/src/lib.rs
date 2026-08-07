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
//! [`ranks_of`] 打破,不看原子的存储下标。这一条有判据守着。
//!
//! # 画不好的地方会说出来
//!
//! 桥环、笼状体系在平面上没有好解,拥挤到一定程度的取代基也排不开。这些如实记在
//! [`Depiction::degraded`] 与 [`Depiction::unresolved`] 里,不假装成功。

#![allow(missing_docs)]

/// 布局用的原子秩:**先按对称等价类,类内再按规范 SMILES 的输出次序**。
///
/// 不是 [`canonical_ranks`](omgkit_io::canon::canonical_ranks) —— 那一个的
/// 深层平局是任取的。
///
/// # `canonical_ranks` 在哪一步失守
///
/// 它自己的模块文档写着「**更深层次的并列仍是任取**」:`break_all_ties` 把每个
/// 还没分开的格里**存储序最靠前**的那个原子劈出去。对大多数分子这不要紧(细化
/// 早就分完了),可一旦剩下真正的对称,秩就跟着写法走。
///
/// 这不是理论上的。**内消旋分子最吃亏**:1-乙炔基-4-苯基环己醇(语料第 573 行)
/// 有一个穿过 C1、C4 的镜面,两条环支路**构造上等价、构型上相反**。1-WL 细化
/// 分不开它们 —— 立体宇称是相对邻居的**等价类**算的,镜面下不变。于是任取一头,
/// 两种写法把楔形画到了相反的方向;两张图都对,但不是同一张。
///
/// **它们并不真的等价** —— 这一点是修法能成立的前提,值得说死。把两种写法各按
/// 自己的 `canonical_ranks` 写成规范风格的串:
///
/// ```text
/// 写法 1: C#C[C@@]1(CC[C@@H](CC1)c1ccccc1)O
/// 写法 2: C#C[C@] 1(CC[C@H] (CC1)c1ccccc1)O
/// ```
///
/// 骨架逐字相同,**每一个立体标记都翻了** —— 两套标号差的是一个**反自同构**
/// (镜面),不是自同构。若真是自同构下等价,两串会完全一样,那么"遍历起点取
/// 字典序最小"也同样分不开,本函数就白写了。正因为不等价,取最小串把这个自由度
/// 消掉了。
///
/// # 为什么是"类在前、序在后",不是直接用输出次序
///
/// 两种秩的**语义不同**,这一点是实测撞出来的:
///
/// - [`symmetry_classes`](omgkit_io::canon::symmetry_classes) 是 1-WL 细化的
///   结果 —— 类编号反映**结构角色**,布局的启发式吃的正是这个。
/// - [`atom_order`](omgkit_io::smiles::Written::atom_order) 是规范 SMILES 的
///   **DFS 输出序** —— 唯一、含立体,但结构上是任意的。
///
/// 直接拿输出次序当秩(试过,全量实测):头号契约降到 2,可**键交叉从 50 涨到
/// 78**,重跑模板生成器也只收回到 70。布局把"秩小"当"重要",而 DFS 序不是那个
/// 意思。
///
/// 分两级就两头都拿到了:主键仍是细化出来的类(结构语义原样保留),**只有类内
/// 那点任取被换成规范次序**。
///
/// # 三个方案的全量对照
///
/// 8831 分子 × 2 规范 × 30 种写法,**三列都用同一张(旧)模板表**,好把"换秩"
/// 这一件事单独看清:
///
/// | | `canonical_ranks` | 纯输出次序 | **类 + 输出次序** |
/// |---|---:|---:|---:|
/// | **写法无关(头号契约)** | 9 | 2 | **3** |
/// | 其中有键交叉 | 50 | 78 | **48** |
/// | 干净 | 16191 | 16169 | **16192** |
/// | 有未解冲突 | 1131 | 1150 | **1130** |
/// | 标签塞不下 | 859 | 846 | **858** |
/// | 键角不过窄 | 180 / 16191 | — | **181 / 16192** |
/// | 硬性质其余八条 | 基准 | 基准 | **一处没动** |
/// | 外部判官 | 496 / 0 | 496 / 0 | **496 / 0** |
///
/// `键角不过窄` 那一格不是新缺陷,是判据的适用范围变大了(一个分子变干净了,
/// 首次进入这条判据)—— 分母也同步 +1,细节见 `harness/README.md`。
///
/// 取第三列:**头号契约降三分之二,而质量指标全部持平或略好**,不是拿别处换的。
///
/// 换秩之后模板表跟着重生成了(否则「重跑逐字节相同」这条验收作废),那一步
/// 另有代价 —— 见 `harness/README.md`。
///
/// # 一定要与 [`hydrogens::with_stereo_hs`] 用同一个
///
/// 补出来的氢按秩排序追加,它若与这里的秩不同源,同一个分子换种写法补出来的
/// 氢就会拿到不同的原子号 —— 后面整条管线跟着变。
#[must_use]
pub fn ranks_of(mol: &omgkit_core::MolBuilder) -> Vec<u32> {
    // 细化到不动点的对称等价类 —— 与写法无关,而且**保住了结构语义**
    let classes = omgkit_io::canon::symmetry_classes(mol);
    // 规范 SMILES 的输出次序 —— 唯一,含立体,只用来打破类内的平局
    let w = omgkit_io::canon::canonical_smiles(mol);
    // **写不全的话,漏掉的原子会静默留在 `pos = 0`**,与类内第一个并列,于是
    // 那一处平局退回存储序 —— 正是这个函数要消掉的东西。实测全语料(含补氢
    // 之后)0 个分子写不全,但这条前提本来没人守。
    debug_assert_eq!(
        w.atom_order.len(),
        mol.num_atoms(),
        "规范 SMILES 没把所有原子写出来,类内平局会退回存储序"
    );
    let mut pos = vec![0u32; mol.num_atoms()];
    for (i, a) in w.atom_order.iter().enumerate() {
        pos[*a as usize] = u32::try_from(i).expect("原子数超出 u32");
    }
    let mut order: Vec<u32> =
        (0..u32::try_from(mol.num_atoms()).expect("原子数超出 u32")).collect();
    order.sort_by_key(|a| (classes[*a as usize], pos[*a as usize]));
    let mut r = vec![0u32; mol.num_atoms()];
    for (i, a) in order.iter().enumerate() {
        r[*a as usize] = u32::try_from(i).expect("原子数超出 u32");
    }
    r
}

pub mod chains;
pub mod geom;
pub mod hydrogens;
pub mod label;
pub mod layout;
pub mod orient;
// 位图输出(PNG / JPEG)。模块自己的 `//!` 已经写清楚了 —— 这里再挂一层 `///`
// 的话,两段文档会合并,而合并后整段的链接是按**外层**(crate 根)的作用域解析的,
// 于是 `[`to_png`]` 这类同模块内的链接全部解析不了,`cargo doc` 直接报错。
/// 桥环的几何摆法。**目前只给离线的模板生成器用**,运行时还走
/// [`rings::relax`] —— 先拿它把表刷好、量清楚,再决定要不要接进在线路径。
#[cfg(test)]
mod arcs;
#[cfg(feature = "raster")]
pub mod raster;
pub mod refine;
pub mod render;
pub mod rings;
pub mod stereo;

pub mod style;
pub mod svg;
/// 桥环骨架的预存坐标表。见模块文档。
pub mod templates;

use std::collections::BTreeMap;

use omgkit_core::MolBuilder;

use geom::Point2;
use rings::Degradation;
use style::Style;

/// 一张 2D 图。
///
/// # 下标是相对**被画的那个分子**的
///
/// 为了画出构型,某些立体中心要补一个显式氢(见 [`hydrogens`])。那时被画的
/// 分子比传进来的多几个原子,而 `coords`、`wedges` 这些逐原子/逐键的向量是按
/// **补完之后**的编号排的 —— 拿 [`Depiction::drawn`] 取回那个分子。
///
/// **前 `mol.num_atoms()` 个原子、前 `mol.num_bonds()` 根键与传入的分子逐项
/// 对应**,所以按原下标索引仍然是对的;多出来的排在后面。
#[derive(Debug, Clone, PartialEq)]
pub struct Depiction {
    /// 逐原子坐标,下标与**被画的那个分子**([`Depiction::drawn`])一致。
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
    /// 为了画出构型补出来的原子/键。空的话画的就是传进来的分子。
    pub added: hydrogens::Augmented,
}

impl Depiction {
    /// 这张图是不是按 `style` 排的。
    ///
    /// 只比布局相关的那部分:换个线宽、换个字体不会让已有坐标失效。
    #[must_use]
    pub fn matches(&self, style: &Style) -> bool {
        self.style_fingerprint == style.layout_fingerprint()
    }

    /// **真正被画的那个分子。** 没补东西时就是传进来的那个,不复制。
    ///
    /// `coords`、`wedges`、`unresolved`、`crossings`、`unwedged`、`degraded`
    /// 的下标全部相对它。渲染与判据都该拿它,而不是拿传进来的分子 —— 否则
    /// 补出来的氢会被静默丢掉,而诊断全绿。
    ///
    /// 返回 [`Cow`](std::borrow::Cow),所以 `&d.drawn(&m)` 在要 `&MolBuilder`
    /// 的地方直接能用(靠 `Deref`)。
    #[must_use]
    pub fn drawn<'a>(&self, mol: &'a MolBuilder) -> std::borrow::Cow<'a, MolBuilder> {
        if self.added.is_empty() {
            std::borrow::Cow::Borrowed(mol)
        } else {
            std::borrow::Cow::Owned(self.added.apply(mol))
        }
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

/// 判据共用的解析 + 净化。
#[cfg(test)]
pub(crate) fn tests_prep(smi: &str) -> MolBuilder {
    let mut m = omgkit_io::smiles::parse(smi).expect("测试用的 SMILES 该能解析");
    omgkit_chem::pipeline::sanitize(&mut m).expect("测试用的分子该能净化");
    m
}

/// 给分子生成 2D 坐标。
///
/// 分子**应当先净化** —— 环感知的结果决定环系统怎么划分。没净化过也能跑,
/// 但环会被当成链画出来。
#[must_use]
pub fn generate(mol: &MolBuilder, style: &Style) -> Depiction {
    generate_with(mol, style, None)
}

/// 同 [`generate`],但可以临时顶替模板表里某一条。**只给离线的模板生成器用。**
///
/// 生成器要问"把这组坐标装进去之后,真实分子画出来好不好",而 `generate` 会查
/// 那张表 —— 表正是它在生成的东西。这个参数把那层循环拆开,见
/// [`templates::Override`]。
pub(crate) fn generate_with(
    mol: &MolBuilder,
    style: &Style,
    over: templates::Override<'_>,
) -> Depiction {
    // **先补显式氢,再做别的。** 有些立体中心三根键全在环上,唯一合法的楔形是
    // C–H —— 那个氢不画出来,构型就只能画到环键上(见 [`hydrogens`])。
    //
    // 补出来的原子接在**末尾**,原有编号一概不变,所以下面整条管线原样跑在补完
    // 的分子上就行:布局、消冲突、规范朝向、楔形指派全都自动把那个氢算进去。
    let added = hydrogens::with_stereo_hs(mol).unwrap_or_default();
    let grown = (!added.is_empty()).then(|| added.apply(mol));
    let mol = grown.as_ref().unwrap_or(mol);

    let ranks = ranks_of(mol);

    // **配位键在几何上就是一根线。** 环感知按化学口径把配位键排除在环外
    // (`omgkit_chem::sssr` 的约定),而布局是几何,照那个口径走的话
    // `N->1CCCCC1` 就被当成一条链,闭环的那根键被拉到 4 个键长长。
    //
    // 所以布局用一份把配位键当单键的副本。拓扑、原子编号都不变,坐标直接对得上;
    // 画的时候仍按原分子的键级走。
    let laid = as_plain_bonds(mol);
    let mol = laid.as_ref().unwrap_or(mol);

    let mut pieces = layout::layout_all(mol, &ranks, style, over);
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
    stereo::fix_cis_trans(mol, &mut flat, &ranks);
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
        added,
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

    /// 完整的图元指纹 —— 连楔形的方向都比。`shape_key` 只比两两距离,
    /// 坐标一模一样而楔形互换的那一类它看不见。
    fn scene_key(smi: &str, style: &Style) -> Vec<String> {
        let m = prep(smi);
        let d = generate(&m, style);
        let q = |p: Point2| format!("{:.3},{:.3}", p.x, p.y);
        let mut v: Vec<String> = render::scene(&m, &d, style)
            .items
            .iter()
            .map(|it| match it {
                render::Primitive::Line { from, to, .. } => {
                    let (x, y) = (q(*from), q(*to));
                    // 线不分方向 —— 谁是起点取决于键的 begin/end
                    if x <= y {
                        format!("L {x} {y}")
                    } else {
                        format!("L {y} {x}")
                    }
                }
                // 楔形分方向:窄端宽端不是一回事
                render::Primitive::Wedge { from, to, .. } => format!("W {} {}", q(*from), q(*to)),
                render::Primitive::Hash { from, to, .. } => format!("H {} {}", q(*from), q(*to)),
                // **文本连内容一起比。** 只比落点的话,`OH` 变成 `HO`、
                // 竖排翻成横排这类退化看不见 —— 而那正是坐标全同、图元不同的
                // 另一大类。
                render::Primitive::Text { at, runs, .. } => format!("T {} {runs:?}", q(*at)),
            })
            .collect();
        v.sort();
        v
    }

    #[test]
    fn a_mirror_symmetric_molecule_gets_the_same_wedges_whichever_way_it_is_written() {
        // **头号契约里最隐蔽的一档。** 这些分子换种写法之后**坐标逐字节相同**,
        // 差的只是楔形与虚楔互换 —— `shape_key` 那条判据完全看不见,得比整份
        // 图元。
        //
        // 根子在 `canonical_ranks` 的深层平局是任取的,而**内消旋分子最吃亏**:
        // 1-乙炔基-4-苯基环己醇有一个穿过 C1、C4 的镜面,两条环支路构造上等价、
        // 构型上相反,1-WL 细化分不开,于是任取一头 —— 两张图都对,但不是同
        // 一张。修法见 [`ranks_of`]。
        //
        // 变异验证:把 `ranks_of` 换回 `canonical_ranks`,这条当场红。
        let groups = [
            // 573:非手性,楔形/虚楔互换
            vec![
                "C(#C)[C@@]1(CC[C@H](C2=CC=CC=C2)CC1)O",
                "C1C[C@H](CC[C@@]1(O)C#C)c1ccccc1",
            ],
            // 2553:笼状胺,坐标多重集相同而线连在不同的点对之间
            vec!["C1CN2CN1CN3CCN(C2)C3", "C1N2CN(CN3CCN(C2)C3)C1"],
        ];
        // **先证明这几对写法真的换了存储序。** 不然改写退化成恒等,这条判据就
        // 静悄悄地空过了 —— `audit.rs` 的搅拌器为这个失效模式专门立过案(旧那个
        // 乘法哈希有 10.85% 的改写原样返回)。
        for ws in &groups {
            let seqs: Vec<Vec<(u8, u8)>> = ws
                .iter()
                .map(|s| {
                    let m = prep(s);
                    (0..m.num_atoms())
                        .map(|i| {
                            let a = m.atoms()[i];
                            (a.atomic_num, a.num_explicit_hs + a.num_implicit_hs)
                        })
                        .collect()
                })
                .collect();
            assert!(
                seqs[1..].iter().any(|s| *s != seqs[0]),
                "{} 与 {} 的存储序一模一样,这一组验不了写法无关",
                ws[0],
                ws[1]
            );
        }
        for style in &Style::ALL {
            for ws in &groups {
                let keys: Vec<Vec<String>> = ws.iter().map(|s| scene_key(s, style)).collect();
                for (w, k) in ws.iter().zip(&keys).skip(1) {
                    assert_eq!(
                        &keys[0], k,
                        "[{}] {w} 与 {} 画出来的图元不同",
                        style.name, ws[0]
                    );
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
