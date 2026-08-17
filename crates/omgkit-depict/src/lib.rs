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
/// 这不是理论上的。**内消旋分子最吃亏**:1-乙炔基-4-苯基环己醇(语料第 574 行)
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
/// # 实现搬到 `omgkit-io` 去了,论证留在这儿
///
/// 三维构象生成([`omgkit-conformer`](https://docs.rs/omgkit-conformer))也要
/// 这个秩,而它**不该依赖绘图 crate**。所以函数体挪到了
/// [`omgkit_io::canon::classed_ranks`],这里只剩一层转发。
///
/// 上面那些数(键交叉、标签塞不下、外部判官)全是**绘图指标**,放在解析/
/// 规范化那一层讲不通,所以留在这里 —— 权威说明是本函数,`classed_ranks`
/// 那边只写契约。
#[must_use]
pub fn ranks_of(mol: &omgkit_core::MolBuilder) -> Vec<u32> {
    omgkit_io::canon::classed_ranks(mol)
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
/// 桥环的几何摆法。**离线的模板生成器与运行时都用它**;运行时排在查表之后、
/// [`rings::relax`] 之前 —— 见模块文档「什么时候轮到它」。
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
///
/// # 顺反还要单独感知一次
///
/// 净化那 12 步里**没有**双键顺反感知(它要用对称等价类,那在净化的上一层,
/// 调不到)。只跑了净化的分子,每根双键的
/// [`stereo`](omgkit_core::BondData::stereo) 都是 `None`,于是顺反校正
/// (`stereo::fix_cis_trans`)整个空转 —— **E/Z 可能画反,而线条本身看着一点
/// 毛病没有**。这是"画错了",不是"没画好",所以 debug 构建下当场拦住:
///
/// ```no_run
/// # use omgkit_core::MolBuilder;
/// # fn demo(mut mol: MolBuilder) {
/// omgkit_chem::pipeline::sanitize(&mut mol).unwrap();
/// omgkit_io::stereo::perceive_bond_stereo(&mut mol); // ← 别漏
/// let d = omgkit_depict::generate(&mol, &omgkit_depict::style::Style::ACS_1996);
/// # let _ = d;
/// # }
/// ```
///
/// 判法见
/// [`directions_not_perceived`](omgkit_io::stereo::directions_not_perceived) ——
/// 只有"双键两端的方向键成对写着、而它自己没有顺反"才报,没写方向的分子一律
/// 放行。release 下不做这个检查。
///
/// (Python 绑定的 `Mol.sanitize()` 把两步合在一起,从 Python 看不出这个区别。)
#[must_use]
pub fn generate(mol: &MolBuilder, style: &Style) -> Depiction {
    debug_assert!(
        !omgkit_io::stereo::directions_not_perceived(mol),
        "这个分子的双键几何**方向键已经写明**、却没有感知过顺反 —— \
         漏了 omgkit_io::stereo::perceive_bond_stereo。这样画不会报错,\
         但顺反校正整个空转,E/Z 可能画反"
    );
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

    // **η5 配位的那 5 根 σ 键会把环感知搅成一团。** 布局只留一根代表键,
    // 见 [`hapto_extras`]。渲染仍拿原分子,一根键不少。
    let hapto = hapto_extras(mol, &ranks);
    let thinned = hapto.as_ref().and_then(|(extras, _)| {
        let mut copy = MolBuilder::with_capacity(mol.num_atoms(), mol.num_bonds());
        for a in mol.atoms() {
            copy.add_atom_data(*a);
        }
        for (bi, b) in mol.bonds().iter().enumerate() {
            if !extras.contains(&bi) {
                copy.add_bond_data(*b).ok()?;
            }
        }
        omgkit_chem::pipeline::sanitize(&mut copy).ok()?;
        Some(copy)
    });
    // 摘细之前的那份 —— 交叉要拿它再算一遍,见下面
    let whole = mol;
    let mol = thinned.as_ref().unwrap_or(mol);

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
    let mut degraded: Vec<Degradation> = pieces.iter().flat_map(|p| p.degraded.clone()).collect();
    // η 配位那几根键在平面上不可能等长,如实记一笔 ——
    // 见 [`Degradation::HaptoCoordination`](rings::Degradation::HaptoCoordination)。
    if let Some((_, told)) = &hapto {
        degraded.extend(told.iter().cloned());
    }

    // 分量并排摆开,再一起消冲突 —— 分量之间也可能撞上。
    //
    // **宽度要按包围盒算,不能只按原子中心。** 原子标签向两侧伸出去,只按中心
    // 排会让相邻分量的标签叠在一起,而消冲突动不了分量(它只翻可旋转键,单原子
    // 分量连键都没有)。实测:NaCl 的两个离子中心间距正好 1.0,而两个标签半径
    // 各约 0.5 —— 正好贴上。
    // **逐原子/逐键的东西一律拿 `whole`。** 摘细是删键,键下标会整体前移
    // —— 拿摘细副本算出来的逐键向量与被画的分子对不上号。
    let radii = refine::radii(whole, style);
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
    stereo::fix_cis_trans(whole, &mut flat, &ranks);
    for (a, q) in pos.iter_mut() {
        *q = flat[*a as usize];
    }

    let mut report = refine::relieve(mol, &mut pos, &ranks, style);

    // **摘细过的话,交叉要拿原分子再算一遍。** 消冲突跑在摘细过的副本上,
    // 看不见被摘掉的那些 η5 键 —— 而它们照样会画出来。二茂铁实测:不补这一步
    // 报的是 0 处交叉,而图上明明有。**画不好就要说出来。**
    //
    // 只在真摘过的时候算(全量语料 2 个分子),所以不给别人添成本。
    if thinned.is_some() {
        report.crossings = refine::crossings(whole, &pos);
    }

    let mut coords = vec![Point2::ORIGIN; mol.num_atoms()];
    for (a, q) in pos {
        coords[a as usize] = q;
    }

    // **摆正要排在楔形指派之前。** 规范朝向里可能含一次镜像,而镜像会把手性
    // 画反;楔形是照最终坐标算的,先摆正再指派,构型自然是对的。反过来做会把
    // 已经画好的楔形悬空 —— 而且线条本身看不出毛病。
    orient::canonicalise(&mut coords, &ranks);

    // 同上:`Depiction::wedges` 的下标必须与**被画的那个分子**一致。
    // 拿摘细副本的话 `wedges.len()` 会比键数少,`dump_molblock` 那种按下标
    // 取的调用方直接越界 —— 实测二茂铁 wedges.len()=12 而键数 20。
    let w = stereo::assign_wedges(whole, &coords, &ranks);

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

/// 环上被 `picked` 选中的那些原子,是不是**首尾相接的一段**(整圈也算)。
///
/// `ring` 是 [`omgkit_chem::sssr::Ring::atoms`],已经按环序排好。
///
/// # 这一条把大环螯合物挡在外面
///
/// 只数"金属打进这个环几根键"是不够的:卟啉、酞菁、环多胺那一类,摘掉金属之后
/// 照样是个环,而金属对它有 4 根键。它们与 η 配位的差别在**位置**——
/// η5 的 Cp 是 5 个**首尾相接**的碳,而 cyclam 的 4 个 N 之间隔着 2~3 个碳。
fn contiguous_on_ring(ring: &[u32], picked: &std::collections::BTreeSet<u32>) -> bool {
    let n = ring.len();
    let k = ring.iter().filter(|a| picked.contains(a)).count();
    if k < 2 {
        return k == 1;
    }
    if k == n {
        return true; // 整圈都选上了
    }
    // 选中的原子之间"断开"了几次;连续的一段只断一次
    let breaks = (0..n)
        .filter(|i| picked.contains(&ring[*i]) && !picked.contains(&ring[(i + 1) % n]))
        .count();
    breaks == 1
}

/// η<sup>n</sup> 配位里**多余的那些键**,给布局用。
///
/// # 二茂铁的 Fe 度数是 10
///
/// SMILES 把 η5 配位写成 5 根独立的 σ 键(`[Fe]23456789` 这样),于是环感知
/// 吐出 9 到 10 个**三元环**(Fe + 环上相邻两个碳),它们全是这个建模方式的
/// 假象。整个体系被当成一个巨大的桥环系统,画出来两个 Cp 环叠在一起 ——
/// 实测 ChemDraw 规范下 8 处键交叉。
///
/// # 认法:摘掉金属的键再看配体自己的环
///
/// 金属 M 的键全摘掉之后做环感知,得到的才是**配体自己的环**。M 与某个这样的
/// 环之间有 ≥3 根键,就是 η<sup>n</sup> 配位。
///
/// 全量语料实测:8831 个分子里**只有 2 个**(都是二茂铁,都是 η5 × 2)。
/// 范围小,所以这里做的也小。
///
/// # 只留一根代表键
///
/// 布局要的是"Cp 是个普通五元环、金属挂在它旁边",所以每个 (金属, 环) 只留
/// **一根**键、其余摘掉。留哪一根按**规范秩**定,不看存储下标 —— 头号契约。
///
/// 摘掉之后二茂铁的 Fe 只剩两根键(每个 Cp 一根),自然成了两个五边形中间的
/// 那个连接原子,也就是夹心式的画法。
///
/// **只动布局。** 渲染仍拿原分子,10 根键一根不少地画出来。
fn hapto_extras(
    mol: &MolBuilder,
    ranks: &[u32],
) -> Option<(std::collections::BTreeSet<usize>, Vec<Degradation>)> {
    /// 会做 π 配位的元素。从宽收 —— 这里只是找候选,`>=3 根键进同一个环`
    /// 那一条才是判据。
    fn is_metal(z: u8) -> bool {
        matches!(z, 3 | 4 | 11..=13 | 19..=32 | 37..=51 | 55..=84 | 87..=118)
    }
    let metals: Vec<u32> = (0..u32::try_from(mol.num_atoms()).ok()?)
        .filter(|a| is_metal(mol.atoms()[*a as usize].atomic_num) && mol.degree(*a) >= 3)
        .collect();
    if metals.is_empty() {
        return None;
    }

    // 配体自己的环:把金属的键全摘掉再感知
    let mut lig = MolBuilder::with_capacity(mol.num_atoms(), mol.num_bonds());
    for a in mol.atoms() {
        lig.add_atom_data(*a);
    }
    for b in mol.bonds() {
        if !metals.contains(&b.begin) && !metals.contains(&b.end) {
            lig.add_bond_data(*b).ok()?;
        }
    }
    omgkit_chem::pipeline::sanitize(&mut lig).ok()?;
    let rings = omgkit_chem::sssr::ring_set(&lig);

    let mut extras = std::collections::BTreeSet::new();
    let mut told = Vec::new();
    for m in &metals {
        for r in &rings {
            // 这个金属打进这个环的那些键
            let mut into: Vec<(u32, usize)> = mol
                .neighbors(*m)
                .filter(|(nb, _)| r.atoms.contains(nb))
                .map(|(nb, bi)| (ranks[nb as usize], bi as usize))
                .collect();
            if into.len() < 3 {
                continue;
            }
            // **必须是环上连续的一段。** 只数键数会把**大环螯合物**一起收进来:
            // 卟啉、酞菁、环多胺那一类,摘掉金属之后照样是个环,而金属对它有
            // 4 根键 —— 全部满足"≥3 根"。
            //
            // 实测 Ni(环四胺) `[Ni]123N4CCN1CCCN2CCN3CCC4`:只数键数的话
            // 键交叉 **0 → 5**,Ni 到四个 N 的距离从 0.95…1.24 变成
            // 1.00 / 3.51 / 4.34 / 5.49 —— **一张本来画对的图被改坏了**,而且
            // 报的退化写着 η 配位,化学上根本不成立(cyclam 是 σ 给体)。
            //
            // 语料里一个大环配合物都没有,所以全量指标看不见这一条 ——
            // 这是**审核拿真实化合物试出来的**,不是量出来的。
            //
            // 真正的 η 配位是金属压在**一整段连续的 π 体系**上(Cp 的 5 个碳
            // 首尾相接),而 cyclam 的 4 个 N 之间隔着 2~3 个碳。
            let bonded: std::collections::BTreeSet<u32> =
                mol.neighbors(*m).map(|(nb, _)| nb).collect();
            if !contiguous_on_ring(&r.atoms, &bonded) {
                continue;
            }
            // 规范秩最小的那个邻居留下,其余摘掉
            into.sort_unstable();
            extras.extend(into.into_iter().skip(1).map(|(_, bi)| bi));
            // **如实报退化。** 那 n 根键在平面上不可能等长,见
            // [`Degradation::HaptoCoordination`]。次序按规范秩,写法无关。
            let mut ring: Vec<u32> = r.atoms.clone();
            ring.sort_by_key(|a| (ranks[*a as usize], *a));
            told.push(Degradation::HaptoCoordination { metal: *m, ring });
        }
    }
    (!extras.is_empty()).then_some((extras, told))
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

    /// 解析 + 净化 + **顺反感知**。第三步不在净化的 12 步里,漏了的话每根双键的
    /// `stereo` 都是 `None`,顺反校正整个空转 —— 顺式反式画成同一张图而判据照样
    /// 绿。[`generate`] 的 `debug_assert!` 现在会当场拦住这种输入。
    fn prep(smi: &str) -> MolBuilder {
        let mut m = omgkit_io::smiles::parse(smi).unwrap();
        omgkit_chem::pipeline::sanitize(&mut m).unwrap();
        omgkit_io::stereo::perceive_bond_stereo(&mut m);
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
            // 6457:**形状全等而姿态差 15.03°**。布局阶段两种写法完全一致
            // (都是 −14.6099°),分岔出在消冲突 —— `far_side` 在"两侧一样多"
            // 时靠 `begin`/`end` 挑边,而那是书写痕迹。绕同一根轴翻这侧还是
            // 那侧,结果差一次整体反射(`orient` 归得掉);可**两次绕不同轴的
            // 反射合成就是任意角度的旋转**,30° 网格归不掉。
            //
            // 只在 ACS 下犯:ChemDraw 的标签小、消冲突根本没动手。
            vec![
                "N1(C(=C(C(=O)OCC)N=N1)CSC2=NC3N(C4C=CC=CC(C(N=N2)=3)=4)C)C5C(=NON=5)N",
                "n1c2c3ccccc3n(c2nc(n1)SCc1n(-c2nonc2N)nnc1C(OCC)=O)C",
            ],
            // 7879:Ni 的四齿配合物,五个环共用同一个金属。**分岔在布局**
            // (前面几个都在消冲突或摆位):拼环时"共用的那两个原子"取自
            // SSSR 的输出序,而 `fuse_on_bond` 选环心是"取远离已放置质心的
            // 那个" —— 质心落在这根键的中垂线上时两个候选等距,`u`/`v` 一
            // 交换就拼到了相反的一侧。这个分子对称度高,四次拼环的 `(u,v)`
            // **每一次都正好反过来**,40/40 个图元全不同。
            vec![
                "O=C1C[N+]23CC[N+]45CC(=O)O[Ni]24(O1)(OC(=O)C3)OC(=O)C5",
                "O1C(=O)C[N+]23CC(=O)O[Ni]1143OC(C[N+]1(CC2)CC(O4)=O)=O",
            ],
            // 1068:唯一一批走**撑开**(`refine` 的第四个算子)的分子。撑开会
            // 改一个键角,是流水线里最晚、也最容易把姿态带偏的一步,所以这一组
            // 守的是"走完这一步之后整张图仍然与写法无关"。
            //
            // **它守不住的那一件事要说清楚:转向怎么命名。** 实测:把撑开的转向
            // 从内蕴的「朝 n / 背 n」改回固定的 `+30° 在先`,这一组照样绿 ——
            // 因为消冲突之前的坐标系今天本来就逐写法一致(起手环落在起始角为
            // **常数**的正多边形上)。那一条由
            // `refine::tests::the_two_splay_directions_are_named_without_looking_at_the_canvas`
            // 钉着:它把坐标整体反射一次直接验等变性,不经过 `orient::canonicalise`。
            vec![
                "C1(/C(NC2=C(N1)C=CC=C2)=N\\C3=CC=C(C(=O)OCC)C=C3)=N/C4=CC=C(C(=O)OCC)C=C4",
                "c1ccc2c(c1)[nH]c(=N/c1ccc(C(OCC)=O)cc1)/c([nH]2)=N\\c1ccc(cc1)C(=O)OCC",
                "c1c(\\N=c2/c([nH]c3ccccc3[nH]2)=N\\c2ccc(cc2)C(=O)OCC)ccc(c1)C(=O)OCC",
            ],
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
    fn a_sandwich_complex_gets_two_proper_rings_and_says_it_is_degraded() {
        // **η5 配位在 SMILES 里是 5 根独立的 σ 键**,于是环感知吐出一堆
        // 三元环(Fe + 环上相邻两个碳),整个体系被当成一个巨大的桥环系统。
        // 实测二茂铁因此画成两个叠在一起的环、8 处交叉,而且骨架指纹算不出来。
        //
        // 布局每个 (金属, 环) 只留一根代表键之后,Cp 成了普通五元环、金属成了
        // 两个五边形之间的连接原子。这条判据钉住三件事:
        //
        // 1. **两个环真的是正五边形** —— 环内每根键都是单位长;
        // 2. **金属在两个环中间** —— 它到两个环心的方向大致相反;
        // 3. **如实报退化** —— 那 10 根 Fe–C 键在平面上不可能等长。
        //
        // 第 3 条最容易漏:少了它,二茂铁会被算成"干净",于是硬性质
        // `键长全等` 开始查它并**当场破 4 处**(实测)。
        for smi in [
            "C12C3=C4C5=C1[Fe]23456789C%10C6=C7C8=C9%10",
            "CN(C)C[C-]12C3=C4C5=C1[Fe++]23456789[C-]%10C6=C7C8=C9%10",
        ] {
            let m = prep(smi);
            for style in &Style::ALL {
                let d = generate(&m, style);
                let hapto: Vec<&rings::Degradation> = d
                    .degraded
                    .iter()
                    .filter(|x| matches!(x, rings::Degradation::HaptoCoordination { .. }))
                    .collect();
                assert_eq!(
                    hapto.len(),
                    2,
                    "[{}] {smi}:该报两处 η 配位退化,实得 {}",
                    style.name,
                    hapto.len()
                );
                // **交叉也要如实报。** 消冲突跑在摘细过的副本上,看不见被摘掉
                // 的那些 η5 键 —— 而它们照样会画出来(金属摆在环外,连到远端
                // 那几个环原子的线必然穿过环)。不补那一步的话报的是 0 处交叉,
                // 而图上明明有 4 处。变异:去掉 `refine::crossings(whole, …)`
                // 那一句 → 这条红。
                assert!(
                    !d.crossings.is_empty(),
                    "[{}] {smi}:扇出去的 η5 键必然穿过环,交叉不能报 0",
                    style.name
                );
                let mut mid = Vec::new();
                for h in &hapto {
                    let rings::Degradation::HaptoCoordination { metal, ring } = h else {
                        unreachable!("上面已经筛过")
                    };
                    assert_eq!(ring.len(), 5, "[{}] {smi}:Cp 该是五元环", style.name);
                    // 一、环内每根键都是单位长(是个正五边形)
                    for w in ring
                        .windows(2)
                        .chain(std::iter::once(&[ring[0], ring[ring.len() - 1]][..]))
                    {
                        let (p, q) = (d.coords[w[0] as usize], d.coords[w[1] as usize]);
                        // `ring` 是按规范秩排的,不是绕环的次序,所以只查
                        // **成键的**那几对
                        if m.neighbors(w[0]).any(|(x, _)| x == w[1]) {
                            let len = p.dist(q);
                            assert!(
                                (len - 1.0).abs() < 1e-6,
                                "[{}] {smi}:Cp 环上的键长 {len:.4},该是 1",
                                style.name
                            );
                        }
                    }
                    let c = ring
                        .iter()
                        .fold(Point2::ORIGIN, |s, a| s + d.coords[*a as usize])
                        * (1.0 / ring.len() as f64);
                    mid.push((d.coords[*metal as usize], c));
                }
                // 二、金属在两个环中间 —— 它指向两个环心的方向大致相反
                let (fe, c0) = mid[0];
                let c1 = mid[1].1;
                let (u, v) = ((c0 - fe).normalized(), (c1 - fe).normalized());
                assert!(
                    u.dot(v) < -0.3,
                    "[{}] {smi}:两个环该分列金属两侧,实得夹角余弦 {:.3}",
                    style.name,
                    u.dot(v)
                );
            }
        }
    }

    #[test]
    fn a_macrocyclic_chelate_is_not_mistaken_for_a_sandwich() {
        // **只数"金属打进这个环几根键"会把大环螯合物一起收进来。** 卟啉、酞菁、
        // 环多胺那一类,摘掉金属之后照样是个环,而金属对它有 4 根键 —— 全部
        // 满足 `>= 3`。
        //
        // Ni(环四胺)是这一类里最干净的例子(式 C₁₀H₂₀N₄Ni,Ni 四配位全接 N,
        // 配体自己是一个含 4 个 N 的 14 元环)。**只数键数的话实测**:
        // 键交叉 0 → 5,Ni 到四个 N 的距离从 0.95…1.24 变成 1.00/3.51/4.34/5.49
        // —— 一张本来画对的图被改坏了。
        //
        // 区别在**位置**:η5 的 Cp 是 5 个首尾相接的碳,而 cyclam 的 4 个 N
        // 之间隔着 2~3 个碳。判据因此要求"环上连续的一段"。
        //
        // **这个反例是代码审核拿真实化合物试出来的,不是从语料量出来的** ——
        // 语料里一个大环配合物都没有(实测 147 个含金属且度 ≥3 的分子里,
        // "金属打进配体环的最大键数"只有 0、1、5 三档)。
        let smi = "[Ni]123N4CCN1CCCN2CCN3CCC4";
        let m = prep(smi);
        let ranks = ranks_of(&m);
        assert!(
            hapto_extras(&m, &ranks).is_none(),
            "{smi} 是 σ 给体的大环螯合,不是 η 配位,一根键都不该摘"
        );
        // 前提要自己成立:金属确实对这个环有 ≥3 根键,只是**不连续**。
        // 少了这一句,把判据换成"从来不摘"也照样绿。
        let ni = (0..u32::try_from(m.num_atoms()).unwrap())
            .find(|a| m.atoms()[*a as usize].atomic_num == 28)
            .expect("该有一个镍");
        assert_eq!(m.degree(ni), 4, "镍该是四配位");
        for style in &Style::ALL {
            let d = generate(&m, style);
            assert!(
                d.crossings.is_empty(),
                "[{}] {smi} 该画得出 0 交叉,实得 {} 处",
                style.name,
                d.crossings.len()
            );
        }
    }

    #[test]
    fn contiguity_is_what_separates_a_sandwich_from_a_chelate() {
        // 纯函数,直接钉住。**这是上一条判据背后的那把尺子。**
        let s = |v: &[u32]| {
            v.iter()
                .copied()
                .collect::<std::collections::BTreeSet<u32>>()
        };
        let ring: Vec<u32> = (0..6).collect();
        assert!(contiguous_on_ring(&ring, &s(&[0, 1, 2])), "连续的三个");
        assert!(
            contiguous_on_ring(&ring, &s(&[4, 5, 0])),
            "跨过接头也算连续"
        );
        assert!(contiguous_on_ring(&ring, &s(&[0, 1, 2, 3, 4, 5])), "整圈");
        assert!(
            !contiguous_on_ring(&ring, &s(&[0, 2, 4])),
            "隔一个的三个不算"
        );
        assert!(!contiguous_on_ring(&ring, &s(&[0, 1, 3])), "两段不算");
        // 五元环整圈 —— 二茂铁就是这一种
        let cp: Vec<u32> = (0..5).collect();
        assert!(contiguous_on_ring(&cp, &s(&[0, 1, 2, 3, 4])));
    }

    #[test]
    fn a_chelate_ring_is_left_alone_because_it_only_exists_through_the_metal() {
        // **这一刀不许伤及无辜。** 螯合物(三乙二胺合钴之类)的环是**经过金属**
        // 才闭合的 —— 把金属的键摘掉之后配体只剩几条链,一个环都没有,于是
        // `hapto_extras` 自然什么也找不到。
        //
        // **这一条守的是"摘键前先摘金属"这个做法本身**,不是那个 `>= 3` 的阈值
        // —— 阈值降到 2 这条也照样绿(变异验过)。真正的区别在于:η 配位的环
        // 在**没有金属时就存在**(Cp⁻ 自己就是个五元环),螯合环不然。
        //
        // 那个 `>= 3` 的阈值**在本语料上没有能区分它的例子**:实测把口径放到
        // `>= 2` 重扫 8831 个分子,"金属与配体环有正好 2 根键"的情形是 **0 个**。
        // 所以它是个设计取舍,不是被数据逼出来的 —— 如实说,不编一个例子来
        // 假装它有判据。
        for smi in [
            "C1CN[Co]23(N1)(NCCN2)NCCN3",
            "[O-]S([O-])(=O)=O.C1CN[Cr+3]23(N1)(NCCN2)NCCN3",
            "CC1=[O+][Co]23([O+]=C(C)C1)([O+]=C(C)CC(=[O+]2)C)[O+]=C(C)CC(=[O+]3)C",
        ] {
            let m = prep(smi);
            let ranks = ranks_of(&m);
            assert!(
                hapto_extras(&m, &ranks).is_none(),
                "{smi} 里没有 η 配位,不该摘任何键"
            );
            // 而且这些分子确实**有**金属、金属确实**有**多根键 —— 否则上面那句
            // 是空过的
            assert!(
                (0..u32::try_from(m.num_atoms()).unwrap())
                    .any(|a| m.atoms()[a as usize].atomic_num > 20 && m.degree(a) >= 4),
                "{smi} 里该有一个多配位的金属,不然这条判据说明不了问题"
            );
        }
    }

    #[test]
    fn which_hapto_bond_survives_does_not_depend_on_how_it_was_written() {
        // 每个 (金属, 环) 只留一根代表键,**留哪一根按规范秩定**。改用存储下标
        // 的话,同一个二茂铁换种写法就会留下另一根 —— 金属挂到环上另一个位置,
        // 整张图跟着转。变异验过:这条当场红,而别的判据全绿。
        // **第二种写法由规范式回写产生,不是手写的。** 自己编 SMILES 会把编错
        // 的风险带进结论 —— 本轮已经踩过一次(为查顺反造的"最小复现"
        // `C/C=C/C` vs `C(/C)=C/C` 其实是 E 与 Z 两个不同分子)。
        //
        // # 这条判据抓不到"代表键改用存储序"
        //
        // 实测:把 `into.sort_unstable()`(按规范秩)换成按存储下标,这条**照样
        // 绿**,全量语料也一处不变。原因是 Cp 环五重对称 —— 金属挂在哪个碳上
        // 画出来都全等,`orient` 再把姿态归一,于是观测不到差别。
        //
        // 那句排序因此是**预防性的**:语料里没有能区分它的例子(带取代基的那个
        // 二茂铁也不行,取代基在另一个环上)。**如实说,不编一个分子来假装它有
        // 判据。** 这条判据守的是别的东西 —— 摘键这件事本身不引入写法依赖。
        for base in [
            "C12C3=C4C5=C1[Fe]23456789C%10C6=C7C8=C9%10",
            "CN(C)C[C-]12C3=C4C5=C1[Fe++]23456789[C-]%10C6=C7C8=C9%10",
        ] {
            let m = prep(base);
            let other = omgkit_io::canon::canonical_smiles(&m).smiles;
            // 先证明它真的换了存储序,否则这条判据是空过的。
            //
            // **要比键表,不能只比元素序** —— 二茂铁除了 Fe 全是碳,两种写法的
            // 元素序碰巧一模一样,拿它当判据是空过的(踩过)。
            let seq = |s: &str| -> Vec<(u32, u32)> {
                prep(s).bonds().iter().map(|b| (b.begin, b.end)).collect()
            };
            assert_ne!(
                seq(base),
                seq(&other),
                "{base} 与它的规范式存储序一样,这条判据验不了东西"
            );
            for style in &Style::ALL {
                assert_eq!(
                    scene_key(base, style),
                    scene_key(&other, style),
                    "[{}] {base}\n  写成 {other} 之后画出来的图元不同",
                    style.name
                );
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
