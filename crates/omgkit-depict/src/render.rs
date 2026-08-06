//! 几何 → 图元。**与输出格式无关的那一半**。
//!
//! 把坐标、键级、标签变成一组画得出来的基本形状(线段、文字、楔形、虚线阶梯),
//! 单位是**磅(pt)**,已按 [`Style::bond_length_pt`] 缩放过。
//!
//! 后端只负责把图元序列化成 SVG / 位图 —— 几何只算这一遍。加新后端不用重写
//! 任何几何,而几何出了错,所有后端一起错、一起被同一批判据抓住。
//!
//! # 三处必须做对的事
//!
//! **y 轴要翻。** 化学图的 y 向上,SVG/位图的 y 向下。不翻的话整张图是上下
//! 镜像的 —— 拓扑、键长全对,只有手性看着反了,而那正是最难发现的一类错。
//!
//! **画双键前要在副本上凯库勒化。** 净化之后芳香键的 `order` 是 `Aromatic`,
//! 直接按它画就成了一堆单线,苯环画成六边形圈。
//!
//! **键要在标签外停住。** 线画进标签里会把字盖掉。停多远由
//! [`Style::margin_width_pt`] 定 —— 这是规范的一部分。

use omgkit_chem::sssr::Ring;
use omgkit_core::{BondFlags, BondOrder, MolBuilder};

use crate::geom::Point2;
use crate::label::{label_for, HSide, Label, LabelDir, Run};
use crate::style::Style;
use crate::Depiction;

/// 一个画得出来的基本形状。坐标单位是**磅**,y 向下。
#[derive(Debug, Clone, PartialEq)]
pub enum Primitive {
    /// 线段
    Line {
        /// 起点
        from: Point2,
        /// 终点
        to: Point2,
        /// 线宽
        width: f64,
    },
    /// 实楔形:`from` 是窄端(立体中心),`to` 是宽端
    Wedge {
        /// 窄端
        from: Point2,
        /// 宽端
        to: Point2,
        /// 宽端的宽度
        wide: f64,
    },
    /// 虚楔形:一叠垂直于键的短横线,从窄到宽
    Hash {
        /// 窄端
        from: Point2,
        /// 宽端
        to: Point2,
        /// 宽端的宽度
        wide: f64,
        /// 横线间距
        spacing: f64,
        /// 线宽
        width: f64,
    },
    /// 一个原子的标签,`at` 是**中心**
    Text {
        /// 中心
        at: Point2,
        /// 文本段
        runs: Vec<Run>,
        /// 正文字号
        size: f64,
    },
}

/// 整张图的图元,连同画布尺寸。
#[derive(Debug, Clone, PartialEq)]
pub struct Scene {
    /// 图元
    pub items: Vec<Primitive>,
    /// 画布宽(pt)
    pub width: f64,
    /// 画布高(pt)
    pub height: f64,
}

/// 画布四周留白,单位是磅。让线条不贴着边缘。
const PAD_PT: f64 = 8.0;

/// 把一张 [`Depiction`] 变成图元。
///
/// # Panics
///
/// `depiction` 与 `mol` 的原子数不符时 panic —— 那说明拿错了图,继续画只会
/// 得到一张张冠李戴的结构式。
#[must_use]
pub fn scene(mol: &MolBuilder, depiction: &Depiction, style: &Style) -> Scene {
    // **画的是补完氢的那个分子。** 有些立体中心三根键全在环上,唯一合法的楔形
    // 是 C–H,那个氢是 `generate` 补出来的(见 [`hydrogens`](crate::hydrogens))。
    // 拿传进来的分子画的话,补出来的氢连同它那根楔形会被**静默丢掉** —— 而
    // `unwedged` 是空的,诊断全绿。
    let grown = depiction.drawn(mol);
    let mol = &*grown;
    assert_eq!(
        depiction.coords.len(),
        mol.num_atoms(),
        "这张图不是这个分子的:坐标 {} 个,原子 {} 个",
        depiction.coords.len(),
        mol.num_atoms()
    );

    // 芳香键要按凯库勒式画,否则苯环只剩一圈单线。只取键级,不动调用方的分子。
    let orders = drawn_orders(mol);

    let scale = style.bond_length_pt;
    let bnd = bounds(&depiction.coords, mol, style);
    let (min_x, min_y, max_x, max_y) = bnd;
    let to_pt = |p: Point2| to_canvas(p, bnd, scale);
    // 原子的画布坐标算一次就够。**双键偏向哪一侧也在这个坐标系里算** ——
    // 见 [`offset_dir`] 的文档:两个坐标系一混,横着的键就会画反。
    let pts: Vec<Point2> = depiction.coords.iter().map(|p| to_pt(*p)).collect();
    let rings = omgkit_chem::sssr::ring_set(mol);

    let labels: Vec<Option<Label>> = (0..mol.num_atoms())
        .map(|i| {
            let a = u32::try_from(i).expect("原子数超出 u32");
            label_at(mol, a, style, &depiction.coords)
        })
        .collect();

    // 标签的字形盒,画布坐标系。**盒心不在原子上** —— 整串朝一侧挪了 `dx`,
    // 好让元素符号落在原子位置上,见 [`Label::dx`](crate::label::Label::dx)。
    let boxes: Vec<Option<(Point2, f64, f64)>> = labels
        .iter()
        .enumerate()
        .map(|(i, l)| {
            l.as_ref().map(|l| {
                (
                    pts[i] + Point2::new(l.dx * scale, 0.0),
                    l.half_w * scale,
                    l.half_h * scale,
                )
            })
        })
        .collect();
    let margin_pt = style.margin() * scale;

    let mut items = Vec::new();

    for (bi, b) in mol.bonds().iter().enumerate() {
        // **楔形的窄端必须在它描述的那个立体中心。** 键的 begin/end 与谁是中心
        // 无关,方向搞反的话楔形从取代基指回中心 —— 读出来是另一个构型,而线条
        // 本身看着没毛病。
        //
        // 窄端由 [`Wedge`](crate::stereo::Wedge) 自己带着,不在这里猜。先前是按
        // "哪头带手性标记"猜的,**两头都是立体中心时(相邻的两个中心共用一根键)
        // 就会猜到前一头去**。
        let wedge = depiction.wedges.get(bi).copied().unwrap_or_default();
        let flip = wedge.narrow() == Some(b.end);
        let (bg, en) = if flip {
            (b.end, b.begin)
        } else {
            (b.begin, b.end)
        };
        let (pa, pb) = (depiction.coords[bg as usize], depiction.coords[en as usize]);
        // 在标签外停住,避免线条盖住字
        let (qa, qb) = trim(pa, pb, &labels[bg as usize], &labels[en as usize], style);
        let (a, bb) = (to_pt(qa), to_pt(qb));
        let w = style.line_width_pt;
        // 从主线**横向平移**出来的那些线(双键的第二条、叁键的两条外侧线),
        // 平移会把端点带回盒里 —— 见 [`escape_boxes`]。主线不走这一步,它归
        // [`trim`] 管,连"塞不下时压缩"那档兜底一起。
        let sidelined = |from: Point2, to: Point2| {
            let (from, to) =
                escape_boxes(from, to, boxes[bg as usize], boxes[en as usize], margin_pt);
            Primitive::Line { from, to, width: w }
        };

        match orders[bi] {
            BondOrder::Double => {
                // `off` 与 `n` 都在画布坐标系里 —— 混用两个坐标系正是先前那个
                // "横着的双键画到环外"的成因,见 [`offset_dir`]。
                let spacing = style.bond_spacing() * scale;
                let off = offset_dir(mol, bi, &pts, &rings, &orders, &labels, spacing);
                let d = (bb - a).normalized();
                let n = Point2::new(-d.y, d.x) * spacing;
                let n = if n.dot(off) < 0.0 { n * -1.0 } else { n };
                if off.norm() < 1e-9 {
                    // 没有偏向的一侧(孤立双键):两条线对称分布
                    let half = n * 0.5;
                    items.push(sidelined(a + half, bb + half));
                    items.push(sidelined(a - half, bb - half));
                } else {
                    items.push(Primitive::Line {
                        from: a,
                        to: bb,
                        width: w,
                    });
                    // 第二条线两端**斜切**到与相邻键接上,见 [`mitre_end`]。
                    //
                    // 有两种端原子不斜切,跟着主线的端点平移就是了 ——
                    // **那一头本来就没有可接的东西**:
                    //
                    // - **带标签的**:斜切点是从原子中心量的角平分线交点,而主线
                    //   为了让开字已经缩回去了,两者不在一个起跑线上;算出来的
                    //   斜切点直接落在字里,还比主线伸得更靠前。RDKit 同规则,
                    //   `doubleBondEnd` 的 `trunc` 取 `!atomLabels_[at]`。
                    // - **度 1 的**:根本没有第二根键。丙烯 `CH₂=CH–CH₃` 末端那
                    //   头就是,两条线齐头收尾才是端烯的通例;RDKit 的
                    //   `doubleBondTerminal` 在那一头也是纯法向平移、不缩进。
                    //
                    // 斜切还取不到(邻居都在另一侧、角平分线退化)时才退回按固定
                    // 比例缩进。
                    let fallback = (bb - a) * 0.12;
                    let end = |e: u32, o: u32, p: Point2, back: Point2| {
                        if labels[e as usize].is_some() || mol.degree(e) == 1 {
                            p + n
                        } else {
                            mitre_end(mol, e, o, &pts, n).unwrap_or(p + n + back)
                        }
                    };
                    let from = end(bg, en, a, fallback);
                    let to = end(en, bg, bb, fallback * -1.0);
                    items.push(sidelined(from, to));
                }
            }
            BondOrder::Triple => {
                let d = (bb - a).normalized();
                let n = Point2::new(-d.y, d.x) * (style.bond_spacing() * scale);
                items.push(Primitive::Line {
                    from: a,
                    to: bb,
                    width: w,
                });
                items.push(sidelined(a + n, bb + n));
                items.push(sidelined(a - n, bb - n));
            }
            _ => match wedge {
                crate::stereo::Wedge::Up { .. } => items.push(Primitive::Wedge {
                    from: a,
                    to: bb,
                    wide: style.bold_width_pt,
                }),
                crate::stereo::Wedge::Down { .. } => items.push(Primitive::Hash {
                    from: a,
                    to: bb,
                    wide: style.bold_width_pt,
                    spacing: style.hash_spacing_pt,
                    width: style.line_width_pt,
                }),
                crate::stereo::Wedge::None => items.push(Primitive::Line {
                    from: a,
                    to: bb,
                    width: w,
                }),
            },
        }
    }

    for (i, l) in labels.iter().enumerate() {
        if let Some(l) = l {
            // 整串朝一侧挪 `dx`,让**元素符号**落在原子位置上
            items.push(Primitive::Text {
                at: to_pt(depiction.coords[i]) + Point2::new(l.dx * scale, 0.0),
                runs: l.runs.clone(),
                size: style.atom_label_pt,
            });
        }
    }

    Scene {
        items,
        width: (max_x - min_x) * scale + 2.0 * PAD_PT,
        height: (max_y - min_y) * scale + 2.0 * PAD_PT,
    }
}

/// 逐键的**画图用键级**:芳香键还原成交替单双键。
///
/// 调用方想知道"这张图里哪根键画成了双键"时也用它 —— 自己再凯库勒化一遍会
/// 挑到另一套单双键,那就不是图上画的那一套了。
///
/// # 为什么要先按规范秩重排一遍
///
/// 凯库勒化要在若干套等价的单双键里挑一套,挑哪一套取决于原子的**存储顺序**。
/// 于是同一个分子换一种 SMILES 写法,苯环上三根双键的位置整个换一圈 ——
/// 坐标逐点相同,画出来的线却不同,而[本 crate 的硬性要求](crate)是
/// **一张图由分子自己决定**。实测:阿司匹林、萘、吡啶都会变。
///
/// 所以先把原子按规范秩重排、键按新编号排序,在这个副本上凯库勒化,再把键级
/// 映回原来的键号。规范秩与写法无关,于是选出的那套单双键也与写法无关。
///
/// 凯库勒化失败(比如芳香体系里有通配原子)时原样返回分子自己的键级 ——
/// 那时芳香环会画成一圈单线,难看但不假装成功。
#[must_use]
pub fn drawn_orders(mol: &MolBuilder) -> Vec<BondOrder> {
    let plain = || mol.bonds().iter().map(|b| b.order).collect::<Vec<_>>();
    let ranks = omgkit_io::canon::canonical_ranks(mol);
    let mut order: Vec<usize> = (0..mol.num_atoms()).collect();
    order.sort_by_key(|i| (ranks[*i], *i));
    let mut pos = vec![0u32; mol.num_atoms()];
    for (new, old) in order.iter().enumerate() {
        pos[*old] = u32::try_from(new).expect("原子数超出 u32");
    }

    let mut copy = MolBuilder::with_capacity(mol.num_atoms(), mol.num_bonds());
    for old in &order {
        copy.add_atom_data(mol.atoms()[*old]);
    }
    // 键也要重排:凯库勒化的搜索顺序同样看键的存储序
    let mut bs: Vec<(u32, u32, usize)> = mol
        .bonds()
        .iter()
        .enumerate()
        .map(|(i, b)| {
            let (x, y) = (pos[b.begin as usize], pos[b.end as usize]);
            (x.min(y), x.max(y), i)
        })
        .collect();
    bs.sort_unstable();
    for (x, y, i) in &bs {
        let mut bd = mol.bonds()[*i];
        bd.begin = *x;
        bd.end = *y;
        // 顺反的参照原子也是原子号,一并映过去 —— 副本只用来取键级,但留着
        // 旧编号的话,谁哪天多读一个字段就会读到别的原子上
        for r in &mut bd.stereo_atoms {
            if *r != omgkit_core::BondData::NO_STEREO_ATOM {
                *r = pos[*r as usize];
            }
        }
        if copy.add_bond_data(bd).is_err() {
            return plain();
        }
    }
    if omgkit_chem::kekulize(&mut copy).is_err() {
        return plain();
    }

    let mut out = plain();
    for (k, (_, _, i)) in bs.iter().enumerate() {
        out[*i] = copy.bonds()[k].order;
    }
    out
}

/// 布局坐标 → 画布坐标。**y 在这里翻**:化学图 y 向上,画布 y 向下。
///
/// `bnd` 是 [`bounds`] 的返回值。单独抽出来是为了让判据能用同一套映射 ——
/// 判据自己再写一遍的话,写错的方式可以和实现一模一样,就守不住了。
pub fn to_canvas(p: Point2, bnd: (f64, f64, f64, f64), scale: f64) -> Point2 {
    Point2::new(
        (p.x - bnd.0) * scale + PAD_PT,
        (bnd.3 - p.y) * scale + PAD_PT,
    )
}

/// 含标签的包围盒,单位是**键长**。
/// 画布的包围盒:`(min_x, min_y, max_x, max_y)`,单位是键长。
///
/// **判据必须用这个函数,不许自己抄一份。** 先前审计里 `canvas_pts` 按"同样的
/// 规则"重算了一遍;等这边改成"用真正的 `h_side`、并把标签的横向偏移 `dx` 算
/// 进去"之后,那份副本没跟着改,判据算出来的原子位置与 `scene` 画出来的线对
/// 不上号 —— `环内双键` 从 0 违例变成 **1403**,全是定位错造成的假阳。
pub fn bounds(coords: &[Point2], mol: &MolBuilder, style: &Style) -> (f64, f64, f64, f64) {
    let (mut x0, mut y0, mut x1, mut y1) = (f64::MAX, f64::MAX, f64::MIN, f64::MIN);
    for (i, p) in coords.iter().enumerate() {
        let a = u32::try_from(i).expect("原子数超出 u32");
        // **要用真正会画出来的那个标签。**
        //
        // 三处先前都错着:写死 `HSide::Right`(氢挂哪边其实看坐标)、完全不看
        // `dx`(整串朝一侧挪了之后**盒心不在原子上**,画布按"盒心在原子上"算
        // 就会短一截),以及只问 `label_for` —— 共线的骨架碳是 [`label_at`] 补
        // 出来的符号,`scene` 画它,画布却没给它留地方。三处各实测踩到过
        // 2、2、4 处,`不出画布` 那条硬性质当场破。
        let (dx, hw, hh) =
            label_at(mol, a, style, coords).map_or((0.0, 0.0, 0.0), |l| (l.dx, l.half_w, l.half_h));
        x0 = x0.min(p.x + dx - hw);
        y0 = y0.min(p.y - hh);
        x1 = x1.max(p.x + dx + hw);
        y1 = y1.max(p.y + hh);
    }
    if x0 > x1 {
        return (0.0, 0.0, 0.0, 0.0);
    }
    (x0, y0, x1, y1)
}

/// 这个原子**真正会画出来**的标签。
///
/// 与 [`label_for`] 的差别是共线的二度原子:骨架碳本来不画符号,但两根键连成
/// 一条直线时顶点处没有拐角,图上根本看不出那里有个原子,所以补一个 ——
/// 见 [`is_collinear`]。
///
/// [`scene`] 画什么、[`bounds`] 按什么留白、判据按什么重建,**必须是同一个
/// 来源**。三处各写各的话,补出来的那个符号就会有的地方算、有的地方不算 ——
/// 实测就是这样破的 `不出画布`:`scene` 给共线的骨架碳画了 `CH`,`bounds`
/// 不知道,半宽 7.22pt 的字直接戳出画布右边(全量 4 处)。
///
/// **有一处没跟过来,要说清楚:**[`refine::radii`](crate::refine) 仍只问
/// `label_for`,共线的骨架碳在那里按裸原子的半径 0.25 个键长算,而它将被画出的
/// `CH` 盒半径是 0.559 —— 布局给它留的地方不到实际的一半。这不是漏改:消冲突
/// 跑的时候坐标还在动,"共线不共线"当时判不了,`label_at` 在那个阶段没有定义。
pub fn label_at(mol: &MolBuilder, atom: u32, style: &Style, coords: &[Point2]) -> Option<Label> {
    let side = h_side(mol, atom, coords);
    label_for(mol, atom, style, side).or_else(|| {
        is_collinear(mol, atom, coords).then(|| crate::label::label_forced(mol, atom, style, side))
    })
}

/// 氢挂哪一侧:挂在**键伸出去的反方向**,免得氢和键叠在一起。
///
/// # 平局必须用容差判,不能直接比 0
///
/// 竖直的键两端 x 本该完全相等,浮点算出来却差个 1e-17 量级的零头。拿
/// `sum > 0.0` 去判,这个零头的符号就决定了写 `OH` 还是 `HO` —— 而它的符号
/// 取决于算到那一步的运算次序。实测:乙醇的羟基因此被写成了 `HO`。
///
/// 这不只是难看,**是确定性隐患**:同一分子的不同写法可能得到两种标签。
/// 所以近似打平时一律取 [`HSide::Right`] —— `OH`、`NH2` 是常规写法,
/// `HO`、`H2N` 只在键从右边来时才用。
///
/// 判据要重建 `scene` 画出来的东西就得知道氢挂哪边,所以这个函数是公开的。
///
/// # 为什么不写成 `label_dir(..).h_side()`
///
/// 竖排还没接进绘制。在那之前把上下两向折成 `East`,等于**丢掉**这批标签本来
/// 定得好好的左右选择 —— 全量语料实测 `—— 有标签在键上塞不下` 从 1271 涨到
/// **1356(+85)**。所以这里就是老口径本身:只看 x 分量。上下两向要等
/// [`crate::label::Label`] 会竖排了才接得上。
///
/// 阈值用的是模块里那一份 `TIE`,**不另抄一个字面量** —— 抄一份就会漂,而漂了
/// 之后没有判据守得住(实测:把这里的 `1e-3` 改成 `0.0`,116 条判据全绿)。
pub fn h_side(mol: &MolBuilder, atom: u32, coords: &[Point2]) -> HSide {
    // **度 0 例外。** 老口径的 x 分量恒为 0,一律落到 `Right`,于是水画成
    // `OH2`、氯化氢画成 `ClH`。这一支照 [`label_dir`]:那批元素写在氢后面。
    // 语料里度 0 原子一个都没有,所以这个例外不动任何已量过的数。
    if mol.degree(atom) == 0 {
        return match label_dir(mol, atom, coords) {
            LabelDir::West => HSide::Left,
            _ => HSide::Right,
        };
    }
    if nbr_sum(mol, atom, coords).0 > TIE {
        HSide::Left // 键都朝右,氢写左边
    } else {
        HSide::Right
    }
}

/// 键向量之和:`Σ (邻居 − 自己)`。标签往它的**反方向**伸。
fn nbr_sum(mol: &MolBuilder, atom: u32, coords: &[Point2]) -> (f64, f64) {
    let here = coords[atom as usize];
    let (mut sx, mut sy) = (0.0, 0.0);
    for (n, _) in mol.neighbors(atom) {
        sx += coords[n as usize].x - here.x;
        sy += coords[n as usize].y - here.y;
    }
    (sx, sy)
}

/// 左右打平的容差。小于这个量的横向偏移一律当作打平 —— 取键长的千分之一,
/// 真正的左右之别至少是半个键长的量级,不可能落进来。
const TIE: f64 = 1e-3;

/// 判"这堆键主要是横的还是竖的"用的斜率阈值:tan 70°。
///
/// # 70° 是抄来的,而且这个数不能随便动
///
/// RDKit 的 `DrawMol::getAtomOrientation` 用的就是它,注释里写明了取值的由来:
/// 吲哚的 `NH` 在它的布局下约 72°,要竖着;而 `c1ccccc1C1CCC(N)(N)CC1` 底部
/// 那两个氨基要横着。70° 把这两种情形分开。
///
/// **换成 30° 会踩上一个结构性的坑。** 本库的图被 [`crate::orient`] 吸附到 30°
/// 栅格上;查阈值的是度 ≥ 2 的原子,它的键向量和是两个以上栅格向量之和,方向
/// 落在 **15° 栅格**上,`|tan|` 的取值是 `{0, 0.268, 0.577, 1, 1.732, 3.732, ∞}`
/// —— **0.577 正是 tan 30°**,于是大批键向量和会**精确落在阈值上**,横竖之分
/// 退化成由浮点末位决定。tan 70° = 2.747 不在这个集合里,离最近的 1.732(60°)
/// 与 3.732(75°)都有余量。全量语料实测:阈值取 70° 时只有 **70** 个原子落在
/// 断点上。
///
/// # 为什么写死一个字面量,不在运行时算
///
/// `tan` 不是正确舍入的,跨平台、跨 libm 版本能差 ulp。本库的头号契约是同一
/// 分子的任何写法画出**逐字节相同**的图元,阈值随平台漂就守不住。写死的这个
/// 是 tan 70° 的最近双精度数(本机 libm 给的比它低 2 ulp,而两簇之间隔着 8 个
/// 数量级,差这 2 ulp 不改任何结论)。[`a_seventy_degree_threshold_is_what_is_written_down`]
/// 这条判据钉住它。
///
/// [`a_seventy_degree_threshold_is_what_is_written_down`]: #
pub const VERT_SLOPE: f64 = 2.747_477_419_454_622; // tan(70°)

/// 判横竖时的容差。
///
/// # 两簇之间空着 12 个数量级
///
/// 判别量取 `|sy| - tan70·|sx|`(不是 `sy/sx`,竖直时那个会溢出)。全量语料
/// 17662 个分子×规范、19694 个带氢标签实测:
///
/// | | 个数 | 判别量 |
/// |---|---:|---|
/// | 精确落在断点上 | 70(其中度 ≥ 3 的 58 个) | ≤ 4.658e-15 |
/// | 真正有横竖之分 | 19624 | ≥ 4.461e-3 |
///
/// 中间空着 12 个数量级,1e-9 落在正中间,两边都够不着。
///
/// **度 ≥ 3 那 58 个是真的破口**:`sum` 是按 `mol.neighbors` 的存储序累加的,
/// 而三项以上的浮点加法不满足结合律 —— 不定死的话,同一分子的两种写法会在
/// 横竖之间摇摆。度 1、度 2 不受影响(单项、两项加法与次序无关)。
pub const DIR_TIE: f64 = 1e-9;

/// 标签往哪个方向伸:**背离键伸出去的方向**,免得压在键上。
///
/// # 端基永远横排
///
/// 度 1 的原子一律给 `East`/`West`。"不竖排"这条照 RDKit(它的注释原话是
/// "atoms of single degree should always be either W or E, never N or S")——
/// 端基只有一根键,横着写读起来就是 `–OH`、`–NH2` 这些常规写法;竖排既没必要,
/// 又会把一串端基排得高低不齐。
///
/// **但左右怎么选,本库与 RDKit 故意不同。** RDKit 把已判成竖的度 1 原子一律
/// 改成 `East`;本库仍按 `sx` 的符号定左右。实测两者在 **805** 个标签上不同
/// (键比 70° 还陡、且 `|sx| > TIE` 的那些,本库给 `HO`/`H2N`,RDKit 给
/// `OH`/`NH2`)。
///
/// # 键向量和整体近似为零时,本库判横、RDKit 判竖
///
/// RDKit 在 `|nbr_sum.x| <= 1e-4` 时把斜率强行设成 1000,于是共线的二度原子、
/// 对称的四度原子都被判成竖排。本库的 [`DIR_TIE`] 让这类落到横排 —— 共线原子
/// 判竖排会把标签盖到键上。实测涉及 **58** 个带氢标签,**故意不同**。
///
/// # 度 3 的特判没有实现
///
/// RDKit 对度 3 还有一支:某根键精确竖直时,按它的方向重定 N/S,免得氢压在
/// 那根键上。本库没实现。实测这条规则会改 **144** 个原子的朝向,其中 **14** 个
/// 有标签、**0** 个标签含氢 —— 现在一个氢都不受影响,所以先不补。等电荷、
/// 同位素也跟着竖排时要重新评估。真要补也**不能照抄**:RDKit 那段用
/// `atan(y/x)` 分不开上下,而且"取第一根竖直的键"是拿存储下标打破平局,
/// 直接违反写法无关。
///
/// # 打平时取 `East`
///
/// 竖直的键两端 x 本该完全相等,浮点算出来却差个 1e-17 量级的零头。拿
/// `sum > 0.0` 去判,这个零头的符号就决定了写 `OH` 还是 `HO` —— 而它的符号
/// 取决于算到那一步的运算次序。实测:乙醇的羟基因此被写成了 `HO`。
/// 所以近似打平时一律取 `East` —— `OH`、`NH2` 是常规写法,`HO`、`H2N` 只在键
/// 从右边来时才用。横竖打平时同理取横,见 [`DIR_TIE`]。
///
/// # `coords` 必须是**布局坐标**(y 向上)
///
/// `Point2` 在本库同时用于布局系与画布系。传画布坐标进来,左右两向照样对,
/// **上下两向会静默反过来**。返回的 `North`/`South` 已经是**页面**口径,
/// `to_canvas` 之后不要再翻第二次。
///
/// 判据要重建 `scene` 画出来的东西就得知道标签朝哪,所以这个函数是公开的。
pub fn label_dir(mol: &MolBuilder, atom: u32, coords: &[Point2]) -> LabelDir {
    // **度 0 单独一支。** 没有键就没有"背离键"这回事,`nbr_sum` 是零向量,
    // 落到通用分支会一律给 `East`,于是水画成 `OH2`、氯化氢画成 `ClH`。
    // 这批元素照 RDKit(`DrawMol.cpp` 的 `getAtomOrientation` else 支)给 `West`,
    // 写出来是 `H2O`、`HCl`。**全量语料查不出这条** —— `large.smi` 里度 0 原子
    // 一个都没有,而水合物、盐的对离子在真实输入里很常见。
    if mol.degree(atom) == 0 {
        return match mol.atoms()[atom as usize].atomic_num {
            8 | 9 | 16 | 17 | 34 | 35 | 52 | 53 | 84 | 85 => LabelDir::West,
            _ => LabelDir::East,
        };
    }
    let (sx, sy) = nbr_sum(mol, atom, coords);
    // 横竖:比 `|sy|` 与 `tan70·|sx|`。写成差而不是商,竖直时才不会溢出。
    let vertical = mol.degree(atom) >= 2 && sy.abs() - VERT_SLOPE * sx.abs() > DIR_TIE;
    if vertical {
        // 键朝上就把标签甩到下面去
        if sy > 0.0 {
            LabelDir::South
        } else {
            LabelDir::North
        }
    } else if sx > TIE {
        LabelDir::West // 键都朝右,氢写左边
    } else {
        LabelDir::East
    }
}

/// 把键的两端各缩短一点,停在标签外。
/// 双键内侧线在 `e` 这一端该收到哪里:与**相邻键的角平分线**求交。
///
/// # 为什么不能按固定比例缩进
///
/// 先前两端各缩固定的 12%。那个数只在某一种夹角下恰好对 —— 苯环上正确的
/// 缩进量是 `s / (2 × 边心距) = 0.18 / (2 × 0.866) = 10.4%`,差 1.6 个百分点,
/// 内侧线的端点因此偏离顶点角平分线 **0.2 pt**,接头处露白。夹角越偏离
/// 120°,差得越多。
///
/// 做法与 RDKit 的 `DrawMol::doubleBondEnd` 一致:取 `e` 的一个邻居 `t`
/// (在偏移那一侧),`e` 处的角平分线是 `‖e→t‖ + ‖e→o‖` 的方向,内侧线是过
/// `e+n`、方向 `e→o` 的直线,两者的交点就是端点。规则多边形上这恰好落在
/// "环心 → 顶点"那条射线上。
///
/// # 挑哪个邻居不能看存储下标
///
/// `e` 可能有好几个邻居在偏移那一侧。挑序只要沾上存储下标,同一个分子换种
/// 写法内侧线就会切在不同的地方 —— 而写法无关是头号契约。这里按**量化后的
/// 投影值**排,再拿量化坐标兜底平局;坐标此时已经规范化过,与写法无关。
fn mitre_end(mol: &MolBuilder, e: u32, o: u32, pts: &[Point2], n: Point2) -> Option<Point2> {
    let pe = pts[e as usize];
    let po = pts[o as usize];
    let to_o = (po - pe).normalized();
    if to_o.norm() < 1e-9 {
        return None;
    }
    // 偏移那一侧的邻居,投影最大的那个
    #[allow(clippy::cast_possible_truncation)]
    let t = mol
        .neighbors(e)
        .map(|(x, _)| x)
        .filter(|x| *x != o)
        .filter_map(|x| {
            let d = pts[x as usize] - pe;
            if d.norm() < 1e-9 {
                return None;
            }
            let dir = d.normalized();
            // 与 n 同侧才算数:异侧的邻居给出的角平分线指向外面
            if dir.dot(n) <= 0.0 {
                return None;
            }
            let p = pts[x as usize];
            Some((
                (dir.dot(n) * 1e9).round() as i64,
                (p.x * 1e9).round() as i64,
                (p.y * 1e9).round() as i64,
                x,
            ))
        })
        .max()?
        .3;

    let to_t = (pts[t as usize] - pe).normalized();
    let bis = to_o + to_t;
    if bis.norm() < 1e-9 {
        return None; // 共线,角平分线退化
    }
    let bis = bis.normalized();
    // 内侧线:过 pe+n,方向 to_o。求 u 使 pe + u·bis 落在这条线上:
    //   (u·bis − n) × to_o = 0  ⟹  u = (n × to_o) / (bis × to_o)
    let denom = bis.cross(to_o);
    if denom.abs() < 1e-9 {
        return None;
    }
    let u = n.cross(to_o) / denom;
    // 切过头的话不如不切 —— 内侧线会反向
    if !u.is_finite() || u <= 0.0 || u > pe.dist(po) {
        return None;
    }
    Some(pe + bis * u)
}

/// 这个二度原子的两根键是不是几乎连成一条直线。
///
/// # 为什么这种原子必须画出符号
///
/// 相邻两根键连成直线时,顶点处没有拐角 —— **图上根本看不出那里有个原子**。
/// 丙二烯 `CH₃CH=C=CHCH₃` 的中心碳是 sp、键角 180°,不画符号的话整张图就是
/// 一条直线加几条平行短线,读起来像顺式二烯,而中心碳无影无踪。
///
/// 判据取几何而不取杂化:`ideal_angle` 已经把 sp 原子摆成 180°,但消冲突之后
/// 别的原子也可能碰巧共线,那时同样看不见。RDKit 的 `isLinearAtom` 是同一个
/// 口径(点积 < −0.95,约 162°),并且同样要求**两根键的键级相同** ——
/// `R—C≡C—R` 的炔碳两侧一单一叁,三条平行线本来就把它标出来了,不用再画符号。
///
/// 判据要数"有多少原子是靠补符号才看得见的",所以这个函数是公开的。
pub fn is_collinear(mol: &MolBuilder, a: u32, coords: &[Point2]) -> bool {
    /// 多共线算共线。取 −0.95 与 RDKit 一致(约 162°)。
    const COS: f64 = -0.95;

    let nbrs: Vec<(u32, u32)> = mol.neighbors(a).collect();
    if nbrs.len() != 2 {
        return false;
    }
    if mol.bonds()[nbrs[0].1 as usize].order != mol.bonds()[nbrs[1].1 as usize].order {
        return false;
    }
    let c = coords[a as usize];
    let u = coords[nbrs[0].0 as usize] - c;
    let v = coords[nbrs[1].0 as usize] - c;
    if u.norm() < 1e-9 || v.norm() < 1e-9 {
        return false;
    }
    u.normalized().dot(v.normalized()) < COS
}

/// 从标签盒的中心沿方向 `d` 走到盒边的距离。`d` 要是单位向量。
///
/// # 为什么不能拿外接圆凑
///
/// 先前切的是 `half_w.hypot(half_h)` —— 标签盒的**外接圆**。圆一定包住盒,所以
/// 线绝不会压到字上,代价是**在盒窄的那个方向上停得太远**:竖直方向去接一个
/// 横向宽的标签,白白空出来的正是 `hypot(w,h) − h`。
///
/// 实测全量语料 129330 个带标签的键端:平均白切 **0.075 个键长**,21.9% 白切
/// 超过 0.1 个键长,最糟的 `[NH2+]` 白切 **0.39 个键长** —— 快四成键长的空白,
/// 一眼就看得出来。
///
/// 标签是 `text-anchor="middle"` 摆的,盒以原子为心,所以这里按居中的轴对齐盒算。
///
/// `from` 是**相对盒心**的位置,要在盒内 —— 盒外的点算出来是到对面那条边的
/// 距离,不是它想要的东西。判据要按同一套口径量净空,所以这个函数是公开的。
pub fn box_reach(from: Point2, d: Point2, half_w: f64, half_h: f64) -> f64 {
    let t = |o: f64, dd: f64, h: f64| {
        if dd.abs() < 1e-12 {
            f64::INFINITY
        } else {
            ((if dd > 0.0 { h } else { -h }) - o) / dd
        }
    };
    t(from.x, d.x, half_w).min(t(from.y, d.y, half_h)).max(0.0)
}

/// 两端要让的净空加起来超过线段的这个比例,就只能按比例压缩。
///
/// [`trim`]、[`escape_boxes`] 与判据的"塞不下"都取这一个数 —— 三处各写一个
/// 字面量的话,改其中一处就会让判据与实现悄悄错位。
const SQUEEZE: f64 = 0.9;

/// 这一端要让出多少净空:从原子沿 `dir` 走出标签盒,再加一个 margin。
///
/// # 判据必须调这个函数
///
/// 净空的口径有两处容易抄错:盒是**偏心**的(整串挪了 `dx`,好让元素符号落在
/// 原子上,见 [`Label::dx`](crate::label::Label::dx)),而且**两端要各算各的**
/// —— 盒不再左右对称,一个值两边通用不成立。
///
/// 审计里 `squeezed_bonds` 就抄过一份**居中盒**的版本:实现改成偏心盒之后判据
/// 没跟上,「有标签在键上塞不下」多报了 547 例(1818 → 1271)。同一个坑更早
/// 还塌过一次(`canvas_pts` 抄 `bounds`,`环内双键` 炸出 1403 处假阳)。
pub fn label_clearance(l: Option<&Label>, dir: Point2, style: &Style) -> f64 {
    l.map_or(0.0, |l| {
        box_reach(Point2::new(-l.dx, 0.0), dir, l.half_w, l.half_h) + style.margin()
    })
}

/// 这根键塞不塞得下两端的标签。
///
/// 塞不下时裁键只能按比例压缩,**端点会落进盒里** —— 那不是切算错了,是
/// ACS 规范下标签本来就占 0.69 个键长,`O⁻—N⁺` 两端要 1.375 个键长的净空,
/// 一个键长装不下。判据据此跳过这类键,口径要与实现完全一致,所以公开。
pub fn is_squeezed(
    pa: Point2,
    pb: Point2,
    la: Option<&Label>,
    lb: Option<&Label>,
    style: &Style,
) -> bool {
    let len = pa.dist(pb);
    if len < 1e-9 {
        return false;
    }
    let d = (pb - pa) * (1.0 / len);
    label_clearance(la, d, style) + label_clearance(lb, d * -1.0, style) >= len * SQUEEZE
}

/// 两端要让的加起来比线段还长时按比例压缩 —— 否则线段会反向,画出一条穿过
/// 标签的短线。剩下的 `1 − SQUEEZE` 那一小截好歹还看得出这里有根键。
fn squeeze(ca: f64, cb: f64, len: f64) -> (f64, f64) {
    if ca + cb >= len * SQUEEZE {
        let k = len * SQUEEZE / (ca + cb).max(1e-9);
        (ca * k, cb * k)
    } else {
        (ca, cb)
    }
}

/// 端点落进标签盒里的话,沿线段自己的方向推出去。盒外的端点原样返回。
///
/// # 为什么 [`trim`] 管不到这些线
///
/// `trim` 切的是**键轴**那条线:从原子中心出发,沿键的方向量到盒边。双键的
/// 第二条线、叁键的两条外侧线都是把它**横向平移**出来的,而横向平移会把端点
/// 带回盒里 —— 盒宽的标签尤其明显:主线从上方出盒时只让开了半高,再横着挪
/// 0.18 个键长,仍然在半宽之内。实测 `HC`(半宽 7.22pt、半高 3.59pt)上
/// 14 处如此。
///
/// 盒是**偏心**的(整串挪了 `dx`),所以按端点相对盒心的位置算,不能按原子算。
///
/// # 盒内/盒外这一刀没有容差,是想清楚了的
///
/// 两个分支差整整一个 margin(ACS 1.6pt、CD 2.0pt),看着像该加个容差。但
/// **容差只会把这个跳变挪个地方,消不掉它** —— "推到盒外"与"不动"之间本来
/// 就差一个 margin。本文件别处的容差(`h_side` 的 `TIE`、`offset_dir` 的
/// `tie`、`mitre_end` 的量化排序)守的是**离散选择**:符号一翻,`OH` 就变成
/// `HO`、第二条线就画到环外,而那个符号取决于求和次序。这里两边给出的都是
/// 合法几何,同一份坐标永远走同一支,写法无关不受影响。
///
/// 实测全量语料上所有平移端点到这条分界的最近距离是 **0.068pt**(ChemDraw 下
/// 芳香 `N⁺` 六元环),比浮点噪声大十几个数量级 —— 目前不在刀锋上。
fn escape_boxes(
    from: Point2,
    to: Point2,
    ba: Option<(Point2, f64, f64)>,
    bb: Option<(Point2, f64, f64)>,
    margin: f64,
) -> (Point2, Point2) {
    let len = from.dist(to);
    if len < 1e-9 {
        return (from, to);
    }
    let d = (to - from) * (1.0 / len);
    let out = |p: Point2, dir: Point2, bx: Option<(Point2, f64, f64)>| -> f64 {
        let Some((c, hw, hh)) = bx else { return 0.0 };
        let r = p - c;
        // **盒外的点不能交给 `box_reach`** —— 它算的是到对面那条边的距离
        if r.x.abs() >= hw || r.y.abs() >= hh {
            return 0.0;
        }
        box_reach(r, dir, hw, hh) + margin
    };
    let (ca, cb) = (out(from, d, ba), out(to, d * -1.0, bb));
    let (ca, cb) = squeeze(ca, cb, len);
    (from + d * ca, to - d * cb)
}

fn trim(
    pa: Point2,
    pb: Point2,
    la: &Option<Label>,
    lb: &Option<Label>,
    style: &Style,
) -> (Point2, Point2) {
    let d = (pb - pa).normalized();
    let (ca, cb) = (
        label_clearance(la.as_ref(), d, style),
        label_clearance(lb.as_ref(), d * -1.0, style),
    );
    let (ca, cb) = squeeze(ca, cb, pa.dist(pb));
    (pa + d * ca, pb - d * cb)
}

/// 双键第二条线该偏向哪一侧。零向量表示"没有偏向",两条线对称画。
///
/// **`pts` 与返回值都是画布坐标**(y 向下)。
///
/// # 坐标系必须说清楚
///
/// 布局坐标 y 向上,画布 y 向下。方向在一个系里算、落点在另一个系里用,
/// 两者的点积化简下来是 `vx² − vy²` 的符号 —— **竖着的键碰巧对,横着的正好
/// 反**,于是"大部分环看着没问题"掩护了整类缺陷。实测:阿司匹林苯环底边
/// 那根双键的第二条线画到了环外。
///
/// # 分几种情形
///
/// | 键 | 第二条线 |
/// |---|---|
/// | 只属于一个环 | 偏向**环内**(射线法判,不是"和环心同侧") |
/// | 稠合处的共用键 | 进**芳香**的那个环;都芳香或都不芳香时取双键最多的。不跨骑 |
/// | 链上、端点共线(累积双键) | 对称跨轴 |
/// | 链上、两端都是端点(乙烯) | 对称跨轴 |
/// | 链上、一端是端点,内侧原子**带标签或还有别的分叉** | 对称跨轴 |
/// | 链上、一端是端点,内侧原子度 2 且无标签 | 偏向内侧原子的另一根键那一侧(端烯、醛) |
/// | 链上、两端都有取代基,票数不为 0 | 偏向取代基多的一侧(顺式那一侧) |
/// | 链上、票数抵消(反式),**两端都带标签** | 对称跨轴 |
/// | 链上、票数抵消(反式),其余 | 按**画布的绝对方向**偏:朝上,竖直键朝左 |
///
/// 后四行与 RDKit `calcDoubleBondLines` / `doubleBondTerminal` 的分档一一对应,
/// 见函数体里的注释。
///
/// 环上那一条**不看取代基**:抗坏血酸环内的 C=C 两端各挂一个 OH,按邻居计数
/// 正好抵消,两条线就骑在环边上,其中一条落到环外去了。
///
/// # 反式那一档的偏侧是**约定**,不是推导
///
/// 两个取代基一边一个时,两侧一样合理。RDKit 按存储序挑(`otherNeighbor(…, 0)`),
/// 那会让同一分子的不同写法画到不同侧 —— 这里改成认**画布的绝对方向**,见
/// 函数体末尾。
fn offset_dir(
    mol: &MolBuilder,
    bi: usize,
    pts: &[Point2],
    rings: &[Ring],
    orders: &[BondOrder],
    labels: &[Option<Label>],
    probe: f64,
) -> Point2 {
    let b = &mol.bonds()[bi];
    let (pa, pb) = (pts[b.begin as usize], pts[b.end as usize]);
    let len = pa.dist(pb);
    if len < f64::EPSILON {
        return Point2::ORIGIN;
    }
    let mid = (pa + pb) * 0.5;
    let axis = (pb - pa) * (1.0 / len);
    let normal = Point2::new(-axis.y, axis.x);
    // 打平的容差按键长取。真正的左右之别至少是半个键长的量级,而共线的邻居
    // 算出来是 1e-15 —— 中间空得很,取哪个数量级都一样。
    let tie = 1e-3 * len;
    // 判"这根键是不是竖着的"的容差。`normal` 是单位向量,所以这是个绝对量。
    // 实测全量语料落进"票数抵消"那一档的 1072 条记录里,`|normal.y|` 只有两簇:
    // 噪声簇 ≤ 1.97e-15、真实簇 ≥ 7.47e-02,**中间空 13 个数量级**,取哪个
    // 数量级都一样。放宽到全部 81488 次调用也一样(5.43e-14 / 6.48e-04)。
    const TIE_DIR: f64 = 1e-9;

    let bond_no = u32::try_from(bi).expect("键数超出 u32");
    let mine: Vec<&Ring> = rings
        .iter()
        .filter(|r| r.bonds.contains(&bond_no))
        .collect();
    if !mine.is_empty() {
        // 第二条线**实际会落在**的两个位置,各问一次"在不在这个环里面"。
        //
        // 先前问的是"环心在哪一侧"。凸环上两者一致,凹环就不一致 —— 卟啉那种
        // 大环、退化布局给出的环,形心可能落在环外,于是判出来的"内侧"其实
        // 是外侧。实测语料里 45 处两者判反。
        let (plus, minus) = (mid + normal * probe, mid - normal * probe);
        // 每个能收下这条线的环记一档,挑一个最该进的
        let mut cands: Vec<(usize, usize, usize, i64, i64, f64)> = Vec::new();
        for r in &mine {
            let poly: Vec<Point2> = r.atoms.iter().map(|a| pts[*a as usize]).collect();
            let side = match (
                crate::geom::point_in_polygon(plus, &poly),
                crate::geom::point_in_polygon(minus, &poly),
            ) {
                (true, false) => 1.0,
                (false, true) => -1.0,
                // 两侧都在或都不在:自交的环、或者窄到放不下这条线的环。弃权。
                _ => continue,
            };
            // **进双键最多的那个环。** 稠合处共用的键两头各是一个环,先前一律
            // 对称画,于是两条线各贴着环边、哪个环都不属于,看着像两根平行的键。
            // 通例是画进"更芳香"的那个环 —— 也就是交替单双键最多的那个,第二条
            // 线于是接上那一圈的交替。
            let doubles = r
                .bonds
                .iter()
                .filter(|b| orders[**b as usize] == BondOrder::Double)
                .count();
            // **芳香的环排在最前。** "双键最多"只是芳香性的代理:kekulé 化之后
            // 芳环的交替双键通常就是最多的那个,所以两条规则几乎总是一致 ——
            // 全量语料 1092 根稠合处的共用双键里只有 3 根分歧。但代理会在
            // "非芳环的双键比隔壁芳环还多"时给出错的答案,而芳香标记是直接的。
            let aromatic = r
                .bonds
                .iter()
                .all(|x| mol.bonds()[*x as usize].flags.contains(BondFlags::AROMATIC));
            let c = poly.iter().fold(Point2::ORIGIN, |s, p| s + *p) * (1.0 / poly.len() as f64);
            #[allow(clippy::cast_possible_truncation)]
            cands.push((
                usize::from(!aromatic),     // 芳环排前
                usize::MAX - doubles,       // 再看双键多的
                r.atoms.len(),              // 并列取小环
                (c.x * 1e6).round() as i64, // 再并列按形心定序 —— 坐标与写法无关
                (c.y * 1e6).round() as i64,
                side,
            ));
        }
        cands.sort_by_key(|a| (a.0, a.1, a.2, a.3, a.4));
        if let Some(best) = cands.first() {
            return normal * best.5;
        }
        // 一个环都收不下(自交的退化环)—— 不猜,交给下面的通用规则
    }

    // **没有"内侧"可言的几种情形,才跨轴对称画。** 与 RDKit
    // `calcDoubleBondLines` / `doubleBondTerminal` 的分档一一对应:
    //
    // - 端点**共线**(累积双键 `C=C=C`):两根键若都把第二条线偏到同一侧,画
    //   出来是一条直线配两条同侧短线 —— 读起来是顺式二烯。RDKit `isLinearAtom`。
    // - **两端都是端基**(乙烯):压根没有别的邻居可偏。RDKit 与上一条同一支。
    if is_collinear(mol, b.begin, pts)
        || is_collinear(mol, b.end, pts)
        || (mol.degree(b.begin) == 1 && mol.degree(b.end) == 1)
    {
        return Point2::ORIGIN;
    }
    // 恰有一端是端基 —— RDKit 的 `doubleBondTerminal`,它**有三支,两支对称**。
    //
    // 先前这里一见度 1 就对称,于是丙烯 `CH₂=CH–CH₃` 的两条线跨着键轴画,而
    // 内侧那个碳上的单键是**收在原子中心**的 —— 两条线谁也不从那一点出发,
    // 顶点合不拢,露出一个豁口。
    //
    // 但**只把这条早退删掉又走过了头**:剩下两支落进下面的投票规则,而"内侧
    // 只有一个别的邻居、票数必然 ±1"对它们同样成立,于是也全变成不对称。实测
    // 多改了 216 根键,方向与 RDKit、与改前都相反:
    //
    // - **内侧带标签**(亚硝基 `R–N=O`、亚胺 `CH₂=N–R`,100 根/95 分子):那一头
    //   所有的键都停在字盒外,压根不存在"单键收在原子中心"这个成因,**没有顶点
    //   要合**。偏一边只丢了对称性,换回来的是 `=` 看着挂在字母底边上。
    // - **内侧还有别的分叉**(丙酮、砜,116 根/68 分子):两侧各被键占着,偏哪边
    //   都压着一根。度 3 时两个邻居分居两侧、票数碰巧抵消,**度 ≥4 抵消不掉** ——
    //   实测砜的两根 S=O 双双朝对方倾斜,四条线挤在中间。
    if mol.degree(b.begin) == 1 || mol.degree(b.end) == 1 {
        let inner = if mol.degree(b.begin) == 1 {
            b.end
        } else {
            b.begin
        };
        if mol.degree(inner) != 2 || labels[inner as usize].is_some() {
            return Point2::ORIGIN;
        }
        // 剩下的是丙烯那一档,落到下面的投票 —— 内侧只有一个别的邻居,
        // 票数必然 ±1(它正好落在键轴上时才是 0,那时本来也无侧可偏)。
    }

    // 两端其它邻居投影到法线上,哪边多就偏哪边。
    //
    // 票数用整数记。用 `f64::signum` 记的话,共线的邻居(累积双键、直线段)
    // 会按 ±0.0 的**符号位**投出一票 —— 而那一位取决于算到那步的运算次序。
    let mut score: i32 = 0;
    for end in [b.begin, b.end] {
        for (n, _) in mol.neighbors(end) {
            if n == b.begin || n == b.end {
                continue;
            }
            let d = (pts[n as usize] - mid).dot(normal);
            if d > tie {
                score += 1;
            } else if d < -tie {
                score -= 1;
            }
        }
    }
    if score != 0 {
        return normal * f64::from(score.signum());
    }

    // 票数抵消 —— **反式双键的通例**,两个取代基一边一个。
    //
    // 对称跨轴是不行的:两端的单键都收在原子中心,而两条线谁也不从那两点出发,
    // **两个顶点都合不拢** —— 与丙烯是同一个毛病,只是这里两头都犯。实测
    // 643 根键 / 574 个分子(单套规范),trans-butene 就是其中一个。
    //
    // **两端都带标签时才对称**:那时键都停在字盒外,压根没有顶点要合(RDKit
    // `calcDoubleBondLines` 里 `atomLabels_[at1] && atomLabels_[at2]` 那一支,
    // 与"端基双键内侧带标签"是同一条道理)。
    if labels[b.begin as usize].is_some() && labels[b.end as usize].is_some() {
        return Point2::ORIGIN;
    }

    // 偏哪一侧**没有道理可讲** —— 两侧各有一个取代基,谁也不比谁更该。RDKit
    // 取 `otherNeighbor(begAt, endAt, 0)`,那是**存储序**:同一个分子换种写法
    // 就画到另一侧去,本库的头号契约不允许照抄。
    //
    // 这里按**画布的绝对方向**定:朝上那一侧,横着打平时朝左。坐标此时已经过
    // `orient::canonicalise`,与写法无关;而且这条规则只认方向、不认 `normal`
    // 的正负,所以**键的 begin/end 换个个儿它也不变**。
    if normal.y < -TIE_DIR {
        normal
    } else if normal.y > TIE_DIR {
        normal * -1.0
    } else if normal.x < 0.0 {
        normal
    } else {
        normal * -1.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::generate;

    fn prep(smi: &str) -> MolBuilder {
        let mut m = omgkit_io::smiles::parse(smi).unwrap();
        omgkit_chem::pipeline::sanitize(&mut m).unwrap();
        m
    }

    fn lines(s: &Scene) -> usize {
        s.items
            .iter()
            .filter(|i| matches!(i, Primitive::Line { .. }))
            .count()
    }
    fn texts(s: &Scene) -> usize {
        s.items
            .iter()
            .filter(|i| matches!(i, Primitive::Text { .. }))
            .count()
    }

    #[test]
    fn an_atom_whose_bonds_all_point_up_puts_its_hydrogen_straight_below() {
        // 对乙酰氨基酚的酰胺氮:一根键指左上、一根指右上,横向分量抵消掉了。
        // 只有左右两向可选时氢会挤到某根键下面,而**正下方整片是空的**。
        let m = prep("CC(=O)Nc1ccc(O)cc1");
        let n = (0..u32::try_from(m.num_atoms()).unwrap())
            .find(|a| m.atoms()[*a as usize].atomic_num == 7)
            .expect("对乙酰氨基酚有一个氮");
        for style in &Style::ALL {
            let d = generate(&m, style);
            // 先确认几何确实是"两根键都朝上"—— 否则这条判据验的不是要验的东西
            let here = d.coords[n as usize];
            for (nb, _) in m.neighbors(n) {
                assert!(
                    d.coords[nb as usize].y > here.y,
                    "[{}] 酰胺氮的邻居 {nb} 没在它上面,这个分子摆得和判据假设的不一样",
                    style.name
                );
            }
            assert_eq!(
                label_dir(&m, n, &d.coords),
                LabelDir::South,
                "[{}] 两根键都朝上,标签该往下伸",
                style.name
            );
        }
    }

    #[test]
    fn a_terminal_atom_never_stacks_its_hydrogen() {
        // 端基只有一根键。竖排既没必要(横着写就是 `–OH`、`–NH2` 这些常规写法),
        // 又会把一串端基排得高低不齐。RDKit 也是特判掉的。
        //
        // 光断言"没有端基竖排"是会空过的 —— 要是这批分子里根本没有键近乎竖直的
        // 端基,怎么改都是绿的。所以先数出**键确实近乎竖直**的端基有几个。
        let mut vertical_terminals = 0usize;
        for smi in [
            "CCO",
            "CC(=O)Oc1ccccc1C(=O)O",
            "NCCO",
            "OCC(O)C(O)C(O)C(O)C=O",
            "CN1C=NC2=C1C(=O)N(C)C(=O)N2C",
            "OC(=O)C(N)Cc1ccccc1",
            "NC(=O)c1ccccc1N",
            "OC(=O)c1ccccc1O",
        ] {
            let m = prep(smi);
            for style in &Style::ALL {
                let d = generate(&m, style);
                for a in 0..u32::try_from(m.num_atoms()).unwrap() {
                    if m.degree(a) != 1 {
                        continue;
                    }
                    let Some(l) = label_at(&m, a, style, &d.coords) else {
                        continue;
                    };
                    if !l.plain().contains('H') {
                        continue;
                    }
                    // **守卫不抄实现的判别式。** 照搬 `|sy| > VERT_SLOPE*|sx|`
                    // 会让守卫跟着实现一起动 —— 常量打错时守卫也跟着错,看不出来。
                    // 这里独立用角度表达同一个前提。
                    let (nb, _) = m.neighbors(a).next().expect("度 1 必有一个邻居");
                    let v = d.coords[nb as usize] - d.coords[a as usize];
                    let deg = v.y.atan2(v.x).to_degrees().abs();
                    if (70.0..110.0).contains(&deg) {
                        vertical_terminals += 1;
                    }
                    assert!(
                        !label_dir(&m, a, &d.coords).is_vertical(),
                        "[{}] {smi} 的端基 {a} 被判成了竖排",
                        style.name
                    );
                }
            }
        }
        assert!(
            vertical_terminals > 0,
            "这批分子里没有键近乎竖直的端基,判据是空过的"
        );
    }

    #[test]
    fn whether_a_label_stacks_does_not_depend_on_how_it_was_written() {
        // 判别量 `|sy| - tan70·|sx|` 是按 `mol.neighbors` 的存储序累加出来的,
        // 而三项以上的浮点加法不满足结合律 —— 度 ≥ 3 的原子上,同一个分子的两种
        // 写法会算出末位不同的判别量。落在断点上时,横竖就由那一位决定了。
        //
        // 取的是**语料里真出现过的**分子,不是自己编的。
        // 头两个是语料里**真的落在 70° 断点上**的分子(度 ≥ 3、判别量 ≤ 4.7e-15)
        // —— 这条判据要守的破口就在它们身上。后两个补一般情形。
        for smi in [
            "CC(=O)Nc1ccc(O)cc1",
            "CN1C=NC2=C1C(=O)N(C)C(=O)N2C",
            "NC(=O)c1ccccc1N",
            "OC(=O)C(N)Cc1ccccc1",
        ] {
            let m = prep(smi);
            let ranks = omgkit_io::canon::canonical_ranks(&m);
            // 规范秩必须是 0..n 的排列,才能拿它单独当排序键。**不留存储下标
            // 兜底** —— 那正是头号契约明令禁止的打破方式,躺在判据里也不行。
            let mut seen: Vec<u32> = ranks.clone();
            seen.sort_unstable();
            assert!(
                seen.iter().enumerate().all(|(i, r)| *r as usize == i),
                "{smi} 的规范秩不是全序,这条判据的前提不成立"
            );
            for style in &Style::ALL {
                let d = generate(&m, style);
                // **按 `scene` 真正画的那个分子算** —— 补过立体氢之后原子多了,
                // 拿原分子算会少算真的画出来的键。
                let g = d.drawn(&m);
                let base = dirs_by_rank(&g, &d.coords, &ranks);

                let mut writings: std::collections::BTreeSet<String> =
                    std::collections::BTreeSet::new();
                let mut compared = 0usize;
                for seed in 0..12u64 {
                    let w =
                        omgkit_io::smiles::write_with_priority(&m, &shuffled(m.num_atoms(), seed));
                    let Ok(mut m2) = omgkit_io::smiles::parse(&w.smiles) else {
                        continue;
                    };
                    if omgkit_chem::pipeline::sanitize(&mut m2).is_err() {
                        continue;
                    }
                    if omgkit_io::canon::canonical_smiles(&m2).smiles
                        != omgkit_io::canon::canonical_smiles(&m).smiles
                    {
                        continue;
                    }
                    let r2 = omgkit_io::canon::canonical_ranks(&m2);
                    let d2 = generate(&m2, style);
                    let g2 = d2.drawn(&m2);
                    let got = dirs_by_rank(&g2, &d2.coords, &r2);
                    compared += 1;
                    writings.insert(w.smiles.clone());
                    assert_eq!(
                        base, got,
                        "[{}] {smi} 写成 {} 之后标签朝向变了",
                        style.name, w.smiles
                    );
                }
                // **12 个种子必须都比上。** `compared > 0` 太松:写出器哪天回归成
                // "写出来读回去不是同一个分子",判据会静悄悄退化成只比一种写法,
                // 而且照样绿 —— 那三个 `continue` 吞的是真 bug。
                assert_eq!(compared, 12, "{smi} 只比上了 {compared} 种写法");
                // 种子数不等于写法数:搅拌出来的写法会重。真正不同的得够多,
                // 否则等于反复比同一种。
                assert!(
                    writings.len() >= 8,
                    "{smi} 12 个种子只产出 {} 种不同写法,判据太弱",
                    writings.len()
                );
            }
        }
    }

    /// 按规范秩排好的方向序列 —— 与原子编号无关的指纹。
    fn dirs_by_rank(m: &MolBuilder, coords: &[Point2], ranks: &[u32]) -> Vec<LabelDir> {
        let mut order: Vec<usize> = (0..m.num_atoms()).collect();
        order.sort_by_key(|i| ranks[*i]);
        order
            .iter()
            .map(|i| label_dir(m, u32::try_from(*i).unwrap(), coords))
            .collect()
    }

    #[test]
    fn a_seventy_degree_threshold_is_what_is_written_down() {
        // `VERT_SLOPE` 是写死的字面量(理由见它的文档:`tan` 不是正确舍入的,
        // 跨平台会漂)。写死就得有东西钉住它 —— 否则打错一位没人发现:实测
        // 换成 tan60 会让 3858 个原子改向、换成 tan80 会让 7587 个改向,而
        // 先前那三条判据在两种变异下**全是绿的**。
        assert!(
            (VERT_SLOPE.atan().to_degrees() - 70.0).abs() < 1e-12,
            "VERT_SLOPE 不是 tan 70°,反解出来是 {}°",
            VERT_SLOPE.atan().to_degrees()
        );
        // 光钉常量还不够 —— 得有真分子的行为把 60/70/80 分开。**两个原子都是
        // 从语料里搜出来的,不是编的**:第一个在 tan60 下会判成竖排,第二个在
        // tan80 下会判成横排,所以两条断言把三个阈值两两分开。
        for (smi, atom, want) in [
            (
                "[C@]1([C@H](C2CC1CC2)C(=O)OC)(C(=O)OC)C",
                2u32,
                LabelDir::East,
            ),
            ("[C@]12(C3CC([C@@H]1COC2=O)CC3)C", 0u32, LabelDir::North),
        ] {
            let m = prep(smi);
            let d = generate(&m, &Style::ACS_1996);
            assert_eq!(
                label_dir(&m, atom, &d.coords),
                want,
                "{smi} 的原子 {atom} 朝向不对 —— 阈值动了?"
            );
        }
    }

    #[test]
    fn a_lone_atom_writes_its_hydrogens_the_way_chemists_do() {
        // 度 0 没有"背离键"这回事。不特判的话通用分支一律给 `East`,
        // 水就画成 `OH2`、氯化氢画成 `ClH`。
        //
        // **全量审计查不出这条** —— `large.smi` 里度 0 原子一个都没有,
        // 而水合物、盐的对离子在真实输入里很常见。
        for (smi, want) in [
            ("O", "H2O"),
            ("Cl", "HCl"),
            ("Br", "HBr"),
            ("S", "H2S"),
            ("N", "NH3"),
            ("C", "CH4"),
        ] {
            let m = prep(smi);
            assert_eq!(m.num_atoms(), 1, "{smi} 不是单原子");
            let d = generate(&m, &Style::ACS_1996);
            let l = label_at(&m, 0, &Style::ACS_1996, &d.coords).expect("单原子必有标签");
            assert_eq!(l.plain(), want, "{smi} 画成了 {}", l.plain());
        }
    }

    /// splitmix64 + Fisher-Yates。**不用 `(i*k+k)%n` 这种仿射映射** —— 它保持
    /// 相邻关系,置换出来的写法太像原来的,判据会空过。
    fn shuffled(n: usize, seed: u64) -> Vec<u32> {
        let mut s = seed.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut next = || {
            s = s.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = s;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        };
        let mut v: Vec<u32> = (0..u32::try_from(n).unwrap()).collect();
        for i in (1..n).rev() {
            let j = usize::try_from(next() % (i as u64 + 1)).unwrap();
            v.swap(i, j);
        }
        v
    }

    #[test]
    fn a_collinear_atom_is_drawn_instead_of_vanishing() {
        // 两根键连成一条直线时顶点处没有拐角,**图上根本看不出那里有个原子**。
        // 丙二烯的中心碳是 sp、键角 180°,不画符号的话整张图是一条直线加几条
        // 平行短线 —— 读起来像顺式二烯,中心碳无影无踪。
        //
        // 这条同时守两件事:符号画出来了,而且两根双键各自**跨轴对称**。少一件
        // 都还是读错:只画符号不对称,两条同侧短线仍然像顺式;只对称不画符号,
        // 中心碳照样看不见。
        for (smi, centre) in [("CC=C=CC", 2u32), ("CC(C)=C=O", 3), ("C=C=C", 1)] {
            for style in &Style::ALL {
                let m = prep(smi);
                let d = generate(&m, style);
                let s = scene(&m, &d, style);
                // `coords`/`wedges` 的下标相对**被画的那个分子** —— 为画出构型补
                // 出来的氢也在里面。`scene` 自己会补,所以它拿的是原分子;这里影子化
                // 之后,下面按下标索引才不会错位甚至越界。
                let m = d.drawn(&m);
                let bnd = bounds(&d.coords, &m, style);
                let here = to_canvas(d.coords[centre as usize], bnd, style.bond_length_pt);

                assert!(
                    is_collinear(&m, centre, &d.coords),
                    "[{}] {smi}:原子 {centre} 该被判成共线",
                    style.name
                );
                let drawn = s.items.iter().any(|it| match it {
                    Primitive::Text { at, .. } => at.dist(here) < 0.5,
                    _ => false,
                });
                assert!(
                    drawn,
                    "[{}] {smi}:共线的原子 {centre} 一个符号都没画,图上看不见它",
                    style.name
                );

                // 中心那两根键的两条线必须对轴对称
                for (nb, bi) in m.neighbors(centre) {
                    if m.bonds()[bi as usize].order != BondOrder::Double {
                        continue;
                    }
                    let pa = to_canvas(d.coords[centre as usize], bnd, style.bond_length_pt);
                    let pb = to_canvas(d.coords[nb as usize], bnd, style.bond_length_pt);
                    let ls = lines_of_bond(&s, pa, pb);
                    assert_eq!(
                        ls.len(),
                        2,
                        "[{}] {smi}:键 {centre}–{nb} 该画两条线,实得 {}",
                        style.name,
                        ls.len()
                    );
                    let axis = (pb - pa).normalized();
                    let n = Point2::new(-axis.y, axis.x);
                    let side = |(u, v): (Point2, Point2)| ((u + v) * 0.5 - (pa + pb) * 0.5).dot(n);
                    let (s0, s1) = (side(ls[0]), side(ls[1]));
                    assert!(
                        s0 * s1 < 0.0 && (s0 + s1).abs() < 0.1 * s0.abs().max(1e-9),
                        "[{}] {smi}:键 {centre}–{nb} 的两条线没跨轴对称(偏移 {s0:.3} 与 {s1:.3})",
                        style.name
                    );
                }
            }
        }
    }

    #[test]
    fn a_bond_stops_at_the_label_box_not_at_its_circumscribed_circle() {
        // 键线在标签外停住是对的,**但停在哪儿有讲究**。先前按标签盒的外接圆
        // 切,圆一定包住盒,所以线绝不会压到字上 —— 代价是在盒窄的那个方向上
        // 停得太远:竖直方向去接一个横向宽的标签,白白空出 `hypot(w,h) − h`。
        //
        // 实测全量语料 129330 个带标签的键端,平均白切 0.075 个键长,最糟的
        // `[NH2+]` 白切 0.39 —— 快四成键长的空白。
        let style = &Style::ACS_1996;
        let m = prep("C=CC(=[NH2+])N");
        let l = label_for(&m, 3, style, HSide::Right).expect("[NH2+] 该有标签");

        // 竖直方向接近:盒边在 half_h,外接圆在 hypot(half_w, half_h)。
        // 盒心横向偏开了 `dx`,而竖直方向不受影响 —— 上下边仍在 ±half_h。
        let up = Point2::new(0.0, 1.0);
        let reach = box_reach(Point2::new(-l.dx, 0.0), up, l.half_w, l.half_h);
        assert!(
            (reach - l.half_h).abs() < 1e-9,
            "竖直方向该切到盒的上边 {},实得 {reach}",
            l.half_h
        );
        // **前提要自己成立才算数。** 这个标签上两种切法差得不够多的话,下面
        // 那条断言换回外接圆也照样绿 —— 判据就是空过的。
        let circle = l.half_w.hypot(l.half_h);
        assert!(
            circle - reach > 0.1,
            "这个标签上外接圆与盒边只差 {:.4} 个键长,拿它当判据说明不了问题",
            circle - reach
        );

        // 真正画出来的那一段必须按盒切。`trim` 两端各留一个 margin。
        let centre = Point2::new(0.0, 0.0);
        let other = Point2::new(0.0, 1.0);
        let (q, _) = trim(centre, other, &Some(l.clone()), &None, style);
        let cut = centre.dist(q);
        assert!(
            (cut - (reach + style.margin())).abs() < 1e-9,
            "竖直接近 [NH2+] 时该切 {},实得 {cut}",
            reach + style.margin()
        );
        assert!(
            cut < circle + style.margin() - 0.1,
            "还是按外接圆切的 —— 白空出 {:.4} 个键长",
            circle - reach
        );

        // 横向接近时两者本来就该一致(盒边正是 half_w),别把这条改坏
        let right = Point2::new(1.0, 0.0);
        assert!(
            (box_reach(Point2::new(-l.dx, 0.0), right, l.half_w, l.half_h) - (l.half_w + l.dx))
                .abs()
                < 1e-9,
            "横向该切到盒的右边"
        );
    }

    #[test]
    fn no_drawn_line_runs_across_an_atom_label() {
        // 切得更近是好事,**但不许近到压着字**。这条守的是另一半。两条判据必须
        // 同时在:只留上面那条的话,把 `trim` 改成"一律不切"能让它更好看,而字
        // 全被划掉。
        //
        // **一处例外要说清楚:两端标签加起来比一个键还长时,`trim` 走压缩兜底,
        // 端点确实会落进盒里。** 那不是切算错了,是 ACS 规范下标签本来就占
        // 0.69 个键长,`O⁻—N⁺` 两端要 1.375 个键长的净空,一个键长塞不下。
        // 换成按盒边切已经把这种键从 2.77% 压到 1.26%(全量 283604 根键),
        // 剩下的要靠逐字形的盒才能再降,不是这条判据管得了的。
        for smi in [
            "C=CC(=[NH2+])N",
            "OC(=O)c1ccccc1OC(C)=O",
            "CC(=O)Nc1ccc(O)cc1",
            "[O-][N+](=O)c1ccccc1S(=O)(=O)O",
            "NCCO",
            // 咖啡因:环内双键的端点带标签。**它在这条上并不吃劲** —— 内侧线
            // 压到 `N` 上时只差 0.14pt 就够不着这条判据留的线宽余量,那份
            // "越过主线伸过去"的难看由
            // [`a_double_bond_inner_line_does_not_reach_past_the_main_line_into_a_label`]
            // 那条量。留在这里是为了这一类形状也走一遍盒的规则。
            "CN1C=NC2=C1C(=O)N(C)C(=O)N2C",
            // 盒**宽**的标签:主线从上方出盒只让开半高,内侧线再横着挪 0.18 个
            // 键长仍在半宽之内 —— 不斜切也不够,还得沿线把端点推出盒去
            "c1cc2ccc[n+]3c2c(c1)SCC3",
        ] {
            for style in &Style::ALL {
                let m = prep(smi);
                let d = generate(&m, style);
                let s = scene(&m, &d, style);
                // `coords`/`wedges` 的下标相对**被画的那个分子** —— 为画出构型补
                // 出来的氢也在里面。`scene` 自己会补,所以它拿的是原分子;这里影子化
                // 之后,下面按下标索引才不会错位甚至越界。
                let m = d.drawn(&m);
                let bnd = bounds(&d.coords, &m, style);
                let scale = style.bond_length_pt;
                let labels: Vec<Option<Label>> = (0..u32::try_from(m.num_atoms()).unwrap())
                    .map(|a| label_at(&m, a, style, &d.coords))
                    .collect();
                // 挨着"塞不下"的键的原子整个跳过 —— 见上面那段。
                // **口径要问实现要,不许在这儿再写一遍**,见 [`is_squeezed`]。
                let squeezed: Vec<bool> = {
                    let mut v = vec![false; m.num_atoms()];
                    for b in m.bonds() {
                        if is_squeezed(
                            d.coords[b.begin as usize],
                            d.coords[b.end as usize],
                            labels[b.begin as usize].as_ref(),
                            labels[b.end as usize].as_ref(),
                            style,
                        ) {
                            v[b.begin as usize] = true;
                            v[b.end as usize] = true;
                        }
                    }
                    v
                };
                for a in 0..u32::try_from(m.num_atoms()).unwrap() {
                    let Some(l) = &labels[a as usize] else {
                        continue;
                    };
                    if squeezed[a as usize] {
                        continue;
                    }
                    let c = to_canvas(d.coords[a as usize], bnd, scale);
                    // 盒是画布坐标系里的轴对齐矩形。**留一点余量** —— 线本身有
                    // 粗细,压边一丝不算划字。
                    let (hw, hh) = (
                        l.half_w * scale - style.line_width_pt,
                        l.half_h * scale - style.line_width_pt,
                    );
                    if hw <= 0.0 || hh <= 0.0 {
                        continue;
                    }
                    for it in &s.items {
                        let Primitive::Line { from, to, .. } = it else {
                            continue;
                        };
                        // **盒心不在原子上** —— 整串朝一侧挪了 `dx`
                        let bc = c + Point2::new(l.dx * scale, 0.0);
                        let inside =
                            |p: &Point2| (p.x - bc.x).abs() < hw && (p.y - bc.y).abs() < hh;
                        assert!(
                            !inside(from) && !inside(to),
                            "[{}] {smi}:原子 {a} 的标签 {} 被一条线的端点压在里面",
                            style.name,
                            l.plain()
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn benzene_is_drawn_with_three_double_bonds_not_six_plain_lines() {
        // 净化之后芳香键的 order 是 Aromatic。直接按它画,苯环就成了一个
        // 六边形圈 —— 拓扑没错、看着也像个环,但少了三条线。
        let m = prep("c1ccccc1");
        let d = generate(&m, &Style::ACS_1996);
        let s = scene(&m, &d, &Style::ACS_1996);
        assert_eq!(lines(&s), 9, "六根键里三根是双键 → 6 + 3 = 9 条线");
        assert_eq!(texts(&s), 0, "苯环全是骨架碳,不该有标签");
    }

    #[test]
    fn a_triple_bond_gets_three_lines() {
        let m = prep("CC#N");
        let d = generate(&m, &Style::ACS_1996);
        let s = scene(&m, &d, &Style::ACS_1996);
        assert_eq!(lines(&s), 4, "一根单键 + 一根三键(3 条线)");
    }

    #[test]
    fn the_y_axis_is_flipped_for_the_canvas() {
        // 化学图 y 向上,画布 y 向下。不翻的话整张图上下镜像 —— 键长拓扑全对,
        // 只有手性看着反了,那是最难发现的一类错。
        let m = prep("CCO");
        let d = generate(&m, &Style::ACS_1996);
        let s = scene(&m, &d, &Style::ACS_1996);
        // 原子坐标里 y 最大的那个,在画布上应当 y 最小(最靠上)
        let hi = (0..d.coords.len())
            .max_by(|a, b| d.coords[*a].y.partial_cmp(&d.coords[*b].y).unwrap())
            .unwrap();
        let lo = (0..d.coords.len())
            .min_by(|a, b| d.coords[*a].y.partial_cmp(&d.coords[*b].y).unwrap())
            .unwrap();
        assert!(d.coords[hi].y > d.coords[lo].y, "测试样本的 y 应当有高低差");

        let scale = Style::ACS_1996.bond_length_pt;
        let (_, _, _, max_y) = bounds(&d.coords, &m, &Style::ACS_1996);
        let y_of = |p: Point2| (max_y - p.y) * scale + PAD_PT;
        assert!(
            y_of(d.coords[hi]) < y_of(d.coords[lo]),
            "y 轴没翻:化学上更高的原子在画布上却更靠下"
        );
        assert!(s.height > 0.0 && s.width > 0.0);
    }

    #[test]
    fn bonds_stop_short_of_atom_labels() {
        // 线画进标签里会把字盖掉。停多远由规范的 margin 定。
        let m = prep("CCO");
        let d = generate(&m, &Style::ACS_1996);
        let style = Style::ACS_1996;
        let o = 2u32; // 氧,有标签 "OH"
        let l = label_for(&m, o, &style, HSide::Right).expect("氧应当有标签");
        let neighbour = m.neighbors(o).next().expect("氧连着一个碳").0;
        let (_, trimmed) = trim(
            d.coords[neighbour as usize],
            d.coords[o as usize],
            &None,
            &Some(l.clone()),
            &style,
        );
        // **不能拿"离原子多远"当判据。** 盒心朝一侧挪了 `dx`(好让元素符号落在
        // 原子位置上),所以离原子远近说明不了有没有盖住字 —— 要量的是**离盒心
        // 多远**,而且要沿着键的方向比盒在那个方向上的半径。
        let centre = d.coords[o as usize] + Point2::new(l.dx, 0.0);
        let dir = (d.coords[o as usize] - d.coords[neighbour as usize]).normalized();
        let from_centre = trimmed - centre;
        // 端点在盒外:两个轴上至少有一个超出
        let outside =
            from_centre.x.abs() >= l.half_w - 1e-9 || from_centre.y.abs() >= l.half_h - 1e-9;
        assert!(
            outside,
            "键停得太近,端点落进标签盒里了:端点离盒心 ({:.4}, {:.4}),盒半宽高 ({:.4}, {:.4})",
            from_centre.x, from_centre.y, l.half_w, l.half_h
        );
        // 而且要**恰好**停在盒边外一个 margin —— 停太远就是先前那种大片空白
        let reach = box_reach(Point2::new(-l.dx, 0.0), dir * -1.0, l.half_w, l.half_h);
        let gap = trimmed.dist(d.coords[o as usize]);
        assert!(
            (gap - (reach + style.margin())).abs() < 1e-9,
            "该停在盒边外一个 margin({}),实得 {gap}",
            reach + style.margin()
        );
    }

    #[test]
    fn the_two_styles_give_different_canvas_sizes() {
        // **规范贯穿到渲染的落点。** ChemDraw 默认的键长是 ACS 的 2.08 倍,
        // 同一个分子的画布也该差不多是那个比例。两边一样大就说明 Style
        // 没有真的传到渲染这一层。
        let m = prep("c1ccc2ccccc2c1");
        let a = scene(&m, &generate(&m, &Style::ACS_1996), &Style::ACS_1996);
        let c = scene(
            &m,
            &generate(&m, &Style::CHEMDRAW_DEFAULT),
            &Style::CHEMDRAW_DEFAULT,
        );
        assert!(
            c.width > a.width * 1.8,
            "画布宽度比只有 {:.2}",
            c.width / a.width
        );
        assert_eq!(lines(&a), lines(&c), "图元数量不该随规范变");
    }

    /// 找出画布上属于某根键的线段:两个端点都贴着这根键的轴线、且落在它的跨度内。
    ///
    /// 相邻的键成 120°(或更小的角),远端离轴线 0.87 个键长开外;环上平行的
    /// 对边至少隔一个键长。都落不进 0.30 这个窗口,所以选出来的只会是这根键
    /// 自己画的那一条或两条线。
    fn lines_of_bond(s: &Scene, pa: Point2, pb: Point2) -> Vec<(Point2, Point2)> {
        let len = pa.dist(pb);
        let axis = (pb - pa) * (1.0 / len);
        let normal = Point2::new(-axis.y, axis.x);
        let mut out = Vec::new();
        for it in &s.items {
            if let Primitive::Line { from, to, .. } = it {
                let inside = |p: Point2| {
                    let v = p - pa;
                    v.dot(normal).abs() < 0.30 * len
                        && v.dot(axis) > -0.10 * len
                        && v.dot(axis) < 1.10 * len
                };
                if inside(*from) && inside(*to) {
                    out.push((*from, *to));
                }
            }
        }
        out
    }

    #[test]
    fn a_fused_double_bond_goes_into_the_aromatic_ring() {
        // 稠合处共用的双键,第二条线该画进**芳香**的那个环 —— 那一圈的交替
        // 单双键要靠它接上。
        //
        // 先前的规则是"双键最多的那个环"。它是芳香性的**代理**:kekulé 化之后
        // 芳环的交替双键通常就是最多的,所以两条规则几乎总是一致 —— 全量语料
        // 1092 根稠合处的共用双键里只有 **3** 根分歧。代理会在"非芳环的双键比
        // 隔壁芳环还多"时给出错的答案。
        //
        // 下面这三个就是量出来的那三根。判据不是我挑的分子,是扫出来的。
        for smi in [
            "[H]/N=c/1\\[nH]c-2c(s1)CSc3c2cccc3",
            "c1cc2ccc3c(c2nc1)N=C[C@H](C3=O)C(=O)[O-]",
            "COCc1[nH]c2c(n1)-c3c(nc(o3)N)[C@H]2c4ccc(cc4)Cl",
        ] {
            for style in &Style::ALL {
                let m = prep(smi);
                let d = generate(&m, style);
                let s = scene(&m, &d, style);
                // `coords`/`wedges` 的下标相对**被画的那个分子** —— 为画出构型补
                // 出来的氢也在里面。`scene` 自己会补,所以它拿的是原分子;这里影子化
                // 之后,下面按下标索引才不会错位甚至越界。
                let m = d.drawn(&m);
                let bnd = bounds(&d.coords, &m, style);
                let pts: Vec<Point2> = d
                    .coords
                    .iter()
                    .map(|p| to_canvas(*p, bnd, style.bond_length_pt))
                    .collect();
                let rings = omgkit_chem::sssr::ring_set(&m);
                let orders = drawn_orders(&m);
                let mut checked = 0usize;
                for (bi, b) in m.bonds().iter().enumerate() {
                    if orders[bi] != BondOrder::Double {
                        continue;
                    }
                    let bn = u32::try_from(bi).unwrap();
                    let mine: Vec<&omgkit_chem::sssr::Ring> =
                        rings.iter().filter(|r| r.bonds.contains(&bn)).collect();
                    if mine.len() < 2 {
                        continue;
                    }
                    let arom = |r: &omgkit_chem::sssr::Ring| {
                        r.bonds
                            .iter()
                            .all(|x| m.bonds()[*x as usize].flags.contains(BondFlags::AROMATIC))
                    };
                    // 只查"一芳一不芳"的共用键 —— 都芳香时这条判据说不出高下
                    let (aro, non): (
                        Vec<&&omgkit_chem::sssr::Ring>,
                        Vec<&&omgkit_chem::sssr::Ring>,
                    ) = mine.iter().partition(|r| arom(r));
                    if aro.len() != 1 || non.is_empty() {
                        continue;
                    }
                    let (pa, pb) = (pts[b.begin as usize], pts[b.end as usize]);
                    let ls = lines_of_bond(&s, pa, pb);
                    let Some(inner) = ls
                        .iter()
                        .find(|(u, v)| u.dist(pa) > 1e-6 || v.dist(pb) > 1e-6)
                        .copied()
                    else {
                        continue;
                    };
                    let mid = (inner.0 + inner.1) * 0.5;
                    let poly: Vec<Point2> = aro[0].atoms.iter().map(|a| pts[*a as usize]).collect();
                    assert!(
                        crate::geom::point_in_polygon(mid, &poly),
                        "[{}] {smi}:键 {}–{} 的第二条线没画进芳环里",
                        style.name,
                        b.begin,
                        b.end
                    );
                    checked += 1;
                }
                assert!(
                    checked > 0,
                    "[{}] {smi}:一根「一芳一不芳」的共用双键都没查到 —— 判据空过了",
                    style.name
                );
            }
        }
    }

    #[test]
    fn a_ring_double_bond_inner_line_meets_the_neighbouring_bonds() {
        // 环内双键的第二条线两端要**斜切**到与相邻键接上,而不是按一个拍脑袋的
        // 定值往里缩。
        //
        // 苯环是规则六边形,这件事有闭式解:内侧线平行于环边、向内偏 `s`,它与
        // 两端顶点的**角平分线**(规则多边形里就是"环心到顶点"那条射线)的交点,
        // 正是一个相似的小六边形的顶点 —— 到环心的距离按 `(边心距 − s)/边心距`
        // 等比缩小。端点落不到那条射线上,画出来就是内侧线比相邻键短一截或者
        // 戳出去一点,接头处露白。
        //
        // 先前是两端各缩固定的 12%。六边形上正确的缩进量是
        // `s / (2 × 边心距) = 0.18 / (2 × 0.866) = 10.4%`,差 1.6 个百分点。
        for style in &Style::ALL {
            let m = prep("c1ccccc1");
            let d = generate(&m, style);
            let s = scene(&m, &d, style);
            // `coords`/`wedges` 的下标相对**被画的那个分子** —— 为画出构型补
            // 出来的氢也在里面。`scene` 自己会补,所以它拿的是原分子;这里影子化
            // 之后,下面按下标索引才不会错位甚至越界。
            let m = d.drawn(&m);
            let bnd = bounds(&d.coords, &m, style);
            let pts: Vec<Point2> = d
                .coords
                .iter()
                .map(|p| to_canvas(*p, bnd, style.bond_length_pt))
                .collect();
            let centre =
                pts.iter().fold(Point2::ORIGIN, |acc, p| acc + *p) * (1.0 / pts.len() as f64);

            let mut checked = 0usize;
            for (bi, b) in m.bonds().iter().enumerate() {
                if drawn_orders(&m)[bi] != BondOrder::Double {
                    continue;
                }
                let (pa, pb) = (pts[b.begin as usize], pts[b.end as usize]);
                let ls = lines_of_bond(&s, pa, pb);
                assert_eq!(ls.len(), 2, "环内双键该画两条线");
                // 主线是与 pa/pb 重合的那条,另一条是内侧线
                let inner = ls
                    .iter()
                    .find(|(u, v)| u.dist(pa) > 1e-6 || v.dist(pb) > 1e-6)
                    .copied()
                    .expect("该有一条内侧线");
                for (end, vertex) in [(inner.0, pa), (inner.1, pb)] {
                    // 端点必须落在"环心 → 顶点"这条射线上
                    let ray = (vertex - centre).normalized();
                    let off = end - centre;
                    let perp = off.x * ray.y - off.y * ray.x;
                    assert!(
                        perp.abs() < 1e-6,
                        "[{}] 苯环内侧线的端点偏离顶点角平分线 {perp:.4} pt —— 接头处对不上",
                        style.name
                    );
                }
                checked += 1;
            }
            assert_eq!(checked, 3, "苯环该有三根凯库勒双键");
        }
    }

    #[test]
    fn a_double_bond_inner_line_does_not_reach_past_the_main_line_into_a_label() {
        // 咖啡因咪唑环上的 C=N。**先前内侧线比主线还伸得靠前**:主线为了让开
        // `N` 缩到离原子中心 5.56pt,内侧线的端点却在 3.20pt —— 那是
        // [`mitre_end`] 从**原子中心**量的角平分线交点,它压根不知道主线已经
        // 缩过了。画出来是一个指着 N 的楔子,而 `N` 的半高只有 3.59pt,内侧线
        // 直接压在字上。
        //
        // 契约:端原子有标签时,内侧线**沿键轴**的落点不许越过主线的落点。
        // 只量沿轴的投影 —— 横向那 0.18 个键长是双键本来就该有的间距。
        //
        // 咖啡因不是随便挑的:此前那条 `no_drawn_line_runs_across_an_atom_label`
        // 列的五个分子里,没有一个的**环内**双键端点带标签,于是这一整类漏了
        // 出去。全量语料上这条错了 3044 处(17.2%)。
        for style in &Style::ALL {
            let m = prep("CN1C=NC2=C1C(=O)N(C)C(=O)N2C");
            let d = generate(&m, style);
            let s = scene(&m, &d, style);
            // `coords`/`wedges` 的下标相对**被画的那个分子** —— 为画出构型补
            // 出来的氢也在里面。`scene` 自己会补,所以它拿的是原分子;这里影子化
            // 之后,下面按下标索引才不会错位甚至越界。
            let m = d.drawn(&m);
            let bnd = bounds(&d.coords, &m, style);
            let scale = style.bond_length_pt;
            let pts: Vec<Point2> = d.coords.iter().map(|p| to_canvas(*p, bnd, scale)).collect();

            let mut checked = 0usize;
            for (bi, b) in m.bonds().iter().enumerate() {
                if drawn_orders(&m)[bi] != BondOrder::Double {
                    continue;
                }
                let (pa, pb) = (pts[b.begin as usize], pts[b.end as usize]);
                let len = pa.dist(pb);
                let axis = (pb - pa) * (1.0 / len);
                let normal = Point2::new(-axis.y, axis.x);
                let ls = lines_of_bond(&s, pa, pb);
                assert_eq!(ls.len(), 2, "[{}] 双键该画两条线", style.name);
                // 主线贴着键轴,内侧线离轴一个间距。都离轴的是对称画法
                // (端基双键),那时没有"主线"可比,跳过。
                let off = |(u, v): &(Point2, Point2)| {
                    (((*u - pa).dot(normal)).abs() + ((*v - pa).dot(normal)).abs()) / 2.0
                };
                let (main, inner) = if off(&ls[0]) < off(&ls[1]) {
                    (ls[0], ls[1])
                } else {
                    (ls[1], ls[0])
                };
                if off(&main) > 0.05 * len {
                    continue;
                }
                for (a, t_of) in [
                    (b.begin, 1.0_f64), // 在 pa 这头,投影越小越靠外
                    (b.end, -1.0),
                ] {
                    if label_at(&m, a, style, &d.coords).is_none() {
                        continue;
                    }
                    // 把两条线各自离这个原子最近的那个端点投影到键轴上
                    let near = |(u, v): &(Point2, Point2)| {
                        let (tu, tv) = ((*u - pa).dot(axis), (*v - pa).dot(axis));
                        if t_of > 0.0 {
                            tu.min(tv)
                        } else {
                            tu.max(tv)
                        }
                    };
                    let (tm, ti) = (near(&main), near(&inner));
                    assert!(
                        (ti - tm) * t_of >= -0.01,
                        "[{}] 咖啡因键 {bi} 的内侧线越过主线伸向带标签的原子 {a}:\
                         主线停在 {tm:.2},内侧线停在 {ti:.2}(键长 {len:.2})",
                        style.name
                    );
                    checked += 1;
                }
            }
            assert!(
                checked > 0,
                "[{}] 咖啡因一处都没查到 —— 判据空过了",
                style.name
            );
        }
    }

    #[test]
    fn a_ring_double_bond_never_puts_a_line_outside_the_ring() {
        // **这条查画出来的线,不查 offset_dir 的返回值。**
        //
        // 先前那条只查后者,于是把"方向算对了、落点却错了"整类缺陷全放过 ——
        // 方向是在布局坐标系(y 向上)里算的,落点用在画布坐标系(y 向下)里,
        // 两边一混,得到的符号是 `vx² − vy²` 的符号:竖着的键碰巧对,横着的
        // 正好反。实测阿司匹林苯环底边那根双键的第二条线画到了环外。
        //
        // 顺带守住另一件事:环上双键偏哪侧**与取代基多少无关**。抗坏血酸环内
        // 那根 C=C 两端各挂一个 OH,按邻居计数正好抵消,两条线就骑在环边上,
        // 其中一条落在环外。
        for smi in [
            "CC(=O)Oc1ccccc1C(=O)O",                     // 阿司匹林:苯环底边是横的
            "OC1=C(O)C(=O)OC1[C@@H](O)CO",               // 抗坏血酸:环内 C=C 两端都有 OH
            "CC1=C(C)CCCC1",                             // 1,2-二甲基环己烯:同上,更干净
            "c1ccccc1",                                  // 苯
            "c1ccc2ccccc2c1",                            // 萘:有共用键
            "c1ccncc1",                                  // 吡啶:环上有标签,键会被截短
            "CN1C=NC2=C1C(=O)N(C)C(=O)N2C",              // 咖啡因
            "C1=CC2=CC=CC3=CC=CC(=C1)C23",               // 苊:环大小不一
            "C1=CC2CCC1CC2", // 双环[2.2.2]辛烯:共用键的两个环在**同一侧**
            "C[C@@]12C(=C[C@@H](O1)C(=O)C23CC3)C(=O)OC", // 螺环 + 凹环:形心判反的那类
        ] {
            for style in &Style::ALL {
                let m = prep(smi);
                let d = generate(&m, style);
                let s = scene(&m, &d, style);
                // `coords`/`wedges` 的下标相对**被画的那个分子** —— 为画出构型补
                // 出来的氢也在里面。`scene` 自己会补,所以它拿的是原分子;这里影子化
                // 之后,下面按下标索引才不会错位甚至越界。
                let m = d.drawn(&m);
                let bnd = bounds(&d.coords, &m, style);
                let pts: Vec<Point2> = d
                    .coords
                    .iter()
                    .map(|p| to_canvas(*p, bnd, style.bond_length_pt))
                    .collect();
                let rings = omgkit_chem::sssr::ring_set(&m);
                // 键级要跟实现取同一个源:实现是在**规范重排后的副本**上凯库勒化的,
                // 判据自己再凯库勒化一遍会挑到另一套单双键,查的就不是同一根键了
                let orders = drawn_orders(&m);

                let mut checked = 0;
                for (bi, b) in m.bonds().iter().enumerate() {
                    let mine: Vec<_> = rings
                        .iter()
                        .filter(|r| r.bonds.contains(&u32::try_from(bi).unwrap()))
                        .collect();
                    if orders[bi] != BondOrder::Double || mine.is_empty() {
                        continue;
                    }
                    let (pa, pb) = (pts[b.begin as usize], pts[b.end as usize]);
                    let len = pa.dist(pb);
                    let mid = (pa + pb) * 0.5;
                    let axis = (pb - pa) * (1.0 / len);
                    let normal = Point2::new(-axis.y, axis.x);

                    let got = lines_of_bond(&s, pa, pb);
                    assert_eq!(
                        got.len(),
                        2,
                        "[{}] {smi}:双键 {bi} 画出了 {} 条线",
                        style.name,
                        got.len()
                    );
                    for (f, t) in got {
                        let side = ((f + t) * 0.5 - mid).dot(normal);
                        // 骑在轴上的那条不偏不倚,放过
                        if side.abs() < 0.02 * len {
                            continue;
                        }
                        // 共用键两侧各是一个环,进哪个都算数;只属于一个环的
                        // 键就只有一个"内"。
                        //
                        // 判"在不在环里"用射线法,不用"和环心同侧" —— 凹环上
                        // 两者不一致,而退化布局给出的环就是凹的、甚至自交的。
                        let ok = mine.iter().any(|r| {
                            let poly: Vec<Point2> =
                                r.atoms.iter().map(|a| pts[*a as usize]).collect();
                            crate::geom::point_in_polygon((f + t) * 0.5, &poly)
                        });
                        assert!(
                            ok,
                            "[{}] {smi}:双键 {bi}({}–{})有一条线画在了环外",
                            style.name, b.begin, b.end
                        );
                    }
                    checked += 1;
                }
                assert!(checked >= 1, "[{}] {smi}:一根环内双键都没查到", style.name);
            }
        }
    }

    #[test]
    fn the_same_molecule_written_differently_draws_the_same_lines() {
        // **坐标相同还不够。** 先前那条判据(lib.rs)比的是坐标,而坐标确实
        // 逐点相同 —— 变的是**键级**:凯库勒化要在几套等价的单双键里挑一套,
        // 挑哪一套取决于原子的存储顺序。同一个分子换个写法,苯环上三根双键的
        // 位置整个换一圈,画出来是另一张图。实测阿司匹林、萘、吡啶都会变。
        let groups = [
            vec!["CC(=O)Oc1ccccc1C(=O)O", "c1cc(OC(C)=O)c(cc1)C(O)=O"],
            vec!["c1ccc2ccccc2c1", "c1cc2c(cc1)cccc2", "c1cc2ccccc2cc1"],
            vec!["c1ccncc1", "c1cccnc1"],
            vec!["Cn1cnc2c1c(=O)n(C)c(=O)n2C", "CN1C=NC2=C1C(=O)N(C)C(=O)N2C"],
            vec!["OC(=O)c1ccccc1", "c1ccccc1C(O)=O"],
            // 反式双键:两个取代基一边一个,第二条线偏哪侧**没有道理可讲**。
            // RDKit 按存储序挑,那正是这条判据要挡的东西 —— 换个写法就画到另
            // 一侧去。我们按画布的绝对方向挑,与写法无关。
            //
            // **必须是反式-2-戊烯,不能用 2-丁烯。** 两条都踩过:
            //
            // - 先前写的 `["C/C=C/C", "C(=C\\C)\\C", "C(/C)=C/C"]` **根本不是
            //   同一个分子** —— 后两个的规范式是 `C/C=C\C`(顺式)。它能绿只
            //   因为这里的 `prep` 不调 `perceive_bond_stereo`,布局看不见 E/Z,
            //   顺反画成了同一张图。**那等于把"顺式和反式必须画得一模一样"
            //   写成了硬性要求**,哪天 prep 对齐了它就红,而红的是判据自己。
            // - 换成三种真反式的 2-丁烯又**抓不住存储序**:它两端等价,
            //   `otherNeighbor(begAt, endAt, 0)` 挑哪一端都落到同一侧。
            //
            // 反式-2-戊烯两端不等价,四种写法的规范式都是 `C/C=C/CC`,而且在
            // 感不感知顺反两套口径下都落进"票数抵消"那一档。
            vec!["C/C=C/CC", "CC/C=C/C", "C(/C)=C\\CC", "C(\\CC)=C/C"],
        ];
        // 图元的多重集:坐标量化后排序。原子编号不同,画出来的线必须一样。
        let key = |smi: &str, style: &Style| -> Vec<String> {
            let m = prep(smi);
            let s = scene(&m, &generate(&m, style), style);
            let q = |p: Point2| format!("{:.3},{:.3}", p.x, p.y);
            let mut v: Vec<String> = s
                .items
                .iter()
                .map(|it| match it {
                    // 线段**不分方向**:`L A B` 与 `L B A` 是同一条线,而谁是
                    // 起点取决于键的 begin/end,那本来就随写法变
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
        };
        for style in &Style::ALL {
            for ws in &groups {
                let a = key(ws[0], style);
                assert!(!a.is_empty(), "{} 什么都没画出来", ws[0]);
                for w in &ws[1..] {
                    assert_eq!(
                        a,
                        key(w, style),
                        "[{}] {w} 与 {} 画出来不一样",
                        style.name,
                        ws[0]
                    );
                }
            }
        }
    }

    #[test]
    fn every_recorded_wedge_actually_reaches_the_canvas() {
        // `Depiction` 里记了楔形,`Scene` 里却没有对应的图元 —— 图上那个立体
        // 中心就是白的,而 `unwedged` 是空的、诊断全绿。
        //
        // 成因是指派可以挑中一根**双键**:双键在这里走的是另一条分支,那条
        // 分支根本不看楔形。四配位的 P(V)、S(VI) 中心会碰到。
        for smi in [
            "CCO[P@@]1(=O)CCCCN1Cc2ccccc2", // 磷(V):中心的一根键是 P=O
            "CCO[P@@]1(=O)CCC[C@@H](C1)C",
            "OC[C@H](O)[C@H]1OC(=O)C(O)=C1O",
            "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O",
        ] {
            for style in &Style::ALL {
                let m = prep(smi);
                let d = generate(&m, style);
                let s = scene(&m, &d, style);
                let recorded = d.wedges.iter().filter(|w| w.narrow().is_some()).count();
                let drawn = s
                    .items
                    .iter()
                    .filter(|it| matches!(it, Primitive::Wedge { .. } | Primitive::Hash { .. }))
                    .count();
                assert_eq!(
                    recorded, drawn,
                    "[{}] {smi}:记了 {recorded} 个楔形,画出来 {drawn} 个",
                    style.name
                );
            }
        }
    }

    #[test]
    fn a_wedge_starts_at_the_stereocentre_it_describes() {
        // 楔形的窄端标出的是"这个构型说的是谁"。**两个立体中心相邻时它们共用
        // 一根键,两头都带手性标记** —— 靠"哪头带标记"去猜窄端就会猜到前一头,
        // 于是这根键描述的构型整个反过来,而图上一点异常都看不出来。
        //
        // 实测:抗坏血酸的两个立体中心正好相邻。
        for smi in [
            "OC[C@H](O)[C@H]1OC(=O)C(O)=C1O",              // 抗坏血酸
            "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O", // 葡萄糖:5 个连成一串
            "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O", // 青霉素 G
        ] {
            for style in &Style::ALL {
                let m = prep(smi);
                let d = generate(&m, style);
                let s = scene(&m, &d, style);
                // `coords`/`wedges` 的下标相对**被画的那个分子** —— 为画出构型补
                // 出来的氢也在里面。`scene` 自己会补,所以它拿的是原分子;这里影子化
                // 之后,下面按下标索引才不会错位甚至越界。
                let m = d.drawn(&m);
                let bnd = bounds(&d.coords, &m, style);
                let pts: Vec<Point2> = d
                    .coords
                    .iter()
                    .map(|p| to_canvas(*p, bnd, style.bond_length_pt))
                    .collect();

                let mut checked = 0;
                for it in &s.items {
                    let (f, t) = match it {
                        Primitive::Wedge { from, to, .. } | Primitive::Hash { from, to, .. } => {
                            (*from, *to)
                        }
                        _ => continue,
                    };
                    // 图元反查是哪根键:楔形的两端**正落在键轴上**(截短只沿轴
                    // 缩),所以要求两个端点离这根键的轴线都不到 2%、且落在跨度
                    // 内。先前用的是"两端配对最近的那根",在拥挤的布局上会配到
                    // 隔壁键去 —— 那样断言比的就是另一根键,可能瞎蒙对。
                    let on = |b: &omgkit_core::BondData| {
                        let (pu, pv) = (pts[b.begin as usize], pts[b.end as usize]);
                        let len = pu.dist(pv);
                        if len < 1e-9 {
                            return false;
                        }
                        let axis = (pv - pu) * (1.0 / len);
                        let normal = Point2::new(-axis.y, axis.x);
                        [f, t].iter().all(|p| {
                            let v = *p - pu;
                            v.dot(normal).abs() < 0.02 * len
                                && v.dot(axis) > -0.05 * len
                                && v.dot(axis) < 1.05 * len
                        })
                    };
                    let hits: Vec<usize> = m
                        .bonds()
                        .iter()
                        .enumerate()
                        .filter(|(_, b)| on(b))
                        .map(|(i, _)| i)
                        .collect();
                    assert_eq!(
                        hits.len(),
                        1,
                        "[{}] {smi}:一个楔形图元落在了 {} 根键上,反查不出唯一的键",
                        style.name,
                        hits.len()
                    );
                    let bi = hits[0];
                    let narrow = d.wedges[bi].narrow().expect("画出了楔形的键必定记着窄端");
                    let b = &m.bonds()[bi];
                    let other = if narrow == b.begin { b.end } else { b.begin };
                    assert!(
                        f.dist(pts[narrow as usize]) < f.dist(pts[other as usize]),
                        "[{}] {smi}:键 {bi}({}–{})的楔形窄端画在了原子 {other} 那头,\
                         而它描述的是原子 {narrow} 的构型",
                        style.name,
                        b.begin,
                        b.end
                    );
                    checked += 1;
                }
                assert!(
                    checked >= 2,
                    "[{}] {smi}:只查到 {checked} 个楔形",
                    style.name
                );
            }
        }
    }

    #[test]
    fn the_recorded_narrow_end_is_what_gets_drawn() {
        // 上面那条走的是完整流水线,而流水线**会尽量避开**两头都是立体中心的键
        // —— 于是"按记录取窄端"和"按哪头带手性标记猜"给出同一个答案,那条判据
        // 分不出这两种实现。
        //
        // 这一条直接把窄端翻到另一头再画,画出来的窄端必须跟着走。渲染若不看
        // 记录、自己去猜,两张图会一模一样。
        let m = prep("N[C@@H](C)O");
        let style = Style::ACS_1996;
        let mut d = generate(&m, &style);
        let bi = d
            .wedges
            .iter()
            .position(|w| w.narrow().is_some())
            .expect("这个分子该有一个楔形");
        let b = &m.bonds()[bi];
        let bnd = bounds(&d.coords, &m, &style);
        let pts: Vec<Point2> = d
            .coords
            .iter()
            .map(|p| to_canvas(*p, bnd, style.bond_length_pt))
            .collect();

        let narrow_of = |s: &Scene| {
            s.items
                .iter()
                .find_map(|it| match it {
                    Primitive::Wedge { from, .. } | Primitive::Hash { from, .. } => Some(*from),
                    _ => None,
                })
                .expect("场景里该有一个楔形图元")
        };

        let before = narrow_of(&scene(&m, &d, &style));
        let (was, other) = if d.wedges[bi].narrow() == Some(b.begin) {
            (b.begin, b.end)
        } else {
            (b.end, b.begin)
        };
        assert!(
            before.dist(pts[was as usize]) < before.dist(pts[other as usize]),
            "改之前窄端就没画在记录的那一头"
        );

        // 只翻窄端,楔形的虚实不动
        d.wedges[bi] = match d.wedges[bi] {
            crate::stereo::Wedge::Up { .. } => crate::stereo::Wedge::Up { narrow: other },
            crate::stereo::Wedge::Down { .. } => crate::stereo::Wedge::Down { narrow: other },
            crate::stereo::Wedge::None => unreachable!("上面刚确认它是个楔形"),
        };
        let after = narrow_of(&scene(&m, &d, &style));
        assert!(
            after.dist(pts[other as usize]) < after.dist(pts[was as usize]),
            "窄端记到原子 {other} 上了,画出来的窄端却还在原子 {was} 那头 —— \
             渲染没看记录,在自己猜"
        );
    }

    #[test]
    fn a_fused_ring_double_bond_goes_inside_one_ring_not_astride_the_bond() {
        // 稠合处共用的那根键两头各是一个环。先前一律对称画,于是两条线各贴着
        // 环边、哪个环都不属于 —— 看着像两根平行的键,而不是一根双键。
        //
        // 通例是画进"更芳香"的那个环。这条查的是**两条线都落在同一个环里**。
        for smi in [
            "CC(=O)C1=CC2=C(C=C1C(C)=O)[C]3(C)CCC[C](C)(C#N)[CH]3CC2=O", // 语料第 2719 行
            "c1ccc2ccccc2c1",                                            // 萘
            "C1=CC2=CC=CC3=CC=CC(=C1)C23",                               // 苊
        ] {
            for style in &Style::ALL {
                let m = prep(smi);
                let d = generate(&m, style);
                let s = scene(&m, &d, style);
                // `coords`/`wedges` 的下标相对**被画的那个分子** —— 为画出构型补
                // 出来的氢也在里面。`scene` 自己会补,所以它拿的是原分子;这里影子化
                // 之后,下面按下标索引才不会错位甚至越界。
                let m = d.drawn(&m);
                let bnd = bounds(&d.coords, &m, style);
                let pts: Vec<Point2> = d
                    .coords
                    .iter()
                    .map(|p| to_canvas(*p, bnd, style.bond_length_pt))
                    .collect();
                let rings = omgkit_chem::sssr::ring_set(&m);
                let orders = drawn_orders(&m);

                for (bi, b) in m.bonds().iter().enumerate() {
                    let no = u32::try_from(bi).expect("键数超出 u32");
                    let mine: Vec<_> = rings.iter().filter(|r| r.bonds.contains(&no)).collect();
                    if orders[bi] != BondOrder::Double || mine.len() < 2 {
                        continue;
                    }
                    let (pa, pb) = (pts[b.begin as usize], pts[b.end as usize]);
                    let got = lines_of_bond(&s, pa, pb);
                    assert_eq!(
                        got.len(),
                        2,
                        "[{}] {smi}:共用键 {bi} 画了 {} 条线",
                        style.name,
                        got.len()
                    );
                    // 两条线的中点必须落在**同一个**环里
                    let home = |p: Point2| -> Vec<usize> {
                        mine.iter()
                            .enumerate()
                            .filter(|(_, r)| {
                                let poly: Vec<Point2> =
                                    r.atoms.iter().map(|a| pts[*a as usize]).collect();
                                crate::geom::point_in_polygon(p, &poly)
                            })
                            .map(|(k, _)| k)
                            .collect()
                    };
                    let h0 = home((got[0].0 + got[0].1) * 0.5);
                    let h1 = home((got[1].0 + got[1].1) * 0.5);
                    assert!(
                        h0.iter().any(|x| h1.contains(x)),
                        "[{}] {smi}:共用键 {bi}({}–{})的两条线分属不同的环 {h0:?} / {h1:?} —— 跨骑在键上了",
                        style.name,
                        b.begin,
                        b.end
                    );
                }
            }
        }
    }

    #[test]
    fn a_double_bond_with_no_inner_side_is_drawn_symmetric() {
        // **只有真的没有内侧可偏时**,两条线才对称跨在键轴两侧:
        //
        // - 累积双键 `O=C=O`:两根键若都把第二条线偏到同一侧,画出来就是一条
        //   直线配两条同侧短线 —— 读起来是顺式二烯;
        // - 内侧原子是分叉点 `CC(C)=C`、`CC(C)=O`:两侧各被一根键占着,偏哪边
        //   都压着一根键;
        // - 两端都是端基 `C=C`:压根没有别的邻居。
        //
        // **"一端是端基"不在此列** —— 那一类见
        // [`a_terminal_double_bond_closes_the_joint_at_the_inner_atom`]。
        //
        // 旧写法拿 `f64::signum` 数邻居,而 `signum` 对 ±0.0 给 ±1 —— 共线的
        // 邻居会按零的符号位投出一票,那一位取决于算到那步的运算次序。
        //
        // **后两个是补进来的**,因为前面那些**挡不住**"把端基那条早退整个删掉"
        // 这个改法:`CC(C)=C`、`CC(C)=O` 的内侧原子度 3、两个邻居分居两侧,
        // 投票碰巧抵消,照样对称 —— 绿得没有道理。
        //
        // - `CN=O`:内侧原子度 2 但**带标签**。那一头的键都停在字盒外,没有
        //   顶点要合;偏一边只会让 `=` 挂在字母底边上。全量 100 根/95 分子。
        // - `Br[CH]1[CH](Br)S(=O)(=O)C2=C1C=CC=C2`:内侧的 S **度 4**,三个别的
        //   邻居,票数是奇数项之和,抵消不掉。删掉早退的话两根 S=O 双双朝对方
        //   倾斜,四条线挤在中间。全量 116 根/68 分子。
        // - `CN=NC`(偶氮):**链上**的双键,两端都有取代基、票数抵消,而
        //   **两端都带标签** —— 那时键都停在字盒外,同样没有顶点要合。这一支
        //   先前一条判据都没有,删掉它全世界不吭声(全量审计也逐字节不变,
        //   因为审计里没有任何一条性质在量偏侧)。删掉之后 `=` 贴着 N 画。
        for smi in [
            "O=C=O",
            "CC(C)=C",
            "C=C",
            "CC(C)=O",
            "CN=O",
            "Br[CH]1[CH](Br)S(=O)(=O)C2=C1C=CC=C2",
            "CN=NC",
        ] {
            let m = prep(smi);
            let style = Style::ACS_1996;
            let d = generate(&m, &style);
            let s = scene(&m, &d, &style);
            let bnd = bounds(&d.coords, &m, &style);
            let pts: Vec<Point2> = d
                .coords
                .iter()
                .map(|p| to_canvas(*p, bnd, style.bond_length_pt))
                .collect();

            let mut checked = 0;
            for b in m.bonds() {
                if b.order != BondOrder::Double {
                    continue;
                }
                let (pa, pb) = (pts[b.begin as usize], pts[b.end as usize]);
                let len = pa.dist(pb);
                let mid = (pa + pb) * 0.5;
                let axis = (pb - pa) * (1.0 / len);
                let normal = Point2::new(-axis.y, axis.x);
                let got = lines_of_bond(&s, pa, pb);
                assert_eq!(got.len(), 2, "{smi}:双键画出了 {} 条线", got.len());
                let sides: Vec<f64> = got
                    .iter()
                    .map(|(f, t)| ((*f + *t) * 0.5 - mid).dot(normal))
                    .collect();
                assert!(
                    sides[0] * sides[1] < 0.0
                        && (sides[0].abs() - sides[1].abs()).abs() < 0.02 * len,
                    "{smi}:两条线没有对称跨在键轴两侧,偏移分别是 {:.3}、{:.3}",
                    sides[0],
                    sides[1]
                );
                checked += 1;
            }
            assert!(checked >= 1, "{smi}:一根双键都没查到");
        }
    }

    #[test]
    fn a_terminal_double_bond_closes_the_joint_at_the_inner_atom() {
        // 丙烯 `CH₂=CH–CH₃`:内侧那个碳上还挂着一根单键,而单键是**收在原子
        // 中心**的。两条线若跨着键轴对称画,谁也不从那一点出发 —— 顶点合不拢,
        // 单键的尖头戳出来,旁边留一个豁口。放大看一眼就明白。
        //
        // 通例(也是 RDKit `doubleBondTerminal` 的第三支):**一条线走键轴**把
        // 顶点接上,另一条偏向内侧原子的**另一根键**那一侧,并在末端那头与主线
        // **齐头**收尾 —— 末端没有第二根键可接,缩进去只会让它看着短一截。
        //
        // 全量语料 404 根键、374 个分子落在这一档。醛 `CC=O` 同理:RDKit 画的
        // 乙醛也是一条线走键轴。
        for smi in ["CC=C", "CC=O", "CCC=C", "C=CC#N"] {
            for style in &Style::ALL {
                let m = prep(smi);
                let d = generate(&m, style);
                let s = scene(&m, &d, style);
                // `coords`/`wedges` 的下标相对**被画的那个分子** —— 为画出构型补
                // 出来的氢也在里面。`scene` 自己会补,所以它拿的是原分子;这里影子化
                // 之后,下面按下标索引才不会错位甚至越界。
                let m = d.drawn(&m);
                let bnd = bounds(&d.coords, &m, style);
                let pts: Vec<Point2> = d
                    .coords
                    .iter()
                    .map(|p| to_canvas(*p, bnd, style.bond_length_pt))
                    .collect();

                let mut checked = 0usize;
                for (bi, b) in m.bonds().iter().enumerate() {
                    if drawn_orders(&m)[bi] != BondOrder::Double {
                        continue;
                    }
                    // 这一档:恰有一端是端基,内侧那端度 2 且不带标签
                    let (da, db) = (m.degree(b.begin), m.degree(b.end));
                    if (da == 1) == (db == 1) {
                        continue;
                    }
                    let (term, inner) = if da == 1 {
                        (b.begin, b.end)
                    } else {
                        (b.end, b.begin)
                    };
                    if m.degree(inner) != 2 || label_at(&m, inner, style, &d.coords).is_some() {
                        continue;
                    }
                    let third = m
                        .neighbors(inner)
                        .map(|(x, _)| x)
                        .find(|x| *x != term)
                        .expect("度 2 的内侧原子必有另一个邻居");

                    let (pi, pt) = (pts[inner as usize], pts[term as usize]);
                    let len = pi.dist(pt);
                    let axis = (pt - pi) * (1.0 / len);
                    let normal = Point2::new(-axis.y, axis.x);
                    let ls = lines_of_bond(&s, pi, pt);
                    assert_eq!(ls.len(), 2, "[{}] {smi}:双键该画两条线", style.name);

                    // 一、有一条线从内侧原子出发 —— 顶点合得拢。
                    //
                    // **这条只等价于"没被画成对称"**:内侧原子不带标签,`trim`
                    // 不动端点,所以主线必然从 `pi` 出发。它管的是 `offset_dir`
                    // 那一档,管不到第二条线怎么收尾 —— 那要靠下面第四条。
                    let starts_at_inner =
                        |(u, v): &(Point2, Point2)| u.dist(pi).min(v.dist(pi)) < 0.02 * len;
                    assert!(
                        ls.iter().any(starts_at_inner),
                        "[{}] {smi}:键 {bi} 两条线都不从内侧原子 {inner} 出发,顶点合不拢",
                        style.name
                    );

                    // 二、另一条偏向内侧原子的另一根键那一侧
                    let off = |(u, v): &(Point2, Point2)| ((*u + *v) * 0.5 - pi).dot(normal);
                    let side = (pts[third as usize] - pi).dot(normal);
                    let outer = ls
                        .iter()
                        .max_by(|x, y| off(x).abs().partial_cmp(&off(y).abs()).expect("坐标非 NaN"))
                        .expect("有两条线");
                    assert!(
                        off(outer) * side > 0.0,
                        "[{}] {smi}:键 {bi} 的第二条线偏到了内侧原子另一根键的反面",
                        style.name
                    );

                    // 三、末端那头齐头 —— 那里没有第二根键可接,不该缩进去
                    let far = |(u, v): &(Point2, Point2)| {
                        let (tu, tv) = ((*u - pi).dot(axis), (*v - pi).dot(axis));
                        tu.max(tv)
                    };
                    let (f0, f1) = (far(&ls[0]), far(&ls[1]));
                    assert!(
                        (f0 - f1).abs() < 0.02 * len,
                        "[{}] {smi}:键 {bi} 两条线在端基那头没齐头,{f0:.2} vs {f1:.2}(键长 {len:.2})",
                        style.name
                    );

                    // 四、第二条线在内侧那头**斜切到角平分线上**,落点有闭式解:
                    // 与内侧原子的距离是 `spacing / sin(θ/2)`,θ 是那里的夹角。
                    //
                    // 这一条是独立的。前三条都只等价于"没画成对称",**把斜切
                    // 整个去掉它们照样全绿** —— 而那时第二条线会穿过相邻的单键
                    // 伸出去(实测丙烯伸出 1.30pt,约 0.09 个键长)。
                    let u = (pts[third as usize] - pi).normalized();
                    let v = (pt - pi).normalized();
                    let bis = (u + v).normalized();
                    let half = u.dot(v).clamp(-1.0, 1.0).acos() / 2.0;
                    let spacing = style.bond_spacing() * style.bond_length_pt;
                    let want = pi + bis * (spacing / half.sin());
                    let inner_end = {
                        let (uu, vv) = *outer;
                        if uu.dist(pi) < vv.dist(pi) {
                            uu
                        } else {
                            vv
                        }
                    };
                    assert!(
                        inner_end.dist(want) < 0.02 * len,
                        "[{}] {smi}:键 {bi} 的第二条线没斜切到角平分线上 —— \
                         落在 ({:.2},{:.2}),闭式解是 ({:.2},{:.2})",
                        style.name,
                        inner_end.x,
                        inner_end.y,
                        want.x,
                        want.y
                    );
                    checked += 1;
                }
                assert!(
                    checked >= 1,
                    "[{}] {smi}:这一档一根键都没查到 —— 判据空过了",
                    style.name
                );
            }
        }
    }

    #[test]
    fn a_trans_double_bond_closes_both_joints() {
        // 反式双键 `CH₃–CH=CH–CH₃`:两个取代基一边一个,按邻居计数正好抵消。
        // 先前抵消就跨轴对称画 —— 于是**两个顶点都合不拢**,与丙烯是同一个
        // 毛病,只是这里两头都犯。画廊里的 trans-butene 就是这样。
        //
        // 全量语料落进这一档的有 **629 根键 / 562 个分子**(单套规范,按审计
        // 那条会感知顺反的管线数;这里的 `prep` 不感知,同一份语料数出 643/574)。
        // 其中 93 根**两端都带标签**,仍旧对称,由
        // [`a_double_bond_with_no_inner_side_is_drawn_symmetric`] 管。真正改了
        // 画法的是 **536 根**。
        //
        // 偏哪一侧没有道理可讲,但**必须偏**:一条线走键轴才接得上两端的单键。
        //
        // **只挑真的落在这一档的分子。** 先前列的 `CC=C(C)CC` 走的是 vote 档
        // (一端一个邻居、另一端两个分居两侧,票数 ±1),在 HEAD 上本来就不
        // 对称,变异下照样绿 —— 对这条判据零贡献。下面按几何显式筛:
        // **两端各恰好一个取代基,且分居键轴两侧**。
        //
        // 断言拿的是**未 trim** 的原子位置,所以列的分子两端都不能带标签 ——
        // 带标签那头主线会缩回字盒外,断言会假红。真实语料里这一档 536 根中
        // 有 246 根是"只有一端带标签",那时只有一个顶点要合,同样该偏,只是
        // 这条判据的写法量不了。
        for smi in ["C/C=C/C", "C/C=C/CC", "C/C=C/C=C/C=C/C"] {
            for style in &Style::ALL {
                let m = prep(smi);
                let d = generate(&m, style);
                let s = scene(&m, &d, style);
                // `coords`/`wedges` 的下标相对**被画的那个分子** —— 为画出构型补
                // 出来的氢也在里面。`scene` 自己会补,所以它拿的是原分子;这里影子化
                // 之后,下面按下标索引才不会错位甚至越界。
                let m = d.drawn(&m);
                let bnd = bounds(&d.coords, &m, style);
                let pts: Vec<Point2> = d
                    .coords
                    .iter()
                    .map(|p| to_canvas(*p, bnd, style.bond_length_pt))
                    .collect();

                let mut checked = 0usize;
                for (bi, b) in m.bonds().iter().enumerate() {
                    if drawn_orders(&m)[bi] != BondOrder::Double {
                        continue;
                    }
                    if m.degree(b.begin) == 1 || m.degree(b.end) == 1 {
                        continue;
                    }
                    let (pa, pb) = (pts[b.begin as usize], pts[b.end as usize]);
                    let len = pa.dist(pb);
                    let mid = (pa + pb) * 0.5;
                    let axis = (pb - pa) * (1.0 / len);
                    let normal = Point2::new(-axis.y, axis.x);
                    // **这一档的定义:两端各恰好一个取代基,且分居键轴两侧。**
                    // 同侧是顺式,走的是另一档(偏向取代基那侧),不在这里。
                    let subs = |e: u32, other: u32| -> Vec<u32> {
                        m.neighbors(e)
                            .map(|(x, _)| x)
                            .filter(|x| *x != other)
                            .collect()
                    };
                    let (sa, sb) = (subs(b.begin, b.end), subs(b.end, b.begin));
                    if sa.len() != 1 || sb.len() != 1 {
                        continue;
                    }
                    let side = |x: u32| (pts[x as usize] - mid).dot(normal);
                    if side(sa[0]) * side(sb[0]) >= 0.0 {
                        continue;
                    }
                    // 两端都带标签的另有规矩(偶氮 `CN=NC`),不在这条里
                    if label_at(&m, b.begin, style, &d.coords).is_some()
                        && label_at(&m, b.end, style, &d.coords).is_some()
                    {
                        continue;
                    }
                    let ls = lines_of_bond(&s, pa, pb);
                    assert_eq!(ls.len(), 2, "[{}] {smi}:双键该画两条线", style.name);
                    // 有一条线整根走在键轴上 —— 两个顶点都接得上
                    assert!(
                        ls.iter().any(|(u, v)| {
                            (u.dist(pa) < 0.02 * len && v.dist(pb) < 0.02 * len)
                                || (u.dist(pb) < 0.02 * len && v.dist(pa) < 0.02 * len)
                        }),
                        "[{}] {smi}:键 {bi} 没有一条线走在键轴上,两个顶点都合不拢",
                        style.name
                    );
                    checked += 1;
                }
                assert!(
                    checked >= 1,
                    "[{}] {smi}:一根都没查到 —— 判据空过了",
                    style.name
                );
            }
        }
    }

    #[test]
    fn nothing_is_drawn_outside_the_canvas() {
        // 画布尺寸只按原子位置和标签算的话,**双键的第二条线与线宽都没算进去** ——
        // 它们会伸到画布外被裁掉。图看着基本正常,只是边上少一截,极容易漏过。
        for smi in [
            "CC(=O)Nc1ccc(O)cc1",
            "c1ccccc1",
            "CC(=O)Oc1ccccc1C(=O)O",
            "CC#N",
            "CN1C=NC2=C1C(=O)N(C)C(=O)N2C",
            "O=C=O",
            // 补出来的符号也要有地方放:这个分子有个共线的骨架碳,`scene` 给它
            // 画了 `CH`,而 [`bounds`] 先前只问 `label_for`,当它不存在 ——
            // 半宽 7.22pt 的字直接戳出画布右边。全量语料上 4 处。
            "c1ccc2c(c1)[C@@H]3CC[C@H]2[n+]4c3cccc4",
        ] {
            for style in &Style::ALL {
                let m = prep(smi);
                let d = generate(&m, style);
                let s = scene(&m, &d, style);
                // `coords`/`wedges` 的下标相对**被画的那个分子** —— 为画出构型补
                // 出来的氢也在里面。`scene` 自己会补,所以它拿的是原分子;这里影子化
                // 之后,下面按下标索引才不会错位甚至越界。
                let m = d.drawn(&m);
                for it in &s.items {
                    // 文字走下面那段 —— **字宽不是 `size/2`**,`CH` 的半宽 7.22pt
                    // 比字号的一半还大 2.22pt,按 `size/2` 量会把出界量少算
                    let pts: Vec<(Point2, f64)> = match it {
                        Primitive::Line { from, to, width } => {
                            vec![(*from, *width / 2.0), (*to, *width / 2.0)]
                        }
                        Primitive::Wedge { from, to, wide }
                        | Primitive::Hash { from, to, wide, .. } => {
                            vec![(*from, *wide / 2.0), (*to, *wide / 2.0)]
                        }
                        Primitive::Text { .. } => continue,
                    };
                    for (p, r) in pts {
                        assert!(
                            p.x - r >= -0.01 && p.x + r <= s.width + 0.01,
                            "[{}] {smi}:图元 x={:.2}(±{r:.2})超出画布宽 {:.2}",
                            style.name,
                            p.x,
                            s.width
                        );
                        assert!(
                            p.y - r >= -0.01 && p.y + r <= s.height + 0.01,
                            "[{}] {smi}:图元 y={:.2}(±{r:.2})超出画布高 {:.2}",
                            style.name,
                            p.y,
                            s.height
                        );
                    }
                }

                let bnd = bounds(&d.coords, &m, style);
                let scale = style.bond_length_pt;
                for a in 0..u32::try_from(m.num_atoms()).unwrap() {
                    let Some(l) = label_at(&m, a, style, &d.coords) else {
                        continue;
                    };
                    let c = to_canvas(d.coords[a as usize], bnd, scale)
                        + Point2::new(l.dx * scale, 0.0);
                    let (hw, hh) = (l.half_w * scale, l.half_h * scale);
                    assert!(
                        c.x - hw >= -0.01
                            && c.x + hw <= s.width + 0.01
                            && c.y - hh >= -0.01
                            && c.y + hh <= s.height + 0.01,
                        "[{}] {smi}:标签 {} 在 ({:.2},{:.2})±({hw:.2},{hh:.2}),画布 {:.2}×{:.2}",
                        style.name,
                        l.plain(),
                        c.x,
                        c.y,
                        s.width,
                        s.height
                    );
                }
            }
        }
    }

    #[test]
    fn a_mismatched_depiction_is_refused_loudly() {
        // 拿错图不该画出一张张冠李戴的结构式
        let a = prep("CCO");
        let b = prep("c1ccccc1");
        let d = generate(&a, &Style::ACS_1996);
        let r = std::panic::catch_unwind(|| scene(&b, &d, &Style::ACS_1996));
        assert!(r.is_err(), "原子数不符却照画不误");
    }
}
