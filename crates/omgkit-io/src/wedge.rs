//! 结构图上的**楔形键**,以及从图上反读一个中心的手性。
//!
//! # 为什么住在 L1(io)而不是画图那一层
//!
//! 楔形不是画图的产物,是**文件格式里的一个字段** —— molblock 的键块第四列
//! 就是它(1 实楔、6 虚楔)。读一个 `.mol` 文件要靠它定手性,画一张图要产生它,
//! 两件事共用同一套语义。
//!
//! 先前这段代码住在 `omgkit-depict`,于是 `omgkit-io` 的 molblock 读取器够不着
//! (depict 依赖 io,反过来不行),只好自己再造一个楔形类型 —— 一个概念两个
//! 类型,迟早分岔。搬下来之后 depict 从这里 `pub use`,两边是同一个东西。
//!
//! # 判据的口径
//!
//! [`chirality_from_wedges`] **只看几何**,不看 `chiral_tag`。它与 depict 的
//! `assign_wedges` 合起来是一次往返;单独看,它回答的是"图上画出来的构型是什么"。
//!
//! # 打标记那几个入口不在这里
//!
//! `assign_*` 四个(二维/三维 × 手性/顺反)全在 [`crate::stereo`]。本模块只管
//! **楔形这个字段怎么读**;"给分子打立体标记"是另一件事,而它有四个口,分散在
//! 两个模块里的话,下一个人只会找到其中一半。

use omgkit_core::{ChiralTag, MolBuilder};

/// 一根键的楔形指派。**窄端在立体中心**,宽端在另一头 —— 这是画结构式的通例:
/// 楔形描述的是"从这个中心看出去",窄端标出了那个中心。
///
/// # 窄端记在类型里,不靠猜
///
/// 两个立体中心相邻时,它们**共用**的那根键两头都是中心,"窄端在哪"就不能
/// 靠"哪头带手性标记"去猜了。猜错的后果是这根键描述的构型整个反过来 ——
/// 而线条本身看着一点毛病没有。实测:抗坏血酸的两个中心正好相邻。
///
/// 所以窄端跟着枚举一起带上,让"有楔形却不知道窄端在哪"这个状态根本表示不出来。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Wedge {
    /// 普通实线
    #[default]
    None,
    /// 实楔形:指向观察者
    Up {
        /// 窄端所在的原子
        narrow: u32,
    },
    /// 虚楔形:背离观察者
    Down {
        /// 窄端所在的原子
        narrow: u32,
    },
}

impl Wedge {
    /// 窄端所在的原子。不是楔形时没有窄端。
    #[must_use]
    pub fn narrow(self) -> Option<u32> {
        match self {
            Wedge::None => None,
            Wedge::Up { narrow } | Wedge::Down { narrow } => Some(narrow),
        }
    }
}

/// 三个画出来的邻居张不出这么大的体积,就判"这张图定不出手性"。
///
/// 单位是**键长的立方**。数值取自 RDKit `Chirality.cpp` 的 `ZERO_VOLUME_TOL`,
/// 照抄是因为**我们导出的 molblock 就是交给它读的** —— 两边用同一把尺,
/// "别人读得回来"才谈得上。
pub(crate) const ZERO_VOLUME_TOL: f64 = 0.1;

enum Ligand {
    Atom(u32, u32),
    ImplicitH,
}

/// 总氢数。两个字段互斥,相加即总数 —— 全仓一致的约定。
pub(crate) fn total_hs(mol: &MolBuilder, atom: u32) -> u8 {
    let a = mol.atoms()[atom as usize];
    a.num_explicit_hs.saturating_add(a.num_implicit_hs)
}

/// 这个原子的第四个配体是不是一对**孤对电子**。
///
/// 三配位的 S、Se、P、As 上那对孤对与隐式氢是同一件事:看不见,但占一个配体
/// 位,构型因此确定。亚砜 `R–S(=O)–R′`、亚膦酸酯都是这样。
///
/// **氮不算。** 三配位氮的孤对翻转极快(氨的翻转垒只有 24 kJ/mol),常温下
/// 两个构型互变,画出楔形是在断言一个不存在的构型。少数被环卡住的(氮杂环丙烷、
/// Tröger 碱)确实稳定,但那要看环张力,不是看元素 —— 这里宁可**漏报**也不
/// 乱画,漏了如实进 `unwedged`。
///
/// 带正电的不算:季铵、锍盐上没有孤对,四个配体全在,走 `(4, 0)` 那一支。
pub(crate) fn has_lone_pair(mol: &MolBuilder, a: u32) -> bool {
    // **表在 `omgkit-core`,这里不再自己写一份。** `omgkit-conf` 抽手性中心时
    // 要问同一个问题(三个邻居也算数),两处各写一份迟早分岔,
    // 而分岔的表现是一半的中心画对了、另一半摆错了。
    let at = mol.atoms()[a as usize];
    omgkit_core::element::has_stereogenic_lone_pair(at.atomic_num, at.formal_charge)
}

/// 从坐标与楔形反读一个中心的手性。判不出来返回 `None`。
///
/// `coords` 逐原子,二维图把 `z` 全给 0 —— 这个函数只用 x、y,z 由楔形给。
///
/// # 判据的口径
///
/// 这个函数**只看几何**,不看 `chiral_tag`。它与 depict 的 `assign_wedges` 合起来是
/// 一次往返;单独看,它是"图上画出来的构型是什么"的答案。
#[must_use]
pub fn chirality_from_wedges(
    mol: &MolBuilder,
    coords: &[[f64; 3]],
    wedges: &[Wedge],
    a: u32,
) -> Option<ChiralTag> {
    let nbrs: Vec<(u32, u32)> = mol.neighbors(a).collect();
    let hs = total_hs(mol, a);

    // 参照序:四个配体。三个邻居加一个隐式氢时,氢占**槽位 1** —— 见模块文档。
    let mut refs: Vec<Ligand> = Vec::with_capacity(4);
    match (nbrs.len(), hs) {
        (4, 0) => {
            for (n, bi) in &nbrs {
                refs.push(Ligand::Atom(*n, *bi));
            }
        }
        // 三根键 + 一个**看不见的第四配体**。两种情形是同一件事:
        //
        // - `(3, 1)` 隐式氢:碳上的常例;
        // - `(3, 0)` **孤对电子**:亚砜的 S、亚砜亚胺、膦氧化物的 P …… 三配位
        //   加一对孤对,构型照样确定,画法也一样 —— 楔形打在三根键之一上,
        //   孤对在它的反面。
        //
        // 先前 `(3, 0)` 直接返回 `None`,于是这些中心一律进 `unwedged`:全量
        // 语料 18 个画不出构型的中心里,**14 个是这一档**。
        (3, 1) | (3, 0) if hs == 1 || has_lone_pair(mol, a) => {
            refs.push(Ligand::Atom(nbrs[0].0, nbrs[0].1));
            refs.push(Ligand::ImplicitH);
            refs.push(Ligand::Atom(nbrs[1].0, nbrs[1].1));
            refs.push(Ligand::Atom(nbrs[2].0, nbrs[2].1));
        }
        _ => return None, // 不是四面体中心,或者配体数说不通
    }

    let centre = coords[a as usize];
    let mut pts: Vec<[f64; 3]> = Vec::with_capacity(4);
    let mut h_slot = None;
    for (k, l) in refs.iter().enumerate() {
        match l {
            Ligand::Atom(n, bi) => {
                let p = coords[*n as usize];
                // 窄端在别人那头的楔形要**反过来读**:那根键说的是"从 n 看出去
                // a 在上面",换成从 a 看,n 就在下面。不翻的话相邻的两个立体
                // 中心里,后读的那个会得到对映体。
                let z = match wedges[*bi as usize] {
                    Wedge::Up { narrow } => {
                        if narrow == a {
                            1.0
                        } else {
                            -1.0
                        }
                    }
                    Wedge::Down { narrow } => {
                        if narrow == a {
                            -1.0
                        } else {
                            1.0
                        }
                    }
                    Wedge::None => 0.0,
                };
                pts.push([p[0] - centre[0], p[1] - centre[1], z]);
            }
            Ligand::ImplicitH => {
                h_slot = Some(k);
                pts.push([0.0, 0.0, 0.0]); // 占位,下面补
            }
        }
    }

    // 隐式氢没有画出来,得把它摆到该在的地方。所有楔形都没有的话,这个中心的
    // 构型在图上根本读不出来 —— 如实返回 None,不猜。
    //
    // 这条只管有隐式氢的情形。四个邻居都画出来的中心不需要摆氢,楔形正负抵消
    // 也照样读得出来 —— 一并早退的话,那种中心会被误判成"读不出",而它的
    // 几何其实一点不含糊。
    //
    // # 面内那一份不能一律取零
    //
    // 先前摆的是 `[0, 0, -zsum]` —— 认为氢**投影到中心上**,只在楔形的反面。
    // 那只在**三个邻居把中心围住**时成立。
    //
    // 四面体的四个键方向之和为零,所以**它们的 2D 投影之和也为零**:三个投影
    // 若全落在中心的同一侧(最大空隙 > 180°),第四个必然落在**对面的空扇区**
    // 里,不可能落在中心上。桥环退化布局里就会出现这种中心。
    //
    // 实测(`large.smi` 全量,量的是**最终交付的图**:导成 molblock 之后按 2D
    // 坐标算每个三邻居中心的最大空隙):273 个楔形窄端三邻居中心里 **2 个**
    // 空隙 > 180°(291.2° 与 234.6°,都是桥环);换成"RDKit 读得出 CIP 码的
    // 中心"这个口径是 271 里 2 个,同一批。
    //
    // (这里原先写着"全量 420 个三邻居中心里 45 个如此,张角只有 120–137°"。
    //  那个数**两种口径都复现不出来**,而它没写清自己数的是哪一批 ——
    //  多半是布局改进之前留下的。数字进注释就要连口径一起写,否则下一个人
    //  没法判断它是不是过期了。)
    //
    // 摆错的后果是**读出对映体**,而线条本身毫无毛病。外部判官
    // (`harness/check_wedge_readback.py`,拿 RDKit 从导出的 molblock 反读)量到
    // **21 个中心画成了对映体而且没进 `unwedged`** —— 全部出在这里。
    //
    // 改成按同一条恒等式摆:面内取三个邻居**单位方向之和的负值**。三个邻居
    // 均匀分布时这个和恰为零,退化成先前的做法;挤在一侧时它指进空扇区。
    // 一个公式管两种情形,中间是连续的。
    if let Some(k) = h_slot {
        let zsum: f64 = pts.iter().map(|p| p[2]).sum();
        if zsum.abs() < 1e-9 {
            return None;
        }
        // # 先问一句:这三根键**定得出**手性吗
        //
        // 三个画出来的邻居里若有两根几乎**共线**,它们张不出面积,楔形再怎么
        // 画也定不出手性 —— 与氢摆在哪无关。判据就是三者的三重积
        // `v1·(v2×v3)`(键长 1 为单位),它等于"底面的 2D 叉积 × 楔形那一维"。
        //
        // 这条口径与 RDKit 的 `Chirality.cpp` 一致:那里同样不摆氢,直接算三个
        // 邻居的三重积,`|vol| ≤ ZERO_VOLUME_TOL = 0.1` 就判 `CHI_UNSPECIFIED`。
        // 常数照抄,因为**我们的 molblock 就是交给它读的**,两边用同一把尺才谈
        // 得上"别人读得回来"。
        //
        // 实测全量语料 404 个三邻居中心里只有 1 个落在这条线下(|vol| = 0.0905,
        // 那两根键相差 184.4°),下一个是 1.24 —— 中间空 13 倍。
        let v: Vec<[f64; 3]> = pts
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != k)
            .map(|(_, p)| *p)
            .collect();
        let cross = [
            v[1][1] * v[2][2] - v[1][2] * v[2][1],
            v[1][2] * v[2][0] - v[1][0] * v[2][2],
            v[1][0] * v[2][1] - v[1][1] * v[2][0],
        ];
        let vol = v[0][0] * cross[0] + v[0][1] * cross[1] + v[0][2] * cross[2];
        if vol.abs() <= ZERO_VOLUME_TOL {
            return None;
        }

        let mut sx = 0.0_f64;
        let mut sy = 0.0_f64;
        for p in &v {
            let n = p[0].hypot(p[1]);
            if n < 1e-9 {
                return None; // 邻居与中心重合,方向无从谈起
            }
            sx += p[0] / n;
            sy += p[1] / n;
        }
        pts[k] = [-sx, -sy, -zsum.signum()];
    }

    // 有向体积:det[p1-p0, p2-p0, p3-p0]。符号的含义由
    // `the_reference_tetrahedron_pins_the_sign` 那条测试钉死。
    let d = |i: usize, j: usize| pts[i][j] - pts[0][j];
    let det = d(1, 0) * (d(2, 1) * d(3, 2) - d(2, 2) * d(3, 1))
        - d(1, 1) * (d(2, 0) * d(3, 2) - d(2, 2) * d(3, 0))
        + d(1, 2) * (d(2, 0) * d(3, 1) - d(2, 1) * d(3, 0));

    if det.abs() < 1e-12 {
        None
    } else if det < 0.0 {
        Some(ChiralTag::Ccw)
    } else {
        Some(ChiralTag::Cw)
    }
}

#[cfg(test)]
mod tests {

    /// 一张画着楔形的二维图,以及把同一张图的某个 `z` 抬起来之后的样子。
    const WEDGED: &str = "\
C[C@H](N)O
     RDKit          2D

  4  3  0  0  0  0  0  0  0  0999 V2000
   -1.2990   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990   -0.7500    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.5000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  1
  2  3  1  0
  2  4  1  0
M  END
";

    fn tagged(block: &str) -> usize {
        let got = crate::molblock::read_v2000(block).expect("读 molblock");
        let mut m = got.mol;
        omgkit_chem::pipeline::sanitize(&mut m).expect("净化");
        crate::stereo::assign_chirality_2d(&mut m, &got.coords, &got.wedges)
    }

    /// 二维那张图读得出中心 —— 这是下一条测试的对照,少了它"三维读不出"
    /// 可能只是因为这张图本来就没有立体。
    #[test]
    fn a_wedge_on_a_flat_drawing_gives_a_centre() {
        assert_eq!(tagged(WEDGED), 1);
    }

    /// 同一张图把一个 `z` 抬起来,就整个不做。
    ///
    /// 三维文件里偶尔留着楔形字段,那时按 xy 投影算出来的体积与分子无关 ——
    /// 空答案可以接受,错答案不行。
    #[test]
    fn the_same_drawing_lifted_into_3d_gives_nothing() {
        let lifted = WEDGED.replace(
            "    0.0000    1.5000    0.0000 O",
            "    0.0000    1.5000    0.9000 O",
        );
        assert_ne!(lifted, WEDGED, "原子块那一行没改到");
        assert_eq!(tagged(&lifted), 0);
    }
}
