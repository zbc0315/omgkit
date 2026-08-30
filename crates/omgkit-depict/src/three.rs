//! **三维分子图** —— 把一组三维坐标画成球棍/空间填充/棍状/线框四种图之一。
//!
//! ```
//! use omgkit_depict::three::{self, Style3D};
//!
//! let mut m = omgkit_io::smiles::parse("CCO").unwrap();
//! let c = omgkit_conf::pipeline::conformer_for(&mut m).unwrap();
//! let d = three::depict(&m, &c.coords, &Style3D::BALL_AND_STICK).unwrap();
//! let svg = omgkit_depict::svg::to_svg(&d.scene, &omgkit_depict::style::Style::ACS_1996);
//! assert!(svg.starts_with("<svg"));
//! ```
//!
//! # 四种样式不是这里发明的
//!
//! [`Style3D`] 的四个常量、连同它们的半径,逐项取自 Jmol 自己文档里的
//! "standard rendering styles"(Jmol Wiki, *Rendering*):
//!
//! | | Jmol 的写法 | 球半径 | 杆半径 |
//! |---|---|---|---|
//! | 空间填充 | `spacefill 100%` | 100% vdW | — |
//! | 球棍 | `wireframe 0.15; spacefill 23%` | 23% vdW | 0.15 Å |
//! | 棍状 | `wireframe 0.3; spacefill off` | = 杆半径 | 0.30 Å |
//! | 线框 | `wireframe 0.01; spacefill 0` | — | 0.01 Å |
//!
//! 范德华半径取 [`omgkit_core`] 元素表里的 `rvdw`(源头是 BODR),不另存一份。
//! 颜色取 [`palette`](crate::palette) 里的 Jmol CPK 表,**键的两半各随自己
//! 那一端的颜色** —— Jmol、PyMOL、VMD、3Dmol.js 都是这么画的。
//!
//! # 与二维那条路的三处不同
//!
//! **一、朝向不许镜像。** 二维的 [`orient`](crate::orient) 在 24 个候选姿态
//! 里挑,其中一半是镜像的 —— 那在二维是安全的,因为构型是靠楔形画出来的,
//! 镜像之后重新指派楔形就行。三维图的构型**就是坐标本身**:镜一下,
//! 手性中心全反,而图上没有任何一处看得出来。所以这里的视角矩阵**行列式恒为
//! +1**,这一条有判据钉着。
//!
//! **二、朝向不做 30° 对齐。** 二维那边对齐是为了保住环上齐整的 120°;
//! 三维图里分子本来就不躺在任何格点上,对齐只会让一个平面分子歪着 15° ——
//! 那个理由不迁移。这里用主轴:方差最大的轴放水平,最小的轴指向观察者
//! (PyMOL 的 `orient`、RDKit 的 `ComputeCanonicalTransform` 都是这么定的)。
//!
//! **三、遮挡靠画家算法。** SVG 没有深度缓冲,只能按深度从远到近画。球按球心
//! 深度排;键劈成两半(各随一端的颜色)之后**再按深度切片**,切到每片的深度
//! 跨度不超过 [`DEPTH_SLICE`] —— 不切的话,一根斜着穿过球的棍会整根压在球上
//! 或者整根被球压住,而两者都不对。
//!
//! # 视角退化会说出来
//!
//! 对称性高的分子(苯、金刚烷、CH₄)主轴不唯一:两个特征值相等时,那两根轴
//! 在它们张成的平面里可以任意转,选哪一根取决于浮点噪声的最后一位。
//! [`View::degenerate`] 如实记这一笔 —— 这类分子换一种 SMILES 写法,画出来
//! **可能不是同一张图**,而图本身看不出毛病。

use omgkit_conf::linalg::symmetric_eigen;
use omgkit_core::{BondOrder, MolBuilder};

use crate::geom::Point2;
use crate::palette::cpk;
use crate::render::{Primitive, Scene, PAD_PT};

/// 一套三维绘图样式。长度一律以**埃**为单位 —— 三维图是比例模型,这一点与
/// 二维的 [`Style`](crate::style::Style) 正相反(那边用磅)。
///
/// # 字段是公开的,改一项就是一次结构体更新
///
/// **并排比较不同样式时必须显式压到同一个比例尺** —— 四套样式各带各的默认值
/// (空间填充 24 磅每埃,其余 36),那是单独出一张图时的合理默认,并排摆就会让
/// 空间填充那格看着小一圈,而分子并没有变。库不会替调用方统一它。
///
/// ```
/// use omgkit_depict::three::Style3D;
///
/// let same_scale = Style3D {
///     scale_pt_per_a: 36.0,
///     ..Style3D::SPACE_FILLING
/// };
/// assert_eq!(same_scale.ball_vdw_frac, 1.0);
/// assert_eq!(Style3D::SPACE_FILLING.scale_pt_per_a, 24.0);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Style3D {
    /// 样式名。
    pub name: &'static str,
    /// 球半径 = 这个比例 × 范德华半径。`0.0` 表示不画球。
    pub ball_vdw_frac: f64,
    /// 键(圆柱)的半径,埃。`0.0` 表示不画键。
    pub stick_radius_a: f64,
    /// 多重键并排画时两根圆柱的中心距,埃。`0.0` 表示一律只画一根。
    pub multiple_bond_spacing_a: f64,
    /// 一埃画多少磅。
    pub scale_pt_per_a: f64,
}

impl Style3D {
    /// **空间填充**(CPK / spacefill)—— 球取满范德华半径,不画键。
    ///
    /// 看的是分子占多大地方。键在这个样式下一根也看不见,所以干脆不画。
    pub const SPACE_FILLING: Style3D = Style3D {
        name: "space-filling",
        ball_vdw_frac: 1.0,
        stick_radius_a: 0.0,
        multiple_bond_spacing_a: 0.0,
        scale_pt_per_a: 24.0,
    };

    /// **球棍**(ball-and-stick)—— 最常用的一种。球 23% vdW、杆 0.15 Å。
    pub const BALL_AND_STICK: Style3D = Style3D {
        name: "ball-and-stick",
        ball_vdw_frac: 0.23,
        stick_radius_a: 0.15,
        multiple_bond_spacing_a: 0.35,
        scale_pt_per_a: 36.0,
    };

    /// **棍状**(stick / licorice)—— 杆 0.30 Å,球与杆同粗所以接头是平滑的。
    ///
    /// 这一档**不并排画多重键**:杆半径 0.30 Å 时两根圆柱的中心距得大于
    /// 0.60 Å 才分得开,而那已经比一根 C–C 键的四成还长,画出来不像双键,
    /// 像两根键。棍状图本来看的就是骨架走向。
    pub const STICK: Style3D = Style3D {
        name: "stick",
        ball_vdw_frac: 0.0,
        stick_radius_a: 0.30,
        multiple_bond_spacing_a: 0.0,
        scale_pt_per_a: 36.0,
    };

    /// **线框**(wireframe)—— 只有细线,原子多的时候不糊成一团。
    pub const WIREFRAME: Style3D = Style3D {
        name: "wireframe",
        ball_vdw_frac: 0.0,
        stick_radius_a: 0.01,
        multiple_bond_spacing_a: 0.25,
        scale_pt_per_a: 36.0,
    };

    /// 本库内置的全部三维样式。
    pub const ALL: [Style3D; 4] = [
        Style3D::SPACE_FILLING,
        Style3D::BALL_AND_STICK,
        Style3D::STICK,
        Style3D::WIREFRAME,
    ];

    /// 棍状样式里球与杆同粗,所以接头看不出来。返回该画多大的球(埃)。
    ///
    /// `ball_vdw_frac` 为 0 而 `stick_radius_a` 非 0 时,球退成**圆头线段的
    /// 那个圆头** —— 半径正好是杆半径,不必单独画。
    #[must_use]
    fn ball_radius(&self, rvdw: f64) -> f64 {
        self.ball_vdw_frac * rvdw
    }
}

/// 元素表里 `rvdw` 为 0 时拿来顶的半径(埃)。
///
/// SMILES 的通配原子 `*`(记作原子序数 0)在元素表里 `rvdw = 0.0`,直接用会
/// 画出半径为 0 的球 —— 也就是**那个原子从图上消失了**,而图看着一点毛病没有。
/// 取碳的 1.7 Å 当替身:大小合理,而颜色仍是刺眼的 deeppink,读图的人看得出
/// 这里有个"不知道是什么"的原子。
const UNKNOWN_RVDW: f64 = 1.7;

/// 键切片之后每片允许的深度跨度(埃)。
///
/// 画家算法按图元的**一个**深度值排序,所以一个图元的深度跨度越小,排出来的
/// 次序越接近真的。取 0.25 Å:比球棍样式的球半径(碳 0.391 Å)小,于是一根
/// 斜穿球体的棍会在球的前后各留下正确的几片。
///
/// 再细下去只是把图元数量堆上去 —— 一根 1.5 Å 的键在最坏情况(正对观察者)
/// 也只切成 6 片,而那时它在投影上根本只有一个点。
pub const DEPTH_SLICE: f64 = 0.25;

/// 判定两个主惯量"相等"的相对容差 —— 见 [`View::degenerate`]。
///
/// # 这个数是量出来的,顺带纠正了一个想当然
///
/// 起初写的理由是"苯、金刚烷这类高对称分子两个特征值差在 1e-12 量级"。
/// **量下来不是**(相邻特征值的相对间隔,构象由 `omgkit-conf` 生成):
///
/// | 分子 | λ₀–λ₁ | λ₁–λ₂ | |
/// |---|---:|---:|---|
/// | 甲烷 `C` | 9.9e-16 | 8.5e-16 | 真简并 |
/// | 四氯化碳 | 5.6e-16 | 5.6e-16 | 真简并 |
/// | 四氟化碳 | 0 | 5.5e-16 | 真简并 |
/// | 氨 `N` | 1.4e-16 | 9.9e-1 | 真简并(两根垂直于 C₃ 轴的) |
/// | 乙炔 `C#C` | 1.0 | **0** | 真简并(线性分子) |
/// | **苯** | **6.1e-3** | 9.9e-1 | **不简并** |
/// | **金刚烷** | **3.3e-2** | 4.5e-2 | **不简并** |
/// | 环己烷 | 5.2e-2 | 5.9e-1 | 不简并 |
/// | 阿司匹林 | 4.8e-1 | 3.7e-1 | 不简并 |
///
/// 苯与金刚烷落在 1e-2 量级,不是 1e-12:**生成出来的构象不是理想的六边形/
/// 笼**,优化器留下的百分之几的不对称把两根轴分开了。它们的视角因此是**定的**
/// (同一组坐标永远同一个视角),只是"定"在一处物理上没有意义的不对称上。
///
/// 取 1e-6:恰好把上表分成两半 —— 对称性**强制**相等的那几个(全在 1e-16
/// 或恰好 0),与仅仅是"接近"的那几个(全在 1e-3 以上)。两边各留三个数量级
/// 的余量,而这条线的位置是上面那张表决定的,不是猜的。
pub const DEGENERATE_TOL: f64 = 1e-6;

/// 比较坐标时的量化精度,与 [`orient`](crate::orient) 同一个数。
///
/// 浮点直接比大小会让"两个图元谁先画"取决于最后一位,而那一位取决于运算次序
/// —— 同一个分子的不同写法就会给出次序不同的 SVG。
const QUANT: f64 = 1e6;

/// 视角:世界坐标 → 屏幕坐标。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct View {
    /// 旋转矩阵,**行主序**:第 `k` 行是屏幕第 `k` 根轴在世界坐标里的方向。
    ///
    /// 轴的含义:0 = 水平向右,1 = 竖直**向上**(不是 SVG 的向下),
    /// 2 = **指向观察者**(值越大越靠前)。
    ///
    /// 行列式恒为 +1 —— 见模块文档「朝向不许镜像」。
    pub rot: [[f64; 3]; 3],
    /// 旋转之前先减掉的中心(世界坐标)。
    pub centre: [f64; 3],
    /// **主轴不唯一** —— 两个主惯量在数值上相等,那两根轴在它们张成的平面里
    /// 可以任意转,选到哪一根取决于坐标最后几位的浮点噪声。
    ///
    /// 落进这一档的是**对称性强制简并**的那些:甲烷、四氯化碳、氨、乙炔
    /// (线性分子垂直于轴的两根)。它们的视角对坐标的微扰不稳 —— 换一份
    /// 构象、甚至同一份构象上加一点噪声,图就转了。
    ///
    /// **这不表示图是错的**,也不表示同一组坐标会画出两张图(那一条由
    /// [`canonical_view`] 的构造保证,与本标志无关)。它表示的是:这张图的
    /// 取向没有承载任何信息,别照着它去比两个分子的姿态。
    ///
    /// 苯、金刚烷这类"看着很对称"的分子**不在**这一档,理由见
    /// [`DEGENERATE_TOL`] 那张实测表。
    pub degenerate: bool,
}

impl View {
    /// 不转、不平移的视角。分子少于两个原子时用它 —— 那时没有主轴可言。
    pub const IDENTITY: View = View {
        rot: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        centre: [0.0, 0.0, 0.0],
        degenerate: false,
    };

    /// 把一个世界坐标点变到屏幕坐标 `(右, 上, 前)`。
    #[must_use]
    pub fn apply(&self, p: [f64; 3]) -> [f64; 3] {
        let d = [
            p[0] - self.centre[0],
            p[1] - self.centre[1],
            p[2] - self.centre[2],
        ];
        let mut out = [0.0; 3];
        for (o, r) in out.iter_mut().zip(&self.rot) {
            *o = r[0] * d[0] + r[1] * d[1] + r[2] * d[2];
        }
        out
    }

    /// 旋转矩阵的行列式。规范视角下恒为 +1。
    #[must_use]
    pub fn determinant(&self) -> f64 {
        let m = &self.rot;
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    }
}

/// 画不出三维图的原因。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Error3D {
    /// 坐标的个数与分子的原子数对不上。**这多半是拿错了构象** —— 继续画只会
    /// 得到一张张冠李戴的图,所以在这里拦住。
    CoordCount {
        /// 分子的原子数
        atoms: usize,
        /// 给进来的坐标数
        coords: usize,
    },
    /// 坐标里有非有限数(NaN / 无穷)。
    NotFinite {
        /// 第一个出问题的原子
        atom: usize,
    },
}

impl std::fmt::Display for Error3D {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error3D::CoordCount { atoms, coords } => {
                write!(f, "分子有 {atoms} 个原子,给进来的坐标有 {coords} 组")
            }
            Error3D::NotFinite { atom } => write!(f, "第 {atom} 个原子的坐标不是有限数"),
        }
    }
}

impl std::error::Error for Error3D {}

/// 一个原子在画布上的落点。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Placed {
    /// 圆心,**画布坐标(磅)**,与 [`Scene`] 里的图元同一个坐标系。
    pub at: Point2,
    /// 深度,**埃**,越大越靠前。画家算法排的就是它。
    ///
    /// 单位与 `at` 不同是故意的:`at` 是画布上的位置(随比例尺变),深度是分子
    /// 自己的量(不随比例尺变),混成一个单位会让"这两个原子差多深"这句话
    /// 没法问。
    pub depth: f64,
    /// 球的半径(磅)。样式不画球时是 `0.0`。
    pub radius: f64,
}

/// 一张三维分子图。
#[derive(Debug, Clone, PartialEq)]
pub struct Depiction3D {
    /// 图元,已经按深度从远到近排好。交给 [`svg`](crate::svg) 就能出图;
    /// 开了 `raster` feature 的话,`raster` 模块出 PNG / JPEG。
    pub scene: Scene,
    /// 每个原子落在画布的哪里,下标与传进来的分子一致。
    ///
    /// 在这里交出去,一是给要在图上加标注的调用方(图元里没有原子号,从
    /// [`Scene`] 反推是猜);二是给判据 —— 少了它,"球心是不是坐标的正交投影"
    /// 这类判据只能去 SVG 里按坐标反认球,认错一个就静默地比了别的东西。
    pub placed: Vec<Placed>,
    /// 用的是哪个视角。
    pub view: View,
    /// 样式名。
    pub style_name: &'static str,
}

impl Depiction3D {
    /// 有没有任何一处没画好。目前只有一档:视角退化。
    #[must_use]
    pub fn is_clean(&self) -> bool {
        !self.view.degenerate
    }
}

/// 给一组三维坐标画一张分子图。
///
/// `coords` 的单位是**埃**,下标与 `mol` 的原子一一对应 —— 通常直接来自
/// [`omgkit_conf::pipeline::conformer_for`],或从三维 molblock 读进来。
///
/// # 氢原子要自己准备好
///
/// 这里**不补氢**。三维图里氢是看得见的实体,补不补是调用方的决定:
/// `conformer_for` 会把立体中心缺的显式氢补进 `mol` 再给坐标,而从 molblock
/// 读进来的分子有什么就是什么。
///
/// # Errors
///
/// 坐标个数与原子数对不上、或坐标里有非有限数。
pub fn depict(
    mol: &MolBuilder,
    coords: &[[f64; 3]],
    style: &Style3D,
) -> Result<Depiction3D, Error3D> {
    if coords.len() != mol.num_atoms() {
        return Err(Error3D::CoordCount {
            atoms: mol.num_atoms(),
            coords: coords.len(),
        });
    }
    for (i, p) in coords.iter().enumerate() {
        if !p.iter().all(|x| x.is_finite()) {
            return Err(Error3D::NotFinite { atom: i });
        }
    }

    let ranks = crate::ranks_of(mol);
    let view = canonical_view(coords, &ranks);
    let screen: Vec<[f64; 3]> = coords.iter().map(|p| view.apply(*p)).collect();

    let mut items = build(mol, &screen, style);
    // **画家算法**:深度小的先画(远),大的后画(近)。平局按图元自身的几何
    // 打破 —— 拿生成次序打破的话,同一个分子换种写法会给出次序不同的 SVG。
    items.sort_by(|a, b| a.key.cmp(&b.key));

    let (scene, map) = to_scene(&items, style);
    let placed = mol
        .atoms()
        .iter()
        .zip(&screen)
        .map(|(a, p)| Placed {
            at: map(p[0], p[1]),
            depth: p[2],
            radius: style.ball_radius(rvdw_of(a.atomic_num)) * style.scale_pt_per_a,
        })
        .collect();

    Ok(Depiction3D {
        scene,
        placed,
        view,
        style_name: style.name,
    })
}

/// 算规范视角:主轴对齐。
///
/// # 为什么不是二维那套候选姿态
///
/// 见模块文档「与二维那条路的三处不同」。这里只补一条实现上的:轴的**符号**
/// 靠规范秩定 —— 特征向量的正负号是特征分解的自由度,不定死的话同一个分子的
/// 两种写法会得到互为 180° 的两张图。定法是"按秩排在最前、且在这根轴上确实
/// 偏离中心的那个原子,坐标为正"。
///
/// 第三根轴不单独定符号,取前两根的**叉积** —— 这样行列式必为 +1,镜像
/// 从构造上就不可能发生。
#[must_use]
pub fn canonical_view(coords: &[[f64; 3]], ranks: &[u32]) -> View {
    if coords.len() < 2 {
        return View::IDENTITY;
    }
    // **累加次序按坐标本身排,不按秩。** 浮点加法不满足结合律,按存储序累加会让
    // 两种编号差最后一位,而那一位足以让接近简并的两根轴对调。
    //
    // 先前这里排的是规范秩,而**秩在深层对称下不唯一**:苯的六个碳同属一个
    // 对称等价类,类内那点平局由规范 SMILES 的输出次序打破,而六种起笔给出的串
    // 一模一样 —— 于是最终还是落回存储序。实测苯换个编号画出来上下翻了个个儿。
    //
    // 按坐标排是**逐项不变**的:重新编号动不了坐标的多重集。平局(坐标逐位
    // 相同的两个原子)再落到秩和下标上,那时两者相加的结果本来就一样。
    let mut order: Vec<usize> = (0..coords.len()).collect();
    order.sort_by(|&i, &j| {
        coords[i][0]
            .total_cmp(&coords[j][0])
            .then(coords[i][1].total_cmp(&coords[j][1]))
            .then(coords[i][2].total_cmp(&coords[j][2]))
            .then(ranks[i].cmp(&ranks[j]))
            .then(i.cmp(&j))
    });

    let n = coords.len() as f64;
    let mut centre = [0.0f64; 3];
    for &i in &order {
        for k in 0..3 {
            centre[k] += coords[i][k];
        }
    }
    for c in &mut centre {
        *c /= n;
    }

    let mut cov = [0.0f64; 9];
    for &i in &order {
        let d = [
            coords[i][0] - centre[0],
            coords[i][1] - centre[1],
            coords[i][2] - centre[2],
        ];
        for a in 0..3 {
            for b in 0..3 {
                cov[a * 3 + b] += d[a] * d[b];
            }
        }
    }

    let Ok(eig) = symmetric_eigen(&cov, 3) else {
        // 分解失败(输入病态)。**不假装有个视角** —— 退回不转,并报退化。
        return View {
            degenerate: true,
            ..View::IDENTITY
        };
    };

    // 特征值降序:最大方差那根放水平,最小的那根指向观察者。
    let mut axes = [[0.0f64; 3]; 3];
    for (k, ax) in axes.iter_mut().enumerate() {
        ax.copy_from_slice(eig.vector(k));
    }

    // **前两根轴的正负号**:特征分解定不了它(`v` 与 `−v` 都是特征向量)。
    // 四种组合各投影一遍,取"投影坐标排序之后字典序最小"的那一种。
    //
    // 为什么是排序过的多重集,而不是"按秩排在最前那个原子的坐标为正":后者
    // 要秩唯一,而秩在苯这类分子上不唯一(见上)。**多重集在重新编号下逐项
    // 不变**,拿它当键就与编号彻底无关了。二维的 [`orient`](crate::orient)
    // 在 24 个候选姿态里挑用的是同一招。
    //
    // 第三根轴不参与挑选,一律取前两根的**叉积** —— 行列式因此必为 +1,
    // 镜像从构造上就不可能发生。
    let base = axes;
    let mut best: Option<Candidate> = None;
    for s0 in [1.0f64, -1.0] {
        for s1 in [1.0f64, -1.0] {
            let a0 = [base[0][0] * s0, base[0][1] * s0, base[0][2] * s0];
            let a1 = [base[1][0] * s1, base[1][1] * s1, base[1][2] * s1];
            let a2 = [
                a0[1] * a1[2] - a0[2] * a1[1],
                a0[2] * a1[0] - a0[0] * a1[2],
                a0[0] * a1[1] - a0[1] * a1[0],
            ];
            let cand = [a0, a1, a2];
            let mut key: Vec<[i64; 3]> = coords
                .iter()
                .map(|p| {
                    let d = [p[0] - centre[0], p[1] - centre[1], p[2] - centre[2]];
                    let mut out = [0i64; 3];
                    for (o, ax) in out.iter_mut().zip(&cand) {
                        *o = q(ax[0] * d[0] + ax[1] * d[1] + ax[2] * d[2]);
                    }
                    out
                })
                .collect();
            key.sort_unstable();
            // MSRV 1.75:`Option::is_none_or` 要 1.82,这里只能用 `map_or`
            #[allow(clippy::unnecessary_map_or)]
            if best.as_ref().map_or(true, |(b, _)| key < *b) {
                best = Some((key, cand));
            }
        }
    }
    let axes = best.expect("四种符号组合总有一个").1;

    // 简并:相邻的两个特征值相等时,那两根轴在它们张成的平面里可以任意转。
    let scale = eig.values[0].abs().max(f64::MIN_POSITIVE);
    let degenerate =
        (0..2).any(|k| (eig.values[k] - eig.values[k + 1]).abs() <= DEGENERATE_TOL * scale);

    View {
        rot: axes,
        centre,
        degenerate,
    }
}

/// 挑视角时的一个候选:(投影坐标的量化多重集, 三根轴)。
type Candidate = (Vec<[i64; 3]>, [[f64; 3]; 3]);

/// 排好序之前的一个图元:几何 + 排序键。
struct Sortable {
    prim: Primitive,
    key: Key,
}

/// 画家算法的排序键。**完全由图元自身的几何决定**,不含生成次序 ——
/// 否则同一分子的两种写法会给出次序不同的 SVG。
///
/// 第一项是深度(量化过的),小的先画。其余各项只用来打破平局。
#[derive(PartialEq, Eq, PartialOrd, Ord)]
struct Key(i64, u8, [i64; 4], [u8; 3]);

fn q(x: f64) -> i64 {
    // 量化到 QUANT 分之一。超出 i64 的坐标在这里饱和 —— 那种输入早在
    // `depict` 的有限性检查之后就只剩天文数字了,饱和好过 `as` 的未定义感。
    let v = (x * QUANT).round();
    if v > i64::MAX as f64 {
        i64::MAX
    } else if v < i64::MIN as f64 {
        i64::MIN
    } else {
        v as i64
    }
}

/// 元素的范德华半径,表里没有(或为 0)时给 [`UNKNOWN_RVDW`]。
fn rvdw_of(atomic_num: u8) -> f64 {
    let r = omgkit_core::element::by_atomic_num(atomic_num).map_or(0.0, |e| f64::from(e.rvdw));
    if r > 0.0 {
        r
    } else {
        UNKNOWN_RVDW
    }
}

/// 造出全部图元(还没排序、还是埃、y 还朝上)。
fn build(mol: &MolBuilder, screen: &[[f64; 3]], style: &Style3D) -> Vec<Sortable> {
    let mut out = Vec::with_capacity(mol.num_atoms() + mol.num_bonds() * 4);

    if style.stick_radius_a > 0.0 {
        for b in mol.bonds() {
            push_bond(
                &mut out,
                mol,
                screen,
                style,
                b.begin as usize,
                b.end as usize,
                b.order,
            );
        }
    }

    for (i, a) in mol.atoms().iter().enumerate() {
        let r = style.ball_radius(rvdw_of(a.atomic_num));
        if r <= 0.0 {
            continue;
        }
        let p = screen[i];
        let at = Point2::new(p[0], p[1]);
        let color = cpk(a.atomic_num);
        out.push(Sortable {
            prim: Primitive::Ball { at, r, color },
            key: Key(q(p[2]), 0, [q(p[0]), q(p[1]), q(r), 0], color),
        });
    }
    out
}

/// 一根键该画几根并排的圆柱。
///
/// 芳香键画**一根**。三维图不做 Kekulé 化 —— 苯环上六根键的键级是一样的,
/// 硬要给其中三根画两条线就得先挑出一套 Kekulé 结构,而那套结构是任取的,
/// 同一个分子换种写法会挑到另一套,于是双键画在了另外三条边上。
/// PyMOL 默认(`valence 0`)也是一根。
fn cylinders(order: BondOrder) -> usize {
    match order {
        BondOrder::Double => 2,
        BondOrder::Triple | BondOrder::Quadruple => 3,
        _ => 1,
    }
}

#[allow(clippy::too_many_arguments)]
fn push_bond(
    out: &mut Vec<Sortable>,
    mol: &MolBuilder,
    screen: &[[f64; 3]],
    style: &Style3D,
    ia: usize,
    ib: usize,
    order: BondOrder,
) {
    let (pa, pb) = (screen[ia], screen[ib]);
    let d = [pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]];

    // 并排的方向:垂直于键在**屏幕平面上的投影**,且躺在屏幕平面里 ——
    // Jmol 的说法是 "twin sticks contained in the viewer's plane"。
    let flat = (d[0] * d[0] + d[1] * d[1]).sqrt();
    let n_cyl = if style.multiple_bond_spacing_a > 0.0 && flat > f64::EPSILON {
        cylinders(order)
    } else {
        // 键正对着观察者时它在屏幕上只有一个点,并排画没有意义,也没有方向可选。
        1
    };
    let perp = if flat > f64::EPSILON {
        [d[1] / flat, -d[0] / flat]
    } else {
        [0.0, 0.0]
    };

    let width = style.stick_radius_a * 2.0;
    let ca = cpk(mol.atoms()[ia].atomic_num);
    let cb = cpk(mol.atoms()[ib].atomic_num);
    let ra = style.ball_radius(rvdw_of(mol.atoms()[ia].atomic_num));
    let rb = style.ball_radius(rvdw_of(mol.atoms()[ib].atomic_num));
    let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();

    for k in 0..n_cyl {
        #[allow(clippy::cast_precision_loss)]
        let off = (k as f64 - (n_cyl as f64 - 1.0) / 2.0) * style.multiple_bond_spacing_a;
        let sa = [pa[0] + perp[0] * off, pa[1] + perp[1] * off, pa[2]];
        let sb = [pb[0] + perp[0] * off, pb[1] + perp[1] * off, pb[2]];
        let mid = [
            (sa[0] + sb[0]) / 2.0,
            (sa[1] + sb[1]) / 2.0,
            (sa[2] + sb[2]) / 2.0,
        ];
        push_half(out, sa, mid, trim(ra, off, len), width, ca);
        push_half(out, sb, mid, trim(rb, off, len), width, cb);
    }
}

/// **半根键有多长一段埋在原子球里** —— 那一段不能画。
///
/// 画了会怎样:靠近原子那几片的深度与球心相当,而键只要朝观察者偏一点,
/// 球外的第一片就比球心深,于是排在球**后面**画……不,是排在球**之后**画,
/// 一个圆头端帽就压在球面上,球上多出一道月牙。实测的第一版正是这样,
/// 每个碳球里都有一道浅灰的弧 —— 图看着"有点脏",而没有任何一处报错。
///
/// 几何是精确的,不是估的:圆柱的轴离球心 `off`(多重键并排时),
/// 它穿出半径 `r` 的球面处离球心 `sqrt(r² − off²)`。`|off| ≥ r` 时圆柱整根
/// 在球外,一点都不用截。
fn trim(r: f64, off: f64, len: f64) -> f64 {
    let inside = (r * r - off * off).max(0.0).sqrt();
    // 截过头就等于整半根都在球里 —— 让调用方一片也别发。
    inside.min(len / 2.0)
}

/// 半根键:从原子那一端到键中点,按深度切成若干片。`skip` 是起点那一头要
/// 让开多长(埋在原子球里的那一段,见 [`trim`])。
///
/// 切片的理由见模块文档。片与片之间**不留缝**:圆头端帽自然接上,而端点是
/// 同一个数算出来的,不会差半个像素。
fn push_half(
    out: &mut Vec<Sortable>,
    from: [f64; 3],
    to: [f64; 3],
    skip: f64,
    width: f64,
    color: [u8; 3],
) {
    let half = dist(from, to);
    if skip >= half {
        return; // 整半根都埋在球里
    }
    let from = if skip > 0.0 {
        lerp(from, to, skip / half)
    } else {
        from
    };
    let dz = (to[2] - from[2]).abs();
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let n = ((dz / DEPTH_SLICE).ceil() as usize).max(1);
    #[allow(clippy::cast_precision_loss)]
    let nf = n as f64;
    for s in 0..n {
        #[allow(clippy::cast_precision_loss)]
        let (t0, t1) = (s as f64 / nf, (s + 1) as f64 / nf);
        let p0 = lerp(from, to, t0);
        let p1 = lerp(from, to, t1);
        let depth = (p0[2] + p1[2]) / 2.0;
        out.push(Sortable {
            prim: Primitive::Stick {
                from: Point2::new(p0[0], p0[1]),
                to: Point2::new(p1[0], p1[1]),
                width,
                color,
            },
            key: Key(q(depth), 1, [q(p0[0]), q(p0[1]), q(p1[0]), q(p1[1])], color),
        });
    }
}

fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    let (x, y, z) = (b[0] - a[0], b[1] - a[1], b[2] - a[2]);
    (x * x + y * y + z * z).sqrt()
}

fn lerp(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    [
        a[0] + (b[0] - a[0]) * t,
        a[1] + (b[1] - a[1]) * t,
        a[2] + (b[2] - a[2]) * t,
    ]
}

/// 埃 → 磅,y 轴翻成 SVG 的朝下,再算画布。
///
/// 一并把那个映射交回去 —— [`Depiction3D::placed`] 要用**同一个**映射,
/// 各算一遍的话原子落点与图元会差半个像素,而那种偏差没有判据看得见。
fn to_scene(items: &[Sortable], style: &Style3D) -> (Scene, impl Fn(f64, f64) -> Point2) {
    let (mut lo, mut hi) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
    let mut note = |x: f64, y: f64, r: f64| {
        lo[0] = lo[0].min(x - r);
        lo[1] = lo[1].min(y - r);
        hi[0] = hi[0].max(x + r);
        hi[1] = hi[1].max(y + r);
    };
    for it in items {
        match &it.prim {
            Primitive::Ball { at, r, .. } => note(at.x, at.y, *r),
            Primitive::Stick {
                from, to, width, ..
            } => {
                note(from.x, from.y, width / 2.0);
                note(to.x, to.y, width / 2.0);
            }
            _ => {}
        }
    }
    if !lo[0].is_finite() {
        // 一个图元都没有(空分子,或者样式把球和杆都关了)
        lo = [0.0, 0.0];
        hi = [0.0, 0.0];
    }

    let s = style.scale_pt_per_a;
    let map = move |x: f64, y: f64| Point2::new((x - lo[0]) * s + PAD_PT, (hi[1] - y) * s + PAD_PT);
    let out = items
        .iter()
        .map(|it| match &it.prim {
            Primitive::Ball { at, r, color } => Primitive::Ball {
                at: map(at.x, at.y),
                r: r * s,
                color: *color,
            },
            Primitive::Stick {
                from,
                to,
                width,
                color,
            } => Primitive::Stick {
                from: map(from.x, from.y),
                to: map(to.x, to.y),
                width: width * s,
                color: *color,
            },
            other => other.clone(),
        })
        .collect();

    (
        Scene {
            items: out,
            width: (hi[0] - lo[0]) * s + PAD_PT * 2.0,
            height: (hi[1] - lo[1]) * s + PAD_PT * 2.0,
        },
        map,
    )
}

#[cfg(test)]
mod tests {
    use super::{
        canonical_view, cylinders, depict, dist, lerp, rvdw_of, trim, Error3D, Placed, Style3D,
        DEPTH_SLICE, UNKNOWN_RVDW,
    };
    use crate::render::Primitive;
    use crate::style::Style;
    use crate::svg::to_svg;
    use omgkit_core::MolBuilder;

    /// 解析 + 生成构象。三维图要的是坐标,而坐标只有这一条路来。
    fn prep(smi: &str) -> (MolBuilder, Vec<[f64; 3]>) {
        let mut m = omgkit_io::smiles::parse(smi).expect("测试用的 SMILES 该能解析");
        let c = omgkit_conf::pipeline::conformer_for(&mut m).expect("测试用的分子该能出构象");
        (m, c.coords)
    }

    /// 判据用的一小批分子:链、芳环、稠环、杂原子、多重键、离子、金刚烷。
    const CORPUS: [&str; 10] = [
        "CCO",
        "c1ccccc1",
        "CC(=O)Oc1ccccc1C(=O)O",
        "C1CC2CCC1CC2",
        "C1C2CC3CC1CC(C2)C3",
        "N#Cc1ccncc1",
        "CC(C)(C)S(=O)(=O)N",
        "[Na+].[Cl-]",
        "C/C=C/C=C/C",
        "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
    ];

    /// **四套样式的数逐项对上 Jmol 自己文档里的 standard rendering styles。**
    ///
    /// 这一条是转录判据,不是设计判据 —— 数是抄来的,抄错了没有任何一处会报错,
    /// 只会画出一张"看着像但不是那个样式"的图。
    #[test]
    fn 样式的半径与_jmol_文档一致() {
        // spacefill 100%
        assert_eq!(Style3D::SPACE_FILLING.ball_vdw_frac, 1.0);
        assert_eq!(Style3D::SPACE_FILLING.stick_radius_a, 0.0);
        // wireframe 0.15; spacefill 23%
        assert_eq!(Style3D::BALL_AND_STICK.ball_vdw_frac, 0.23);
        assert_eq!(Style3D::BALL_AND_STICK.stick_radius_a, 0.15);
        // wireframe 0.3; spacefill off
        assert_eq!(Style3D::STICK.ball_vdw_frac, 0.0);
        assert_eq!(Style3D::STICK.stick_radius_a, 0.30);
        // wireframe 0.01; spacefill 0
        assert_eq!(Style3D::WIREFRAME.ball_vdw_frac, 0.0);
        assert_eq!(Style3D::WIREFRAME.stick_radius_a, 0.01);
    }

    /// **视角矩阵必须是真旋转** —— 正交、且行列式为 +1。
    ///
    /// 行列式为 −1 的话整个分子被镜像了,而三维图的构型**就是坐标本身**:
    /// 画出来的每一个手性中心都反了,图上没有任何一处看得出来。二维那条路
    /// 允许镜像(楔形是镜像之后重新指派的),三维不允许 —— 这是两条路最要紧的
    /// 一处不同。
    #[test]
    fn 视角是真旋转不是镜像() {
        for smi in CORPUS {
            let (m, c) = prep(smi);
            let v = canonical_view(&c, &crate::ranks_of(&m));
            assert!(
                (v.determinant() - 1.0).abs() < 1e-9,
                "{smi}:行列式 {},不是 +1 —— 分子被镜像了",
                v.determinant()
            );
            for i in 0..3 {
                for j in 0..3 {
                    let dot: f64 = (0..3).map(|k| v.rot[i][k] * v.rot[j][k]).sum();
                    let want = f64::from(u8::from(i == j));
                    assert!(
                        (dot - want).abs() < 1e-9,
                        "{smi}:第 {i} 行与第 {j} 行的内积是 {dot},该是 {want}"
                    );
                }
            }
        }
    }

    /// 把原子重新编号,画出来**逐字节相同**。
    ///
    /// 这是本仓的头号契约在三维上的样子。视角靠特征分解定,而特征分解吃的是
    /// 二阶矩 —— 按存储序累加的话,两种编号会差最后一位,接近简并的两根轴
    /// 就此对调,整张图转了 90°。质心与二阶矩因此都按规范秩累加。
    ///
    /// 取的分子都**没有立体中心**:重新编号会改变每个原子的邻居次序,而手性
    /// 标记是相对邻居次序存的,一起改才对得上 —— 那是另一件事,不该混进这条
    /// 判据。带立体的分子由全语料判官 `harness/check_depict3d.py` 覆盖。
    #[test]
    fn 换个原子编号画出来逐字节相同() {
        for smi in [
            "CC(=O)Oc1ccccc1C(=O)O",
            "c1ccccc1",
            "C1CC2CCC1CC2",
            "N#Cc1ccncc1",
            "OCCOCCO",
        ] {
            let (m, c) = prep(smi);
            let (m2, c2) = renumbered(&m, &c);
            assert_eq!(
                omgkit_io::canon::canonical_smiles(&m).smiles,
                omgkit_io::canon::canonical_smiles(&m2).smiles,
                "{smi}:重新编号之后不是同一个分子,判据本身坏了"
            );
            for style in &Style3D::ALL {
                let a = to_svg(&depict(&m, &c, style).unwrap().scene, &Style::ACS_1996);
                let b = to_svg(&depict(&m2, &c2, style).unwrap().scene, &Style::ACS_1996);
                assert_eq!(a, b, "{smi} / {}:重新编号之后画出来不一样", style.name);
            }
        }
    }

    /// 原子倒序重排,坐标跟着走。
    fn renumbered(mol: &MolBuilder, coords: &[[f64; 3]]) -> (MolBuilder, Vec<[f64; 3]>) {
        let n = mol.num_atoms();
        let mut out = MolBuilder::with_capacity(n, mol.num_bonds());
        for a in mol.atoms().iter().rev() {
            out.add_atom_data(*a);
        }
        #[allow(clippy::cast_possible_truncation)]
        let map = |old: u32| n as u32 - 1 - old;
        for b in mol.bonds() {
            let mut nb = *b;
            nb.begin = map(b.begin);
            nb.end = map(b.end);
            out.add_bond_data(nb).unwrap();
        }
        let mut c = vec![[0.0; 3]; n];
        for (i, p) in coords.iter().enumerate() {
            c[n - 1 - i] = *p;
        }
        (out, c)
    }

    /// 图元一个都不许出画布。
    #[test]
    fn 图元不出画布() {
        for smi in CORPUS {
            let (m, c) = prep(smi);
            for style in &Style3D::ALL {
                let d = depict(&m, &c, style).unwrap();
                let s = &d.scene;
                for it in &s.items {
                    let pts: Vec<(crate::geom::Point2, f64)> = match it {
                        Primitive::Ball { at, r, .. } => vec![(*at, *r)],
                        Primitive::Stick {
                            from, to, width, ..
                        } => vec![(*from, *width / 2.0), (*to, *width / 2.0)],
                        _ => panic!("三维图里出现了二维图元"),
                    };
                    for (p, r) in pts {
                        assert!(
                            p.x - r >= -0.01 && p.x + r <= s.width + 0.01,
                            "{smi} / {}:图元 x={:.2}(±{r:.2})出了画布宽 {:.2}",
                            style.name,
                            p.x,
                            s.width
                        );
                        assert!(
                            p.y - r >= -0.01 && p.y + r <= s.height + 0.01,
                            "{smi} / {}:图元 y={:.2}(±{r:.2})出了画布高 {:.2}",
                            style.name,
                            p.y,
                            s.height
                        );
                    }
                }
            }
        }
    }

    /// **投影上重叠的两个球,深的先画。** 画家算法唯一要守的事。
    ///
    /// 反过来的话,后面的球盖在前面的球上 —— 图上读到的立体关系是反的,
    /// 而线条本身一点毛病没有。
    #[test]
    fn 重叠的球近的后画() {
        for smi in CORPUS {
            let (m, c) = prep(smi);
            for style in [&Style3D::SPACE_FILLING, &Style3D::BALL_AND_STICK] {
                let d = depict(&m, &c, style).unwrap();
                let order = ball_order(&d.scene, &d.placed);
                for i in 0..d.placed.len() {
                    for j in 0..i {
                        let (a, b) = (&d.placed[i], &d.placed[j]);
                        if a.at.dist(b.at) >= a.radius + b.radius {
                            continue; // 投影上不重叠,谁先画都行
                        }
                        let (near, far) = if a.depth > b.depth { (i, j) } else { (j, i) };
                        if (a.depth - b.depth).abs() < 1e-12 {
                            continue; // 一样深,次序由几何平局规则定
                        }
                        assert!(
                            order[near] > order[far],
                            "{smi} / {}:原子 {near}(深度 {:.3})比 {far}(深度 {:.3})靠前,\
                             却画在了前面",
                            style.name,
                            d.placed[near].depth,
                            d.placed[far].depth
                        );
                    }
                }
            }
        }
    }

    /// 每个原子的球在 `scene.items` 里排第几。按落点认球。
    fn ball_order(scene: &crate::render::Scene, placed: &[Placed]) -> Vec<usize> {
        placed
            .iter()
            .map(|p| {
                scene
                    .items
                    .iter()
                    .position(|it| {
                        matches!(it, Primitive::Ball { at, .. }
                                 if (at.x - p.at.x).abs() < 1e-9 && (at.y - p.at.y).abs() < 1e-9)
                    })
                    .expect("每个原子都该有一个球")
            })
            .collect()
    }

    /// **每根键的半棍确实从球面上起步,而不是从球心。**
    ///
    /// 埋在球里的那一段不截掉的话,朝观察者偏的键会在球面上留下一道月牙 ——
    /// 实测第一版每个碳球里都有一道浅灰的弧,图看着"有点脏",而没有任何一处
    /// 报错。
    ///
    /// 先前这条判据写成"任何棍的端点都不许落进任何球的投影里",那是**错的**:
    /// 一根跟这个原子无关的键从它前面或后面横过去,投影上本来就会压进那个圆。
    /// 判据要问的是每根键**自己**那一头,所以按键逐根构造期望的起点再去图里找。
    #[test]
    fn 半棍从球面起步而不是从球心() {
        let style = &Style3D::BALL_AND_STICK;
        for smi in CORPUS {
            let (m, c) = prep(smi);
            let d = depict(&m, &c, style).unwrap();
            let s = style.scale_pt_per_a;
            for b in m.bonds() {
                if cylinders(b.order) != 1 {
                    continue; // 并排的多重键起点各自偏开,另算,不在这条判据里
                }
                for (near, far) in [(b.begin, b.end), (b.end, b.begin)] {
                    let p = d.view.apply(c[near as usize]);
                    let qv = d.view.apply(c[far as usize]);
                    // **期望值不许调 `trim`** —— 那是被测的东西,两边一起改就
                    // 永远打不红(实测:把 `trim` 整个乘 0,这条判据一声不吭)。
                    // 单根圆柱走的是轴过球心那一支,交点离球心正好一个半径。
                    let r = style.ball_vdw_frac
                        * f64::from(
                            omgkit_core::element::by_atomic_num(
                                m.atoms()[near as usize].atomic_num,
                            )
                            .unwrap()
                            .rvdw,
                        );
                    let len = dist(p, qv);
                    let start = lerp(p, qv, r.min(len / 2.0) / len);
                    // 画布坐标:与 `to_scene` 的映射同构 —— 平移由原子落点给出,
                    // 尺度是比例尺,y 翻向。
                    let at = d.placed[near as usize].at;
                    let want = crate::geom::Point2::new(
                        at.x + (start[0] - p[0]) * s,
                        at.y - (start[1] - p[1]) * s,
                    );
                    assert!(
                        d.scene.items.iter().any(|it| matches!(
                            it,
                            Primitive::Stick { from, .. } if from.dist(want) < 0.01
                        )),
                        "{smi}:键 {near}-{far} 该从 ({:.2},{:.2}) 起步,图里没有这一段",
                        want.x,
                        want.y
                    );
                }
            }
        }
    }

    /// **把整个分子刚体旋转一下,画出来还是同一张图。**
    ///
    /// 视角的活儿就是把输入摆正,所以输入本来怎么摆不该有影响 —— 同一份构象
    /// 从 SMILES 生成、还是从别人写的 molblock 读进来(那是另一个坐标系),
    /// 该给出同一张图。
    ///
    /// 这一条同时是**符号选取规则的判据**。特征向量的正负号是特征分解的自由度,
    /// 分解器换个扫描次序就可能整根反过来;不定一条几何上的规则,文档里那批
    /// 配图会在某次升级之后集体上下颠倒。规则是"投影坐标的多重集字典序最小",
    /// 而这条判据正是它守住的东西 —— 实测把规则换成"永远取第一种组合",
    /// 这条判据当场红。
    ///
    /// 比的是原子落点而不是 SVG 字节:旋转会经过一遍三角函数,末位噪声必然有,
    /// 而那点噪声在 0.01 磅(300 dpi 下 0.04 像素)之下看不出来。
    #[test]
    fn 输入怎么摆都画出同一张图() {
        // 一个随手取的真旋转:绕 (1,1,1)/√3 转 1 弧度。
        let (c1, s1) = (1.0f64.cos(), 1.0f64.sin());
        let u = [
            std::f64::consts::FRAC_1_SQRT_2 * 0.816_496_580_927_726,
            0.577_350_269_189_626,
            0.577_350_269_189_626,
        ];
        let u = {
            let n = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
            [u[0] / n, u[1] / n, u[2] / n]
        };
        let rot = |p: [f64; 3]| {
            // 罗德里格斯公式
            let dot = u[0] * p[0] + u[1] * p[1] + u[2] * p[2];
            let cross = [
                u[1] * p[2] - u[2] * p[1],
                u[2] * p[0] - u[0] * p[2],
                u[0] * p[1] - u[1] * p[0],
            ];
            let mut out = [0.0; 3];
            for k in 0..3 {
                out[k] = p[k] * c1 + cross[k] * s1 + u[k] * dot * (1.0 - c1) + 3.7;
            }
            out
        };
        for smi in CORPUS {
            let (m, c) = prep(smi);
            let moved: Vec<[f64; 3]> = c.iter().map(|p| rot(*p)).collect();
            for style in &Style3D::ALL {
                let a = depict(&m, &c, style).unwrap();
                let b = depict(&m, &moved, style).unwrap();
                assert!(
                    (a.scene.width - b.scene.width).abs() < 0.01
                        && (a.scene.height - b.scene.height).abs() < 0.01,
                    "{smi} / {}:换个摆法画布尺寸都变了({:.2}×{:.2} vs {:.2}×{:.2})",
                    style.name,
                    a.scene.width,
                    a.scene.height,
                    b.scene.width,
                    b.scene.height
                );
                for (i, (x, y)) in a.placed.iter().zip(&b.placed).enumerate() {
                    assert!(
                        x.at.dist(y.at) < 0.01 && (x.depth - y.depth).abs() < 0.001,
                        "{smi} / {}:原子 {i} 换个摆法落到了别处 \
                         (({:.3},{:.3}) vs ({:.3},{:.3}))",
                        style.name,
                        x.at.x,
                        x.at.y,
                        y.at.x,
                        y.at.y
                    );
                }
            }
        }
    }

    /// **视角只由坐标定,与样式无关。**
    ///
    /// 文档里那张"四套样式并排"的图靠的正是这一条 —— 说的是"同一组坐标、
    /// 同一个视角,差别只在怎么画"。视角要是跟着样式变,那句话就是假的,
    /// 而四张图并排看还是像那么回事。
    #[test]
    fn 视角与样式无关() {
        for smi in CORPUS {
            let (m, c) = prep(smi);
            let first = depict(&m, &c, &Style3D::ALL[0]).unwrap().view;
            for style in &Style3D::ALL[1..] {
                let v = depict(&m, &c, style).unwrap().view;
                assert_eq!(
                    v.rot, first.rot,
                    "{smi}:{} 的视角与别的样式不同",
                    style.name
                );
                assert_eq!(v.centre, first.centre, "{smi}:{} 的中心不同", style.name);
            }
        }
    }

    /// 样式说不画球就一个球都没有,说不画棍就一根棍都没有。
    #[test]
    fn 样式关掉的东西一个都不发() {
        let (m, c) = prep("CC(=O)Oc1ccccc1C(=O)O");
        let sf = depict(&m, &c, &Style3D::SPACE_FILLING).unwrap();
        assert!(
            sf.scene
                .items
                .iter()
                .all(|it| matches!(it, Primitive::Ball { .. })),
            "空间填充不该有棍"
        );
        assert_eq!(
            sf.scene.items.len(),
            m.num_atoms(),
            "空间填充该一个原子一个球"
        );
        for style in [&Style3D::STICK, &Style3D::WIREFRAME] {
            let d = depict(&m, &c, style).unwrap();
            assert!(
                d.scene
                    .items
                    .iter()
                    .all(|it| matches!(it, Primitive::Stick { .. })),
                "{} 不该有球",
                style.name
            );
        }
    }

    /// 三维图里没有文字图元,所以**拿哪套二维规范去渲染都一样**。
    ///
    /// 这一条要判死:`svg::to_svg` 收一个 `&Style`,而三维那条路只能随便给一个。
    /// 哪天渲染开始读规范里的别的字段,这条判据当场红,而不是让调用方在某个
    /// 分子上撞见两张不同的图。
    #[test]
    fn 三维图与二维规范无关() {
        let (m, c) = prep("CC(=O)Oc1ccccc1C(=O)O");
        for style in &Style3D::ALL {
            let d = depict(&m, &c, style).unwrap();
            assert_eq!(
                to_svg(&d.scene, &Style::ACS_1996),
                to_svg(&d.scene, &Style::CHEMDRAW_DEFAULT),
                "{}:换套二维规范画出来不一样了",
                style.name
            );
        }
    }

    /// 对称性高的分子必须报视角退化,而不对称的必须不报。
    ///
    /// **两头都要判**:只判"苯报退化"的话,把这个标志写死成 `true` 也照样绿。
    #[test]
    fn 对称强制简并的分子报视角退化() {
        // 对称性**强制**两个主惯量相等的:正四面体、C₃ᵥ、线性
        for smi in ["C", "C(Cl)(Cl)(Cl)Cl", "FC(F)(F)F", "N", "C#C"] {
            let (m, c) = prep(smi);
            let v = canonical_view(&c, &crate::ranks_of(&m));
            assert!(v.degenerate, "{smi} 的主轴不唯一,该报退化");
        }
        // 反着也要判:标志写死成 true 的话下面这几条会红。
        //
        // **苯与金刚烷故意放在这一侧** —— 直觉会把它们归到"高对称"里,而
        // 实测它们相邻特征值的相对间隔是 6.1e-3 与 3.3e-2(见 DEGENERATE_TOL
        // 那张表),离简并差三个数量级。判据照实测写,不照直觉写。
        for smi in [
            "CC(=O)Oc1ccccc1C(=O)O",
            "CCCCCCO",
            "N#Cc1ccncc1",
            "c1ccccc1",
            "C1C2CC3CC1CC(C2)C3",
            "O",
        ] {
            let (m, c) = prep(smi);
            let v = canonical_view(&c, &crate::ranks_of(&m));
            assert!(!v.degenerate, "{smi} 的主轴是唯一的,不该报退化");
        }
    }

    /// 坐标个数对不上、或者有 NaN,当场报错 —— 不许画出一张张冠李戴的图。
    #[test]
    fn 坏输入报错而不是画出来() {
        let (m, c) = prep("CCO");
        assert_eq!(
            depict(&m, &c[..2], &Style3D::BALL_AND_STICK),
            Err(Error3D::CoordCount {
                atoms: m.num_atoms(),
                coords: 2
            })
        );
        let mut bad = c.clone();
        bad[1][2] = f64::NAN;
        assert_eq!(
            depict(&m, &bad, &Style3D::BALL_AND_STICK),
            Err(Error3D::NotFinite { atom: 1 })
        );
        bad[1][2] = f64::INFINITY;
        assert_eq!(
            depict(&m, &bad, &Style3D::BALL_AND_STICK),
            Err(Error3D::NotFinite { atom: 1 })
        );
    }

    /// 通配原子 `*` 在元素表里 `rvdw = 0`,不顶一个半径就会从图上消失。
    #[test]
    fn 没有范德华半径的原子仍然画得出来() {
        assert_eq!(rvdw_of(0), UNKNOWN_RVDW);
        assert!((rvdw_of(6) - 1.7).abs() < 1e-6, "碳该是 1.7");
        let (m, c) = prep("*CC");
        let d = depict(&m, &c, &Style3D::SPACE_FILLING).unwrap();
        assert_eq!(d.scene.items.len(), m.num_atoms(), "通配原子也该有一个球");
        assert!(d.placed[0].radius > 0.0, "通配原子的球半径不许是 0");
    }

    /// 截长的几何是精确的,不是估的。
    #[test]
    fn 截长按球与圆柱轴的交点算() {
        // 轴过球心:截到球面,正好是半径
        assert!((trim(0.4, 0.0, 1.5) - 0.4).abs() < 1e-12);
        // 轴偏离球心 0.3:交点在 sqrt(0.4² − 0.3²) = 0.264...
        assert!((trim(0.4, 0.3, 1.5) - (0.16f64 - 0.09).sqrt()).abs() < 1e-12);
        // 轴整根在球外:一点都不截
        assert_eq!(trim(0.4, 0.5, 1.5), 0.0);
        // 球比半根键还大:最多截到中点,不许截成负的
        assert_eq!(trim(2.0, 0.0, 1.5), 0.75);
    }

    /// **键按深度切得够细** —— 每一片的深度跨度不超过 [`DEPTH_SLICE`]。
    ///
    /// 这是画家算法的前提:一个图元的深度跨度越大,拿单个深度值给它排序就越
    /// 可能排错。判据逐键算出该有多少片,再与实际发出的片数比 —— 只看"总数
    /// 大于零"的话,把切片整个关掉也照样绿。
    ///
    /// 用棍状样式:它一根键一根圆柱(不并排)、又不画球(不用截),片数因此
    /// 完全由深度跨度决定,算得出闭式。
    #[test]
    fn 键按深度切得够细() {
        for smi in CORPUS {
            let (m, c) = prep(smi);
            let d = depict(&m, &c, &Style3D::STICK).unwrap();
            let want: usize = m
                .bonds()
                .iter()
                .map(|b| {
                    let dz = (d.placed[b.begin as usize].depth - d.placed[b.end as usize].depth)
                        .abs()
                        / 2.0;
                    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
                    let per_half = ((dz / DEPTH_SLICE).ceil() as usize).max(1);
                    2 * per_half
                })
                .sum();
            assert_eq!(
                d.scene.items.len(),
                want,
                "{smi}:{} 根键该切成 {want} 片",
                m.num_bonds()
            );
            // 分母闸:一根键都没有的分子上,上面那条恒真。
            assert!(want >= 2 * m.num_bonds(), "期望片数算错了");
        }
    }
}
