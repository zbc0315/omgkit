//! **端到端:分子进,一组三维坐标出。**
//!
//! ```text
//! 界矩阵 → 三角光滑化 → 取上限矩阵 U 当参考距离表 → 度量矩阵嵌入
//!        → 全局手性定向(离散,一次)→ L-BFGS 精修
//! ```
//!
//! **全程没有随机数。** 同一个分子永远同一个答案,不需要 seed。
//!
//! # 每一步为什么在这儿
//!
//! | 步 | 为什么 |
//! |---|---|
//! | 光滑化 | 之后的上限矩阵 `U` 按构造满足三角不等式,是一张画得出来的距离表 |
//! | 用 `U` 而不是随机取 | RDKit 逐对独立随机取,取出来的表常常摆不出来,它的应对是作废重掷 |
//! | **全局手性定向** | 反射不在 `SO(3)` 连通分支里,**下降法翻不过去**,只能离散地定一次 |
//! | L-BFGS | 目标只有 `C¹`,线搜索必须上强 Wolfe(实测 Armijo 会退化成最速下降) |
//!
//! # 失败是什么意思
//!
//! 只有**界矩阵自相矛盾**(光滑化判不可行)才叫失败 —— 那时连一张自洽的距离表
//! 都没有。实测全语料 8831 个分子里 **1 个(0.01%)**,而 RDKit ETKDGv3 2025.09.2
//! 在同一份语料上失败 36 个(0.41%)—— 见 `harness/baseline_rdkit_etkdg.py`。
//!
//! 嵌入摆不进三维**不算失败**(压掉那一维,精修去救),精修没收敛也**不算失败**
//! (给出当前最好的坐标,并如实报出残差)。理由是:这一步的产物是给力场优化用的
//! **起点**,起点差一点可以修,没有起点才是灾难。

use crate::bounds;
use crate::chiral::{self, Center};
use crate::embed::{self, reference_distances};
use crate::field::Field;
use crate::optimize::{minimize, Options};
use crate::smooth::{triangle_smooth, Bounds, SmoothError};
use omgkit_core::MolBuilder;

/// 生成构型时失败的原因。
// **不再 `Copy`**:`Sanitize` 那一档要把净化的具体原因带上(哪个原子超价),
// 而那是有分配的。翻成一句笼统的"生成失败"等于把排查线索丢掉。
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ConformerError {
    /// 界矩阵自相矛盾 —— 连一张自洽的距离表都拿不出来。
    Infeasible {
        /// 卡住的那一对原子。
        pair: (usize, usize),
    },
    /// 嵌入那一步的输入坏了(非有限数、特征分解不收敛)。这不该发生。
    Embed(crate::embed::EmbedError),
    /// 净化过不去 —— 只有 [`conformer_for`] 会给出它。
    Sanitize(omgkit_chem::SanitizeError),
}

impl core::fmt::Display for ConformerError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Infeasible { pair } => write!(
                f,
                "界矩阵自相矛盾:原子 {} 与 {} 之间的上下界交空",
                pair.0, pair.1
            ),
            Self::Embed(e) => write!(f, "嵌入失败:{e:?}"),
            Self::Sanitize(e) => write!(f, "净化失败:{e}"),
        }
    }
}

/// 一个构型,外加"它有多好"的账。
#[derive(Debug, Clone)]
pub struct Conformer {
    /// 每个原子的坐标。
    pub coords: Vec<[f64; 3]>,
    /// 精修**之前**的误差函数值。
    pub energy_before: f64,
    /// 精修**之后**的误差函数值。
    pub energy: f64,
    /// 精修迭代了多少次。
    pub iterations: usize,
    /// 精修**有没有收敛**(梯度的无穷范数降到 `grad_tol` 以下)。
    ///
    /// # 这一项先前不报,于是没有任何判据看得见"没收敛"
    ///
    /// 只报 `energy` 是不够的:端到端判据报的是**全语料平均**误差
    /// (5.9e1 → 6.4e-4),一个卡在 4e-2 的分子淹在平均里。实测甲醇 `CO`
    /// 就是这样 —— C–O 键出来 1.211 Å 而界是 [1.374, 1.394],精修 10 步之后
    /// 线搜索步长缩到底、当场返回,没收敛。
    pub converged: bool,
    /// 精修停下来时梯度的无穷范数。没收敛时它说明差多远。
    pub grad_norm: f64,
    /// 全局手性定向那一步有没有把结构翻过来。
    pub reflected: bool,
    /// 破对称动了几个原子(嵌入给出重合坐标的那些)。
    ///
    /// **这个数不该常年是 0** —— 对称分子本来就会撞上简并。它一直是 0
    /// 反倒说明这一步没接上,或者判据的样本里没有对称分子。
    pub spread: usize,
    /// 手性中心数,以及精修之后号正确的个数。
    pub chiral_total: usize,
    /// 见 [`Conformer::chiral_total`]。
    pub chiral_ok: usize,
}

/// 残差高于它就走[重试阶梯](conformer)。
///
/// 取"没到 ~0"这条线,而不是某个"够好了"的阈值:规则简单,而且落进局部极小的
/// 分子残差普遍在 1e-3 以上。实测语料里 `large` 10.85%、`hard` 25.00% 的分子
/// 越过这条线,所以阶梯的代价是按这个比例摊的。
const RETRY_RESIDUAL: f64 = 1e-6;

/// 重试阶梯最多试几级。
///
/// 实测 12 个扰动起点里 11 个能收到 ~0,所以 4 级足够;再多是给代价买不到东西。
const RETRY_STEPS: usize = 4;

/// 每一级的扰动幅度(Å),第 `k` 级是 `k` 倍。
///
/// 键长量级 1.4 Å,0.1 Å 的位移足以跳出驻点又不至于把结构打散。
const RETRY_AMPLITUDE: f64 = 0.1;

/// 第 `step` 级给第 `idx` 个坐标分量的扰动,取值在 `[-1, 1)`。
///
/// **整数哈希,不用 `sin`。** libm 的最后几位在不同平台上未必一样,而这个项目
/// 要的是"同一个分子在哪儿跑都给同一组坐标"。splitmix64 的两轮混合足够把
/// 相邻下标打散。
fn jitter(step: usize, idx: usize) -> f64 {
    let mut z = (step as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15)
        ^ (idx as u64).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^= z >> 31;
    #[allow(clippy::cast_precision_loss)]
    let unit = (z >> 11) as f64 / (1u64 << 52) as f64;
    unit.mul_add(2.0, -1.0)
}

/// 一组扁平坐标上有几个手性中心的号是对的。
fn chiral_ok_of(x: &[f64], centers: &[Center]) -> usize {
    let coords: Vec<[f64; 3]> = x.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    chiral::correct_count(&coords, centers)
}

/// 挑选阶梯里"哪一次更好"时,一根键算不算**断了**的判据(Å)。
///
/// 与 `examples/feasibility.rs` 那条判据同一个口径 —— 那不是迎合判据:
/// 键是最硬的一档,参数表给的区间宽只有 0.02 Å,越界 0.1 Å 就是化学上错了。
const BROKEN_BOND_TOL: f64 = 0.1;

/// 一组扁平坐标上有几根键越界超过 [`BROKEN_BOND_TOL`]。
fn broken_bonds_of(x: &[f64], mol: &MolBuilder, b: &Bounds) -> usize {
    mol.bonds()
        .iter()
        .filter(|bd| {
            let (i, j) = (bd.begin as usize, bd.end as usize);
            let d = (0..3)
                .map(|k| (x[3 * i + k] - x[3 * j + k]).powi(2))
                .sum::<f64>()
                .sqrt();
            (b.lower(i, j) - d).max(d - b.upper(i, j)) > BROKEN_BOND_TOL
        })
        .count()
}

/// 精修的迭代上限。
///
/// 400 与 RDKit 第一段极小化的 `field->minimize(400, ...)` 同量级。
/// 到了上限不算失败 —— 给出当前坐标并把残差报出来。
pub const MAX_REFINE_ITER: usize = 400;

/// 从一个**刚解析出来**的分子直接拿到构型。
///
/// # 这五步的顺序不能换
///
/// 1. **净化** —— 后面全部依赖价键、环、芳香性;
/// 2. **感知双键顺反** —— 不做的话 `bounds::stereo_path_torsion` 一次都不发力,
///    顺反整档退回自由旋转(语料里 954 根方向键、404 根双键立体);
/// 3. **补显式氢** —— 氢参与界矩阵,少了它摆出来的是另一个分子。插入顺序按
///    规范秩定,保证同一个分子在哪儿跑都给同一组坐标;
/// 4. **抽手性中心** —— 必须在补氢**之后**,三配位中心的第四格才落得下;
/// 5. 生成([`conformer`])。
///
/// # 为什么它得在库里,而不是每个调用方各抄一遍
///
/// 这个配方先前在 `examples/feasibility.rs`、`examples/dump_conformers.rs`、
/// `examples/bench_conformers.rs`、`tests/small_molecule_geometry.rs` 里各写了
/// 一份,而 Python 绑定还要再写一份 —— 绑定那一层的规矩是"只做翻译,不做化学",
/// 一份抄在那里的配方是整个项目里唯一没有 Rust 判据覆盖的一块。
///
/// **`mol` 会被就地改掉**(净化、补氢),返回的坐标对应改完之后的原子表。
///
/// # Errors
///
/// 净化失败、界矩阵自相矛盾,或嵌入的输入坏掉。
pub fn conformer_for(mol: &mut MolBuilder) -> Result<Conformer, ConformerError> {
    let centers = prepare(mol)?;
    conformer(mol, &centers)
}

/// [`conformer_for`] 的前四步:净化 → 感知双键顺反 → 补显式氢 → 抽手性中心。
///
/// 单独拿出来是给**要在中间插一脚**的调用方用的(比如"只导带立体标记的分子"
/// 那种筛选)。顺序与理由见 [`conformer_for`]。
///
/// # Errors
///
/// 净化失败。
pub fn prepare(mol: &mut MolBuilder) -> Result<Vec<Center>, ConformerError> {
    omgkit_chem::pipeline::sanitize(mol).map_err(ConformerError::Sanitize)?;
    omgkit_io::stereo::perceive_bond_stereo(mol);
    let ranks = omgkit_io::canon::classed_ranks(mol);
    omgkit_chem::add_explicit_hs(mol, &ranks);
    Ok(chiral::centers(mol))
}

/// 给一个分子生成**一个**三维构型。
///
/// # 调用方要先做的两件事
///
/// 1. **补氢**(界矩阵与手性都按显式氢算)。
/// 2. **把 SMILES 的 `/` `\` 折算成双键自己的 `BondStereo`** ——
///    `omgkit_io::stereo::perceive_bond_stereo`。顺反记在**相邻单键**的
///    `direction` 上,不折算的话 `bounds::stereo_path_torsion` 一次都不发力,
///    双键的 1-4 扭转整档退回"顺式到反式的全程",交付的几何会有一半站错边。
///
/// 第 2 条整条流水线先前**压根没做**:实测全语料 405 条双键受影响,
/// 外部判据(RDKit 从坐标读回)上 10 个分子交付的是错的几何。
/// 这不是能靠文档解决的 —— 所以它同时被
/// `omgkit_io::stereo::directions_not_perceived` 这个谓词看着
/// (那个谓词与感知**由构造保证一致**),`examples/feasibility.rs` 拿它当闸。
///
/// `centers` 由调用方给 —— 见 [`crate::chiral::centers`]。
///
/// # Errors
///
/// 界矩阵自相矛盾,或嵌入的输入坏掉。
pub fn conformer(mol: &MolBuilder, centers: &[Center]) -> Result<Conformer, ConformerError> {
    let n = mol.num_atoms();
    let (mut b, _) = bounds::build(mol);
    if let Err(SmoothError::Infeasible { pair }) = triangle_smooth(&mut b) {
        return Err(ConformerError::Infeasible { pair });
    }
    let e = embed::embed(&reference_distances(&b), n).map_err(ConformerError::Embed)?;
    let mut coords = e.coords;

    // **破对称必须在优化器之前。** 对称分子的 Gram 矩阵有重特征值,等价原子会拿到
    // 逐位相同的坐标 —— 而完全重合的两个原子**梯度恰好为零**(方向向量是零向量),
    // 优化器永远分不开它们。实测语料里 0.50% 的分子这样,全语料 44 个,
    // 而且是静默的:坐标照样返回,只是废的。见 `crate::spread`。
    let spread = crate::spread::break_coincidence(&mut coords);

    // **全局手性定向:离散,一次,必须在精修之前。**
    // 反射不在 SO(3) 的连通分支里 —— 连续下降要走到镜像必须把整个分子压平,
    // 下降法不会付这个势垒,所以精修救不了整体定向。
    let reflected = chiral::needs_reflection(&coords, centers);
    if reflected {
        chiral::reflect(&mut coords);
    }

    let field = Field::new(&b, centers);
    let mut x: Vec<f64> = coords.iter().flat_map(|p| p.iter().copied()).collect();
    let mut g = vec![0.0; x.len()];
    let energy_before = {
        use crate::optimize::Objective;
        field.value_and_grad(&x, &mut g)
    };
    let opts = Options {
        max_iter: MAX_REFINE_ITER,
        grad_tol: 1e-6,
        memory: 8,
    };
    let report = minimize(&field, &mut x, &opts);

    // ---- 确定性重试阶梯 ----
    //
    // 精修会落进**局部极小**,而那不是"差一点"—— 甲醇 `CO` 就停在残差
    // 4.04e-2 上,C–O 键出来 1.211 Å 而界是 [1.374, 1.394],短了 0.163。
    // 从产物原地再起一次走 0 步(梯度 4.4e-7,确实是驻点),换参考距离表
    // (下界↔上界之间插值)也一样 —— 键的界宽只有 0.02 Å,插值改不动形状。
    // 而把**起点扰动**一下,12 次里 11 次收到 ~0。甲硫醇、高氯酸、叠氮甲烷同理。
    //
    // 所以这里按阶梯重试:残差没到 ~0 就换个扰动过的起点再跑,取最好的一次。
    // 扰动用整数哈希而不是 `sin` —— libm 的最后几位在不同平台上未必一样,
    // 而这个项目要的是"同一个分子在哪儿跑都给同一组坐标",全程无随机数。
    //
    // 取"最好"时按三把尺子,**顺序不能换**:
    //
    // 1. 手性对了几个 —— 硬不变量(判据要求交付坐标上 100% 正确),
    //    残差小一点换来一个对映体是笔坏买卖;
    // 2. **断了几根键** —— 键是最硬的一档,区间宽只有 0.02 Å;
    // 3. 残差。
    //
    // 第 2 把尺子是量出来才加的:只按(手性, 残差)挑,`hard.smi` 上那个卟啉
    // 分子的越界键从 1 根变成 4 根 —— 阶梯挑了个"总残差更小"的解,而它把误差
    // 摊到了更多键上。总残差小不等于化学上更对。
    let mut best_report = report;
    let mut best_score = (
        chiral_ok_of(&x, centers),
        -isize::try_from(broken_bonds_of(&x, mol, &b)).unwrap_or(isize::MIN),
        -best_report.value,
    );
    let base = x.clone();
    let mut step = 1;
    while best_report.value > RETRY_RESIDUAL && step <= RETRY_STEPS {
        let mut xk = base.clone();
        let amp = RETRY_AMPLITUDE * f64::from(u32::try_from(step).unwrap_or(u32::MAX));
        for (idx, v) in xk.iter_mut().enumerate() {
            *v += amp * jitter(step, idx);
        }
        let r = minimize(&field, &mut xk, &opts);
        let score = (
            chiral_ok_of(&xk, centers),
            -isize::try_from(broken_bonds_of(&xk, mol, &b)).unwrap_or(isize::MIN),
            -r.value,
        );
        if score > best_score {
            best_score = score;
            best_report = r;
            x.copy_from_slice(&xk);
        }
        step += 1;
    }
    let report = best_report;

    for (i, p) in coords.iter_mut().enumerate() {
        *p = [x[3 * i], x[3 * i + 1], x[3 * i + 2]];
    }

    Ok(Conformer {
        chiral_total: centers.len(),
        chiral_ok: chiral::correct_count(&coords, centers),
        coords,
        energy_before,
        spread,
        energy: report.value,
        iterations: report.iterations,
        converged: report.converged,
        grad_norm: report.grad_norm,
        reflected,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn prep(smi: &str) -> MolBuilder {
        let mut m = omgkit_io::smiles::parse(smi).expect("SMILES 该解析得了");
        omgkit_chem::pipeline::sanitize(&mut m).expect("该 sanitize 得了");
        // **这一步先前没有,整条流水线都没有。** 见 `conformer` 的前置条件那一节。
        omgkit_io::stereo::perceive_bond_stereo(&mut m);
        let r = omgkit_io::canon::classed_ranks(&m);
        omgkit_chem::add_explicit_hs(&mut m, &r);
        m
    }

    /// 一条 1-2-3-4 路径在给定坐标下的扭转角(度)。
    fn torsion(p: &[[f64; 3]], a: usize, b: usize, c: usize, d: usize) -> f64 {
        let sub = |u: [f64; 3], v: [f64; 3]| [u[0] - v[0], u[1] - v[1], u[2] - v[2]];
        let dot = |u: [f64; 3], v: [f64; 3]| u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
        let cross = |u: [f64; 3], v: [f64; 3]| {
            [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ]
        };
        let (b0, b1, b2) = (sub(p[a], p[b]), sub(p[c], p[b]), sub(p[d], p[c]));
        let n = dot(b1, b1).sqrt();
        let u = [b1[0] / n, b1[1] / n, b1[2] / n];
        let proj = |v: [f64; 3]| {
            let t = dot(v, u);
            [v[0] - t * u[0], v[1] - t * u[1], v[2] - t * u[2]]
        };
        let (v, w) = (proj(b0), proj(b2));
        // `y.atan2(x)`,**别把两个参数写反** —— 反了得到的是 90° − τ,
        // 而那种错在"离 0 还是离 180"上极难看出来(诊断时踩过)。
        dot(cross(u, v), w).atan2(dot(v, w)).to_degrees()
    }

    #[test]
    fn 双键顺反必须落到正确的一侧() {
        // **期望值写死在这里,不从 `bd.stereo` 读。**
        // 从 `bd.stereo` 读的话,"整批顺反反了"这种变异会让期望和几何一起翻,
        // 永远自洽 —— 实测把 `perceive_bond_stereo` 的 Cis/Trans 对调,
        // 那种写法照样全绿。成对给出顺式与反式也救不了,原因同上。
        //
        // 第二列是**几何事实**:`/F ... /F` 两个 F 在双键两侧(反式,|τ| > 90°)。
        for (smi, cis) in [
            ("F/C=C/F", false),
            ("F/C=C\\F", true),
            ("C/C=C/C", false),
            ("C/C=C\\C", true),
            ("Cl/C(C)=C(C)/Br", false),
            ("Cl/C(C)=C(C)\\Br", true),
            // 环上的那一档 —— 先前交付的几何就是在这里站错边的
            ("[H]/N=C/1\\N[C@]2(CSC(=[NH+]2)N)CS1", true),
            ("[H]/N=C/1\\N=C([C@H](S1)CC(=O)[O-])O", true),
            ("CCOC(=O)[C@@H]1C(=N/C(=N/CC=C)/S1)C", false),
        ] {
            let m = prep(smi);
            // 先确认标记真的被折算出来了 —— 否则这条测试测了个寂寞
            assert!(
                !omgkit_io::stereo::directions_not_perceived(&m),
                "{smi}:有方向键没折算,`prep` 漏了 perceive_bond_stereo"
            );
            let marked: Vec<_> = m
                .bonds()
                .iter()
                .filter(|b| b.stereo != omgkit_core::BondStereo::None)
                .copied()
                .collect();
            assert!(!marked.is_empty(), "{smi}:一根带立体标记的双键都没有");
            let centers = chiral::centers(&m);
            let c = conformer(&m, &centers).unwrap_or_else(|e| panic!("{smi} 失败:{e:?}"));
            // 每个分子的第一根带标记的双键 —— 就是 SMILES 里写明的那一根,
            // 拿**写死的**期望去比
            let bd = marked[0];
            let (i, j) = (bd.stereo_atoms[0] as usize, bd.stereo_atoms[1] as usize);
            let t = torsion(&c.coords, i, bd.begin as usize, bd.end as usize, j);
            assert_eq!(
                t.abs() < 90.0,
                cis,
                "{smi} 键 {}={}({:?}):参照 {i}/{j} 的扭转 {t:.1}°,应当在 {} 一侧",
                bd.begin,
                bd.end,
                bd.stereo,
                if cis {
                    "顺式(|τ|<90°)"
                } else {
                    "反式(|τ|>90°)"
                }
            );
        }
    }

    #[test]
    fn 三配位立体中心抽得出来且两个对映体互为镜像() {
        // 亚砜/亚磺酰胺的 S、膦的 P:三根键 + 一对孤对,构型照样确定。
        // 先前 `centers()` 凑不够四个邻居就整个丢掉 —— 语料 13 个分子、16 个中心
        // 的构型因此是掷硬币。
        // 每一对都必须是**真的对映体**。头一版里 `C[S@](=O)C`(DMSO,两个甲基一样)
        // 与 `c1ccc(cc1)[P@@]2CCCCC2`(环两臂等价)**根本不是立体中心** ——
        // RDKit 直接把标记清掉。拿它们当夹具的话,"号必须相反"是由
        // `sign = match chiral_tag` 直接保证的恒真式,测了个寂寞。
        // 下面每一条都用 RDKit 规范化确认过标记不会被清掉。
        for (a, b) in [
            ("C[S@](=O)CC", "C[S@@](=O)CC"),
            ("C[S@](=O)c1ccccc1", "C[S@@](=O)c1ccccc1"),
            (
                "C[C@@H]1CO[S@@](=O)N1c2ccccc2",
                "C[C@@H]1CO[S@](=O)N1c2ccccc2",
            ),
            ("C[P@H]CC", "C[P@@H]CC"),
            ("CC[P@](C)c1ccccc1", "CC[P@@](C)c1ccccc1"),
            ("C[C@@H]1CC[P@](c2ccccc2)C1", "C[C@@H]1CC[P@@](c2ccccc2)C1"),
        ] {
            let (ma, mb) = (prep(a), prep(b));
            let (ca, cb) = (chiral::centers(&ma), chiral::centers(&mb));
            // 抽得出来,而且确实是三配位那一档
            let three_a = ca.iter().filter(|c| c.is_three_coordinate()).count();
            assert!(three_a > 0, "{a}:一个三配位中心都没抽出来");
            assert_eq!(
                three_a,
                cb.iter().filter(|c| c.is_three_coordinate()).count()
            );
            // **只比三配位那些中心** —— 上面几对里有的只翻了 S / P,
            // 分子里别的手性中心两边一样,号本来就该相同。
            // (这条断言的第一版没分,当场被 `C[C@@H]1CO[S@@]…` 那一对打红,
            //  它的碳在两边完全一致。)
            //
            // 断的是"两个对映体给出相反的号"。它抓不住"整批号取反"
            // (那种错法两边一起翻,仍然相反),但抓得住"根本没读标记"
            // —— 那时两边会给同一个号。绝对约定由外部判据 `verify_stereo.py` 守,
            // 而三配位 **P** 那一档判官够不着(RDKit 自己都读不回),
            // 所以 P 的绝对号目前只有这条必要条件钉着,别当成已验证。
            for (x, y) in ca.iter().zip(cb.iter()) {
                assert_eq!(x.atom, y.atom);
                if !x.is_three_coordinate() {
                    continue;
                }
                assert_eq!(
                    x.sign, -y.sign,
                    "{a} / {b}:对映体在 {} 号上应当相反,实得 {} vs {}",
                    x.atom, x.sign, y.sign
                );
            }
            // 交付的坐标要真的照着摆:两边的号都得对上自己的目标
            let (fa, fb) = (
                conformer(&ma, &ca).unwrap_or_else(|e| panic!("{a}:{e:?}")),
                conformer(&mb, &cb).unwrap_or_else(|e| panic!("{b}:{e:?}")),
            );
            assert_eq!(fa.chiral_ok, fa.chiral_total, "{a}:交付的号不对");
            assert_eq!(fb.chiral_ok, fb.chiral_total, "{b}:交付的号不对");
            // 而且两边的**中心基点体积**必须反号 —— 这才说明几何真的镜像了,
            // 不是两边都朝同一个方向摆然后判据也跟着错。
            for (x, y) in ca.iter().zip(cb.iter()) {
                if !x.is_three_coordinate() {
                    continue;
                }
                let (va, vb) = (
                    chiral::center_volume(&fa.coords, x),
                    chiral::center_volume(&fb.coords, y),
                );
                assert!(
                    va * vb < 0.0,
                    "{a} / {b}:中心 {} 的体积没反号({va:+.3} vs {vb:+.3})",
                    x.atom
                );
            }
        }
    }

    #[test]
    fn 三配位但没有孤对的不算立体中心() {
        // `[S+]` 三配位是平面的,`[N+]`/`[C]` 三配位是 sp² —— 都不是四面体中心。
        // 这条守的是"别把任何凑不够四邻居的带标记原子都当成三配位中心"。
        for smi in ["C[S+](C)C", "C[N+](C)C"] {
            let mut m = omgkit_io::smiles::parse(smi).expect("解析");
            omgkit_chem::pipeline::sanitize(&mut m).expect("净化");
            // 强行按四面体标记,看 `centers` 会不会上钩
            for i in 0..m.num_atoms() as u32 {
                let z = m.atoms()[i as usize].atomic_num;
                if z == 16 || z == 7 {
                    if let Some(a) = m.atom_mut(i) {
                        a.chiral_tag = omgkit_core::ChiralTag::Ccw;
                    }
                }
            }
            let r = omgkit_io::canon::classed_ranks(&m);
            omgkit_chem::add_explicit_hs(&mut m, &r);
            for c in chiral::centers(&m) {
                assert!(
                    !c.is_three_coordinate(),
                    "{smi}:原子 {} 不该被当成三配位立体中心",
                    c.atom
                );
            }
        }
    }

    #[test]
    fn 漏了顺反折算会被谓词看见() {
        // `directions_not_perceived` 是把"前置条件"变成机器可查的那一条。
        // 它一旦恒为 false,`feasibility` 那道闸就什么都没守住。
        let mut m = omgkit_io::smiles::parse("F/C=C/F").expect("解析");
        omgkit_chem::pipeline::sanitize(&mut m).expect("净化");
        // **故意不调** perceive_bond_stereo
        let r = omgkit_io::canon::classed_ranks(&m);
        omgkit_chem::add_explicit_hs(&mut m, &r);
        assert!(
            omgkit_io::stereo::directions_not_perceived(&m),
            "漏了折算却没被谓词看见 —— 那道前置条件闸是瞎的"
        );
    }

    #[test]
    fn 常见分子都给得出构型() {
        for smi in [
            "CCO",
            "c1ccccc1",
            "C1CCCCC1",
            "CC(=O)Nc1ccc(O)cc1",
            "C1CC2CCC1CC2",
            "CC(C)(C)OC(=O)N1CCC(CC1)N",
            "FS(F)(F)(F)(F)F",
            "C=C=C",
        ] {
            let m = prep(smi);
            let c = conformer(&m, &[]).unwrap_or_else(|e| panic!("{smi} 失败:{e:?}"));
            assert_eq!(c.coords.len(), m.num_atoms(), "{smi} 坐标数不对");
            for (i, p) in c.coords.iter().enumerate() {
                assert!(
                    p.iter().all(|v| v.is_finite()),
                    "{smi} 第 {i} 个原子坐标不是有限数:{p:?}"
                );
            }
            // **精修必须真的降下去**,不能原地不动
            assert!(
                c.energy <= c.energy_before,
                "{smi} 精修之后反而更差:{} → {}",
                c.energy_before,
                c.energy
            );
        }
    }

    #[test]
    fn 精修确实在干活() {
        // 起点是嵌入给的,残差不小;精修应当把它压掉一大截。
        let m = prep("CC(=O)Nc1ccc(O)cc1");
        let c = conformer(&m, &[]).unwrap();
        assert!(
            c.energy_before > 0.1,
            "起点残差 {} 太小,测不到东西",
            c.energy_before
        );
        assert!(
            c.energy < c.energy_before * 0.2,
            "只压掉了 {:.1}%:{} → {}",
            100.0 * (1.0 - c.energy / c.energy_before),
            c.energy_before,
            c.energy
        );
    }

    #[test]
    fn 同一个分子两次给逐位相同的坐标() {
        // 全程无随机数 —— 这条一旦红,说明哪里混进了不确定性
        // (HashMap 迭代序、并行归约、未定义的排序 tie-break)。
        let m = prep("CC(C)(C)OC(=O)N1CCC(CC1)N");
        let a = conformer(&m, &[]).unwrap();
        let b = conformer(&m, &[]).unwrap();
        assert_eq!(a.coords, b.coords, "两次跑出来的坐标不逐位相同");
        assert_eq!(a.energy, b.energy);
        assert_eq!(a.iterations, b.iterations);
    }

    #[test]
    fn 界不可行时如实报失败() {
        // 手工造一个自相矛盾的界。这里直接调底层,确认错误一路传得上来。
        use crate::smooth::Bounds;
        let mut b = Bounds::new(3, 0.0, 10.0);
        b.set_lower(0, 1, 5.0);
        b.set_upper(0, 1, 5.0);
        b.set_lower(1, 2, 5.0);
        b.set_upper(1, 2, 5.0);
        b.set_lower(0, 2, 50.0);
        b.set_upper(0, 2, 50.0);
        assert!(triangle_smooth(&mut b).is_err(), "这组界本来就该判不可行");
    }
}
