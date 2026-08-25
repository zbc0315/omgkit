//! **距离几何的误差函数** —— 精修阶段要最小化的那个目标。
//!
//! 原子对的距离违反,以及手性中心的有符号体积违反(**体积那边是两条**)。
//!
//! # 距离项:沿用 RDKit 的形状,换掉参与的对
//!
//! 逐对:落在 `[l, u]` 里罚 0,出界才罚。
//!
//! ```text
//! d > u :  val = d²/u² − 1
//! d < l :  val = 2l²/(l² + d²) − 1
//! E    += w · val²
//! ```
//!
//! **形状照抄是有理由的**,不是省事:
//!
//! - 罚的是**相对**越界,不是绝对。同样偏 0.1 Å,1.5 Å 的键上比 10 Å 的长程对
//!   贵约 7 倍 —— 而嵌入恰恰是拿键长换长程的(经典 MDS 拟合平方距离,
//!   大距离天然占主导),所以这个加权正是需要的。
//! - 下限那一支在 `d → 0` 时**饱和**到 1,梯度也有界(`preFactor/d → −8/l²`)。
//!   看着像缺陷,其实是故意的:两个原子完全重叠时罚项不会把优化器拽爆。
//!   代价是"穿过去"的势垒有限(最多 `w`),所以立体化学与自穿要靠**别的**判据守,
//!   不能指望这一项(见 [`crate::chiral`] 与 [`crate::threading`])。
//!
//! **参与的对是全部 `N²`,不是 RDKit 默认的 `u − l ≤ 5.0`。** 实测那条过滤在语料上
//! 丢掉 16.8% 的对(柔性大分子最高 51%),而那正是自穿会发生的地方。
//! 这也不是我们的发明:`Embedder.cpp:866` 里 RDKit 走随机坐标时自己就把
//! `basinThresh` 设成 `1e8`(等于全部对),注释理由一字不差。
//!
//! # 手性项:有符号体积落进区间。**是两条,不是一条**
//!
//! `E += w·(V − lo)²`(`V < lo` 时)或 `w·(V − hi)²`(`V > hi` 时),而 `V` 有两个:
//!
//! | 体积 | 式子 | 它管什么 | 它看不见什么 |
//! |---|---|---|---|
//! | 四配体 | `det[l₁−l₀, l₂−l₀, l₃−l₀]` | 整体定向、别把四个配体压平 | **中心原子在哪** |
//! | 中心基点 | `det[l₀−c, l₁−c, l₂−c]` | 别伞形翻转 | 第四个配体 |
//!
//! RDKit 只有前者(`ChiralViolationContrib`),于是中心被挤到配体四面体外面时
//! 它没有任何回复力。
//!
//! # 第二条现在**一次都没在罚**,留着是因为代价是零、堵的洞是真的
//!
//! 实测三栏(`smoke.chirality.jsonl`,247 个中心):
//!
//! | | 精修步数 | 1-2 键 | 1-3 角 | 真值口径手性 |
//! |---|---|---|---|---|
//! | 手性项**全关** | 167 | 0.0% | 0.4% | **64.8%** |
//! | 只有四配体项(`WEIGHT_UMBRELLA = 0`) | 169 | 0.0% | 0.1% | **100%** |
//! | **两条都有** | 168 | 0.0% | 0.1% | 100% |
//!
//! **第二行与第三行分不出来。** 也就是说中心基点这一项在当前语料上
//! 一个中心都没救到 —— 交付坐标里中心原子在配体四面体外的是 **0 个**。
//! 距离项已经间接挡住了大半(中心跑出去就得拉长某条键或压缩某个 1-3)。
//!
//! 头一版这张表的第二行标的是"只有四配体项 … 64.8%",**那是标错的**:
//! 当时代码里只剩中心基点一项,把权重清零关掉的是**全部**手性项。
//! 独立审核复现时拆穿了这一点,这里照实改。
//!
//! 留着它的理由不是"它修好了什么",是:四配体行列式对伞形翻转**在数学上就是瞎的**,
//! 而这一档真出事就是灾难 —— 给 300 个分子各翻一个中心的伞,外部判据
//! `verify_stereo.py` 从 290/301 掉到 **2/301**。代价实测为零(168 vs 169 步),
//! 且有一条判据钉着它(`手性项挡得住伞形翻转`,把 `WEIGHT_UMBRELLA` 清零当场红)。
//!
//! # 只用中心基点项**不行**,那是扫过参数才确认的
//!
//! | `UMBRELLA_LO` | 0.2 | 0.3 | 0.8 | 1.5 |
//! |---|---|---|---|---|
//! | 1-3 角越界 | 4.1% | 4.0% | 4.1% | 4.7% |
//!
//! 权重也从 0.5 扫到 10,1-3 角一直在 4.1–4.8% —— 不是调参能解决的,
//! 所以两条并存而不是替换。(另外拿 `VOL_LO = 5` 直接配中心基点公式跑过一次,
//! 键长越界 19.8%:两个体积的尺度差 3.2 倍,常数不能照搬。)
//!
//! # 对齐号约定要翻**区间**,不能翻向量
//!
//! 两条的号相反(中心取四配体质心时 `V_配体 = −4·V_中心`)。对齐的办法是给目标
//! 区间取反号 —— 给向量取负会让梯度与能量对不上号,而这个坑我踩过一次。
//! 注意**数值差分校验只在夹具非退化时才看得见它**:见
//! `两项一起的梯度也一致` 与 `手性项的梯度与能量一致` 里那两条护栏。
//!
//! **梯度是自己推的,不是从 RDKit 抄的。** 那边 `ChiralViolationContribs` 的
//! 能量是 `w(V−lo)²` 而 `preFactor` 是 `w(V−lo)` —— **差一个因子 2**
//! (第四维项同样如此)。RDKit 的极小化是梯度驱动的,所以它那两项的**实际权重
//! 是名义值的一半**;而 L-BFGS 对"能量与梯度不是同一个函数"零容忍。
//! 所以这里写数学上一致的梯度,并用 [`crate::optimize::max_grad_error`] 钉住。
//!
//! # 还没有的
//!
//! 第四维项**不做** —— 理由是量出来的:一次离散的全局反射就把手性正确率
//! 从 53% 拉到 86%,剩下的是个别中心错,那一档三维精修够得着。见 [`crate::chiral`]。

use crate::chiral::Center;
use crate::optimize::Objective;
use crate::smooth::Bounds;

/// **四配体项**的体积目标下限(绝对值)。与 RDKit 的 `volLowerBound` 同值。
///
/// **它同时是一道"别压平"的闸**:只要求符号对的话,把四个配体压成近乎共面、
/// 体积 `+1e−6`,符号照样对而分子是废的。要求 `|V| ≥ 5` 就堵住了这条路。
///
/// 实测 247 个真实构象上的中心:`|V_配体|` 最小 5.882、中位 8.374、最大 16.446,
/// 所以 5.0 在实测最小值下面留了 15% 余量。
pub const VOL_LO: f64 = 5.0;

/// 四配体项的体积目标上限(绝对值)。与 RDKit 的 `volUpperBound` 同值。
pub const VOL_HI: f64 = 100.0;

/// **中心基点项**的体积目标下限(绝对值)—— 专治伞形翻转。
///
/// # 为什么必须单立一项
///
/// 四配体行列式 `det[l₁−l₀, l₂−l₀, l₃−l₀]` **完全不看中心原子在哪**。
/// 中心被挤到配体四面体外面时它一点变化都没有,而真实立体化学已经翻了 ——
/// 于是优化器对这件事**没有任何回复力**。这一档当前语料上没有在发生
/// (交付坐标里 484 个中心,0 个在四面体外),留着是因为代价为零;
/// 详见 [`crate::chiral::center_volume`] 与本模块文档那一段。
///
/// # 阈值为什么取得这么小
///
/// 这一项的活**只是把号守住**,不是去塑形几何。实测扫过 0.2 → 1.5,
/// 真值手性一律 100%,而 1-3 角越界 4.1 / 4.0 / 4.1 / 4.7%(0.2 / 0.3 / 0.8 / 1.5)
/// —— 不单调,但大了只有坏处。
///
/// 另外它不能照 `VOL_LO / 3.234` 折算:那是拿 RDKit 的配体序算出来的比值,
/// 而**换一组三配体,`|V_中心|` 能差三成**(三条键长的乘积不同)。
///
/// 取 **0.3**。余量要按**我们自己的配体序、在我们交付的坐标上**量,不能拿
/// RDKit 那份 MMFF 优化过的基准(那上面最小 1.925,会算出 6 倍余量,是错的):
/// 实测全语料 484 个中心最小 **0.7948**(`C1COC2(O1)[C@@]3(C[C@@]3(C(=[NH+]2)N)C#N)C#N`
/// 的螺环丙烷中心),余量 **2.65 倍**。这一项在收敛点上一个中心都没在罚 ——
/// 它只在号已经错了(`V` 跑到区间对面)的时候才发力。
pub const UMBRELLA_LO: f64 = 0.3;

/// 中心基点项的体积目标上限(绝对值)。实测 `|V_中心|` 最大 4.043,
/// 取 **30** 等于不设上限 —— 这一项不负责压制过大的体积,那是四配体项的活。
pub const UMBRELLA_HI: f64 = 30.0;

/// 中心基点项相对四配体项的权重倍数。
///
/// 由端到端判据定,不是推出来的 —— 见 [`WEIGHT_CHIRAL`] 的说明。
pub const WEIGHT_UMBRELLA: f64 = 1.0;

/// 手性项的默认权重。
///
/// RDKit 名义上用 1.0,但它的手性梯度少了因子 2,**实际权重是 0.5**。
/// 这里的梯度是对的,所以要复现它的行为该取 0.5 —— 不过这是个自由参数,
/// 该取多少由端到端的判据说了算,不是照抄一个数。
pub const WEIGHT_CHIRAL: f64 = 1.0;

/// 距离几何的误差函数。
#[derive(Debug, Clone)]
pub struct Field {
    n: usize,
    /// 逐对下限,行主序 `n×n`,只用上三角。
    lower: Vec<f64>,
    /// 逐对上限,同上。
    upper: Vec<f64>,
    centers: Vec<Center>,
    /// 手性项的权重。
    pub weight_chiral: f64,
}

impl Field {
    /// 从界矩阵与手性中心建一个误差函数。
    ///
    /// **界矩阵应当是光滑化之后的** —— 没光滑过的界里有大量废话上限
    /// (`MAX_UPPER = 1000`),那些对等于没有约束。
    ///
    /// # Panics
    ///
    /// 任何一个 [`Center`] 的中心原子或配体下标 ≥ 原子数时 panic。
    ///
    /// **在这里查,不留到热循环里去。** `Center::atom` 自从手性项改用中心基点
    /// 之后成了热路径索引,越界会变成 `value_and_grad` 深处一句没有上下文的
    /// "index out of bounds";而 `Center` 的字段全是 `pub`,手写调用方够得着。
    #[must_use]
    pub fn new(b: &Bounds, centers: &[Center]) -> Self {
        let n = b.len();
        for c in centers {
            assert!(
                (c.atom as usize) < n,
                "手性中心的中心原子下标 {} 越界(原子数 {n})",
                c.atom
            );
            for &l in c.real_ligands() {
                assert!(
                    (l as usize) < n,
                    "手性中心 {} 的配体下标 {l} 越界(原子数 {n})",
                    c.atom
                );
            }
        }
        let mut lower = vec![0.0; n * n];
        let mut upper = vec![0.0; n * n];
        for i in 0..n {
            for j in (i + 1)..n {
                lower[i * n + j] = b.lower(i, j);
                upper[i * n + j] = b.upper(i, j);
            }
        }
        Self {
            n,
            lower,
            upper,
            centers: centers.to_vec(),
            weight_chiral: WEIGHT_CHIRAL,
        }
    }

    /// 原子数。
    #[must_use]
    pub fn len(&self) -> usize {
        self.n
    }

    /// 有没有原子。
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    /// 只算距离项(诊断与单元测试用)。
    fn distance_term(&self, x: &[f64], g: &mut [f64]) -> f64 {
        let n = self.n;
        let mut e = 0.0;
        // **固定的双重循环,不许并行归约** —— 浮点加法不结合,顺序一变答案就变,
        // 而这个 crate 承诺同一输入同一输出。
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = [
                    x[3 * i] - x[3 * j],
                    x[3 * i + 1] - x[3 * j + 1],
                    x[3 * i + 2] - x[3 * j + 2],
                ];
                let d2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
                let (lo, hi) = (self.lower[i * n + j], self.upper[i * n + j]);
                let (lo2, hi2) = (lo * lo, hi * hi);

                // 分三支:超上限、低于下限、落在区间里(罚 0)
                let dedd_over_d = if d2 > hi2 && hi2 > 0.0 {
                    let val = d2 / hi2 - 1.0;
                    e += val * val;
                    // E = val²,dE/d(d²) = 2·val/hi²,而 dE/dx = dE/d(d²)·2·dx
                    4.0 * val / hi2
                } else if d2 < lo2 {
                    let s = lo2 + d2;
                    let val = 2.0 * lo2 / s - 1.0;
                    e += val * val;
                    // dval/d(d²) = −2lo²/s²,dE/d(d²) = 2·val·(−2lo²/s²)
                    -8.0 * val * lo2 / (s * s)
                } else {
                    continue;
                };
                for t in 0..3 {
                    let f = dedd_over_d * dx[t];
                    g[3 * i + t] += f;
                    g[3 * j + t] -= f;
                }
            }
        }
        e
    }

    /// 只算手性项(诊断与单元测试用)。
    fn chiral_term(&self, x: &[f64], g: &mut [f64]) -> f64 {
        let mut e = 0.0;
        for c in &self.centers {
            let p: Vec<[f64; 3]> = c
                .real_ligands()
                .iter()
                .map(|&a| {
                    let a = a as usize;
                    [x[3 * a], x[3 * a + 1], x[3 * a + 2]]
                })
                .collect();
            let sub = |u: [f64; 3], v: [f64; 3]| [u[0] - v[0], u[1] - v[1], u[2] - v[2]];
            let cross = |u: [f64; 3], v: [f64; 3]| {
                [
                    u[1] * v[2] - u[2] * v[1],
                    u[2] * v[0] - u[0] * v[2],
                    u[0] * v[1] - u[1] * v[0],
                ]
            };
            let ctr = c.atom as usize;
            let o = [x[3 * ctr], x[3 * ctr + 1], x[3 * ctr + 2]];
            // 一个通用的"体积落在 [lo, hi] 之内"罚项。
            //
            // `v = a·(b×cc)`,其中 `a, b, cc` 是三个**差向量**,各自的基点原子是
            // `base`;`idx` 给出这三个差向量的头端原子。
            // `want` 是这个体积**该有的号**。两项的号约定相反,所以由调用方给,
            // **不许**靠给向量取负来对齐 —— 那会让梯度与能量对不上号
            // (`term` 按 `v = a·(b×cc)` 推导数,传进去 `−a` 就差一个号)。
            let term = |a: [f64; 3],
                        b: [f64; 3],
                        cc: [f64; 3],
                        idx: [usize; 3],
                        base: usize,
                        want: f64,
                        lo_abs: f64,
                        hi_abs: f64,
                        w: f64,
                        g: &mut [f64]| {
                let bxc = cross(b, cc);
                let v = a[0] * bxc[0] + a[1] * bxc[1] + a[2] * bxc[2];
                let (lo, hi) = if want < 0.0 {
                    (-hi_abs, -lo_abs)
                } else {
                    (lo_abs, hi_abs)
                };
                let dev = if v < lo {
                    v - lo
                } else if v > hi {
                    v - hi
                } else {
                    return 0.0;
                };
                // dE/dV = 2·w·dev。**因子 2 不能少** —— RDKit 那边就少了它,
                // 于是它的手性项实际权重只有名义值的一半。
                let k = 2.0 * w * dev;
                let cxa = cross(cc, a);
                let axb = cross(a, b);
                for t in 0..3 {
                    g[3 * idx[0] + t] += k * bxc[t];
                    g[3 * idx[1] + t] += k * cxa[t];
                    g[3 * idx[2] + t] += k * axb[t];
                    // 基点原子拿走三者之和的负数 —— 平移不变性要求梯度和为零
                    g[3 * base + t] -= k * (bxc[t] + cxa[t] + axb[t]);
                }
                w * dev * dev
            };

            // ---- 一、四配体项:管"整体定向"与"别压平" ----
            //
            // **三配位中心走不到这一项** —— 第四个"配体"是孤对电子,没有坐标。
            // 那一档只剩中心基点项,而它本来就只用三个配体,够用:
            // 号由它定,"别压平"也由它的 `UMBRELLA_LO` 兜着
            // (中心落到三个配体所在平面上时 `V → 0`,正是要挡的那件事)。
            if !c.is_three_coordinate() {
                // 这一项与 RDKit 的 `ChiralViolationContrib` 同形,几何代价最低
                // (实测端到端 1-3 角越界 0.1%)。但它**完全不看中心原子在哪**。
                let (a, b, cc) = (sub(p[1], p[0]), sub(p[2], p[0]), sub(p[3], p[0]));
                let idx = [
                    c.ligands[1] as usize,
                    c.ligands[2] as usize,
                    c.ligands[3] as usize,
                ];
                // `Center::sign` 是按**中心基点**标定的,而这一项算的是四配体行列式,
                // 两者反号(正四面体上 `V_配体 = −4·V_中心`)—— 所以目标号取 `−sign`。
                e += term(
                    a,
                    b,
                    cc,
                    idx,
                    c.ligands[0] as usize,
                    -c.sign,
                    VOL_LO,
                    VOL_HI,
                    self.weight_chiral,
                    g,
                );
            }

            // ---- 二、中心基点项:只管"别翻伞" ----
            // 中心原子被挤到配体四面体**外面**时,四配体行列式号不变而真实立体
            // 化学已经翻了 —— 第一项对这件事没有任何回复力。见 `chiral::center_volume`。
            //
            // 阈值取得**小**:它的活只是把号守住,不是去塑形几何。
            let (a2, b2, c2) = (sub(p[0], o), sub(p[1], o), sub(p[2], o));
            let idx2 = [
                c.ligands[0] as usize,
                c.ligands[1] as usize,
                c.ligands[2] as usize,
            ];
            e += term(
                a2,
                b2,
                c2,
                idx2,
                ctr,
                c.sign,
                UMBRELLA_LO,
                UMBRELLA_HI,
                self.weight_chiral * WEIGHT_UMBRELLA,
                g,
            );
        }
        e
    }
}

impl Objective for Field {
    fn value_and_grad(&self, x: &[f64], grad: &mut [f64]) -> f64 {
        // **非有限坐标必须当场变成 NaN,不能悄悄罚 0。**
        //
        // 距离项分三支比大小(`d2 > hi2` / `d2 < lo2` / 落在区间里),而 NaN
        // 与任何数相比都是 false —— 带 NaN 的原子对于是走第三支,**罚 0、
        // 梯度 0**。实测四个原子里放一个 NaN:误差 0、梯度全 0,
        // `minimize` 当场报 `converged=true, grad_norm=0, value=0` ——
        // 一个废掉的结构拿到满分,而且是最高分。
        //
        // 这一遍是 O(n),外面那两层循环是 O(n²)。实测全语料
        // (`dump_conformers -- large.smi`,各跑三次):有守卫 0.70–0.90 s,
        // 拿掉守卫 0.71–1.05 s —— 两个区间叠在一起,代价在跑次波动之下。
        if !x.iter().all(|v| v.is_finite()) {
            for v in grad.iter_mut() {
                *v = f64::NAN;
            }
            return f64::NAN;
        }
        for v in grad.iter_mut() {
            *v = 0.0;
        }
        self.distance_term(x, grad) + self.chiral_term(x, grad)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::optimize::{max_grad_error, minimize, Options};

    /// 造一个 `n` 个原子的界:界照一组给定坐标定成 `真实距离 ± half`,
    /// 于是那组坐标本身就是零点。
    fn field_from(points: &[[f64; 3]], half: f64, centers: &[Center]) -> Field {
        let n = points.len();
        let mut b = Bounds::new(n, 0.0, 1000.0);
        for i in 0..n {
            for j in (i + 1)..n {
                let d = ((points[i][0] - points[j][0]).powi(2)
                    + (points[i][1] - points[j][1]).powi(2)
                    + (points[i][2] - points[j][2]).powi(2))
                .sqrt();
                b.set_lower(i, j, (d - half).max(0.05));
                b.set_upper(i, j, d + half);
            }
        }
        Field::new(&b, centers)
    }

    /// **非有限坐标不许被判成满分。**
    ///
    /// 修之前四个入口**互相独立地**把 NaN 洗成最好的分数,实测(四个原子
    /// 里放一个 NaN):
    ///
    /// | 入口 | 修之前 | 应有 |
    /// |---|---|---|
    /// | `Field::value_and_grad` | 误差 **0**、梯度全 **0** | NaN |
    /// | `minimize` | `converged=true, grad_norm=0, value=0` | 不收敛 |
    /// | `max_grad_error` | 偏差 **0** | NaN |
    ///
    /// 根因两条:距离项分三支比大小,而 NaN 与任何数比都是 false,于是带 NaN
    /// 的原子对走"落在区间里"那一支、罚 0;`f64::max` 遇 NaN 返回另一个操作数,
    /// 于是"最坏偏差"的归约把 NaN 洗成 0(见
    /// [`crate::linalg::max_nan_wins`])。
    ///
    /// `+inf` 也要走同一条路 —— 它经一次减法就变成 NaN。
    #[test]
    fn 非有限坐标不许被判成满分() {
        let pts = [
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [0.0, 1.5, 0.0],
            [0.0, 0.0, 1.5],
        ];
        let f = field_from(&pts, 0.1, &[]);
        let clean: Vec<f64> = pts.iter().flatten().copied().collect();

        // **先证明干净坐标上这份场是满分的**,否则下面每一条都可能是空断言:
        // 一个恒返回 NaN 的实现也能让它们全绿。
        let mut g = vec![0.0; clean.len()];
        assert_eq!(
            f.value_and_grad(&clean, &mut g),
            0.0,
            "界就是照这组点定的,它本身该是零点"
        );
        assert!(g.iter().all(|v| *v == 0.0), "零点上梯度该是零");
        assert!(
            max_grad_error(&f, &clean, 1e-5) < 1e-6,
            "干净坐标上梯度该对"
        );
        let r = minimize(&f, &mut clean.clone(), &Options::default());
        assert!(r.converged && r.value.is_finite(), "干净坐标上该收敛");

        // 三种非有限值 × 三个不同位置(第一个原子、中间、最后一个分量)
        for (tag, bad) in [
            ("NaN", f64::NAN),
            ("+inf", f64::INFINITY),
            ("-inf", f64::NEG_INFINITY),
        ] {
            for poison in [0usize, 5, 11] {
                let mut x = clean.clone();
                x[poison] = bad;
                let mut g = vec![0.0; x.len()];
                let e = f.value_and_grad(&x, &mut g);
                assert!(e.is_nan(), "{tag}@{poison}:误差该是 NaN,实得 {e}");
                assert!(
                    g.iter().all(|v| v.is_nan()),
                    "{tag}@{poison}:梯度该整条 NaN,实得 {g:?}"
                );
                assert!(
                    max_grad_error(&f, &x, 1e-5).is_nan(),
                    "{tag}@{poison}:梯度校验不许报出一个有限偏差"
                );
                let mut xx = x.clone();
                let r = minimize(&f, &mut xx, &Options::default());
                assert!(!r.converged, "{tag}@{poison}:不许报收敛");
                assert!(
                    !r.grad_norm.is_finite(),
                    "{tag}@{poison}:梯度范数不许是有限数,实得 {}",
                    r.grad_norm
                );
                assert!(
                    !r.value.is_finite(),
                    "{tag}@{poison}:目标值不许是有限数,实得 {}",
                    r.value
                );
            }
        }
    }

    fn lcg(state: &mut u64) -> f64 {
        *state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1);
        ((*state >> 11) as f64) / ((1u64 << 53) as f64) * 2.0 - 1.0
    }

    fn random_points(n: usize, seed: u64, scale: f64) -> Vec<[f64; 3]> {
        let mut st = seed;
        (0..n)
            .map(|_| {
                [
                    lcg(&mut st) * scale,
                    lcg(&mut st) * scale,
                    lcg(&mut st) * scale,
                ]
            })
            .collect()
    }

    fn flat(p: &[[f64; 3]]) -> Vec<f64> {
        p.iter().flat_map(|q| q.iter().copied()).collect()
    }

    #[test]
    fn 距离项的梯度与能量一致() {
        // **这一条是这个模块的地基。** L-BFGS 对能量/梯度不一致零容忍,
        // 而 RDKit 的手性项与第四维项恰好差一个因子 2 —— 照抄就会栽在这里。
        let pts = random_points(10, 5, 3.0);
        let f = field_from(&pts, 0.3, &[]);
        // 起点要打乱,让上下两支都被激活
        let mut st = 77u64;
        let x: Vec<f64> = flat(&pts).iter().map(|v| v + lcg(&mut st) * 2.0).collect();
        let e = max_grad_error(&f, &x, 1e-6);
        assert!(e < 1e-6, "距离项梯度对不上:{e:.3e}");
    }

    /// 断言**两条手性项都真的在罚**。
    ///
    /// 只问"手性项合起来非零"是不够的:四配体项在罚就够让那个断言绿,
    /// 而伞形项照样可以落在区间内、能量与梯度都是 0 —— 于是伞形项的梯度
    /// 从来没被数值差分碰过。实测过:那样的话给伞形项的向量取负(号错)
    /// 全套测试照样绿。
    fn assert_both_chiral_active(pts: &[[f64; 3]], c: &Center) {
        let inside = |v: f64, lo: f64, hi: f64, want: f64| {
            if want < 0.0 {
                v <= -lo && v >= -hi
            } else {
                v >= lo && v <= hi
            }
        };
        let l = c.ligands.map(|k| pts[k as usize]);
        let vl = crate::chiral::signed_volume(l[0], l[1], l[2], l[3]);
        let vc = crate::chiral::center_volume(pts, c);
        assert!(
            !inside(vl, VOL_LO, VOL_HI, -c.sign),
            "sign={} 四配体项没在罚(V={vl}),这一半测了个寂寞",
            c.sign
        );
        assert!(
            !inside(vc, UMBRELLA_LO, UMBRELLA_HI, c.sign),
            "sign={} 伞形项没在罚(V={vc}),这一半测了个寂寞",
            c.sign
        );
    }

    #[test]
    fn 手性项的梯度与能量一致() {
        // 下标 0 是**中心原子**,1..=4 是配体。先前这里写的是
        // `atom: 0, ligands: [0,1,2,3]` —— 中心就是第一个配体,于是伞形项的
        // `a₂ = l₀ − c ≡ 0`,那一项给出**常数能量、恒零梯度**,数值差分
        // 永远看不见它。这是全套里唯一一条专门差分手性项的测试,
        // 它退化就等于伞形项的梯度完全没人验。
        // 夹具是**算出来的**:要让两项对**两个号都**违反,只需把两个体积都做小 ——
        // `|V| < lo` 时,`want = +1` 那支 `V < lo` 违反,`want = −1` 那支
        // `V > −lo` 也违反。所以:四个配体近共面(`|V_配体| ≈ 0.3 < VOL_LO`),
        // 中心近落在 `l₀l₁l₂` 平面上(`|V_中心| ≈ 0.05 < UMBRELLA_LO`)。
        // 位置故意不对称,免得下标搞混了还能碰巧对上。
        let pts: Vec<[f64; 3]> = vec![
            [0.03, -0.02, 0.02], // 中心
            [1.10, 0.00, 0.00],
            [0.00, 0.95, 0.00],
            [-1.05, 0.00, 0.00],
            [0.00, -1.00, 0.15],
        ];
        for sign in [-1.0, 1.0] {
            let c = Center {
                atom: 0,
                ligands: [1, 2, 3, 4],
                sign,
            };
            // 界给得极宽,只留手性项
            let mut b = Bounds::new(5, 0.01, 1000.0);
            for i in 0..5 {
                for j in (i + 1)..5 {
                    b.set_lower(i, j, 0.01);
                    b.set_upper(i, j, 1000.0);
                }
            }
            let f = Field::new(&b, &[c]);
            let x = flat(&pts);
            // **两条都必须真的在罚**,分开断言 —— 只问"合起来非零"的话,
            // 四配体项在罚就够了,伞形项照样可以是零,而那正是要验的那一条。
            let mut g = vec![0.0; 15];
            let e0 = f.value_and_grad(&x, &mut g);
            assert!(e0 > 0.0, "sign={sign} 这个构型没违反手性,测不到东西");
            assert_both_chiral_active(&pts, &c);
            let e = max_grad_error(&f, &x, 1e-6);
            assert!(e < 1e-6, "sign={sign} 手性项梯度对不上:{e:.3e}");
        }
    }

    #[test]
    fn 两项一起的梯度也一致() {
        let pts = random_points(8, 21, 2.5);
        // 下标 0 是**中心原子**,1..=4 是配体。先前这里写的是
        // `atom: 0, ligands: [0,1,2,3]` —— 中心就是第一个配体,那是个不存在的构型。
        let c = Center {
            atom: 0,
            ligands: [1, 2, 3, 4],
            sign: -1.0,
        };
        let f = field_from(&pts, 0.2, &[c]);
        let mut st = 33u64;
        let x: Vec<f64> = flat(&pts).iter().map(|v| v + lcg(&mut st) * 1.5).collect();
        // **护栏:两条手性项都必须真的在罚**,否则这条名为"两项一起"的测试
        // 其实只测了距离项(或者只测了两条手性项里的一条)。
        // 注意断言要用**扰动之后**的 `x`,那才是 `max_grad_error` 实际吃的构型。
        let px: Vec<[f64; 3]> = (0..pts.len())
            .map(|i| [x[3 * i], x[3 * i + 1], x[3 * i + 2]])
            .collect();
        assert_both_chiral_active(&px, &c);
        let e = max_grad_error(&f, &x, 1e-6);
        assert!(e < 1e-6, "合起来梯度对不上:{e:.3e}");
    }

    #[test]
    fn 落在界内时罚为零() {
        let pts = random_points(12, 3, 4.0);
        let f = field_from(&pts, 0.3, &[]);
        let x = flat(&pts);
        let mut g = vec![0.0; x.len()];
        let e = f.value_and_grad(&x, &mut g);
        assert_eq!(e, 0.0, "界是照这组点定的,它本身应当零罚");
        for (k, v) in g.iter().enumerate() {
            assert_eq!(*v, 0.0, "第 {k} 个分量的梯度应当是 0");
        }
    }

    #[test]
    fn 精修能把违反压下去() {
        // 起点塌到原点附近,让几乎每一对都低于下限
        let pts = random_points(15, 13, 4.0);
        let f = field_from(&pts, 0.25, &[]);
        let mut st = 5u64;
        let mut x: Vec<f64> = (0..45).map(|_| lcg(&mut st) * 0.3).collect();
        let mut g = vec![0.0; 45];
        let e0 = f.value_and_grad(&x, &mut g);
        let r = minimize(
            &f,
            &mut x,
            &Options {
                max_iter: 2000,
                grad_tol: 1e-10,
                memory: 8,
            },
        );
        assert!(e0 > 1.0, "起点应当很差,实际 {e0:.3e}");
        assert!(
            r.value < 1e-10,
            "没压下去:{:.3e}(起点 {e0:.3e},迭代 {})",
            r.value,
            r.iterations
        );
    }

    #[test]
    fn 宽区间的对也必须进力场() {
        // **这一条守的是"全部 N² 对"那个决定。** RDKit 默认只放
        // `u − l ≤ basinThresh(5.0)` 的对,于是宽区间的对(拓扑上远的那些)
        // 在力场里**没有任何约束** —— 而那正是自穿会发生的地方。
        //
        // 前面几条测试都用 `真实距离 ± 0.3` 造界,区间宽只有 0.6,
        // 那道过滤一条都滤不掉 —— 所以它们**抓不住**这个退化(实测确实逃脱)。
        // 这里专门给一个宽到 19 的区间。
        let mut b = Bounds::new(2, 0.0, 1000.0);
        b.set_lower(0, 1, 1.0);
        b.set_upper(0, 1, 20.0);
        assert!(
            b.upper(0, 1) - b.lower(0, 1) > 5.0,
            "构造失效:区间必须宽过 basinThresh 才测得到东西"
        );
        let f = Field::new(&b, &[]);
        let mut g = vec![0.0; 6];
        // 拉到 30 Å —— 超了上限 10 Å
        let e = f.value_and_grad(&[0.0, 0.0, 0.0, 30.0, 0.0, 0.0], &mut g);
        assert!(
            e > 0.0,
            "宽区间的对被漏掉了 —— 那是 RDKit basinThresh 的行为,不是我们的"
        );
        assert!(g[0] != 0.0, "梯度也该有,实际是 {}", g[0]);
    }

    #[test]
    fn 下限那一支是饱和的() {
        // **这是个已知的性质,不是缺陷** —— 但必须钉住,因为它决定了
        // "两个原子穿过彼此"的势垒上限,而立体化学要靠别的判据守。
        //
        // d → 0 时 val → 2l²/l² − 1 = 1,所以单对的罚**最多 1**(权重 1 时)。
        let mut b = Bounds::new(2, 0.0, 1000.0);
        b.set_lower(0, 1, 2.0);
        b.set_upper(0, 1, 3.0);
        let f = Field::new(&b, &[]);
        let mut g = vec![0.0; 6];
        // 完全重合
        let e_zero = f.value_and_grad(&[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], &mut g);
        assert!(
            (e_zero - 1.0).abs() < 1e-12,
            "完全重叠时的罚应当恰好是 1,实得 {e_zero}"
        );
        // 而超上限那一支是**无界**的:拉到 30 Å
        let e_far = f.value_and_grad(&[0.0, 0.0, 0.0, 30.0, 0.0, 0.0], &mut g);
        assert!(e_far > 50.0, "超上限应当罚得很重,实得 {e_far}");
    }

    #[test]
    fn 手性项能把号翻回来() {
        // 给一个手性反了的四面体,让优化器把它翻正 —— 这正是
        // "个别中心错,三维精修够得着"那句话的判据。
        // 下标 0 是**中心原子**,1..=4 是四个配体 —— 先前这里写的是
        // `atom: 0, ligands: [0,1,2,3]`(中心就是第一个配体),那是个不存在的
        // 构型:`a = l₀ − c = 0`,梯度整体恒零,测试实际上什么都没测。
        let pts: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.2],
            [1.13, 0.0, -0.4],
            [-0.57, 0.98, -0.4],
            [-0.57, -0.98, -0.4],
        ];
        let c = Center {
            atom: 0,
            ligands: [1, 2, 3, 4],
            sign: -1.0,
        };
        // 这组点的中心基点体积是正的;`sign = -1` 要求它变成负的
        let v0 = crate::chiral::center_volume(&pts, &c);
        assert!(v0 > 0.0, "构造有问题:起始体积应当为正,实得 {v0}");
        let mut b = Bounds::new(5, 0.5, 1000.0);
        for i in 1..5 {
            // 中心到配体:键长那一档
            b.set_lower(0, i, 1.0);
            b.set_upper(0, i, 1.4);
            for j in (i + 1)..5 {
                b.set_lower(i, j, 1.6);
                b.set_upper(i, j, 4.0);
            }
        }
        let f = Field::new(&b, &[c]);
        let mut x = flat(&pts);
        let r = minimize(
            &f,
            &mut x,
            &Options {
                max_iter: 3000,
                grad_tol: 1e-9,
                memory: 8,
            },
        );
        let p: Vec<[f64; 3]> = (0..5)
            .map(|i| [x[3 * i], x[3 * i + 1], x[3 * i + 2]])
            .collect();
        let v1 = crate::chiral::center_volume(&p, &c);
        // 中心基点体积比 `UMBRELLA_LO`,四配体体积比 `VOL_LO` —— **两个尺度不一样**,
        // 别拿其中一个的阈值去量另一个(实测中心基点中位 2.573、四配体 8.374)。
        assert!(
            v1 <= -UMBRELLA_LO + 1e-6,
            "手性没翻过来:{v0:.3} → {v1:.3}(目标 ≤ −{UMBRELLA_LO},残值 {:.3e})",
            r.value
        );
        let vl = crate::chiral::signed_volume(p[1], p[2], p[3], p[4]);
        assert!(
            vl >= VOL_LO - 1e-6,
            "四配体那一项也该到位(号与中心基点相反):{vl:.3},目标 ≥ {VOL_LO}"
        );
    }

    #[test]
    fn 手性项挡得住伞形翻转() {
        // **这一条是先前完全没有的那档回复力。** 四个配体摆成正四面体不动,
        // 把中心原子推到四面体**外面** —— 四配体行列式看不见这件事,
        // 所以旧的手性项对它没有任何回复力。
        // 先摆一个正常的中心,再把中心原子**沿配体 1,2,3 所在平面镜像**过去。
        // `V` 正比于中心到那张平面的有号距离,所以镜像**精确变号** —— 这是
        // 算出来的构造,不是试出来的。
        let mut pts: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.2],
            [1.13, 0.0, -0.4],
            [-0.57, 0.98, -0.4],
            [-0.57, -0.98, -0.4],
        ];
        let c = Center {
            atom: 0,
            ligands: [1, 2, 3, 4],
            sign: 1.0,
        };
        let good = crate::chiral::center_volume(&pts, &c);
        {
            let (p, q, r) = (pts[1], pts[2], pts[3]);
            let sub = |u: [f64; 3], v: [f64; 3]| [u[0] - v[0], u[1] - v[1], u[2] - v[2]];
            let (a, b) = (sub(q, p), sub(r, p));
            let n = [
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0],
            ];
            let nn = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
            let d = sub(pts[0], p);
            let t = (d[0] * n[0] + d[1] * n[1] + d[2] * n[2]) / nn;
            for k in 0..3 {
                pts[0][k] -= 2.0 * t * n[k];
            }
        }
        let v0 = crate::chiral::center_volume(&pts, &c);
        assert!(good > 0.0, "起始那个正常构型该是正的,实得 {good}");
        assert!(v0 < 0.0, "构造有问题:翻伞之后体积该为负,实得 {v0}");
        // 只给距离上的活动余地,不直接钉中心的位置 —— 要让**手性项**把它拉回去
        let mut b = Bounds::new(5, 0.5, 1000.0);
        for i in 1..5 {
            b.set_lower(0, i, 1.0);
            b.set_upper(0, i, 1.6);
            for j in (i + 1)..5 {
                b.set_lower(i, j, 1.6);
                b.set_upper(i, j, 4.0);
            }
        }
        let f = Field::new(&b, &[c]);
        let mut x = flat(&pts);
        minimize(
            &f,
            &mut x,
            &Options {
                max_iter: 3000,
                grad_tol: 1e-9,
                memory: 8,
            },
        );
        let p: Vec<[f64; 3]> = (0..5)
            .map(|i| [x[3 * i], x[3 * i + 1], x[3 * i + 2]])
            .collect();
        let v1 = crate::chiral::center_volume(&p, &c);
        assert!(v1 > 0.0, "伞形翻转没被拉回来:{v0:.3} → {v1:.3}");
    }
}
