//! **距离几何的误差函数** —— 精修阶段要最小化的那个目标。
//!
//! 两项:原子对的距离违反,以及手性中心的有符号体积违反。
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
//! # 手性项:有符号体积落进区间
//!
//! `E += w·(V − lo)²`(`V < lo` 时)或 `w·(V − hi)²`(`V > hi` 时)。
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

/// 手性体积的目标下限(绝对值)。与 RDKit 的 `volLowerBound` 同值。
///
/// **它同时是一道"别压平"的闸**:只要求符号对的话,把中心压成近乎共面、
/// 体积 `+1e−6`,符号照样对而分子是废的。要求 `|V| ≥ 5` 就堵住了这条路。
pub const VOL_LO: f64 = 5.0;

/// 手性体积的目标上限(绝对值)。与 RDKit 的 `volUpperBound` 同值。
pub const VOL_HI: f64 = 100.0;

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
    #[must_use]
    pub fn new(b: &Bounds, centers: &[Center]) -> Self {
        let n = b.len();
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
                .ligands
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
            let (a, b, cc) = (sub(p[1], p[0]), sub(p[2], p[0]), sub(p[3], p[0]));
            let bxc = cross(b, cc);
            let v = a[0] * bxc[0] + a[1] * bxc[1] + a[2] * bxc[2];

            // 目标区间:号由 `Center::sign` 定,大小是 [VOL_LO, VOL_HI]
            let (lo, hi) = if c.sign < 0.0 {
                (-VOL_HI, -VOL_LO)
            } else {
                (VOL_LO, VOL_HI)
            };
            let dev = if v < lo {
                v - lo
            } else if v > hi {
                v - hi
            } else {
                continue;
            };
            e += self.weight_chiral * dev * dev;
            // dE/dV = 2·w·dev。**因子 2 不能少** —— RDKit 那边就少了它,
            // 于是它的手性项实际权重只有名义值的一半。
            let k = 2.0 * self.weight_chiral * dev;
            // V = a·(b×c) 对四个点的导数
            let cxa = cross(cc, a);
            let axb = cross(a, b);
            for t in 0..3 {
                let (l1, l2, l3, l0) = (
                    c.ligands[1] as usize,
                    c.ligands[2] as usize,
                    c.ligands[3] as usize,
                    c.ligands[0] as usize,
                );
                g[3 * l1 + t] += k * bxc[t];
                g[3 * l2 + t] += k * cxa[t];
                g[3 * l3 + t] += k * axb[t];
                g[3 * l0 + t] -= k * (bxc[t] + cxa[t] + axb[t]);
            }
        }
        e
    }
}

impl Objective for Field {
    fn value_and_grad(&self, x: &[f64], grad: &mut [f64]) -> f64 {
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

    #[test]
    fn 手性项的梯度与能量一致() {
        let pts = random_points(4, 9, 1.5);
        for sign in [-1.0, 1.0] {
            let c = Center {
                atom: 0,
                ligands: [0, 1, 2, 3],
                sign,
            };
            // 界给得极宽,只留手性项
            let mut b = Bounds::new(4, 0.01, 1000.0);
            for i in 0..4 {
                for j in (i + 1)..4 {
                    b.set_lower(i, j, 0.01);
                    b.set_upper(i, j, 1000.0);
                }
            }
            let f = Field::new(&b, &[c]);
            let x = flat(&pts);
            // 先确认这个构型真的违反了(否则罚项为 0,梯度恒 0,测了个寂寞)
            let mut g = vec![0.0; 12];
            let e0 = f.value_and_grad(&x, &mut g);
            assert!(e0 > 0.0, "sign={sign} 这个构型没违反手性,测不到东西");
            let e = max_grad_error(&f, &x, 1e-6);
            assert!(e < 1e-6, "sign={sign} 手性项梯度对不上:{e:.3e}");
        }
    }

    #[test]
    fn 两项一起的梯度也一致() {
        let pts = random_points(8, 21, 2.5);
        let c = Center {
            atom: 0,
            ligands: [0, 1, 2, 3],
            sign: 1.0,
        };
        let f = field_from(&pts, 0.2, &[c]);
        let mut st = 33u64;
        let x: Vec<f64> = flat(&pts).iter().map(|v| v + lcg(&mut st) * 1.5).collect();
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
        let pts: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 1.2],
            [1.13, 0.0, -0.4],
            [-0.57, 0.98, -0.4],
            [-0.57, -0.98, -0.4],
        ];
        // 这组点的有符号体积是负的(Ccw);要求它变成正的(Cw)
        let v0 = crate::chiral::signed_volume(pts[0], pts[1], pts[2], pts[3]);
        assert!(v0 < 0.0, "构造有问题:起始体积应当为负,实得 {v0}");
        let c = Center {
            atom: 0,
            ligands: [0, 1, 2, 3],
            sign: 1.0,
        };
        let mut b = Bounds::new(4, 0.5, 1000.0);
        for i in 0..4 {
            for j in (i + 1)..4 {
                b.set_lower(i, j, 1.0);
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
        let p: Vec<[f64; 3]> = (0..4)
            .map(|i| [x[3 * i], x[3 * i + 1], x[3 * i + 2]])
            .collect();
        let v1 = crate::chiral::signed_volume(p[0], p[1], p[2], p[3]);
        assert!(
            v1 >= VOL_LO - 1e-6,
            "手性没翻过来:{v0:.3} → {v1:.3}(目标 ≥ {VOL_LO},残值 {:.3e})",
            r.value
        );
    }
}
