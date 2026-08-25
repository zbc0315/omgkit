//! **L-BFGS** —— 有限内存的拟牛顿下降,手写,零外部依赖。
//!
//! # 这个目标函数有三个坑,不处理就不收敛
//!
//! 距离几何的误差函数是**平底罚项**:原子对落在 `[l, u]` 里罚 0、出界才罚,
//! 于是它只有 `C¹` —— 每有一对原子跨过自己的边界,二阶导就跳变一次。
//! 三个后果,逐条对付:
//!
//! 1. **曲率信息会被污染。** L-BFGS 靠 `(s, y)` 对近似 Hessian,而 `y·s ≤ 0` 的对
//!    意味着"沿这个方向曲率为负",用它更新出来的方向不是下降方向。
//! 2. **所以线搜索必须上强 Wolfe,不能只用 Armijo 回溯。**
//!    Armijo 只要求"降得够多",不管曲率,于是它会接受**跨过谷底**的步长 ——
//!    那一步梯度反号、`y·s < 0`,曲率对被丢掉;丢光之后 L-BFGS 退化成最速下降。
//!
//!    **这是实测出来的,不是推的。** 头一版就写的 Armijo + 丢弃,8 维 Rosenbrock 上:
//!
//!    | | Armijo + 丢弃 | 强 Wolfe |
//!    |---|---|---|
//!    | Rosenbrock(8 维) | **5000 步没收敛**,丢 4901 个对、回溯 73882 次 | **71 步**,丢 0,回溯 163 |
//!    | Powell | 1112 步,丢 1037 | **152 步**,丢 72 |
//!
//!    `y·s` 从第 3 步起一直是负的、步长冻在 3e-6 原地打转。
//!    强 Wolfe 的曲率条件 `|g₊·d| ≤ c₂|g·d|` **保证** `y·s > 0`:
//!    `y·s = α(g₊−g)·d ≥ α(c₂−1)(g·d) > 0`。
//!
//!    **谨慎更新因此从主力退成保险。** 实测把它整个关掉:Rosenbrock 与距离盒
//!    **一步都不变**,Powell 反而从 152 步降到 103 步 —— 也就是说它在这几个目标上
//!    从不该触发,偶尔触发时还是在拖后腿。
//!
//!    留着它是**拿一点速度换一层保险**:强 Wolfe 的保证是在精确算术下成立的,
//!    而 `y·s` 与 `‖s‖²` 同时掉到舍入量级时(Hessian 奇异那类),
//!    一个数值上的垃圾对照样能混进历史。这笔账是量过的,不是想当然。
//!
//!    **别写判据去逼它触发** —— 那等于断言一件不该发生的事。同理"方向不下降就
//!    重置历史"那一段:实测也从不触发。
//! 3. **求和顺序必须固定。** 目标函数是 `N²` 项求和,浮点加法不结合,
//!    顺序一变结果就变。所以**不许并行归约** —— 这个 crate 承诺同一输入同一输出。
//!
//! # 判官
//!
//! 标准测试函数(Rosenbrock、Beale、Powell)有解析最优解,不依赖任何外部实现。
//! **但它们都是光滑的,验不出上面那三条。** 所以另配一个**分段二次的平底罚项**
//! 测试函数 —— 那才是真正要跑的那一类。见本文件末尾的测试。
//!
//! 另有 [`max_grad_error`]:中心差分数值梯度。**每一个进优化器的项都必须过这一条**,
//! 因为 L-BFGS 对"能量与梯度不是同一个函数"零容忍(线搜索用能量、方向用梯度,
//! 两者不一致时曲率对全是垃圾)。RDKit 的手性项与第四维项就**差一个因子 2**
//! (能量 `w(v−lo)²`、梯度 `w(v−lo)`),照抄会静默改掉权重比。

/// 一个可微的目标函数。
pub trait Objective {
    /// 同时给出能量与梯度。
    ///
    /// **必须是同一个函数的能量与它的导数。** 不一致的话 L-BFGS 的线搜索与
    /// 方向计算各说各话,曲率对全是垃圾 —— 表现为"就是不收敛",而且极难查。
    /// 用 [`max_grad_error`] 把这一条钉进单元测试。
    fn value_and_grad(&self, x: &[f64], grad: &mut [f64]) -> f64;
}

/// 优化的旋钮。
#[derive(Debug, Clone, Copy)]
pub struct Options {
    /// 最多迭代多少次。
    pub max_iter: usize,
    /// 梯度的无穷范数降到这个数以下就算收敛。
    pub grad_tol: f64,
    /// 记忆多少个 `(s, y)` 对。5~10 是常用范围。
    pub memory: usize,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            max_iter: 400,
            grad_tol: 1e-6,
            memory: 8,
        }
    }
}

/// 优化跑完之后的账。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Report {
    /// 迭代了多少次。
    pub iterations: usize,
    /// 最终的目标值。
    pub value: f64,
    /// 最终梯度的无穷范数。
    pub grad_norm: f64,
    /// 有没有达到 [`Options::grad_tol`]。
    pub converged: bool,
    /// **谨慎更新丢弃了多少个 `(s, y)` 对。**
    ///
    /// 这个数是**线搜索是否健康**的信号,不是目标难度的信号:
    /// 强 Wolfe 保证 `y·s > 0`,所以它正常情况下应当是 0 或很小。
    /// **它一大就说明线搜索出了问题** —— 头一版用 Armijo 回溯时,
    /// 8 维 Rosenbrock 上这个数是 4901/5000,而那正是它跑不动的原因。
    pub discarded: usize,
    /// 线搜索总共回溯了多少次。
    pub backtracks: usize,
}

/// 中心差分数值梯度与解析梯度的**最大相对偏差**。
///
/// 给单元测试用 —— 每一个进优化器的项都必须过这一条,理由见模块文档。
///
/// 相对化的分母取 `max(|数值|, |解析|, 1)`,免得在梯度接近 0 的分量上炸出假警报。
#[must_use]
pub fn max_grad_error(obj: &dyn Objective, x: &[f64], h: f64) -> f64 {
    let n = x.len();
    let mut g = vec![0.0; n];
    obj.value_and_grad(x, &mut g);
    let mut probe = x.to_vec();
    let mut worst: f64 = 0.0;
    let mut scratch = vec![0.0; n];
    for i in 0..n {
        let orig = probe[i];
        probe[i] = orig + h;
        let f1 = obj.value_and_grad(&probe, &mut scratch);
        probe[i] = orig - h;
        let f0 = obj.value_and_grad(&probe, &mut scratch);
        probe[i] = orig;
        let num = (f1 - f0) / (2.0 * h);
        let denom = num.abs().max(g[i].abs()).max(1.0);
        let rel = (num - g[i]).abs() / denom;
        // 归约必须让 NaN 赢,否则"解析梯度是 NaN"会被判成偏差 0 ——
        // 见 [`crate::linalg::max_nan_wins`]。
        worst = crate::linalg::max_nan_wins(worst, rel);
    }
    worst
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    // **固定顺序求和。** 不许换成并行归约 —— 浮点加法不结合,顺序一变答案就变。
    let mut s = 0.0;
    for k in 0..a.len() {
        s += a[k] * b[k];
    }
    s
}

/// 无穷范数。**任何一个分量是 NaN,结果就是 NaN。**
///
/// 不能写成 `fold(0.0, f64::max)`:`f64::max` 碰上 NaN 时**返回另一个操作数**
/// (两个方向都如此,实测 `0.0f64.max(NAN)` 与 `NAN.max(0.0)` 都是 `0`),
/// 于是整条梯度全是 NaN 时这里给出 **0** —— 恰好就是"梯度为零、已经收敛"的
/// 意思。收敛判据 `grad_norm <= grad_tol` 于是当场成立,一个废掉的结构
/// 拿到满分。传播 NaN 之后那个比较恒为 false,不收敛就是不收敛。
fn inf_norm(v: &[f64]) -> f64 {
    v.iter()
        .fold(0.0_f64, |m, x| crate::linalg::max_nan_wins(m, x.abs()))
}

/// 强 Wolfe 线搜索的结果。
struct LineSearch {
    ok: bool,
    /// 接受点的目标值
    f: f64,
    /// 一共评估了多少次(诊断用)
    evals: usize,
}

/// **强 Wolfe 线搜索**(Nocedal & Wright 算法 3.5 / 3.6:先扩张括号,再 zoom)。
///
/// 两个条件:
/// - 充分下降(Armijo):`f(α) ≤ f₀ + c₁·α·(g·d)`
/// - 曲率(强):`|g(α)·d| ≤ c₂·|g·d|`
///
/// 第二条正是 Armijo 缺的那条,也是 `y·s > 0` 的保证 —— 理由写在调用处。
///
/// zoom 段用**二分**而不是三次插值:二分完全确定、没有除零分支,
/// 而这个 crate 承诺同一输入同一输出。代价是多几次函数求值,可以忍。
fn wolfe_search(
    obj: &dyn Objective,
    x: &[f64],
    dir: &[f64],
    f0: f64,
    slope0: f64,
    x_out: &mut [f64],
    g_out: &mut [f64],
) -> LineSearch {
    const C1: f64 = 1e-4;
    const C2: f64 = 0.9;
    const MAX_EVAL: usize = 60;
    let n = x.len();
    let mut evals = 0;

    let eval = |a: f64, xo: &mut [f64], go: &mut [f64]| -> (f64, f64) {
        for j in 0..n {
            xo[j] = x[j] + a * dir[j];
        }
        let fa = obj.value_and_grad(xo, go);
        (fa, dot(go, dir))
    };

    let (mut a_prev, mut f_prev) = (0.0_f64, f0);
    let mut a = 1.0_f64;
    // 括号阶段
    let (mut lo, mut hi) = (0.0_f64, -1.0_f64); // hi < 0 表示还没括起来
    for i in 0..MAX_EVAL {
        let (fa, sa) = eval(a, x_out, g_out);
        evals += 1;
        if fa > f0 + C1 * a * slope0 || (i > 0 && fa >= f_prev) {
            lo = a_prev;
            hi = a;
            break;
        }
        if sa.abs() <= -C2 * slope0 {
            return LineSearch {
                ok: true,
                f: fa,
                evals,
            }; // 两条都满足,直接收
        }
        if sa >= 0.0 {
            lo = a;
            hi = a_prev;
            break;
        }
        a_prev = a;
        f_prev = fa;
        a *= 2.0;
        if a > 1e10 {
            // 一路降下去没有极小 —— 目标在这个方向上无界
            return LineSearch {
                ok: false,
                f: fa,
                evals,
            };
        }
    }
    if hi < 0.0 {
        // 括号阶段用满了还没括住
        return LineSearch {
            ok: false,
            f: f0,
            evals,
        };
    }

    // zoom 阶段:二分
    let mut f_lo = {
        let (fl, _) = eval(lo, x_out, g_out);
        evals += 1;
        fl
    };
    while evals < MAX_EVAL {
        let am = 0.5 * (lo + hi);
        let (fa, sa) = eval(am, x_out, g_out);
        evals += 1;
        if fa > f0 + C1 * am * slope0 || fa >= f_lo {
            hi = am;
        } else {
            if sa.abs() <= -C2 * slope0 {
                return LineSearch {
                    ok: true,
                    f: fa,
                    evals,
                };
            }
            if sa * (hi - lo) >= 0.0 {
                hi = lo;
            }
            lo = am;
            f_lo = fa;
        }
        if (hi - lo).abs() < 1e-16 {
            break;
        }
    }
    // 用满了:退回区间左端(它满足充分下降),只是曲率条件没达标
    let (fa, _) = eval(lo, x_out, g_out);
    evals += 1;
    LineSearch {
        ok: lo > 0.0 && fa < f0,
        f: fa,
        evals,
    }
}

/// 最小化 `obj`,`x` 既是初值也是出口。
///
/// # Panics
///
/// `x` 为空时 panic。
pub fn minimize(obj: &dyn Objective, x: &mut [f64], opts: &Options) -> Report {
    assert!(!x.is_empty(), "没有变量可优化");
    let n = x.len();
    let m = opts.memory.max(1);

    let mut g = vec![0.0; n];
    let mut f = obj.value_and_grad(x, &mut g);
    let mut report = Report {
        iterations: 0,
        value: f,
        grad_norm: inf_norm(&g),
        converged: inf_norm(&g) <= opts.grad_tol,
        discarded: 0,
        backtracks: 0,
    };
    if report.converged {
        return report;
    }

    // 环形存 (s, y) 对
    let mut s_hist: Vec<Vec<f64>> = Vec::with_capacity(m);
    let mut y_hist: Vec<Vec<f64>> = Vec::with_capacity(m);
    let mut rho: Vec<f64> = Vec::with_capacity(m);

    let mut dir = vec![0.0; n];
    let mut x_new = vec![0.0; n];
    let mut g_new = vec![0.0; n];
    let mut alpha = vec![0.0; m];

    for iter in 1..=opts.max_iter {
        // ---- 两循环递推,算出下降方向 ----
        dir.copy_from_slice(&g);
        let k = s_hist.len();
        for i in (0..k).rev() {
            let a = rho[i] * dot(&s_hist[i], &dir);
            alpha[i] = a;
            for j in 0..n {
                dir[j] -= a * y_hist[i][j];
            }
        }
        // 初始 Hessian 的缩放 γ = (sᵀy)/(yᵀy);第一步没有历史,取 1
        let gamma = if k > 0 {
            let ys = dot(&s_hist[k - 1], &y_hist[k - 1]);
            let yy = dot(&y_hist[k - 1], &y_hist[k - 1]);
            if yy > 0.0 {
                ys / yy
            } else {
                1.0
            }
        } else {
            1.0
        };
        for d in dir.iter_mut() {
            *d *= gamma;
        }
        for i in 0..k {
            let beta = rho[i] * dot(&y_hist[i], &dir);
            let coef = alpha[i] - beta;
            for j in 0..n {
                dir[j] += coef * s_hist[i][j];
            }
        }
        for d in dir.iter_mut() {
            *d = -*d;
        }

        // 方向必须真的下降。不下降就退回最速下降 —— 这是兜底,不是常态。
        let mut slope = dot(&g, &dir);
        if slope >= 0.0 {
            for (d, gi) in dir.iter_mut().zip(&g) {
                *d = -gi;
            }
            slope = -dot(&g, &g);
            s_hist.clear();
            y_hist.clear();
            rho.clear();
        }

        // ---- 强 Wolfe 线搜索 ----
        //
        // **这里不能只用 Armijo 回溯。** Armijo 只要求"降得够多",不管曲率,
        // 于是它会接受**跨过谷底**的步长 —— 那一步上梯度反号,`y·s < 0`,
        // 曲率对被谨慎更新丢掉;丢光之后 L-BFGS 退化成最速下降,而最速下降
        // 的首试步长 1.0 在 Rosenbrock 这种病态谷里又要回溯十几次。
        //
        // 实测(8 维 Rosenbrock):Armijo 版跑满 5000 步没收敛,丢弃 4901 个曲率对、
        // 回溯 73882 次,`y·s` 从第 3 步起一直是负的、步长冻在 3e-6 原地打转。
        //
        // 强 Wolfe 的曲率条件 `|g₊·d| ≤ c₂|g·d|` **保证** `y·s > 0`:
        // `y·s = α(g₊−g)·d ≥ α(c₂−1)(g·d) > 0`(因为 `g·d < 0`、`c₂ < 1`)。
        // 也就是说曲率对天然合格,谨慎更新退化成一道保险而不是主力。
        let ls = wolfe_search(obj, x, &dir, f, slope, &mut x_new, &mut g_new);
        report.backtracks += ls.evals;
        let ok = ls.ok;
        if ok {
            f = ls.f;
        }
        if !ok {
            // 步长缩到底还降不下去:到极小了,或者目标不连续。停,并如实报出来。
            report.iterations = iter;
            report.value = f;
            report.grad_norm = inf_norm(&g);
            report.converged = report.grad_norm <= opts.grad_tol;
            return report;
        }

        // ---- 谨慎更新:曲率不正的 (s, y) 对**丢掉** ----
        let mut s_new = vec![0.0; n];
        let mut y_new = vec![0.0; n];
        for j in 0..n {
            s_new[j] = x_new[j] - x[j];
            y_new[j] = g_new[j] - g[j];
        }
        let ys = dot(&s_new, &y_new);
        let ss = dot(&s_new, &s_new);
        // Li–Fukushima 的谨慎条件。平底罚项上这一条会**经常**触发,那是正常的 ——
        // 边界跳变造出来的曲率信息本来就不该用。
        const CAUTION: f64 = 1e-8;
        if ys > CAUTION * ss && ss > 0.0 {
            if s_hist.len() == m {
                s_hist.remove(0);
                y_hist.remove(0);
                rho.remove(0);
            }
            rho.push(1.0 / ys);
            s_hist.push(s_new);
            y_hist.push(y_new);
        } else {
            report.discarded += 1;
        }

        x.copy_from_slice(&x_new);
        g.copy_from_slice(&g_new);
        report.iterations = iter;
        report.value = f;
        report.grad_norm = inf_norm(&g);
        if report.grad_norm <= opts.grad_tol {
            report.converged = true;
            return report;
        }
    }
    report
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Rosenbrock;
    impl Objective for Rosenbrock {
        fn value_and_grad(&self, x: &[f64], g: &mut [f64]) -> f64 {
            let mut f = 0.0;
            for v in g.iter_mut() {
                *v = 0.0;
            }
            for i in 0..(x.len() - 1) {
                let (a, b) = (x[i], x[i + 1]);
                let t1 = b - a * a;
                let t2 = 1.0 - a;
                f += 100.0 * t1 * t1 + t2 * t2;
                g[i] += -400.0 * a * t1 - 2.0 * t2;
                g[i + 1] += 200.0 * t1;
            }
            f
        }
    }

    struct Beale;
    impl Objective for Beale {
        fn value_and_grad(&self, x: &[f64], g: &mut [f64]) -> f64 {
            let (a, b) = (x[0], x[1]);
            let f1 = 1.5 - a + a * b;
            let f2 = 2.25 - a + a * b * b;
            let f3 = 2.625 - a + a * b * b * b;
            g[0] = 2.0 * (f1 * (b - 1.0) + f2 * (b * b - 1.0) + f3 * (b * b * b - 1.0));
            g[1] = 2.0 * (f1 * a + f2 * 2.0 * a * b + f3 * 3.0 * a * b * b);
            f1 * f1 + f2 * f2 + f3 * f3
        }
    }

    struct Powell;
    impl Objective for Powell {
        fn value_and_grad(&self, x: &[f64], g: &mut [f64]) -> f64 {
            let (a, b, c, d) = (x[0], x[1], x[2], x[3]);
            let t1 = a + 10.0 * b;
            let t2 = c - d;
            let t3 = b - 2.0 * c;
            let t4 = a - d;
            g[0] = 2.0 * t1 + 40.0 * t4 * t4 * t4;
            g[1] = 20.0 * t1 + 4.0 * t3 * t3 * t3;
            g[2] = 10.0 * t2 - 8.0 * t3 * t3 * t3;
            g[3] = -10.0 * t2 - 40.0 * t4 * t4 * t4;
            t1 * t1 + 5.0 * t2 * t2 + t3 * t3 * t3 * t3 + 10.0 * t4 * t4 * t4 * t4
        }
    }

    /// **平底罚项** —— 这才是真正要跑的那一类。
    ///
    /// `Σ max(0, |xᵢ − cᵢ| − w)²`:落在盒子里罚 0、出界才罚。只有 `C¹`,
    /// 每跨一次边界二阶导就跳变一次 —— 前面三个光滑函数验不出这件事。
    struct FlatBottom {
        center: Vec<f64>,
        half_width: f64,
    }
    impl Objective for FlatBottom {
        fn value_and_grad(&self, x: &[f64], g: &mut [f64]) -> f64 {
            let mut f = 0.0;
            for i in 0..x.len() {
                let d = x[i] - self.center[i];
                let over = d.abs() - self.half_width;
                if over > 0.0 {
                    f += over * over;
                    g[i] = 2.0 * over * d.signum();
                } else {
                    g[i] = 0.0;
                }
            }
            f
        }
    }

    #[test]
    fn 解析梯度与数值梯度一致() {
        // **每一个目标都必须过这一条。** L-BFGS 对"能量与梯度不是同一个函数"
        // 零容忍 —— 而 RDKit 的手性项与第四维项恰好差一个因子 2,
        // 照抄那种写法就会栽在这里。
        let flat = FlatBottom {
            center: vec![0.0, 0.0, 0.0],
            half_width: 1.0,
        };
        let cases: Vec<(&dyn Objective, Vec<f64>)> = vec![
            (&Rosenbrock, vec![-1.2, 1.0, 0.7, -0.3]),
            (&Beale, vec![1.0, 0.5]),
            (&Powell, vec![3.0, -1.0, 0.0, 1.0]),
            (&flat, vec![2.5, -3.0, 0.25]),
        ];
        for (obj, x) in cases {
            let e = max_grad_error(obj, &x, 1e-6);
            assert!(e < 1e-6, "梯度与能量对不上,最大相对偏差 {e:.3e}");
        }
    }

    #[test]
    fn rosenbrock_收敛到已知最优() {
        for n in [2usize, 10] {
            let mut x = vec![-1.2; n];
            for (i, v) in x.iter_mut().enumerate() {
                if i % 2 == 1 {
                    *v = 1.0;
                }
            }
            let r = minimize(
                &Rosenbrock,
                &mut x,
                &Options {
                    max_iter: 2000,
                    grad_tol: 1e-8,
                    memory: 8,
                },
            );
            assert!(r.value < 1e-12, "n={n} 目标值 {:.3e}", r.value);
            for (i, v) in x.iter().enumerate() {
                assert!((v - 1.0).abs() < 1e-5, "n={n} x[{i}] = {v}");
            }
        }
    }

    #[test]
    fn beale_与_powell_收敛到已知最优() {
        let mut x = vec![1.0, 1.0];
        let r = minimize(&Beale, &mut x, &Options::default());
        assert!(r.value < 1e-12, "Beale 目标值 {:.3e}", r.value);
        assert!(
            (x[0] - 3.0).abs() < 1e-4 && (x[1] - 0.5).abs() < 1e-4,
            "{x:?}"
        );

        // Powell 的 Hessian 在最优处奇异,收敛慢 —— 判据放在目标值上,不放在坐标上
        let mut x = vec![3.0, -1.0, 0.0, 1.0];
        let r = minimize(
            &Powell,
            &mut x,
            &Options {
                max_iter: 5000,
                grad_tol: 1e-12,
                memory: 8,
            },
        );
        assert!(r.value < 1e-8, "Powell 目标值 {:.3e}", r.value);
    }

    /// **耦合的**平底罚项:一堆点,逐对距离要落进 `[lo, hi]`。
    ///
    /// 这就是距离几何误差函数的缩小版 —— 上面那个可分离的 `FlatBottom` 太容易了:
    /// 每个坐标各走各的,进了盒子梯度就是 0,**从来不会造出坏曲率**,
    /// 于是谨慎更新那一段一次都触发不到(实测 `discarded` 恒为 0)。
    /// 真正的目标是**耦合**的:动一个原子会同时改变它与所有其他原子的距离,
    /// 边界跨来跨去,坏曲率才出得来。
    struct DistanceBox {
        n: usize,
        lo: Vec<f64>,
        hi: Vec<f64>,
    }
    impl DistanceBox {
        fn idx(&self, i: usize, j: usize) -> usize {
            i * self.n + j
        }
    }
    impl Objective for DistanceBox {
        fn value_and_grad(&self, x: &[f64], g: &mut [f64]) -> f64 {
            for v in g.iter_mut() {
                *v = 0.0;
            }
            let mut f = 0.0;
            for i in 0..self.n {
                for j in (i + 1)..self.n {
                    let d3 = [
                        x[3 * i] - x[3 * j],
                        x[3 * i + 1] - x[3 * j + 1],
                        x[3 * i + 2] - x[3 * j + 2],
                    ];
                    let d = (d3[0] * d3[0] + d3[1] * d3[1] + d3[2] * d3[2]).sqrt();
                    if d < 1e-12 {
                        continue;
                    }
                    let k = self.idx(i, j);
                    let over = d - self.hi[k];
                    let under = self.lo[k] - d;
                    let (pen, sign) = if over > 0.0 {
                        (over, 1.0)
                    } else if under > 0.0 {
                        (under, -1.0)
                    } else {
                        continue;
                    };
                    f += pen * pen;
                    let c = 2.0 * pen * sign / d;
                    for t in 0..3 {
                        g[3 * i + t] += c * d3[t];
                        g[3 * j + t] -= c * d3[t];
                    }
                }
            }
            f
        }
    }

    /// 造一个**有解**的实例:先随机取点,再把界定成"真实距离 ± 半宽"。
    /// 原来那组点就是一个零点,所以最优值确定是 0。
    fn distance_box(n: usize, half: f64, seed: u64) -> (DistanceBox, Vec<f64>) {
        let mut st = seed;
        let mut lcg = || {
            st = st.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            ((st >> 11) as f64) / ((1u64 << 53) as f64) * 2.0 - 1.0
        };
        let pts: Vec<f64> = (0..3 * n).map(|_| lcg() * 5.0).collect();
        let mut lo = vec![0.0; n * n];
        let mut hi = vec![0.0; n * n];
        for i in 0..n {
            for j in (i + 1)..n {
                let d = ((pts[3 * i] - pts[3 * j]).powi(2)
                    + (pts[3 * i + 1] - pts[3 * j + 1]).powi(2)
                    + (pts[3 * i + 2] - pts[3 * j + 2]).powi(2))
                .sqrt();
                lo[i * n + j] = (d - half).max(0.1);
                hi[i * n + j] = d + half;
            }
        }
        (DistanceBox { n, lo, hi }, pts)
    }

    #[test]
    fn 耦合平底罚项_梯度一致且线搜索健康() {
        let (obj, pts) = distance_box(12, 0.25, 7);
        // 先钉梯度 —— 起点故意打乱,让一堆对处在出界状态
        let mut st = 99u64;
        let mut lcg = || {
            st = st.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            ((st >> 11) as f64) / ((1u64 << 53) as f64) * 2.0 - 1.0
        };
        // **起点要真的差**:把所有点塌到原点附近。轻微扰动的起点收敛太快
        // (实测 12 步就到零),边界几乎不跨,测不出坏曲率。
        let start: Vec<f64> = (0..pts.len()).map(|_| lcg() * 0.3).collect();
        let e = max_grad_error(&obj, &start, 1e-6);
        assert!(e < 1e-6, "耦合平底罚项的梯度对不上:{e:.3e}");

        let mut x = start.clone();
        let r = minimize(
            &obj,
            &mut x,
            &Options {
                max_iter: 3000,
                grad_tol: 1e-9,
                memory: 8,
            },
        );
        // 界是照一组真实点定的,所以零点存在
        assert!(
            r.value < 1e-12,
            "没压下去:{:.3e}(迭代 {})",
            r.value,
            r.iterations
        );
        // **反过来断言**:线搜索健康时,丢弃应当很少。
        //
        // 头一版这里断言的是 `discarded > 0`("谨慎更新必须起作用"),那是错的 ——
        // 它把保险当成了主力。强 Wolfe 保证 `y·s > 0`,所以丢弃**本来就该接近 0**;
        // 这个数一大恰恰说明线搜索坏了(Armijo 版在 Rosenbrock 上是 4901/5000)。
        assert!(
            r.discarded * 4 < r.iterations.max(4),
            "丢了 {} 个曲率对 / 迭代 {} —— 线搜索可能退化了",
            r.discarded,
            r.iterations
        );
    }

    #[test]
    fn 平底罚项能压到零() {
        // 起点远在盒子外,最优是"进盒子",目标值 0。
        let n = 20;
        let obj = FlatBottom {
            center: (0..n).map(|i| (i as f64) * 0.1).collect(),
            half_width: 0.5,
        };
        let mut x: Vec<f64> = (0..n).map(|i| 10.0 + (i as f64)).collect();
        let r = minimize(
            &obj,
            &mut x,
            &Options {
                max_iter: 2000,
                grad_tol: 1e-10,
                memory: 8,
            },
        );
        assert!(r.value < 1e-16, "没压到零:{:.3e}", r.value);
        for (i, v) in x.iter().enumerate() {
            let d = (v - obj.center[i]).abs();
            assert!(d <= obj.half_width + 1e-6, "x[{i}] 没进盒子:偏 {d}");
        }
        // 注意:**可分离**的平底目标上 `discarded` 恒为 0,那是正常的 ——
        // 每个坐标各走各的,造不出坏曲率。谨慎更新那一段由上面那条
        // **耦合**的测试来守。
    }

    #[test]
    #[ignore]
    fn 探针_各目标丢弃了多少曲率对() {
        let mut x = vec![-1.2, 1.0, 0.7, -0.3, 2.0, -1.5, 0.2, 3.0];
        let r = minimize(
            &Rosenbrock,
            &mut x,
            &Options {
                max_iter: 5000,
                grad_tol: 1e-10,
                memory: 8,
            },
        );
        println!(
            "Rosenbrock(8维): 迭代 {} 丢弃 {} 回溯 {}",
            r.iterations, r.discarded, r.backtracks
        );
        let mut x = vec![3.0, -1.0, 0.0, 1.0];
        let r = minimize(
            &Powell,
            &mut x,
            &Options {
                max_iter: 5000,
                grad_tol: 1e-12,
                memory: 8,
            },
        );
        println!(
            "Powell:        迭代 {} 丢弃 {} 回溯 {}",
            r.iterations, r.discarded, r.backtracks
        );
        for (n, half, seed) in [(12usize, 0.25, 7u64), (30, 0.1, 11), (40, 0.05, 3)] {
            let (obj, _) = distance_box(n, half, seed);
            let mut st = seed ^ 0xabcd;
            let mut lcg = || {
                st = st.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
                ((st >> 11) as f64) / ((1u64 << 53) as f64) * 2.0 - 1.0
            };
            let mut x: Vec<f64> = (0..3 * n).map(|_| lcg() * 0.3).collect();
            let r = minimize(
                &obj,
                &mut x,
                &Options {
                    max_iter: 5000,
                    grad_tol: 1e-10,
                    memory: 8,
                },
            );
            println!(
                "距离盒 n={n} 半宽{half}: 迭代 {} 丢弃 {} 回溯 {} 值 {:.3e}",
                r.iterations, r.discarded, r.backtracks, r.value
            );
        }
    }
    #[test]
    fn 迭代次数必须配得上_lbfgs() {
        // **只断言"收敛到已知最优"是不够的。** 把 L-BFGS 退化成最速下降之后,
        // 这些目标照样收敛,只是慢一个数量级 —— 于是所有"收敛"判据全绿。
        // 实测:五条变异(曲率条件恒真、两循环丢掉 beta、初始缩放恒取 1、
        // 谨慎更新关掉、方向不下降不重置)**全部逃脱**。
        //
        // 所以这里断**迭代次数的上界**。上界是安全的:算法变好只会更少,
        // 不会让判据变假(与"断言随能力长进变成假话"那类不同)。
        //
        // # 闸设在哪儿,是量出来的
        //
        // | 变异 | Rosenbrock 8 维 | Powell | 距离盒 |
        // |---|---|---|---|
        // | 健康 | **69** | 152 | 24 |
        // | 两循环递推丢掉 beta | **110** | 145 | 30 |
        // | 初始 Hessian 缩放恒取 1 | **177** | 80 | 31 |
        // | 谨慎更新关掉 | 69 | 103 | 24 |
        // | 方向不下降时不重置 | 69 | 152 | 24 |
        //
        // 只有 Rosenbrock 那一列分得开,所以闸压到 **100**(健康值 69,余量 45%),
        // 它抓得住前两条。Powell 与距离盒那两条只当回归护栏(退化反而更快,
        // 分辨不了),闸留宽。
        //
        // **后两条变异抓不住,而且那是对的**:实测它们一步都不改变 ——
        // 强 Wolfe 已经保证 `y·s > 0`,所以谨慎更新与"方向不下降就重置"
        // 在这几个目标上**从来不触发**。它们是保险,不是主力;
        // 写判据去逼它们触发,等于断言一件不该发生的事。
        //
        // 起点用**交错的 −1.2 / 1.0**(Rosenbrock 的标准起点)。别随手写一组:
        // 链式 Rosenbrock 在 n ≥ 4 时有局部极小,实测
        // `[-1.2,1,0.7,-0.3,2,-1.5,0.2,3]` 会停在 f = 3.986 —— 优化器没错,
        // 那真是个局部极小,但拿它比迭代次数就没有意义了。
        let mut x: Vec<f64> = (0..8)
            .map(|i| if i % 2 == 0 { -1.2 } else { 1.0 })
            .collect();
        let r = minimize(
            &Rosenbrock,
            &mut x,
            &Options {
                max_iter: 5000,
                grad_tol: 1e-10,
                memory: 8,
            },
        );
        // 判**目标值**而不是 `converged`:8 维 Rosenbrock 在 f64 下压到 1e-10 的
        // 梯度范数已经到舍入极限,线搜索会先一步找不到合格点而正常退出 ——
        // 那不是"没收敛",是"到底了"。拿 `converged` 当判据等于在断言浮点精度。
        assert!(r.value < 1e-14, "8 维 Rosenbrock 目标值 {:.3e}", r.value);
        assert!(
            r.iterations < 100,
            "8 维 Rosenbrock 用了 {} 步(健康值 69)—— 拟牛顿那一套退化了",
            r.iterations
        );

        let mut x = vec![3.0, -1.0, 0.0, 1.0];
        let r = minimize(
            &Powell,
            &mut x,
            &Options {
                max_iter: 5000,
                grad_tol: 1e-12,
                memory: 8,
            },
        );
        assert!(
            r.iterations < 400,
            "Powell 用了 {} 步(健康值 152;这一条只当回归护栏)",
            r.iterations
        );

        let (obj, _) = distance_box(30, 0.1, 11);
        let mut st = 11u64 ^ 0xabcd;
        let mut lcg = || {
            st = st.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            ((st >> 11) as f64) / ((1u64 << 53) as f64) * 2.0 - 1.0
        };
        let mut x: Vec<f64> = (0..90).map(|_| lcg() * 0.3).collect();
        let r = minimize(
            &obj,
            &mut x,
            &Options {
                max_iter: 5000,
                grad_tol: 1e-10,
                memory: 8,
            },
        );
        assert!(r.value < 1e-12, "距离盒没压到零:{:.3e}", r.value);
        assert!(
            r.iterations < 60,
            "距离盒用了 {} 步(健康值 24;这一条只当回归护栏)",
            r.iterations
        );
    }

    #[test]
    fn 同样的输入两次给同样的答案() {
        // 求和顺序一旦被换成并行归约,这条就会红。
        let run = || {
            let mut x = vec![-1.2, 1.0, 0.7, -0.3, 2.0];
            let r = minimize(&Rosenbrock, &mut x, &Options::default());
            (x, r)
        };
        let (x1, r1) = run();
        let (x2, r2) = run();
        assert_eq!(x1, x2, "两次跑出来的坐标不逐位相同");
        assert_eq!(r1, r2, "两次跑出来的账不一样");
    }

    #[test]
    fn 起点就是最优时不乱动() {
        let mut x = vec![1.0, 1.0, 1.0];
        let r = minimize(&Rosenbrock, &mut x, &Options::default());
        assert!(r.converged);
        assert_eq!(r.iterations, 0, "起点已是最优,不该迭代");
        for v in &x {
            assert!((v - 1.0).abs() < 1e-15);
        }
    }
}
