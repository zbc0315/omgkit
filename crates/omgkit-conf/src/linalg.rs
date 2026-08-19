//! **对称矩阵的特征分解** —— 手写循环 Jacobi 旋转。
//!
//! 嵌入那一步要把度量矩阵(Gram 矩阵)分解成 `V^T Λ V`,取最大的三个特征对当坐标。
//! 这个模块只做线性代数,不认识分子。
//!
//! # 为什么不用 RDKit 那套幂迭代
//!
//! RDKit 走的是 `PowerEigenSolver`(带收缩的幂法),取前 `dim` 个特征对。
//! 它有四处结构性的毛病:
//!
//! 1. **幂法收敛到的是绝对值最大的特征值,不是最大的。** Gram 矩阵本该半正定,
//!    但参考距离表摆不进任何维度时它就不是 —— 这时一个**大负特征值**会被当成
//!    "第一特征对"取回来,而真正的最大正特征值排在后面。
//! 2. **要一个随机起始向量**(`v.setToRandom(seed)`)。RDKit 用
//!    `(int)(sumSqD2 * N)` 当种子,所以它其实是确定的,但这个种子是距离表的
//!    一个**有损函数** —— 两张不同的表可能撞出同一个起手向量,也可能因为一个
//!    末位比特的差别换一个。
//! 3. **不收敛时它静默返回垃圾。** `powerEigenSolver` 返回 `bool`,
//!    某一个特征对迭代满 1000 次没收敛就 `break`,后面的特征值原地留着初值;
//!    而 `computeInitialCoords` **根本没接这个返回值**
//!    (`DistGeomUtils.cpp:126`)。特征值相近时幂法收敛得任意慢,而分子的
//!    Gram 矩阵在近似对称的结构上正好有成对的重特征值。
//! 4. **收缩会累积误差**:每取一个特征对就从矩阵里减掉 `λ v v^T`,
//!    第三个特征对身上带着前两个的误差。而它的收敛判据是 `1e-3` —— 很松。
//!
//! Jacobi 没有这四条:它**一次算出全部** `n` 个特征对,不需要起始向量,
//! 扫描顺序固定所以完全确定,精度到机器 eps,而且对小特征值有**高相对精度**
//! (这一条正是判断"这张距离表能不能装进三维"所需要的 —— 要看的就是那些
//! 接近零的特征值究竟是正是负)。
//!
//! 代价是 `O(n³)` 每一轮、要扫若干轮,比只取三个特征对的幂法贵。
//! **但本算法一个分子只跑一次嵌入**(参考距离表是确定的,没有重掷),
//! 而 RDKit 每个分子最多跑 `10×N` 次 —— 单次贵一点换掉几十次重试是划算的。
//!
//! # 判官
//!
//! 见 `examples/eigen_oracle.rs`:拿真实分子的 Gram 矩阵与 numpy 的
//! `linalg.eigh`(LAPACK)逐个特征值比。另外还有一条不依赖任何外部实现的硬判据 ——
//! **真实三维构象的 Gram 矩阵秩必须恰好是 3**,多一个非零特征值就是错的。

/// 最多扫多少轮。Jacobi 典型 6~10 轮收敛;扫到这个数还不收敛说明输入病态,
/// 这时候要报错,不能返回一个"扫累了"的近似解。
const MAX_SWEEPS: usize = 60;

/// 对称矩阵特征分解的结果。
#[derive(Debug, Clone, PartialEq)]
pub struct Eigen {
    /// 特征值,**降序**。长度 `n`。
    pub values: Vec<f64>,
    /// 特征向量,行主序 `n×n`,**第 `k` 行是第 `k` 个特征向量**,已归一化。
    ///
    /// 按行存(而不是照 Jacobi 内部那样按列)是因为取用的一方总是
    /// "把第 `k` 个特征向量整根拿走",按行存这一步是连续内存。
    pub vectors: Vec<f64>,
    /// 实际扫了几轮。诊断用 —— 轮数异常高就是矩阵病态的信号。
    pub sweeps: usize,
}

impl Eigen {
    /// 矩阵阶数。
    #[must_use]
    pub fn len(&self) -> usize {
        self.values.len()
    }

    /// 是不是空矩阵。
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    /// 第 `k` 个特征向量(降序,`k = 0` 是最大特征值那个)。
    #[must_use]
    pub fn vector(&self, k: usize) -> &[f64] {
        let n = self.len();
        &self.vectors[k * n..(k + 1) * n]
    }
}

/// 特征分解失败的原因。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EigenError {
    /// 输入里有非有限数(NaN / 无穷)。
    ///
    /// **必须挡在门口**:收敛判据是"上三角非对角元的绝对值之和等于 0",
    /// 而 `NaN == 0.0` 恒假 —— 放进来就是扫满 `MAX_SWEEPS` 轮、
    /// 返回一堆 NaN,而调用方拿到的是 `Ok`。
    NotFinite {
        /// 第一个出问题的元素在行主序里的下标。
        index: usize,
    },
    /// 扫满 `MAX_SWEEPS` 轮仍未收敛。
    NoConverge {
        /// 扫了几轮。
        sweeps: usize,
    },
}

/// 判断 `g` 相对于 `x` 是否已经小到可以忽略。
///
/// 这是 Jacobi 的经典写法 `|x| + g == |x|` 的显式版本 —— 那种浮点等号写法
/// 意图要靠猜,这里直接写成相对比较。
#[inline]
fn is_negligible(g: f64, x: f64) -> bool {
    g <= f64::EPSILON * x.abs()
}

/// Jacobi 的成对旋转:把 `m[a]` 与 `m[b]` 一起转过去(`a`/`b` 是行主序的平铺下标)。
#[inline]
fn rot(m: &mut [f64], a: usize, b: usize, s: f64, tau: f64) {
    let g = m[a];
    let h = m[b];
    m[a] = g - s * (h + g * tau);
    m[b] = h + s * (g - h * tau);
}

/// 求实对称矩阵 `a`(行主序 `n×n`)的全部特征对。
///
/// **只读上三角与对角**,下三角当它不存在 —— 但 debug 构建下会断言矩阵确实对称,
/// 免得调用方只填了一半还以为读的是自己填的那一半。
///
/// # Errors
///
/// 输入含非有限数,或扫满 `MAX_SWEEPS` 轮未收敛。
///
/// # Panics
///
/// `a.len() != n * n` 时 panic —— 这是调用方的编码错误,不是运行时状况。
pub fn symmetric_eigen(a: &[f64], n: usize) -> Result<Eigen, EigenError> {
    assert_eq!(a.len(), n * n, "矩阵不是 {n}×{n}");
    for (i, &x) in a.iter().enumerate() {
        if !x.is_finite() {
            return Err(EigenError::NotFinite { index: i });
        }
    }
    debug_assert!(
        (0..n).all(|i| (0..i).all(|j| {
            let (x, y) = (a[i * n + j], a[j * n + i]);
            (x - y).abs() <= 1e-12 * (1.0 + x.abs().max(y.abs()))
        })),
        "矩阵不对称:只读上三角,下三角与它对不上说明调用方填错了一半"
    );

    let mut m = a.to_vec();
    // `v` 的第 k **列**是第 k 个特征向量(Jacobi 内部的自然布局),最后转置出去
    let mut v = vec![0.0; n * n];
    for i in 0..n {
        v[i * n + i] = 1.0;
    }
    // `d` 是当前对角;`b`/`z` 是 Jacobi 的累加技巧:一轮里的增量先攒在 `z` 上,
    // 轮末一次性加到 `b` —— 直接往 `d` 上累加会把一堆小量加到大量上,精度掉得快。
    let mut d: Vec<f64> = (0..n).map(|i| m[i * n + i]).collect();
    let mut b = d.clone();
    let mut z = vec![0.0; n];

    for sweep in 1..=MAX_SWEEPS {
        // 上三角非对角元绝对值之和。Jacobi 的每一次旋转都让它严格下降,
        // 降到 0 就是收敛(不是"足够小"就停 —— 是真的到 0)。
        let mut sm = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                sm += m[i * n + j].abs();
            }
        }
        if sm == 0.0 {
            return Ok(finish(n, &d, &v, sweep - 1));
        }
        // 前三轮设一道门槛,只转"大"的元素。小元素这时候转纯属做无用功 ——
        // 后面的大旋转会把它们再搅乱一遍。
        let tresh = if sweep < 4 {
            0.2 * sm / (n * n) as f64
        } else {
            0.0
        };

        for p in 0..n {
            for q in (p + 1)..n {
                let apq = m[p * n + q];
                let g = 100.0 * apq.abs();
                // 这一项相对两个对角元已经小到加上去都不改变它们:转它只会引入
                // 舍入噪声,直接置零。前四轮不做这个判断 —— 那时对角还没稳定。
                if sweep > 4 && is_negligible(g, d[p]) && is_negligible(g, d[q]) {
                    m[p * n + q] = 0.0;
                    continue;
                }
                if apq.abs() <= tresh {
                    continue;
                }
                // 解 2×2 子块的旋转角。`t = tan θ`,取绝对值较小的那个根 ——
                // 它对应 |θ| < π/4,是数值上稳的那一支。
                let h = d[q] - d[p];
                let t = if is_negligible(g, h) {
                    // 这一支只在 |h| 远大于 g 时走到,而 g > 0(否则上面就 continue 了),
                    // 所以 h 必不为零,除法安全。
                    apq / h
                } else {
                    let theta = 0.5 * h / apq;
                    let t = 1.0 / (theta.abs() + (1.0 + theta * theta).sqrt());
                    if theta < 0.0 {
                        -t
                    } else {
                        t
                    }
                };
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = t * c;
                let tau = s / (1.0 + c);
                let dh = t * apq;
                z[p] -= dh;
                z[q] += dh;
                d[p] -= dh;
                d[q] += dh;
                m[p * n + q] = 0.0;

                // 转 m 的第 p、q 行列。三段拆开是为了始终只碰上三角。
                for j in 0..p {
                    rot(&mut m, j * n + p, j * n + q, s, tau);
                }
                for j in (p + 1)..q {
                    rot(&mut m, p * n + j, j * n + q, s, tau);
                }
                for j in (q + 1)..n {
                    rot(&mut m, p * n + j, q * n + j, s, tau);
                }
                // 转 v 的第 p、q 列
                for j in 0..n {
                    rot(&mut v, j * n + p, j * n + q, s, tau);
                }
            }
        }
        for i in 0..n {
            b[i] += z[i];
            d[i] = b[i];
            z[i] = 0.0;
        }
    }
    Err(EigenError::NoConverge { sweeps: MAX_SWEEPS })
}

/// 收尾:按特征值降序排,转置成按行存,再把符号钉死。
fn finish(n: usize, d: &[f64], v: &[f64], sweeps: usize) -> Eigen {
    let mut order: Vec<usize> = (0..n).collect();
    // 降序。特征值相等时 `sort_by` 是稳定排序,保持 Jacobi 的原顺序,
    // 于是同样的输入永远给同样的输出。
    order.sort_by(|&x, &y| {
        d[y].partial_cmp(&d[x])
            .expect("非有限数已在入口挡掉,这里不可能有 NaN")
    });
    let mut values = Vec::with_capacity(n);
    let mut vectors = vec![0.0; n * n];
    for (k, &src) in order.iter().enumerate() {
        values.push(d[src]);
        for i in 0..n {
            vectors[k * n + i] = v[i * n + src];
        }
        canonical_sign(&mut vectors[k * n..(k + 1) * n]);
    }
    Eigen {
        values,
        vectors,
        sweeps,
    }
}

/// 把特征向量的符号钉死:**绝对值最大的分量取正**,并列时取下标最小的那个。
///
/// 不钉的话同一个分子可能拿到**镜像**的坐标:坐标是 `sqrt(λ_j)·v_j[i]`,
/// 把 `v_j` 整根翻号就是把第 `j` 根坐标轴翻过去 —— 那是一次反射,手性跟着翻。
///
/// **注意这条约定钉不住重特征值。** 特征值相等时,那个特征子空间里的基
/// 只定到一个旋转,符号约定管不着。旋转不改变分子(整体转一下而已),
/// 所以无害;但也因此,**整体手性不能指望嵌入这一步给对** —— 那是后面
/// 精修阶段带手性体积项去管的事。
fn canonical_sign(v: &mut [f64]) {
    let mut best = 0usize;
    for i in 1..v.len() {
        if v[i].abs() > v[best].abs() {
            best = i;
        }
    }
    if v.get(best).is_some_and(|x| *x < 0.0) {
        for x in v.iter_mut() {
            *x = -*x;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 重建残差 `‖A − Σ λ_k v_k v_k^T‖_max`。
    ///
    /// **这一条不是自洽性检查**:一个什么都不做的实现(返回对角线当特征值、
    /// 单位阵当特征向量)在非对角矩阵上会立刻被它抓住。
    fn recon_err(a: &[f64], n: usize, e: &Eigen) -> f64 {
        let mut worst: f64 = 0.0;
        for i in 0..n {
            for j in 0..n {
                let mut s = 0.0;
                for k in 0..n {
                    s += e.values[k] * e.vectors[k * n + i] * e.vectors[k * n + j];
                }
                worst = worst.max((s - a[i * n + j]).abs());
            }
        }
        worst
    }

    /// 正交性偏差 `‖V V^T − I‖_max`。
    fn orth_err(n: usize, e: &Eigen) -> f64 {
        let mut worst: f64 = 0.0;
        for k in 0..n {
            for l in 0..n {
                let mut s = 0.0;
                for i in 0..n {
                    s += e.vectors[k * n + i] * e.vectors[l * n + i];
                }
                worst = worst.max((s - f64::from(u8::from(k == l))).abs());
            }
        }
        worst
    }

    /// 确定性的伪随机,只为造测试矩阵 —— 不进生产路径。
    fn lcg(state: &mut u64) -> f64 {
        *state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1);
        ((*state >> 11) as f64) / ((1u64 << 53) as f64) * 2.0 - 1.0
    }

    fn sym_random(n: usize, seed: u64) -> Vec<f64> {
        let mut st = seed;
        let mut a = vec![0.0; n * n];
        for i in 0..n {
            for j in i..n {
                let x = lcg(&mut st);
                a[i * n + j] = x;
                a[j * n + i] = x;
            }
        }
        a
    }

    #[test]
    fn 二阶的解析解() {
        // [[2,1],[1,2]] 的特征值是 3 与 1,特征向量 (1,1)/√2 与 (1,−1)/√2
        let a = [2.0, 1.0, 1.0, 2.0];
        let e = symmetric_eigen(&a, 2).unwrap();
        assert!((e.values[0] - 3.0).abs() < 1e-14, "{:?}", e.values);
        assert!((e.values[1] - 1.0).abs() < 1e-14, "{:?}", e.values);
        let r = std::f64::consts::FRAC_1_SQRT_2;
        assert!((e.vector(0)[0] - r).abs() < 1e-14);
        assert!((e.vector(0)[1] - r).abs() < 1e-14);
        // 符号约定:两个分量绝对值相等 → 取下标最小的那个为正
        assert!((e.vector(1)[0] - r).abs() < 1e-14, "{:?}", e.vector(1));
        assert!((e.vector(1)[1] + r).abs() < 1e-14, "{:?}", e.vector(1));
    }

    #[test]
    fn 已经是对角阵时不动它() {
        let a = [3.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 7.0];
        let e = symmetric_eigen(&a, 3).unwrap();
        assert_eq!(e.values, vec![7.0, 3.0, 1.0]);
        assert_eq!(e.sweeps, 0, "对角阵一轮都不该扫");
        assert!(recon_err(&a, 3, &e) < 1e-15);
    }

    #[test]
    fn 退化特征值() {
        // 单位阵:三重特征值 1。特征向量任取一组正交基都对,只要求重建与正交。
        let mut a = vec![0.0; 9];
        for i in 0..3 {
            a[i * 3 + i] = 1.0;
        }
        let e = symmetric_eigen(&a, 3).unwrap();
        for v in &e.values {
            assert!((v - 1.0).abs() < 1e-15);
        }
        assert!(orth_err(3, &e) < 1e-14);
        assert!(recon_err(&a, 3, &e) < 1e-14);
    }

    #[test]
    fn 负特征值也要拿对() {
        // [[0,1],[1,0]] → +1 与 −1。幂法在这里会先撞上哪个取决于起手向量;
        // Jacobi 必须两个都给对,而且**降序**。
        let a = [0.0, 1.0, 1.0, 0.0];
        let e = symmetric_eigen(&a, 2).unwrap();
        assert!((e.values[0] - 1.0).abs() < 1e-14, "{:?}", e.values);
        assert!((e.values[1] + 1.0).abs() < 1e-14, "{:?}", e.values);
    }

    #[test]
    fn 随机对称阵的重建与正交() {
        for n in [1usize, 2, 3, 5, 8, 13, 21, 40] {
            for seed in [1u64, 0xf00d, 0xdead_beef] {
                let a = sym_random(n, seed);
                let e = symmetric_eigen(&a, n).unwrap();
                assert!(
                    recon_err(&a, n, &e) < 1e-12,
                    "n={n} seed={seed} 重建残差 {}",
                    recon_err(&a, n, &e)
                );
                assert!(orth_err(n, &e) < 1e-12, "n={n} seed={seed} 正交性");
                // 降序
                for k in 1..n {
                    assert!(e.values[k - 1] >= e.values[k], "n={n} 没按降序排");
                }
            }
        }
    }

    #[test]
    fn 空矩阵与一阶() {
        let e = symmetric_eigen(&[], 0).unwrap();
        assert!(e.is_empty());
        let e = symmetric_eigen(&[-2.5], 1).unwrap();
        assert_eq!(e.values, vec![-2.5]);
        assert_eq!(e.vectors, vec![1.0]);
    }

    #[test]
    fn 三维点集的_gram_矩阵秩恰好是三() {
        // 这一条是嵌入那一步真正依赖的性质,也是最能抓错的一条:
        // 由**真实三维坐标**造出来的 Gram 矩阵,非零特征值必须恰好三个。
        // 少一个说明丢了信息,多一个说明算错了 —— 这不是自洽性,是几何事实。
        let mut st = 42u64;
        let n = 12;
        let pts: Vec<[f64; 3]> = (0..n)
            .map(|_| [lcg(&mut st) * 5.0, lcg(&mut st) * 5.0, lcg(&mut st) * 5.0])
            .collect();
        // 先移到质心,Gram 矩阵才是 x_i·x_j
        let mut c = [0.0; 3];
        for p in &pts {
            for k in 0..3 {
                c[k] += p[k] / n as f64;
            }
        }
        let mut g = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                g[i * n + j] = (0..3)
                    .map(|k| (pts[i][k] - c[k]) * (pts[j][k] - c[k]))
                    .sum();
            }
        }
        let e = symmetric_eigen(&g, n).unwrap();
        let scale = e.values[0].abs();
        assert!(e.values[2] > 1e-6 * scale, "前三个特征值必须都显著为正");
        for k in 3..n {
            assert!(
                e.values[k].abs() < 1e-10 * scale,
                "第 {k} 个特征值应当是零,实际 {}",
                e.values[k]
            );
        }
    }

    #[test]
    fn 非有限数挡在门口() {
        let a = [1.0, 0.0, 0.0, f64::NAN];
        assert_eq!(
            symmetric_eigen(&a, 2),
            Err(EigenError::NotFinite { index: 3 })
        );
        let a = [f64::INFINITY, 0.0, 0.0, 1.0];
        assert_eq!(
            symmetric_eigen(&a, 2),
            Err(EigenError::NotFinite { index: 0 })
        );
    }

    #[test]
    fn 符号约定是确定的() {
        // 同一个矩阵跑两次必须逐位相同,而且约定确实生效:
        // 每根特征向量绝对值最大的分量为正。
        let a = sym_random(7, 0x1234);
        let e1 = symmetric_eigen(&a, 7).unwrap();
        let e2 = symmetric_eigen(&a, 7).unwrap();
        assert_eq!(e1, e2, "同样的输入两次给的结果不一样");
        for k in 0..7 {
            let v = e1.vector(k);
            let best = (0..7)
                .max_by(|&x, &y| v[x].abs().total_cmp(&v[y].abs()))
                .unwrap();
            assert!(v[best] > 0.0, "第 {k} 根特征向量的最大分量是负的");
        }
    }

    #[test]
    fn 病态的近重特征值() {
        // 相邻特征值只差 1e-10 —— 幂法在这里要迭代到天荒地老(收敛率是
        // λ₂/λ₁ 的幂),Jacobi 照样几轮出结果且两个值都对。
        let a = [1.0, 0.0, 0.0, 1.0 + 1e-10];
        let e = symmetric_eigen(&a, 2).unwrap();
        assert!(
            (e.values[0] - (1.0 + 1e-10)).abs() < 1e-15,
            "{:?}",
            e.values
        );
        assert!((e.values[1] - 1.0).abs() < 1e-15, "{:?}", e.values);
        assert!(e.sweeps <= 4, "扫了 {} 轮,太多", e.sweeps);
    }
}
