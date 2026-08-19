//! **度量矩阵嵌入** —— 把一张距离表变成一组三维坐标。
//!
//! 这是经典多维标度(classical MDS / Crippen–Havel 嵌入):距离表 → 关于质心的
//! Gram 矩阵 → 特征分解 → 取最大的三个特征对当坐标。
//!
//! # 公式与它为什么长这样
//!
//! 手上只有两两距离 `d_ij`,想要坐标。先把坐标的原点定在质心上,那么每个原子到
//! 质心的平方距离可以只由距离表算出来:
//!
//! ```text
//! d_0i² = (1/n) Σ_j d_ij² − (1/n²) Σ_{j<k} d_jk²
//! ```
//!
//! 有了它,内积就出来了 —— 由余弦定理 `d_ij² = d_0i² + d_0j² − 2 x_i·x_j`:
//!
//! ```text
//! T_ij = x_i · x_j = ½ (d_0i² + d_0j² − d_ij²)
//! ```
//!
//! `T` 就是关于质心的 Gram 矩阵。它必然对称;**如果这张距离表真能在三维里摆出来**,
//! 那么 `T = X Xᵀ` 且秩恰好是 3,于是取三个最大特征对 `x_i[k] = √λ_k · v_k[i]`
//! 就精确还原了坐标(相差一个刚体变换)。
//!
//! # 参考距离表取哪一张:**取光滑化之后的上限矩阵 `U`**
//!
//! 这是本算法与 RDKit 唯一的实质分岔(见 [`crate::smooth`])。RDKit 在区间里
//! **逐对独立随机取**,取出来的表常常不是度量,`d_0i²` 于是可能为负 ——
//! 它一旦发现负值就**作废整次尝试**(`DistGeomUtils.cpp:114`),而病因是结构性的
//! 时候重掷 `10×N` 次会以同样的方式全部失败。
//!
//! `U` 是 Floyd–Warshall 的产物,**按构造满足三角不等式**,所以它是一张画得出来的
//! 距离表。而且它有个附带的好处:所有距离都取到上限,意味着初始结构是**最舒展**的
//! 那一个 —— 链是伸开的、没有自穿、没有打结。对后面的精修是个好起点。
//!
//! # 不合格时怎么办:**降级,不是失败**
//!
//! `U` 是度量不等于它能装进**三维**(正三角形的三个顶点是度量,却要二维;
//! `n` 个两两等距的点要 `n−1` 维)。所以前三个特征值仍可能有非正的。
//!
//! RDKit 遇到这个就 `return false` 让整次尝试作废。**这里不失败** ——
//! 非正的那一维坐标记为 0,同时把"离三维有多远"如实报回去
//! ([`Embedding::fit3`]、[`Embedding::negative_share`])。
//! 嵌入本来就只是给优化器一个起点,起点差一点是可以修的,直接没有起点才是灾难。
//!
//! # 判官
//!
//! 见 `examples/eigen_oracle.rs`。另有一条不依赖任何外部实现的硬判据,写在下面的
//! 单元测试里:**拿真实构象的精确距离表回嵌,还原出来的距离必须逐对相同** ——
//! 这一条把上面那两个公式完全钉死,写错任何一项都会当场红。

use crate::linalg::{symmetric_eigen, EigenError};
use crate::smooth::Bounds;

/// 三维。这个常数只是为了让下面的 `3` 都有名字,不是可调项 ——
/// 改成别的值整个模块的语义都不成立。
const DIM: usize = 3;

/// 特征值小到什么程度就当它是零。
///
/// 判的是**相对于最大特征值**的比例,不是绝对值:分子大小不同,谱的量级差很多,
/// 拿绝对阈值卡会在大分子上把真实的维度判掉。
///
/// (RDKit 这里用的是绝对阈值 `EIGVAL_TOL = 0.001`,而且同一个常数还兼任
/// `d_0i²` 的下限判据 —— 两个量纲不一样的东西共用一个数。)
const EIGVAL_REL_TOL: f64 = 1e-12;

/// 嵌入的结果。
#[derive(Debug, Clone, PartialEq)]
pub struct Embedding {
    /// 每个原子的坐标,**已经以质心为原点**(Gram 矩阵是关于质心建的,天然如此)。
    pub coords: Vec<[f64; DIM]>,
    /// 用上的三个特征值(降序)。非正的那些已经被压成 0。
    pub eigenvalues: [f64; DIM],
    /// 前三个特征值占**全谱正质量**的比例。1.0 表示这张距离表精确地是三维的。
    ///
    /// 这是"起点有多好"的度量,不是判据 —— 低了只说明后面的精修要多干点活。
    pub fit3: f64,
    /// 负特征值的绝对值占全谱绝对质量的比例。0.0 表示距离表是半正定的。
    ///
    /// 这个数直接量出"这张表离能摆出来有多远"。实测:RDKit 随机取距离是 0.268,
    /// 直接用光滑化后的 `U` 是 0.039。
    pub negative_share: f64,
    /// 前三个特征值里有几个非正 —— 有几个就有几根坐标轴被压成了 0。
    pub degenerate_axes: usize,
    /// 有几个原子的"到质心平方距离"算出来是负的。
    ///
    /// 这正是 RDKit 直接判死整次尝试的那个条件(`sqD0i[i] < EIGVAL_TOL`)。
    /// 这里只记账不判死,留给判官看。
    pub negative_centroid_sq: usize,
}

/// 嵌入失败的原因。**只有输入本身坏掉才会失败**,几何上摆不出来是降级不是失败。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EmbedError {
    /// 距离表里有非有限数,或者特征分解没收敛。
    Eigen(EigenError),
    /// 距离表不是 `n×n`。
    BadShape,
}

impl From<EigenError> for EmbedError {
    fn from(e: EigenError) -> Self {
        Self::Eigen(e)
    }
}

/// 从光滑化之后的界矩阵取参考距离表:**整张取上限**。
///
/// 理由见模块文档 —— `U` 按构造满足三角不等式,而下限一定 `≤` 上限,
/// 所以这张表同时也满足全部下限约束。
#[must_use]
pub fn reference_distances(b: &Bounds) -> Vec<f64> {
    let n = b.len();
    let mut d = vec![0.0; n * n];
    for i in 0..n {
        for j in (i + 1)..n {
            let u = b.upper(i, j);
            d[i * n + j] = u;
            d[j * n + i] = u;
        }
    }
    d
}

/// 由距离表算关于质心的 Gram 矩阵。
///
/// 返回 `(T, 到质心平方距离为负的原子数)`。
///
/// # Panics
///
/// `dist.len() != n * n` 时 panic。
#[must_use]
pub fn metric_matrix(dist: &[f64], n: usize) -> (Vec<f64>, usize) {
    assert_eq!(dist.len(), n * n, "距离表不是 {n}×{n}");
    if n == 0 {
        return (Vec::new(), 0);
    }
    let nf = n as f64;

    // Σ_{j<k} d_jk² / n² —— 注意求和只走**上三角**,不是整张表。
    // 走整张就等于把每一对数了两遍,`d_0i²` 会整体偏小 n 倍量级的一项。
    let mut sum_sq = 0.0;
    for i in 0..n {
        for j in (i + 1)..n {
            sum_sq += dist[i * n + j] * dist[i * n + j];
        }
    }
    sum_sq /= nf * nf;

    // d_0i²:第 i 个原子到质心的平方距离
    let mut sq0 = vec![0.0; n];
    let mut negative = 0;
    for i in 0..n {
        let row: f64 = (0..n).map(|j| dist[i * n + j] * dist[i * n + j]).sum();
        sq0[i] = row / nf - sum_sq;
        if sq0[i] < 0.0 {
            negative += 1;
        }
    }

    let mut t = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            t[i * n + j] = 0.5 * (sq0[i] + sq0[j] - dist[i * n + j] * dist[i * n + j]);
        }
    }
    (t, negative)
}

/// 把一张距离表嵌进三维。
///
/// # Errors
///
/// 距离表形状不对,或者含非有限数 / 特征分解不收敛。**几何上摆不出三维不算失败** ——
/// 那种情况把对应的坐标轴压成 0,并在 [`Embedding::fit3`] 里如实报出来。
pub fn embed(dist: &[f64], n: usize) -> Result<Embedding, EmbedError> {
    if dist.len() != n * n {
        return Err(EmbedError::BadShape);
    }
    let (t, negative_centroid_sq) = metric_matrix(dist, n);
    let eig = symmetric_eigen(&t, n)?;

    // 谱质量:正的一半用来算 fit3,绝对值那一半用来算负份额
    let mut pos_mass = 0.0;
    let mut abs_mass = 0.0;
    let mut neg_mass = 0.0;
    for &v in &eig.values {
        abs_mass += v.abs();
        if v > 0.0 {
            pos_mass += v;
        } else {
            neg_mass += -v;
        }
    }

    let scale = eig.values.first().map_or(0.0, |v| v.abs());
    let cut = EIGVAL_REL_TOL * scale;
    let mut eigenvalues = [0.0; DIM];
    let mut degenerate_axes = 0;
    for (k, lam) in eigenvalues.iter_mut().enumerate() {
        let v = eig.values.get(k).copied().unwrap_or(0.0);
        if v > cut {
            *lam = v;
        } else {
            // 这一根轴摆不出来:坐标压成 0,同时记账。
            // RDKit 在这里 `return false` 作废整次尝试;见模块文档。
            degenerate_axes += 1;
        }
    }

    let mut coords = vec![[0.0; DIM]; n];
    for (k, &lam) in eigenvalues.iter().enumerate() {
        // 退化轴的特征值上面已经压成了 0,`√0 · v = 0` 自然就把那一维写成 0,
        // 不需要再加一道 `if <= 0 { continue }` —— 那道判断永远与直接算等价。
        //
        // **写成 0 是有意的,不是凑合。** 真正平面的分子(苯环、酰胺)第三个
        // 特征值本来就该是 0,给它随手撒一点扰动等于把平面结构掰弯。
        // RDKit 那个 `randNegEig` 就是往这一维填随机数,而它填的时机是
        // "整次尝试本来要作废了" —— 与这里不是一回事。
        // 对称结构导致优化器卡住是**精修**阶段的事,在那里按需破对称。
        let s = lam.sqrt();
        // n < DIM 时特征向量根本不存在,跳过 —— 那一维本来就没有信息
        let Some(vk) = (k < eig.len()).then(|| eig.vector(k)) else {
            continue;
        };
        for (i, c) in coords.iter_mut().enumerate() {
            c[k] = s * vk[i];
        }
    }

    Ok(Embedding {
        coords,
        eigenvalues,
        fit3: if pos_mass > 0.0 {
            eigenvalues.iter().sum::<f64>() / pos_mass
        } else {
            0.0
        },
        negative_share: if abs_mass > 0.0 {
            neg_mass / abs_mass
        } else {
            0.0
        },
        degenerate_axes,
        negative_centroid_sq,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 确定性伪随机,只为造测试点集。
    fn lcg(state: &mut u64) -> f64 {
        *state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1);
        ((*state >> 11) as f64) / ((1u64 << 53) as f64) * 2.0 - 1.0
    }

    fn random_points(n: usize, seed: u64) -> Vec<[f64; 3]> {
        let mut st = seed;
        (0..n)
            .map(|_| [lcg(&mut st) * 6.0, lcg(&mut st) * 6.0, lcg(&mut st) * 6.0])
            .collect()
    }

    fn exact_distances(p: &[[f64; 3]]) -> Vec<f64> {
        let n = p.len();
        let mut d = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                d[i * n + j] = (0..3)
                    .map(|k| (p[i][k] - p[j][k]).powi(2))
                    .sum::<f64>()
                    .sqrt();
            }
        }
        d
    }

    /// 距离表之间的最大逐对偏差。距离在刚体变换下不变,所以这是比坐标更该比的量。
    fn max_dist_dev(a: &[f64], b: &[f64]) -> f64 {
        a.iter()
            .zip(b)
            .fold(0.0_f64, |w, (x, y)| w.max((x - y).abs()))
    }

    #[test]
    fn 真实三维点集精确回嵌() {
        // **这是整个模块最硬的一条判据。** 由真实三维坐标算出精确距离表,
        // 嵌回去必须逐对重合到机器精度。上面那两个公式写错任何一项都过不了。
        for n in [4usize, 5, 9, 17, 33] {
            for seed in [1u64, 0xf00d, 0xdead_beef] {
                let pts = random_points(n, seed);
                let d = exact_distances(&pts);
                let e = embed(&d, n).unwrap();
                let back = exact_distances(&e.coords);
                let dev = max_dist_dev(&d, &back);
                assert!(dev < 1e-9, "n={n} seed={seed} 距离偏差 {dev:.3e}");
                assert_eq!(e.degenerate_axes, 0, "n={n} seed={seed} 不该有退化轴");
                assert!(e.fit3 > 1.0 - 1e-12, "n={n} seed={seed} fit3={}", e.fit3);
                assert!(
                    e.negative_share < 1e-12,
                    "n={n} 负份额 {}",
                    e.negative_share
                );
                assert_eq!(e.negative_centroid_sq, 0);
            }
        }
    }

    #[test]
    fn 嵌出来的坐标以质心为原点() {
        let pts = random_points(11, 7);
        let e = embed(&exact_distances(&pts), 11).unwrap();
        for k in 0..3 {
            let c: f64 = e.coords.iter().map(|p| p[k]).sum::<f64>() / 11.0;
            assert!(c.abs() < 1e-12, "第 {k} 轴质心 {c:.3e} 不在原点");
        }
    }

    #[test]
    fn 平面点集只用掉两根轴() {
        // 全在 z=0 上的点:第三个特征值必须是 0,而且距离仍要精确还原。
        let mut st = 3u64;
        let pts: Vec<[f64; 3]> = (0..8)
            .map(|_| [lcg(&mut st) * 4.0, lcg(&mut st) * 4.0, 0.0])
            .collect();
        let d = exact_distances(&pts);
        let e = embed(&d, 8).unwrap();
        assert_eq!(e.degenerate_axes, 1, "特征值 {:?}", e.eigenvalues);
        assert!(max_dist_dev(&d, &exact_distances(&e.coords)) < 1e-9);
    }

    #[test]
    fn 接近平面但不是平面的结构不许被压平() {
        // 芳环轻微皱折、酰胺略微扭出平面 —— 第三维是**真的有**,只是很小。
        // 零特征值的阈值定得太松就会把这一维判掉,而结构被压平之后距离全错。
        //
        // 这里 z 的幅度是 x/y 的 1/200,于是 λ₃/λ₁ ≈ 2.5e-5:
        // 它落在"真零"(1e-12)与"松阈值"(比如 1e-3)之间,正好把阈值钉住。
        let mut st = 5u64;
        let pts: Vec<[f64; 3]> = (0..10)
            .map(|_| [lcg(&mut st) * 4.0, lcg(&mut st) * 4.0, lcg(&mut st) * 0.02])
            .collect();
        let d = exact_distances(&pts);
        let e = embed(&d, 10).unwrap();
        let ratio = e.eigenvalues[2] / e.eigenvalues[0];
        assert!(
            (1e-12..1e-3).contains(&ratio),
            "构造失效:λ₃/λ₁ = {ratio:.3e} 没落在要钉的区间里"
        );
        assert_eq!(e.degenerate_axes, 0, "第三维是真的,不许判成退化");
        let dev = max_dist_dev(&d, &exact_distances(&e.coords));
        assert!(dev < 1e-9, "距离偏差 {dev:.3e} —— 结构被压平了");
    }

    #[test]
    fn 摆不进三维时降级而不是失败() {
        // 五个两两等距的点要四维才摆得下。RDKit 在这里判死整次尝试;
        // 这里必须**照样给出坐标**,并且把差距如实报出来。
        let n = 5;
        let mut d = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                if i != j {
                    d[i * n + j] = 1.0;
                }
            }
        }
        let e = embed(&d, n).expect("必须给出坐标,不能失败");
        assert_eq!(e.coords.len(), n);
        assert!(e.fit3 < 1.0, "fit3={} 应当小于 1", e.fit3);
        // 正单形的 Gram 矩阵是半正定的,所以负份额是 0 —— 它只是维数不够
        assert!(e.negative_share < 1e-12, "负份额 {}", e.negative_share);
        assert_eq!(e.degenerate_axes, 0, "四维单形的前三个特征值都是正的");
    }

    #[test]
    fn 自相矛盾的距离表会露出负特征值() {
        // "A 离 B 一米、B 离 C 一米、A 离 C 十米" —— 违反三角不等式,
        // 任何维度都摆不出来。必须体现为负特征值,而且**仍然要返回坐标**。
        let n = 3;
        let mut d = vec![0.0; n * n];
        let set = |d: &mut Vec<f64>, i: usize, j: usize, v: f64| {
            d[i * n + j] = v;
            d[j * n + i] = v;
        };
        set(&mut d, 0, 1, 1.0);
        set(&mut d, 1, 2, 1.0);
        set(&mut d, 0, 2, 10.0);
        let e = embed(&d, n).expect("坏表也要给坐标");
        assert!(e.negative_share > 0.01, "负份额只有 {}", e.negative_share);
        // **坐标必须永远是有限数。** 负特征值一旦被当成有效轴,`√负数` 就是 NaN,
        // 而 NaN 会一路淌进后面的精修,那时再查就晚了。零阈值卡在正数一侧
        // (`v > cut`,`cut > 0`)保证了这一点,这条断言是它的闸。
        for (i, c) in e.coords.iter().enumerate() {
            assert!(
                c.iter().all(|x| x.is_finite()),
                "第 {i} 个原子的坐标不是有限数:{c:?}"
            );
        }
    }

    #[test]
    fn 参考距离表整张取上限() {
        let mut b = Bounds::new(3, 1.0, 5.0);
        b.set_upper(0, 1, 2.5);
        b.set_lower(0, 1, 2.0);
        let d = reference_distances(&b);
        let at = |i: usize, j: usize| d[i * 3 + j];
        assert_eq!(at(0, 1), 2.5);
        assert_eq!(at(1, 0), 2.5, "必须对称");
        assert_eq!(at(0, 0), 0.0, "对角必须是 0");
        assert_eq!(at(0, 2), 5.0);
    }

    #[test]
    fn 空表与退化尺寸() {
        let e = embed(&[], 0).unwrap();
        assert!(e.coords.is_empty());
        // 一个原子:坐标在原点,三根轴全退化
        let e = embed(&[0.0], 1).unwrap();
        assert_eq!(e.coords, vec![[0.0; 3]]);
        assert_eq!(e.degenerate_axes, 3);
        // 两个原子:只有一根轴有信息
        let d = vec![0.0, 1.5, 1.5, 0.0];
        let e = embed(&d, 2).unwrap();
        assert_eq!(e.degenerate_axes, 2);
        let back = exact_distances(&e.coords);
        assert!((back[1] - 1.5).abs() < 1e-12, "两原子距离 {}", back[1]);
    }

    #[test]
    fn 形状不对要报错() {
        assert_eq!(embed(&[1.0, 2.0], 3), Err(EmbedError::BadShape));
    }

    #[test]
    fn 平移旋转不改变结果的距离() {
        // 嵌入只认距离,所以把输入点整体搬走 / 转过去,嵌出来的距离表必须不变。
        let pts = random_points(9, 99);
        let d1 = exact_distances(&pts);
        let (c, s) = (0.6_f64, 0.8_f64);
        let moved: Vec<[f64; 3]> = pts
            .iter()
            .map(|p| {
                [
                    c * p[0] - s * p[1] + 100.0,
                    s * p[0] + c * p[1] - 50.0,
                    p[2] + 7.0,
                ]
            })
            .collect();
        let d2 = exact_distances(&moved);
        assert!(max_dist_dev(&d1, &d2) < 1e-12, "构造有问题:距离本该不变");
        let e1 = embed(&d1, 9).unwrap();
        let e2 = embed(&d2, 9).unwrap();
        let dev = max_dist_dev(&exact_distances(&e1.coords), &exact_distances(&e2.coords));
        assert!(dev < 1e-9, "偏差 {dev:.3e}");
    }
}
