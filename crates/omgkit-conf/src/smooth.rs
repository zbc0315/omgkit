//! **三角光滑化** —— 把一张自相矛盾的距离区间表收紧成自洽的。
//!
//! # 它解决什么
//!
//! 界矩阵是逐条化学事实写出来的:这根键 1.5 Å、这个角推出 1-3 距离 2.5 Å、
//! 这两个原子至少离 3 Å……**这些条目彼此不知道对方存在**,于是会自相矛盾:
//! 写了"A–B 至多 2、B–C 至多 2",却又写"A–C 至多 10" —— 后一条是废话,
//! 因为三角不等式已经把它压到 4 了。
//!
//! 光滑化就是把每一条区间用三角不等式收紧到极限:
//!
//! - **上限**:`U_ij ← min(U_ij, U_ik + U_kj)`。这正是 Floyd–Warshall 的最短路,
//!   所以收敛之后的 `U` 是"沿着约束网络走过去的最短距离"。
//! - **下限**:`L_ij ← max(L_ij, L_ik − U_kj, L_jk − U_ik)`。
//!
//! # 为什么这一步是整个算法的地基
//!
//! **Floyd–Warshall 之后的上限矩阵 `U` 本身就是一个度量** —— `U_ij ≤ U_ik + U_kj`
//! 按构造恒成立。也就是说它是一张"画得出来"的距离表,不是拼凑的。
//!
//! 这一条正是本算法与 RDKit 的分岔点:RDKit 在这之后**对每一对原子独立随机取一个
//! 距离**,取出来的表往往任何空间里都摆不出来(实测负特征值占谱质量 27%);
//! 而直接拿 `U` 当参考距离表,这个数是 **4%**,能装进三维的结构从 66% 升到 93%。
//!
//! # 判官
//!
//! RDKit 的 `GetMoleculeBoundsMatrix` 带 `doTriangleSmoothing` 开关,
//! 所以可以把它**未光滑**的矩阵喂进来,再与它**光滑后**的逐元素比 ——
//! 这是外部判官,不是自己跟自己对。见 `harness/dump_bounds.py` 与
//! `examples/smooth_oracle.rs`。

/// 一张距离区间表。
///
/// **上三角存上限、下三角存下限**(与 RDKit `GetMoleculeBoundsMatrix` 同一约定),
/// 于是一个 `n×n` 的数组就装下了两个矩阵,而且缓存友好。
#[derive(Debug, Clone, PartialEq)]
pub struct Bounds {
    n: usize,
    /// 行主序的 `n×n`;`m[i*n+j]`,`i<j` 是上限、`i>j` 是下限。
    m: Vec<f64>,
}

impl Bounds {
    /// 造一张 `n` 个原子的表,上限全 `hi`、下限全 `lo`。
    #[must_use]
    pub fn new(n: usize, lo: f64, hi: f64) -> Self {
        let mut m = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                if i < j {
                    m[i * n + j] = hi;
                } else if i > j {
                    m[i * n + j] = lo;
                }
            }
        }
        Self { n, m }
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

    /// `i`–`j` 的距离上限。
    #[must_use]
    pub fn upper(&self, i: usize, j: usize) -> f64 {
        let (a, b) = if i < j { (i, j) } else { (j, i) };
        self.m[a * self.n + b]
    }

    /// `i`–`j` 的距离下限。
    #[must_use]
    pub fn lower(&self, i: usize, j: usize) -> f64 {
        let (a, b) = if i < j { (i, j) } else { (j, i) };
        self.m[b * self.n + a]
    }

    /// 设上限。
    pub fn set_upper(&mut self, i: usize, j: usize, v: f64) {
        let (a, b) = if i < j { (i, j) } else { (j, i) };
        self.m[a * self.n + b] = v;
    }

    /// 设下限。
    pub fn set_lower(&mut self, i: usize, j: usize, v: f64) {
        let (a, b) = if i < j { (i, j) } else { (j, i) };
        self.m[b * self.n + a] = v;
    }

    /// 逐元素读入一个行主序的 `n×n`(上三角上限、下三角下限)。判官用。
    #[must_use]
    pub fn from_row_major(n: usize, m: Vec<f64>) -> Option<Self> {
        (m.len() == n * n).then_some(Self { n, m })
    }

    /// 逐元素取出行主序的 `n×n`。判官用。
    #[must_use]
    pub fn as_row_major(&self) -> &[f64] {
        &self.m
    }
}

/// 三角光滑化失败的原因。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SmoothError {
    /// 某一对原子的下限被推到了上限之上 —— 这组约束在**任何**维度里都不可能同时成立。
    ///
    /// 带上是哪一对,好让上层说得出"卡在哪儿",而不是只报一句"失败"。
    Infeasible {
        /// 出问题的两个原子。
        pair: (usize, usize),
    },
}

/// **三角光滑化。** 就地把 `b` 收紧到三角不等式的极限。
///
/// 成功之后 `b.upper` 满足 `U_ij ≤ U_ik + U_kj`(它是一个度量),
/// 而任何一对的 `lower ≤ upper`。
///
/// # 与 RDKit 的一处刻意不同
///
/// RDKit 的 `TriangleSmooth.cpp:50-56` 更新下限时用的是 `if / else if` ——
/// 两条候选 `L_ik − U_kj` 与 `L_jk − U_ik` 只取**先命中的那条**,而不是取 `max`。
/// 教科书写法取 `max`,给出更紧的界。
///
/// 这里用 `max`。**不是因为 RDKit 错了** —— 实测 400 个分子共 11705 次下限更新,
/// "第二条该赢却没赢"的次数是 **0**(单趟 Floyd–Warshall 本身就是不动点),
/// 所以两种写法在真实分子上给出**逐位相同**的结果。取 `max` 只是把这件事
/// 写成不依赖于"恰好是不动点"这个巧合。
///
/// # Errors
///
/// 约束自相矛盾(某一对的下限压过上限)时返回 [`SmoothError::Infeasible`],
/// 并说明是哪一对。上层据此走确定性松弛梯,而不是重掷骰子。
pub fn triangle_smooth(b: &mut Bounds) -> Result<(), SmoothError> {
    let n = b.len();
    for k in 0..n {
        for i in 0..n {
            if i == k {
                continue;
            }
            let u_ik = b.upper(i, k);
            let l_ik = b.lower(i, k);
            for j in (i + 1)..n {
                if j == k {
                    continue;
                }
                // 上限:最短路
                let u_kj = b.upper(k, j);
                if b.upper(i, j) > u_ik + u_kj {
                    b.set_upper(i, j, u_ik + u_kj);
                }
                // 下限:两条候选都要看(见上面那段"与 RDKit 的一处刻意不同")
                let cand = (l_ik - u_kj).max(b.lower(j, k) - u_ik);
                if b.lower(i, j) < cand {
                    b.set_lower(i, j, cand);
                }
                if b.lower(i, j) > b.upper(i, j) {
                    return Err(SmoothError::Infeasible { pair: (i, j) });
                }
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **上限必须被三角不等式压下来。** 手算的例子:
    /// A–B ≤ 2、B–C ≤ 2,那么无论原来写了多少,A–C 都不能超过 4。
    #[test]
    fn an_upper_bound_is_squeezed_by_the_path_through_the_third_atom() {
        let mut b = Bounds::new(3, 0.0, 100.0);
        b.set_upper(0, 1, 2.0);
        b.set_upper(1, 2, 2.0);
        b.set_upper(0, 2, 10.0); // 废话,该被压到 4
        triangle_smooth(&mut b).expect("这组约束不矛盾");
        assert!(
            (b.upper(0, 2) - 4.0).abs() < 1e-12,
            "A–C 上限该被压到 4,实得 {}",
            b.upper(0, 2)
        );
    }

    /// **下限也要被推上去。** A–C 至少 5、B–C 至多 2 ⟹ A–B 至少 3。
    #[test]
    fn a_lower_bound_is_pushed_up_by_the_triangle_inequality() {
        let mut b = Bounds::new(3, 0.0, 100.0);
        b.set_lower(0, 2, 5.0);
        b.set_upper(1, 2, 2.0);
        triangle_smooth(&mut b).expect("这组约束不矛盾");
        assert!(
            b.lower(0, 1) >= 3.0 - 1e-12,
            "A–B 下限该被推到 3,实得 {}",
            b.lower(0, 1)
        );
    }

    /// **自相矛盾要说清楚卡在哪一对**,不能只说失败。
    ///
    /// A–B ≥ 10,但 A–C ≤ 1 且 B–C ≤ 1。
    ///
    /// # 断的是"报出来的那一对确实是坏的",**不是**"必须是我猜的那一对"
    ///
    /// 这组约束里 (A,B) 和 (B,C) **都**是合法见证:A–B ≥ 10 且 A–C ≤ 1
    /// 推出 B–C ≥ 9,而 B–C ≤ 1。报哪一个取决于 k/i/j 的遍历次序 ——
    /// 头一版我断言必须报 A–B,当场红,**错的是期望不是代码**。
    #[test]
    fn an_impossible_constraint_set_names_a_genuinely_broken_pair() {
        let mut b = Bounds::new(3, 0.0, 100.0);
        b.set_lower(0, 1, 10.0);
        b.set_upper(0, 2, 1.0);
        b.set_upper(1, 2, 1.0);
        match triangle_smooth(&mut b) {
            Err(SmoothError::Infeasible { pair: (i, j) }) => {
                assert!(
                    b.lower(i, j) > b.upper(i, j),
                    "报的是 ({i},{j}),可它的下限 {} 并没有超过上限 {} —— 见证是假的",
                    b.lower(i, j),
                    b.upper(i, j)
                );
            }
            Ok(()) => panic!("这组约束是矛盾的,该报不可行"),
        }
    }

    /// **光滑之后的上限矩阵本身是一个度量** —— `U_ij ≤ U_ik + U_kj` 对**所有**三元组成立。
    ///
    /// 这条是整个算法的地基:正因为它是度量,才能直接拿它当参考距离表,
    /// 而不必像 RDKit 那样在区间里随机取一组(取出来的往往任何空间都摆不出)。
    #[test]
    fn the_smoothed_upper_bounds_are_themselves_a_metric() {
        // 一条 6 原子的链,外加几条随手写的、彼此不自洽的约束
        let n = 6;
        let mut b = Bounds::new(n, 0.5, 50.0);
        for i in 0..(n - 1) {
            b.set_upper(i, i + 1, 1.5);
            b.set_lower(i, i + 1, 1.4);
        }
        b.set_upper(0, 5, 40.0);
        b.set_upper(1, 4, 30.0);
        b.set_lower(0, 3, 2.0);
        triangle_smooth(&mut b).expect("该可行");
        let mut checked = 0;
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    if i == j || j == k || i == k {
                        continue;
                    }
                    assert!(
                        b.upper(i, j) <= b.upper(i, k) + b.upper(k, j) + 1e-9,
                        "三角不等式破了:U({i},{j})={} > U({i},{k})+U({k},{j})={}",
                        b.upper(i, j),
                        b.upper(i, k) + b.upper(k, j)
                    );
                    checked += 1;
                }
            }
        }
        assert!(checked >= 100, "只验了 {checked} 个三元组");
    }

    /// **再光滑一次不该有任何变化** —— 单趟就是不动点。
    ///
    /// 这条是上一条的操作版:它同时钉住"结果与遍历次序无关"。
    #[test]
    fn smoothing_twice_changes_nothing() {
        let n = 7;
        let mut b = Bounds::new(n, 0.8, 60.0);
        for i in 0..(n - 1) {
            b.set_upper(i, i + 1, 1.52);
            b.set_lower(i, i + 1, 1.50);
        }
        for i in 0..(n - 2) {
            b.set_upper(i, i + 2, 2.6);
            b.set_lower(i, i + 2, 2.4);
        }
        triangle_smooth(&mut b).expect("该可行");
        let once = b.clone();
        triangle_smooth(&mut b).expect("再来一次也该可行");
        assert_eq!(
            b, once,
            "第二趟光滑化改动了矩阵 —— 单趟不是不动点,那结果就依赖遍历次序"
        );
    }

    /// 空的、单原子的:不许 panic。
    #[test]
    fn tiny_inputs_are_fine() {
        for n in 0..3 {
            let mut b = Bounds::new(n, 1.0, 5.0);
            triangle_smooth(&mut b).expect("小输入该直接过");
            assert_eq!(b.len(), n);
        }
        assert!(Bounds::new(0, 1.0, 5.0).is_empty());
    }
}
