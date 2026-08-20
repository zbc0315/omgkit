//! **确定性地把重合的原子分开** —— 必须在优化器**之前**跑。
//!
//! # 为什么非做不可
//!
//! 对称分子的 Gram 矩阵有**重特征值**,而对称等价的原子在那些特征子空间里
//! 拿到的分量可以完全相同 —— 于是嵌入给出**逐位相同的坐标**。
//!
//! 那不是"挤得有点近",是精确重合,而精确重合的两个原子**梯度恰好为零**:
//! 距离罚项的梯度是 `dE/dd · (rᵢ − rⱼ)/d`,分子是零向量。优化器再跑一万步也
//! 分不开它们。RDKit 靠随机数天然避开了这件事,我们全程无随机数,躲不掉。
//!
//! **实测这不是理论顾虑:**
//!
//! | 语料 | 分子 | 原子完全重合(< 1e-6 Å)的分子 |
//! |---|---|---|
//! | `hard.smi` | 68 | **8(11.8%)** |
//! | `large.smi` | 8830 | **44(0.50%)** |
//!
//! 0.50% 已经与 RDKit 的整体失败率同量级,而且是**静默**的 —— 坐标照样返回,
//! 只是废的。典型案例全是四个相同取代基的四面体中心:
//! `C[Si](C)(C)C`、`C[Ge](C)(C)C`、`C[Sn](C)(C)C`、`C[P+](C)(C)C`。
//!
//! # 判据先前为什么没抓到
//!
//! 端到端判官跑的是 `smoke.chirality.jsonl` —— 那是**带手性中心**的分子,
//! 而手性中心本身就破了对称。**判据的输入分布系统性地排除了要测的那一档。**
//! 这与 `basinThresh` 那次是同一个病:判据看着在守,实际什么都没守。
//!
//! # 怎么破:van der Corput,不是随机数
//!
//! 每个原子按**自己的下标**取一个确定的位移。用 van der Corput 基数反转序列
//! (基 2/3/5 各给一个分量)—— 它确定、低差异、下标不同则位移不同,
//! 正好满足"把简并的对称拆开"又"同一分子永远同一答案"。
//!
//! 位移只给**需要的那些原子**,幅度也只要够让梯度有方向 —— 剩下的交给精修,
//! 它会把结构拉回界内。

/// 认为两个原子"重合"的阈值(Å)。
///
/// 取 1e-6:真正的简并给出的是**逐位相同**的坐标(实测就是 0.000000),
/// 而正常结构里最近的非键原子对也在 1 Å 以上 —— 中间空了六个数量级,
/// 阈值放在哪儿都一样,不存在误判。
pub const COINCIDENT_TOL: f64 = 1e-6;

/// 位移的幅度(Å)。
///
/// 只要够让梯度有个方向就行,剩下的交给精修。取 0.1 Å:
/// 比键长小一个数量级,不会把结构推乱;又远大于 [`COINCIDENT_TOL`],
/// 一次就分得开。
pub const AMPLITUDE: f64 = 0.1;

/// van der Corput 基数反转:把 `n` 的 `base` 进制展开翻到小数点后。
///
/// 结果落在 `[0, 1)`,低差异 —— 相邻下标给出的值分得很开,
/// 这正是"把简并拆开"要的性质。**完全确定**,没有随机数。
#[must_use]
pub fn van_der_corput(mut n: usize, base: usize) -> f64 {
    debug_assert!(base >= 2, "基数至少是 2");
    let mut q = 0.0;
    let mut bk = 1.0 / base as f64;
    while n > 0 {
        q += (n % base) as f64 * bk;
        n /= base;
        bk /= base as f64;
    }
    q
}

/// 把重合的原子确定性地分开,返回动了几个原子。
///
/// 只动**确实与别人重合**的那些 —— 没有重合就一个都不碰,
/// 这样绝大多数分子走这一步是零开销、零影响。
pub fn break_coincidence(coords: &mut [[f64; 3]]) -> usize {
    let n = coords.len();
    // 先标记:谁与谁重合。O(n²),而 n 是原子数,可忽略。
    let mut needs = vec![false; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let d2 = (0..3)
                .map(|t| (coords[i][t] - coords[j][t]).powi(2))
                .sum::<f64>();
            if d2 <= COINCIDENT_TOL * COINCIDENT_TOL {
                needs[i] = true;
                needs[j] = true;
            }
        }
    }
    let mut moved = 0;
    for (i, flag) in needs.iter().enumerate() {
        if !flag {
            continue;
        }
        // 三个分量取不同的基,位移就散在三维里而不是一条线上。
        //
        // **但共线其实也能分开** —— 变异验证过:三个分量都用基 2 时判据仍然全绿,
        // 因为不同下标的幅度本来就不同,沿同一方向推不同的距离照样分得开。
        // 用三个基是为了扰动的条件数更好,不是正确性必需 —— 别把它说成后者。
        let off = [
            van_der_corput(i + 1, 2) - 0.5,
            van_der_corput(i + 1, 3) - 0.5,
            van_der_corput(i + 1, 5) - 0.5,
        ];
        for t in 0..3 {
            coords[i][t] += AMPLITUDE * off[t];
        }
        moved += 1;
    }
    moved
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn van_der_corput_的前几项() {
        // 基 2:1/2, 1/4, 3/4, 1/8, 5/8, 3/8, 7/8 …
        let want = [0.5, 0.25, 0.75, 0.125, 0.625, 0.375, 0.875];
        for (k, w) in want.iter().enumerate() {
            let got = van_der_corput(k + 1, 2);
            assert!((got - w).abs() < 1e-15, "第 {} 项:{got} vs {w}", k + 1);
        }
        // 基 3:1/3, 2/3, 1/9 …
        assert!((van_der_corput(1, 3) - 1.0 / 3.0).abs() < 1e-15);
        assert!((van_der_corput(2, 3) - 2.0 / 3.0).abs() < 1e-15);
        assert!((van_der_corput(3, 3) - 1.0 / 9.0).abs() < 1e-15);
        // 0 落在 0
        assert_eq!(van_der_corput(0, 2), 0.0);
    }

    #[test]
    fn 相邻下标的位移互不相同() {
        // 这是"能把简并拆开"的前提:下标不同 → 位移不同。
        let mut seen: Vec<[f64; 3]> = Vec::new();
        for i in 1..=64 {
            let o = [
                van_der_corput(i, 2),
                van_der_corput(i, 3),
                van_der_corput(i, 5),
            ];
            for prev in &seen {
                let same = (0..3).all(|t| (prev[t] - o[t]).abs() < 1e-12);
                assert!(!same, "下标 {i} 与前面某个撞了:{o:?}");
            }
            seen.push(o);
        }
    }

    #[test]
    fn 不重合就一个都不动() {
        let mut c = vec![[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [0.0, 1.5, 0.0]];
        let before = c.clone();
        assert_eq!(break_coincidence(&mut c), 0);
        assert_eq!(c, before, "没有重合时不许碰坐标");
    }

    #[test]
    fn 重合的会被分开() {
        // 四个原子挤在一个点上 —— 正是四取代四面体那种简并的样子
        let mut c = vec![[0.0; 3]; 4];
        c.push([5.0, 0.0, 0.0]); // 一个正常的,不该被动
        let moved = break_coincidence(&mut c);
        assert_eq!(moved, 4, "四个重合的都该动");
        assert_eq!(c[4], [5.0, 0.0, 0.0], "不重合的那个不许动");
        // 分开之后两两都不再重合
        for i in 0..4 {
            for j in (i + 1)..4 {
                let d = (0..3)
                    .map(|t| (c[i][t] - c[j][t]).powi(2))
                    .sum::<f64>()
                    .sqrt();
                assert!(d > COINCIDENT_TOL, "{i}/{j} 还是重合:{d}");
            }
        }
    }

    #[test]
    fn 位移是确定的() {
        let build = || {
            let mut c = vec![[0.0; 3]; 6];
            break_coincidence(&mut c);
            c
        };
        assert_eq!(build(), build(), "两次跑出来的位移必须逐位相同");
    }

    #[test]
    fn 幅度受控() {
        // 位移不能大到把结构推乱 —— 每个分量最多 AMPLITUDE/2
        let mut c = vec![[0.0; 3]; 32];
        break_coincidence(&mut c);
        for (i, p) in c.iter().enumerate() {
            for (t, v) in p.iter().enumerate() {
                assert!(
                    v.abs() <= AMPLITUDE / 2.0 + 1e-12,
                    "第 {i} 个原子第 {t} 分量位移 {v} 超过 {}",
                    AMPLITUDE / 2.0
                );
            }
        }
    }
}
