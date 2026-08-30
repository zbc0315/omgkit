//! 置换的宇称 —— **一个概念,一处实现**。
//!
//! 四面体手性、丙二烯轴手性、反应产物的重定基、SMARTS 的参照系换算,做的都是
//! 同一件事:同一组配体的两种排列,差的是奇数次对换还是偶数次。奇数就要把
//! `@` 与 `@@` 对调。
//!
//! 这份先前在四个地方各写了一遍(`omgkit-io` 的 `smiles`、`smarts::mol`,
//! `omgkit-match` 的 `matcher`、`react`)。当时四份给出的答案一致,可**"当前一致"
//! 不是性质,是巧合**:边界条件本来就已经分岔了 —— 有一份长度对不上时返回
//! `false`(按偶处理,不翻),另外三份返回 `None`(说不出来)。四份里改一份,
//! 分子会在一条路上是自己、在另一条路上是对映体,而两边各自自洽、谁都不报错。
//!
//! **判据里那一份不算重复。** `omgkit-io/tests/canonical_invariance.rs` 里有一份
//! 独立实现,那是故意的:判据拿被测代码当真值就什么都没验。

/// `from` → `to` 置换的宇称。两者不是同一多重集时返回 `None`。
///
/// `None` 的意思是"说不出来",**不是"偶"** —— 调用方拿到它时该放弃翻转,
/// 而不是当作"不用翻"。两者的区别在于:前者会让上层知道这个标记表达不出来,
/// 后者会安静地交出一个可能是对映体的分子。
///
/// 配体数 n ≤ 6,O(n²) 完全够用。
#[must_use]
pub fn permutation_is_odd<T: PartialEq + Clone>(from: &[T], to: &[T]) -> Option<bool> {
    if from.len() != to.len() {
        return None;
    }
    let mut cur = from.to_vec();
    let mut swaps = 0usize;
    for (i, want) in to.iter().enumerate() {
        if cur[i] == *want {
            continue;
        }
        let j = (i + 1..cur.len()).find(|&j| cur[j] == *want)?;
        cur.swap(i, j);
        swaps += 1;
    }
    Some(swaps % 2 == 1)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 断的是**性质**:恒等是偶、一次对换是奇、奇偶随对换次数交替。
    #[test]
    fn parity_alternates_with_each_transposition() {
        let base = [10u32, 20, 30, 40];
        assert_eq!(permutation_is_odd(&base, &base), Some(false));
        assert_eq!(permutation_is_odd(&base, &[20, 10, 30, 40]), Some(true));
        assert_eq!(permutation_is_odd(&base, &[20, 30, 10, 40]), Some(false));
        assert_eq!(permutation_is_odd(&base, &[40, 30, 20, 10]), Some(false));
        assert_eq!(permutation_is_odd(&base, &[10, 20, 40, 30]), Some(true));
    }

    /// **不是同一多重集时给 `None`,不给"偶"。**
    ///
    /// 先前四份实现里有一份在这一档返回 `false`(按偶处理)。那一档正是
    /// "这个标记表达不出来",安静地当成"不用翻"会交出可能是对映体的分子。
    #[test]
    fn a_mismatched_multiset_is_unanswerable_not_even() {
        assert_eq!(permutation_is_odd(&[1u32, 2, 3], &[1, 2]), None, "长度不同");
        assert_eq!(
            permutation_is_odd(&[1u32, 2, 3], &[1, 2, 4]),
            None,
            "元素不同"
        );
        assert_eq!(
            permutation_is_odd(&[1u32, 1, 2], &[1, 2, 2]),
            None,
            "重数不同"
        );
        // 空的两边是同一个多重集,恒等,偶
        assert_eq!(permutation_is_odd::<u32>(&[], &[]), Some(false));
    }

    /// 重复元素:同样的两个东西对调是自同构,宇称仍要说得出来。
    #[test]
    fn repeated_elements_still_have_a_parity() {
        assert_eq!(
            permutation_is_odd(&[1u32, 1, 2, 3], &[1, 1, 3, 2]),
            Some(true)
        );
        assert_eq!(
            permutation_is_odd(&[1u32, 1, 2, 3], &[1, 1, 2, 3]),
            Some(false)
        );
    }
}
