//! 配位几何(`@SP` / `@TB` / `@OH`)的排列表与换参照系的换算。
//!
//! # 序号是相对"配体按什么顺序列出"的
//!
//! `[Co@OH25](N)(O)(S)(P)Cl` 里的 25 说的是"按这个列出顺序,六个配体各占八面体的
//! 哪个顶点"。列出顺序一变,同一个分子的序号就变 —— 而解析(书写序 → 存储序)与
//! 写出(存储序 → 输出序)两处都会变顺序。不换算就等于换了个分子。
//!
//! # 表是从参照实现穷举量出来的,不是照着规范条文写的
//!
//! 对每一类:取**互不相同**的配体(相同的配体会让不同序号指向同一个分子,
//! 把规则盖住),穷举"每个序号 × 每种列出顺序"的全部写法,交给
//! RDKit 2025.09.2 规范化,看哪些写法落到同一个分子:
//!
//! | 类别 | 配体数 | 序号数 | 写法数 | 归成几组 | 转动群阶 |
//! |---|---|---|---|---|---|
//! | `@SP` | 4 | 3 | 3 × 24 = **72** | 3 | 8 |
//! | `@TB` | 5 | 20 | 20 × 120 = **2400** | 20 | 6 |
//! | `@OH` | 6 | 30 | 30 × 720 = **21600** | 30 | 24 |
//!
//! 组数恰好等于序号数,每组大小恰好等于顺序数 —— 序号与立体异构体一一对应。
//! 转动群的阶也对得上几何:方形的对称群 D4 是 8(平面四方没有手性,镜像还是
//! 它自己),三角双锥是 6,八面体是 24。
//!
//! # 模型:把"序号 1 在恒等顺序下"的那些位置定义成多面体的顶点
//!
//! 于是两样东西都直接从测量里读出来,不用手写几何:
//!
//! - **转动群 `R`** = 序号 1 的稳定子(哪些顺序重排之后还是同一个分子);
//! - **每个序号的顶点对应 `p_i`**,满足 `isomer(i, σ) = isomer(1, σ∘p_i)`。
//!
//! 模型在**全部写法**上验过:72 / 2400 / 21600 条,反例 0。外加陪集判据
//! `isomer(1,α) = isomer(1,β) ⟺ ∃r∈R: α = β∘r`,分别在 576 / 14400 / 518400 对上
//! 验过,反例 0(其中真的同构的分别是 192 / 720 / 17280 对 —— 不是空过)。
//!
//! 这两条性质在 `tests/` 里有本地版本(不需要 RDKit),端到端那条在
//! `harness/check_stereo_perm.py`(需要 RDKit,在 CI 里跑)。
//!
//! # 换算怎么做
//!
//! 要把序号 `i` 从顺序 `from` 换到顺序 `to`:isomer 相同即
//! `to∘p_j = (from∘p_i)∘r` 对某个 `r ∈ R` 成立。左乘 `to⁻¹` 得
//! `p_j = q∘r`,其中 `q = to⁻¹∘from∘p_i`。因为那些 `p_j` 恰好是各个右陪集的
//! 代表元,满足条件的 `j` 有且只有一个。

use crate::ChiralTag;

/// 一类配位几何的排列表。
struct Table {
    /// 转动群,每个元素是顶点下标的一个置换
    rotations: &'static [&'static [u8]],
    /// 序号 → "顶点 → 列出位置"的对应。下标 0 对应序号 1。
    slots: &'static [&'static [u8]],
}

const SP_ROTATIONS: &[&[u8]] = &[
    &[0, 1, 2, 3],
    &[0, 3, 2, 1],
    &[1, 0, 3, 2],
    &[1, 2, 3, 0],
    &[2, 1, 0, 3],
    &[2, 3, 0, 1],
    &[3, 0, 1, 2],
    &[3, 2, 1, 0],
];

const SP_SLOTS: &[&[u8]] = &[&[0, 1, 2, 3], &[0, 2, 1, 3], &[0, 1, 3, 2]];

const TB_ROTATIONS: &[&[u8]] = &[
    &[0, 1, 2, 3, 4],
    &[0, 2, 3, 1, 4],
    &[0, 3, 1, 2, 4],
    &[4, 1, 3, 2, 0],
    &[4, 2, 1, 3, 0],
    &[4, 3, 2, 1, 0],
];

const TB_SLOTS: &[&[u8]] = &[
    &[0, 1, 2, 3, 4],
    &[0, 1, 3, 2, 4],
    &[0, 1, 2, 4, 3],
    &[0, 1, 4, 2, 3],
    &[0, 1, 3, 4, 2],
    &[0, 1, 4, 3, 2],
    &[0, 2, 3, 4, 1],
    &[0, 2, 4, 3, 1],
    &[1, 0, 2, 3, 4],
    &[1, 0, 2, 4, 3],
    &[1, 0, 3, 2, 4],
    &[1, 0, 4, 2, 3],
    &[1, 0, 3, 4, 2],
    &[1, 0, 4, 3, 2],
    &[2, 0, 1, 3, 4],
    &[2, 0, 1, 4, 3],
    &[3, 0, 1, 2, 4],
    &[3, 0, 2, 1, 4],
    &[2, 0, 4, 1, 3],
    &[2, 0, 3, 1, 4],
];

const OH_ROTATIONS: &[&[u8]] = &[
    &[0, 1, 2, 3, 4, 5],
    &[0, 2, 3, 4, 1, 5],
    &[0, 3, 4, 1, 2, 5],
    &[0, 4, 1, 2, 3, 5],
    &[1, 0, 4, 5, 2, 3],
    &[1, 2, 0, 4, 5, 3],
    &[1, 4, 5, 2, 0, 3],
    &[1, 5, 2, 0, 4, 3],
    &[2, 0, 1, 5, 3, 4],
    &[2, 1, 5, 3, 0, 4],
    &[2, 3, 0, 1, 5, 4],
    &[2, 5, 3, 0, 1, 4],
    &[3, 0, 2, 5, 4, 1],
    &[3, 2, 5, 4, 0, 1],
    &[3, 4, 0, 2, 5, 1],
    &[3, 5, 4, 0, 2, 1],
    &[4, 0, 3, 5, 1, 2],
    &[4, 1, 0, 3, 5, 2],
    &[4, 3, 5, 1, 0, 2],
    &[4, 5, 1, 0, 3, 2],
    &[5, 1, 4, 3, 2, 0],
    &[5, 2, 1, 4, 3, 0],
    &[5, 3, 2, 1, 4, 0],
    &[5, 4, 3, 2, 1, 0],
];

const OH_SLOTS: &[&[u8]] = &[
    &[0, 1, 2, 3, 4, 5],
    &[0, 1, 4, 3, 2, 5],
    &[0, 1, 2, 3, 5, 4],
    &[0, 1, 2, 4, 3, 5],
    &[0, 1, 2, 5, 3, 4],
    &[0, 1, 2, 4, 5, 3],
    &[0, 1, 2, 5, 4, 3],
    &[0, 1, 3, 2, 4, 5],
    &[0, 1, 3, 2, 5, 4],
    &[0, 1, 4, 2, 3, 5],
    &[0, 1, 5, 2, 3, 4],
    &[0, 1, 4, 2, 5, 3],
    &[0, 1, 5, 2, 4, 3],
    &[0, 1, 3, 4, 2, 5],
    &[0, 1, 3, 5, 2, 4],
    &[0, 1, 5, 3, 2, 4],
    &[0, 1, 4, 5, 2, 3],
    &[0, 1, 5, 4, 2, 3],
    &[0, 1, 3, 4, 5, 2],
    &[0, 1, 3, 5, 4, 2],
    &[0, 1, 4, 3, 5, 2],
    &[0, 1, 5, 3, 4, 2],
    &[0, 1, 4, 5, 3, 2],
    &[0, 1, 5, 4, 3, 2],
    &[0, 2, 3, 4, 5, 1],
    &[0, 2, 3, 5, 4, 1],
    &[0, 2, 4, 3, 5, 1],
    &[0, 2, 5, 3, 4, 1],
    &[0, 2, 4, 5, 3, 1],
    &[0, 2, 5, 4, 3, 1],
];

fn table_for(tag: ChiralTag) -> Option<Table> {
    match tag {
        ChiralTag::SquarePlanar => Some(Table {
            rotations: SP_ROTATIONS,
            slots: SP_SLOTS,
        }),
        ChiralTag::TrigonalBipyramidal => Some(Table {
            rotations: TB_ROTATIONS,
            slots: TB_SLOTS,
        }),
        ChiralTag::Octahedral => Some(Table {
            rotations: OH_ROTATIONS,
            slots: OH_SLOTS,
        }),
        _ => None,
    }
}

/// 这类配位几何要几个配体。四面体与其它类别没有排列表,给 `None`。
#[must_use]
pub fn ligand_count(tag: ChiralTag) -> Option<usize> {
    table_for(tag).map(|t| t.slots[0].len())
}

/// 把配位几何的排列序号从一种配体顺序换算到另一种。
///
/// `from` / `to` 是同一组配体的两种排列(用什么标识都行 —— 键下标、原子下标,
/// 只要在这一组里唯一)。返回 `to` 那个顺序下表达同一个分子的序号。
///
/// 多面体的某个顶点不是键的时候(方括号里的氢、一个空的配位位置),给它一个
/// 占位标识、跟着一起传进来即可。它在两侧各排在哪一位由调用方定,见
/// `omgkit-io` 的 `smiles::coordination_ligands`。
///
/// # 返回 `None` 的情形
///
/// 类别不在这三种配位几何之内(四面体、丙二烯轴手性、未标注 —— 前两类各有
/// 各的机制)、序号不在 `1..=序号数`、
/// 两侧配体数与该几何对不上、或者两侧不是同一组配体。调用方应当把这几种
/// 都当作"这个标记表达不出来",而不是猜一个值 —— 猜出来的是另一个分子。
#[must_use]
pub fn renumber(tag: ChiralTag, perm: u8, from: &[u32], to: &[u32]) -> Option<u8> {
    let t = table_for(tag)?;
    let n = t.slots[0].len();
    if from.len() != n || to.len() != n {
        return None;
    }
    let p_i = *t.slots.get(usize::from(perm).checked_sub(1)?)?;

    // q = to⁻¹ ∘ from ∘ p_i:第 k 个顶点上的配体,在 `to` 里排第几
    let mut q = [0u8; 6];
    for k in 0..n {
        let ligand = *from.get(usize::from(p_i[k]))?;
        let pos = to.iter().position(|&x| x == ligand)?;
        q[k] = u8::try_from(pos).ok()?;
    }

    // 找那个落在同一个右陪集里的序号:p_j = q∘r
    for (j, p_j) in t.slots.iter().enumerate() {
        if t.rotations
            .iter()
            .any(|r| (0..n).all(|k| usize::from(p_j[k]) == usize::from(q[usize::from(r[k])])))
        {
            return u8::try_from(j + 1).ok();
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `0..n` 的全部排列。
    fn permutations(n: usize) -> Vec<Vec<u32>> {
        let mut out = Vec::new();
        let mut cur: Vec<u32> = (0..n as u32).collect();
        fn go(k: usize, cur: &mut Vec<u32>, out: &mut Vec<Vec<u32>>) {
            if k == cur.len() {
                out.push(cur.clone());
                return;
            }
            for i in k..cur.len() {
                cur.swap(k, i);
                go(k + 1, cur, out);
                cur.swap(k, i);
            }
        }
        go(0, &mut cur, &mut out);
        out
    }

    fn factorial(n: usize) -> usize {
        (1..=n).product()
    }

    const ALL: [ChiralTag; 3] = [
        ChiralTag::SquarePlanar,
        ChiralTag::TrigonalBipyramidal,
        ChiralTag::Octahedral,
    ];

    /// **表的结构就是判据**:序号数 × 转动群阶 = 顺序数,而且各个序号落在
    /// 互不相同的右陪集里。
    ///
    /// 这两条一起把表钉死:少一个序号、多一个转动、某两行落进同一个陪集,
    /// 都当场红。它们也正是"序号与立体异构体一一对应"这句话的形式化。
    #[test]
    fn 排列表是转动群的一组陪集代表元() {
        for tag in ALL {
            let t = table_for(tag).expect("有表");
            let n = t.slots[0].len();
            assert_eq!(
                t.slots.len() * t.rotations.len(),
                factorial(n),
                "{tag:?}:序号数 {} × 转动群阶 {} ≠ {n}! = {}",
                t.slots.len(),
                t.rotations.len(),
                factorial(n)
            );
            // 任意两个序号不在同一个右陪集里
            for (a, pa) in t.slots.iter().enumerate() {
                for (b, pb) in t.slots.iter().enumerate().skip(a + 1) {
                    let same = t
                        .rotations
                        .iter()
                        .any(|r| (0..n).all(|k| pa[k] == pb[usize::from(r[k])]));
                    assert!(
                        !same,
                        "{tag:?}:序号 {} 与 {} 落在同一个陪集里",
                        a + 1,
                        b + 1
                    );
                }
            }
        }
    }

    /// 转动群本身要是个群:含恒等、对复合封闭、每个元素有逆。
    ///
    /// 表是量出来的,而"量出来的东西是不是群"没有任何东西保证 —— 这里补上。
    #[test]
    fn 转动群是个群() {
        for tag in ALL {
            let t = table_for(tag).expect("有表");
            let n = t.slots[0].len();
            let ident: Vec<u8> = (0..n as u8).collect();
            assert!(
                t.rotations.iter().any(|r| r == &&ident[..]),
                "{tag:?}:转动群里没有恒等"
            );
            for a in t.rotations {
                for b in t.rotations {
                    let ab: Vec<u8> = (0..n).map(|k| a[usize::from(b[k])]).collect();
                    assert!(
                        t.rotations.iter().any(|r| r == &&ab[..]),
                        "{tag:?}:转动群对复合不封闭"
                    );
                }
            }
        }
    }

    /// 换算是可逆的,而且对固定的目标顺序是个双射。
    ///
    /// 双射这一条挡的是"把好几个序号映到同一个"那类退化 —— 一个恒返回 1 的
    /// 实现能过"可逆"(不能),但过不了双射。
    #[test]
    fn 换算可逆且是双射() {
        for tag in ALL {
            let t = table_for(tag).expect("有表");
            let n = t.slots[0].len();
            let ident: Vec<u32> = (0..n as u32).collect();
            for to in permutations(n) {
                let mut seen = vec![false; t.slots.len()];
                for i in 1..=t.slots.len() as u8 {
                    let j = renumber(tag, i, &ident, &to)
                        .unwrap_or_else(|| panic!("{tag:?}:序号 {i} 换到 {to:?} 换不出来"));
                    assert!(
                        !seen[usize::from(j) - 1],
                        "{tag:?}:两个序号都换成了 {j},换算不是双射"
                    );
                    seen[usize::from(j) - 1] = true;
                    assert_eq!(
                        renumber(tag, j, &to, &ident),
                        Some(i),
                        "{tag:?}:序号 {i} 换过去再换回来不是自己"
                    );
                }
            }
        }
    }

    /// 平面四方那张表要与**人读得懂的那条规则**对得上:序号说的是
    /// "按列出顺序哪两对配体互为反位"。
    ///
    /// 规则只写在这里(判据侧),实现走的是通用的陪集机制 —— 两条路互相对账。
    /// 这也是当初量出通用表之前先量出来的那条规则:`@SP1` 反位 (1,3)(2,4)、
    /// `@SP2` (1,2)(3,4)、`@SP3` (1,4)(2,3),与规范里 U / 4 / Z 三种形状对应。
    #[test]
    fn 平面四方的表与反位配对规则对得上() {
        const TRANS: [[(usize, usize); 2]; 3] =
            [[(0, 2), (1, 3)], [(0, 1), (2, 3)], [(0, 3), (1, 2)]];
        fn pairing(perm: u8, ligands: &[u32]) -> Vec<[u32; 2]> {
            let mut v: Vec<[u32; 2]> = TRANS[usize::from(perm) - 1]
                .iter()
                .map(|&(i, j)| {
                    let (a, b) = (ligands[i], ligands[j]);
                    if a <= b {
                        [a, b]
                    } else {
                        [b, a]
                    }
                })
                .collect();
            v.sort_unstable();
            v
        }
        let ident: Vec<u32> = (0..4).collect();
        for to in permutations(4) {
            for i in 1..=3u8 {
                let j = renumber(ChiralTag::SquarePlanar, i, &ident, &to).expect("换得出来");
                assert_eq!(
                    pairing(i, &ident),
                    pairing(j, &to),
                    "序号 {i} 换到顺序 {to:?} 得到 {j},反位配对却变了"
                );
            }
        }
    }

    /// 表达不出来的输入一律给 `None` —— 绝不猜一个值。
    #[test]
    fn 换算不了时不猜() {
        let l: Vec<u32> = (0..4).collect();
        assert_eq!(renumber(ChiralTag::SquarePlanar, 0, &l, &l), None, "序号 0");
        assert_eq!(
            renumber(ChiralTag::SquarePlanar, 4, &l, &l),
            None,
            "序号越界"
        );
        assert_eq!(renumber(ChiralTag::Cw, 1, &l, &l), None, "四面体没有排列表");
        assert_eq!(
            renumber(ChiralTag::Allene, 1, &l, &l),
            None,
            "@AL 没有排列表"
        );
        assert_eq!(
            renumber(ChiralTag::SquarePlanar, 1, &l[..3], &l[..3]),
            None,
            "配体数与几何对不上"
        );
        assert_eq!(
            renumber(ChiralTag::SquarePlanar, 1, &l, &[0, 1, 2, 9]),
            None,
            "两侧不是同一组配体"
        );
        assert_eq!(ligand_count(ChiralTag::Cw), None);
        assert_eq!(ligand_count(ChiralTag::Octahedral), Some(6));
    }
}
