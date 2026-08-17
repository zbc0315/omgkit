//! 把隐式氢补成真正的原子。三维构象要用。
//!
//! # 与 [`remove_hs`](mod@crate::remove_hs) **不对称**,别照着抄
//!
//! 名字看着是一对,行为不是:
//!
//! | | `remove_hs` | 本模块 |
//! |---|---|---|
//! | 做法 | 删原子 ⇒ 必须**重建**整个分子 | 纯**追加** |
//! | 原有下标 | **全部失效** | **一个不变** |
//! | 邻居相对序 | 保持(重建时按原键序) | 保持,新氢排在**最后** |
//!
//! "原有下标一个不变"是三维这边最需要的性质:坐标、秩、判据全都按下标对应,
//! 补个氢就重排下标的话,上下游得整体重算。所以这里**不许**为了对称把它写成
//! 就地重建。
//!
//! # 补的次序按秩,不按存储序
//!
//! 同一个分子换种写法,存储序会变;新氢若按存储序追加,它们拿到的原子号就跟着
//! 变,后面整条管线(坐标、指纹、判据)全跟着变。所以中心按**传进来的秩**排序,
//! 秩相同(不可能,秩是双射,但防一手)再按下标。
//!
//! 秩是**参数**而不是在这里算的:`omgkit-chem` 不依赖 `omgkit-io`,而规范秩住在
//! 那边。把依赖写在签名上,与 [`adjust_hs`](fn@crate::adjust_hs) 把
//! `implicit_before` 写成参数是同一个道理 —— 隐式依赖迟早被人绕过。
//!
//! # 手性标记的参照顺序会变,**这一版不管**
//!
//! 四面体手性是相对"邻居的排列顺序"定义的,而新氢一律排在最后。原先那个隐式氢
//! 在 SMILES 里的位置(紧跟中心原子,或中心是首原子时排第一)与"最后"一般不是
//! 同一个位置,**位置一挪,奇偶可能翻,`@` 与 `@@` 就反了**。
//!
//! 本函数**不去修正这个**,也不假装修正了:它只补原子与键,不碰
//! [`chiral_tag`](omgkit_core::AtomData::chiral_tag)。谁要在补氢之后用手性,
//! 必须自己把这一层想清楚 —— 有判据把"补氢不碰 tag"这件事钉住,免得将来有人
//! 以为它管了。

use omgkit_core::{BondOrder, MolBuilder};

/// 把隐式氢(以及方括号里记着的显式氢计数)补成真正的氢原子。
///
/// `order` 是与写法无关的原子秩,长度必须与原子数一致 —— 补出来的氢按它排序,
/// 见模块文档。返回补了几个氢。
///
/// **在 `sanitize` 之后调用。** 净化会重算隐式氢计数,补完再净化会把补出来的
/// 氢又算一遍(中心的计数已经清零,所以不会重复补,但价键检查会看到不同的分子)。
///
/// # 原有下标一个不变
///
/// 新氢追加在原子表末尾,新键追加在键表末尾。调用前拿到的任何原子/键下标在
/// 返回之后仍然有效 —— 有判据钉住。
pub fn add_explicit_hs(mol: &mut MolBuilder, order: &[u32]) -> usize {
    let n = mol.num_atoms();
    debug_assert_eq!(order.len(), n, "秩与原子数必须一一对应");
    if order.len() != n {
        return 0;
    }

    // 先把"谁要补几个"收集齐,再按秩排 —— 边遍历边加原子会把 `num_atoms` 搅乱
    let mut todo: Vec<(u32, u8)> = Vec::new();
    for a in 0..u32::try_from(n).unwrap_or(u32::MAX) {
        let at = mol.atoms()[a as usize];
        let k = at.num_implicit_hs.saturating_add(at.num_explicit_hs);
        if k > 0 {
            todo.push((a, k));
        }
    }
    todo.sort_by_key(|(a, _)| (order[*a as usize], *a));

    let mut added = 0usize;
    for (a, k) in todo {
        for _ in 0..k {
            let h = mol.add_atom(1);
            // 加键失败只可能是下标越界,而 `h` 刚由本函数造出来 —— 不会发生。
            // 真发生了也不 panic:少一根键,后面的判据会看见。
            if mol.add_bond(a, h, BondOrder::Single).is_ok() {
                added += 1;
            }
        }
        // 氢已经是真原子了,计数清零,否则它们会被数两遍
        if let Some(at) = mol.atom_mut(a) {
            at.num_implicit_hs = 0;
            at.num_explicit_hs = 0;
        }
    }
    added
}

#[cfg(test)]
mod tests {
    use super::*;

    fn prep(smi: &str) -> MolBuilder {
        let mut m = omgkit_io::smiles::parse(smi).expect("测试用的 SMILES 该能解析");
        crate::pipeline::sanitize(&mut m).expect("测试用的分子该能净化");
        m
    }

    fn ranks(m: &MolBuilder) -> Vec<u32> {
        omgkit_io::canon::classed_ranks(m)
    }

    #[test]
    fn the_hydrogen_count_comes_out_right() {
        for (smi, want) in [
            ("C", 4),
            ("CC", 6),
            ("CCO", 6),
            ("c1ccccc1", 6),
            // 吡咯是 C4H5N:四个环碳各 1 个氢 + 氮上 1 个 = **5**,不是 6。
            // 我按 6 写,判据当场红 —— 是我算错。
            ("[nH]1cccc1", 5),
            ("CC(=O)OC", 6),
            ("[Na+].[Cl-]", 0),
        ] {
            let mut m = prep(smi);
            let r = ranks(&m);
            let before = m.num_atoms();
            let got = add_explicit_hs(&mut m, &r);
            assert_eq!(got, want, "{smi} 补出来的氢数");
            assert_eq!(m.num_atoms(), before + want, "{smi} 的原子总数");
            // 补出来的确实都是氢
            for i in before..m.num_atoms() {
                assert_eq!(m.atoms()[i].atomic_num, 1, "{smi} 第 {i} 个新原子不是氢");
            }
        }
    }

    #[test]
    fn every_original_index_still_means_the_same_atom() {
        // **这是本模块存在的全部理由。** 三维那边坐标、秩、判据都按下标对应,
        // 补个氢就重排下标的话上下游得整体重算(`remove_hs` 就是那样,所以它
        // 明确写着"下标全部失效")。
        let mut m = prep("CC(=O)OCC");
        let r = ranks(&m);
        let before: Vec<(u8, i8)> = m
            .atoms()
            .iter()
            .map(|a| (a.atomic_num, a.formal_charge))
            .collect();
        let bonds_before: Vec<(u32, u32, BondOrder)> = m
            .bonds()
            .iter()
            .map(|b| (b.begin, b.end, b.order))
            .collect();

        add_explicit_hs(&mut m, &r);

        for (i, want) in before.iter().enumerate() {
            let got = (m.atoms()[i].atomic_num, m.atoms()[i].formal_charge);
            assert_eq!(got, *want, "第 {i} 个原子变了");
        }
        for (i, want) in bonds_before.iter().enumerate() {
            let b = m.bonds()[i];
            assert_eq!((b.begin, b.end, b.order), *want, "第 {i} 根键变了");
        }
    }

    #[test]
    fn the_new_hydrogens_are_ordered_by_rank_not_by_storage() {
        // 换种写法,补出来的氢必须拿到"同样的"原子号 —— 否则后面整条管线
        // (坐标、指纹、判据)全跟着写法变。
        //
        // 判法:补完之后,按秩排好的"每个新氢的父亲的秩"序列必须一致。
        let mut seen: Option<Vec<u32>> = None;
        for smi in ["CCO", "OCC", "C(O)C"] {
            let mut m = prep(smi);
            let r = ranks(&m);
            let before = m.num_atoms();
            add_explicit_hs(&mut m, &r);
            // 新氢按追加顺序排列,取它们各自父亲的秩
            let seq: Vec<u32> = (before..m.num_atoms())
                .map(|h| {
                    let p = m
                        .neighbors(u32::try_from(h).unwrap())
                        .next()
                        .expect("氢总有一个父亲")
                        .0;
                    r[p as usize]
                })
                .collect();
            match &seen {
                None => seen = Some(seq),
                Some(first) => assert_eq!(first, &seq, "{smi} 补氢的次序与别的写法不同"),
            }
        }
    }

    #[test]
    fn adding_twice_adds_nothing_the_second_time() {
        // 计数清零了才不会重复补。没清零的话第二次会再补一遍,原子数翻倍。
        let mut m = prep("CCO");
        let r = ranks(&m);
        let first = add_explicit_hs(&mut m, &r);
        assert_eq!(first, 6);
        let r2 = ranks(&m);
        let second = add_explicit_hs(&mut m, &r2);
        assert_eq!(second, 0, "第二次又补了 {second} 个氢");
    }

    #[test]
    fn the_chiral_tag_is_left_alone_on_purpose() {
        // 补氢**不修正**手性标记的参照顺序(新氢排在最后,而原先那个隐式氢在
        // SMILES 里的位置一般不是最后)。模块文档写明了这一点;这条判据把
        // "确实没碰" 钉住,免得将来有人以为它管了。
        let mut m = prep("N[C@@H](C)C(=O)O");
        let r = ranks(&m);
        let tags: Vec<_> = m.atoms().iter().map(|a| a.chiral_tag).collect();
        add_explicit_hs(&mut m, &r);
        for (i, want) in tags.iter().enumerate() {
            assert_eq!(
                m.atoms()[i].chiral_tag,
                *want,
                "第 {i} 个原子的手性标记被改了"
            );
        }
    }
}
