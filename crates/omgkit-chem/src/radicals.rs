//! 自由基电子数(净化第 6 步)。
//!
//! # 为什么排在 kekulize 之后、芳香性感知之前
//!
//! 自由基数由**具体的 Kekulé 结构**决定,芳香写法下这个信息根本不存在:
//! `[N]1C=CC=C1` 里的氮带一个自由基,而同一个分子写成 `[n]1cccc1` 就无从
//! 判断。所以必须等 kekulize 定下键级之后再算。
//!
//! 而芳香性感知又要用到自由基数,所以本步排在它前面。
//!
//! 这与[隐式氢的推断](crate::valence)正好相反 —— 后者是在 kekulize
//! **之前**就预判了结果。
//!
//! # 只处理方括号里写出来的原子
//!
//! 判据是 [`AtomFlags::NO_IMPLICIT`](omgkit_core::AtomFlags::NO_IMPLICIT),
//! 即原子写在 `[]` 中。理由很直接:没写方括号的原子会自动补隐式氢把价填满,
//! 不可能剩下未成对电子;写了方括号就等于声明"氢数我说了算",剩下的空位
//! 才是自由基。通配原子(`*`,Z = 0)也跳过 —— 身份未知,不能猜。
//!
//! # 三个容易写错的地方
//!
//! **1. 总价用截断,不是四舍五入。**
//! 这里是 `(accum + 0.1) as i32`,而[显式价](crate::valence)那边是
//! `round(accum + 0.1)`。两者对 `accum = 1.5` 分别给 1 和 2。本步跑在
//! kekulize 之后,正常不该再有 1.5 的键级,但边角输入上这个区别是可观察的。
//!
//! **2. 价表取自**真实**原子序数,不是"有效原子序数"。**
//! 与价键计算不同,这里 **不**先扣掉形式电荷再查表;电荷作为独立项进算式。
//! 搞混会让带电原子的自由基数系统性偏移。
//!
//! **3. 两条路线取小,但只在第二条非负时。**
//! 按八隅体反推与按外层电子数反推,分别适用于电负性较大和较小的元素;
//! 取二者较小值,**前提是后者非负**。

use omgkit_core::{element, AtomFlags, MolBuilder};

/// 计算并就地写回每个原子的自由基电子数。
///
/// 本步不会失败 —— 它只读结构、只写 `num_radical_electrons`。
pub fn assign_radicals(mol: &mut MolBuilder) {
    for i in 0..mol.num_atoms() as u32 {
        if let Some(n) = radicals_of(mol, i) {
            if let Some(a) = mol.atom_mut(i) {
                a.num_radical_electrons = n;
            }
        }
    }
}

/// 单个原子的自由基电子数;返回 `None` 表示该原子不在处理范围内(保持原值)。
fn radicals_of(mol: &MolBuilder, idx: u32) -> Option<u8> {
    let atom = mol.atoms()[idx as usize];
    let z = atom.atomic_num;

    // 只处理方括号原子,且跳过通配原子
    if !atom.flags.contains(AtomFlags::NO_IMPLICIT) || z == 0 {
        return None;
    }

    let valences = element::by_atomic_num(z).map_or(&[-1i8][..], |e| e.valences);
    let n_outer = i32::from(element::by_atomic_num(z).map_or(0, |e| e.outer_electrons));
    let chg = i32::from(atom.formal_charge);

    // 该元素有明确的价约束吗?
    let constrained = valences.len() != 1 || valences[0] != -1;
    if !constrained {
        // 过渡金属一类:没有可靠的价信息,不猜。
        // 只要成了键就不赋自由基。
        if mol.degree(idx) > 0 {
            return Some(0);
        }
        // 孤立离子:按外层电子数的奇偶给 0 或 1。
        // 电荷离谱时(算出负数)退回 0,不继续算。
        let n_valence = (n_outer - chg).max(0);
        return Some((n_valence % 2) as u8);
    }

    // 注意:截断而非四舍五入 —— 见模块文档第 1 点
    let accum: f32 = mol
        .neighbors(idx)
        .map(|(_, bi)| mol.bonds()[bi as usize].valence_contribution_to(idx))
        .sum::<f32>()
        + f32::from(atom.num_explicit_hs);
    let total_valence = (accum + 0.1) as i32;

    let base_count = if matches!(z, 1 | 2) { 2 } else { 8 };

    // 路线一:按满壳层反推还差几个电子(适用于电负性较大的元素)
    let mut n_radicals = base_count - n_outer - total_valence + chg;
    if n_radicals < 0 {
        n_radicals = 0;
        // 该元素可能超价:挑第一个够用的价态
        if valences.len() > 1 {
            for &val in valences {
                let val = i32::from(val);
                if val - total_valence + chg >= 0 {
                    n_radicals = val - total_valence + chg;
                    break;
                }
            }
        }
    }

    // 路线二:按外层电子数反推(适用于电负性较小的元素);非负时取较小者
    let n_radicals2 = n_outer - total_valence - chg;
    if n_radicals2 >= 0 {
        n_radicals = n_radicals.min(n_radicals2);
    }

    u8::try_from(n_radicals.max(0)).ok()
}

#[cfg(test)]
mod tests {
    use omgkit_io::smiles;

    use super::*;
    use crate::{clean_up, kekulize, perceive_rings, update_property_cache};

    /// 跑完整的第 1/3/4/5/6 步前缀
    fn radicals(smi: &str) -> Vec<u8> {
        let mut m = smiles::parse(smi).unwrap_or_else(|e| panic!("{}", e.render()));
        clean_up(&mut m);
        update_property_cache(&mut m).expect("价键校验应通过");
        let _ = perceive_rings(&mut m);
        kekulize(&mut m).expect("应能 kekulize");
        assign_radicals(&mut m);
        m.atoms().iter().map(|a| a.num_radical_electrons).collect()
    }

    /// 没写方括号的原子一律 0 —— 它们靠隐式氢把价填满,不可能剩下未成对电子。
    #[test]
    fn atoms_without_brackets_never_get_radicals() {
        for smi in ["CCO", "c1ccccc1", "CC(=O)O", "C1CCCCC1", "N#Cc1ccccc1"] {
            assert!(
                radicals(smi).iter().all(|&r| r == 0),
                "{smi}: 非方括号原子不应有自由基"
            );
        }
    }

    /// 经典自由基:甲基、亚甲基卡宾、氮氧自由基。
    #[test]
    fn classic_radicals() {
        assert_eq!(radicals("[CH3]"), vec![1], "甲基自由基:一个未成对电子");
        assert_eq!(radicals("[CH2]"), vec![2], "亚甲基:两个");
        assert_eq!(radicals("[CH]"), vec![3]);
        assert_eq!(radicals("[C]"), vec![4]);
        assert_eq!(radicals("[OH]"), vec![1], "羟基自由基");
        assert_eq!(radicals("[NH2]"), vec![1], "氨基自由基");
    }

    /// `[N]1C=CC=C1` 的氮必须带自由基。这个信息在芳香写法 `[n]1cccc1` 下
    /// 根本不存在 —— 正是本步必须排在 kekulize 之后的理由。
    #[test]
    fn pyrrolyl_nitrogen_is_a_radical() {
        let r = radicals("[N]1C=CC=C1");
        assert_eq!(r[0], 1, "N 应带一个自由基电子,实际 {r:?}");
        assert!(r[1..].iter().all(|&x| x == 0), "碳不应带自由基");
    }

    /// 满价的方括号原子不该被判成自由基。
    #[test]
    fn saturated_bracket_atoms_have_none() {
        assert_eq!(radicals("[CH4]"), vec![0]);
        assert_eq!(radicals("[NH3]"), vec![0]);
        assert_eq!(radicals("[OH2]"), vec![0]);
        assert_eq!(radicals("[NH4+]"), vec![0]);
        assert_eq!(radicals("[OH-]"), vec![0]);
    }

    /// 电荷参与运算,而且**不**是先扣掉电荷再查价表 —— 见模块文档第 2 点。
    #[test]
    fn charge_enters_the_formula_not_the_lookup() {
        assert_eq!(radicals("[CH3+]"), vec![0], "甲基正离子:六电子,无未成对");
        assert_eq!(radicals("[CH3-]"), vec![0], "甲基负离子:孤对,非自由基");
        assert_eq!(radicals("[NH2-]"), vec![0]);
        assert_eq!(radicals("[NH3+]"), vec![1], "氨自由基阳离子");
    }

    /// 通配原子身份未知,不能猜。
    #[test]
    fn dummy_atoms_are_left_alone() {
        assert_eq!(radicals("[*]"), vec![0]);
        assert_eq!(radicals("C[*]"), vec![0, 0]);
    }

    /// 无价约束的元素(过渡金属):成键就判 0,孤立离子按外层电子奇偶。
    #[test]
    fn unconstrained_elements() {
        assert_eq!(radicals("[Fe]"), vec![0], "Fe 外层 8 个电子,偶数");
        assert_eq!(radicals("[Cu]"), vec![1], "Cu 外层 11 个电子,奇数");
        assert_eq!(radicals("[Cu]Cl")[0], 0, "一旦成键就不赋自由基");
    }

    /// 幂等:跑两遍与跑一遍结果相同。
    #[test]
    fn is_idempotent() {
        for smi in ["[CH3]", "[N]1C=CC=C1", "CCO", "[Fe]", "[OH]"] {
            let mut m = smiles::parse(smi).unwrap();
            clean_up(&mut m);
            update_property_cache(&mut m).unwrap();
            let _ = perceive_rings(&mut m);
            kekulize(&mut m).unwrap();

            assign_radicals(&mut m);
            let once: Vec<u8> = m.atoms().iter().map(|a| a.num_radical_electrons).collect();
            assign_radicals(&mut m);
            let twice: Vec<u8> = m.atoms().iter().map(|a| a.num_radical_electrons).collect();
            assert_eq!(once, twice, "{smi}: 不幂等");
        }
    }
}
