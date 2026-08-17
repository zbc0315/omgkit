//! 净化第 12 步:氢的隐式/显式表示调整。
//!
//! # 消失的氢
//!
//! 芳香性感知会让一部分原子的**隐式氢算不出来了**。吡咯的氮是最典型的例子:
//!
//! | 阶段 | 键 | 键级和 | 默认价 | 隐式氢 |
//! |---|---|---|---|---|
//! | kekulize 之后 | 两条单键 | 2 | 3 | **1** |
//! | 芳香化之后 | 两条芳香键 | 1.5 + 1.5 = 3 | 3 | **0** |
//!
//! 氢并没有真的消失 —— 分子还是吡咯,氮上仍有一个氢。丢的是**表示**:
//! 芳香键各算 1.5 之后,价已经被键占满,推断不出还剩多少氢。
//!
//! 本步骤把这一份差额从"隐式"挪到"显式",总氢数不变。写出时也就还原成
//! `[nH]` 而不是 `n` —— 后者是个连 kekulize 都做不了的假分子。
//!
//! # 这一步为什么容易被漏掉
//!
//! 它不改变任何**语义量**:总氢数、价、芳香性全都不变,变的只是同一个氢记在
//! 哪个字段里。所以只比"隐式氢"或只比"显式氢"的测试都发现不了它缺失 ——
//! 缺了这一步,吡咯氮的隐式氢是 0(与参照一致),显式氢也是 0(与参照不一致,
//! 但没人在比)。判据必须是**总氢数**。
//!
//! # 它依赖一个"过期"的值,所以那个值要显式传进来
//!
//! 判据是"重算出来的隐式氢比原来少了多少",而"原来"指的是**芳香化之前**算的
//! 那一份。芳香化之后若先跑一遍价键计算,那份值就被覆盖了,本步骤再去比就是
//! 拿新值比新值,永远相等,什么也不做 —— 而且不报错。
//!
//! 常见做法是让价键计算"惰性"、靠标脏来保住旧值。那样这条依赖是隐式的:
//! 谁在中间多调一次价键计算,这一步就静默失效。所以这里把旧值做成**参数**,
//! 依赖写在签名上。

use omgkit_core::MolBuilder;

use crate::valence::{explicit_valence_nonstrict, implicit_hs_nonstrict};

/// 把因芳香化而推断不出来的隐式氢改记为显式氢(第 12 步)。
///
/// `implicit_before` 是**芳香化之前**每个原子的隐式氢数 —— 见模块文档,
/// 这个依赖必须显式给出。
///
/// 必须排在芳香性感知**之后**。在那之前隐式氢还算得出来,这一步无事可做。
///
/// 返回被调整的原子数。触发面很窄(只有芳香环上带氢的杂原子),所以调用方
/// 若要断言"它确实开了火",拿这个数比零。
///
/// # Panics
/// `implicit_before` 长度与原子数不符时 panic —— 那是调用方的编程错误。
pub fn adjust_hs(mol: &mut MolBuilder, implicit_before: &[u8]) -> usize {
    assert_eq!(
        implicit_before.len(),
        mol.num_atoms(),
        "芳香化前的隐式氢数组长度与原子数不符"
    );
    let mut changed = 0;
    for i in 0..mol.num_atoms() as u32 {
        let orig_implicit = i32::from(implicit_before[i as usize]);
        let ev = explicit_valence_nonstrict(mol, i);
        let new_implicit = i32::from(implicit_hs_nonstrict(mol, i, ev));

        // 只处理"变少了"。变多说明键级减少了,那是别的步骤该管的事,
        // 在这里补显式氢会凭空造出氢来。
        if new_implicit >= orig_implicit {
            continue;
        }
        let gained = orig_implicit - new_implicit;
        if let Some(a) = mol.atom_mut(i) {
            a.num_explicit_hs = a
                .num_explicit_hs
                .saturating_add(u8::try_from(gained).unwrap_or(u8::MAX));
            a.num_implicit_hs = u8::try_from(new_implicit).unwrap_or(0);
        }
        changed += 1;
    }
    changed
}

#[cfg(test)]
mod tests {
    use omgkit_core::AtomFlags;
    use omgkit_io::smiles;

    use super::*;
    use crate::{
        assign_radicals, clean_up, kekulize, perceive_rings, set_aromaticity,
        valence::update_property_cache,
    };

    /// 跑到第 9 步为止(不含本步),并返回芳香化之前的隐式氢快照
    fn upto_step9(smi: &str) -> (MolBuilder, Vec<u8>) {
        let mut m = smiles::parse(smi).unwrap_or_else(|e| panic!("{}", e.render()));
        clean_up(&mut m);
        update_property_cache(&mut m).expect("价键");
        let _ = perceive_rings(&mut m);
        kekulize(&mut m).expect("kekulize");
        assign_radicals(&mut m);
        set_aromaticity(&mut m);
        let before: Vec<u8> = m.atoms().iter().map(|a| a.num_implicit_hs).collect();
        update_property_cache(&mut m).expect("收尾价键");
        (m, before)
    }

    fn total_hs(m: &MolBuilder, i: usize) -> u32 {
        let a = m.atoms()[i];
        u32::from(a.num_explicit_hs) + u32::from(a.num_implicit_hs)
    }

    /// 吡咯的氮:芳香化之后隐式氢算成了 0,本步骤要把它挪回显式。
    #[test]
    fn pyrrole_nitrogen_keeps_its_hydrogen() {
        let (mut m, before) = upto_step9("c1cc[nH]c1");
        assert_eq!(total_hs(&m, 3), 0, "第 9 步之后氢确实'消失'了");

        assert_eq!(adjust_hs(&mut m, &before), 1, "应当只调整了这一个原子");
        assert_eq!(m.atoms()[3].num_explicit_hs, 1);
        assert_eq!(m.atoms()[3].num_implicit_hs, 0);
        assert_eq!(total_hs(&m, 3), 1, "总氢数恢复");
        assert!(
            m.atoms()[3].flags.contains(AtomFlags::AROMATIC),
            "仍是芳香的"
        );
    }

    /// 苯上没有杂原子,本步骤不该动任何东西。
    #[test]
    fn benzene_is_untouched() {
        let (mut m, implicit_before) = upto_step9("c1ccccc1");
        let atoms_before: Vec<_> = m.atoms().to_vec();
        assert_eq!(adjust_hs(&mut m, &implicit_before), 0);
        assert_eq!(m.atoms(), &atoms_before[..]);
    }

    /// N-甲基吡咯的氮没有氢,不该凭空补出一个。
    #[test]
    fn substituted_nitrogen_gains_nothing() {
        let (mut m, before) = upto_step9("Cn1cccc1");
        adjust_hs(&mut m, &before);
        assert_eq!(total_hs(&m, 1), 0, "N 上是甲基,不是氢");
    }

    /// 幂等:再跑一次不该继续加氢。
    ///
    /// 这一步会改写它自己下次要读的字段,写成"每次都补差额"的话,
    /// 反复调用会让氢无限增长。
    #[test]
    fn is_idempotent() {
        for smi in ["c1cc[nH]c1", "c1ccc2[nH]ccc2c1", "c1ccccc1", "CCO"] {
            let (mut m, before) = upto_step9(smi);
            adjust_hs(&mut m, &before);
            let once: Vec<_> = m.atoms().to_vec();
            // 第二次要用**当前**的隐式氢当基准 —— 拿旧快照再比一次等于重放,
            // 那是调用方用错,不是不幂等
            let now: Vec<u8> = m.atoms().iter().map(|a| a.num_implicit_hs).collect();
            assert_eq!(adjust_hs(&mut m, &now), 0, "{smi}:第二次不该再有调整");
            assert_eq!(m.atoms(), &once[..], "{smi}:不幂等");
        }
    }
}
