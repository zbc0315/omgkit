//! 净化第 2 步:把超价的"非金属—金属"单键改成配位键。
//!
//! # 为什么需要这一步
//!
//! 配体与金属的配位键按 IUPAC 的建议画成**单键**。这与价键模型冲突:
//! 二茂铁的环戊二烯基碳按单键算是五价,直接被判超价而拒绝净化。
//!
//! 把那条键改成配位键就解决了 —— 配位键的电子对由配体提供,不计入金属侧的
//! 价,也不让配体超价(见 `BondData::valence_contribution_to`)。
//!
//! # 判据
//!
//! | 条件 | 说明 |
//! |---|---|
//! | 原子是**非金属** | 金属定义见 [`is_metal`],是排除法 |
//! | 原子**超价** | 显式价超过"按形式电荷调整后的等效元素"的最大允许价 |
//! | 不是 H / He / F / Ne | 这几个元素不接受配位键 |
//! | 与金属之间是**单键** | 已经是配位键或双键的不动 |
//!
//! 芳香原子有一条特例:显式价**恰好等于**最大允许价、且总连接数为 4 时也算
//! 超价。环戊二烯基–金属体系正需要这条 —— 那里碳的显式价算出来正好是 4。
//!
//! # 改动是有方向的
//!
//! 新键的 `begin` 是**非金属**(给电子的一端),`end` 是金属。端点朝向在配位键
//! 上有语义,不能随手写。这也是净化里唯一会交换键端点的一步。
//!
//! # 多个金属时怎么挑
//!
//! 一个超价原子同时连着**多个**金属时得挑一个来改。判据必须与输入编号无关,
//! 否则同一个分子换个编号会改到不同的键。规范秩是最自然的次级判据,但它属于
//! L3,这一层拿不到,所以这里用局部不变量(元素、电荷、度数)代替。
//!
//! 触发面窄:8839 条语料里这一步只改动 **2 条**分子,两条的超价原子都只连着
//! 一个金属,平局分支一次都没走到。等 L3 的规范秩能被这一层用到时再收紧;
//! 在那之前,不同的排法会挑中不同的金属,但都是合法的画法。

use omgkit_core::{AtomFlags, BondOrder, MolBuilder};

use crate::valence::explicit_valence_nonstrict;

/// 把超价的非金属—金属单键改成配位键(第 2 步)。
///
/// 必须排在价键计算**之前** —— 它的作用正是让那些原子不再超价,
/// 放在后面的话价键计算已经先一步拒绝了整个分子。
///
/// 返回被改动的键数。触发面极窄(8839 条语料里 2 条),调用方若要断言
/// "它确实开了火",拿这个数比零。
pub fn cleanup_organometallics(mol: &mut MolBuilder) -> usize {
    // 先扫一遍看有没有活干。绝大多数分子没有金属,这一趟就直接返回了。
    if !(0..mol.num_atoms() as u32).any(|a| needs_fixing(mol, a)) {
        return 0;
    }

    // 处理顺序影响"多个非金属抢同一个金属"时的结果。用局部不变量排序而不是
    // 原子下标 —— 下标随输入编号而变,同一个分子换个编号就可能得到不同的键。
    let mut order: Vec<u32> = (0..mol.num_atoms() as u32).collect();
    order.sort_by_key(|&a| local_key(mol, a));

    let mut changed = 0;
    for a in order {
        if !needs_fixing(mol, a) {
            continue;
        }
        let Some(bond) = pick_metal_bond(mol, a) else {
            continue;
        };
        if let Some(mut b) = mol.bond_mut(bond) {
            b.set_order(BondOrder::Dative);
        }
        // 配位键的 begin 必须是给电子的一端
        if mol.bonds()[bond as usize].begin != a {
            mol.swap_bond_ends(bond).expect("键下标来自遍历");
        }
        changed += 1;
    }
    changed
}

/// 该原子是不是"超价的非金属,且连着金属单键"。
fn needs_fixing(mol: &MolBuilder, a: u32) -> bool {
    if !is_hypervalent_non_metal(mol, a) || !accepts_dative(mol.atoms()[a as usize].atomic_num) {
        return false;
    }
    mol.neighbors(a).any(|(other, bi)| {
        mol.bonds()[bi as usize].order == BondOrder::Single
            && is_metal(mol.atoms()[other as usize].atomic_num)
    })
}

/// 挑一条要改的键:配位键最少的那个金属;并列时用局部不变量。
fn pick_metal_bond(mol: &MolBuilder, a: u32) -> Option<u32> {
    mol.neighbors(a)
        .filter(|&(other, bi)| {
            mol.bonds()[bi as usize].order == BondOrder::Single
                && is_metal(mol.atoms()[other as usize].atomic_num)
        })
        .min_by_key(|&(other, _)| (dative_count(mol, other), local_key(mol, other)))
        .map(|(_, bi)| bi)
}

fn dative_count(mol: &MolBuilder, a: u32) -> usize {
    mol.neighbors(a)
        .filter(|&(_, bi)| mol.bonds()[bi as usize].order == BondOrder::Dative)
        .count()
}

/// 与输入编号无关的局部排序键。见模块文档里那条已知差异。
fn local_key(mol: &MolBuilder, a: u32) -> (u8, i8, usize, i32) {
    let at = mol.atoms()[a as usize];
    (
        at.atomic_num,
        at.formal_charge,
        mol.degree(a),
        explicit_valence_nonstrict(mol, a),
    )
}

/// 该元素是否接受配位键。氢、氦、氟、氖不接受。
fn accepts_dative(z: u8) -> bool {
    !matches!(z, 1 | 2 | 9 | 10)
}

/// 该元素是否是金属。
///
/// 定义是**排除法** —— 列出不是金属的那些,其余都算。通配原子(0)不算。
#[must_use]
pub fn is_metal(z: u8) -> bool {
    !matches!(
        z,
        0 | 1
            | 2
            | 5
            | 6
            | 7
            | 8
            | 9
            | 10
            | 14
            | 15
            | 16
            | 17
            | 18
            | 33
            | 34
            | 35
            | 36
            | 52
            | 53
            | 54
            | 85
            | 86
    )
}

/// 非金属原子是否超价。
///
/// 判据用的是**按形式电荷调整后**的等效元素:`N+` 按碳算,`O+` 按氮算。
/// 这样 `c1cccc[n+]1-[Fe]` 不会被误改成配位键 —— 毕竟 `c1cccc[n+]1-C`
/// 本来就是合法的。
fn is_hypervalent_non_metal(mol: &MolBuilder, a: u32) -> bool {
    let at = mol.atoms()[a as usize];
    if is_metal(at.atomic_num) {
        return false;
    }
    let eff = i32::from(at.atomic_num) - i32::from(at.formal_charge);
    if eff <= 0 {
        return false;
    }
    let Ok(eff) = u8::try_from(eff) else {
        return false;
    };
    let Some(elem) = omgkit_core::element::by_atomic_num(eff) else {
        return false;
    };
    let Some(&max_v) = elem.valences.last() else {
        return false;
    };
    if max_v <= 0 {
        return false;
    }
    let max_v = i32::from(max_v);
    let ev = explicit_valence_nonstrict(mol, a);
    let total_degree =
        mol.degree(a) + usize::from(at.num_explicit_hs) + usize::from(at.num_implicit_hs);

    // 芳香 + 四连接的特例:环戊二烯基–金属体系里碳的显式价正好等于最大价,
    // 少了这条就修不了二茂铁
    ev > max_v || (ev == max_v && at.flags.contains(AtomFlags::AROMATIC) && total_degree == 4)
}

#[cfg(test)]
mod tests {
    use omgkit_io::smiles;

    use super::*;
    use crate::clean_up;

    fn after(smi: &str) -> MolBuilder {
        let mut m = smiles::parse(smi).unwrap_or_else(|e| panic!("{smi}: {}", e.render()));
        clean_up(&mut m);
        cleanup_organometallics(&mut m);
        m
    }

    /// 二茂铁:那条让整个分子净化不了的键要变成配位键。
    #[test]
    fn ferrocene_carbon_metal_bond_becomes_dative() {
        let smi = "CN(C)C[C-]12C3=C4C5=C1[Fe++]23456789[C-]%10C6=C7C8=C9%10";
        let m = after(smi);
        let dative: Vec<_> = m
            .bonds()
            .iter()
            .filter(|b| b.order == BondOrder::Dative)
            .collect();
        assert_eq!(dative.len(), 1, "应当只有一条键被改成配位键");
        // 给电子的一端必须在 begin
        let d = dative[0];
        assert!(
            !is_metal(m.atoms()[d.begin as usize].atomic_num),
            "begin 应当是给电子的非金属"
        );
        assert!(
            is_metal(m.atoms()[d.end as usize].atomic_num),
            "end 应当是金属"
        );
    }

    /// 没有金属的分子一条键都不该动。
    #[test]
    fn plain_organics_are_untouched() {
        for smi in ["CCO", "c1ccccc1", "CC(=O)O", "N[C@@H](C)C(=O)O"] {
            let mut m = smiles::parse(smi).unwrap();
            clean_up(&mut m);
            let before: Vec<_> = m.bonds().to_vec();
            assert_eq!(cleanup_organometallics(&mut m), 0, "{smi}");
            assert_eq!(m.bonds(), &before[..], "{smi}");
        }
    }

    /// 金属存在但非金属**没有超价**时不动。
    ///
    /// 这条把判据的"超价"那一半钉住 —— 少了它,所有金属—碳键都会被改成
    /// 配位键,而那是错的:格氏试剂的 C—Mg 本来就是单键。
    #[test]
    fn non_hypervalent_metal_bonds_are_untouched() {
        for smi in ["C[Mg]Br", "[Li]CCCC", "C[Pt](C)(C)C"] {
            let m = after(smi);
            assert!(
                m.bonds().iter().all(|b| b.order != BondOrder::Dative),
                "{smi}:碳没有超价,不该改"
            );
        }
    }

    /// 金属的定义是排除法。
    #[test]
    fn metal_predicate() {
        for z in [1, 6, 7, 8, 9, 15, 16, 17, 35, 53] {
            assert!(!is_metal(z), "原子序数 {z} 不该算金属");
        }
        for z in [3, 12, 26, 29, 78, 92] {
            assert!(is_metal(z), "原子序数 {z} 应当算金属");
        }
        assert!(!is_metal(0), "通配原子不算金属");
    }

    /// 幂等:改过之后那条键已经是配位键,第二次不该再动。
    #[test]
    fn is_idempotent() {
        let smi = "CN(C)C[C-]12C3=C4C5=C1[Fe++]23456789[C-]%10C6=C7C8=C9%10";
        let mut m = smiles::parse(smi).unwrap();
        clean_up(&mut m);
        assert_eq!(cleanup_organometallics(&mut m), 1);
        let once: Vec<_> = m.bonds().to_vec();
        assert_eq!(cleanup_organometallics(&mut m), 0, "第二次不该再改");
        assert_eq!(m.bonds(), &once[..]);
    }
}
