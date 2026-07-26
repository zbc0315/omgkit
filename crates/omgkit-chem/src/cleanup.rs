//! 非标准画法修正(净化第 1 步)。
//!
//! 把化学上等价、但价键不合法的画法改写成两性离子形式:
//!
//! ```text
//! CN(=O)=O      →  C[N+](=O)[O-]      硝基
//! C-N=N#N       →  C-N=[N+]=[N-]      叠氮
//! C=P(=O)X      →  C=[P+]([O-])X      磷酰
//! X(=O)(=O)O    →  [X+2]([O-])([O-])O 高卤酸(X = Cl/Br/I)
//! ```
//!
//! # 为什么必须排在价键校验之前
//!
//! 上面这些写法本身就是超价的,严格校验会直接把分子拒掉。修正之后它们才能
//! 通过。触发面很窄:只在这几种非标准画法上动手,别的结构原样放过。
//!
//! # 全程使用非严格价键计算
//!
//! 本步跑在校验之前,分子里本就存在 `N(=O)=O` 这类超价写法,严格模式会把
//! 它们拒掉,而本步的职责恰恰是修它们。故一律用
//! `valence::explicit_valence_nonstrict`。
//!
//! # 三个易错细节
//!
//! 1. **氮**取**第一个**匹配的氧,**磷**取**最后一个**。
//! 2. **卤素只处理 Cl/Br/I,不含 F**,而且不检查度数。
//! 3. 叠氮那一遍只遍历硝基那一遍**标记过**的原子,不是全部氮。

use omgkit_core::{BondOrder, MolBuilder};

use crate::valence::explicit_valence_nonstrict;

/// 对分子就地施加第 1 步的全部修正。
///
/// 本步不会失败 —— 它使用非严格价键计算,遇到看不懂的结构就原样放过。
pub fn clean_up(mol: &mut MolBuilder) {
    nitrogens_cleanup(mol);
    for i in 0..mol.num_atoms() as u32 {
        match mol.atoms()[i as usize].atomic_num {
            15 => phosphorus_cleanup(mol, i),
            17 | 35 | 53 => halogen_cleanup(mol, i),
            _ => {}
        }
    }
}

/// 找出 `atom` 的第一个满足条件的关联键,返回 `(键下标, 邻居下标)`。
///
/// 遍历顺序即邻接索引的顺序,也就是键的插入顺序。
fn find_bond<F>(mol: &MolBuilder, atom: u32, mut pred: F) -> Option<(u32, u32)>
where
    F: FnMut(&omgkit_core::BondData, u32) -> bool,
{
    mol.neighbors(atom).find_map(|(other, bi)| {
        let b = &mol.bonds()[bi as usize];
        pred(b, other).then_some((bi, other))
    })
}

/// 中性五价氮 → 两性离子形式。
fn nitrogens_cleanup(mol: &mut MolBuilder) {
    // 第一遍:硝基 `N(=O)` → `[N+]-[O-]`
    let mut considered: Vec<u32> = Vec::new();
    for i in 0..mol.num_atoms() as u32 {
        let a = mol.atoms()[i as usize];
        if a.atomic_num != 7 || a.formal_charge != 0 {
            continue;
        }
        // 逐个原子现算 —— 前面的修改会影响后面原子的价
        if explicit_valence_nonstrict(mol, i) != 5 {
            continue;
        }
        considered.push(i);

        let hit = find_bond(mol, i, |b, other| {
            let o = mol.atoms()[other as usize];
            o.atomic_num == 8 && o.formal_charge == 0 && b.order == BondOrder::Double
        });
        if let Some((bi, other)) = hit {
            mol.bond_mut(bi)
                .expect("下标来自遍历")
                .set_order(BondOrder::Single);
            mol.atom_mut(i).expect("下标来自遍历").formal_charge = 1;
            mol.atom_mut(other).expect("邻居下标合法").formal_charge = -1;
        }
    }

    // 第二遍:叠氮 `N#N` → `[N+]=[N-]`。只看第一遍标记过的氮。
    for &i in &considered {
        let hit = find_bond(mol, i, |b, other| {
            let n = mol.atoms()[other as usize];
            n.atomic_num == 7 && n.formal_charge == 0 && b.order == BondOrder::Triple
        });
        if let Some((bi, other)) = hit {
            mol.bond_mut(bi)
                .expect("下标来自遍历")
                .set_order(BondOrder::Double);
            mol.atom_mut(i).expect("下标来自遍历").formal_charge = 1;
            mol.atom_mut(other).expect("邻居下标合法").formal_charge = -1;
        }
    }
}

/// `C=P(=O)X` → `C=[P+]([O-])X`。
fn phosphorus_cleanup(mol: &mut MolBuilder, idx: u32) {
    if mol.atoms()[idx as usize].formal_charge != 0 {
        return;
    }
    if explicit_valence_nonstrict(mol, idx) != 5 || mol.degree(idx) != 3 {
        return;
    }

    // 注意:这里刻意**不提前退出**:`dbl_to_o` 取的是**最后一个**匹配的氧
    let mut dbl_to_o: Option<(u32, u32)> = None;
    let mut has_double_to_c_or_n = false;
    for (other, bi) in mol.neighbors(idx) {
        let b = mol.bonds()[bi as usize];
        let nbr = mol.atoms()[other as usize];
        if nbr.atomic_num == 8 && nbr.formal_charge == 0 && b.order == BondOrder::Double {
            dbl_to_o = Some((bi, other));
        } else if (nbr.atomic_num == 6 || nbr.atomic_num == 7)
            && b.order == BondOrder::Double
            && mol.degree(other) >= 2
        {
            has_double_to_c_or_n = true;
        }
    }

    if let (true, Some((bi, o))) = (has_double_to_c_or_n, dbl_to_o) {
        mol.atom_mut(o).expect("邻居下标合法").formal_charge = -1;
        mol.bond_mut(bi)
            .expect("下标来自遍历")
            .set_order(BondOrder::Single);
        mol.atom_mut(idx).expect("下标来自遍历").formal_charge = 1;
    }
}

/// `X(=O)(=O)O` → `[X+2]([O-])([O-])O`,X = Cl / Br / I。
fn halogen_cleanup(mol: &mut MolBuilder, idx: u32) {
    if mol.atoms()[idx as usize].formal_charge != 0 {
        return;
    }
    let ev = explicit_valence_nonstrict(mol, idx);
    if !matches!(ev, 3 | 5 | 7) {
        return;
    }
    // 全部邻居都必须是氧
    let all_o = mol
        .neighbors(idx)
        .all(|(o, _)| mol.atoms()[o as usize].atomic_num == 8);
    if !all_o {
        return;
    }

    let doubles: Vec<(u32, u32)> = mol
        .neighbors(idx)
        .filter(|&(_, bi)| mol.bonds()[bi as usize].order == BondOrder::Double)
        .map(|(other, bi)| (bi, other))
        .collect();

    let charge = i8::try_from(doubles.len()).unwrap_or(i8::MAX);
    for (bi, other) in doubles {
        mol.bond_mut(bi)
            .expect("下标来自遍历")
            .set_order(BondOrder::Single);
        mol.atom_mut(other).expect("邻居下标合法").formal_charge = -1;
    }
    mol.atom_mut(idx).expect("下标来自遍历").formal_charge = charge;
}

#[cfg(test)]
mod tests {
    use omgkit_io::smiles;

    use super::*;
    use crate::valence::update_property_cache;

    /// 返回 (形式电荷, 键级) 便于逐项断言
    fn cleaned(smi: &str) -> (Vec<i8>, Vec<BondOrder>) {
        let mut m = smiles::parse(smi).unwrap_or_else(|e| panic!("{}", e.render()));
        clean_up(&mut m);
        (
            m.atoms().iter().map(|a| a.formal_charge).collect(),
            m.bonds().iter().map(|b| b.order).collect(),
        )
    }

    #[test]
    fn nitro_group_becomes_zwitterion() {
        // CN(=O)=O -> C[N+](=O)[O-]:第一个双键氧变成单键负氧
        let (q, o) = cleaned("CN(=O)=O");
        assert_eq!(q, vec![0, 1, -1, 0], "N 应带正电,第一个 O 带负电");
        assert_eq!(o[1], BondOrder::Single, "第一条 N=O 变单键");
        assert_eq!(o[2], BondOrder::Double, "第二条保持双键");
    }

    #[test]
    fn azide_becomes_zwitterion() {
        // C-N=N#N -> C-N=[N+]=[N-]
        let (q, o) = cleaned("CN=N#N");
        assert_eq!(q, vec![0, 0, 1, -1]);
        assert_eq!(o[2], BondOrder::Double, "N#N 变 N=N");
    }

    #[test]
    fn charged_nitrogen_is_left_alone() {
        // 已经带电的氮不能再动
        let before = cleaned("O=[n+]1occcc1");
        let mut m = smiles::parse("O=[n+]1occcc1").unwrap();
        let q0: Vec<i8> = m.atoms().iter().map(|a| a.formal_charge).collect();
        clean_up(&mut m);
        let q1: Vec<i8> = m.atoms().iter().map(|a| a.formal_charge).collect();
        assert_eq!(q0, q1, "带电氮不应被修改");
        assert_eq!(before.0, q1);
    }

    #[test]
    fn neutral_low_valence_nitrogen_untouched() {
        // 普通三价氮不满足 valence == 5,不动
        let (q, o) = cleaned("CN(C)C");
        assert_eq!(q, vec![0, 0, 0, 0]);
        assert!(o.iter().all(|&x| x == BondOrder::Single));
    }

    #[test]
    fn phosphorus_zwitterion() {
        // CC=P(=O)CC -> CC=[P+]([O-])CC
        let (q, o) = cleaned("CC=P(=O)CC");
        assert_eq!(q, vec![0, 0, 1, -1, 0, 0]);
        assert_eq!(o[1], BondOrder::Double, "C=P 保持双键");
        assert_eq!(o[2], BondOrder::Single, "P=O 变单键");
    }

    /// 双键对象的**度数必须 ≥ 2** —— 这条限制很容易漏掉。
    /// `C=P(=O)CC` 里那个 C 只连着 P(度数 1),故不做修改。
    #[test]
    fn phosphorus_requires_neighbor_degree_at_least_two() {
        let (q, o) = cleaned("C=P(=O)CC");
        assert!(q.iter().all(|&x| x == 0), "度数不足时不应修改电荷");
        assert_eq!(o[0], BondOrder::Double);
        assert_eq!(o[1], BondOrder::Double, "P=O 应保持双键");
    }

    #[test]
    fn phosphorus_without_double_to_c_or_n_untouched() {
        // 磷酸:只有到 O 的双键,没有到 C/N 的双键 → 不动
        let (q, o) = cleaned("OP(=O)(O)O");
        assert!(q.iter().all(|&x| x == 0));
        assert_eq!(o[1], BondOrder::Double);
    }

    #[test]
    fn perchloric_acid() {
        // OCl(=O)(=O)=O -> [Cl+3]([O-])([O-])([O-])O
        let (q, o) = cleaned("OCl(=O)(=O)=O");
        assert_eq!(q[1], 3, "Cl 带 +3");
        assert_eq!(q.iter().filter(|&&x| x == -1).count(), 3, "三个负氧");
        assert_eq!(o.iter().filter(|&&x| x == BondOrder::Double).count(), 0);
    }

    #[test]
    fn fluorine_is_not_a_halogen_here() {
        // 只覆盖 Cl/Br/I,不含 F
        let (q, o) = cleaned("OF(=O)(=O)=O");
        assert!(q.iter().all(|&x| x == 0), "F 不参与卤素修正");
        assert_eq!(o.iter().filter(|&&x| x == BondOrder::Double).count(), 3);
    }

    #[test]
    fn halogen_with_non_oxygen_neighbor_untouched() {
        let (q, _) = cleaned("CCl(=O)=O");
        assert!(q.iter().all(|&x| x == 0), "邻居里有碳则不动");
    }

    /// 本步的实际意义:让原本超价被拒的分子通过价键校验。
    #[test]
    fn cleanup_makes_overvalent_nitro_pass_valence_check() {
        let smi = "CCCCOC(=O)c1c(cccc1N(=O)=O)C(=O)O";

        let mut without = smiles::parse(smi).unwrap();
        assert!(
            update_property_cache(&mut without).is_err(),
            "不做第 1 步时,五价中性氮应当被价键校验拒绝"
        );

        let mut with = smiles::parse(smi).unwrap();
        clean_up(&mut with);
        assert!(
            update_property_cache(&mut with).is_ok(),
            "做过第 1 步后应当通过"
        );
    }

    #[test]
    fn is_idempotent() {
        // 跑两遍与跑一遍结果相同 —— 修正过的结构不该被反复改写
        for smi in ["CN(=O)=O", "CN=N#N", "C=P(=O)CC", "OCl(=O)(=O)=O", "CCO"] {
            let mut once = smiles::parse(smi).unwrap();
            clean_up(&mut once);
            let mut twice = smiles::parse(smi).unwrap();
            clean_up(&mut twice);
            clean_up(&mut twice);
            assert_eq!(once.atoms(), twice.atoms(), "{smi}: 原子不幂等");
            assert_eq!(once.bonds(), twice.bonds(), "{smi}: 键不幂等");
        }
    }
}
