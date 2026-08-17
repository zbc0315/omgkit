//! 显式氢的合并。
//!
//! # 判据:合并之后必须还是同一个分子
//!
//! 每条用例给一对 SMILES:一条把氢写成图里的原子,一条把同一个分子的氢写在
//! 括号里。合并之后两者的**规范 SMILES 必须逐字节相同** —— 规范化对输入编号
//! 不敏感,所以这条判据不受"氢写在第几位"的影响,正好用来验参照系换算。
//!
//! 手性那一档尤其要紧:氢在邻居里的位置一变,标记就相对另一个参照系了。
//! 错了不会报错,只是分子变成镜像 —— 原子数、键集合、连通性全对,纯拓扑
//! 比对永远发现不了。

use omgkit_chem::{is_removable, remove_hs, sanitize};
use omgkit_core::MolBuilder;
use omgkit_io::{canon, smiles};

fn parsed(smi: &str) -> MolBuilder {
    smiles::parse(smi).unwrap_or_else(|e| panic!("{smi}: {}", e.render()))
}

fn canonical(smi: &str) -> String {
    let mut m = parsed(smi);
    sanitize(&mut m).unwrap_or_else(|e| panic!("{smi}: {e}"));
    canon::canonical_smiles(&m).smiles
}

/// 合并显式氢之后的规范 SMILES。
fn canonical_after_removal(smi: &str) -> String {
    let mut m = parsed(smi);
    remove_hs(&mut m);
    sanitize(&mut m).unwrap_or_else(|e| panic!("{smi} 合并氢之后净化失败: {e}"));
    canon::canonical_smiles(&m).smiles
}

/// 最基本的:氢没了,重原子一个不少。
#[test]
fn explicit_hydrogens_are_merged_into_the_count() {
    for (smi, heavy) in [
        ("[H]C", 1),
        ("[H]C([H])([H])[H]", 1),
        ("[H]OC([H])([H])C", 3),
    ] {
        let mut m = parsed(smi);
        let n = remove_hs(&mut m);
        assert!(n > 0, "{smi}: 一个氢都没合并");
        assert_eq!(m.num_atoms(), heavy, "{smi}: 剩下的重原子数不对");
    }
}

/// 合并不改变总氢数 —— 氢只是从图里的节点搬进了计数。
#[test]
fn merging_preserves_the_total_hydrogen_count() {
    for smi in [
        "[H]C",
        "[H]OC([H])([H])C",
        "[H]N([H])C",
        "c1cc[nH]c1",
        "[H]c1ccccc1",
    ] {
        let mut a = parsed(smi);
        sanitize(&mut a).unwrap_or_else(|e| panic!("{smi}: {e}"));
        let before: u32 = a
            .atoms()
            .iter()
            .map(|x| u32::from(x.num_explicit_hs) + u32::from(x.num_implicit_hs))
            .sum::<u32>()
            + a.atoms().iter().filter(|x| x.atomic_num == 1).count() as u32;

        let mut b = parsed(smi);
        remove_hs(&mut b);
        sanitize(&mut b).unwrap_or_else(|e| panic!("{smi} 合并后: {e}"));
        let after: u32 = b
            .atoms()
            .iter()
            .map(|x| u32::from(x.num_explicit_hs) + u32::from(x.num_implicit_hs))
            .sum::<u32>()
            + b.atoms().iter().filter(|x| x.atomic_num == 1).count() as u32;

        assert_eq!(before, after, "{smi}: 合并前后总氢数不同");
    }
}

/// 四面体手性:氢在邻居里的**每一个位置**都要换对参照系。
///
/// 右边那条是同一个分子不写显式氢的形式。四个位置各测一遍 —— 只测一个位置的话,
/// "从不翻"和"总是翻"里总有一个能蒙混过去。
#[test]
fn chirality_survives_the_merge_from_every_neighbour_position() {
    for (explicit, same_molecule) in [
        // 氢在第 0 位(手性原子是串首)
        ("[C@]([H])(N)(O)F", "N[C@@H](O)F"),
        ("[C@@]([H])(N)(O)F", "N[C@H](O)F"),
        ("[H][C@](N)(O)F", "N[C@@H](O)F"),
        ("[H][C@@](N)(O)F", "N[C@H](O)F"),
        // 第 1 位
        ("N[C@]([H])(O)F", "N[C@H](O)F"),
        ("N[C@@]([H])(O)F", "N[C@@H](O)F"),
        ("C[C@]([H])(O)N", "C[C@@H](N)O"),
        // 第 2 位
        ("N[C@](O)([H])F", "N[C@@H](O)F"),
        ("N[C@@](O)([H])F", "N[C@H](O)F"),
        ("C[C@](O)([H])N", "C[C@H](N)O"),
        // 第 3 位
        ("N[C@](O)(F)[H]", "N[C@H](O)F"),
        ("N[C@@](O)(F)[H]", "N[C@@H](O)F"),
        ("C[C@](O)(N)[H]", "C[C@@H](N)O"),
    ] {
        assert_eq!(
            canonical_after_removal(explicit),
            canonical(same_molecule),
            "{explicit} 合并氢之后应当与 {same_molecule} 是同一个分子 —— \
             不同就是参照系没换对,分子成了镜像"
        );
    }
}

/// 合并**不改变**拓扑上等价的那些分子的判定:没有手性的照样对得上。
#[test]
fn achiral_molecules_round_trip_too() {
    for (explicit, same_molecule) in [
        ("[H]OC([H])([H])C", "CCO"),
        ("[H]c1ccccc1", "c1ccccc1"),
        ("[H]N([H])c1ccccc1", "Nc1ccccc1"),
        ("[H]C(=O)O", "O=CO"),
    ] {
        assert_eq!(
            canonical_after_removal(explicit),
            canonical(same_molecule),
            "{explicit} 合并氢之后应当等于 {same_molecule}"
        );
    }
}

/// 这些氢**不能**删 —— 每一条挡的都是一类会丢信息的删除。
#[test]
fn information_carrying_hydrogens_are_kept() {
    for (smi, why) in [
        ("[2H]C", "氘是另一种核素"),
        ("[3H]C", "氚同上"),
        ("[H+]", "质子是独立物种,而且没有邻居可并"),
        ("[H][H]", "氢分子:两个都删就什么都不剩"),
        ("[H-]", "带电荷"),
        ("[CH3:1][H:2]", "带映射号,反应模板按号引用"),
    ] {
        let mut m = parsed(smi);
        let before = m.num_atoms();
        let n = remove_hs(&mut m);
        assert_eq!(n, 0, "{smi}: 不该删({why}),实际删了 {n} 个");
        assert_eq!(m.num_atoms(), before, "{smi}: 原子数不该变");
    }
}

/// 承载双键方向的那根键上的氢要留着 —— 删了顺反就没了。
#[test]
fn hydrogens_carrying_a_bond_direction_are_kept() {
    let mut m = parsed("[H]/C=C/[H]");
    let n = remove_hs(&mut m);
    assert_eq!(n, 0, "带方向键的氢不该删,否则顺反信息丢失");
    // 反面:同一个骨架但氢不带方向,就该删
    let mut m2 = parsed("[H]C([H])=C([H])[H]");
    assert!(remove_hs(&mut m2) > 0, "不带方向的氢应当照常合并");
}

/// 桥氢与孤立氢不动。
#[test]
fn bridging_and_lone_hydrogens_are_kept() {
    // 硼烷型桥氢:一个氢连着两个硼
    let mut bridge = MolBuilder::new();
    let b1 = bridge.add_atom(5);
    let b2 = bridge.add_atom(5);
    let h = bridge.add_atom(1);
    bridge
        .add_bond(h, b1, omgkit_core::BondOrder::Single)
        .expect("建键");
    bridge
        .add_bond(h, b2, omgkit_core::BondOrder::Single)
        .expect("建键");
    assert!(!is_removable(&bridge, h), "桥氢并给谁都是猜,不该删");

    let mut lone = MolBuilder::new();
    let only = lone.add_atom(1);
    assert!(!is_removable(&lone, only), "孤立的氢并不进谁");
}

/// 幂等:合并过一次之后再合并不再有东西可动。
#[test]
fn merging_is_idempotent() {
    for smi in ["[H]OC([H])([H])C", "N[C@]([H])(O)F", "[H]c1ccccc1"] {
        let mut m = parsed(smi);
        assert!(remove_hs(&mut m) > 0, "{smi}: 第一次应当有得删");
        let before = canon::canonical_smiles(&m).smiles;
        assert_eq!(remove_hs(&mut m), 0, "{smi}: 第二次不该再删出东西");
        assert_eq!(
            canon::canonical_smiles(&m).smiles,
            before,
            "{smi}: 第二次调用改变了分子"
        );
    }
}
