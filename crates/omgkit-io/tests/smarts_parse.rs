//! SMARTS 整体解析的测试。
//!
//! 判据分两层:
//!
//! - **拓扑**:原子数、键数、连接关系 —— 与 SMILES 同一套约定,可以直接比
//! - **语义**:每个原子/键的查询树 —— 逐个比对树的结构
//!
//! 只比拓扑是不够的:`[C]` 和 `[N]` 拓扑完全一样。只比查询树也不够:
//! 环闭合接错了原子,树还是对的。

use omgkit_io::smarts::{self, AtomExpr, AtomPrim, BondExpr, BondPrim};

fn q(s: &str) -> smarts::QueryMol {
    smarts::parse(s).unwrap_or_else(|e| panic!("{s}:\n{}", e.render()))
}

fn elem(z: u8, aromatic: bool) -> AtomExpr {
    AtomExpr::Prim(AtomPrim::Element {
        z,
        aromatic: Some(aromatic),
    })
}

/// 拓扑与 SMILES 同构:分支、环闭合、片段都照搬。
#[test]
fn topology_matches_smiles_shape() {
    for (pat, na, nb) in [
        ("C", 1, 0),
        ("CC", 2, 1),
        ("CCO", 3, 2),
        ("C(=O)O", 3, 2),
        ("c1ccccc1", 6, 6),
        ("C1CC1.CC", 5, 4),
        ("[#6][#7]", 2, 1),
        ("[C,N][O,S]", 2, 1),
    ] {
        let m = q(pat);
        assert_eq!((m.num_atoms(), m.num_bonds()), (na, nb), "{pat}");
        assert!(m.is_consistent(), "{pat}:查询树与拓扑长度不一致");
    }
}

/// 方括号外的裸原子在 SMARTS 里是**查询**,不是具体元素。
#[test]
fn bare_atoms_are_queries() {
    assert_eq!(q("C").atoms[0], elem(6, false), "C 是脂肪碳");
    assert_eq!(q("c").atoms[0], elem(6, true), "c 是芳香碳");
    assert_eq!(q("*").atoms[0], AtomExpr::Prim(AtomPrim::Any));
    assert_eq!(q("a").atoms[0], AtomExpr::Prim(AtomPrim::Aromatic));
    assert_eq!(q("A").atoms[0], AtomExpr::Prim(AtomPrim::Aliphatic));
    assert_eq!(q("Cl").atoms[0], elem(17, false), "两字符有机子集");
}

/// 没写键符号时的默认是"单键**或**芳香键"。
///
/// 这与 SMILES 不同 —— 那里靠两端原子的芳香性当场定死,而 SMARTS 的两端
/// 可能是通配,定不下来,所以默认值本身就是个析取式。
#[test]
fn omitted_bond_is_single_or_aromatic() {
    assert_eq!(q("CC").bonds[0], BondExpr::default_bond());
    assert_eq!(q("**").bonds[0], BondExpr::default_bond());
    assert_eq!(
        q("C=C").bonds[0],
        BondExpr::Prim(BondPrim::Double),
        "写了就按写的来"
    );
}

/// 环闭合键同样延后追加到键表末尾 —— 与 SMILES 同一套约定。
#[test]
fn ring_bonds_are_appended_last() {
    let m = q("C1CC1");
    assert_eq!(m.num_bonds(), 3);
    let b = m.topology.bonds();
    assert_eq!((b[0].begin, b[0].end), (0, 1));
    assert_eq!((b[1].begin, b[1].end), (1, 2));
    assert_eq!(
        (b[2].begin, b[2].end),
        (2, 0),
        "环键在最后,闭合原子在 begin"
    );
}

/// 环闭合上也能写键表达式。
#[test]
fn ring_closure_carries_a_bond_expression() {
    let m = q("C=1CCCCC1");
    assert_eq!(m.bonds[m.num_bonds() - 1], BondExpr::Prim(BondPrim::Double));
    let m = q("C1CCCCC=1");
    assert_eq!(m.bonds[m.num_bonds() - 1], BondExpr::Prim(BondPrim::Double));
    // 两端写了互相矛盾的键表达式要报错
    assert!(smarts::parse("C=1CCCCC#1").is_err());
}

/// `[H]` 特例表:落在表里的是氢元素,其余是氢计数。
#[test]
fn hydrogen_special_case_table() {
    let h_elem = AtomExpr::Prim(AtomPrim::Element {
        z: 1,
        aromatic: None,
    });
    // 在表里 —— 氢元素
    assert_eq!(q("[H]").atoms[0], h_elem, "[H] 是氢元素");
    assert_eq!(
        q("[H+]").atoms[0],
        AtomExpr::And(vec![h_elem.clone(), AtomExpr::Prim(AtomPrim::Charge(1))]),
        "[H+] 在表里"
    );
    assert_eq!(
        q("[2H]").atoms[0],
        AtomExpr::And(vec![AtomExpr::Prim(AtomPrim::Isotope(2)), h_elem.clone()]),
        "[2H] 是氘"
    );
    assert_eq!(
        q("[H:1]").atoms[0],
        AtomExpr::And(vec![h_elem, AtomExpr::Prim(AtomPrim::AtomMap(1))]),
        "[H:1] 在表里"
    );

    // 不在表里 —— 氢计数
    assert_eq!(
        q("[H1]").atoms[0],
        AtomExpr::Prim(AtomPrim::TotalHs(1)),
        "[H1] 的 H 是计数"
    );
    assert_eq!(
        q("[CH]").atoms[0],
        AtomExpr::And(vec![elem(6, false), AtomExpr::Prim(AtomPrim::TotalHs(1))]),
        "[CH] 的 H 是计数"
    );
    assert_eq!(
        q("[HH]").atoms[0],
        AtomExpr::And(vec![
            AtomExpr::Prim(AtomPrim::TotalHs(1)),
            AtomExpr::Prim(AtomPrim::TotalHs(1))
        ]),
        "[HH] 两个都是计数"
    );
    assert_eq!(
        q("[H,C]").atoms[0],
        AtomExpr::Or(vec![AtomExpr::Prim(AtomPrim::TotalHs(1)), elem(6, false)]),
        "[H,C] 的 H 是计数"
    );
}

/// 优先级在整条模式里也要成立 —— 不只是孤立的表达式解析器里。
#[test]
fn precedence_holds_inside_a_pattern() {
    let m = q("[C,N;H1][O,S]");
    assert_eq!(
        m.atoms[0],
        AtomExpr::And(vec![
            AtomExpr::Or(vec![elem(6, false), elem(7, false)]),
            AtomExpr::Prim(AtomPrim::TotalHs(1)),
        ])
    );
    assert_eq!(
        m.atoms[1],
        AtomExpr::Or(vec![elem(8, false), elem(16, false)])
    );
}

/// 键表达式里的逻辑运算。
#[test]
fn bond_expressions_in_a_pattern() {
    assert_eq!(
        q("C-,=C").bonds[0],
        BondExpr::Or(vec![
            BondExpr::Prim(BondPrim::Single),
            BondExpr::Prim(BondPrim::Double)
        ])
    );
    assert_eq!(
        q("C!@C").bonds[0],
        BondExpr::Not(Box::new(BondExpr::Prim(BondPrim::InRing))),
        "非环键"
    );
    assert_eq!(q("C~C").bonds[0], BondExpr::Prim(BondPrim::Any));
}

/// 递归 SMARTS `$(...)`:括号里是一条完整的模式,整体是**一个**原子的条件。
#[test]
fn recursive_smarts() {
    let m = q("[$(CC)]");
    assert_eq!(m.num_atoms(), 1, "递归整体只占一个原子");
    let AtomExpr::Prim(AtomPrim::Recursive(sub)) = &m.atoms[0] else {
        panic!("应当解析成递归基元,实际 {:?}", m.atoms[0]);
    };
    assert_eq!(sub.num_atoms(), 2, "子模式 CC 有两个原子");
    assert_eq!(sub.atoms[0], elem(6, false));

    // 嵌套方括号:外层的 `]` 不能被内层的骗走
    let m = q("[$([CH3])]");
    assert_eq!(m.num_atoms(), 1);
    let AtomExpr::Prim(AtomPrim::Recursive(sub)) = &m.atoms[0] else {
        panic!("应当解析成递归基元");
    };
    assert_eq!(sub.num_atoms(), 1);
    assert_eq!(
        sub.atoms[0],
        AtomExpr::And(vec![elem(6, false), AtomExpr::Prim(AtomPrim::TotalHs(3))])
    );

    // 递归里套递归
    let m = q("[$([#6]-[$([#7])])]");
    assert_eq!(m.num_atoms(), 1);

    // 与其它基元组合,以及出现在或分支里 —— 语料里最常见的形状
    let m = q("[$([#8]),$([#7;!R])]");
    assert_eq!(m.num_atoms(), 1);
    let AtomExpr::Or(parts) = &m.atoms[0] else {
        panic!("应当是或,实际 {:?}", m.atoms[0]);
    };
    assert_eq!(parts.len(), 2);

    // `$` 后面没有括号是语法错误
    assert!(smarts::parse("[$C]").is_err());
    assert!(smarts::parse("[$(CC]").is_err(), "递归括号未闭合");
}

/// 键表达式里的并置就是与。
///
/// 键符号是单字符,但并置不歧义 —— 键表达式夹在两个原子之间,后面没有放
/// 第二条键的位置。`=!@`(双键且非环键)在真实语料里很常见。
#[test]
fn bond_juxtaposition_is_conjunction() {
    assert_eq!(
        q("C=!@C").bonds[0],
        BondExpr::And(vec![
            BondExpr::Prim(BondPrim::Double),
            BondExpr::Not(Box::new(BondExpr::Prim(BondPrim::InRing)))
        ]),
        "双键且非环键"
    );
    assert_eq!(q("C-@C").bonds[0], q("C-&@C").bonds[0]);
}

/// 语法错误要带位置。
#[test]
fn syntax_errors_have_positions() {
    assert!(smarts::parse("").is_err());
    assert!(smarts::parse("C1CC").is_err(), "环未闭合");
    assert!(smarts::parse("C(").is_err(), "括号不匹配");
    assert!(smarts::parse("[C").is_err(), "方括号未闭合");
    assert!(smarts::parse("C=").is_err(), "悬空键");
}
