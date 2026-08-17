//! `MolBuilder` 邻接索引的正确性验证,跑在**真实解析结果**上。
//!
//! 索引是随建边增量维护的(半边链表),不是缓存。增量维护的东西一旦写错不会
//! 报错,只会让某个原子悄悄少一个邻居 —— 然后价键、成环、kekulize 全部跟着错,
//! 而且错得像化学问题。所以这里对每个分子把索引与**暴力扫描键表**的结果逐项
//! 比对,核心断言有二:
//!
//! 1. 邻居集合一致,不多不少
//! 2. 邻居**顺序** = 键的插入顺序 —— 手性判定依赖这个顺序
//!
//! 单元测试用的是手搭分子,拓扑规整;真实 SMILES 里有环闭合(键在末尾统一
//! 追加,端点顺序还可能颠倒)、断开的片段、孤立原子 —— 这些才是容易出事的地方。

use std::path::PathBuf;

use omgkit_core::MolBuilder;
use omgkit_io::smiles;

fn corpus(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../harness/corpus")
        .join(name)
}

/// 暴力重算:扫全部键,挑出关联 `atom` 的,按键下标升序。
fn brute_force_neighbors(mol: &MolBuilder, atom: u32) -> Vec<(u32, u32)> {
    mol.bonds()
        .iter()
        .enumerate()
        .filter_map(|(bi, b)| b.other_end(atom).map(|other| (other, bi as u32)))
        .collect()
}

/// 返回 (比对分子数, 比对原子数, 出错描述)
fn check_corpus(name: &str) -> (usize, usize, Vec<String>) {
    let path = corpus(name);
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("读不到语料 {}: {e}", path.display()));

    let (mut n_mols, mut n_atoms) = (0usize, 0usize);
    let mut bad: Vec<String> = Vec::new();

    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let Some(smi) = line.split_whitespace().next() else {
            continue;
        };
        // 解析失败的由 L1 差分测试守着,这里只管解析成功的
        let Ok(mol) = smiles::parse(smi) else {
            continue;
        };
        n_mols += 1;

        for a in 0..mol.num_atoms() as u32 {
            n_atoms += 1;
            let expected = brute_force_neighbors(&mol, a);
            let actual: Vec<(u32, u32)> = mol.neighbors(a).collect();
            if expected != actual {
                bad.push(format!(
                    "{smi}:原子 #{a} 邻居不一致\n    暴力={expected:?}\n    索引={actual:?}"
                ));
            }
            if mol.degree(a) != expected.len() {
                bad.push(format!(
                    "{smi}:原子 #{a} 度数不一致 暴力={} 索引={}",
                    expected.len(),
                    mol.degree(a)
                ));
            }
        }

        // 每条键都必须能被 bond_between 从任一端找回来
        for (bi, b) in mol.bonds().iter().enumerate() {
            let bi = bi as u32;
            for (x, y) in [(b.begin, b.end), (b.end, b.begin)] {
                match mol.bond_between(x, y) {
                    Some(found) if found == bi => {}
                    other => bad.push(format!("{smi}:bond_between({x},{y}) = {other:?},应为 {bi}")),
                }
            }
        }

        if !mol.adjacency_index_is_consistent() {
            bad.push(format!("{smi}:自检 adjacency_index_is_consistent 失败"));
        }
    }

    (n_mols, n_atoms, bad)
}

fn report(name: &str, bad: &[String]) -> String {
    let mut out = format!("\n{name}:邻接索引有 {} 处不一致\n", bad.len());
    for m in bad.iter().take(20) {
        out.push_str("  ");
        out.push_str(m);
        out.push('\n');
    }
    if bad.len() > 20 {
        out.push_str(&format!("  ...(另有 {} 处)\n", bad.len() - 20));
    }
    out
}

#[test]
fn adjacency_index_matches_brute_force_on_smoke() {
    let (mols, atoms, bad) = check_corpus("smoke.smi");
    assert!(mols > 0, "一条都没解析成功");
    assert!(bad.is_empty(), "{}", report("冒烟语料", &bad));
    println!("邻接索引自检通过:{mols} 条分子,{atoms} 个原子");
}

#[test]
#[ignore = "大语料较慢;用 cargo test -- --ignored 运行"]
fn adjacency_index_matches_brute_force_on_large() {
    let (mols, atoms, bad) = check_corpus("large.smi");
    assert!(mols > 1000, "语料不完整:只解析成功 {mols} 条");
    assert!(bad.is_empty(), "{}", report("大语料", &bad));
    println!("邻接索引自检通过:{mols} 条分子,{atoms} 个原子");
}

/// 环闭合键在解析末尾统一追加,而且**按环标号排序**而非按出现先后 ——
/// 这是最容易让"邻居顺序 = 插入顺序"这条不变量出问题的地方。
#[test]
fn ring_closure_bonds_keep_insertion_order() {
    for smi in [
        "c1ccccc1",            // 单环
        "C1CC2CCC1CC2",        // 桥环
        "c1ccc2ccccc2c1",      // 萘,两个环共享一条键
        "C%10CC%10",           // 两位环标号
        "C1CC1C2CC2",          // 两个独立环,标号复用先后
        "N[C@@H](C)C(=O)O",    // 带手性 —— 顺序错了手性就错
        "C1CCCCC1.C1CCCCC1",   // 两个片段
        "[Na+].[Cl-]",         // 全是孤立原子
        "c1cc2cc3cccc3cc2cc1", // 多环稠合
    ] {
        let mol = smiles::parse(smi).unwrap_or_else(|e| panic!("{smi}: {}", e.render()));
        for a in 0..mol.num_atoms() as u32 {
            assert_eq!(
                mol.neighbors(a).collect::<Vec<_>>(),
                brute_force_neighbors(&mol, a),
                "{smi}:原子 #{a} 的邻居顺序不是键的插入顺序"
            );
        }
    }
}
