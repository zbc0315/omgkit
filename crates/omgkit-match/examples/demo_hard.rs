//! 难例:对称多位点、成环/断环、递归 SMARTS 选择性、立体依赖。
//!
//! 除了产物本身,还打印**产物集个数** —— 底物上有几条反应路径就有几组,
//! 内容可能重复。个数错了是很容易漏掉的一类缺陷:产物看着都对,只是多了
//! 或少了一条路径。
use omgkit_chem::sanitize;
use omgkit_core::MolBuilder;
use omgkit_io::{canon, smarts, smiles, stereo};
use omgkit_match::{run_reactants, MolProps};

fn prep(smi: &str) -> (MolBuilder, MolProps) {
    let mut m = smiles::parse(smi).unwrap_or_else(|e| panic!("{smi}: {}", e.render()));
    sanitize(&mut m).unwrap_or_else(|e| panic!("{smi}: {e}"));
    stereo::perceive_bond_stereo(&mut m);
    let p = MolProps::compute(&m);
    (m, p)
}

fn main() {
    let cases: &[(&str, &str, &[&str])] = &[
        (
            "对称二元醇:两个羟基等价,应有两条路径",
            "[C:1][OH]>>[C:1]Cl",
            &["OCCCCO"],
        ),
        (
            "对二甲苯:四个等价芳环位点",
            "[cH:1]>>[c:1]Br",
            &["Cc1ccc(C)cc1"],
        ),
        (
            "成环:分子内酯化(羟基酸 → 内酯)",
            "[C:1](=[O:2])[OH].[OH][C:3]>>[C:1](=[O:2])O[C:3]",
            &["OC(=O)CCCCO", "OC(=O)CCCCO"],
        ),
        (
            "断环:哌啶开环(遍历不能绕回另一个模板原子)",
            "[C:1][N:2]>>[C:1].[N:2]",
            &["C1CCNCC1"],
        ),
        (
            "递归 SMARTS 选择性:只酰化非酰胺氮",
            "[NX3;H2;!$(N-C=O):1]>>[N:1]C(C)=O",
            &["NCCNC(C)=O"],
        ),
        (
            "递归 SMARTS 选择性:只氧化苄位",
            "[CH2;$(C-c1ccccc1):1]>>[C:1]=O",
            &["c1ccccc1CCC"],
        ),
        (
            "内消旋酒石酸:两个手性中心相关",
            "[C:1][OH]>>[C:1]Cl",
            &["O[C@@H](C(=O)O)[C@H](O)C(=O)O"],
        ),
        (
            "相互依赖的手性中心(1,4-二取代环己烷)",
            "[C:1][OH]>>[C:1]Cl",
            &["O[C@H]1CC[C@@H](N)CC1"],
        ),
        (
            "共轭多烯:只动一端的羧基,其余双键立体要全保住",
            "[C:1](=[O:2])[OH]>>[C:1](=[O:2])Cl",
            &["C/C=C/C=C/C(=O)O"],
        ),
        (
            "芳香环去芳构化(苯 → 环己二烯)",
            "[c:1]1[c:2][c:3][c:4][c:5][c:6]1>>[C:1]1[C:2]=[C:3][C:4][C:5]=[C:6]1",
            &["c1ccccc1"],
        ),
        (
            "杂环上的电荷改变(吡咯型氮不该丢氢)",
            "[n:1]>>[n+:1]",
            &["c1cc[nH]c1"],
        ),
        (
            "取代基被替换,手性中心度数不变",
            "[C:1][NH2]>>[C:1]O",
            &["N[C@@H](C)C(=O)O"],
        ),
    ];

    for (name, rxn_smarts, subs) in cases {
        println!("\n【{name}】");
        println!("  模板  {rxn_smarts}");
        println!("  底物  {}", subs.join("  +  "));
        let Ok(rxn) = smarts::parse_reaction(rxn_smarts) else {
            println!("  !! 模板解析失败");
            continue;
        };
        let inputs: Vec<_> = subs.iter().map(|s| prep(s)).collect();
        let sets = run_reactants(&rxn, &inputs, 50, false);
        println!("  路径  {} 条", sets.len());
        let mut seen: Vec<String> = Vec::new();
        for outcome in sets {
            let row: Vec<String> = outcome
                .products
                .into_iter()
                .map(|mut p| match sanitize(&mut p) {
                    Ok(()) => canon::canonical_smiles(&p).smiles,
                    Err(e) => format!("<净化失败: {e}>"),
                })
                .collect();
            let j = row.join(" + ");
            if !seen.contains(&j) {
                seen.push(j);
            }
        }
        if seen.is_empty() {
            println!("  产物  (无)");
        }
        for p in &seen {
            println!("  产物  {p}");
        }
    }
}
