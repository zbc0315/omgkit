//! 最难的一档:多组分、模板自带立体要求、大环成环、匹配数爆炸。
use std::time::Instant;

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

fn show(name: &str, rxn_smarts: &str, subs: &[&str], cap: usize) {
    println!("\n【{name}】");
    println!("  模板  {rxn_smarts}");
    println!("  底物  {}", subs.join("  +  "));
    let Ok(rxn) = smarts::parse_reaction(rxn_smarts) else {
        println!("  !! 模板解析失败");
        return;
    };
    let inputs: Vec<_> = subs.iter().map(|s| prep(s)).collect();
    let t = Instant::now();
    let sets = run_reactants(&rxn, &inputs, cap, false);
    let dt = t.elapsed();
    println!("  路径  {} 条   耗时 {dt:?}", sets.len());
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
    println!("  不同产物 {} 种", seen.len());
    for p in seen.iter().take(6) {
        println!("    {p}");
    }
    if seen.len() > 6 {
        println!("    ...(另有 {} 种)", seen.len() - 6);
    }
}

fn main() {
    // 一、多组分:Ugi 四组分反应
    show(
        "Ugi 四组分(醛 + 胺 + 羧酸 + 异腈)",
        "[C:1]=[O:2].[N;H2:3][C:4].[C:5](=[O:6])[OH].[C-:7]#[N+:8][C:9]\
         >>[C:1]([N:3][C:4])[C:7](=[O:2])[N:9].[C:5](=[O:6])[OH].[N:8]",
        &["CC=O", "NCc1ccccc1", "CC(=O)O", "[C-]#[N+]C(C)(C)C"],
        50,
    );

    // 二、模板自带立体要求:只对特定构型反应
    show(
        "模板要求 @ 构型(只该命中一半底物)",
        "[C@H:1][OH]>>[C@H:1]Cl",
        &["C[C@H](O)CC"],
        20,
    );
    show(
        "同一模板作用在相反构型上",
        "[C@H:1][OH]>>[C@H:1]Cl",
        &["C[C@@H](O)CC"],
        20,
    );

    // 三、大环成环:双端二胺 + 双端二酸,分子内成大环
    show(
        "大环内酰胺化(分子内,12 元环)",
        "[C:1](=[O:2])[OH].[N;H2:3]>>[C:1](=[O:2])[N:3]",
        &["NCCCCCCCCCCC(=O)O", "NCCCCCCCCCCC(=O)O"],
        20,
    );

    // 四、匹配数爆炸:高度对称的底物 + 宽松模板
    for n in [1usize, 2, 3, 4] {
        let sub = format!("c1ccccc1{}", "Cc1ccccc1".repeat(n));
        show(
            &format!("宽松模板 × {} 个苯环(匹配数爆炸)", n + 1),
            "[cH:1]>>[c:1]F",
            &[&sub],
            10_000,
        );
    }

    // 五、双模板笛卡尔积:两侧各有多个位点
    show(
        "两个反应物各有多个位点(笛卡尔积)",
        "[C:1][OH].[C:2][NH2]>>[C:1][N:2]",
        &["OCCCCO", "NCCCCN"],
        100,
    );
}
