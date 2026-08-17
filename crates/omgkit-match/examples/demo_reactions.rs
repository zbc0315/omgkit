//! 用有名字的实际有机反应跑一遍产物生成。
//!
//! 输出 `反应名 / 底物 / 产物`,产物已净化并写成规范 SMILES。
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
            "酰胺缩合(羧酸 + 胺)",
            "[C:1](=[O:2])[OH].[N;H2,H1;!$(N-C=O):3]>>[C:1](=[O:2])[N:3]",
            &["CC(=O)O", "NCc1ccccc1"],
        ),
        (
            "费歇尔酯化(羧酸 + 醇)",
            "[C:1](=[O:2])[OH].[OX2H1:3][C:4]>>[C:1](=[O:2])[O:3][C:4]",
            &["OC(=O)c1ccccc1", "CCO"],
        ),
        (
            "SN2(卤代烃 + 醇盐)",
            "[C:1][Br].[O-:2][C:3]>>[C:1][O:2][C:3]",
            &["CCCCBr", "[O-]CC"],
        ),
        (
            "还原胺化(醛 + 仲胺)",
            "[C:1]=[O:2].[N;H1:3]>>[C:1][N:3]",
            &["CCC=O", "C1CCNCC1"],
        ),
        (
            "铃木偶联(芳基溴 + 硼酸)",
            "[c:1][Br].[c:2][B]([OH])[OH]>>[c:1][c:2]",
            &["Brc1ccc(C)cc1", "OB(O)c1ccccc1"],
        ),
        (
            "狄尔斯–阿尔德(丁二烯 + 亲双烯体)",
            "[C:1]=[C:2][C:3]=[C:4].[C:5]=[C:6]>>[C:1]1[C:2]=[C:3][C:4][C:6][C:5]1",
            &["C=CC=C", "C=CC(=O)OC"],
        ),
        (
            "酯水解",
            "[C:1](=[O:2])[O:3][C:4]>>[C:1](=[O:2])[OH].[OH][C:4]",
            &["CC(=O)OCC"],
        ),
        (
            "醇氧化成酮",
            "[C:1][C:2]([OH])[C:3]>>[C:1][C:2](=O)[C:3]",
            &["CC(O)C"],
        ),
        (
            "保留立体的取代(手性中心不该翻)",
            "[C:1][OH]>>[C:1]Cl",
            &["N[C@@H](C)CO"],
        ),
        (
            "保留双键立体",
            "[C:1](=[O:2])[OH]>>[C:1](=[O:2])OC",
            &["C/C=C/C(=O)O"],
        ),
    ];

    for (name, rxn_smarts, subs) in cases {
        println!("\n【{name}】");
        println!("  模板  {rxn_smarts}");
        println!("  底物  {}", subs.join("  +  "));
        let rxn = match smarts::parse_reaction(rxn_smarts) {
            Ok(r) => r,
            Err(e) => {
                println!("  !! 模板解析失败:{}", e.render());
                continue;
            }
        };
        let inputs: Vec<_> = subs.iter().map(|s| prep(s)).collect();
        let sets = run_reactants(&rxn, &inputs, 20, false);
        if sets.is_empty() {
            println!("  产物  (无 —— 底物匹配不上模板)");
            continue;
        }
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
            let joined = row.join(" + ");
            if !seen.contains(&joined) {
                seen.push(joined);
            }
        }
        for (i, p) in seen.iter().enumerate() {
            println!(
                "  产物{}  {p}",
                if seen.len() > 1 {
                    format!("{}", i + 1)
                } else {
                    String::from(" ")
                }
            );
        }
    }

    demo_atom_mapping();
}

/// 开启原子映射号之后,两侧写出来就是一条完整的映射反应。
///
/// 号是运行时发的:模板只覆盖分子的一小块,其余原子的对应关系来自搬运,
/// 模板里的号在这里派不上用场。
fn demo_atom_mapping() {
    println!("\n===== 原子映射号 =====");
    let cases: &[(&str, &str, &[&str])] = &[
        ("羟基变氯", "[C:1][OH:2]>>[C:1][Cl:2]", &["CC(C)CO"]),
        (
            "酰胺缩合",
            "[C:1](=[O:2])[OH].[N;H2:3]>>[C:1](=[O:2])[N:3]",
            &["CC(=O)O", "NCc1ccccc1"],
        ),
        ("断 C—N(两个产物)", "[C:1][N:2]>>[C:1].[N:2]", &["CCNC"]),
    ];

    for (name, rxn_smarts, subs) in cases {
        println!("\n[{name}]");
        let Ok(rxn) = smarts::parse_reaction(rxn_smarts) else {
            println!("  !! 模板解析失败");
            continue;
        };
        let inputs: Vec<_> = subs.iter().map(|s| prep(s)).collect();
        let outs = run_reactants(&rxn, &inputs, 1, true);
        let Some(o) = outs.into_iter().next() else {
            println!("  (底物匹配不上模板)");
            continue;
        };
        // 不净化不规范化 —— 两者都会重排原子,看号贴在哪就得按存储顺序写
        let side = |ms: &[MolBuilder]| {
            ms.iter()
                .map(|m| smiles::write(m).smiles)
                .collect::<Vec<_>>()
                .join(".")
        };
        println!("  {} >> {}", side(&o.reactants), side(&o.products));
    }
}
