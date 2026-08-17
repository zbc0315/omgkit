//! 对「反应语料 × 分子语料」跑产物生成,输出有产物的那些。
//!
//! 输出 `反应序号<TAB>分子序号<TAB>产物`,产物形如 `CCCl|CCO.CC`
//! (同一组的多个产物用 `.` 连,组间用 `|`,组已排序)。净化失败的产物
//! 记作 `<invalid>`。
//!
//! # 写的是本实现的 SMILES,不是规范形式
//!
//! 判官是外部实现:两边的产物都交给它规范化再比。**不能**拿本实现的规范
//! SMILES 去比外部实现的规范 SMILES —— 那是两套不同的规范化算法,同一个分子
//! 的字符串本来就不一样,比出来的"分歧"全是噪声。
//!
//! # 第一行记的是"问了多少个分子"
//!
//! 判官那边要以两侧的并集为准,才能查出"本实现少了一条"。可它读的是整份分子
//! 语料,而这里可以只跑前 N 个 —— 两个数对不上时,第 N 个之后的分子在判官
//! 眼里全是"只有基准有产物",凭空多出成千上万条假分歧,而且数量大到一眼看去
//! 不像实现问题。
//!
//! 所以覆盖范围写进输出的第一行 `#mols<TAB>N`,判官读不到就直接报错 ——
//! 靠文档提醒对齐守不住,靠数据自己携带才守得住。

use std::io::{BufWriter, Write};

use omgkit_chem::sanitize;
use omgkit_core::MolBuilder;
use omgkit_io::{smarts, smiles};
use omgkit_match::{run_reactants, MolProps};

/// 与 harness/oracle 一致的上限
const MAX_PRODUCT_SETS: usize = 100;

fn read_corpus(path: &str) -> Vec<String> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("读不到 {path}: {e}"))
        .lines()
        .filter_map(|l| l.split_whitespace().next())
        .filter(|t| !t.starts_with('#'))
        .map(String::from)
        .collect()
}

fn main() {
    let mut args = std::env::args().skip(1);
    let (Some(rxn_path), Some(mol_path)) = (args.next(), args.next()) else {
        eprintln!("用法: dump_reactions <反应.txt> <分子.smi> [分子数上限]");
        std::process::exit(2);
    };
    let limit: usize = args
        .next()
        .and_then(|s| s.parse().ok())
        .unwrap_or(usize::MAX);

    let rxn_src = read_corpus(&rxn_path);
    let smis: Vec<String> = read_corpus(&mol_path).into_iter().take(limit).collect();
    let rxns: Vec<Option<smarts::Reaction>> = rxn_src
        .iter()
        .map(|s| smarts::parse_reaction(s).ok())
        .collect();

    let stdout = std::io::stdout();
    let mut out = BufWriter::new(stdout.lock());
    // 覆盖范围必须随数据走,见模块文档
    let _ = writeln!(out, "#mols\t{}", smis.len());
    let (mut pairs, mut with_products) = (0usize, 0usize);

    for (mi, smi) in smis.iter().enumerate() {
        let Ok(mut mol) = smiles::parse(smi) else {
            continue;
        };
        if sanitize(&mut mol).is_err() {
            continue;
        }
        // 先把方向键感知成双键自己的顺反 —— 反应可能正好断掉承载方向的那根键,
        // 那时只有存成双键属性的立体信息才活得下来
        omgkit_io::stereo::perceive_bond_stereo(&mut mol);
        let props = MolProps::compute(&mol);
        let inputs = [(mol, props)];
        for (ri, rxn) in rxns.iter().enumerate() {
            let Some(rxn) = rxn else { continue };
            if rxn.reactants.len() != 1 {
                continue;
            }
            pairs += 1;
            let sets = run_reactants(rxn, &inputs, MAX_PRODUCT_SETS, false);
            if sets.is_empty() {
                continue;
            }
            with_products += 1;
            let mut groups: Vec<String> = sets
                .into_iter()
                .map(|outcome| {
                    outcome
                        .products
                        .into_iter()
                        .map(write_product)
                        .collect::<Vec<_>>()
                        .join(".")
                })
                .collect();
            groups.sort();
            let _ = writeln!(out, "{ri}\t{mi}\t{}", groups.join("|"));
        }
    }
    let _ = out.flush();
    eprintln!("组合 {pairs},有产物 {with_products}");
}

fn write_product(mut p: MolBuilder) -> String {
    match sanitize(&mut p) {
        Ok(()) => smiles::write(&p).smiles,
        Err(_) => "<invalid>".to_string(),
    }
}
