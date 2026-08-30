//! 对「分子语料 × SMARTS 语料」的每个组合跑一次匹配,输出有命中的那些。
//!
//! 输出格式(每行一条):`分子序号<TAB>模式序号<TAB>命中的原子集合`,
//! 集合形如 `0,1,2|3,4,5`(组内升序,组间按字典序)。
//!
//! 只输出有命中的组合 —— 绝大多数组合是零命中,全写出来的话文件里
//! 99% 是噪声。
//!
//! 用法:`dump_matches <分子.smi> <模式.txt> [分子数上限] [模式数上限]`。

use std::io::{BufWriter, Write};

use omgkit_chem::sanitize;
use omgkit_core::MolBuilder;
use omgkit_io::{smarts, smiles};
use omgkit_match::{substructure_matches, MatchOptions, MolProps};

fn sanitized(smi: &str) -> Option<MolBuilder> {
    let mut m = smiles::parse(smi).ok()?;
    sanitize(&mut m).ok()?;
    Some(m)
}

fn main() {
    let mut args = std::env::args().skip(1);
    let (Some(mols_path), Some(pats_path)) = (args.next(), args.next()) else {
        eprintln!("用法: dump_matches <分子.smi> <模式.txt> [分子数上限] [模式数上限]");
        std::process::exit(2);
    };
    let limit_mols: usize = args
        .next()
        .and_then(|s| s.parse().ok())
        .unwrap_or(usize::MAX);
    let limit_pats: usize = args
        .next()
        .and_then(|s| s.parse().ok())
        .unwrap_or(usize::MAX);

    let mol_text = std::fs::read_to_string(&mols_path).expect("读不到分子语料");
    let pat_text = std::fs::read_to_string(&pats_path).expect("读不到模式语料");

    // 语料里解析失败的条目要**保留占位**,否则序号会错位,
    // 而错位的比对会把全部结果报成分歧
    let mols: Vec<Option<(MolBuilder, MolProps)>> = mol_text
        .lines()
        .filter_map(|l| l.split_whitespace().next())
        .filter(|s| !s.is_empty() && !s.starts_with('#'))
        .take(limit_mols)
        .map(|s| {
            let m = sanitized(s)?;
            let p = MolProps::compute(&m);
            Some((m, p))
        })
        .collect();

    let pats: Vec<Option<smarts::QueryMol>> = pat_text
        .lines()
        .map(str::trim)
        .filter(|s| !s.is_empty() && !s.starts_with('#'))
        .take(limit_pats)
        .map(|s| smarts::parse(s).ok())
        .collect();

    let stdout = std::io::stdout();
    let mut out = BufWriter::new(stdout.lock());
    let opts = MatchOptions {
        max_matches: 1000,
        uniquify: true,
        use_chirality: true,
    };
    let mut pairs = 0usize;
    let mut with_hits = 0usize;

    for (mi, mol) in mols.iter().enumerate() {
        let Some((m, props)) = mol else { continue };
        for (pi, pat) in pats.iter().enumerate() {
            let Some(q) = pat else { continue };
            pairs += 1;
            let hits = substructure_matches(q, m, props, opts);
            if hits.is_empty() {
                continue;
            }
            with_hits += 1;
            let mut sets: Vec<String> = hits
                .iter()
                .map(|h| {
                    let mut v = h.clone();
                    v.sort_unstable();
                    v.iter().map(u32::to_string).collect::<Vec<_>>().join(",")
                })
                .collect();
            sets.sort();
            let _ = writeln!(out, "{mi}\t{pi}\t{}", sets.join("|"));
        }
    }
    let _ = out.flush();
    eprintln!("组合 {pairs} 个,有命中 {with_hits} 个");
}
