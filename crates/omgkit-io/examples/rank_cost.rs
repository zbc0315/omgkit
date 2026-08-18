//! 量一遍两种原子秩的代价:`canonical_ranks` vs `classed_ranks`。
//!
//! ```shell
//! cargo run -p omgkit-io --release --example rank_cost -- harness/corpus/large.smi
//! ```
//!
//! # 为什么要有这个例子
//!
//! `classed_ranks`「比 `canonical_ranks` 贵 3.7 倍」这句话被三处文档引用
//! (`omgkit_depict::hydrogens` 的注释、`classed_ranks` 自己、
//! `omgkit-conf` 的方案),而它最初只是某一次改动里随手记下的一个数。
//! 引用一个没人复核的数,与编一个没有区别 —— 这个例子把它变成可重跑的。
//!
//! 口径:只解析,**不净化**。两个函数吃的是同一个分子,比的是相对代价,
//! 净化与否不改变比值的量级(而且省掉一个 dev-dependency)。

use std::time::Instant;

fn main() {
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "harness/corpus/large.smi".to_string());
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| {
        eprintln!("读不了 {path}:{e}");
        std::process::exit(1);
    });

    let mols: Vec<_> = text
        .lines()
        .filter_map(|l| l.split('\t').next())
        .filter(|s| !s.is_empty())
        .filter_map(|s| omgkit_io::smiles::parse(s).ok())
        .collect();
    let atoms: usize = mols.iter().map(omgkit_core::MolBuilder::num_atoms).sum();
    println!("{} 个分子,{atoms} 个原子(未补氢)", mols.len());

    // 各跑两遍,取第二遍 —— 第一遍会把缓存预热,量出来的是冷启动
    let mut cheap = 0.0;
    let mut dear = 0.0;
    for round in 0..2 {
        let t = Instant::now();
        let mut sink = 0u64;
        for m in &mols {
            sink += u64::from(omgkit_io::canon::canonical_ranks(m).len() as u32);
        }
        let a = t.elapsed().as_secs_f64();

        let t = Instant::now();
        for m in &mols {
            sink += u64::from(omgkit_io::canon::classed_ranks(m).len() as u32);
        }
        let b = t.elapsed().as_secs_f64();

        if round == 1 {
            cheap = a;
            dear = b;
        }
        // 用一下结果,免得优化器把整个循环删掉 —— 那会量出一个漂亮的 0
        assert!(sink > 0);
    }

    println!("canonical_ranks  {:>8.1} ms", cheap * 1e3);
    println!("classed_ranks    {:>8.1} ms", dear * 1e3);
    println!("倍数             {:>8.2}×", dear / cheap);
}
