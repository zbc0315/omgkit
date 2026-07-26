//! 量整条管线各阶段的占比 —— 优化之前先知道瓶颈在哪。
//!
//! 2026-07-26 实测(8839 条语料):解析 21%、净化 68%、写出 12%。
//!
//! 留着这个示例是因为它挡下过一次没必要的优化:写出里的方向判定占写出耗时的
//! 10%,看着值得动手,一量才发现那是整条管线的 **1.2%** —— 而净化占 68%。
//! **对着一个阶段的百分比调优,很容易优化掉一个根本不重要的东西。**
use std::time::Instant;

fn main() {
    let text = std::fs::read_to_string("harness/corpus/large.smi").expect("读不到语料");
    let toks: Vec<&str> = text
        .lines()
        .filter_map(|l| l.split_whitespace().next())
        .filter(|t| !t.starts_with('#'))
        .collect();

    let t = Instant::now();
    let mols: Vec<_> = toks
        .iter()
        .filter_map(|t| omgkit_io::smiles::parse(t).ok())
        .collect();
    let parse = t.elapsed().as_secs_f64() * 1e3;

    let t = Instant::now();
    let sane: Vec<_> = mols
        .iter()
        .filter_map(|m| {
            let mut m = m.clone();
            omgkit_chem::sanitize(&mut m).ok().map(|()| m)
        })
        .collect();
    let sanitize = t.elapsed().as_secs_f64() * 1e3;

    let t = Instant::now();
    for m in &sane {
        std::hint::black_box(omgkit_io::smiles::write(m));
    }
    let write = t.elapsed().as_secs_f64() * 1e3;

    let t = Instant::now();
    for m in &sane {
        std::hint::black_box(omgkit_io::canon::canonical_smiles(m));
    }
    let canonical = t.elapsed().as_secs_f64() * 1e3;

    let total = parse + sanitize + write;
    println!("语料 {} 条 / 净化通过 {} 条", mols.len(), sane.len());
    println!(
        "  解析      {parse:8.1}ms  ({:4.1}%)",
        100.0 * parse / total
    );
    println!(
        "  净化      {sanitize:8.1}ms  ({:4.1}%)",
        100.0 * sanitize / total
    );
    println!(
        "  写出      {write:8.1}ms  ({:4.1}%)",
        100.0 * write / total
    );
    println!("  ——");
    println!(
        "  规范写出  {canonical:8.1}ms  (另计,是写出的 {:.1} 倍)",
        canonical / write
    );
}
