//! **测量,不是判据**:在同一份语料上量我方生成构型的失败率与耗时,
//! 好与 `harness/baseline_rdkit_etkdg.py` 量出来的参照放在一起看。
//!
//! ```shell
//! cargo run -q -p omgkit-conf --release --example bench_conformers -- harness/corpus/large.smi
//! python3 harness/baseline_rdkit_etkdg.py harness/corpus/large.smi
//! ```
//!
//! # 口径必须对齐,否则那个比值是假的
//!
//! 两侧的计时都**只包住生成本身**:参照那边 `MolFromSmiles` 与 `AddHs` 在计时
//! 之外,这边解析、净化、感知顺反、补显式氢也都在计时之外。两侧都是单线程、
//! 同一份 `large.smi`、同一台机器。
//!
//! **仍然不是同一件工作**,引用比值时要连这句一起写:参照跑的是 ETKDGv3
//! (距离几何 + 它自己那套扭转/手性项),这边跑的是"建界 → 三角光滑 → 嵌入 →
//! 破对称 → 全局定向 → 精修"。比的是"从 SMILES 拿到一组三维坐标要多久",
//! 不是同一个算法的实现快慢。
//!
//! # 为什么不做成闸
//!
//! 墙钟会抖。本仓库为此栽过两次(一次把改别的 crate 的提交打红了),复杂度
//! 那几条判据早已换成数工作量。这里是**报数**,不设阈值。
use std::time::Instant;

fn main() {
    let path = std::env::args().nth(1).expect("语料");
    let text = std::fs::read_to_string(&path).expect("读语料");

    let mut times: Vec<f64> = Vec::new();
    let mut worst: Vec<(f64, String)> = Vec::new();
    let (mut parse_fail, mut fail) = (0usize, 0usize);
    let t0 = Instant::now();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let smi = line.split('\t').next().unwrap_or("").trim();
        let Ok(mut mol) = omgkit_io::smiles::parse(smi) else {
            parse_fail += 1;
            continue;
        };
        if omgkit_chem::pipeline::sanitize(&mut mol).is_err() {
            parse_fail += 1;
            continue;
        }
        omgkit_io::stereo::perceive_bond_stereo(&mut mol);
        let ranks = omgkit_io::canon::classed_ranks(&mol);
        omgkit_chem::add_explicit_hs(&mut mol, &ranks);
        let centers = omgkit_conf::chiral::centers(&mol);

        let t = Instant::now();
        let got = omgkit_conf::pipeline::conformer(&mol, &centers);
        let dt = t.elapsed().as_secs_f64();

        if got.is_err() {
            fail += 1;
        }
        times.push(dt);
        worst.push((dt, smi.to_string()));
    }
    let total = t0.elapsed().as_secs_f64();

    assert!(!times.is_empty(), "一个分子都没跑到 —— 语料是空的?");
    times.sort_by(f64::total_cmp);
    let n = times.len();
    #[allow(clippy::cast_precision_loss)]
    let mean_ms = 1000.0 * times.iter().sum::<f64>() / n as f64;
    #[allow(clippy::cast_precision_loss)]
    let fail_pct = 100.0 * fail as f64 / n as f64;
    #[allow(
        clippy::cast_precision_loss,
        clippy::cast_sign_loss,
        clippy::cast_possible_truncation
    )]
    let p99 = times[(0.99 * n as f64) as usize];
    println!("语料:{path}");
    println!(
        "墙钟合计 {total:.1} s;生成 {n} 个;平均 {mean_ms:.1} ms;中位 {:.1} ms",
        1000.0 * times[n / 2]
    );
    println!("p99 {:.0} ms;最慢 {:.2} s", 1000.0 * p99, times[n - 1]);
    println!("解析/净化失败 {parse_fail};生成失败 {fail}({fail_pct:.2}%)");
    worst.sort_by(|a, b| b.0.total_cmp(&a.0));
    println!("最慢的 5 个:");
    for (dt, smi) in worst.iter().take(5) {
        println!("   {dt:.2} s  {}", &smi[..smi.len().min(80)]);
    }
}
