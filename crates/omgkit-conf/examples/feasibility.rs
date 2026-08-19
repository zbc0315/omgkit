//! **全语料的界可行率** —— 一个分子在进嵌入之前就死掉的比例。
//!
//! # 为什么单独立一个
//!
//! `bounds_oracle` 只跑 400 个分子(受 RDKit 导出基准的规模所限),而这个数
//! 是整个项目的**头号指标**:要赢的正是 RDKit 那 0.52% 的失败率。
//! 界矩阵自相矛盾的分子连嵌入都进不去,它直接计入失败。
//!
//! 这个判据不需要任何外部基准,只要一份 SMILES 语料:
//!
//! ```shell
//! cargo run -p omgkit-conf --release --example feasibility -- harness/corpus/large.smi
//! ```

use omgkit_conf::{bounds, smooth::triangle_smooth};

fn main() {
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "harness/corpus/large.smi".into());
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| {
        eprintln!("读不了语料 {path}:{e}");
        std::process::exit(1);
    });

    let (mut n, mut n_parse_fail, mut n_empty, mut n_infeasible) = (0u64, 0u64, 0u64, 0u64);
    let mut empty_cases: Vec<String> = Vec::new();
    let mut infeasible_cases: Vec<String> = Vec::new();

    for line in text.lines() {
        let smi = line.split('\t').next().unwrap_or("").trim();
        if smi.is_empty() {
            continue;
        }
        let Ok(mut mol) = omgkit_io::smiles::parse(smi) else {
            n_parse_fail += 1;
            continue;
        };
        if omgkit_chem::pipeline::sanitize(&mut mol).is_err() {
            n_parse_fail += 1;
            continue;
        }
        // 补氢要给一个与写法无关的秩;这里只关心界可不可行,用恒等秩即可
        let order: Vec<u32> = (0..mol.num_atoms() as u32).collect();
        omgkit_chem::explicit_hs::add_explicit_hs(&mut mol, &order);
        n += 1;
        let (mut b, _) = bounds::build(&mol);
        // 建完界先看有没有区间当场就空 —— 那是**表自相矛盾**,与几何无关
        let nat = b.len();
        let mut empty = false;
        for i in 0..nat {
            for j in (i + 1)..nat {
                if b.lower(i, j) > b.upper(i, j) {
                    empty = true;
                }
            }
        }
        if empty {
            n_empty += 1;
            if empty_cases.len() < 8 {
                empty_cases.push(smi.to_string());
            }
        }
        if triangle_smooth(&mut b).is_err() {
            n_infeasible += 1;
            if infeasible_cases.len() < 8 {
                infeasible_cases.push(smi.to_string());
            }
        }
    }

    #[allow(clippy::cast_precision_loss)]
    let (pe, pi) = (
        100.0 * n_empty as f64 / n.max(1) as f64,
        100.0 * n_infeasible as f64 / n.max(1) as f64,
    );
    println!("语料 {path}:建界 {n} 个分子(解析/净化失败 {n_parse_fail})");
    println!("  建界即空区间   {n_empty}({pe:.2}%)");
    println!("  光滑化判不可行 {n_infeasible}({pi:.2}%)  ← 这些分子连嵌入都进不去");
    println!("  RDKit ETKDG 的失败率是 0.52% —— 这一行必须低于它,否则算法还没开始就输了");
    if !empty_cases.is_empty() {
        println!("  空区间的例子:{}", empty_cases.join("  "));
    }
    if !infeasible_cases.is_empty() {
        println!("  不可行的例子:{}", infeasible_cases.join("  "));
    }
}
