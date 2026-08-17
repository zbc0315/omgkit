//! 量标签朝向的分布:横竖之分的判别量落在哪、竖排会有多少。
//!
//! ```shell
//! cargo run -p omgkit-depict --release --example label_dirs -- harness/corpus/large.smi
//! ```
//!
//! # 为什么要有这个程序
//!
//! `render::VERT_SLOPE` 的阈值、`render::DIR_TIE` 的容差,理由全靠"两簇之间隔
//! 着多少个数量级"这类实测数字撑着。数字写进文档就会过期,而**下一个人没法
//! 复核**。所以把测量本身入库,跟 `harness/gen_elements.py` 是同一个路子。
//!
//! # 口径
//!
//! 按 `Depiction::drawn` 给出的**真正画出来的**那个分子算 —— 补过立体氢之后
//! 原子多了,拿原分子算会少算真的画出来的键。每个分子跑 `Style::ALL` 两套规范。

use omgkit_depict::{
    generate,
    label::LabelDir,
    render::{label_at, label_dir, VERT_SLOPE},
    style::Style,
};

fn main() {
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "harness/corpus/large.smi".into());
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("读不了 {path}:{e}"));

    let mut atoms = 0u64; // 全部原子
    let mut labelled = 0u64; // 有标签的
    let mut with_h = 0u64; // 标签里有氢的
    let mut vertical = 0u64; // 判成竖排的(带氢标签里)
    let mut ties = 0u64; // 判别量落在容差以内的
    let mut ties_deg3 = 0u64; // 其中度 ≥ 3 的 —— 求和次序会变,是真破口
    let mut tie_max = 0.0f64; // 平局簇的上沿
    let mut real_min = f64::MAX; // 真实簇的下沿
    let mut blocked_by_degree_one = 0u64; // 度 1 特判挡下的
    let mut by_deg = [0u64; 8];
    let mut mols = 0u64;

    for line in text.lines() {
        let smi = line.split_whitespace().next().unwrap_or("");
        if smi.is_empty() {
            continue;
        }
        let Ok(mut m) = omgkit_io::smiles::parse(smi) else {
            continue;
        };
        if omgkit_chem::pipeline::sanitize(&mut m).is_err() {
            continue;
        }
        omgkit_io::stereo::perceive_bond_stereo(&mut m);
        mols += 1;
        for style in &Style::ALL {
            let d = generate(&m, style);
            let grown = d.drawn(&m);
            let mol = &*grown;
            for a in 0..u32::try_from(mol.num_atoms()).unwrap_or(0) {
                atoms += 1;
                let Some(l) = label_at(mol, a, style, &d.coords) else {
                    continue;
                };
                labelled += 1;
                if !l.plain().contains('H') {
                    continue;
                }
                with_h += 1;
                let deg = mol.degree(a);
                if deg == 0 {
                    continue;
                }
                let here = d.coords[a as usize];
                let (mut sx, mut sy) = (0.0f64, 0.0f64);
                for (n, _) in mol.neighbors(a) {
                    sx += d.coords[n as usize].x - here.x;
                    sy += d.coords[n as usize].y - here.y;
                }
                let c = sy.abs() - VERT_SLOPE * sx.abs();
                if c.abs() <= 1e-9 {
                    ties += 1;
                    tie_max = tie_max.max(c.abs());
                    if deg >= 3 {
                        ties_deg3 += 1;
                    }
                    continue;
                }
                real_min = real_min.min(c.abs());
                if c > 0.0 {
                    if deg == 1 {
                        blocked_by_degree_one += 1;
                    } else {
                        vertical += 1;
                        by_deg[(deg as usize).min(7)] += 1;
                        debug_assert!(label_dir(mol, a, &d.coords).is_vertical());
                    }
                }
            }
        }
    }

    println!("语料 {path}:{mols} 个分子 × {} 套规范", Style::ALL.len());
    println!("  原子(按真正画出来的分子算)  {atoms}");
    println!("  其中有标签                    {labelled}");
    println!("  其中标签里有氢                {with_h}");
    println!();
    println!("横竖判别量 |sy| − tan70·|sx|(只统计带氢标签、度 ≥ 1):");
    println!("  落在 ±1e-9 以内(平局)       {ties}   其中度 ≥ 3 的 {ties_deg3}");
    println!("    平局簇上沿                  {tie_max:.3e}");
    println!("  其余                          {}", with_h - ties);
    println!("    真实簇下沿                  {real_min:.3e}");
    println!();
    println!(
        "  判成竖排(度 ≥ 2)             {vertical}  ({:.1}%)",
        100.0 * vertical as f64 / with_h as f64
    );
    println!("    按度数 2..5                 {:?}", &by_deg[2..6]);
    println!("  度 1 特判挡下的               {blocked_by_degree_one}");
    println!();
    println!("  一个标签朝向的例子:{:?}", LabelDir::South);
}
