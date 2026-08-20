//! **把我们生成的构型导出成 jsonl**,交给 `harness/verify_stereo.py` 让 RDKit
//! 读回立体化学。
//!
//! 这是**唯一真正外部**的立体化学判据:它不看我们的任何公式,只把坐标交出去,
//! 问"RDKit 从这组坐标读出来的立体化学,与输入 SMILES 指定的一致吗"。
//!
//! ```shell
//! cargo run -p omgkit-conf --release --example dump_conformers -- harness/corpus/large.smi > /tmp/ours.jsonl
//! .venv/bin/python harness/verify_stereo.py /tmp/ours.jsonl
//! ```
//!
//! **进不了 CI**(CI 机器上没有 RDKit),所以它不是闸,是**闸的设计验证**:
//! 能进 CI 的那一条是 `conformer_oracle` 里的"真值口径"——
//! 中心与配体序取自基准、在我们交付的坐标上复算体积。这里跑一遍是为了确认
//! 那一条确实盯住了该盯的东西。
//!
//! 实测(`large.smi` 里 301 个带立体标记的分子,2026-08-20):
//!
//! | | 一致 | 说明 |
//! |---|---|---|
//! | 中心基点手性项落地**前** | 288 / 301(95.68%) | |
//! | 中心基点手性项落地**后** | **290 / 301(96.35%)** | 修好 2 个,弄坏 0 个 |
//!
//! 剩下 11 个:**10 个是环上双键的 E/Z**、1 个是三配位硫(亚磺酰胺)——
//! 两档都是已知的、与四面体手性无关的缺口。四面体手性一个都没错。

fn main() {
    let path = std::env::args().nth(1).expect("语料");
    let limit: usize = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(usize::MAX);
    let text = std::fs::read_to_string(&path).expect("读语料");
    let mut n = 0usize;
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') || n >= limit {
            continue;
        }
        let smi = line.split('\t').next().unwrap_or("").trim();
        // 只看带立体标记的分子
        if !smi.contains('@') {
            continue;
        }
        let Ok(mut mol) = omgkit_io::smiles::parse(smi) else {
            continue;
        };
        if omgkit_chem::pipeline::sanitize(&mut mol).is_err() {
            continue;
        }
        let r = omgkit_io::canon::classed_ranks(&mol);
        omgkit_chem::add_explicit_hs(&mut mol, &r);
        let centers = omgkit_conf::chiral::centers(&mol);
        if centers.is_empty() {
            continue;
        }
        let Ok(conf) = omgkit_conf::pipeline::conformer(&mol, &centers) else {
            continue;
        };
        let z: Vec<u8> = mol.atoms().iter().map(|a| a.atomic_num).collect();
        let chg: Vec<i8> = mol.atoms().iter().map(|a| a.formal_charge).collect();
        let bonds: Vec<[u32; 3]> = mol
            .bonds()
            .iter()
            .map(|b| {
                let o = match b.order {
                    omgkit_core::BondOrder::Double => 2,
                    omgkit_core::BondOrder::Triple => 3,
                    omgkit_core::BondOrder::Aromatic => 4,
                    _ => 1,
                };
                [b.begin, b.end, o]
            })
            .collect();
        println!(
            "{}",
            serde_json::json!({
                "smiles": smi,
                "z": z,
                "charge": chg,
                "bonds": bonds,
                "xyz": conf.coords,
                // 中心与配体序也导出来 —— 下游要能在**我们交付的坐标**上
                // 独立复算体积,而不是只能比 SMILES。
                "centers": centers.iter().map(|c| serde_json::json!({
                    "atom": c.atom, "ligands": c.ligands, "sign": c.sign,
                })).collect::<Vec<_>>(),
            })
        );
        n += 1;
    }
    eprintln!("导出 {n} 个分子");
}
