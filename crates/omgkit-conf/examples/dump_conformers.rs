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
//! **它在 CI 里**(`external` 那个 job 装了 RDKit),`harness/gates.sh` 与
//! `.github/workflows/ci.yml` 都跑它 —— 这段先前写着"进不了 CI(CI 机器上没有
//! RDKit)",那是接进去之前的话。
//!
//! 与它互补的是 `conformer_oracle` 里的"真值口径"(中心与配体序取自基准、
//! 在我们交付的坐标上复算体积):那一条不需要 RDKit,连 `gates` 这个 job 也跑。
//!
//! 实测(`large.smi` 里带立体标记的分子,2026-08-20):
//!
//! | | 覆盖 | 一致 | 判官够不着 |
//! |---|---|---|---|
//! | 只收带 `@` 的、中心基点手性项落地前 | 301 | 288(95.68%) | — |
//! | 同上,落地后 | 301 | 290(96.35%) | — |
//! | 补上双键顺反折算 + 也收 `/` `\` | 632 | 631(99.84%) | — |
//! | 补上三配位立体中心 + 判官自校准 | **642** | **640 / 640(100%)** | **2** |
//!
//! 那 2 个是**三配位磷**:RDKit 的 `AssignStereochemistryFrom3D` 不给它赋手性,
//! 连 RDKit 自己嵌出来的构象都读不回 —— 判官够不着,不是我们摆错。

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
        // 只看带立体标记的分子。
        //
        // **别只收带 `@` 的。** 语料里 331 个分子有 `/` `\` 却没有 `@`,
        // 只收 `@` 的话双键顺反那一档从来进不了外部判据 —— 而那正是这次修的一档。
        // 判据的输入分布排除掉要测的那一档,是这个仓库反复踩的同一个坑。
        if !smi.contains('@') && !smi.contains('/') && !smi.contains('\\') {
            continue;
        }
        let Ok(mut mol) = omgkit_io::smiles::parse(smi) else {
            continue;
        };
        if omgkit_chem::pipeline::sanitize(&mut mol).is_err() {
            continue;
        }
        // 不折算的话双键顺反整档丢掉,见 `pipeline::conformer` 的前置条件那一节
        omgkit_io::stereo::perceive_bond_stereo(&mut mol);
        let r = omgkit_io::canon::classed_ranks(&mol);
        omgkit_chem::add_explicit_hs(&mut mol, &r);
        let centers = omgkit_conf::chiral::centers(&mol);
        // **有手性中心 or 有双键顺反,两者之一就收。** 先前这里只看 `centers`,
        // 于是上面刚放宽的那 331 个"只有 `/` `\`"的分子又在这儿被滤掉了 ——
        // 放宽一道闸而下一道还卡着,等于没放宽。
        let has_ez = mol
            .bonds()
            .iter()
            .any(|b| b.stereo != omgkit_core::BondStereo::None);
        if centers.is_empty() && !has_ez {
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
                //
                // **三配位中心的第四格不是原子下标**,导 `ligands` 全四格的话
                // `4294967295` 会被下游当成原子号:要么越界 panic,
                // 要么在 Python 里静默取错。只导真正落在原子上的那些。
                "centers": centers.iter().map(|c| serde_json::json!({
                    "atom": c.atom,
                    "ligands": c.real_ligands(),
                    "three_coordinate": c.is_three_coordinate(),
                    "sign": c.sign,
                })).collect::<Vec<_>>(),
            })
        );
        n += 1;
    }
    eprintln!("导出 {n} 个分子");
}
