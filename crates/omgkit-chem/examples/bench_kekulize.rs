//! 逐分子测 kekulize 耗时,用于观察规模增长曲线。
use std::time::Instant;

fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("用法: bench_kekulize <corpus.smi>");
    let text = std::fs::read_to_string(path).unwrap();
    let per_mol = std::env::args().nth(2).is_some();

    let mut ready = Vec::new();
    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let smi = line.split_whitespace().next().unwrap();
        let Ok(mut m) = omgkit_io::smiles::parse(smi) else {
            continue;
        };
        omgkit_chem::clean_up(&mut m);
        if omgkit_chem::update_property_cache(&mut m).is_err() {
            continue;
        }
        let _ = omgkit_chem::perceive_rings(&mut m);
        ready.push((smi.to_string(), m));
    }

    if per_mol {
        // kekulize 会就地改分子,所以每次都得先克隆。克隆本身是 O(V+E),
        // 在大分子上并不便宜 —— 单独测一列减掉,否则读到的是"克隆 + kekulize"
        // 的和,会把一个纯线性的开销误当成算法没做好。
        println!(
            "{:>6}  {:>12}  {:>12}  {:>12}",
            "原子数", "克隆", "合计", "kekulize"
        );
        for (_, m) in &ready {
            let reps = 200usize;

            let t = Instant::now();
            for _ in 0..reps {
                std::hint::black_box(m.clone());
            }
            let clone_us = t.elapsed().as_micros() as f64 / reps as f64;

            let t = Instant::now();
            for _ in 0..reps {
                let mut c = m.clone();
                let _ = omgkit_chem::kekulize(&mut c);
            }
            let total_us = t.elapsed().as_micros() as f64 / reps as f64;

            println!(
                "{:>6}  {:>9.1} µs  {:>9.1} µs  {:>9.1} µs",
                m.num_atoms(),
                clone_us,
                total_us,
                (total_us - clone_us).max(0.0)
            );
        }
    } else {
        let t = Instant::now();
        let mut ok = 0;
        for (_, m) in &ready {
            let mut c = m.clone();
            if omgkit_chem::kekulize(&mut c).is_ok() {
                ok += 1
            }
        }
        println!(
            "{} 分子, kekulize 成功 {}, 合计 {:.1?} ({:.1} µs/分子)",
            ready.len(),
            ok,
            t.elapsed(),
            t.elapsed().as_micros() as f64 / ready.len() as f64
        );
    }
}
