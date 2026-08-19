//! **端到端判官:整条流水线跑完之后,构型到底好不好。**
//!
//! 前面几条判官各守一段(界矩阵、特征分解、嵌入、手性、自穿),这一条守的是
//! **产物**:分子进去、坐标出来,那组坐标满不满足化学。
//!
//! 量四件事,**精修前后各量一遍** —— 只报"之后"看不出精修有没有在干活:
//!
//! | 量 | 为什么它是硬的 |
//! |---|---|
//! | 键长/键角越界(按拓扑档分) | 长程对占 87%,只看总数会被稀释 |
//! | 自穿(键交叉 + 环穿刺) | 距离判据看不见它 —— 穿过去时每一对距离都可以合法 |
//! | 手性号 | 立体化学错了分子就是错的,不是"差一点" |
//! | 耗时 | 要与 RDKit 比的那条线 |
//!
//! ```shell
//! cargo run -p omgkit-conf --release --example conformer_oracle -- harness/baseline/smoke.chirality.jsonl
//! ```

use omgkit_conf::chiral::Center;
use omgkit_conf::{bounds, chiral, embed, pipeline, smooth, threading};
use omgkit_core::{BondOrder, ChiralTag, MolBuilder};

/// 精修之后,**键长**(1-2)越界超过这个数的对占比上限。
///
/// 键是最硬的一档:参数表给的区间宽只有 `2×DIST12_TOL = 0.02 Å`,
/// 而精修的误差函数罚**相对**越界,键上最贵 —— 压不下去说明优化器没干活。
const MAX_BOND_VIOL_FRAC: f64 = 0.02;

/// 精修之后仍有环穿刺的分子占比上限(小样本给下限,同 `threading_oracle`)。
const MAX_PIERCE_FRAC: f64 = 0.05;
const MIN_PIERCE_ALLOWANCE: u64 = 2;

/// 精修之后手性号正确的比例下限。
///
/// 嵌入 + 全局反射之后是 86.2%,剩下 13.8% 是个别中心错 ——
/// **那一档正是三维精修该救回来的**,所以这条闸要比 86.2% 高。
const MIN_CHIRAL_OK: f64 = 0.90;

/// 逐档统计越界:`[1-2, 1-3, 1-4, 长程]` 的 `(越界数, 总数)`。
fn viol_by_class(coords: &[[f64; 3]], b: &smooth::Bounds, topo: &[u8]) -> [(u64, u64); 5] {
    let n = b.len();
    let mut out = [(0u64, 0u64); 5];
    for i in 0..n {
        for j in (i + 1)..n {
            let d = ((coords[i][0] - coords[j][0]).powi(2)
                + (coords[i][1] - coords[j][1]).powi(2)
                + (coords[i][2] - coords[j][2]).powi(2))
            .sqrt();
            let over = (b.lower(i, j) - d).max(d - b.upper(i, j)).max(0.0);
            let c = topo[i * n + j] as usize;
            out[c].1 += 1;
            if over > 0.1 {
                out[c].0 += 1;
            }
        }
    }
    out
}

/// 拓扑距离,封顶 4。
fn topo_dist(mol: &MolBuilder, n: usize) -> Vec<u8> {
    let mut topo = vec![4u8; n * n];
    for start in 0..n {
        let mut d = vec![u8::MAX; n];
        d[start] = 0;
        let mut q = std::collections::VecDeque::from([start]);
        while let Some(x) = q.pop_front() {
            if d[x] >= 3 {
                continue;
            }
            let Ok(xu) = u32::try_from(x) else { continue };
            for (y, _) in mol.neighbors(xu) {
                let y = y as usize;
                if y < n && d[y] == u8::MAX {
                    d[y] = d[x] + 1;
                    q.push_back(y);
                }
            }
        }
        for j in 0..n {
            topo[start * n + j] = d[j].min(4);
        }
    }
    topo
}

fn main() {
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "harness/baseline/smoke.chirality.jsonl".into());
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| {
        eprintln!("读不了判官基准 {path}:{e}");
        std::process::exit(1);
    });

    let (mut n_mol, mut n_fail) = (0u64, 0u64);
    let mut before = [(0u64, 0u64); 5];
    let mut after = [(0u64, 0u64); 5];
    let (mut cross_b, mut cross_a) = (0u64, 0u64);
    let (mut pierce_mol_b, mut pierce_mol_a) = (0u64, 0u64);
    let (mut chi_total, mut chi_before, mut chi_after) = (0u64, 0u64, 0u64);
    let (mut iters, mut e_before, mut e_after) = (0u64, 0.0f64, 0.0f64);
    let start = std::time::Instant::now();

    for line in text.lines() {
        let Ok(v) = serde_json::from_str::<serde_json::Value>(line) else {
            continue;
        };
        let nat = v["n"].as_u64().unwrap_or(0) as usize;
        let z: Vec<u8> = v["z"]
            .as_array()
            .map(|a| {
                a.iter()
                    .filter_map(|x| x.as_u64().map(|y| y as u8))
                    .collect()
            })
            .unwrap_or_default();
        let chg: Vec<i8> = v["charge"]
            .as_array()
            .map(|a| {
                a.iter()
                    .filter_map(|x| x.as_i64().map(|y| y as i8))
                    .collect()
            })
            .unwrap_or_default();
        if z.len() != nat {
            continue;
        }
        let mut m = MolBuilder::new();
        for (k, &a) in z.iter().enumerate() {
            let mut ad = omgkit_core::AtomData::new(a);
            ad.formal_charge = chg.get(k).copied().unwrap_or(0);
            m.add_atom_data(ad);
        }
        let mut ok = true;
        for e in v["bonds"].as_array().into_iter().flatten() {
            let Some(t) = e.as_array() else { continue };
            let (Some(i), Some(j), Some(o)) = (
                t.first().and_then(serde_json::Value::as_u64),
                t.get(1).and_then(serde_json::Value::as_u64),
                t.get(2).and_then(serde_json::Value::as_u64),
            ) else {
                continue;
            };
            let ord = match o {
                2 => BondOrder::Double,
                3 => BondOrder::Triple,
                4 => BondOrder::Aromatic,
                _ => BondOrder::Single,
            };
            if m.add_bond(i as u32, j as u32, ord).is_err() {
                ok = false;
            }
        }
        if !ok || m.num_atoms() != nat {
            continue;
        }
        for c in v["centers"].as_array().into_iter().flatten() {
            let (Some(atom), Some(tag)) = (c["atom"].as_u64(), c["tag"].as_i64()) else {
                continue;
            };
            if let Some(a) = m.atom_mut(atom as u32) {
                a.chiral_tag = if tag == 2 {
                    ChiralTag::Ccw
                } else {
                    ChiralTag::Cw
                };
            }
        }
        let centers: Vec<Center> = chiral::centers(&m);

        // ---- 精修前:嵌入 + 全局反射 ----
        let (mut b, _) = bounds::build(&m);
        if smooth::triangle_smooth(&mut b).is_err() {
            n_fail += 1;
            continue;
        }
        let Ok(emb) = embed::embed(&embed::reference_distances(&b), nat) else {
            n_fail += 1;
            continue;
        };
        let mut pre = emb.coords;
        if chiral::needs_reflection(&pre, &centers) {
            chiral::reflect(&mut pre);
        }

        // ---- 精修后 ----
        let Ok(conf) = pipeline::conformer(&m, &centers) else {
            n_fail += 1;
            continue;
        };
        n_mol += 1;
        iters += conf.iterations as u64;
        e_before += conf.energy_before;
        e_after += conf.energy;

        let topo = topo_dist(&m, nat);
        for (k, (bad, tot)) in viol_by_class(&pre, &b, &topo).iter().enumerate() {
            before[k].0 += bad;
            before[k].1 += tot;
        }
        for (k, (bad, tot)) in viol_by_class(&conf.coords, &b, &topo).iter().enumerate() {
            after[k].0 += bad;
            after[k].1 += tot;
        }
        let tb = threading::detect(&m, &pre);
        let ta = threading::detect(&m, &conf.coords);
        cross_b += tb.crossings as u64;
        cross_a += ta.crossings as u64;
        if tb.pierces > 0 {
            pierce_mol_b += 1;
        }
        if ta.pierces > 0 {
            pierce_mol_a += 1;
        }
        chi_total += centers.len() as u64;
        chi_before += chiral::correct_count(&pre, &centers) as u64;
        chi_after += conf.chiral_ok as u64;
    }
    let elapsed = start.elapsed();

    #[allow(clippy::cast_precision_loss)]
    let pct = |a: u64, b: u64| 100.0 * a as f64 / b.max(1) as f64;
    println!("判官:端到端构型,分子 {n_mol} 个(失败 {n_fail})");
    println!(
        "  精修:平均 {} 步;误差 {:.3e} → {:.3e}(降 {:.1}%)",
        iters / n_mol.max(1),
        e_before / n_mol.max(1) as f64,
        e_after / n_mol.max(1) as f64,
        100.0 * (1.0 - e_after / e_before.max(1e-30))
    );
    println!("  ── 越界 >0.1 Å 的对,按拓扑档(精修前 → 后)──");
    for (c, name) in [
        (1usize, "1-2 键"),
        (2, "1-3 角"),
        (3, "1-4 扭转"),
        (4, "长程"),
    ] {
        println!(
            "    {name:8} 对数 {:6}   {:5.1}% → {:5.1}%",
            after[c].1,
            pct(before[c].0, before[c].1),
            pct(after[c].0, after[c].1)
        );
    }
    println!("  ── 自穿 ──");
    println!("    键交叉 {cross_b} → {cross_a};有环穿刺的分子 {pierce_mol_b} → {pierce_mol_a}(共 {n_mol})");
    println!("  ── 手性 ──");
    println!(
        "    中心 {chi_total} 个:{:.1}% → {:.1}%",
        pct(chi_before, chi_total),
        pct(chi_after, chi_total)
    );
    println!(
        "  ── 耗时 ── 合计 {:.2} 秒,每分子 {:.2} ms",
        elapsed.as_secs_f64(),
        1000.0 * elapsed.as_secs_f64() / n_mol.max(1) as f64
    );

    let mut fatal = false;
    if n_mol == 0 {
        eprintln!("\n一个分子都没跑成");
        fatal = true;
    }
    let bond_frac = pct(after[1].0, after[1].1) / 100.0;
    if bond_frac > MAX_BOND_VIOL_FRAC {
        eprintln!(
            "\n精修之后键长越界的对占 {:.1}% > {:.0}% —— 最硬的一档没压下去",
            100.0 * bond_frac,
            100.0 * MAX_BOND_VIOL_FRAC
        );
        fatal = true;
    }
    #[allow(
        clippy::cast_precision_loss,
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss
    )]
    let allowed = ((n_mol as f64 * MAX_PIERCE_FRAC).ceil() as u64).max(MIN_PIERCE_ALLOWANCE);
    if pierce_mol_a > allowed {
        eprintln!("\n精修之后仍有 {pierce_mol_a} 个分子环穿刺 > 允许的 {allowed} 个");
        fatal = true;
    }
    if chi_total > 0 && pct(chi_after, chi_total) / 100.0 < MIN_CHIRAL_OK {
        eprintln!(
            "\n精修之后手性正确率 {:.1}% < {:.0}% —— 个别中心那一档本该由精修救回来",
            pct(chi_after, chi_total),
            100.0 * MIN_CHIRAL_OK
        );
        fatal = true;
    }
    if fatal {
        std::process::exit(1);
    }
    println!("\n端到端判据全过。");
}
