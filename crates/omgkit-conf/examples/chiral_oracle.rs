//! **外部判官:手性中心抽得对不对。**
//!
//! # 为什么这一条非要外部判官
//!
//! 手性抽取里"四个配体按什么槽位排"这个约定,排错的话**错法是整批一致**的 ——
//! 于是"符号正确率"要么 0% 要么 100%,两个数看起来都像"约定定死了",
//! 实际一个全对一个全错。**这种错推不出来,只能拿真值比。**
//!
//! 真值取自 `harness/dump_chirality.py`:每个中心的有符号体积**在真实构象上的
//! 实际符号**。注意真值不是"标记推出来的号"—— 那正是待验的东西,拿它当真值
//! 就成了自证。导出时另外自检过一遍:标记推的号与真实构象算的号处处一致。
//!
//! ```shell
//! .venv/bin/python harness/dump_chirality.py harness/corpus/large.smi harness/baseline/smoke.chirality.jsonl 150
//! cargo run -p omgkit-conf --release --example chiral_oracle
//! ```

use omgkit_conf::bounds;
use omgkit_conf::chiral;
use omgkit_conf::embed::{embed, reference_distances};
use omgkit_conf::smooth::triangle_smooth;
use omgkit_core::{BondOrder, ChiralTag, MolBuilder};

/// 符号预测错的中心数上限。**这一条必须是 0** ——
/// 它不是统计量:每个中心的号要么对要么错,没有"大部分对"这回事。
const MAX_WRONG: u64 = 0;

/// 抽漏的中心数上限。真值里有而我们没抽出来的,同样是 0。
///
/// **单看"符号正确率"是不够的**:一个什么都不抽的实现,正确率是 0/0,
/// 而 0 个错误 —— 单向的闸又会奖励"什么都不做"。
const MAX_MISSED: u64 = 0;

/// 做完一次全局反射之后,手性号正确的中心占比**下限**。
///
/// # 这个数决定了要不要上四维
///
/// 实测(247 个中心):嵌完直接对的 **53.0%** —— 基本就是掷硬币,
/// 与"嵌入给出的坐标系定向是任意的"完全吻合。做一次全局反射之后 **86.2%**。
///
/// 剩下的 13.8% 是**个别中心相对多数错**,全局反射按定义救不了它们。
/// 但那一档**三维精修救得了**:翻一个中心只要它自己的手性体积过零,
/// 是局部的、有限的势垒;而全局反射要求**所有**中心同时压平 ——
/// 两者的势垒不是一个量级。所以四维先不做,等精修落地之后量剩下多少再说。
///
/// 闸设在 0.80,是贴着现值 0.862 的棘轮:全局反射那段逻辑一旦退化,
/// 这个数会掉回 0.53 附近,当场就红。
const MIN_CHIRAL_AFTER_REFLECT: f64 = 0.80;

fn floats3(v: &serde_json::Value) -> Vec<[f64; 3]> {
    v.as_array()
        .map(|a| {
            a.iter()
                .filter_map(|p| {
                    let q = p.as_array()?;
                    Some([
                        q.first()?.as_f64()?,
                        q.get(1)?.as_f64()?,
                        q.get(2)?.as_f64()?,
                    ])
                })
                .collect()
        })
        .unwrap_or_default()
}

fn main() {
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "harness/baseline/smoke.chirality.jsonl".into());
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| {
        eprintln!("读不了判官基准 {path}:{e}");
        eprintln!(
            "先跑:.venv/bin/python harness/dump_chirality.py harness/corpus/large.smi {path} 150"
        );
        std::process::exit(1);
    });

    let (mut n_mol, mut n_build_fail) = (0u64, 0u64);
    let (mut n_truth, mut n_found, mut n_wrong, mut n_missed) = (0u64, 0u64, 0u64, 0u64);
    let (mut n_vol_bad, mut worst_vol) = (0u64, 0.0f64);
    // 嵌入之后手性对不对 —— 这一组数决定"要不要上四维"
    let (mut n_emb_total, mut n_emb_before, mut n_emb_after) = (0u64, 0u64, 0u64);
    let (mut n_reflected, mut n_mol_embedded, mut n_mol_all_right) = (0u64, 0u64, 0u64);
    let mut wrong_cases: Vec<String> = Vec::new();

    for line in text.lines() {
        let Ok(v) = serde_json::from_str::<serde_json::Value>(line) else {
            continue;
        };
        let smi = v["smiles"].as_str().unwrap_or("").to_string();
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
        let coords = floats3(&v["coords"]);
        if z.len() != nat || coords.len() != nat {
            n_build_fail += 1;
            continue;
        }

        // 按导出的连接表建分子 —— 下标天生对齐,而且**四配位齐全、没有隐式氢**,
        // 正好满足 `chiral::centers` 的两条前置条件。
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
        if !ok {
            n_build_fail += 1;
            continue;
        }
        // 立体标记按真值给(RDKit 的 tag:1 = Cw/@@、2 = Ccw/@)
        let truth: Vec<(u32, i64, i64)> = v["centers"]
            .as_array()
            .map(|a| {
                a.iter()
                    .filter_map(|c| {
                        Some((
                            c["atom"].as_u64()? as u32,
                            c["tag"].as_i64()?,
                            c["sign"].as_i64()?,
                        ))
                    })
                    .collect()
            })
            .unwrap_or_default();
        for &(atom, tag, _) in &truth {
            if let Some(a) = m.atom_mut(atom) {
                a.chiral_tag = if tag == 2 {
                    ChiralTag::Ccw
                } else {
                    ChiralTag::Cw
                };
            }
        }
        n_mol += 1;
        n_truth += truth.len() as u64;

        let got = chiral::centers(&m);
        n_found += got.len() as u64;

        // ---- 顺带量一件决定性的事:**嵌入给出的手性对不对,全局反射能救回多少** ----
        //
        // 这个数直接决定"要不要上四维"。RDKit 一有手性中心就嵌到四维,
        // 为的是让手性能连续翻(三维里翻不了,反射不在 SO(3) 连通分支里)。
        // 如果嵌入 + 一次离散的全局反射就能把绝大多数中心弄对,四维就不必做。
        if !got.is_empty() {
            let (mut bm, _) = bounds::build(&m);
            if triangle_smooth(&mut bm).is_ok() {
                if let Ok(e) = embed(&reference_distances(&bm), nat) {
                    let before = chiral::correct_count(&e.coords, &got);
                    let mut c2 = e.coords.clone();
                    if chiral::needs_reflection(&c2, &got) {
                        chiral::reflect(&mut c2);
                        n_reflected += 1;
                    }
                    let after = chiral::correct_count(&c2, &got);
                    n_emb_total += got.len() as u64;
                    n_emb_before += before as u64;
                    n_emb_after += after as u64;
                    if after == got.len() {
                        n_mol_all_right += 1;
                    }
                    n_mol_embedded += 1;
                }
            }
        }
        for &(atom, _, sign) in &truth {
            let Some(c) = got.iter().find(|c| c.atom == atom) else {
                n_missed += 1;
                continue;
            };
            // 我们**预测**的号是 c.sign;真值是真实构象上量出来的 sign
            #[allow(clippy::cast_precision_loss)]
            let want = sign as f64;
            if (c.sign - want).abs() > 1e-9 {
                n_wrong += 1;
                if wrong_cases.len() < 6 {
                    wrong_cases.push(format!("{smi} 原子 {atom}"));
                }
            }
            // 再独立验一遍几何:拿真实坐标算体积,号必须与真值一致
            let vol = chiral::center_volume(&coords, c);
            if vol == 0.0 || vol.signum() != want {
                n_vol_bad += 1;
            }
            worst_vol = worst_vol.max(vol.abs());
        }
    }

    println!("判官:手性中心,分子 {n_mol} 个(建不出来 {n_build_fail})");
    println!("  真值里的中心 {n_truth} 个;我们抽出 {n_found} 个,漏 {n_missed}(上限 {MAX_MISSED})");
    println!("  符号预测错 {n_wrong} 个(上限 {MAX_WRONG})");
    println!("  拿真实坐标复算体积、号对不上的 {n_vol_bad} 个");
    #[allow(clippy::cast_precision_loss)]
    let pct = |a: u64, b: u64| 100.0 * a as f64 / b.max(1) as f64;
    println!("  ── 嵌入之后的手性(决定要不要上四维)──");
    println!(
        "    中心 {n_emb_total} 个:嵌完直接对的 {n_emb_before}({:.1}%),\
         做一次全局反射之后 {n_emb_after}({:.1}%)",
        pct(n_emb_before, n_emb_total),
        pct(n_emb_after, n_emb_total)
    );
    println!(
        "    分子 {n_mol_embedded} 个:翻了 {n_reflected} 个;**全部中心都对**的 {n_mol_all_right}({:.1}%)",
        pct(n_mol_all_right, n_mol_embedded)
    );
    if !wrong_cases.is_empty() {
        println!("  错的例子:{}", wrong_cases.join("  "));
    }

    let mut fatal = false;
    // 全局反射不是可有可无的一步:没有它,手性正确率就是掷硬币(实测 53%)。
    if pct(n_emb_after, n_emb_total) / 100.0 < MIN_CHIRAL_AFTER_REFLECT {
        eprintln!(
            "\n全局反射之后手性正确率只有 {:.1}% < {:.0}% —— 反射那段退化了",
            pct(n_emb_after, n_emb_total),
            100.0 * MIN_CHIRAL_AFTER_REFLECT
        );
        fatal = true;
    }
    if n_truth == 0 {
        eprintln!("\n真值里一个手性中心都没有 —— 基准文件不对");
        fatal = true;
    }
    if n_missed > MAX_MISSED {
        eprintln!("\n漏抽 {n_missed} 个中心 —— 抽不出来的中心,后面没人管它的手性");
        fatal = true;
    }
    if n_wrong > MAX_WRONG {
        eprintln!(
            "\n符号预测错 {n_wrong} 个 —— 槽位约定错了。注意这种错是**整批一致**的,\
             正确率不是 0% 就是 100%"
        );
        fatal = true;
    }
    if n_vol_bad > 0 {
        eprintln!("\n拿真实坐标复算,{n_vol_bad} 个中心的号对不上 —— 几何那一层就错了");
        fatal = true;
    }
    if fatal {
        std::process::exit(1);
    }
    println!("\n手性判据全过。");
}
