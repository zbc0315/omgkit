//! **外部判官:自穿。** 链有没有从环里穿过去、两根键有没有交叉。
//!
//! # 两步,顺序不能反
//!
//! 1. **先拿真实构象校准检测器。** MMFF 优化过的真实分子里不存在自穿,
//!    所以在它们身上跑出来的 `crossings` / `pierces` **必须是 0** ——
//!    不是 0 就说明检测器在乱报,后面量什么都没意义。
//!    顺便量出真实分子里"不共享原子的键对"最近能到多少,
//!    [`CROSS_TOL`](omgkit_conf::threading::CROSS_TOL) 的取值就靠这个数撑着。
//! 2. **再量我们自己嵌出来的坐标。**
//!
//! 反过来做(先量自己的、看着数不大就收工)是自证:检测器要是根本报不出东西,
//! 那个 0 只说明它没在看。
//!
//! # 这条判据守的是哪个决定
//!
//! 精修阶段"力场里放全部 N² 对"这个改动,目的就是防自穿。**目的必须有判据** ——
//! 否则收益证明不了,也防不住后续改动把它弄丢。
//!
//! ```shell
//! cargo run -p omgkit-conf --release --example threading_oracle -- harness/baseline/smoke.bounds.jsonl
//! ```

use omgkit_conf::bounds;
use omgkit_conf::embed::{embed, reference_distances};
use omgkit_conf::smooth::triangle_smooth;
use omgkit_conf::threading::{self, CROSS_TOL};
use omgkit_core::{BondOrder, MolBuilder};

/// 真实构象上允许报出多少次自穿。**必须是 0** —— 真实分子里没有自穿,
/// 报出来就是检测器在乱报,那样后面量我们自己的坐标毫无意义。
const MAX_FALSE_POSITIVE: u64 = 0;

/// 我们嵌出来的坐标里,**环穿刺**的分子占比上限。
///
/// 环穿刺是"链穿过环",没有中间状态 —— 这一档是硬伤,精修很难救回来
/// (要把链抽出来得先把它塞回去,中途穿过环面,能量墙很高)。
const MAX_PIERCE_FRAC: f64 = 0.05;

/// 但比例闸在小样本上是**触发即抖**的:冒烟档 27 个分子里出 1 个就是 3.7%,
/// 出 2 个就 7.4% 直接红,而全量档 400 个里同样是那 1 个分子(0.2%)。
/// 所以给一个下限:允许的个数是 `max(该数, 比例算出来的)`。
///
/// 这不是放水 —— 全量档上 5% 等于 20 个,这个下限只在小样本上起作用。
const MIN_PIERCE_ALLOWANCE: u64 = 2;

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

fn quantile(v: &[f64], f: f64) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let i = ((v.len() as f64 - 1.0) * f).round() as usize;
    v[i]
}

fn main() {
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "harness/baseline/smoke.bounds.jsonl".into());
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| {
        eprintln!("读不了判官基准 {path}:{e}");
        std::process::exit(1);
    });

    let (mut n_mol, mut n_emb) = (0u64, 0u64);
    // 真实构象上的
    let (mut real_cross, mut real_pierce) = (0u64, 0u64);
    let mut real_gaps: Vec<f64> = Vec::new();
    // 我们嵌出来的
    let (mut our_cross, mut our_pierce, mut our_pierce_mol) = (0u64, 0u64, 0u64);
    let mut our_gaps: Vec<f64> = Vec::new();
    let mut worst: Vec<String> = Vec::new();

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
        n_mol += 1;

        // ---- 一、真实构象:检测器的校准 ----
        if coords.len() == nat {
            let t = threading::detect(&m, &coords);
            real_cross += t.crossings as u64;
            real_pierce += t.pierces as u64;
            if t.pairs > 0 && t.min_gap.is_finite() {
                real_gaps.push(t.min_gap);
            }
        }

        // ---- 二、我们嵌出来的坐标 ----
        let (mut b, _) = bounds::build(&m);
        if triangle_smooth(&mut b).is_err() {
            continue;
        }
        let Ok(e) = embed(&reference_distances(&b), nat) else {
            continue;
        };
        n_emb += 1;
        let t = threading::detect(&m, &e.coords);
        our_cross += t.crossings as u64;
        our_pierce += t.pierces as u64;
        if t.pierces > 0 {
            our_pierce_mol += 1;
            if worst.len() < 6 {
                worst.push(format!("{smi}(穿刺 {})", t.pierces));
            }
        }
        if t.pairs > 0 && t.min_gap.is_finite() {
            our_gaps.push(t.min_gap);
        }
    }

    real_gaps.sort_by(f64::total_cmp);
    our_gaps.sort_by(f64::total_cmp);
    #[allow(clippy::cast_precision_loss)]
    let pierce_frac = our_pierce_mol as f64 / n_emb.max(1) as f64;

    println!("判官:自穿,分子 {n_mol} 个(嵌出来 {n_emb} 个)");
    println!("  ── 一、真实构象(检测器的校准)──");
    println!("    键交叉 {real_cross}、环穿刺 {real_pierce}(都必须是 {MAX_FALSE_POSITIVE})");
    println!(
        "    不共享原子的键对,最近距离:最小 {:.3} / p05 {:.3} / 中位 {:.3} Å",
        real_gaps.first().copied().unwrap_or(f64::NAN),
        quantile(&real_gaps, 0.05),
        quantile(&real_gaps, 0.5)
    );
    println!("    ↑ CROSS_TOL = {CROSS_TOL} 必须低于上面那个**最小值**,否则真实分子会被误报");
    println!("  ── 二、我们嵌出来的坐标 ──");
    println!("    键交叉 {our_cross}、环穿刺 {our_pierce};有穿刺的分子 {our_pierce_mol}/{n_emb}({:.1}%,上限 {:.0}%)",
        100.0 * pierce_frac, 100.0 * MAX_PIERCE_FRAC);
    println!(
        "    最近距离:最小 {:.3} / p05 {:.3} / 中位 {:.3} Å",
        our_gaps.first().copied().unwrap_or(f64::NAN),
        quantile(&our_gaps, 0.05),
        quantile(&our_gaps, 0.5)
    );
    if !worst.is_empty() {
        println!("    穿刺的例子:{}", worst.join("  "));
    }

    let mut fatal = false;
    if n_mol == 0 {
        eprintln!("\n一个分子都没读到");
        fatal = true;
    }
    // **校准这一条最重要**:检测器在真实分子上乱报的话,第二段的数全不作数
    if real_cross > MAX_FALSE_POSITIVE || real_pierce > MAX_FALSE_POSITIVE {
        eprintln!(
            "\n真实构象上报出了自穿(键交叉 {real_cross}、环穿刺 {real_pierce})—— \
             检测器在乱报,或者 CROSS_TOL 定高了。这一条不过,下面的数没有意义"
        );
        fatal = true;
    }
    if let Some(&m) = real_gaps.first() {
        if m <= CROSS_TOL {
            eprintln!(
                "\n真实分子里键对最近能到 {m:.3} Å,而 CROSS_TOL = {CROSS_TOL} —— \
                 阈值没有留在真实分布之下,迟早误报"
            );
            fatal = true;
        }
    }
    #[allow(
        clippy::cast_precision_loss,
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss
    )]
    let allowed = ((n_emb as f64 * MAX_PIERCE_FRAC).ceil() as u64).max(MIN_PIERCE_ALLOWANCE);
    if our_pierce_mol > allowed {
        eprintln!(
            "\n有环穿刺的分子 {our_pierce_mol} 个 > 允许的 {allowed} 个({:.1}%,\
             闸 {:.0}% 且至少 {MIN_PIERCE_ALLOWANCE} 个)—— 链穿过环,精修很难救回来",
            100.0 * pierce_frac,
            100.0 * MAX_PIERCE_FRAC
        );
        fatal = true;
    }
    if fatal {
        std::process::exit(1);
    }
    println!("\n自穿判据全过。");
}
