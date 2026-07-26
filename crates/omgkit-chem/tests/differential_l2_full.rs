//! 整条净化管线的差分测试,基准是**全部 12 步跑完**的状态。
//!
//! # 为什么必须有这一档
//!
//! 逐步骤的差分测试各自停在某一步,谁也看不到最后一步的产出。第 12 步
//! (氢的隐式/显式调整)不改变任何逐步骤测试在比的量,少了它没有任何一档会红。
//!
//! # 判据是**总氢数**,不是隐式氢
//!
//! 吡咯的氮:
//!
//! | | 显式氢 | 隐式氢 | 总氢 |
//! |---|---|---|---|
//! | 基准 | 1 | 0 | **1** |
//! | 缺第 12 步时 | 0 | 0 | **0** |
//!
//! 只比隐式氢的话,两边都是 0,测试全绿而分子已经错了。氢记在哪个字段里是
//! **表示**(kekulize 会来回搬),总氢数才是语义量。所以这里两个都比:
//! 总氢数保证不丢氢,显式/隐式分别比保证表示也对齐。

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use omgkit_chem::sanitize;
use omgkit_core::{AtomFlags, BondFlags};
use omgkit_io::smiles;

/// 原子列下标
const A_EXPL_H: usize = 3;
const A_AROMATIC: usize = 6;
const A_IMPLICIT_H: usize = 8;
const A_EXPLICIT_VALENCE: usize = 9;
const A_HYBRID: usize = 10;
const A_IN_RING: usize = 11;
const A_RADICALS: usize = 13;
/// 键列下标
const B_AROMATIC: usize = 4;
const B_IN_RING: usize = 5;
const B_CONJUGATED: usize = 6;

/// l2 原子行的列数,见 harness/README.md 的列规范。
const A_COLS: usize = 15;

/// 已定位根因的分歧。**不是豁免名单** —— 根因修好之后条目不删,测试同样会红。
///
/// 当前为空。根因修好后测试会主动报错要求删除对应条目,登记表因此不会
/// 沉淀成没人再看的豁免名单。
const KNOWN_DIVERGENCES: &[&str] = &[];

struct Mismatch {
    smi: String,
    field: String,
    baseline: String,
    omgkit: String,
}

#[derive(Default)]
struct Stats {
    n: usize,
    compared: usize,
    /// 至少有一个原子的显式氢非零 —— 第 12 步的产出主要落在这里
    with_explicit_hs: usize,
    /// 命中的已登记分歧数
    known_hit: usize,
}

fn baseline(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../harness/baseline")
        .join(name)
}

fn diff_against(path: &Path) -> (Stats, Vec<Mismatch>) {
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("读不到全管线基准 {}: {e}", path.display()));

    let mut bad: Vec<Mismatch> = Vec::new();
    let mut stats = Stats::default();
    let mut hit_known: std::collections::BTreeSet<&str> = std::collections::BTreeSet::new();

    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let rec: serde_json::Value = serde_json::from_str(line).expect("基准 JSONL 格式错误");
        let smi = rec["smi"].as_str().expect("缺 smi").to_string();
        stats.n += 1;

        let Ok(mut mol) = smiles::parse(&smi) else {
            continue;
        };
        stats.compared += 1;

        let rd_ok = rec["ok"].as_bool().unwrap_or(false);
        let ours = sanitize(&mut mol);

        let known = KNOWN_DIVERGENCES.iter().find(|&&k| k == smi);
        let mut push = |field: String, baseline: String, omgkit: String| {
            if let Some(k) = known {
                hit_known.insert(*k);
                return;
            }
            bad.push(Mismatch {
                smi: smi.clone(),
                field,
                baseline,
                omgkit,
            });
        };

        match (rd_ok, &ours) {
            (false, Ok(())) => {
                push(
                    "净化结果".into(),
                    format!("失败({})", rec["err"].as_str().unwrap_or("?")),
                    "通过".into(),
                );
                continue;
            }
            (true, Err(e)) => {
                push("净化结果".into(), "通过".into(), format!("失败({e})"));
                continue;
            }
            (false, Err(_)) => continue,
            (true, Ok(())) => {}
        }

        if mol.atoms().iter().any(|a| a.num_explicit_hs > 0) {
            stats.with_explicit_hs += 1;
        }

        for (i, row) in rec["atoms"].as_array().unwrap().iter().enumerate() {
            let v: Vec<i64> = row
                .as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_i64().unwrap())
                .collect();
            assert_eq!(
                v.len(),
                A_COLS,
                "{smi}:基准的原子列数是 {},本文件按 {A_COLS} 列解读 —— \
                 基准过期或列号未同步",
                v.len()
            );
            let a = mol.atoms()[i];
            let mut cmp = |name: &str, exp: i64, got: i64| {
                if exp != got {
                    push(
                        format!("原子[{i}].{name}"),
                        exp.to_string(),
                        got.to_string(),
                    );
                }
            };
            // **总氢数**是语义量,必须比。显式/隐式的分配是表示,也比,
            // 但那是为了对齐表示,不是为了守住氢有没有丢。
            cmp(
                "总氢数",
                v[A_EXPL_H] + v[A_IMPLICIT_H],
                i64::from(a.num_explicit_hs) + i64::from(a.num_implicit_hs),
            );
            cmp("显式氢", v[A_EXPL_H], i64::from(a.num_explicit_hs));
            cmp("隐式氢", v[A_IMPLICIT_H], i64::from(a.num_implicit_hs));
            cmp("杂化", v[A_HYBRID], i64::from(a.hybridization as u8));
            cmp(
                "芳香",
                v[A_AROMATIC],
                i64::from(a.flags.contains(AtomFlags::AROMATIC)),
            );
            cmp(
                "在环中",
                v[A_IN_RING],
                i64::from(a.flags.contains(AtomFlags::IN_RING)),
            );
            cmp("自由基", v[A_RADICALS], i64::from(a.num_radical_electrons));
            let _ = A_EXPLICIT_VALENCE; // 显式价由总氢数与键级共同决定,不单独比
        }

        for (i, row) in rec["bonds"].as_array().unwrap().iter().enumerate() {
            let v: Vec<i64> = row
                .as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_i64().unwrap())
                .collect();
            let b = mol.bonds()[i];
            let mut cmp = |name: &str, exp: i64, got: i64| {
                if exp != got {
                    push(format!("键[{i}].{name}"), exp.to_string(), got.to_string());
                }
            };
            cmp(
                "芳香",
                v[B_AROMATIC],
                i64::from(b.flags.contains(BondFlags::AROMATIC)),
            );
            cmp(
                "在环中",
                v[B_IN_RING],
                i64::from(b.flags.contains(BondFlags::IN_RING)),
            );
            cmp(
                "共轭",
                v[B_CONJUGATED],
                i64::from(b.flags.contains(BondFlags::CONJUGATED)),
            );
        }
    }
    // 已登记但**已不再分歧**的条目要报错并要求删除 —— 一个永远躺着没人动的
    // 豁免名单本身就是暗坑
    let stale: Vec<&&str> = KNOWN_DIVERGENCES
        .iter()
        .filter(|k| !hit_known.contains(**k) && text.contains(**k))
        .collect();
    assert!(
        stale.is_empty(),
        "以下已登记的分歧现在不再分歧了,请从 KNOWN_DIVERGENCES 删除:{stale:?}"
    );
    stats.known_hit = hit_known.len();
    (stats, bad)
}

fn report(stats: &Stats, bad: &[Mismatch], limit: usize) -> String {
    let mut by_field: BTreeMap<&str, usize> = BTreeMap::new();
    for m in bad {
        *by_field
            .entry(m.field.split('.').next_back().unwrap_or(&m.field))
            .or_default() += 1;
    }
    let mut by_smi: BTreeMap<&str, Vec<&Mismatch>> = BTreeMap::new();
    for m in bad {
        by_smi.entry(&m.smi).or_default().push(m);
    }
    let mut out = format!(
        "\nL2 全管线差分失败:基准 {} 条,比对 {} 条,{} 条有分歧,共 {} 处\n\n分歧字段分布:\n",
        stats.n,
        stats.compared,
        by_smi.len(),
        bad.len()
    );
    for (f, c) in &by_field {
        out.push_str(&format!("  {f:<16} {c}\n"));
    }
    out.push_str("\n前若干条:\n");
    for (smi, ms) in by_smi.iter().take(limit) {
        out.push_str(&format!("  {smi}\n"));
        for m in ms.iter().take(6) {
            out.push_str(&format!(
                "      {:<22} 基准={:<8} 本实现={}\n",
                m.field, m.baseline, m.omgkit
            ));
        }
    }
    if by_smi.len() > limit {
        out.push_str(&format!("  ...(另有 {} 条)\n", by_smi.len() - limit));
    }
    out
}

#[test]
fn l2_full_smoke() {
    let (stats, bad) = diff_against(&baseline("smoke.l2.jsonl"));
    assert!(stats.compared > 0, "一条都没比对上");
    assert!(bad.is_empty(), "{}", report(&stats, &bad, 20));
    assert!(
        stats.with_explicit_hs > 0,
        "语料里没有带显式氢的分子,第 12 步是空过的"
    );
    println!(
        "L2 全管线冒烟差分通过:比对 {} 条,{} 条含显式氢",
        stats.compared, stats.with_explicit_hs
    );
}

#[test]
#[ignore = "需要先生成大语料基准;用 cargo test -- --ignored 运行"]
fn l2_full_large() {
    let (stats, bad) = diff_against(&baseline("large.l2.jsonl"));
    assert!(
        stats.compared > 1000,
        "基准不完整:只比对了 {} 条",
        stats.compared
    );
    assert!(bad.is_empty(), "{}", report(&stats, &bad, 15));
    println!(
        "L2 全管线大语料差分通过:比对 {} 条,{} 条含显式氢,已登记分歧命中 {} 条",
        stats.compared, stats.with_explicit_hs, stats.known_hit
    );
}
