//! 价键计算与隐式氢推断(净化第 3 步)的差分测试。
//!
//! # 为什么必须用"只跑第 3 步"的基准
//!
//! 完整的 l2 基准反映的是 **12 步全跑完**的最终状态:芳香性感知会改写芳香
//! 标志,收尾还会再算一次价键。凯库勒式写法的分子在"只跑第 3 步"和"全跑完"
//! 两种状态下,隐式氢/显式价是**不同**的。拿全量基准验证单步会得到一堆看不懂
//! 的分歧。
//!
//! 所以基准生成器支持 `--sanitize-ops`,能只跑管线的任意子集:
//!
//! ```sh
//! python3 harness/oracle_pipeline.py --input harness/corpus/smoke.smi \
//!     --stage l2 --sanitize-ops PROPERTIES --out harness/baseline/smoke.l2-properties.jsonl
//! ```
//!
//! # 比对的内容
//!
//! | 列 | 来源 |
//! |---|---|
//! | 原子 `隐式氢` | l2 原子列 8 |
//! | 原子 `显式价` | l2 原子列 9 |
//! | **能否通过第 3 步** | 基准的成败标记 vs 本实现返回的 `Err` |
//!
//! 第三项同样重要:它验证严格模式下超价判定的边界条件 —— 语料里有若干分子
//! 本就该在这一步失败,漏判或误判都会立刻暴露。

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use omgkit_chem::update_property_cache;
use omgkit_io::smiles;

/// l2 原子列下标
const A_IMPLICIT_H: usize = 8;
const A_EXPLICIT_VALENCE: usize = 9;

/// l2 原子行的列数。基准与本文件的列号必须同步 —— 对不上时立即炸,
/// 而不是让错位比对变成一堆无从解释的"化学分歧"。新列一律追加到行尾,
/// 见 harness/README.md 的列规范。
const A_COLS: usize = 15;

/// 已知且已定位根因的分歧。与环感知测试同一套机制:登记项若不再分歧,
/// 测试会失败并要求删除 —— 防止名单沉淀成暗坑。
const KNOWN_DIVERGENCES: &[(&str, &str)] = &[];

fn is_known(smi: &str) -> bool {
    KNOWN_DIVERGENCES.iter().any(|(s, _)| *s == smi)
}

struct Mismatch {
    smi: String,
    field: String,
    baseline: String,
    omgkit: String,
}

struct DiffResult {
    n: usize,
    compared: usize,
    unexpected: Vec<Mismatch>,
    hit_known: Vec<String>,
    known_present: Vec<String>,
}

fn baseline(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../harness/baseline")
        .join(name)
}

fn diff_against(path: &Path) -> DiffResult {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| {
        panic!(
            "读不到第 3 步基准 {}: {e}\n生成:python3 harness/oracle_pipeline.py \
             --stage l2 --sanitize-ops PROPERTIES --input <corpus> --out <此文件>",
            path.display()
        )
    });

    let mut bad: Vec<Mismatch> = Vec::new();
    let mut known_present = Vec::new();
    let (mut n, mut compared) = (0usize, 0usize);

    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let rec: serde_json::Value = serde_json::from_str(line).expect("基准 JSONL 格式错误");
        let smi = rec["smi"].as_str().expect("缺 smi").to_string();
        n += 1;
        if is_known(&smi) {
            known_present.push(smi.clone());
        }

        // 解析失败的由 L1 的差分测试守着,此处不重复
        let Ok(mut mol) = smiles::parse(&smi) else {
            continue;
        };
        compared += 1;

        let rd_ok = rec["ok"].as_bool().unwrap_or(false);
        let ours = update_property_cache(&mut mol);

        let mut push = |field: &str, baseline: String, omgkit: String| {
            bad.push(Mismatch {
                smi: smi.clone(),
                field: field.to_string(),
                baseline,
                omgkit,
            });
        };

        // -- 通过 / 失败必须一致 --
        match (rd_ok, &ours) {
            (false, Ok(_)) => {
                push(
                    "第3步结果",
                    format!("失败({})", rec["err"].as_str().unwrap_or("?")),
                    "通过".into(),
                );
                continue;
            }
            (true, Err(e)) => {
                push("第3步结果", "通过".into(), format!("失败({e})"));
                continue;
            }
            (false, Err(_)) => continue, // 双方都失败,无需逐列比对
            (true, Ok(_)) => {}
        }
        let r = ours.expect("上面的 match 已保证是 Ok");

        for (i, row) in rec["atoms"].as_array().unwrap().iter().enumerate() {
            let vals: Vec<i64> = row
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_i64().unwrap())
                .collect();
            assert_eq!(
                vals.len(),
                A_COLS,
                "{smi}:基准的原子列数是 {},本文件按 {A_COLS} 列解读 —— \
                 基准过期或列号未同步,重新生成基准(见 harness/README.md)",
                vals.len()
            );
            let ih = i64::from(r.implicit_hs[i]);
            if vals[A_IMPLICIT_H] != ih {
                push(
                    &format!("原子[{i}].隐式氢"),
                    vals[A_IMPLICIT_H].to_string(),
                    ih.to_string(),
                );
            }
            let ev = i64::from(r.explicit_valence[i]);
            if vals[A_EXPLICIT_VALENCE] != ev {
                push(
                    &format!("原子[{i}].显式价"),
                    vals[A_EXPLICIT_VALENCE].to_string(),
                    ev.to_string(),
                );
            }
        }
    }

    let mut hit_known: Vec<String> = bad
        .iter()
        .filter(|m| is_known(&m.smi))
        .map(|m| m.smi.clone())
        .collect();
    hit_known.sort();
    hit_known.dedup();
    bad.retain(|m| !is_known(&m.smi));

    DiffResult {
        n,
        compared,
        unexpected: bad,
        hit_known,
        known_present,
    }
}

fn report(r: &DiffResult, limit: usize) -> String {
    let mut by_field: BTreeMap<&str, usize> = BTreeMap::new();
    for m in &r.unexpected {
        *by_field
            .entry(m.field.split('.').next_back().unwrap_or(&m.field))
            .or_default() += 1;
    }
    let mut by_smi: BTreeMap<&str, Vec<&Mismatch>> = BTreeMap::new();
    for m in &r.unexpected {
        by_smi.entry(&m.smi).or_default().push(m);
    }

    let mut out = format!(
        "\nL2 第 3 步差分失败:基准 {} 条,比对 {} 条,其中 {} 条有分歧,共 {} 处\n\n\
         分歧字段分布:\n",
        r.n,
        r.compared,
        by_smi.len(),
        r.unexpected.len()
    );
    for (field, count) in &by_field {
        out.push_str(&format!("  {field:<14} {count}\n"));
    }
    out.push_str("\n前若干条:\n");
    for (smi, ms) in by_smi.iter().take(limit) {
        out.push_str(&format!("  {smi}\n"));
        for m in ms.iter().take(6) {
            out.push_str(&format!(
                "      {:<22} 基准={:<28} omgkit={}\n",
                m.field, m.baseline, m.omgkit
            ));
        }
    }
    if by_smi.len() > limit {
        out.push_str(&format!("  ...(另有 {} 条)\n", by_smi.len() - limit));
    }
    out
}

fn assert_clean(r: &DiffResult, limit: usize) {
    assert!(r.unexpected.is_empty(), "{}", report(r, limit));

    let mut stale: Vec<&String> = r
        .known_present
        .iter()
        .filter(|s| !r.hit_known.contains(s))
        .collect();
    stale.sort();
    stale.dedup();
    assert!(
        stale.is_empty(),
        "KNOWN_DIVERGENCES 中有条目已不再产生分歧,请删除:\n{}",
        stale
            .iter()
            .map(|s| format!("  {s}"))
            .collect::<Vec<_>>()
            .join("\n")
    );
}

#[test]
fn l2_valence_smoke() {
    let r = diff_against(&baseline("smoke.l2-properties.jsonl"));
    assert!(r.compared > 0, "一条都没比对上");
    assert_clean(&r, 20);
    println!("L2 第 3 步冒烟差分通过:比对 {} 条,零未登记分歧", r.compared);
}

/// 大语料(~8800 条)。生成基准见 `harness/README.md`。
#[test]
#[ignore = "需要先生成大语料基准;用 cargo test -- --ignored 运行"]
fn l2_valence_large() {
    let r = diff_against(&baseline("large.l2-properties.jsonl"));
    assert!(r.compared > 1000, "基准不完整:只比对了 {} 条", r.compared);
    assert_clean(&r, 15);
    println!(
        "L2 第 3 步大语料差分通过:比对 {} 条,零未登记分歧",
        r.compared
    );
}
