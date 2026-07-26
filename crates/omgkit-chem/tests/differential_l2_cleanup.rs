//! 非标准画法修正(净化第 1 步)的差分测试。
//!
//! 分两档:
//!
//! | 测试 | 比对内容 |
//! |---|---|
//! | `l2_cleanup_*` | 形式电荷、键级 |
//! | `l2_cleanup_properties_*` | 上述 + 隐式氢、显式价、能否通过校验 |
//!
//! 第二档是真正的集成测试:本步改写结构、价键计算据此重算,两步串起来才
//! 对应实际的净化前缀。
//!
//! 本步的触发面很窄,只改动语料里的少数分子(全是 `N(=O)=O` 硝基)。
//! 因此单看"分歧为零"意义有限,必须同时确认它**确实改动了**那些分子,
//! 否则一个什么都不做的空实现也能通过。这一点由 `cleanup_actually_fires` 守着。

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use omgkit_chem::{clean_up, update_property_cache};
use omgkit_io::smiles;

/// 原子列下标
const A_CHARGE: usize = 1;
const A_IMPLICIT_H: usize = 8;
const A_EXPLICIT_VALENCE: usize = 9;
/// 键列下标
const B_ORDER: usize = 2;

/// l2 原子行的列数。基准与本文件的列号必须同步 —— 对不上时立即炸,
/// 而不是让错位比对变成一堆无从解释的"化学分歧"。新列一律追加到行尾,
/// 见 harness/README.md 的列规范。
const A_COLS: usize = 15;

struct Mismatch {
    smi: String,
    field: String,
    baseline: String,
    omgkit: String,
}

struct DiffResult {
    n: usize,
    compared: usize,
    /// omgkit 相对 L1 解析结果确实发生了改动的分子数
    changed: usize,
    bad: Vec<Mismatch>,
}

fn baseline(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../harness/baseline")
        .join(name)
}

/// `with_properties` 为真时额外跑第 3 步并比对价键相关列。
fn diff_against(path: &Path, with_properties: bool) -> DiffResult {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| {
        panic!(
            "读不到基准 {}: {e}\n生成方式见 harness/README.md",
            path.display()
        )
    });

    let mut bad: Vec<Mismatch> = Vec::new();
    let (mut n, mut compared, mut changed) = (0usize, 0usize, 0usize);

    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let rec: serde_json::Value = serde_json::from_str(line).expect("基准 JSONL 格式错误");
        let smi = rec["smi"].as_str().expect("缺 smi").to_string();
        n += 1;

        // 解析失败的由 L1 的差分测试守着
        let Ok(mut mol) = smiles::parse(&smi) else {
            continue;
        };
        compared += 1;

        let before = mol.clone();
        clean_up(&mut mol);
        if mol.atoms() != before.atoms() || mol.bonds() != before.bonds() {
            changed += 1;
        }

        let mut push = |field: String, baseline: String, omgkit: String| {
            bad.push(Mismatch {
                smi: smi.clone(),
                field,
                baseline,
                omgkit,
            });
        };

        let rd_ok = rec["ok"].as_bool().unwrap_or(false);
        let valence = if with_properties {
            let r = update_property_cache(&mut mol);
            match (rd_ok, &r) {
                (false, Ok(_)) => {
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
                (true, Ok(_)) => {}
            }
            Some(r.expect("上面的 match 已保证是 Ok"))
        } else {
            assert!(rd_ok, "只跑 CLEANUP 不应失败,但基准里 {smi} 标为失败");
            None
        };

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
            let a = mol.atoms()[i];

            if vals[A_CHARGE] != i64::from(a.formal_charge) {
                push(
                    format!("原子[{i}].电荷"),
                    vals[A_CHARGE].to_string(),
                    a.formal_charge.to_string(),
                );
            }
            if let Some(v) = &valence {
                let ih = i64::from(v.implicit_hs[i]);
                if vals[A_IMPLICIT_H] != ih {
                    push(
                        format!("原子[{i}].隐式氢"),
                        vals[A_IMPLICIT_H].to_string(),
                        ih.to_string(),
                    );
                }
                let ev = i64::from(v.explicit_valence[i]);
                if vals[A_EXPLICIT_VALENCE] != ev {
                    push(
                        format!("原子[{i}].显式价"),
                        vals[A_EXPLICIT_VALENCE].to_string(),
                        ev.to_string(),
                    );
                }
            }
        }

        for (i, row) in rec["bonds"].as_array().unwrap().iter().enumerate() {
            let vals: Vec<i64> = row
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_i64().unwrap())
                .collect();
            let order = mol.bonds()[i].order as i64;
            if vals[B_ORDER] != order {
                push(
                    format!("键[{i}].键级"),
                    vals[B_ORDER].to_string(),
                    order.to_string(),
                );
            }
        }
    }

    DiffResult {
        n,
        compared,
        changed,
        bad,
    }
}

fn report(r: &DiffResult, limit: usize) -> String {
    let mut by_field: BTreeMap<&str, usize> = BTreeMap::new();
    for m in &r.bad {
        *by_field
            .entry(m.field.split('.').next_back().unwrap_or(&m.field))
            .or_default() += 1;
    }
    let mut by_smi: BTreeMap<&str, Vec<&Mismatch>> = BTreeMap::new();
    for m in &r.bad {
        by_smi.entry(&m.smi).or_default().push(m);
    }

    let mut out = format!(
        "\nL2 第 1 步差分失败:基准 {} 条,比对 {} 条,{} 条有分歧,共 {} 处\n\n\
         分歧字段分布:\n",
        r.n,
        r.compared,
        by_smi.len(),
        r.bad.len()
    );
    for (field, count) in &by_field {
        out.push_str(&format!("  {field:<14} {count}\n"));
    }
    out.push_str("\n前若干条:\n");
    for (smi, ms) in by_smi.iter().take(limit) {
        out.push_str(&format!("  {smi}\n"));
        for m in ms.iter().take(6) {
            out.push_str(&format!(
                "      {:<20} 基准={:<26} omgkit={}\n",
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
fn l2_cleanup_smoke() {
    let r = diff_against(&baseline("smoke.l2-cleanup.jsonl"), false);
    assert!(r.compared > 0, "一条都没比对上");
    assert!(r.bad.is_empty(), "{}", report(&r, 20));
    println!(
        "L2 第 1 步冒烟差分通过:比对 {} 条,改动 {} 条",
        r.compared, r.changed
    );
}

#[test]
fn l2_cleanup_properties_smoke() {
    let r = diff_against(&baseline("smoke.l2-cleanup-properties.jsonl"), true);
    assert!(r.compared > 0, "一条都没比对上");
    assert!(r.bad.is_empty(), "{}", report(&r, 20));
    println!(
        "L2 第 1+3 步冒烟差分通过:比对 {} 条,改动 {} 条",
        r.compared, r.changed
    );
}

#[test]
#[ignore = "需要先生成大语料基准;用 cargo test -- --ignored 运行"]
fn l2_cleanup_large() {
    let r = diff_against(&baseline("large.l2-cleanup.jsonl"), false);
    assert!(r.compared > 1000, "基准不完整:只比对了 {} 条", r.compared);
    assert!(r.bad.is_empty(), "{}", report(&r, 15));
    println!(
        "L2 第 1 步大语料差分通过:比对 {} 条,改动 {} 条",
        r.compared, r.changed
    );
}

#[test]
#[ignore = "需要先生成大语料基准;用 cargo test -- --ignored 运行"]
fn l2_cleanup_properties_large() {
    let r = diff_against(&baseline("large.l2-cleanup-properties.jsonl"), true);
    assert!(r.compared > 1000, "基准不完整:只比对了 {} 条", r.compared);
    assert!(r.bad.is_empty(), "{}", report(&r, 15));
    println!(
        "L2 第 1+3 步大语料差分通过:比对 {} 条,改动 {} 条",
        r.compared, r.changed
    );
}

/// 防止"空实现也能通过"。
///
/// 第 1 步触发面很窄:8839 条语料里只改动 15 条。如果只断言"零分歧",
/// 一个什么都不做的 `clean_up` 同样能通过 —— 因为剩下 8824 条本来就不该变。
/// 所以必须显式确认它**确实改动了预期数量**的分子。
#[test]
#[ignore = "需要先生成大语料基准;用 cargo test -- --ignored 运行"]
fn cleanup_actually_fires() {
    let r = diff_against(&baseline("large.l2-cleanup.jsonl"), false);
    assert_eq!(
        r.changed, 15,
        "第 1 步应当改动 15 条分子。\
         数字变化意味着实现或语料发生了变动,请先查清原因再更新此断言。"
    );
}
