//! 自由基电子数(净化第 6 步)的差分测试。
//!
//! # 触发面极窄,所以必须同时确认"确实生效"
//!
//! 语料里只有极少数分子含非零自由基,一个什么都不做的空实现同样能拿到
//! "零分歧"。由 [`radicals_actually_fire`] 守着。
//!
//! 冒烟语料因此专门补了自由基分子,逐条覆盖各个分支:满壳层反推、超价回退、
//! 无价约束元素的成键与孤立两条路。
//!
//! # 比对哪些列
//!
//! 除自由基电子数外,**价键量也一并比对**:本步之后隐式氢推断会读到非零
//! 自由基,如果两者的耦合接错了,这里就会炸出来。

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use omgkit_chem::{
    assign_radicals, clean_up, cleanup_organometallics, kekulize, perceive_rings,
    update_property_cache, ValenceResult,
};
use omgkit_core::MolBuilder;
use omgkit_io::smiles;

/// 原子列下标
const A_EXPL_H: usize = 3;
const A_NO_IMPLICIT: usize = 7;
const A_IMPLICIT_H: usize = 8;
const A_EXPLICIT_VALENCE: usize = 9;
const A_RADICALS: usize = 13;

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
    /// 含至少一个非零自由基的分子数
    with_radicals: usize,
    bad: Vec<Mismatch>,
}

fn baseline(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../harness/baseline")
        .join(name)
}

/// 净化前缀:第 1、3、4、5、6 步。
///
/// 收尾再算一次价键 —— kekulize 改过键级,而且第 6 步之后自由基电子数才非零,
/// **这一次重算才会真正读到它**。第 6 步与隐式氢推断的耦合正是靠这里暴露的。
fn run_pipeline(mol: &mut MolBuilder) -> Result<ValenceResult, String> {
    clean_up(mol);
    // 第 2 步必须在价键计算之前 —— 基准的 ops 也含它
    cleanup_organometallics(mol);
    update_property_cache(mol).map_err(|e| format!("第3步: {e}"))?;
    let _ = perceive_rings(mol);
    kekulize(mol).map_err(|e| format!("第5步: {e}"))?;
    assign_radicals(mol);
    update_property_cache(mol).map_err(|e| format!("收尾第3步: {e}"))
}

fn diff_against(path: &Path) -> DiffResult {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| {
        panic!(
            "读不到第 6 步基准 {}: {e}\n生成:python3 harness/oracle_pipeline.py --stage l2 \
             --sanitize-ops CLEANUP,PROPERTIES,SYMMRINGS,KEKULIZE,FINDRADICALS ...",
            path.display()
        )
    });

    let mut bad: Vec<Mismatch> = Vec::new();
    let (mut n, mut compared, mut with_radicals) = (0usize, 0usize, 0usize);

    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let rec: serde_json::Value = serde_json::from_str(line).expect("基准 JSONL 格式错误");
        let smi = rec["smi"].as_str().expect("缺 smi").to_string();
        n += 1;

        let Ok(mut mol) = smiles::parse(&smi) else {
            continue;
        };
        compared += 1;

        let rd_ok = rec["ok"].as_bool().unwrap_or(false);
        let ours = run_pipeline(&mut mol);

        let mut push = |field: String, baseline: String, omgkit: String| {
            bad.push(Mismatch {
                smi: smi.clone(),
                field,
                baseline,
                omgkit,
            });
        };

        let valence = match (rd_ok, &ours) {
            (false, Ok(_)) => {
                push(
                    "管线结果".into(),
                    format!("失败({})", rec["err"].as_str().unwrap_or("?")),
                    "通过".into(),
                );
                continue;
            }
            (true, Err(e)) => {
                push("管线结果".into(), "通过".into(), format!("失败({e})"));
                continue;
            }
            (false, Err(_)) => continue,
            (true, Ok(_)) => ours.expect("上面的 match 已保证是 Ok"),
        };

        if mol.atoms().iter().any(|a| a.num_radical_electrons > 0) {
            with_radicals += 1;
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
                 基准过期或列号未同步,重新生成基准(见 harness/README.md)",
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
            cmp("自由基", v[A_RADICALS], i64::from(a.num_radical_electrons));
            // 第 6 步之后隐式氢推断会读自由基;耦合接错了这里会炸
            cmp("隐式氢", v[A_IMPLICIT_H], i64::from(a.num_implicit_hs));
            cmp("显式氢", v[A_EXPL_H], i64::from(a.num_explicit_hs));
            cmp(
                "不推断隐式氢",
                v[A_NO_IMPLICIT],
                i64::from(a.flags.contains(omgkit_core::AtomFlags::NO_IMPLICIT)),
            );
            cmp(
                "显式价",
                v[A_EXPLICIT_VALENCE],
                i64::from(valence.explicit_valence[i]),
            );
        }
    }

    DiffResult {
        n,
        compared,
        with_radicals,
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
        "\nL2 第 6 步差分失败:基准 {} 条,比对 {} 条,{} 条有分歧,共 {} 处\n\n\
         分歧字段分布:\n",
        r.n,
        r.compared,
        by_smi.len(),
        r.bad.len()
    );
    for (f, c) in &by_field {
        out.push_str(&format!("  {f:<16} {c}\n"));
    }
    out.push_str("\n前若干条:\n");
    for (smi, ms) in by_smi.iter().take(limit) {
        out.push_str(&format!("  {smi}\n"));
        for m in ms.iter().take(6) {
            out.push_str(&format!(
                "      {:<22} 基准={:<20} omgkit={}\n",
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
fn l2_radicals_smoke() {
    let r = diff_against(&baseline("smoke.l2-radicals.jsonl"));
    assert!(r.compared > 0, "一条都没比对上");
    assert!(r.bad.is_empty(), "{}", report(&r, 20));
    assert!(
        r.with_radicals >= 10,
        "冒烟语料里只有 {} 条分子带自由基,该档接近空过 —— \
         语料应当覆盖 assignRadicals 的各条分支",
        r.with_radicals
    );
    println!(
        "L2 第 6 步冒烟差分通过:比对 {} 条,其中 {} 条带自由基",
        r.compared, r.with_radicals
    );
}

#[test]
#[ignore = "需要先生成大语料基准;用 cargo test -- --ignored 运行"]
fn l2_radicals_large() {
    let r = diff_against(&baseline("large.l2-radicals.jsonl"));
    assert!(r.compared > 1000, "基准不完整:只比对了 {} 条", r.compared);
    assert!(r.bad.is_empty(), "{}", report(&r, 15));
    println!(
        "L2 第 6 步大语料差分通过:比对 {} 条,其中 {} 条带自由基",
        r.compared, r.with_radicals
    );
}

/// 防止"空实现也能通过"。
///
/// 8839 条语料里只有 9 条带非零自由基。只断言"零分歧"的话,一个什么都不做的
/// `assign_radicals` 同样能过 —— 因为剩下 8830 条本来就该是 0。
#[test]
#[ignore = "需要先生成大语料基准;用 cargo test -- --ignored 运行"]
fn radicals_actually_fire() {
    let r = diff_against(&baseline("large.l2-radicals.jsonl"));
    assert_eq!(
        r.with_radicals, 9,
        "第 6 步应当在 9 条分子上产生非零自由基。\
         数字变化意味着实现或语料发生了变动,请先查清原因再更新此断言。"
    );
}
