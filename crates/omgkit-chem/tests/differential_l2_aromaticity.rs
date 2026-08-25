//! 芳香性感知(净化第 7 步)的差分测试。
//!
//! 比对原子与键的芳香标志,以及键级 —— 判定为芳香的键会从单/双键改成芳香键,
//! 键级对不上说明标记范围错了。
//!
//! 同时比对第 3 步的价键量:芳香性感知会改键级,而隐式氢是据键级算的,
//! 两者若脱节这里会暴露。
//!
//! 芳香标志是**感知**出来的,不是从输入读的:凯库勒式写法与芳香写法必须
//! 得到同样的结果,这一点由语料本身覆盖(两种写法都有)。

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use omgkit_chem::{
    assign_radicals, clean_up, cleanup_organometallics, kekulize, perceive_rings, set_aromaticity,
    update_property_cache,
};
use omgkit_core::{AtomFlags, BondFlags, MolBuilder};
use omgkit_io::smiles;

/// 原子列下标
const A_AROMATIC: usize = 6;
const A_IMPLICIT_H: usize = 8;
const A_EXPLICIT_VALENCE: usize = 9;
const A_RADICALS: usize = 13;
/// 键列下标
const B_ORDER: usize = 2;
const B_AROMATIC: usize = 4;

/// 键行的列数。与原子那侧同一个道理:列号必须与基准同步,对不上时立即炸,
/// 而不是让错位比对变成一堆无从解释的"化学分歧"。
///
/// 原子那侧一直有这道闸,**键这侧九个读取方一个都没有** —— 键元组从 6 列长到
/// 7 列(末尾的"共轭")的时候没人看得见,而下标一旦是插在中间加的,
/// 每个读取方都会静默地比错列。
const B_COLS: usize = 7;

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
    /// 含至少一个芳香原子的分子数
    with_aromatic: usize,
    bad: Vec<Mismatch>,
}

fn baseline(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../harness/baseline")
        .join(name)
}

/// 净化前缀:第 1、3、4、5、6、7 步,收尾再算一次价键。
fn run_pipeline(mol: &mut MolBuilder) -> Result<(), String> {
    clean_up(mol);
    // 第 2 步必须在价键计算之前 —— 基准的 ops 也含它
    cleanup_organometallics(mol);
    update_property_cache(mol).map_err(|e| format!("第3步: {e}"))?;
    let _ = perceive_rings(mol);
    kekulize(mol).map_err(|e| format!("第5步: {e}"))?;
    assign_radicals(mol);
    set_aromaticity(mol);
    update_property_cache(mol).map_err(|e| format!("收尾第3步: {e}"))?;
    Ok(())
}

fn diff_against(path: &Path) -> DiffResult {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| {
        panic!(
            "读不到第 7 步基准 {}: {e}\n生成:python3 harness/oracle_pipeline.py --stage l2 \
             --sanitize-ops CLEANUP,PROPERTIES,SYMMRINGS,KEKULIZE,FINDRADICALS,SETAROMATICITY ...",
            path.display()
        )
    });

    let mut bad: Vec<Mismatch> = Vec::new();
    let (mut n, mut compared, mut with_aromatic) = (0usize, 0usize, 0usize);

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

        match (rd_ok, &ours) {
            (false, Ok(())) => {
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
            (true, Ok(())) => {}
        }

        if mol
            .atoms()
            .iter()
            .any(|a| a.flags.contains(AtomFlags::AROMATIC))
        {
            with_aromatic += 1;
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
            cmp(
                "芳香",
                v[A_AROMATIC],
                i64::from(a.flags.contains(AtomFlags::AROMATIC)),
            );
            cmp("隐式氢", v[A_IMPLICIT_H], i64::from(a.num_implicit_hs));
            cmp("自由基", v[A_RADICALS], i64::from(a.num_radical_electrons));
            let _ = A_EXPLICIT_VALENCE;
        }

        for (i, row) in rec["bonds"].as_array().unwrap().iter().enumerate() {
            let v: Vec<i64> = row
                .as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_i64().unwrap())
                .collect();
            assert_eq!(
                v.len(),
                B_COLS,
                "{smi}:基准的键列数是 {},本文件按 {B_COLS} 列解读 —— \
                 基准过期或列号未同步,重新生成基准(见 harness/README.md)",
                v.len()
            );
            let b = mol.bonds()[i];
            let arom = i64::from(b.flags.contains(BondFlags::AROMATIC));
            if v[B_AROMATIC] != arom {
                push(
                    format!("键[{i}].芳香"),
                    v[B_AROMATIC].to_string(),
                    arom.to_string(),
                );
            }
            if v[B_ORDER] != b.order as i64 {
                push(
                    format!("键[{i}].键级"),
                    v[B_ORDER].to_string(),
                    (b.order as i64).to_string(),
                );
            }
        }
    }

    DiffResult {
        n,
        compared,
        with_aromatic,
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
        "\nL2 第 7 步差分失败:基准 {} 条,比对 {} 条,{} 条有分歧,共 {} 处\n\n\
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
fn l2_aromaticity_smoke() {
    let r = diff_against(&baseline("smoke.l2-aromaticity.jsonl"));
    assert!(r.compared > 0, "一条都没比对上");
    assert!(r.bad.is_empty(), "{}", report(&r, 20));
    assert!(
        r.with_aromatic > 0,
        "冒烟语料里没有任何芳香分子,该档是空过的"
    );
    println!(
        "L2 第 7 步冒烟差分通过:比对 {} 条,其中 {} 条含芳香原子",
        r.compared, r.with_aromatic
    );
}

#[test]
#[ignore = "需要先生成大语料基准;用 cargo test -- --ignored 运行"]
fn l2_aromaticity_large() {
    let r = diff_against(&baseline("large.l2-aromaticity.jsonl"));
    assert!(r.compared > 1000, "基准不完整:只比对了 {} 条", r.compared);
    assert!(r.bad.is_empty(), "{}", report(&r, 15));
    println!(
        "L2 第 7 步大语料差分通过:比对 {} 条,其中 {} 条含芳香原子",
        r.compared, r.with_aromatic
    );
}

/// 防止"什么都不标也能过"。
///
/// 语料里绝大多数分子含芳香环,一个空实现会在这里立刻暴露 —— 但仍显式断言,
/// 与其他触发面窄的步骤保持同一套防线。
#[test]
#[ignore = "需要先生成大语料基准;用 cargo test -- --ignored 运行"]
fn aromaticity_actually_fires() {
    let r = diff_against(&baseline("large.l2-aromaticity.jsonl"));
    assert_eq!(
        r.with_aromatic, 6660,
        "应有 6660 条分子含芳香原子。数字变化意味着实现或语料发生了变动,\
         请先查清原因再更新此断言。"
    );
}
