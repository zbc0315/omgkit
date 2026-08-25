//! kekulize(净化第 5 步)的差分测试:与参考实现逐字段比对。
//!
//! # 比对哪些量
//!
//! 一个芳香体系通常有多个同样合法的 Kekulé 式(萘有两个),选中哪一个由
//! 遍历顺序决定,没有化学含义。所以只比对**化学上确定**的量:
//!
//! | 量 | 是否比对 | 理由 |
//! |---|---|---|
//! | 芳香标志已清除 | 比 | 化学决定 |
//! | 显式价 / 隐式氢 | 比 | 化学决定,且跨不同 Kekulé 式不变 |
//! | 显式氢 / 不推断隐式氢 | 比 | 吡咯型 `[nH]` 的处理 |
//! | 每个原子的双键数 | 比 | 完美匹配覆盖同一批顶点,故也不变 |
//! | **逐键下标的键级** | 不比 | 取决于遍历顺序,不同解都合法 |
//!
//! 拿逐键键级当判据会产生大量假失败:化学量全对,只是选了另一个等价解。
//!
//! # 不依赖参考实现的性质
//!
//! - **确定性**:同一输入必须给出同一结果,否则规范化输出无从稳定
//! - **完备性**:存在合法结构就必须找到,且"搜索未尽"与"确实无解"要分开报

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use omgkit_chem::{
    clean_up, cleanup_organometallics, kekulize, perceive_rings, update_property_cache,
    ValenceResult,
};
use omgkit_core::{AtomFlags, BondFlags, MolBuilder};
use omgkit_io::smiles;

/// 原子列下标
const A_EXPL_H: usize = 3;
const A_AROMATIC: usize = 6;
const A_NO_IMPLICIT: usize = 7;
const A_IMPLICIT_H: usize = 8;
const A_EXPLICIT_VALENCE: usize = 9;
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
    /// 键级确实被 kekulize 改动过的分子数
    changed: usize,
    bad: Vec<Mismatch>,
}

fn baseline(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../harness/baseline")
        .join(name)
}

/// omgkit 侧的四步前缀。返回 `Err` 表示某一步失败。
fn run_pipeline(mol: &mut MolBuilder) -> Result<(ValenceResult, bool), String> {
    clean_up(mol);
    // 第 2 步必须在价键计算之前 —— 基准的 ops 也含它
    cleanup_organometallics(mol);
    update_property_cache(mol).map_err(|e| format!("第3步: {e}"))?;
    let _ = perceive_rings(mol);

    // 只测**第 5 步本身**改动了什么 —— 从 kekulize 之前取快照。
    // 从解析结果取快照会把第 1 步 clean_up 的改动也算进来,数字就对不上了。
    let before: Vec<_> = mol.bonds().iter().map(|b| b.order).collect();
    kekulize(mol).map_err(|e| format!("第5步: {e}"))?;
    let kekulize_changed = mol.bonds().iter().map(|b| b.order).ne(before);

    // kekulize 改了键级,价键要重算后才能与基准比对
    let v = update_property_cache(mol).map_err(|e| format!("收尾第3步: {e}"))?;
    Ok((v, kekulize_changed))
}

fn diff_against(path: &Path) -> DiffResult {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| {
        panic!(
            "读不到 kekulize 基准 {}: {e}\n生成:python3 harness/oracle_pipeline.py \
             --stage l2 --sanitize-ops CLEANUP,PROPERTIES,SYMMRINGS,KEKULIZE ...",
            path.display()
        )
    });

    let mut bad: Vec<Mismatch> = Vec::new();
    let (mut n, mut compared, mut changed) = (0usize, 0usize, 0usize);

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
        let (valence, kekulize_changed) = valence;
        if kekulize_changed {
            changed += 1;
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
            cmp("显式氢", v[A_EXPL_H], i64::from(a.num_explicit_hs));
            cmp(
                "不推断隐式氢",
                v[A_NO_IMPLICIT],
                i64::from(a.flags.contains(AtomFlags::NO_IMPLICIT)),
            );
            cmp("隐式氢", v[A_IMPLICIT_H], i64::from(valence.implicit_hs[i]));
            cmp(
                "显式价",
                v[A_EXPLICIT_VALENCE],
                valence.explicit_valence[i].into(),
            );
        }

        // 键:只比对化学上确定的量。
        // 逐键下标的**键级**刻意不比对 —— 见文件头。
        let n_atoms = mol.num_atoms();
        let mut rd_dbl = vec![0i64; n_atoms];
        let mut our_dbl = vec![0i64; n_atoms];
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

            // 芳香标志必须全部清除
            let arom = i64::from(b.flags.contains(BondFlags::AROMATIC));
            if v[B_AROMATIC] != arom {
                push(
                    format!("键[{i}].芳香"),
                    v[B_AROMATIC].to_string(),
                    arom.to_string(),
                );
            }
            if arom != 0 || v[B_AROMATIC] != 0 {
                push(format!("键[{i}].残留芳香"), "0".into(), "非0".into());
            }
            // 键级本身不比对,但双键的**分布**要比 —— 完美匹配覆盖同一批原子
            if v[B_ORDER] == 2 {
                rd_dbl[v[0] as usize] += 1;
                rd_dbl[v[1] as usize] += 1;
            }
            if b.order == omgkit_core::BondOrder::Double {
                our_dbl[b.begin as usize] += 1;
                our_dbl[b.end as usize] += 1;
            }
        }
        for i in 0..n_atoms {
            if rd_dbl[i] != our_dbl[i] {
                push(
                    format!("原子[{i}].双键数"),
                    rd_dbl[i].to_string(),
                    our_dbl[i].to_string(),
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
        "\nL2 第 5 步差分失败:基准 {} 条,比对 {} 条,{} 条有分歧,共 {} 处\n\n\
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
                "      {:<22} 基准={:<26} omgkit={}\n",
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
fn l2_kekulize_smoke() {
    let r = diff_against(&baseline("smoke.l2-kekulize.jsonl"));
    assert!(r.compared > 0, "一条都没比对上");
    assert!(r.bad.is_empty(), "{}", report(&r, 20));
    assert!(
        r.changed > 0,
        "冒烟语料里没有分子被 kekulize 改动,该档是空过的"
    );
    println!(
        "L2 第 5 步冒烟差分通过:比对 {} 条,改动 {} 条",
        r.compared, r.changed
    );
}

#[test]
#[ignore = "需要先生成大语料基准;用 cargo test -- --ignored 运行"]
fn l2_kekulize_large() {
    let r = diff_against(&baseline("large.l2-kekulize.jsonl"));
    assert!(r.compared > 1000, "基准不完整:只比对了 {} 条", r.compared);
    assert!(r.bad.is_empty(), "{}", report(&r, 15));
    assert_eq!(
        r.changed, 933,
        "kekulize 应当改动 933 条分子。\
         数字变化意味着实现或语料发生了变动,请先查清原因。"
    );
    println!(
        "L2 第 5 步大语料差分通过:比对 {} 条,改动 {} 条",
        r.compared, r.changed
    );
}
