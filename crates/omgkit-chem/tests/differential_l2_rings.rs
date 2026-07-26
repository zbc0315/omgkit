//! 环感知(净化第 4 步)的差分测试。
//!
//! 只比对环相关的三列:原子"在环中"、原子"最小环大小"、键"在环中"。
//!
//! 前提是**净化不改变连通性**(已在两份语料上验证:规模改变 0、无向边集
//! 改变 0),所以可以在解析结果的图上直接做环感知,再与净化后的基准比对。
//! 唯一的例外是第 2 步会交换某些键的端点并改成配位键,但那不影响无向边集,
//! 也不改变键的下标顺序。

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use omgkit_chem::{clean_up, cleanup_organometallics, perceive_rings};
use omgkit_io::smiles;

/// l2 原子列下标
const A_IN_RING: usize = 11;
const A_MIN_RING: usize = 12;
/// l2 键列下标
const B_IN_RING: usize = 5;

/// l2 原子行的列数。基准与本文件的列号必须同步 —— 对不上时立即炸,
/// 而不是让错位比对变成一堆无从解释的"化学分歧"。新列一律追加到行尾,
/// 见 harness/README.md 的列规范。
const A_COLS: usize = 15;

/// 已知且**已定位根因**的分歧。
///
/// 这不是豁免名单:测试会检查每个条目是否**仍然**产生分歧,若某条目已不再
/// 分歧就直接失败并要求删除它。一个永远躺着没人动的豁免名单本身就是暗坑。
///
/// 当前为空。登记条目时要留意一类假分歧:基准的 `--sanitize-ops` 若与本测试
/// 实际跑的步骤对不上,差出来的是步骤而不是实现 —— 环感知依赖的图正是被
/// 第 2 步(有机金属键改配位键)改过的。
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

fn baseline(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../harness/baseline")
        .join(name)
}

/// 一次比对的产出。
struct DiffResult {
    /// 基准记录数
    n: usize,
    /// 实际比对的分子数
    compared: usize,
    /// 未登记的分歧 —— 必须为空
    unexpected: Vec<Mismatch>,
    /// 命中 [`KNOWN_DIVERGENCES`] 的分子(去重)
    hit_known: Vec<String>,
    /// 基准中出现过的、已登记的 SMILES(用于检测过期条目)
    known_present: Vec<String>,
}

fn diff_against(path: &Path) -> DiffResult {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| {
        panic!(
            "读不到 L2 基准 {}: {e}\n生成:python3 harness/oracle_pipeline.py \
             --stage l2 --input <corpus> --out {}",
            path.display(),
            path.display()
        )
    });

    let mut bad = Vec::new();
    let mut known_present: Vec<String> = Vec::new();
    let (mut n, mut compared) = (0usize, 0usize);

    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let rec: serde_json::Value = serde_json::from_str(line).expect("基准 JSONL 格式错误");
        let smi = rec["smi"].as_str().expect("缺 smi").to_string();
        n += 1;
        if is_known(&smi) {
            known_present.push(smi.clone());
        }

        // 基准里净化失败的分子无从比对
        if !rec["ok"].as_bool().unwrap_or(false) {
            continue;
        }
        let Ok(mut mol) = smiles::parse(&smi) else {
            // L1 已有独立的差分测试守着可解析性,此处不重复报
            continue;
        };

        let na = rec["na"].as_u64().unwrap() as usize;
        let nb = rec["nb"].as_u64().unwrap() as usize;
        if mol.num_atoms() != na || mol.num_bonds() != nb {
            bad.push(Mismatch {
                smi: smi.clone(),
                field: "规模".into(),
                baseline: format!("{na}原子/{nb}键"),
                omgkit: format!("{}原子/{}键", mol.num_atoms(), mol.num_bonds()),
            });
            continue;
        }

        // 基准是**全管线**跑完的状态,而环感知只依赖图。净化里只有第 1、2 步
        // 会改动图(第 2 步把某些键改成配位键并交换端点),所以补这两步就够 ——
        // 少了第 2 步,有机金属分子的图与基准不是同一个图。
        clean_up(&mut mol);
        cleanup_organometallics(&mut mol);
        let r = perceive_rings(&mut mol);
        compared += 1;

        let mut push = |field: String, baseline: String, omgkit: String| {
            bad.push(Mismatch {
                smi: smi.clone(),
                field,
                baseline,
                omgkit,
            });
        };

        for (i, row) in rec["atoms"].as_array().unwrap().iter().enumerate() {
            let r_row: Vec<i64> = row
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_i64().unwrap())
                .collect();
            assert_eq!(
                r_row.len(),
                A_COLS,
                "{smi}:基准的原子列数是 {},本文件按 {A_COLS} 列解读 —— \
                 基准过期或列号未同步,重新生成基准(见 harness/README.md)",
                r_row.len()
            );
            let in_ring = i64::from(r.atom_in_ring[i]);
            if r_row[A_IN_RING] != in_ring {
                push(
                    format!("原子[{i}].在环中"),
                    r_row[A_IN_RING].to_string(),
                    in_ring.to_string(),
                );
            }
            let size = i64::from(r.atom_min_ring_size[i]);
            if r_row[A_MIN_RING] != size {
                push(
                    format!("原子[{i}].最小环"),
                    r_row[A_MIN_RING].to_string(),
                    size.to_string(),
                );
            }
        }

        for (i, row) in rec["bonds"].as_array().unwrap().iter().enumerate() {
            let r_row: Vec<i64> = row
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_i64().unwrap())
                .collect();
            let in_ring = i64::from(r.bond_in_ring[i]);
            if r_row[B_IN_RING] != in_ring {
                push(
                    format!("键[{i}].在环中"),
                    r_row[B_IN_RING].to_string(),
                    in_ring.to_string(),
                );
            }
        }
    }

    let mut hit_known: Vec<String> = bad
        .iter()
        .filter(|m: &&Mismatch| is_known(&m.smi))
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

fn report(n: usize, compared: usize, bad: &[Mismatch], limit: usize) -> String {
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
        "\nL2 环感知差分失败:基准 {n} 条,实际比对 {compared} 条,\
         其中 {} 条有分歧,共 {} 处\n\n分歧字段分布:\n",
        by_smi.len(),
        bad.len()
    );
    for (field, count) in &by_field {
        out.push_str(&format!("  {field:<12} {count}\n"));
    }
    out.push_str("\n前若干条:\n");
    for (smi, ms) in by_smi.iter().take(limit) {
        out.push_str(&format!("  {smi}\n"));
        for m in ms.iter().take(6) {
            out.push_str(&format!(
                "      {:<20} 基准={:<6} omgkit={}\n",
                m.field, m.baseline, m.omgkit
            ));
        }
    }
    if by_smi.len() > limit {
        out.push_str(&format!("  ...(另有 {} 条)\n", by_smi.len() - limit));
    }
    out
}

/// 断言:未登记的分歧为零,且每个出现在本语料中的已登记条目**仍然**分歧。
///
/// 后一半是关键 —— 它保证 [`KNOWN_DIVERGENCES`] 不会在根因修好后
/// 悄悄留存下来变成暗坑。
fn assert_clean(r: &DiffResult, limit: usize) {
    assert!(
        r.unexpected.is_empty(),
        "{}",
        report(r.n, r.compared, &r.unexpected, limit)
    );

    let mut stale: Vec<&String> = r
        .known_present
        .iter()
        .filter(|s| !r.hit_known.contains(s))
        .collect();
    stale.sort();
    stale.dedup();
    assert!(
        stale.is_empty(),
        "KNOWN_DIVERGENCES 中有条目已不再产生分歧,请从名单里删除:\n{}",
        stale
            .iter()
            .map(|s| format!("  {s}"))
            .collect::<Vec<_>>()
            .join("\n")
    );
}

#[test]
fn l2_rings_smoke() {
    let r = diff_against(&baseline("smoke.l2.jsonl"));
    assert!(r.compared > 0, "一条都没比对上");
    assert_clean(&r, 20);
    println!(
        "L2 环感知冒烟差分通过:比对 {} 条,未登记分歧 0(已登记 {} 条)",
        r.compared,
        r.hit_known.len()
    );
}

/// 大语料(~8800 条)。生成基准见 `harness/README.md`。
#[test]
#[ignore = "需要先生成大语料基准;用 cargo test -- --ignored 运行"]
fn l2_rings_large() {
    let r = diff_against(&baseline("large.l2.jsonl"));
    assert!(
        r.compared > 1000,
        "大语料基准看起来不完整:只比对了 {} 条",
        r.compared
    );
    assert_clean(&r, 15);
    println!(
        "L2 环感知大语料差分通过:比对 {} 条,未登记分歧 0(已登记 {} 条:{:?})",
        r.compared,
        r.hit_known.len(),
        r.hit_known
    );
}
