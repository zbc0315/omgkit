//! 共轭标记与杂化状态(净化第 8、9 步)的差分测试。
//!
//! 基准跑到第 9 步为止。比对:
//!
//! | 量 | 列 |
//! |---|---|
//! | 原子杂化 | 原子列 10 |
//! | 键共轭 | 键列 6 |
//! | 芳香标志、隐式氢、自由基 | 前置步骤的结果,一并守着 |
//!
//! 两步耦合:杂化在轨道数为 4 时要看有没有共轭键才决定 sp³ 还是 sp²,
//! 所以必须一起比对 —— 只比其中一个,耦合接错了也发现不了。

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use omgkit_chem::{
    assign_radicals, clean_up, cleanup_organometallics, kekulize, perceive_rings, set_aromaticity,
    set_conjugation, set_hybridization, update_property_cache,
};
use omgkit_core::{AtomFlags, BondFlags, MolBuilder};
use omgkit_io::smiles;

/// 原子列下标
const A_AROMATIC: usize = 6;
const A_IMPLICIT_H: usize = 8;
const A_HYBRID: usize = 10;
const A_RADICALS: usize = 13;
/// 键列下标
const B_AROMATIC: usize = 4;
const B_CONJUGATED: usize = 6;

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
    /// 含至少一条共轭键的分子数
    with_conjugated: usize,
    /// 出现过的杂化取值(确认各分支都被走到)
    hybrids_seen: std::collections::BTreeSet<u8>,
    bad: Vec<Mismatch>,
}

fn baseline(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../harness/baseline")
        .join(name)
}

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
    set_conjugation(mol);
    set_hybridization(mol);
    Ok(())
}

fn diff_against(path: &Path) -> DiffResult {
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("读不到第 8/9 步基准 {}: {e}", path.display()));

    let mut bad: Vec<Mismatch> = Vec::new();
    let (mut n, mut compared, mut with_conjugated) = (0usize, 0usize, 0usize);
    let mut hybrids_seen = std::collections::BTreeSet::new();

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
            .bonds()
            .iter()
            .any(|b| b.flags.contains(BondFlags::CONJUGATED))
        {
            with_conjugated += 1;
        }
        for a in mol.atoms() {
            hybrids_seen.insert(a.hybridization as u8);
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
            cmp("杂化", v[A_HYBRID], i64::from(a.hybridization as u8));
            cmp(
                "芳香",
                v[A_AROMATIC],
                i64::from(a.flags.contains(AtomFlags::AROMATIC)),
            );
            cmp("隐式氢", v[A_IMPLICIT_H], i64::from(a.num_implicit_hs));
            cmp("自由基", v[A_RADICALS], i64::from(a.num_radical_electrons));
        }

        for (i, row) in rec["bonds"].as_array().unwrap().iter().enumerate() {
            let v: Vec<i64> = row
                .as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_i64().unwrap())
                .collect();
            let b = mol.bonds()[i];
            let conj = i64::from(b.flags.contains(BondFlags::CONJUGATED));
            if v[B_CONJUGATED] != conj {
                push(
                    format!("键[{i}].共轭"),
                    v[B_CONJUGATED].to_string(),
                    conj.to_string(),
                );
            }
            let arom = i64::from(b.flags.contains(BondFlags::AROMATIC));
            if v[B_AROMATIC] != arom {
                push(
                    format!("键[{i}].芳香"),
                    v[B_AROMATIC].to_string(),
                    arom.to_string(),
                );
            }
        }
    }

    DiffResult {
        n,
        compared,
        with_conjugated,
        hybrids_seen,
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
        "\nL2 第 8/9 步差分失败:基准 {} 条,比对 {} 条,{} 条有分歧,共 {} 处\n\n\
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
fn l2_hybrid_smoke() {
    let r = diff_against(&baseline("smoke.l2-hybrid.jsonl"));
    assert!(r.compared > 0, "一条都没比对上");
    assert!(r.bad.is_empty(), "{}", report(&r, 20));
    assert!(r.with_conjugated > 0, "冒烟语料里没有共轭键,该档是空过的");

    // 配位几何的三种杂化只能由立体标记给出,电子计数推不出来。语料里各有
    // 专门的用例(`@SP` / `@TB` / `@OH`),取值缺失即说明那条分支空过。
    for (h, what) in [(5u8, "sp²d"), (6, "sp³d"), (7, "sp³d²")] {
        assert!(
            r.hybrids_seen.contains(&h),
            "杂化 {what} 在冒烟语料里没出现,配位几何那条分支是空过的。实际:{:?}",
            r.hybrids_seen
        );
    }

    println!(
        "L2 第 8/9 步冒烟差分通过:比对 {} 条,{} 条含共轭键,杂化取值 {:?}",
        r.compared, r.with_conjugated, r.hybrids_seen
    );
}

#[test]
#[ignore = "需要先生成大语料基准;用 cargo test -- --ignored 运行"]
fn l2_hybrid_large() {
    let r = diff_against(&baseline("large.l2-hybrid.jsonl"));
    assert!(r.compared > 1000, "基准不完整:只比对了 {} 条", r.compared);
    assert!(r.bad.is_empty(), "{}", report(&r, 15));
    println!(
        "L2 第 8/9 步大语料差分通过:比对 {} 条,{} 条含共轭键,杂化取值 {:?}",
        r.compared, r.with_conjugated, r.hybrids_seen
    );
}

/// 防止"只覆盖了少数几个分支也能过"。
///
/// 杂化有 8 种取值,若语料只走到其中两三种,大量分支就没有任何用例守着。
/// 这条断言把覆盖面钉死。
#[test]
#[ignore = "需要先生成大语料基准;用 cargo test -- --ignored 运行"]
fn hybridization_branches_are_covered() {
    let r = diff_against(&baseline("large.l2-hybrid.jsonl"));
    // 0=未定 1=s 2=sp 3=sp2 4=sp3 6=sp3d 7=sp3d2
    for expect in [0u8, 1, 2, 3, 4] {
        assert!(
            r.hybrids_seen.contains(&expect),
            "杂化取值 {expect} 在语料里从未出现,该分支无人守护。实际出现:{:?}",
            r.hybrids_seen
        );
    }
    assert!(
        r.with_conjugated > 1000,
        "含共轭键的分子只有 {} 条,覆盖不足",
        r.with_conjugated
    );
}
