//! SMILES 解析的差分测试。
//!
//! 基准由 `harness/oracle_pipeline.py --stage l1` 生成。两档:
//!
//! - [`l1_smoke`] —— 冒烟语料(含非法输入),基准已入库,**始终运行**。
//!   一个会静默跳过的差分测试本身就是暗坑。
//! - [`l1_large`] —— 大语料(~8800 条),基准体积大不入库,
//!   标记为 `#[ignore]`,用 `cargo test -- --ignored` 显式运行。
//!
//! 比对的列见 `harness/README.md` 的列规范表。**键的"在环中"列不比对** ——
//! 环感知属于净化,不属于解析。

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use omgkit_core::{AtomFlags, BondFlags};
use omgkit_io::smiles;

/// l1 原子列:[Z, 电荷, 同位素, 显式氢, 映射号, 立体类别, 芳香, 不推断隐式氢,
/// 立体排列序号]
const A_Z: usize = 0;
const A_CHARGE: usize = 1;
const A_ISOTOPE: usize = 2;
const A_EXPL_H: usize = 3;
const A_MAP: usize = 4;
const A_CHIRAL: usize = 5;
const A_AROMATIC: usize = 6;
const A_NO_IMPLICIT: usize = 7;
/// 第 8 列是立体排列序号,**不比对** —— 见下方 `cmp` 处的说明
/// l1 原子行的列数。基准与本文件的列号必须同步 —— 对不上时立即炸,
/// 而不是让越界或错位比对变成一堆无从解释的分歧。
const A_COLS: usize = 9;

/// l1 键列:[起点, 终点, 键级, 方向, 芳香, 在环中, 共轭] —— **7 列**。
///
/// 末两列 L1 都不比:环感知与共轭感知都属于净化(L2),这里比它们
/// 只会把净化的缺陷记到解析头上。到 L2 由 `B_IN_RING` / `B_CONJUGATED` 比。
/// (这一行先前写的是 6 列,而基准一直是 7 列 —— 见 `B_COLS`。)
const B_BEGIN: usize = 0;
const B_END: usize = 1;
const B_ORDER: usize = 2;
const B_DIR: usize = 3;
const B_AROMATIC: usize = 4;

/// 键行的列数。与原子那侧同一个道理:列号必须与基准同步,对不上时立即炸,
/// 而不是让错位比对变成一堆无从解释的"化学分歧"。
///
/// 原子那侧一直有这道闸,**键这侧九个读取方一个都没有** —— 键元组从 6 列长到
/// 7 列(末尾的"共轭")的时候没人看得见,而下标一旦是插在中间加的,
/// 每个读取方都会静默地比错列。
const B_COLS: usize = 7;

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

/// 语法特性的覆盖计数。
///
/// 差分测试全绿有两种可能:实现对,或者语料压根没走到那条路径。这些计数把
/// 两者分开 —— 它们本身也是断言的对象。
#[derive(Default)]
struct Coverage {
    /// 配位键条数
    dative_bonds: usize,
    /// 非四面体立体标记的原子数(`@SP` / `@TB` / `@OH` / `@AL`)
    nontetrahedral_atoms: usize,
    /// 排列序号非零的原子数
    stereo_perms: usize,
    /// 四面体立体标记的原子数
    tetrahedral_atoms: usize,
}

/// 对一份基准跑完整比对,返回 (记录数, 分歧列表, 覆盖计数)。
fn diff_against(path: &Path) -> (usize, Vec<Mismatch>, Coverage) {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| {
        panic!(
            "读不到 L1 基准 {}: {e}\n生成方式见 harness/README.md",
            path.display()
        )
    });

    let mut bad = Vec::new();
    let mut n = 0usize;
    let mut cov = Coverage::default();

    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let rec: serde_json::Value = serde_json::from_str(line).expect("基准 JSONL 格式错误");
        let smi = rec["smi"].as_str().expect("缺 smi 字段").to_string();
        n += 1;

        let mut push = |field: String, baseline: String, omgkit: String| {
            bad.push(Mismatch {
                smi: smi.clone(),
                field,
                baseline,
                omgkit,
            });
        };

        let rd_ok = rec["ok"].as_bool().unwrap_or(false);
        let parsed = smiles::parse(&smi);

        if !rd_ok {
            if parsed.is_ok() {
                push("可解析性".into(), "失败".into(), "成功".into());
            }
            continue;
        }
        let mol = match parsed {
            Ok(m) => m,
            Err(e) => {
                push("可解析性".into(), "成功".into(), format!("失败: {e}"));
                continue;
            }
        };

        // -- 规模 --
        let na = rec["na"].as_u64().unwrap() as usize;
        let nb = rec["nb"].as_u64().unwrap() as usize;
        if mol.num_atoms() != na {
            push("原子数".into(), na.to_string(), mol.num_atoms().to_string());
            continue; // 规模不符,逐项比对无意义
        }
        if mol.num_bonds() != nb {
            push("键数".into(), nb.to_string(), mol.num_bonds().to_string());
            continue;
        }

        cov.dative_bonds += mol
            .bonds()
            .iter()
            .filter(|b| b.order == omgkit_core::BondOrder::Dative)
            .count();
        for a in mol.atoms() {
            if a.chiral_tag.is_tetrahedral() {
                cov.tetrahedral_atoms += 1;
            } else if a.chiral_tag != omgkit_core::ChiralTag::Unspecified {
                cov.nontetrahedral_atoms += 1;
            }
            if a.stereo_perm != 0 {
                cov.stereo_perms += 1;
            }
        }

        // -- 原子列 --
        for (i, row) in rec["atoms"].as_array().unwrap().iter().enumerate() {
            let r: Vec<i64> = row
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_i64().unwrap())
                .collect();
            assert_eq!(
                r.len(),
                A_COLS,
                "{smi}:基准的原子列数是 {},本文件按 {A_COLS} 列解读 —— \
                 基准过期或列号未同步,重新生成基准(见 harness/README.md)",
                r.len()
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
            cmp("元素", r[A_Z], i64::from(a.atomic_num));
            cmp("电荷", r[A_CHARGE], i64::from(a.formal_charge));
            cmp("同位素", r[A_ISOTOPE], i64::from(a.isotope));
            cmp("显式氢", r[A_EXPL_H], i64::from(a.num_explicit_hs));
            cmp("映射号", r[A_MAP], i64::from(a.atom_map));
            cmp("立体类别", r[A_CHIRAL], a.chiral_tag as i64);
            // 第 8 列(立体排列序号)**刻意不比对**:两边存的不是同一个量。
            //
            // 基准存的是归一到分子键序之后的序号(缺配体还会补位),omgkit 存的
            // 是书写时的字面值 —— 见 `AtomData::stereo_perm`。二者只在"没有
            // 重排且配体齐全"时才碰巧相等,拿来比会得到一堆并非缺陷的分歧。
            // 归一算法属于 L6,做完之后这一列才有比对的意义。
            cmp(
                "芳香",
                r[A_AROMATIC],
                i64::from(a.flags.contains(AtomFlags::AROMATIC)),
            );
            cmp(
                "不推断隐式氢",
                r[A_NO_IMPLICIT],
                i64::from(a.flags.contains(AtomFlags::NO_IMPLICIT)),
            );
        }

        // -- 键列 --
        for (i, row) in rec["bonds"].as_array().unwrap().iter().enumerate() {
            let r: Vec<i64> = row
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_i64().unwrap())
                .collect();
            assert_eq!(
                r.len(),
                B_COLS,
                "{smi}:基准的键列数是 {},本文件按 {B_COLS} 列解读 —— \
                 基准过期或列号未同步,重新生成基准(见 harness/README.md)",
                r.len()
            );
            let b = mol.bonds()[i];
            let mut cmp = |name: &str, exp: i64, got: i64| {
                if exp != got {
                    push(format!("键[{i}].{name}"), exp.to_string(), got.to_string());
                }
            };
            cmp("起点", r[B_BEGIN], i64::from(b.begin));
            cmp("终点", r[B_END], i64::from(b.end));
            cmp("键级", r[B_ORDER], b.order as i64);
            cmp("方向", r[B_DIR], b.direction as i64);
            cmp(
                "芳香",
                r[B_AROMATIC],
                i64::from(b.flags.contains(BondFlags::AROMATIC)),
            );
            // 第 5 列"在环中"刻意不比对 —— 见文件头说明
        }
    }

    (n, bad, cov)
}

/// 分歧报告。按分子聚合,并统计分歧字段的分布 —— 大规模跑时,
/// 分布比逐条列表更能指出问题的性质。
fn report(n: usize, bad: &[Mismatch], limit: usize) -> String {
    let mut by_field: BTreeMap<&str, usize> = BTreeMap::new();
    for m in bad {
        let key = m.field.split('.').next_back().unwrap_or(&m.field);
        *by_field.entry(key).or_default() += 1;
    }

    let mut by_smi: BTreeMap<&str, Vec<&Mismatch>> = BTreeMap::new();
    for m in bad {
        by_smi.entry(&m.smi).or_default().push(m);
    }

    let mut out = format!(
        "\nL1 差分失败:{} 条分子中 {} 条有分歧,共 {} 处\n\n分歧字段分布:\n",
        n,
        by_smi.len(),
        bad.len()
    );
    for (field, count) in &by_field {
        out.push_str(&format!("  {field:<16} {count}\n"));
    }
    out.push_str("\n前若干条:\n");
    for (smi, ms) in by_smi.iter().take(limit) {
        out.push_str(&format!("  {smi}\n"));
        for m in ms.iter().take(6) {
            out.push_str(&format!(
                "      {:<22} 基准={:<10} omgkit={}\n",
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
fn l1_smoke() {
    let (n, bad, cov) = diff_against(&baseline("smoke.l1.jsonl"));
    assert!(n > 0, "基准是空的");
    assert!(bad.is_empty(), "{}", report(n, &bad, 20));

    // 全绿必须是"实现对",不能是"语料没走到"。这几条路径在语料里各有专门的
    // 用例,数为零就说明用例掉了或者被静默跳过了。
    assert!(
        cov.dative_bonds > 0,
        "冒烟语料里一条配位键都没有,该路径是空过的"
    );
    assert!(
        cov.nontetrahedral_atoms > 0,
        "冒烟语料里没有非四面体立体标记,`@SP`/`@TB`/`@OH` 是空过的"
    );
    assert!(
        cov.stereo_perms > 0,
        "冒烟语料里排列序号全为零,该列是空过的"
    );
    assert!(cov.tetrahedral_atoms > 0, "冒烟语料里没有四面体手性");

    println!(
        "L1 冒烟差分通过:{n} 条分子,零分歧;配位键 {},非四面体立体 {},\
         排列序号 {},四面体手性 {}",
        cov.dative_bonds, cov.nontetrahedral_atoms, cov.stereo_perms, cov.tetrahedral_atoms
    );
}

/// 大语料(~8800 条,取自公开的 NCI / ZINC 子集)。
///
/// 语料与基准的生成方式见 `harness/README.md`。
#[test]
#[ignore = "需要先生成大语料基准,见函数文档;用 cargo test -- --ignored 运行"]
fn l1_large() {
    let (n, bad, cov) = diff_against(&baseline("large.l1.jsonl"));
    assert!(n > 1000, "大语料基准看起来不完整:只有 {n} 条");
    assert!(bad.is_empty(), "{}", report(n, &bad, 15));
    println!(
        "L1 大语料差分通过:{n} 条分子,零分歧;配位键 {},非四面体立体 {},\
         四面体手性 {}",
        cov.dative_bonds, cov.nontetrahedral_atoms, cov.tetrahedral_atoms
    );
}
