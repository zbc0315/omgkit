//! SMARTS 解析的差分测试。
//!
//! 基准由 `harness/oracle_smarts.py` 生成,语料是真实项目里在用的 SMARTS
//! (PAINS 过滤器、官能团层级、RLewis 库),不是手挑的 —— 手挑会不自觉地
//! 只挑自己已经实现的写法。
//!
//! # 比三件事,缺一不可
//!
//! | 量 | 少了它会怎样 |
//! |---|---|
//! | **可解析性** | 静默拒绝一批合法写法,或静默接受非法写法 |
//! | **原子数** | 表达式边界切错(比如把 `=!@` 当成两条键) |
//! | **键数** | 环闭合配错,或分支处理错 |
//!
//! 只比可解析性最危险:一个把所有内容都当通配符的实现能 100% 通过。

use std::path::PathBuf;

use omgkit_io::smarts;

fn baseline() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../harness/baseline/smarts.jsonl")
}

struct Mismatch {
    smarts: String,
    what: String,
}

#[derive(Default)]
struct Stats {
    n: usize,
    /// 两边都能解析
    both_ok: usize,
    /// 两边都拒绝
    both_err: usize,
    /// 用到递归 SMARTS `$(...)` 的条数
    with_recursive: usize,
    /// 用到键表达式逻辑运算的条数
    with_bond_logic: usize,
}

fn diff() -> (Stats, Vec<Mismatch>) {
    let path = baseline();
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| {
        panic!(
            "读不到 SMARTS 基准 {}: {e}\n生成方式见 harness/README.md",
            path.display()
        )
    });

    let mut stats = Stats::default();
    let mut bad = Vec::new();

    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let rec: serde_json::Value = serde_json::from_str(line).expect("基准 JSONL 格式错误");
        let s = rec["smarts"].as_str().expect("缺 smarts").to_string();
        stats.n += 1;
        if s.contains("$(") {
            stats.with_recursive += 1;
        }
        if s.contains("!@") || s.contains(",:") || s.contains("-,") {
            stats.with_bond_logic += 1;
        }

        let rd_ok = rec["ok"].as_bool().unwrap_or(false);
        let ours = smarts::parse(&s);

        match (rd_ok, &ours) {
            (false, Ok(_)) => bad.push(Mismatch {
                smarts: s,
                what: "基准拒绝,本实现接受".into(),
            }),
            (true, Err(e)) => bad.push(Mismatch {
                smarts: s,
                what: format!("基准接受,本实现拒绝:{}", e.kind),
            }),
            (false, Err(_)) => stats.both_err += 1,
            (true, Ok(q)) => {
                stats.both_ok += 1;
                let (na, nb) = (
                    rec["na"].as_u64().unwrap_or(0) as usize,
                    rec["nb"].as_u64().unwrap_or(0) as usize,
                );
                if q.num_atoms() != na || q.num_bonds() != nb {
                    bad.push(Mismatch {
                        smarts: s,
                        what: format!(
                            "规模不同:基准 {na} 原子 {nb} 键,本实现 {} 原子 {} 键",
                            q.num_atoms(),
                            q.num_bonds()
                        ),
                    });
                } else if !q.is_consistent() {
                    bad.push(Mismatch {
                        smarts: s,
                        what: "查询树与拓扑长度不一致".into(),
                    });
                }
            }
        }
    }
    (stats, bad)
}

fn report(bad: &[Mismatch], limit: usize) -> String {
    let mut out = format!("\nSMARTS 差分失败 {} 条:\n\n", bad.len());
    for m in bad.iter().take(limit) {
        out.push_str(&format!("  {}\n      {}\n", m.smarts, m.what));
    }
    if bad.len() > limit {
        out.push_str(&format!("  ...(另有 {} 条)\n", bad.len() - limit));
    }
    out
}

#[test]
fn smarts_parse_matches_baseline() {
    let (stats, bad) = diff();
    assert!(bad.is_empty(), "{}", report(&bad, 15));
    assert!(stats.n > 500, "语料只有 {} 条", stats.n);

    // 全绿必须是"实现对",不能是"语料没走到"
    assert!(
        stats.with_recursive > 0,
        "语料里没有递归 SMARTS,`$(...)` 是空过的"
    );
    assert!(
        stats.with_bond_logic > 0,
        "语料里没有键表达式的逻辑运算,那条路径是空过的"
    );
    assert!(
        stats.both_err > 0,
        "语料里没有双方都拒绝的用例 —— 拒绝行为没有守护,一个来者不拒的实现照样能全绿"
    );

    println!(
        "SMARTS 差分通过:{} 条,两边都成功 {},都拒绝 {};含递归 {},含键逻辑 {}",
        stats.n, stats.both_ok, stats.both_err, stats.with_recursive, stats.with_bond_logic
    );
}
