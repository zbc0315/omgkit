//! 环集的差分测试:相关环集合与参考实现逐环比对。
//!
//! 与 `differential_l2_rings` 比的不是一回事:那个比"键在不在环里""过某原子
//! 的最短环多大",是纯图论量,与选哪组环无关;这里比的是**具体是哪些环**,
//! 是芳香性感知的直接输入。
//!
//! # 判据:按**原子集合**比对,不按顺序
//!
//! 环内原子的排列是遍历产物,环的**原子集合**才是语义量。两边都归一成
//! "排序后的原子下标集合"再比。
//!
//! # 触发面很窄
//!
//! 本实现的环集与"任取一组最小环基"的差别只在少数分子上显现,且恒为多出一个环。
//! 只断言零分歧的话,一个只求最小环基的实现在绝大多数分子上都是对的 ——
//! 所以另有 [`symmetrization_actually_adds_rings`] 确认那些分子确实多出了环。
//!
//! # 已登记的分歧
//!
//! 基准的环集补全在高对称的笼状/带状/螺环体系上并不完全,所以两边在这类
//! 分子上必然不同 —— 本实现给出的环**更多**。[`KNOWN_DIVERGENCES`] 逐条
//! 登记这些分子。
//!
//! 这不是豁免名单:已登记但**已不再分歧**的条目会让测试失败并要求删除。
//!
//! 分歧被 [`ring_set_divergence_does_not_reach_aromaticity`] 限定在环集内 ——
//! 那条测试确认这些分子的芳香输出两边完全一致。

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use omgkit_chem::{
    clean_up, cleanup_organometallics, perceive_rings, ring_set, update_property_cache,
};
use omgkit_io::smiles;

/// 已登记的环集分歧 —— 参考实现的环集补全不完全,笼状体系上两边必然不同。
///
/// 不是豁免名单:测试会检查每个条目是否**仍然**产生分歧,已不再分歧的条目
/// 会让测试失败并要求删除。
const KNOWN_DIVERGENCES: &[(&str, &str)] = &[
    (
        "C12C3C4C1C3C24",
        "K₃,₃ 骨架:6 原子 9 键,圈秩 4。参考实现给 7 个四元环,本实现给 9 个 —— \
         9 个四元面地位完全对等,少给哪两个都没有道理。",
    ),
    (
        "C1C23CC14CC(C2)(C3)C4",
        "三个螺环丁烷首尾相接:参考实现只给 3 个四元环,本实现另给出 6 个六元环。",
    ),
    (
        "C1C2CC3CC4CC1C1CC2CC3CC4C1",
        "[4]cyclacene 分子带:参考实现 4 个六元环,本实现另给 2 个八元环。",
    ),
    (
        "c12c3c4c1c3c24",
        "K₃,₃ 骨架的 sp² 版本(可 kekulize):同上,7 对 9。芳香输出两边一致。",
    ),
];

fn is_known(smi: &str) -> bool {
    KNOWN_DIVERGENCES.iter().any(|(s, _)| *s == smi)
}

struct Mismatch {
    smi: String,
    baseline: Vec<Vec<u32>>,
    omgkit: Vec<Vec<u32>>,
}

struct DiffResult {
    n: usize,
    compared: usize,
    /// 含至少一个环的分子数
    with_rings: usize,
    /// 环总数
    total_rings: usize,
    /// 环数**超过圈秩**的分子数,即本实现比最小环基多出的部分
    beyond_cyclomatic: usize,
    /// 实际产生分歧的已登记分子
    known_hit: std::collections::BTreeSet<String>,
    bad: Vec<Mismatch>,
}

fn baseline(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../harness/baseline")
        .join(name)
}

fn diff_against(path: &Path) -> DiffResult {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| {
        panic!(
            "读不到环集基准 {}: {e}\n生成:python3 harness/oracle_pipeline.py --stage l2 \
             --sanitize-ops CLEANUP,PROPERTIES,SYMMRINGS ...",
            path.display()
        )
    });

    let mut bad: Vec<Mismatch> = Vec::new();
    let (mut n, mut compared, mut with_rings) = (0usize, 0usize, 0usize);
    let (mut total_rings, mut beyond_cyclomatic) = (0usize, 0usize);
    let mut known_hit: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();

    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let rec: serde_json::Value = serde_json::from_str(line).expect("基准 JSONL 格式错误");
        let smi = rec["smi"].as_str().expect("缺 smi").to_string();
        n += 1;

        if !rec["ok"].as_bool().unwrap_or(false) {
            continue;
        }
        let Some(rd_rings) = rec["rings"].as_array() else {
            continue; // 环感知未运行
        };
        let Ok(mut mol) = smiles::parse(&smi) else {
            continue;
        };
        clean_up(&mut mol);
        // 第 2 步必须在价键计算之前 —— 基准的 ops 也含它
        cleanup_organometallics(&mut mol);
        if update_property_cache(&mut mol).is_err() {
            continue;
        }
        let _ = perceive_rings(&mut mol);
        compared += 1;

        let expected: BTreeSet<Vec<u32>> = rd_rings
            .iter()
            .map(|r| {
                let mut v: Vec<u32> = r
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|x| x.as_u64().unwrap() as u32)
                    .collect();
                v.sort_unstable();
                v
            })
            .collect();

        let got: BTreeSet<Vec<u32>> = ring_set(&mol)
            .into_iter()
            .map(|r| {
                let mut v = r.atoms;
                v.sort_unstable();
                v
            })
            .collect();

        if !got.is_empty() {
            with_rings += 1;
            total_rings += got.len();
            // 圈秩 = 键数 - 原子数 + 连通分量数(只算非配位键)
            let n_bonds = mol
                .bonds()
                .iter()
                .filter(|b| b.order != omgkit_core::BondOrder::Dative)
                .count();
            let n_comp = count_components(&mol);
            let cyclomatic = n_bonds + n_comp - mol.num_atoms();
            if got.len() > cyclomatic {
                beyond_cyclomatic += 1;
            }
        }

        if expected != got {
            if is_known(&smi) {
                known_hit.insert(smi.clone());
            } else {
                bad.push(Mismatch {
                    smi: smi.clone(),
                    baseline: expected.into_iter().collect(),
                    omgkit: got.into_iter().collect(),
                });
            }
        }
    }

    DiffResult {
        n,
        compared,
        with_rings,
        total_rings,
        beyond_cyclomatic,
        known_hit,
        bad,
    }
}

/// 连通分量数(把配位键也算作连接 —— 只用于圈秩的粗略统计)
fn count_components(mol: &omgkit_core::MolBuilder) -> usize {
    let n = mol.num_atoms();
    let mut seen = vec![false; n];
    let mut k = 0;
    let mut stack: Vec<u32> = Vec::new();
    for s in 0..n as u32 {
        if seen[s as usize] {
            continue;
        }
        k += 1;
        seen[s as usize] = true;
        stack.push(s);
        while let Some(x) = stack.pop() {
            for (y, bi) in mol.neighbors(x) {
                if mol.bonds()[bi as usize].order == omgkit_core::BondOrder::Dative {
                    continue;
                }
                if !seen[y as usize] {
                    seen[y as usize] = true;
                    stack.push(y);
                }
            }
        }
    }
    k
}

fn report(r: &DiffResult, limit: usize) -> String {
    let mut by_shape: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    for m in &r.bad {
        *by_shape
            .entry((m.baseline.len(), m.omgkit.len()))
            .or_default() += 1;
    }
    let mut out = format!(
        "\n环集差分失败:基准 {} 条,比对 {} 条,{} 条有分歧\n\n\
         (基准环数, 本实现环数) 分布:\n",
        r.n,
        r.compared,
        r.bad.len()
    );
    for ((a, b), c) in &by_shape {
        out.push_str(&format!("  基准 {a:>3} 环 / 本实现 {b:>3} 环   {c} 条\n"));
    }
    out.push_str("\n前若干条:\n");
    for m in r.bad.iter().take(limit) {
        out.push_str(&format!("  {}\n", m.smi));
        let only_rd: Vec<_> = m
            .baseline
            .iter()
            .filter(|r| !m.omgkit.contains(r))
            .collect();
        let only_om: Vec<_> = m
            .omgkit
            .iter()
            .filter(|r| !m.baseline.contains(r))
            .collect();
        if !only_rd.is_empty() {
            out.push_str(&format!("      仅基准有: {only_rd:?}\n"));
        }
        if !only_om.is_empty() {
            out.push_str(&format!("      仅本实现有: {only_om:?}\n"));
        }
    }
    if r.bad.len() > limit {
        out.push_str(&format!("  ...(另有 {} 条)\n", r.bad.len() - limit));
    }
    out
}

#[test]
fn ringset_smoke() {
    let r = diff_against(&baseline("smoke.l2-rings.jsonl"));
    assert!(r.compared > 0, "一条都没比对上");
    assert!(r.bad.is_empty(), "{}", report(&r, 20));
    assert!(r.with_rings > 0, "冒烟语料里一个环都没有,该档是空过的");
    assert_stale_divergences_removed(&r);
    println!(
        "环集冒烟差分通过:比对 {} 条,{} 条含环,共 {} 个环,已登记分歧 {} 条命中",
        r.compared,
        r.with_rings,
        r.total_rings,
        r.known_hit.len()
    );
}

/// 已登记但**已不再分歧**的条目必须删除。
///
/// 一个永远躺着没人动的豁免名单本身就是暗坑;让它在根因消失后主动报警,
/// 才不会沉淀。
fn assert_stale_divergences_removed(r: &DiffResult) {
    let stale: Vec<&str> = KNOWN_DIVERGENCES
        .iter()
        .map(|(s, _)| *s)
        .filter(|s| !r.known_hit.contains(*s))
        .collect();
    assert!(
        stale.is_empty(),
        "以下分子已登记为已知分歧,但本次比对**没有**产生分歧 —— \
         说明根因已消失,请从 KNOWN_DIVERGENCES 中删除这些条目:\n{stale:#?}"
    );
}

/// 把环集分歧**限制在环集内**。
///
/// 已登记的分歧分子上,环集两边不同,但芳香原子数与芳香键数必须完全一致 ——
/// 否则分歧就传导到了下游,性质完全不同。
#[test]
fn ring_set_divergence_does_not_reach_aromaticity() {
    use omgkit_chem::{assign_radicals, kekulize, set_aromaticity};
    use omgkit_core::{AtomFlags, BondFlags};

    // 数值取自参考实现
    for (smi, expect_atoms, expect_bonds) in [
        ("C12C3C4C1C3C24", 0usize, 0usize),
        ("C1C23CC14CC(C2)(C3)C4", 0, 0),
        ("C1C2CC3CC4CC1C1CC2CC3CC4C1", 0, 0),
        ("c12c3c4c1c3c24", 6, 9),
    ] {
        let mut m = smiles::parse(smi).unwrap_or_else(|e| panic!("{smi}: {}", e.render()));
        clean_up(&mut m);
        update_property_cache(&mut m).unwrap_or_else(|e| panic!("{smi}: {e}"));
        let _ = perceive_rings(&mut m);
        kekulize(&mut m).unwrap_or_else(|e| panic!("{smi}: {e}"));
        assign_radicals(&mut m);
        set_aromaticity(&mut m);

        let na = m
            .atoms()
            .iter()
            .filter(|a| a.flags.contains(AtomFlags::AROMATIC))
            .count();
        let nb = m
            .bonds()
            .iter()
            .filter(|b| b.flags.contains(BondFlags::AROMATIC))
            .count();
        assert_eq!(na, expect_atoms, "{smi}: 芳香原子数");
        assert_eq!(nb, expect_bonds, "{smi}: 芳香键数");
    }
}

#[test]
#[ignore = "需要先生成大语料基准;用 cargo test -- --ignored 运行"]
fn ringset_large() {
    let r = diff_against(&baseline("large.l2-rings.jsonl"));
    assert!(r.compared > 1000, "基准不完整:只比对了 {} 条", r.compared);
    assert!(r.bad.is_empty(), "{}", report(&r, 15));
    println!(
        "环集大语料差分通过:比对 {} 条,{} 条含环,共 {} 个环,\
         其中 {} 条的环数超过圈秩",
        r.compared, r.with_rings, r.total_rings, r.beyond_cyclomatic
    );
}

/// 防止"只求一组最小环基也能过"。
///
/// 相关环与最小环基的差别只在很少的分子上显现,只断言"零分歧"的话,一个只求
/// 最小环基的实现在绝大多数分子上都是对的 —— 必须显式确认那些分子**确实**
/// 多出了环。
#[test]
#[ignore = "需要先生成大语料基准;用 cargo test -- --ignored 运行"]
fn symmetrization_actually_adds_rings() {
    let r = diff_against(&baseline("large.l2-rings.jsonl"));
    assert_eq!(
        r.beyond_cyclomatic, 49,
        "应有 49 条分子的环数超过圈秩。\
         数字变化意味着实现或语料发生了变动,请先查清原因再更新此断言。"
    );
}
