//! 子结构匹配的差分测试:「分子语料 × SMARTS 语料」的全部命中。
//!
//! 基准由 `harness/oracle_matches.py` 生成,格式与
//! `cargo run -p omgkit-match --example dump_matches` 的输出一致。
//!
//! # 比的是**匹配原子集合**,不是匹配数
//!
//! 只比数量的话,"命中了别的原子但个数凑巧相同"发现不了。集合逐条比对才能
//! 钉住"命中的是哪些原子"。
//!
//! # 序号错位是这类测试最阴的失效方式
//!
//! 基准里存的是**行号**,两侧对语料的行过滤规则必须完全一致。少跳一类行
//! (例如 `#` 注释行)就会让序号整体偏移,报出来的每一条分歧都指向错误的
//! 分子 —— 内容看着像化学问题,根因却在读文件那一步。
//!
//! 所以这里**不重新读语料**:基准里已经有分子和模式的序号,测试自己按同一套
//! 规则读一次,并断言条数与基准里出现过的最大序号相容。

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use omgkit_chem::sanitize;
use omgkit_io::{smarts, smiles};
use omgkit_match::{substructure_matches, MatchOptions, MolProps};

/// 与基准生成器保持一致。两边不同的话,高度对称的分子上会出现
/// "只是截断点不同"的假分歧。
const MAX_MATCHES: usize = 1000;

/// **基准里有命中、而本实现净化不了**的分子。
///
/// 语料里本来就有一批净化不了的分子(超价、不可 kekulize),那些两边都失败,
/// 基准里也没有它们的命中,无需登记。只有"基准能匹配、我却连净化都过不去"
/// 才是真缺口 —— 名单登记的就是这一类。
///
/// **不是豁免名单**:名单里的分子若能净化了,测试会报错要求删除条目。
///
/// 当前为空。
const KNOWN_UNSANITIZABLE: &[&str] = &[];

fn corpus(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../harness/corpus")
        .join(name)
}

fn baseline(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../harness/baseline")
        .join(name)
}

/// 语料的行过滤规则。**必须与 `harness/oracle_matches.py` 的 `read_lines`
/// 一字不差** —— 基准里存的是行号。
fn read_corpus(path: &Path) -> Vec<String> {
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("读不到语料 {}: {e}", path.display()));
    text.lines()
        .filter_map(|l| l.split_whitespace().next())
        .filter(|t| !t.starts_with('#'))
        .map(String::from)
        .collect()
}

/// 基准里的一行:(分子序号, 模式序号) → 命中的原子集合(已归一)
/// 基准的内容,连同它**覆盖了多少个分子**。
///
/// 覆盖范围随数据走。基准可能只跑了前 N 个分子,而这里读的是整份语料 ——
/// 两个数对不上时,第 N 个之后的分子全变成"只有本实现有命中",凭空多出成千
/// 上万条假分歧。凡是 dump 侧可以只覆盖一部分输入的判据都要这么防。
struct Baseline {
    hits: BTreeMap<(usize, usize), String>,
    n_mols: usize,
}

fn read_baseline(path: &Path) -> Baseline {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| {
        panic!(
            "读不到匹配基准 {}: {e}\n生成方式见 harness/README.md",
            path.display()
        )
    });
    let mut n_mols = None;
    let mut hits = BTreeMap::new();
    for l in text.lines().filter(|l| !l.trim().is_empty()) {
        if let Some(rest) = l.strip_prefix("#mols\t") {
            n_mols = Some(rest.parse().expect("#mols 首行不是数字"));
            continue;
        }
        let mut it = l.split('\t');
        let mi: usize = it.next().expect("缺分子序号").parse().expect("分子序号");
        let pi: usize = it.next().expect("缺模式序号").parse().expect("模式序号");
        let sets = it.next().expect("缺命中集合").to_string();
        hits.insert((mi, pi), sets);
    }
    let n_mols = n_mols.unwrap_or_else(|| {
        panic!(
            "基准 {} 缺 `#mols<TAB>N` 首行,无从知道它覆盖了多少个分子。\n\
             用当前版本的 harness/oracle_matches.py 重新生成",
            path.display()
        )
    });
    Baseline { hits, n_mols }
}

/// 把一次匹配的结果归一成基准的字符串形式。
fn encode(hits: &[Vec<u32>]) -> String {
    let mut sets: Vec<String> = hits
        .iter()
        .map(|h| {
            let mut v = h.clone();
            v.sort_unstable();
            v.iter().map(u32::to_string).collect::<Vec<_>>().join(",")
        })
        .collect();
    sets.sort();
    sets.join("|")
}

struct Report {
    only_baseline: Vec<(usize, usize)>,
    only_ours: Vec<(usize, usize)>,
    different: Vec<(usize, usize)>,
    /// 两边都有且完全一致的条数
    agreed: usize,
    /// 因净化失败而跳过的分子
    skipped_molecules: BTreeSet<String>,
}

/// 返回 (报告, 分子串, 模式串, **基准覆盖到的分子数**, **语料总条数**)。
///
/// 后两个是**分母**:基准可以只覆盖语料的一段,而那一段有多大必须报出来 ——
/// 见 `matches_large` 的说明(它读的基准只覆盖 8839 条里的前 2000 条)。
fn run(mol_corpus: &str, baseline_name: &str) -> (Report, Vec<String>, Vec<String>, usize, usize) {
    let pats_raw = read_corpus(&corpus("smarts.txt"));
    let base = read_baseline(&baseline(baseline_name));
    let expect = base.hits;
    // 只比基准覆盖到的那一段
    let mut smis = read_corpus(&corpus(mol_corpus));
    let total_mols = smis.len();
    assert!(
        smis.len() >= base.n_mols,
        "基准覆盖 {} 个分子,而语料只有 {} 个 —— 语料对不上",
        base.n_mols,
        smis.len()
    );
    smis.truncate(base.n_mols);

    let queries: Vec<Option<smarts::QueryMol>> =
        pats_raw.iter().map(|s| smarts::parse(s).ok()).collect();

    let opts = MatchOptions {
        max_matches: MAX_MATCHES,
        uniquify: true,
        use_chirality: true,
    };
    let mut got: BTreeMap<(usize, usize), String> = BTreeMap::new();
    let mut skipped = BTreeSet::new();

    for (mi, smi) in smis.iter().enumerate() {
        let Ok(mut mol) = smiles::parse(smi) else {
            continue;
        };
        if sanitize(&mut mol).is_err() {
            skipped.insert(smi.clone());
            continue;
        }
        let props = MolProps::compute(&mol);
        for (pi, q) in queries.iter().enumerate() {
            let Some(q) = q else { continue };
            let hits = substructure_matches(q, &mol, &props, opts);
            if !hits.is_empty() {
                got.insert((mi, pi), encode(&hits));
            }
        }
    }

    // 净化失败**且基准里有命中**的分子才需要登记 —— 两边都失败的无需登记
    let skipped_idx: BTreeSet<usize> = smis
        .iter()
        .enumerate()
        .filter(|(_, s)| skipped.contains(*s))
        .map(|(i, _)| i)
        .collect();
    let lost: BTreeSet<&str> = expect
        .keys()
        .filter(|(mi, _)| skipped_idx.contains(mi))
        .map(|(mi, _)| smis[*mi].as_str())
        .collect();
    let unregistered: Vec<&&str> = lost
        .iter()
        .filter(|s| !KNOWN_UNSANITIZABLE.contains(*s))
        .collect();
    assert!(
        unregistered.is_empty(),
        "以下分子基准里有命中,本实现却净化不了,且未登记:{unregistered:#?}"
    );
    let stale: Vec<&&str> = KNOWN_UNSANITIZABLE
        .iter()
        .filter(|k| smis.iter().any(|s| s == *k) && !lost.contains(*k))
        .collect();
    assert!(
        stale.is_empty(),
        "以下已登记的分子现在不再丢失匹配了,请从 KNOWN_UNSANITIZABLE 删除:{stale:?}"
    );
    let keep = |k: &(usize, usize)| !skipped_idx.contains(&k.0);

    let mut rep = Report {
        only_baseline: Vec::new(),
        only_ours: Vec::new(),
        different: Vec::new(),
        agreed: 0,
        skipped_molecules: skipped,
    };
    for (k, v) in &expect {
        if !keep(k) {
            continue;
        }
        match got.get(k) {
            None => rep.only_baseline.push(*k),
            Some(g) if g != v => rep.different.push(*k),
            Some(_) => rep.agreed += 1,
        }
    }
    for k in got.keys() {
        if keep(k) && !expect.contains_key(k) {
            rep.only_ours.push(*k);
        }
    }
    (rep, smis, pats_raw, base.n_mols, total_mols)
}

fn describe(rep: &Report, smis: &[String], pats: &[String], limit: usize) -> String {
    let mut out = format!(
        "\n匹配差分失败:一致 {},只有基准 {},只有本实现 {},集合不同 {}\n\n",
        rep.agreed,
        rep.only_baseline.len(),
        rep.only_ours.len(),
        rep.different.len()
    );
    let show = |out: &mut String, tag: &str, ks: &[(usize, usize)]| {
        for &(mi, pi) in ks.iter().take(limit) {
            out.push_str(&format!(
                "  [{tag}] 分子#{mi} {}\n          模式#{pi} {}\n",
                &smis[mi][..smis[mi].len().min(90)],
                &pats[pi][..pats[pi].len().min(90)]
            ));
        }
    };
    show(&mut out, "只有基准", &rep.only_baseline);
    show(&mut out, "只有本实现", &rep.only_ours);
    show(&mut out, "集合不同", &rep.different);
    out
}

fn assert_clean(rep: &Report, smis: &[String], pats: &[String]) {
    assert!(
        rep.only_baseline.is_empty() && rep.only_ours.is_empty() && rep.different.is_empty(),
        "{}",
        describe(rep, smis, pats, 8)
    );
}

#[test]
fn matches_smoke() {
    let (rep, smis, pats, base_mols, total_mols) = run("smoke.smi", "smoke.matches.tsv");
    assert_clean(&rep, &smis, &pats);
    assert!(rep.agreed > 20, "只对上了 {} 条,语料太弱", rep.agreed);
    println!(
        "匹配差分(冒烟):基准覆盖语料前 {base_mols} 条(语料共 {total_mols} 条),\
         一致 {} 条;跳过(已登记){} 条分子",
        rep.agreed,
        rep.skipped_molecules.len()
    );
}

#[test]
// **不再 `#[ignore]`。** 它读的 `harness/baseline/matches.tsv` 先前没有入库,
// 所以这一条连本地 `--ignored` 都要先手工生成基准 —— 现在那份基准跟着入库了
// (164 651 字节),理由见 `.gitignore`:默认 `cargo test` 会跑到的基准必须入库。
//
// # 名字叫 `large`,实际只覆盖前 2000 条 —— 这一点必须写在脸上
//
// `matches.tsv` 首行是 `#mols 2000`,而 `large.smi` 有 8839 条 —— **覆盖 22.6%**。
// 独立审核指出:测试名叫 `matches_large`、读的是 `large.smi`、打印"大语料",
// 而 2000 这个数**一处都没露过面**。所以下面把 `base.n_mols` 打出来。
//
// 审核还用锁里钉的 RDKit 2025.09.2 重导了**全量**基准(8839 条,678 kB),
// 报出 6 条本实现命中而 RDKit 不命中的方向键模式(都在同一个分子上,小环内的
// C=C 被我们给了方向语义)。**那一条我没能自己复核**(复现要重导全量基准),
// 已按原样记进任务,不当作已确证的结论。
//
// 也就是说:**这不只是"覆盖率数字不好看",截断很可能挡住了活的分歧。**
// 扩到全量是另一件事(要连着处理那 6 条),见任务清单。
fn matches_large() {
    let (rep, smis, pats, base_mols, total_mols) = run("large.smi", "matches.tsv");
    assert_clean(&rep, &smis, &pats);
    assert!(rep.agreed > 10_000, "只对上了 {} 条", rep.agreed);
    // **跳过的分子要有上限。** 先前只印不判 —— 那是这个仓库刚清过一批的
    // "单向过滤器":它只会把分歧变成不计数。实测现值 0。
    assert!(
        rep.skipped_molecules.len() <= 5,
        "跳过了 {} 个分子(上限 5)—— 跳过的越多,上面那个'一致 N 条'越不说明问题",
        rep.skipped_molecules.len()
    );
    println!(
        "匹配差分:基准覆盖语料前 {} 条(语料共 {} 条),一致 {} 条;跳过(已登记){} 条分子",
        base_mols,
        total_mols,
        rep.agreed,
        rep.skipped_molecules.len()
    );
}
