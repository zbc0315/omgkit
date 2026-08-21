//! SMARTS 写出必须**保住语义**:写出去再读回来,匹配到的东西要一模一样。
//!
//! # 为什么幂等 + 规模守恒不够
//!
//! `omgkit-io/tests/roundtrip_smarts.rs` 判的是"写出幂等 + 原子数键数不变",
//! 而它的文档自己写着"语义等价要外部实现当判官,那是另一档"。
//! 问题是那一档的判官(`harness/check_smarts_write.py`)**不在 CI 里**,
//! 于是这一档实际上没人守 —— 而幂等这条**挡不住把谓词写成常量**:
//!
//! ```text
//! 把 write.rs 里 `AtomPrim::Charge(c) => charge_string(*c)` 改成
//! `charge_string(0)`,于是:
//!
//!     [N+]   -> [N&+0]      写出幂等 ✓  规模守恒 ✓
//!     [O-]   -> [O&+0]      阴离子氧写成了中性氧
//!     [NX3+] -> [N&X3&+0]
//!
//! `cargo test --release` 与 `cargo test --workspace` **全绿**(实测)。
//! ```
//!
//! **776 条语料 SMARTS 里 700 条带形式电荷。**
//!
//! # 这一档不需要外部实现
//!
//! "我们的 SMARTS 语义对不对"确实要外部判官;但"**我们的写出器保不保得住
//! 我们自己的语义**"不用 —— 拿同一个匹配器跑两遍就够了:写丢了电荷,
//! `parse(write(q))` 在我们自己的匹配器下就会匹配到不同的原子。
//!
//! # 判据自己也要有分母
//!
//! "两边匹配集合相同"对**两个都匹配不到东西**的查询同样成立。所以还要断言
//! 真正命中过的模式数与总命中次数(`MIN_PATTERNS_HIT` / `MIN_HITS`)。
//!
//! 这条闸当场就发挥了作用:我头一版拿**冒烟语料**当探针,而它前几十条是乙醇、
//! 甲烷、环己烷这类小分子 —— 756 条 SMARTS 总共只命中 **69** 次。判据几乎是空的,
//! 而语义比对那部分照样"全过"。**分母闸不只防回归,也防写判据的人把它写空。**

use std::collections::BTreeSet;
use std::path::PathBuf;

use omgkit_chem::sanitize;
use omgkit_io::{smarts, smiles};
use omgkit_match::{substructure_matches, MatchOptions, MolProps};

/// 探针分子:`large.smi` 里前 N 个**解析且净化得了**的。**不用全部** ——
/// 756 条 SMARTS(语料 795 行,其余是注释/空行/解析不了的)× 全部分子是平方级,
/// 而判据要的是"覆盖到各类谓词",不是"跑遍语料"。
///
/// **用药物样语料,不用冒烟语料。** 冒烟档前几十条是乙醇、甲烷、环己烷这类
/// 小分子,756 条 SMARTS 在它们上面总共只命中 **69** 次 —— 判据几乎是空的。
/// 这一条是被下面那个 `MIN_HITS` 当场抓出来的:我头一版按"应该上万"拍了个
/// 下限,跑出来 69,于是才发现探针选错了。**分母闸不只防回归,也防我自己
/// 把判据写空。**
const N_PROBES: usize = 2000;

/// 至少要有多少条 SMARTS **命中过**,以及总命中多少次。见模块文档"判据也要有分母"。
///
/// 逐档实测(`large.smi` 前 N 条作探针,756 条 SMARTS):
///
/// | 探针分子 | 命中过的模式 | 总命中 | 测试档耗时 |
/// |---|---|---|---|
/// | 40 | 71 | 435 | 0.2 s |
/// | 200 | 108 | 1 675 | 0.6 s |
/// | 600 | 163 | 4 928 | 1.8 s |
/// | **2000** | **212** | **15 092** | **6.4 s** |
///
/// 取 2000:再往上加分子的边际覆盖很低(600→2000 只多 49 条模式),
/// 而 CI 要跑两遍(release 与测试档)。加上手工补的那 17 条谓词之后,
/// 现值是 **773 条里 226 条命中过、总命中 58 863 次**。
///
/// **如实记下来:773 条里只有 226 条被真正验证过**(29%)。其余那些在这批
/// 药物样分子上一次都匹配不到 —— 它们的写出结果没有被验证。要提上去只能换
/// 更杂的探针语料,那是另一件事。
///
/// 下限设在实测值的 76~84%,不贴着现值 —— 语料微调不该让判据变红。
/// (`MIN_PATTERNS_HIT` 190/226 = 84%,`MIN_HITS` 45 000/58 863 = 76%。)
const MIN_PATTERNS_HIT: usize = 190;
/// 见 [`MIN_PATTERNS_HIT`]。
const MIN_HITS: usize = 45_000;

/// 每条查询在每个分子上最多找几个匹配。**只为封顶耗时** ——
/// 两边用同一个上限,所以"集合相同"这条判据不受它影响。
const MAX_MATCHES: usize = 8;

fn corpus(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../harness/corpus")
        .join(name)
}

fn probes() -> Vec<(String, omgkit_core::MolBuilder, MolProps)> {
    let text = std::fs::read_to_string(corpus("large.smi")).expect("读得到探针语料");
    let mut out = Vec::new();
    for line in text.lines() {
        if out.len() >= N_PROBES {
            break;
        }
        let Some(smi) = line.split_whitespace().next() else {
            continue;
        };
        if smi.is_empty() || smi.starts_with('#') {
            continue;
        }
        let Ok(mut m) = smiles::parse(smi) else {
            continue;
        };
        if sanitize(&mut m).is_err() {
            continue;
        }
        let props = MolProps::compute(&m);
        out.push((smi.to_string(), m, props));
    }
    out
}

/// 一条查询在一个分子上匹配到的**原子集合的集合**。
///
/// **不能比映射元组。** `smarts::write` 会重排查询原子的编号(生成树的遍历
/// 顺序与原编号无关,`roundtrip_smarts` 的文档就写着这一点),于是回读之后
/// 同一个匹配的元组是**置换过的**。头一版比元组,当场报出两条"语义变了",
/// 逐个查下去两边匹配到的是同一批分子原子、只是次序不同 —— 那是判据自己的
/// 假红,不是写出器的缺陷。
///
/// 比"匹配到哪些原子"就与编号无关了。代价是看不见"匹配到同一批原子但对应
/// 关系变了"这一档。
///
/// (先前这里写着"那一档由 `roundtrip_smarts` 的规模守恒 + 幂等兜着" ——
///  **假的**。独立审核实测:把写出的原子映射号一律换成 `*`,规模守恒、幂等,
///  `roundtrip_smarts` 与本判据**双双全绿**。那一档现在无人守,已立任务。)
fn hits(
    q: &omgkit_io::smarts::QueryMol,
    m: &omgkit_core::MolBuilder,
    p: &MolProps,
) -> BTreeSet<BTreeSet<u32>> {
    substructure_matches(
        q,
        m,
        p,
        MatchOptions {
            max_matches: MAX_MATCHES,
            uniquify: true,
            use_chirality: true,
        },
    )
    .into_iter()
    .map(|mp| mp.into_iter().collect::<BTreeSet<u32>>())
    .collect()
}

#[test]
fn smarts_写出去再读回来匹配到的东西不变() {
    let probes = probes();
    // `probes()` 是"一直往后读到凑满 `N_PROBES` 个",所以正常情况下正好拿满;
    // 拿不满只可能是语料没喂进来或整批净化不了。
    assert!(
        probes.len() >= N_PROBES * 9 / 10,
        "只准备出 {} 个探针分子(要 {N_PROBES} 个),判据的覆盖面不够",
        probes.len()
    );

    let corpus_text = std::fs::read_to_string(corpus("smarts.txt")).expect("读得到 SMARTS 语料");
    // **语料盖不到的谓词要手工补上,而"盖不盖得到"要用变异量,不能 grep。**
    //
    // 我头一版是 `grep -cE` 数的,数出来"电荷 707 条、`R<n>` 0 条" —— 两个都错:
    // `-` 在 SMARTS 里还是单键符号,`R` 还是环键基元。**光看字符串数不出谓词。**
    //
    // 换成变异量:把写出器里某一档谓词写成常量,看判据红不红。
    // (`large.smi` 前 2000 条作探针)
    //
    // | 写出器变异 | 只用语料 | 语料 + 下面这 17 条 |
    // |---|---|---|
    // | 电荷写成 `+0` | **红** | 红 |
    // | `r<n>` 写成 `r` | 绿 | **红** |
    // | `R<n>` 写成 `R` | 绿 | **红** |
    // | `v<n>` 写成 `v4` | 绿 | **红** |
    // | `x<n>` 写成 `x` | 绿 | **红** |
    // | `h<n>` 写成 `h0` | 绿 | **红** |
    // | 同位素丢掉 | 绿 | **红** |
    // | 原子映射号丢掉 | 绿 | **红**(但见下) |
    //
    // 八档里七档靠这 17 条撑着 —— **判据的覆盖面被语料的谓词覆盖面卡死**,
    // 而这一点光看判据代码看不出来。
    //
    // **原子映射号那一档要打个折扣**:映射号不参与匹配,所以这条判据看不见
    // "映射号写错了",只看得见"写出的串坏到解析不了"。真要守它得比
    // `write_reaction` 两边的映射对应关系 —— 那一档目前无人守,已立任务。
    let handmade = "\
[C;r5]\n[C;r6]\n[C;R1]\n[C;R2]\n[c;x2]\n[c;x3]\n[C;v4]\n[N;v3]\n\
[C;h0]\n[C;h1]\n[13C]\n[C:1]\n[C@H](F)(Cl)Br\n[C@@H](F)(Cl)Br\n\
[C;R1;r6]\n[N;X3;v3]\n[c;x2;R1]\n";
    let text = format!("{corpus_text}\n{handmade}");
    let mut n_pat = 0usize;
    let mut n_hits = 0usize;
    let mut n_pat_hit = 0usize;
    let mut bad: Vec<String> = Vec::new();

    for line in text.lines() {
        let src = line.trim();
        if src.is_empty() || src.starts_with('#') {
            continue;
        }
        let Ok(q1) = smarts::parse(src) else { continue };
        let w = smarts::write(&q1);
        let q2 = match smarts::parse(&w) {
            Ok(q) => q,
            Err(e) => {
                bad.push(format!("  {src}\n    写出 {w} 之后解析失败:{}", e.render()));
                continue;
            }
        };
        n_pat += 1;
        let mut this_hits = 0usize;
        for (smi, m, p) in &probes {
            let (a, b) = (hits(&q1, m, p), hits(&q2, m, p));
            n_hits += a.len();
            this_hits += a.len();
            if a != b {
                if bad.len() < 10 {
                    bad.push(format!(
                        "  {src}\n    写出 {w}\n    在 {smi} 上:原查询 {} 个匹配,回读 {} 个",
                        a.len(),
                        b.len()
                    ));
                }
                break;
            }
        }
        if this_hits > 0 {
            n_pat_hit += 1;
        }
    }
    println!("比了 {n_pat} 条 SMARTS,其中 {n_pat_hit} 条至少命中一次;总命中 {n_hits} 次");

    // 列表封顶在 10 条,所以**报的是"至少这么多"** —— 别把它当成总数。
    assert!(
        bad.is_empty(),
        "至少 {} 条 SMARTS 写出去之后匹配到的东西变了(只列前 {}，共比了 {n_pat} 条):\n{}",
        bad.len(),
        bad.len().min(10),
        bad.join("\n")
    );
    // **分母。** 两边都匹配不到东西时"集合相同"恒成立 —— 见模块文档。
    // 真正该量的是**多少条 SMARTS 至少命中过一次**:总命中次数会被少数几条
    // 高频模式撑起来,而那说明不了别的模式被测到了。
    assert!(
        n_pat_hit >= MIN_PATTERNS_HIT,
        "只有 {n_pat_hit} 条 SMARTS 命中过(下限 {MIN_PATTERNS_HIT},共比了 {n_pat} 条)\
         —— 其余那些的写出结果没有被真正验证过"
    );
    assert!(
        n_hits >= MIN_HITS,
        "总命中只有 {n_hits} 次(下限 {MIN_HITS})—— 判据可能什么也没测到"
    );
    assert!(n_pat >= 700, "只比了 {n_pat} 条 SMARTS,语料没喂进来");
}
