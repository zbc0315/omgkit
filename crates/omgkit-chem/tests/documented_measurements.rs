//! 复核源码注释里那些**可测量的**定量声明。
//!
//! # 为什么要有这一档
//!
//! 注释里写"实测 N 条"的地方,时间一长很容易与实现脱节 —— 尤其是被**别处的
//! 改动间接影响**的那些。它们不会报错,只会静静地变成假话,然后把后来的人
//! 引向错误的排查方向。
//!
//! 典型形状:canon.rs 那句"打破对称会影响几条分子"取决于写出器写出多少信息,
//! 而写出器在 write.rs。改一边不会让另一边变红,没有任何东西把两者连起来。
//!
//! # 口径必须与被验证的那句话一致
//!
//! 复核方与被复核的那句话口径不同时,报出来的"对不上"是假的:
//!
//! | 声明 | 口径要求 |
//! |---|---|
//! | 含环分子数 | 与 `differential_l2_ringset` 同口径,只数净化得了的那些;数全部 8839 条会多出 6 条 |
//! | 打破对称影响数 | 必须走 `canonical_smiles` 的立体预处理,否则量的不是同一个分子 |
//!
//! **复核工具本身量错口径,比不复核更糟** —— 它会让人去改本来正确的注释。
//! 所以每条断言都把口径写在断言旁边。
//!
//! # 数字变了怎么办
//!
//! 先判断是实现退步还是**能力增长**(写出得越细,能区分起点的分子就越多,
//! 打破对称那个数只增不减)。确认之后**同时**改注释和本文件的期望值 ——
//! 只改一处就等于把这道闸门关掉了。

use std::path::PathBuf;

use omgkit_chem::{clean_up, cleanup_organometallics, perceive_rings, sanitize};
use omgkit_io::{canon, smiles};

/// 语料总条数。差分测试比对的是其中**两侧都净化得了**的 8831 条。
const CORPUS_TOTAL: usize = 8839;

/// `organometallics.rs` 与 `builder.rs` 声称:第 2 步在语料上只改动 2 条分子。
const CLAIM_ORGANOMETALLIC_MOLECULES: usize = 2;

/// `sssr.rs` 声称:7553 条含环分子。口径同 `differential_l2_ringset`。
const CLAIM_MOLECULES_WITH_RINGS: usize = 7553;

/// `canon.rs` 声称:7 条分子取不同起点会写出不同的串。
const CLAIM_TIE_BREAK_MATTERS: usize = 7;

fn corpus() -> Vec<String> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../harness/corpus/large.smi");
    std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("读不到语料 {}: {e}", path.display()))
        .lines()
        .filter_map(|l| l.split_whitespace().next())
        .filter(|t| !t.starts_with('#'))
        .map(String::from)
        .collect()
}

fn expect(name: &str, got: usize, want: usize, where_: &str) {
    assert_eq!(
        got, want,
        "\n{name}:实测 {got},而 {where_} 的注释写着 {want}。\n\
         先判断是实现退步还是能力增长,再**同时**改注释和本文件的期望值 —— \
         只改一处等于把这道闸门关掉。"
    );
}

#[test]
#[ignore = "要跑全量语料;用 cargo test -- --ignored 运行"]
fn documented_numbers_still_hold() {
    let smis = corpus();
    assert_eq!(smis.len(), CORPUS_TOTAL, "语料条数变了,下面的数都要重量");

    // 一、第 2 步的触发面
    let mut organo = 0usize;
    for s in &smis {
        let Ok(mut m) = smiles::parse(s) else {
            continue;
        };
        clean_up(&mut m);
        if cleanup_organometallics(&mut m) > 0 {
            organo += 1;
        }
    }
    expect(
        "第 2 步改动的分子数",
        organo,
        CLAIM_ORGANOMETALLIC_MOLECULES,
        "organometallics.rs / builder.rs",
    );

    // 二、含环分子数
    //
    // 口径:只数**净化得了**的那些 —— 与 differential_l2_ringset 一致。
    // 少了这道过滤会数出 7559,然后被误读成"注释过期"。
    let mut with_rings = 0usize;
    for s in &smis {
        let Ok(m) = smiles::parse(s) else { continue };
        let mut probe = m.clone();
        if sanitize(&mut probe).is_err() {
            continue;
        }
        let mut m = m;
        clean_up(&mut m);
        cleanup_organometallics(&mut m);
        if perceive_rings(&mut m).bond_in_ring.iter().any(|&x| x) {
            with_rings += 1;
        }
    }
    expect(
        "含环分子数",
        with_rings,
        CLAIM_MOLECULES_WITH_RINGS,
        "sssr.rs",
    );

    // 三、打破对称真正影响结果的分子数
    //
    // 口径:`tie_break_matters` 与 `canonical_smiles` 共用同一条立体预处理,
    // 否则量的不是同一个分子。
    let mut tie_break = 0usize;
    for s in &smis {
        let Ok(mut m) = smiles::parse(s) else {
            continue;
        };
        if sanitize(&mut m).is_err() {
            continue;
        }
        if canon::tie_break_matters(&m) {
            tie_break += 1;
        }
    }
    expect(
        "打破对称真正影响结果的分子数",
        tie_break,
        CLAIM_TIE_BREAK_MATTERS,
        "canon.rs",
    );

    println!(
        "注释里的定量声明全部复核通过:第 2 步 {organo} 条、含环 {with_rings} 条、\
         打破对称 {tie_break} 条"
    );
}
