//! 入库基准的规模是**契约**,不是随手生成的副产品。
//!
//! # 判官不可能自己发现基准被截短了
//!
//! `harness/baseline/` 下的判官基准都是逐行一个分子。判官读多少行就比多少个
//! 分子,而**它无从知道自己本该收到多少行** —— 于是基准一旦变小,每一条判据
//! 都会在更少的分子上算出更好看的数,而且全部退 0。实测(把基准截成 1 行):
//!
//! | 判官 | 27 行 | 截成 1 行 |
//! |---|---|---|
//! | `smooth_oracle` | 过 | **过** |
//! | `eigen_oracle` | 过 | **过** |
//! | `threading_oracle` | 过 | **过** |
//! | `conformer_oracle`(150 → 1) | 过 | **过** |
//!
//! 判官内部的分母闸(`MAX_UNREAD` / `MAX_FAIL` / 两个基准行数相等)挡的是
//! "读进来了却没比到";这一条挡的是另一半:**读进来的本身就少了**。
//!
//! # 为什么用一条测试而不是判官的参数
//!
//! 让每个判官收一个"至少该有几行"的参数,那个数就要在 `harness/gates.sh` 与
//! `.github/workflows/ci.yml` 里各手抄一遍 —— 本仓库为手抄的数栽过好几次
//! (`gates.sh` 的 `N/14`、`ci.yml` 的"四道闸门"、`verify_stereo.py` 的
//! "第 14 道闸")。写成一条测试,数只有一份,而且改基准时会**当场红**,
//! 逼着改的人确认那是有意的。
//!
//! # 改了基准怎么办
//!
//! 重新生成基准之后这条测试会红。**先确认新规模是有意的**(比如换了语料、
//! 换了抽样步长),再改下面的数,并在提交信息里说明为什么。
//! 绝不要"跑红了就把数改成实际值" —— 那正是这条测试要拦的事。

use std::path::PathBuf;

/// 每个入库基准该有多少行。**改这里之前先读模块文档。**
///
/// 这里只列 **omgkit-conf 的五条判官真的会读**的那四份。其余 13 份实测过:
/// `smoke.l1` / `smoke.l2*` / `smoke.matches.tsv` / `smarts.jsonl` 截短之后
/// `cargo test --release` 都会红(它们的读取方逐行比对,少一行就少一批断言,
/// 而且有特征覆盖闸),所以不用在这里重复钉。
///
/// `matches.tsv` 更进一步:它首行写着自己覆盖了多少个分子,`matches_large`
/// 断言那个数等于语料条数 —— 截短与"重导时加了 `--limit-mols`"都当场红。
///
/// **唯一的例外是 `smoke.l3.jsonl`**:149 行,截成 1 行全套测试照样绿 ——
/// 全仓库只有 `harness/README.md` 的生成命令提到它,**没有任何读取方**。
/// 那是一份死基准,不是截短风险;要么给它补个读取方,要么删掉。见任务清单。
const EXPECTED: &[(&str, usize)] = &[
    // omgkit-conf 的五条判官(CI 第 6–8、11–14 步)
    ("smoke.bounds.jsonl", 27),
    ("smoke.gram_eigs.jsonl", 27),
    ("smoke.chirality.jsonl", 150),
    ("smoke.lonepair.jsonl", 15),
];

fn baseline(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../harness/baseline")
        .join(name)
}

#[test]
fn 入库基准的行数是契约() {
    let mut bad = Vec::new();
    for &(name, want) in EXPECTED {
        let path = baseline(name);
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("读不了入库基准 {path:?}:{e}"));
        let got = text.lines().filter(|l| !l.trim().is_empty()).count();
        if got != want {
            bad.push(format!("  {name}:{got} 行,契约写着 {want} 行"));
        }
    }
    assert!(
        bad.is_empty(),
        "入库基准的规模变了:\n{}\n\n\
         判官读多少行就比多少个分子,它无从知道自己本该收到多少 —— 基准一变小,\n\
         每一条判据都会在更少的分子上算出更好看的数,而且全部退 0。\n\
         **先确认新规模是有意的**,再改 `EXPECTED`,并在提交信息里说明。",
        bad.join("\n")
    );
}

/// 每一行都必须是合法 JSON。
///
/// 行数对得上、内容坏掉,是这套闸门里曾经唯一还完全敞着的缝:非法行既不进
/// 判官的分母(先前 `n_lines += 1` 排在解析之后),也不改变这里的行数
/// (数的是非空行)。**两层闸同时失明。** 判官那边已经把行计数挪到解析之前,
/// 这里再钉一次 —— 那边只在被喂到这份基准时才看得见,这条测试每次都跑。
#[test]
fn 入库基准里没有坏行() {
    let mut bad = Vec::new();
    for &(name, _) in EXPECTED {
        let path = baseline(name);
        let text = std::fs::read_to_string(&path).expect("读得到");
        for (i, line) in text.lines().enumerate() {
            if line.trim().is_empty() {
                continue;
            }
            if serde_json::from_str::<serde_json::Value>(line).is_err() {
                bad.push(format!("  {name} 第 {} 行不是合法 JSON", i + 1));
            }
        }
    }
    assert!(bad.is_empty(), "入库基准里有坏行:\n{}", bad.join("\n"));
}

/// 两个特征分解基准必须是**同一批分子**导出来的。
///
/// `eigen_oracle` 用 `zip` 配对两个文件,而 `zip` 取短的那个 —— 行数不等时它
/// 静默截断,判据照样"跑得通",只是在比更少的分子。判官里已经有一条同样的
/// 前置检查,这里再钉一次:那边挡的是运行时,这里挡的是**入库的文件本身**。
#[test]
fn 界矩阵与特征分解两份基准同批() {
    let a = std::fs::read_to_string(baseline("smoke.bounds.jsonl")).expect("读得到");
    let b = std::fs::read_to_string(baseline("smoke.gram_eigs.jsonl")).expect("读得到");
    let (na, nb) = (a.lines().count(), b.lines().count());
    assert_eq!(
        na, nb,
        "smoke.bounds.jsonl 有 {na} 行而 smoke.gram_eigs.jsonl 有 {nb} 行 —— \
         两者必须是同一批分子导出来的"
    );
}
