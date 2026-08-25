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
/// 这里列的是 **omgkit-conf 的五条判官真的会读**的那四份,外加 `smoke.l3.jsonl`
/// (它的读取方在 `omgkit-io/tests/differential_l3.rs`,不在本 crate ——
/// 数只该有一份,所以仍旧钉在这里)。其余 12 份实测过:
/// `smoke.l1` / `smoke.l2*` / `smoke.matches.tsv` / `smarts.jsonl` 截短之后
/// `cargo test --release` 都会红(它们的读取方逐行比对,少一行就少一批断言,
/// 而且有特征覆盖闸),所以不用在这里重复钉。
///
/// `matches.tsv` 更进一步:它首行写着自己覆盖了多少个分子,`matches_large`
/// 断言那个数等于语料条数 —— 截短与"重导时加了 `--limit-mols`"都当场红。
///
/// `smoke.l3.jsonl` 一度是**唯一的例外**:149 行,截成 1 行全套测试照样绿 ——
/// 全仓库只有 `harness/README.md` 的生成命令提到它,没有任何读取方。那不是截短
/// 风险,是一份死基准。现在它有读取方了(`omgkit-io/tests/differential_l3.rs`,
/// 拿 RDKit 的规范串当"同一个分子的另一种写法",接上的第一次运行就抓出 5 条
/// 我方缺陷),所以行数也进了下面这张表。
///
/// 那条判据里钉死的 11 条例外顺带也是一道截短闸:截掉任何一条,"钉住的例外
/// 不见了"当场红。但那只覆盖那 11 行,行数仍旧要在这里钉。
const EXPECTED: &[(&str, usize)] = &[
    // omgkit-conf 的五条判官(CI 第 6–8、11–14 步)
    ("smoke.bounds.jsonl", 27),
    ("smoke.gram_eigs.jsonl", 27),
    ("smoke.chirality.jsonl", 150),
    ("smoke.lonepair.jsonl", 15),
    // 读取方在 omgkit-io:tests/differential_l3.rs
    ("smoke.l3.jsonl", 149),
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

/// 行数是契约,**行里装了多少东西也是**。
///
/// 上面那条数的是行。手性基准每一行装着若干个中心,而中心数可以在**行数不变**
/// 的情况下变少 —— 判官照样读满 150 行,只是每行比得少了,而它无从知道。
///
/// 这不是假想。提交 `61b8d58` 教会了 `dump_chirality.py` 收**三配位立体中心**
/// (亚砜的 S、膦的 P),却没有重导 `smoke.chirality.jsonl`:入库的那份里
/// 247 个中心**全是四配位**,而当时的脚本导出来是 248 个、其中 8 个三配位。
/// 行数一样(150),判官全绿,那个提交声称落地的那一档在主手性判官眼里
/// 根本不存在 —— 四个月没人看得见。
///
/// `harness/check_baseline_schema.py` 挡的是"脚本长了新字段、基准没重导"
/// (它比结构);这一条挡的是另一半:**结构没变而内容变少**。
///
/// 顺反那一列同理:`dump_chirality.py` 补上它之前,这 150 个分子里 23 个带
/// E/Z 的分子在两条手性判官眼里是**没有顺反**的分子(界矩阵少了解 1-4 顺反
/// 析取的依据)。数出来钉死,删掉那一列当场红。
#[test]
fn 手性基准装了多少中心与顺反也是契约() {
    // (基准, 中心总数, 其中三配位, 带顺反的双键)
    const WANT: &[(&str, usize, usize, usize)] = &[
        ("smoke.chirality.jsonl", 248, 8, 28),
        ("smoke.lonepair.jsonl", 21, 17, 0),
    ];
    let mut bad = Vec::new();
    for &(name, n_c, n_three, n_stereo) in WANT {
        let text = std::fs::read_to_string(baseline(name)).expect("读得到");
        let (mut c, mut three, mut st) = (0usize, 0usize, 0usize);
        for line in text.lines().filter(|l| !l.trim().is_empty()) {
            let v: serde_json::Value = serde_json::from_str(line).expect("合法 JSON");
            for x in v["centers"].as_array().into_iter().flatten() {
                c += 1;
                // 按**配体个数**数,不按 `three_coordinate` 那个标注 ——
                // 判官消费的是 `nbrs`(三配位那一档就是只有三个配体),
                // 标注只是导出时顺手写的。拿标注当契约会绕回自证。
                let three_c = x["nbrs"].as_array().map_or(0, Vec::len) == 3;
                assert_eq!(
                    three_c,
                    x["three_coordinate"].as_bool() == Some(true),
                    "{name} 里有个中心的 `three_coordinate` 标注与 `nbrs` 的长度不符"
                );
                if three_c {
                    three += 1;
                }
            }
            for b in v["bonds"].as_array().into_iter().flatten() {
                // 第 4 列是 RDKit 的 `BondStereo` 号(2 Z、3 E、4 cis、5 trans),
                // 第 5–6 列是两个参照原子。三样齐全才算这根键真的带顺反。
                let has = matches!(b.get(3).and_then(serde_json::Value::as_i64), Some(2..=5))
                    && b.get(4)
                        .and_then(serde_json::Value::as_i64)
                        .is_some_and(|x| x >= 0)
                    && b.get(5)
                        .and_then(serde_json::Value::as_i64)
                        .is_some_and(|x| x >= 0);
                if has {
                    st += 1;
                }
            }
        }
        if (c, three, st) != (n_c, n_three, n_stereo) {
            bad.push(format!(
                "  {name}:中心 {c}(契约 {n_c})、其中三配位 {three}(契约 {n_three})、\
                 带顺反的双键 {st}(契约 {n_stereo})"
            ));
        }
    }
    assert!(
        bad.is_empty(),
        "手性基准装的内容变了:\n{}\n\n\
         行数不变而中心/顺反变少,判官全都看不见 —— 它读满了行,只是每行比得少了。\n\
         **先确认新的数是有意的**,再改这里的契约,并在提交信息里说明。",
        bad.join("\n")
    );
}
