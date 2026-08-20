//! 写出去的串,别人读回来**再净化**,分子式不许变。
//!
//! # 为什么这条测试必须在 omgkit-chem 里
//!
//! `omgkit-io` 自己的往返测试(`roundtrip_smiles`)天生看不见这一类缺陷:
//! **io 的解析器不推隐式氢**(`num_implicit_hs` 恒为 0),氢要到净化才出现。
//! 于是"解析 → 写出 → 再解析"在氢数这件事上是**共谋**的:写出侧漏掉的氢,
//! 解析侧也不会补,两边一致,往返恒等照样成立。
//!
//! 氢是在 `omgkit-chem` 的净化里补的,所以只有跨到这一层,才问得出真正要问的
//! 那个问题:**别人照着我们写的串读回来,分子还是不是原来那个。**
//!
//! # 抓到过什么
//!
//! `smiles::write` 的 `hs_survive_without_brackets`(判断能不能省掉方括号)
//! 先前只看元素的**首位**默认价,而补氢取的是"第一个不小于已用价的允许价"。
//! 多价元素上两者不是一回事:
//!
//! | | 键级和 | 价表 | 读者补的氢 | 写出侧当时的判断 |
//! |---|---|---|---|---|
//! | `Cl[I]Cl` 的 I | 2 | `[1, 3, 5]` | 1(补到 3 价) | `2 >= 1` → 去掉方括号 |
//! | `NC[S](=O)=O` 的 S | 5 | `[2, 4, 6]` | 1(补到 6 价) | `5 >= 2` → 去掉方括号 |
//!
//! 于是写出 `ClICl` / `NCS(=O)=O`,读回来再净化就多一个氢 —— **是另一个分子**。
//! 全语料 9 条如此,而当时全部单元测试与往返测试都是绿的;发现它的是外部判据
//! `harness/check_write.py --canonical`,而那条判据当时不在 CI 里。
//!
//! # 两条写出路都要走
//!
//! 缺陷只在**不净化**那条路上显形,而那正是外部判据走的路
//! (`examples/write_smiles.rs` 只解析、不净化)。净化之后再写出是另一条路:
//! 净化会把隐式氢挪进 `num_explicit_hs`,判断因此走另一个分支。
//!
//! 两条都测,而且**比对一律在净化之后做** —— 只有净化才会补氢,不净化的话
//! 两边的 `num_implicit_hs` 都是 0,分子式当然相同,判据形同虚设。
//! (这一点是实测的:不净化比对时,`Cl[I]Cl` 与 `ClICl` 的分子式逐字相同。)
//!
//! 变异实测,三种变异各自被哪条抓住:
//!
//! | 变异 | 不净化_全量 | 多价元素的方括号不许省 | 净化后_全量 | 不净化_冒烟 |
//! |---|---|---|---|---|
//! | **A** 恒返回 `false`(一律留框) | 绿 | **红** `[SH2]` 该省没省 | 绿 | 绿 |
//! | **B** 删掉芳香分支 | 绿 | **红** `c1cc[sH]c1`→`c1cccs1` | 绿 | 绿 |
//! | **C** 退回只看首位默认价 | **红** 8 条分子式变了 | **红** `Cl[I]Cl`→`ClICl` | 绿 | 绿 |
//!
//! 几件要如实记下来、别当成"都验过了"的事:
//!
//! - **A 和 B 只有逐条那个测试抓得住。** A 改了全语料 **1044 行**输出、B 改了 0 行,
//!   而两条外部判据(`check_write.py` 两个方向)对两者**全绿** —— 多写方括号
//!   语义不变,外部判官看不见。
//! - **净化后那两条在三种变异下全绿。** 净化会把隐式氢挪进 `num_explicit_hs`,
//!   判断走的是另一个分支,那个分支眼下恰好是对的。它们守的是那条路将来
//!   单独回归的可能,不是这次的缺陷。
//! - 全量比冒烟多抓 8 条,是因为 `Cl[I]Cl` 那一档分子**只在大语料里有**。
//!
//! (A 这条变异是独立审核做出来的。先前 `bare` 那一组用的是 `CS(=O)(=O)C` 之类
//!  **氢数为 0** 的原子,`hs_survive_without_brackets` 根本不会被调用 —— 五条测试
//!  在 A 下全绿。现在换成带氢的 `[SH2]`/`[CH4]`/`[NH3]`/`[OH2]`,断言也换成全等。)
//!
//! # 判据为什么比分子式
//!
//! 规范写出会**重排原子**,所以逐原子比对要先建映射,而建映射本身就要用到
//! 规范秩 —— 那又是被测代码。分子式(逐元素计数,氢单独算)与顺序无关,
//! 多一个氢当场就看得见,不需要任何来自被测代码的辅助。
//!
//! 它当然不是全能的:重排原子而不改分子式的缺陷它看不见。
//!
//! **别指望 `roundtrip_smiles` 补这一段** —— 它只调 `smiles::write` 与
//! `write_with_priority`,两者都是 `WriteStyle::Faithful`,而 Faithful 的
//! `needs_brackets` 压根不会走到 `hs_survive_without_brackets`。它**一次都没跑过
//! 规范写出这条路**。(先前这里写着"那一类由 roundtrip_smiles 守着",
//! 是独立审核查出来的假话。)
//!
//! 规范路上"重排 / 改键级而分子式不变"这一类,眼下唯一的守卫是外部判据
//! `harness/check_write.py --canonical`(两边各自规范化再比串),
//! 而它现在进了 CI。

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use omgkit_core::MolBuilder;
use omgkit_io::{canon, smiles};

fn corpus(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../harness/corpus")
        .join(name)
}

/// 逐元素原子数;氢单独一档(显式氢原子与隐式氢数合并计)。
fn formula(mol: &MolBuilder) -> BTreeMap<u8, u32> {
    let mut f = BTreeMap::new();
    for a in mol.atoms() {
        *f.entry(a.atomic_num).or_insert(0) += 1;
        let hs = u32::from(a.num_explicit_hs) + u32::from(a.num_implicit_hs);
        if hs != 0 {
            *f.entry(1).or_insert(0) += hs;
        }
    }
    f
}

fn render(f: &BTreeMap<u8, u32>) -> String {
    f.iter()
        .map(|(z, n)| {
            let sym = omgkit_core::element::by_atomic_num(*z).map_or("?", |e| e.symbol);
            format!("{sym}{n}")
        })
        .collect::<Vec<_>>()
        .join("")
}

/// 解析 + 净化,拿到"下游真正会看到的那个分子"。净化不过就返回 `None`。
fn sanitized(smi: &str) -> Option<MolBuilder> {
    let mut m = smiles::parse(smi).ok()?;
    omgkit_chem::pipeline::sanitize(&mut m).ok()?;
    Some(m)
}

struct Failure {
    smi: String,
    written: String,
    before: String,
    after: String,
}

/// 跑一份语料。`sanitize_first` 决定写出的是净化前还是净化后的分子 ——
/// 见模块文档"两条写出路都要走"。
///
/// 返回 (进入比对的分子数, 失败列表)。
fn check(path: &Path, sanitize_first: bool) -> (usize, Vec<Failure>) {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("读不了 {path:?}:{e}"));
    let mut n = 0;
    let mut bad = Vec::new();
    for line in text.lines() {
        let smi = line.split_whitespace().next().unwrap_or("");
        if smi.is_empty() || smi.starts_with('#') {
            continue;
        }
        // 参照永远是"输入净化之后"的分子:下游看到的就是它。
        // 净化不过的分子谈不上"分子式变没变",跳过。
        let Some(want) = sanitized(smi) else {
            continue;
        };
        let source = if sanitize_first {
            want.clone()
        } else {
            let Ok(m) = smiles::parse(smi) else { continue };
            m
        };
        n += 1;
        let written = canon::canonical_smiles(&source).smiles;
        let f0 = formula(&want);
        let Some(back) = sanitized(&written) else {
            bad.push(Failure {
                smi: smi.to_string(),
                written,
                before: render(&f0),
                after: "写出来的串读不回/净化不过".into(),
            });
            continue;
        };
        let f1 = formula(&back);
        if f0 != f1 {
            bad.push(Failure {
                smi: smi.to_string(),
                written,
                before: render(&f0),
                after: render(&f1),
            });
        }
    }
    (n, bad)
}

fn report(bad: &[Failure], limit: usize) -> String {
    let mut s = String::new();
    for f in bad.iter().take(limit) {
        s.push_str(&format!(
            "\n  {}\n    写出 {}\n    分子式 {} -> {}",
            f.smi, f.written, f.before, f.after
        ));
    }
    if bad.len() > limit {
        s.push_str(&format!("\n  …… 另有 {} 条", bad.len() - limit));
    }
    s
}

/// **分母一起断言。** 语料读不到 / 全军覆没时,"零失败"是假的。
fn assert_clean(n: usize, floor: usize, bad: &[Failure], what: &str) {
    assert!(
        n >= floor,
        "{what}:只有 {n} 个分子进了比对,分母不对(至少该有 {floor})"
    );
    assert!(
        bad.is_empty(),
        "{what}:{} 条分子式变了:{}",
        bad.len(),
        report(bad, 10)
    );
}

#[test]
fn 不净化就写出_读回来分子式不变_冒烟() {
    let (n, bad) = check(&corpus("smoke.smi"), false);
    assert_clean(n, 135, &bad, "冒烟 / 不净化写出");
}

/// **`Cl[I]Cl` 那一类只在大语料里出现**,冒烟档一条都没有 —— 所以这条不标
/// `#[ignore]`:它是唯一见得到那一档的测试,标了等于没写。
/// 五条测试合计 release 下 0.33 秒、debug 下 0.45 秒。
#[test]
fn 不净化就写出_读回来分子式不变_全量() {
    let (n, bad) = check(&corpus("large.smi"), false);
    assert_clean(n, 8800, &bad, "全量 / 不净化写出");
}

#[test]
fn 净化后写出_读回来分子式不变_冒烟() {
    let (n, bad) = check(&corpus("smoke.smi"), true);
    assert_clean(n, 135, &bad, "冒烟 / 净化后写出");
}

#[test]
fn 净化后写出_读回来分子式不变_全量() {
    let (n, bad) = check(&corpus("large.smi"), true);
    assert_clean(n, 8800, &bad, "全量 / 净化后写出");
}

/// 逐条钉死几个**多价 / 超价**的边界。上面那几条全量测试真正在守的就是这一档,
/// 单独写出来是为了失败时一眼看得出是哪一类,而不是从八千条里翻。
///
/// 每条都跑**两条写出路**(净化不过的分子只跑不净化那条)。
#[test]
fn 多价元素的方括号不许省() {
    /// 两条路各写一遍;净化不过的只有不净化那一条。
    fn forms(smi: &str) -> Vec<String> {
        let raw = smiles::parse(smi).expect("测试用的 SMILES 该能解析");
        let mut v = vec![canon::canonical_smiles(&raw).smiles];
        if let Some(m) = sanitized(smi) {
            v.push(canon::canonical_smiles(&m).smiles);
        }
        v
    }

    // (输入, 写出里必须出现的片段)
    let keep = [
        // I 价表 [1,3,5]:两根键会被补到 3 价 = 多一个氢
        ("Cl[I]Cl", "[I]"),
        // S 价表 [2,4,6]:五价会被补到 6 价 = 多一个氢
        ("NC[S](=O)=O", "[S]"),
        // P 价表 [3,5]:六价**超过全部允许价**,读者补几个氢取决于它那份表 → 留框。
        // (这个分子过不了本仓的严格净化,所以它只跑不净化那条路。)
        ("F[P](F)(F)(F)(F)F", "[P]"),
        // 中性碳零键:去框写成 `C` 读回来补四个氢
        ("[C]", "[C]"),
        // **芳香那一档。** 噻吩的 S:两根芳香键 = 键级和 3,已经超过首位默认价 2。
        // 非芳香规则会算"补到 4 价 = 1 个氢,与实际相符 → 去框",写出 `c1cccs1`,
        // 而读者按芳香分支读回来是 0 个氢 —— 氢丢了。
        // 变异实测:删掉写出侧的芳香分支,全语料写出**一行都不变**、五条测试与
        // 两条外部判据全绿,只有这一条会红。它是整块改动里唯一靠它守住的逻辑。
        ("c1cc[sH]c1", "[sH]"),
    ];
    for (smi, want) in keep {
        for w in forms(smi) {
            assert!(w.contains(want), "{smi} 写成了 {w},里面该有 {want}");
        }
    }

    // 反过来:**该省的还得省**。
    //
    // 这一组先前用的是 `CS(=O)(=O)C` / `CSC` / `CI` / `COP(=O)(OC)OC`,而那四个
    // 原子**氢数都是 0** —— `needs_brackets` 里 `author_fixed_hs` 为假,
    // `hs_survive_without_brackets` **根本不会被调用**。于是那一组是空断言:
    // 独立审核实测,把该函数改成恒返回 `false`(= 一律留框),
    // 全语料写出改了 **1044 行**,而五条测试与两条外部判据**全部照旧全绿**。
    //
    // 换成**带氢**的原子才走得到那条分支:它们的键级和 + 氢数恰好等于首位默认价,
    // 裸写读回来补出的正是这些氢。断言也换成**全等**,不再用"不含某片段" ——
    // 那种否定断言遇到 `[SH0]`、`[S+0]` 这类写法照样能空过。
    let bare = [
        ("[SH2]", "S"), // 硫化氢:一价都没有 + 2 个氢 = 首位默认价 2
        ("[CH4]", "C"), // 甲烷
        ("[NH3]", "N"), // 氨
        ("[OH2]", "O"), // 水
    ];
    for (smi, want) in bare {
        for w in forms(smi) {
            assert_eq!(w, want, "{smi} 写成了 {w},该省的方括号没省");
        }
    }
}
