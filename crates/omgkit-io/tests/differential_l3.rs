//! l3 差分:**规范 SMILES 的收敛判据**,参照是 RDKit 写出的“同一个分子的另一种写法”。
//!
//! 基准由 `harness/oracle_pipeline.py --stage l3 --remove-hs` 生成
//! (RDKit 2025.09.2,与 `harness/requirements.lock` 同版本),每行带
//! `can` 一列 —— RDKit 自己的规范串。判据是:
//!
//! ```text
//! canonical(解析(语料里的写法))  ==  canonical(解析(RDKit 的规范串))
//! ```
//!
//! 外加三列直接比:能不能读(解析+净化成不成功)、去氢后的原子数、键数。
//!
//! # 与已有的三条规范化判据都不重叠
//!
//! | 判据 | 换掉的是什么 | 第二种写法从哪来 |
//! |---|---|---|
//! | `canonical_smiles_is_invariant_under_renumbering_*` | 原子编号与键序 | 自己重排 |
//! | `canonical_smiles_is_a_fixed_point`(经 wheel) | 写出再读回 | **自己的写出器** |
//! | `different_spellings_converge_to_one_canonical_form` | 起点/分支/环闭合 | **手挑的 5 组** |
//! | 本文件 | 以上全部,外加芳香/凯库勒、氢的记法、方向符号的落点 | **RDKit**,149 + 8839 条真语料 |
//!
//! 前两条的第二种写法都出自我方:自家写出器漏掉的东西,读回来一样漏,两边共谋着
//! “收敛”。第三条的写法是手挑的 —— 只走得到写的人想得到的形态。这一条的另一种
//! 写法由外部实现给出,没人商量过。
//!
//! # 这份基准死了很久
//!
//! `smoke.l3.jsonl` 入库以来**没有任何读取方**,截成一行全套测试照样绿
//! (`omgkit-conf/tests/baseline_sizes.rs` 的模块文档里记着这件事)。接上读取方的
//! 第一次运行就抓出 8 条不收敛,其中 5 条是我方的缺陷 —— 配位键旁边的原子该不该
//! 留方括号,取决于氢是写在括号里还是推断出来的,于是同一个分子有两串
//! (见 `smiles/write.rs::hs_survive_without_brackets`)。另外 3 条在参照那一侧,
//! 逐条钉在 [`NOT_CONVERGENT`]。
//!
//! 顺带发现基准与生成它的脚本也脱了钩:`--remove-hs` 当时在**解析**那一步做,
//! 而生成器后来改成了“不净化地解析 + 显式跑净化”。两者一撞,RDKit 会把方括号里的
//! 氢数、`noImplicit` 标志和手性标记一并抹掉,三条手性用例被换成了另一个分子。
//! 已修(去氢挪到净化之后),基准已重导。

use std::path::{Path, PathBuf};

use omgkit_core::MolBuilder;
use omgkit_io::{canon, smiles};

/// 不收敛、且**已查明原因在参照那一侧**的条目(冒烟语料)。
///
/// 这张表是双向的:多一条新的不收敛红,少一条也红 —— 后者逼着改的人回来确认
/// “消失”是对的,而不是把判据悄悄放松了一格。
const NOT_CONVERGENT: &[(&str, &str)] = &[
    (
        "[C@@H]1CCCCC1O",
        "RDKit 的写出器不写这个原子的手性(2 根重原子键 + 方括号里 1 个氢 + 1 个\
         自由基电子),分子对象里标记还在(CHI_TETRAHEDRAL_CCW),写出来就没了。\
         于是 chiral-ring-open-cw 与 -ccw 两条语料塌成同一串 OC1[CH]CCCC1。\
         我方两串不同 —— 本仓有一条判据专门守“不同的分子不能塌到同一个规范串”。",
    ),
    ("[C@H]1CCCCC1O", "同上,是上一条的对映体。"),
    (
        "OI(=O)O",
        "**RDKit 自己的规范串不是自己的不动点。** 它把 OI(=O)O 写成 O=[IH](O)O,\
         而 O=[IH](O)O 再读回去,净化第 1 步的卤素修正就会触发(显式写出的那个氢\
         把显式价从 4 顶到 5),得到 [O-][IH+](O)O。我方逐字节复现了这个行为 ——\
         两次都与 RDKit 2025.09.2 一致,分歧只在参照的两端之间。",
    ),
];

/// 规范串里**有没有立体标记**,逐条与参照比;分歧全部钉死,写明是哪一侧写不出。
///
/// 收敛判据有一个天然的盲区:**两边一起把立体标记丢干净,照样收敛**。
/// 语料里那几条非四面体的用例正是这样 —— 我方两次都写不出,两串相同,
/// 收敛判据全绿,而分子的立体信息一次都没落到纸上。所以这一列单独比、单独钉。
///
/// 平面四方那三条**已经修好了**(`@SP` 写得出来了),表里只剩三角双锥与八面体。
/// 写出器补上 `@TB`/`@OH` 之后这张表还会红,那时再删掉对应行
/// (`harness/check_write.py` 的 `NON_TETRAHEDRAL_GAP` 是同一件事的另一半)。
///
/// 修 `@SP` 的时候顺带纠正了这里先前写着的一句错话:`[Pt@SP1](Cl)(Cl)(N)N` 与
/// `[Pt@SP3](Cl)(Cl)(N)N` **不是**顺铂与反铂,而是**同一个异构体的两种写法** ——
/// 序号说的是“按列出顺序哪两对互为反位”,而两个氯彼此可换,于是这两个序号在
/// 这个分子上表达同一件事。判据是 RDKit 自己给的:把它重排原子再规范化,RDKit
/// 会在 `@SP1` 与 `@SP3` 两串之间跳(`@SP2` 只有一串)。我方两者收敛到同一串。
const STEREO_MISMATCH: &[(&str, &str)] = &[
    (
        "[C@@H]1CCCCC1O",
        "参照写不出(理由见 NOT_CONVERGENT 的同名条目)",
    ),
    ("[C@H]1CCCCC1O", "同上"),
    ("F[P@TB15](Cl)(Br)(I)S", "我方写不出:三角双锥 `@TB`"),
    ("C[Co@OH25](N)(O)(S)(P)Cl", "我方写不出:八面体 `@OH`"),
    (
        "[Co@OH5]1(N)(O)(S)(P)CCC1",
        "我方写不出:八面体 `@OH`(带环的)",
    ),
];

fn baseline(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../harness/baseline")
        .join(name)
}

/// **产品那条路** —— 与 `harness/oracle_pipeline.py` 的 l3 阶段逐步对应:
/// 解析 → 净化 → 去氢 → 感知顺反 → 规范化。
///
/// 去氢排在净化之后,不是解析里顺手做掉:生成器那侧同样如此,而且是被迫的 ——
/// 放在解析那一步 RDKit 会把方括号里的氢连同手性标记一起抹掉(见模块文档)。
fn product_path(smi: &str) -> Result<MolBuilder, String> {
    let mut m = smiles::parse(smi).map_err(|e| format!("解析:{}", e.render()))?;
    omgkit_chem::sanitize(&mut m).map_err(|e| format!("净化:{e}"))?;
    omgkit_chem::remove_hs(&mut m);
    omgkit_io::stereo::perceive_bond_stereo(&mut m);
    Ok(m)
}

/// 串里有几个带立体标记的原子。`@@` 算一个,`@SP1`/`@TB15`/`@OH2` 也各算一个。
///
/// 不能直接数 `@` 字符:`@@` 是两个字符一个中心,而两侧的规范串原子次序本就不同,
/// 同一个中心在一边写 `@`、另一边写 `@@` 是正常的(参照系不同),数字符会假红。
fn n_stereo(s: &str) -> usize {
    let b = s.as_bytes();
    let mut i = 0;
    let mut n = 0;
    while i < b.len() {
        if b[i] == b'@' {
            n += 1;
            i += if i + 1 < b.len() && b[i + 1] == b'@' {
                2
            } else {
                1
            };
        } else {
            i += 1;
        }
    }
    n
}

struct Bad {
    smi: String,
    what: String,
    reference: String,
    ours: String,
}

/// 一趟比对的产物。两张“分歧清单”是给钉死的例外表对账用的。
#[derive(Default)]
struct Stats {
    /// 基准行数(**在解析之前数**,坏行不能悄悄从分母里消失)
    lines: usize,
    /// 两边都读得进来的分子数
    ok: usize,
    /// 参照的规范串里带立体标记的分子数
    stereo_ref: usize,
    /// 我方的规范串里带立体标记的分子数
    stereo_ours: usize,
    ez: usize,
    multi: usize,
    dative: usize,
    /// `can_renumbered` 与 `can` 真的不同的分子数 —— 那一列还带不带信息
    renumbered_distinct: usize,
    /// 不收敛的条目(不论有没有被钉住)
    diverged: Vec<String>,
    /// 立体标记有无与参照不同的条目(同上)
    stereo_mismatch: Vec<String>,
}

fn diff_against(
    path: &Path,
    pinned: &[(&str, &str)],
    pinned_stereo: &[(&str, &str)],
) -> (Stats, Vec<Bad>) {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| {
        panic!(
            "读不到 l3 基准 {}: {e}\n生成方式见 harness/README.md",
            path.display()
        )
    });

    let mut st = Stats::default();
    let mut bad: Vec<Bad> = Vec::new();

    for line in text.lines() {
        if line.trim().is_empty() {
            continue;
        }
        st.lines += 1;
        let v: serde_json::Value = serde_json::from_str(line)
            .unwrap_or_else(|e| panic!("第 {} 行不是合法 JSON:{e}", st.lines));
        let smi = v["smi"].as_str().expect("每行都有 smi").to_string();
        let rd_ok = v["ok"].as_bool().expect("每行都有 ok");
        let ours = product_path(&smi);

        // 一、能不能读:两边必须同时成功或同时失败
        if ours.is_ok() != rd_ok {
            bad.push(Bad {
                smi,
                what: "能否解析+净化".into(),
                reference: if rd_ok { "成功" } else { "失败" }.into(),
                ours: match &ours {
                    Ok(_) => "成功".into(),
                    Err(e) => format!("失败({e})"),
                },
            });
            continue;
        }
        let Ok(m) = ours else { continue };
        st.ok += 1;

        // 二、去氢之后的原子数与键数
        for (field, got, want) in [
            (
                "原子数",
                m.atoms().len() as u64,
                v["na"].as_u64().expect("有 na"),
            ),
            (
                "键数",
                m.bonds().len() as u64,
                v["nb"].as_u64().expect("有 nb"),
            ),
        ] {
            if got != want {
                bad.push(Bad {
                    smi: smi.clone(),
                    what: field.into(),
                    reference: want.to_string(),
                    ours: got.to_string(),
                });
            }
        }

        let rd_can = v["can"].as_str().expect("l3 每条成功记录都有 can");
        let ours_from_smi = canon::canonical_smiles(&m).smiles;

        // 覆盖计数:全绿必须是“实现对”,不能是“语料里压根没有那个形态”
        if rd_can.contains('@') {
            st.stereo_ref += 1;
        }
        if ours_from_smi.contains('@') {
            st.stereo_ours += 1;
        }
        if rd_can.contains('/') || rd_can.contains('\\') {
            st.ez += 1;
        }
        if rd_can.contains('.') {
            st.multi += 1;
        }
        if rd_can.contains("->") || rd_can.contains("<-") {
            st.dative += 1;
        }

        // 三、带立体标记的原子数 —— 收敛判据对“两边一起丢”是瞎的。
        // 比的是**个数**不是有无:三个中心写出两个,只比有无是看不见的。
        if n_stereo(rd_can) != n_stereo(&ours_from_smi) {
            st.stereo_mismatch.push(smi.clone());
            if !pinned_stereo.iter().any(|&(s, _)| s == smi) {
                bad.push(Bad {
                    smi: smi.clone(),
                    what: format!(
                        "带立体标记的原子数 {} vs {}",
                        n_stereo(rd_can),
                        n_stereo(&ours_from_smi)
                    ),
                    reference: rd_can.into(),
                    ours: ours_from_smi.clone(),
                });
            }
        }

        // 四、收敛:RDKit 的规范串读回来,规范化之后必须是同一串。
        //
        // `can_renumbered` 是 RDKit 把原子顺序反转后再写一遍的结果,与 `can`
        // 不同时就是**第三种写法**,一并喂进来。语料里恰好有一条:
        // `[Co@OH5]1(N)(O)(S)(P)CCC1` —— RDKit 自己的规范化在八面体立体上不是
        // 重排不变的(`[Co@OH2]` vs `[Co@OH30]`)。这一列先前没有任何读取方。
        let rd_ren = v["can_renumbered"].as_str().unwrap_or(rd_can);
        if rd_ren != rd_can {
            st.renumbered_distinct += 1;
        }
        let mut references = vec![rd_can];
        if rd_ren != rd_can {
            references.push(rd_ren);
        }
        for reference in references {
            match product_path(reference) {
                Err(e) => bad.push(Bad {
                    smi: smi.clone(),
                    what: "读 RDKit 的规范串".into(),
                    reference: reference.into(),
                    ours: e,
                }),
                Ok(m2) => {
                    let ours_from_can = canon::canonical_smiles(&m2).smiles;
                    if ours_from_can != ours_from_smi {
                        if !st.diverged.contains(&smi) {
                            st.diverged.push(smi.clone());
                        }
                        if !pinned.iter().any(|&(s, _)| s == smi) {
                            bad.push(Bad {
                                smi: smi.clone(),
                                what: format!("不收敛(RDKit 的规范串是 {reference})"),
                                reference: ours_from_can,
                                ours: ours_from_smi.clone(),
                            });
                        }
                    }
                }
            }
        }
    }

    (st, bad)
}

/// 钉住的例外必须**条条还在**。某条不再分歧了是好消息,但表要跟着改 ——
/// 留着就等于把判据放松了一格,而且没人会注意到。
fn assert_still_pinned(pinned: &[(&str, &str)], seen: &[String], what: &str) {
    let gone: Vec<&str> = pinned
        .iter()
        .map(|&(s, _)| s)
        .filter(|s| !seen.iter().any(|d| d == s))
        .collect();
    assert!(
        gone.is_empty(),
        "{what} 里这几条已经不分歧了:{gone:?}\n\
         这是好消息,但表要跟着删 —— 留着就是把判据放松了一格。"
    );
}

fn report(st: &Stats, bad: &[Bad]) -> String {
    let mut out = format!(
        "\nl3 差分失败:{} 行({} 条读得进来),{} 处分歧\n\n",
        st.lines,
        st.ok,
        bad.len()
    );
    for b in bad.iter().take(20) {
        out.push_str(&format!(
            "  {}\n      {:<30} 参照={}\n      {:<30} 我方={}\n",
            b.smi, b.what, b.reference, "", b.ours
        ));
    }
    if bad.len() > 20 {
        out.push_str(&format!("  ...(另有 {} 处)\n", bad.len() - 20));
    }
    out.push_str(
        "\n判据是“同一个分子的两种写法必须收敛到同一个规范串”,第二种写法由 RDKit 写出。\n\
         查明原因在参照那一侧时,把条目加进 NOT_CONVERGENT 并写清为什么。\n",
    );
    out
}

#[test]
fn l3_冒烟语料的规范串与_rdkit_的写法收敛() {
    let (st, bad) = diff_against(&baseline("smoke.l3.jsonl"), NOT_CONVERGENT, STEREO_MISMATCH);
    assert!(st.lines > 0, "l3 基准是空的");
    assert_still_pinned(NOT_CONVERGENT, &st.diverged, "NOT_CONVERGENT");
    assert_still_pinned(STEREO_MISMATCH, &st.stereo_mismatch, "STEREO_MISMATCH");
    assert!(bad.is_empty(), "{}", report(&st, &bad));

    assert!(st.stereo_ref > 0, "语料里一条立体标记都没有,那一路是空过的");
    assert!(st.ez > 0, "语料里一条顺反都没有,那一路是空过的");
    assert!(st.multi > 0, "语料里没有多组分分子,那一路是空过的");
    assert!(st.dative > 0, "语料里没有配位键,那一路是空过的");
    assert!(
        st.renumbered_distinct > 0,
        "没有一条的 can_renumbered 与 can 不同 —— 那一列不再带信息了。\n\
         语料里本来有一条(RDKit 的八面体立体不是重排不变的)。\
         若是参照那侧修好了,把这条断言删掉即可。"
    );

    println!(
        "l3 冒烟收敛判据通过:{} 行,{} 条读得进来;带立体标记的分子 {}/{}(参照/我方,\
         分歧 {} 条已钉),顺反 {},多组分 {},配位键 {};不收敛的例外 {} 条",
        st.lines,
        st.ok,
        st.stereo_ref,
        st.stereo_ours,
        STEREO_MISMATCH.len(),
        st.ez,
        st.multi,
        st.dative,
        NOT_CONVERGENT.len(),
    );
}

/// 大语料(8839 条,取自公开的 NCI / ZINC 子集)。**一条例外都没有。**
///
/// 冒烟语料那 11 条钉死的例外全是刻意造的边角(自由基上的手性、非四面体立体、
/// 卤素修正的触发边界),真实语料里一条都没出现:8831 个读得进来的分子上,
/// 收敛、立体标记有无、原子数、键数、能否解析,五列全是零分歧。
///
/// 所以这里两张例外表都传空的 —— 大语料上出现**任何**一条,都必须先解释清楚。
///
/// 基准体积 1.5 MB,按 `.gitignore` 的规矩不入库,生成命令见 `harness/README.md`。
#[test]
#[ignore = "需要先生成大语料基准,见函数文档;用 cargo test -- --ignored 运行"]
fn l3_大语料的规范串与_rdkit_的写法收敛() {
    let (st, bad) = diff_against(&baseline("large.l3.jsonl"), &[], &[]);
    assert!(
        st.lines > 1000,
        "大语料基准看起来不完整:只有 {} 行",
        st.lines
    );
    assert!(bad.is_empty(), "{}", report(&st, &bad));
    assert!(
        st.stereo_ref > 0,
        "大语料里一条立体标记都没有,那一路是空过的"
    );
    assert!(st.ez > 0, "大语料里一条顺反都没有,那一路是空过的");
    println!(
        "l3 大语料收敛判据通过:{} 行,{} 条读得进来,零分歧;\
         带立体标记的分子 {}/{}(参照/我方),顺反 {},多组分 {}",
        st.lines, st.ok, st.stereo_ref, st.stereo_ours, st.ez, st.multi
    );
}
