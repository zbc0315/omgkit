//! 子结构匹配的基本正确性。
//!
//! 这一层测的是"匹配算得对不对",判据是手工可验的小分子。大规模的正确性
//! 由 `differential.rs` 对着外部实现比。

use omgkit_chem::sanitize;
use omgkit_core::MolBuilder;
use omgkit_io::{canon, smarts, smiles};
use omgkit_match::{substructure_matches, MatchOptions, MolProps};

/// 解析并跑完整的净化 —— 匹配要用到隐式氢、芳香标志、环信息。
fn sanitized(smi: &str) -> MolBuilder {
    let mut m = smiles::parse(smi).unwrap_or_else(|e| panic!("{smi}: {}", e.render()));
    sanitize(&mut m).unwrap_or_else(|e| panic!("{smi}: {e}"));
    m
}

fn count(pat: &str, smi: &str) -> usize {
    hits(pat, smi, MatchOptions::default()).len()
}

fn hits(pat: &str, smi: &str, opts: MatchOptions) -> Vec<Vec<u32>> {
    let q = smarts::parse(pat).unwrap_or_else(|e| panic!("{pat}: {}", e.render()));
    let m = sanitized(smi);
    let props = MolProps::compute(&m);
    substructure_matches(&q, &m, &props, opts)
}

#[test]
fn single_atom_patterns() {
    assert_eq!(count("[OH]", "CCO"), 1);
    assert_eq!(count("C", "CCO"), 2, "两个脂肪碳");
    assert_eq!(count("[#6]", "CCO"), 2);
    assert_eq!(count("c", "c1ccccc1"), 6);
    assert_eq!(count("C", "c1ccccc1"), 0, "芳香碳不是脂肪碳");
    assert_eq!(count("*", "CCO"), 3);
}

#[test]
fn two_atom_patterns() {
    assert_eq!(count("CO", "CCO"), 1);
    assert_eq!(count("CC", "CCO"), 1);
    assert_eq!(count("CC", "CCC"), 2, "丙烷有两条 C-C");
    assert_eq!(count("C=O", "CC=O"), 1);
    assert_eq!(count("C=O", "CCO"), 0, "单键不是双键");
}

/// 环上的匹配 —— 苯环里的 `cc` 有 6 条键。
#[test]
fn ring_patterns() {
    assert_eq!(count("cc", "c1ccccc1"), 6);
    assert_eq!(count("c1ccccc1", "c1ccccc1"), 1, "整个环,去重后只剩一种");
    assert_eq!(count("[R1]", "c1ccccc1"), 6);
    assert_eq!(count("[R2]", "c1ccc2ccccc2c1"), 2, "萘的两个稠合碳");
    assert_eq!(
        count("[r5]", "C1CC2CCCCC2C1"),
        5,
        "氢化茚:最小环为 5 的原子"
    );
    assert_eq!(count("[r6]", "C1CC2CCCCC2C1"), 4);
}

/// 去重开关的差别:同一组原子的多个排列。
#[test]
fn uniquify_changes_the_count() {
    let uniq = hits(
        "c1ccccc1",
        "c1ccccc1",
        MatchOptions {
            uniquify: true,
            ..Default::default()
        },
    );
    let all = hits(
        "c1ccccc1",
        "c1ccccc1",
        MatchOptions {
            uniquify: false,
            ..Default::default()
        },
    );
    assert_eq!(uniq.len(), 1, "按原子集合去重后只剩一种");
    assert_eq!(all.len(), 12, "苯的对称群阶为 12(6 个旋转 × 2 个方向)");
}

/// `max_matches` 早停。高度对称的分子上匹配数会爆炸,
/// 只问"有没有"时不该把它们全枚举出来。
#[test]
fn max_matches_stops_early() {
    let one = hits(
        "cc",
        "c1ccccc1",
        MatchOptions {
            max_matches: 1,
            uniquify: false,
            use_chirality: true,
        },
    );
    assert_eq!(one.len(), 1);
}

/// 逻辑运算与查询基元。
#[test]
fn query_expressions() {
    assert_eq!(count("[C,N]", "CCN"), 3);
    assert_eq!(count("[C;H3]", "CCO"), 1, "只有一个甲基");
    assert_eq!(count("[!C]", "CCO"), 1, "氧");
    assert_eq!(count("[CX4]", "CC=O"), 1, "只有甲基是四连接");
    assert_eq!(count("[$(CO)]", "CCO"), 1, "递归:连着氧的碳");
    assert_eq!(count("[$(CC)]", "CCO"), 2, "两个碳都连着碳");
}

/// 映射给出的原子下标必须真的对上。
#[test]
fn mapping_points_at_the_right_atoms() {
    // CCO:原子 0=C 1=C 2=O
    let h = hits("CO", "CCO", MatchOptions::default());
    assert_eq!(h.len(), 1);
    assert_eq!(h[0], vec![1, 2], "查询原子 0(C)→1,原子 1(O)→2");
}

/// 断开的模式(多片段)要能匹配到分子的不同部分。
#[test]
fn disconnected_patterns() {
    assert_eq!(count("C.O", "CCO"), 2, "两个碳各配一次氧");
    assert_eq!(count("[OH].[OH]", "OCCO"), 1, "两个羟基,去重后一种");
}

/// 模式比分子大时直接无解,不该崩。
#[test]
fn pattern_larger_than_molecule() {
    assert_eq!(count("CCCC", "CC"), 0);
}

/// 配位键的方向有语义:`->` 与 `<-` 匹配的是相反的朝向。
#[test]
fn dative_direction_matters() {
    assert_eq!(count("N->[Cu]", "N->[Cu]"), 1);
    assert_eq!(count("[Cu]<-N", "N->[Cu]"), 1, "同一条键,反着写也对");
    assert_eq!(count("[Cu]->N", "N->[Cu]"), 0, "方向反了就不匹配");
}

/// 手性匹配要换参照系,不能比原始标记。
///
/// 标记相对**各自分子的邻居存储顺序**,查询与底物的顺序不同。直接比原始标记
/// 是拿两个参照系里的值去比,得到的构型可以正好相反 —— 而这类错误只在写全
/// 邻居的查询上显形,欠定查询照样"通过",很容易漏掉。
///
/// 期望值对应的规则:查询原子度 ≥ 3 才判构型,更少时只要求底物有手性。
#[test]
fn chirality_matching_rebases_the_reference_frame() {
    let targets = ["C[C@H](O)CC", "C[C@@H](O)CC", "CC(O)CC"];
    for (query, want) in [
        // 欠定(度 1):只要求"有手性"
        ("[C@H:1][OH]", [1usize, 1, 0]),
        // 写全(度 3):判构型,而且两个构型互补
        ("[C@:1]([OH])([CH3])[CH2]", [0, 1, 0]),
        ("[C@@:1]([OH])([CH3])[CH2]", [1, 0, 0]),
    ] {
        let q = smarts::parse(query).unwrap_or_else(|e| panic!("{query}: {}", e.render()));
        for (i, smi) in targets.iter().enumerate() {
            let mut m = smiles::parse(smi).unwrap();
            sanitize(&mut m).unwrap();
            let props = MolProps::compute(&m);
            let got = substructure_matches(
                &q,
                &m,
                &props,
                MatchOptions {
                    max_matches: 0,
                    uniquify: true,
                    use_chirality: true,
                },
            )
            .len();
            assert_eq!(got, want[i], "{query} 对 {smi}");
        }
    }
}

/// **首原子的括号氢排第几** —— SMARTS 必须与 SMILES 读成同一件事。
///
/// 规范:`[C@H]` 的括号氢占四元组里"紧跟前一个原子"那一位。首原子没有前一个
/// 原子,氢因此落到**第一位**;而 `C[C@H](N)O` 里第一位是前面那个碳、氢排第二。
/// 两者差一次对换,所以 `[C@@H](C)(N)O` 与 `C[C@H](N)O` 是**同一个**构型。
///
/// # 为什么值得单独钉一条
///
/// 上面那条 `chirality_matching_rebases_the_reference_frame` 的查询里手性原子
/// 虽然也在首位,但**没写括号氢**(`[C@:1](...)`)—— 正好绕开了这一档。
/// 于是"首原子 + 括号氢"这一支在 Rust 侧一条判据都没有,只有
/// `harness/check_smarts_chirality.py` 看得见,而它要 wheel、先前也不在 CI 里。
///
/// # 期望值从哪来:RDKit **自己跟自己**对不上
///
/// 同一批查询,RDKit 2022.09.5 与 2025.09.2 给出**相反**的匹配
/// (`[C@@H](C)(N)O` 在 2025 上匹配 `C[C@H](N)O`,在 2022 上匹配它的对映体);
/// 而两版对**同一串当 SMILES 读**的结果完全一致(都规范成 `C[C@H](N)O`)。
/// 也就是说 2022 的 SMARTS 与它自己的 SMILES 读法自相矛盾,2025 修好了。
/// 本仓库钉 2025.09.2,本实现与它一致。
///
/// 这条判据因此**不诉诸任何外部实现**:它断的是"SMARTS 与 SMILES 同一个读法",
/// 而 SMILES 那一侧由全量差分判据独立守着。
#[test]
fn a_leading_bracket_h_is_read_the_same_in_smarts_and_smiles() {
    // (同一串, 它当 SMILES 读等于谁, 谁是对映体)
    for (query, same, mirror) in [
        ("[C@@H](C)(N)O", "C[C@H](N)O", "C[C@@H](N)O"),
        ("[C@H](C)(N)O", "C[C@@H](N)O", "C[C@H](N)O"),
    ] {
        // 一、先把 SMILES 那一侧钉死 —— 它是这条判据的立足点,
        // 不钉的话下面两条断言可以在"两边一起错"的情况下双双通过
        let a = canon::canonical_smiles(&sanitized(query)).smiles;
        let b = canon::canonical_smiles(&sanitized(same)).smiles;
        let c = canon::canonical_smiles(&sanitized(mirror)).smiles;
        assert_eq!(a, b, "{query} 当 SMILES 读该等于 {same}");
        assert_ne!(b, c, "{same} 与 {mirror} 该是对映体");

        // 二、SMARTS 那一侧必须给出同一个读法
        let opts = MatchOptions {
            max_matches: 0,
            uniquify: true,
            use_chirality: true,
        };
        assert_eq!(
            hits(query, same, opts).len(),
            1,
            "{query} 该匹配 {same} —— 首原子的括号氢排第一位"
        );
        assert_eq!(
            hits(query, mirror, opts).len(),
            0,
            "{query} 不该匹配对映体 {mirror}"
        );
    }
}

/// 双键顺反匹配也要换参照系 —— 与手性同源的一处缺陷。
///
/// `/` 相对键自己的 `begin → end` 朝向,查询与底物的朝向不同。早先直接比原始
/// 方向,于是**答案取决于分子怎么写的**:`F/C=C/F` 与 `C(\F)=C/F` 是同一个
/// 分子,前者匹配、后者不匹配。
///
/// 所以每个构型都用**两种写法**测 —— 只用一种写法测,这个 bug 照样全绿。
#[test]
fn cis_trans_matching_is_independent_of_how_the_molecule_was_written() {
    // 前两个都是反式,中间两个都是顺式,最后一个未指定
    let targets = ["F/C=C/F", "C(\\F)=C/F", "F/C=C\\F", "C(/F)=C/F", "FC=CF"];
    for (query, want) in [
        ("[F]/[C]=[C]/[F]", [1usize, 1, 0, 0, 0]),
        ("[F]/[C]=[C]\\[F]", [0, 0, 1, 1, 0]),
        // 只写一侧 —— 说明不了相对位置,只要求"有方向"
        ("[F]/[C]=[C][F]", [1, 1, 1, 1, 1]),
        // 取代基对不上,一个都不该中
        ("[C]/[C]=[C]/[C]", [0, 0, 0, 0, 0]),
    ] {
        let q = smarts::parse(query).unwrap_or_else(|e| panic!("{query}: {}", e.render()));
        for (i, smi) in targets.iter().enumerate() {
            let mut m = smiles::parse(smi).unwrap();
            sanitize(&mut m).unwrap();
            let props = MolProps::compute(&m);
            let got = substructure_matches(
                &q,
                &m,
                &props,
                MatchOptions {
                    max_matches: 0,
                    uniquify: true,
                    use_chirality: true,
                },
            )
            .len();
            assert_eq!(got, want[i], "{query} 对 {smi}");
        }
    }
}

/// 方向键的**键级**要求就是默认键那一个:单键或芳香键。
///
/// 写成"仅单键"会比不写方向还严,于是模板把方向落在一根**芳香**键上时一处也
/// 匹配不上 —— 而稠环模板里这很常见。取自语料 US08841334B2(2-亚氨基苯并噻唑
/// 的 N-酰化,逆向):模板在环内的 S—c 键上写了 `/`,底物那根键是芳香键。
///
/// 判据卡两头:放宽之后要匹配得上,而**几何约束不能跟着丢**。
#[test]
fn a_direction_written_on_an_aromatic_bond_still_matches() {
    // 底物的方向写在 n—c 上,模板写在 s—c 上 —— 参照原子分处双键两侧
    let substrate = "COCCn1/c(=N/C(=O)CC2CCCC2)sc2cc(F)ccc21";

    assert_eq!(
        count("[#16;a]/[c]", substrate),
        2,
        "芳香键上写方向,键级不该被否掉"
    );
    assert_eq!(
        count("[#16;a][c]", substrate),
        2,
        "不写方向的同一个查询 —— 两者在键级上必须同义"
    );
    assert_eq!(
        count("[#16;a]-[c]", substrate),
        0,
        "写死单键才该落空 —— `-` 与 `/` 不是一回事"
    );

    // 语料原模板的反应物侧:换算过参照系之后几何是一致的
    let tpl = "[#16;a:4]/[c:5](=[N;H0;D2;+0:6]\\[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3]):[#7;a:7]";
    assert_eq!(count(tpl, substrate), 1, "参照系换算过来该匹配");

    // 把 N 那一端的方向翻过来 —— 要的几何变了,一处也不该中
    let flipped = "[#16;a:4]/[c:5](=[N;H0;D2;+0:6]/[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3]):[#7;a:7]";
    assert_eq!(
        count(flipped, substrate),
        0,
        "键级放宽了,几何约束不能跟着丢"
    );
}

/// **总价只有一处实现。** 稠环的稠合位是这条的压力点。
///
/// `omgkit-match` 先前自己写了一份"键级和四舍五入加氢数"。它对**两根**芳香键
/// 的原子与 `omgkit-core` 那份一致,而稠合位有**三根**芳香键:4.5 进位成 5,于是
/// 萘的两个稠合碳被判成五价 —— 任何用 `[v4]` 挑碳的 SMARTS 在稠环上都漏掉稠合位。
///
/// 断的是**性质**,不是当时那份实现的输出:芳香碳一律四价,五价的一个都不该有。
#[test]
fn a_fused_aromatic_carbon_is_tetravalent_like_every_other_carbon() {
    for (smi, carbons) in [
        ("c1ccccc1", 6),                // 苯:全是两根芳香键
        ("c1ccc2ccccc2c1", 10),         // 萘:两个稠合位有三根
        ("c1ccc2c(c1)ccc1ccccc12", 14), // 蒽/菲骨架:四个稠合位
        ("c1ccc2[nH]ccc2c1", 8),        // 吲哚:杂环也一样
    ] {
        assert_eq!(count("[v5;#6]", smi), 0, "{smi}: 有芳香碳被判成五价");
        assert_eq!(
            count("[v4;#6]", smi),
            carbons,
            "{smi}: 四价碳的个数不对 —— 稠合位多半又被算成五价了"
        );
    }
}
