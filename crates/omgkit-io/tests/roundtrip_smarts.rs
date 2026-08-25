//! SMARTS 写出的判据:**写出幂等且规模守恒**。
//!
//! # 为什么不比表达式树
//!
//! 写出会重排原子编号(生成树的遍历顺序与原编号无关),再解析得到的
//! `QueryMol` 拓扑与原来同构但编号不同。直接比树会把"编号变了"报成"写错了"。
//!
//! 所以判据是:
//!
//! - **规模守恒** —— 原子数与键数不变
//! - **映射号守恒** —— 原子映射号(`:n`)的多重集不变
//! - **写出幂等** —— `write(parse(write(q))) == write(q)`
//!
//! 幂等这条比看上去强:写出的结果一旦与它自己的再解析对不上,说明这一趟里
//! 有信息被丢掉或被凭空添上。
//!
//! # 但幂等对"整批丢掉"是**瞎的**
//!
//! 映射号被写出器整个丢掉的话:第一趟写出没有 `:n`,再解析也没有,第二趟写出
//! 还是没有 —— **幂等成立、规模守恒成立**,而反应模板已经废了(产物构建全靠
//! 映射号对位)。所以映射号要单独数。
//!
//! 同一个洞外部判据也堵不上:`harness/check_smarts_write.py` 比的是**匹配语义**,
//! 而映射号不参与匹配;`harness/corpus/smarts.txt` 的 776 条模式里**一条映射号
//! 都没有**。真正带映射号的语料是 `harness/corpus/reactions.txt`,
//! 见 `reaction_templates_keep_their_atom_maps`。
//!
//! 语义等价(两个 SMARTS 匹配到同样的东西)要外部实现当判官,那是另一档。

use omgkit_io::smarts::{self, AtomExpr, AtomPrim, QueryMol};

/// 一个查询里所有原子映射号,**排好序的多重集**。
///
/// 写出会重排原子编号,所以只能比多重集,不能逐位比。递归 SMARTS
/// (`$(...)`)里面的也要收 —— 那里同样写得下 `:n`。
fn atom_maps(q: &QueryMol) -> Vec<u16> {
    fn walk(e: &AtomExpr, out: &mut Vec<u16>) {
        match e {
            AtomExpr::Prim(AtomPrim::AtomMap(m)) => out.push(*m),
            AtomExpr::Prim(AtomPrim::Recursive(sub)) => {
                for a in &sub.atoms {
                    walk(a, out);
                }
            }
            AtomExpr::Prim(_) => {}
            AtomExpr::Not(inner) => walk(inner, out),
            AtomExpr::And(parts) | AtomExpr::Or(parts) => {
                for p in parts {
                    walk(p, out);
                }
            }
        }
    }
    let mut out = Vec::new();
    for a in &q.atoms {
        walk(a, &mut out);
    }
    out.sort_unstable();
    out
}

fn check(src: &str) -> Result<(), String> {
    let q = smarts::parse(src).map_err(|e| format!("解析失败: {}", e.render()))?;
    let w1 = smarts::write(&q);
    let q2 = smarts::parse(&w1).map_err(|e| format!("写出 {w1} 之后解析失败: {}", e.render()))?;
    let w2 = smarts::write(&q2);
    if q2.num_atoms() != q.num_atoms() || q2.num_bonds() != q.num_bonds() {
        return Err(format!(
            "规模变了:{}原子/{}键 -> {}原子/{}键(写出 {w1})",
            q.num_atoms(),
            q.num_bonds(),
            q2.num_atoms(),
            q2.num_bonds()
        ));
    }
    // **映射号要单独数。** 整批丢掉时幂等与规模守恒双双成立,见模块文档。
    let (m1, m2) = (atom_maps(&q), atom_maps(&q2));
    if m1 != m2 {
        return Err(format!("原子映射号变了:{m1:?} -> {m2:?}(写出 {w1})"));
    }
    if w1 != w2 {
        return Err(format!("写出不幂等:{w1} -> {w2}"));
    }
    Ok(())
}

/// 手工用例,覆盖各档语法。
#[test]
fn handpicked_patterns_round_trip() {
    for src in [
        "[C]",
        "[CH3]",
        // 三档优先级都要走到 —— `;` 与 `&` 选错会改语义
        "[C,N;H1]",
        "[C,N&H1]",
        "[!C;!N]",
        "[#6]",
        "[c]",
        "[a]",
        "[A]",
        "CCO",
        "c1ccccc1",
        "[C:1][OH:2]",
        "[$(C=O)]",
        "[!$(C=O)]",
        "[C;R2]",
        "[C]~[N]",
        "C1CC1",
        "C1CCC2CCCCC2C1",
        "[C](=O)[OH]",
        "[N+]",
        "[N-2]",
        "[13C]",
        "[C@H]",
        "*",
        "CC.NN",
        "[C]->[O]",
        "[C]<-[O]",
        "[C]@[C]",
        "[C]!@[C]",
    ] {
        check(src).unwrap_or_else(|e| panic!("{src}: {e}"));
    }
}

/// `;` 与 `&` 的选择是**语义相关**的 —— 原子表达式里没有括号可以分组。
///
/// `And[Or[C,N], H1]` 写成 `&` 会变成 `C,(N&H1)`,那是另一个查询。
#[test]
fn and_operator_choice_preserves_grouping() {
    let a = smarts::parse("[C,N;H1]").unwrap();
    let b = smarts::parse("[C,N&H1]").unwrap();
    assert_ne!(a.atoms[0], b.atoms[0], "两者本来就该不同");
    assert_eq!(smarts::write(&a), "[C,N;H1]");
    assert_eq!(smarts::write(&b), "[C,N&H1]");
}

/// 默认键(单键或芳香键)不写符号 —— 写成 `-,:` 语义没错但每根键多两个字符。
#[test]
fn default_bonds_are_omitted() {
    let q = smarts::parse("CCO").unwrap();
    let w = smarts::write(&q);
    assert!(!w.contains(','), "默认键被展开成了析取式:{w}");
    assert_eq!(w, "[C][C][O]");
}

/// 反应写成固定的三段式,试剂段为空也要留两个 `>`。
#[test]
fn reactions_keep_three_sections() {
    for src in [
        "[C:1][OH:2]>>[C:1][Cl:2]",
        "[C:1][OH:2].[N:3]>>[C:1][N:3]",
        "[C:1]=[O:2]>[Pd]>[C:1][O:2]",
    ] {
        let r = smarts::parse_reaction(src).unwrap();
        let w = smarts::write_reaction(&r);
        assert_eq!(w.matches('>').count(), 2, "{src} 写成了 {w},段数不对");
        let r2 = smarts::parse_reaction(&w)
            .unwrap_or_else(|e| panic!("{src} -> {w} 解析不了: {}", e.render()));
        assert_eq!(r2.reactants.len(), r.reactants.len(), "{src} -> {w}");
        assert_eq!(r2.agents.len(), r.agents.len(), "{src} -> {w}");
        assert_eq!(r2.products.len(), r.products.len(), "{src} -> {w}");
        assert_eq!(smarts::write_reaction(&r2), w, "反应写出不幂等");
    }
}

/// **反应模板的映射号必须原样活下来** —— 语料是 `reactions.txt`。
///
/// # 为什么单开一条,不并进上面的语料
///
/// `harness/corpus/smarts.txt` 的 776 条模式里**一条映射号都没有**(实测)。
/// 也就是说上面 `corpus_round_trips` 里那条映射号判据在那份语料上是**空过**的:
/// 数出来两边都是空多重集,当然相等。真正带映射号的语料只有反应模板。
///
/// 而这一档丢了就是灾难:产物构建全靠映射号对位,丢掉之后模板还"写得出来、
/// 解析得回、幂等",只是它描述的不再是同一个反应。
///
/// # 三段都要走
///
/// 反应写成 `反应物 > 试剂 > 产物`,三段各自是一串 `.` 分开的模式。
/// 只测反应物那一段的话,写出器在产物段上的回归照样全绿。
#[test]
fn reaction_templates_keep_their_atom_maps() {
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../harness/corpus/reactions.txt");
    let text =
        std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("读不到 {}: {e}", path.display()));

    let mut n_rxn = 0usize;
    let mut n_maps = 0usize;
    let mut bad: Vec<String> = Vec::new();
    for line in text.lines() {
        let src = line.trim();
        if src.is_empty() || src.starts_with('#') {
            continue;
        }
        let r = match smarts::parse_reaction(src) {
            Ok(r) => r,
            Err(e) => {
                bad.push(format!("  {src}\n     解析不了:{}", e.render()));
                continue;
            }
        };
        n_rxn += 1;
        let before: Vec<Vec<u16>> = [&r.reactants, &r.agents, &r.products]
            .into_iter()
            .map(|part| {
                let mut v: Vec<u16> = part.iter().flat_map(atom_maps).collect();
                v.sort_unstable();
                v
            })
            .collect();
        n_maps += before.iter().map(Vec::len).sum::<usize>();

        let w = smarts::write_reaction(&r);
        let Ok(r2) = smarts::parse_reaction(&w) else {
            bad.push(format!("  {src}\n     写成 {w} 之后解析不了"));
            continue;
        };
        let after: Vec<Vec<u16>> = [&r2.reactants, &r2.agents, &r2.products]
            .into_iter()
            .map(|part| {
                let mut v: Vec<u16> = part.iter().flat_map(atom_maps).collect();
                v.sort_unstable();
                v
            })
            .collect();
        // **逐段比。** 合起来比的话,把一个映射号从反应物挪到产物照样"守恒"。
        for (i, seg) in ["反应物", "试剂", "产物"].into_iter().enumerate() {
            if before[i] != after[i] {
                bad.push(format!(
                    "  {src}\n     写成 {w}\n     {seg}段的映射号 {:?} -> {:?}",
                    before[i], after[i]
                ));
            }
        }
    }

    assert!(
        bad.is_empty(),
        "反应模板的映射号没活下来:\n{}",
        bad.join("\n")
    );
    // **区分力闸。** 语料换成一份没有映射号的,上面每一条断言都会在空多重集上
    // 成立 —— 那是这条判据最会骗人的一种绿。数出来钉死。
    assert_eq!(n_rxn, 20, "反应语料的条数变了 —— 先确认是有意的再改这个数");
    assert_eq!(
        n_maps, 96,
        "语料里的映射号总数变了 —— 先确认是有意的再改这个数;\
         这个数掉到 0 时上面每一条断言都会空过"
    );
}

/// 全量 SMARTS 语料。
#[test]
fn corpus_round_trips() {
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../harness/corpus/smarts.txt");
    let text =
        std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("读不到 {}: {e}", path.display()));
    let mut n = 0usize;
    let mut bad: Vec<String> = Vec::new();
    for line in text.lines() {
        let src = line.trim();
        if src.is_empty() || src.starts_with('#') {
            continue;
        }
        // 语料里有故意的非法输入,解析器另有测试守着,这里只管解析得了的
        if smarts::parse(src).is_err() {
            continue;
        }
        n += 1;
        if let Err(e) = check(src) {
            bad.push(format!("  {src}\n     {e}"));
        }
    }
    assert!(n > 300, "语料只走到 {n} 条,判据几乎是空过的");
    assert!(
        bad.is_empty(),
        "SMARTS 往返失败 {} 条:\n{}",
        bad.len(),
        bad.iter().take(15).cloned().collect::<Vec<_>>().join("\n")
    );
    println!("SMARTS 写出往返通过:{n} 条");
}

/// 手性标记的参照系:解析归一到存储序,写出换回书写序,两者必须互逆。
///
/// # 为什么要有这一档
///
/// 查询里的 `@` 相对**书写顺序**,而下游(匹配、产物构建)一律按存储序读它。
/// 中间差着两项:环闭合键在串里紧跟原子、存储时却追加到末尾;括号里的氢在
/// 原子是片段首原子时排第一位、否则排第二位。少补任何一项,查询就会匹配到
/// 镜像分子 —— 而原子数、键数、连通性全都对。
///
/// 判据是**写出幂等**:`写出(解析(s))` 与 `写出(解析(写出(解析(s))))` 必须
/// 逐字节相同。不能拿"两次解析的树相等"当判据 —— 写出会重排原子,拓扑本就
/// 不同,那样比永远不过。
///
/// 幂等抓得住这个缺陷:写出若不做逆向换算,一来一回就连翻两次,第二次写出的
/// 标记与第一次相反。
#[test]
fn chirality_survives_parse_write_parse() {
    for s in [
        // 首原子 + 括号氢,0/1/2 个环闭合
        "[C@H](N)(O)F",
        "[C@@H](N)(O)F",
        "[C@H]1(C)CCC[C@@H]1O",
        "[C@@H]1(C)CCC[C@@H]1O",
        "[C@H]12CC[C@H]3CCCC[C@H]1CC23",
        // 非首原子
        "N[C@H](O)F",
        "C[C@H]1CC[C@@H](O)C1",
        "O[C@H]1CCCC[C@@H]1N",
        // 无括号氢
        "[C@](N)(O)(F)C",
        "C[C@]1(N)CC[C@@H](O)C1",
    ] {
        let w1 = smarts::write(&smarts::parse(s).unwrap_or_else(|e| panic!("{s}: {}", e.render())));
        let w2 = smarts::write(
            &smarts::parse(&w1).unwrap_or_else(|e| panic!("{s} 写成 {w1} 后读不回:{}", e.render())),
        );
        assert_eq!(
            w1, w2,
            "{s}:写出不幂等,{w1} 再走一遍成了 {w2} —— \
             解析与写出的宇称换算没有互逆,查询会匹配到镜像分子"
        );
    }
}
