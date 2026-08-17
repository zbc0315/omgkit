//! 写出的判据:**解析 → 写出 → 再解析,得到同一个分子**。
//!
//! 这条判据不依赖任何外部参照 —— 它只说"写出与解析互为逆",而这正是写出
//! 唯一需要满足的性质。逐字节比字面是**规范** SMILES 的判据,要等规范化排序
//! 做完才谈得上;在那之前拿字面去比,比的是排序而不是写出。
//!
//! # 比对要用输出顺序换算
//!
//! 写出会重排原子:输出里的第 `i` 个原子是原分子的 `atom_order[i]`。不换算
//! 就只能比"是否同构",那要跑图同构 —— 既慢,又会把写出的错误和匹配的错误
//! 混成一团。
//!
//! # 键的端点顺序不比
//!
//! 除配位键外,键是无向的,`(a,b)` 与 `(b,a)` 是同一条键。环闭合键的端点
//! 朝向还取决于键级符号写在开环端还是闭合端,是**书写痕迹**而非分子的性质。
//! 配位键例外 —— 它的 `begin` 是给电子的一端,朝向有语义,必须逐条比。

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

use omgkit_core::{AtomFlags, BondOrder, MolBuilder};
use omgkit_io::smiles;

fn corpus(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../harness/corpus")
        .join(name)
}

/// 一条往返失败。
struct Failure {
    smi: String,
    written: String,
    why: String,
}

#[derive(Default)]
struct Stats {
    /// 原语料里能解析的条数
    parsed: usize,
    /// 往返成功的条数
    round_tripped: usize,
    /// 走到过的语法特性,确认判据不是在一堆平凡分子上空过
    with_rings: usize,
    with_aromatic: usize,
    with_fragments: usize,
    with_brackets: usize,
    with_dative: usize,
    /// 需要 `%NN` 形式环标号的分子数
    with_high_ring_labels: usize,
}

/// 逐原子比对(按输出顺序换算回原分子)。
fn compare_atoms(before: &MolBuilder, after: &MolBuilder, order: &[u32]) -> Option<String> {
    for (i, &orig) in order.iter().enumerate() {
        let a = before.atoms()[orig as usize];
        let b = after.atoms()[i];
        let same = a.atomic_num == b.atomic_num
            && a.formal_charge == b.formal_charge
            && a.isotope == b.isotope
            && a.num_explicit_hs == b.num_explicit_hs
            && a.atom_map == b.atom_map
            && a.num_radical_electrons == b.num_radical_electrons
            && a.flags.contains(AtomFlags::AROMATIC) == b.flags.contains(AtomFlags::AROMATIC)
            && a.flags.contains(AtomFlags::NO_IMPLICIT) == b.flags.contains(AtomFlags::NO_IMPLICIT)
            // 四面体手性要往返;配位几何尚不写出,故只在能写的那类上比
            && (!a.chiral_tag.is_tetrahedral() || a.chiral_tag == b.chiral_tag);
        if !same {
            return Some(format!(
                "原子 {i}(原下标 {orig}):写出前 {a:?}\n            写出后 {b:?}"
            ));
        }
    }
    None
}

/// 键集合比对。普通键按无序端点比,配位键按有序端点比。
fn compare_bonds(before: &MolBuilder, after: &MolBuilder, order: &[u32]) -> Option<String> {
    // order[i] = 输出位置 i 对应的原下标;反过来查
    let mut to_out = vec![u32::MAX; before.num_atoms()];
    for (i, &orig) in order.iter().enumerate() {
        to_out[orig as usize] = i as u32;
    }

    let key = |begin: u32, end: u32, order_: BondOrder| {
        // 配位键的 begin 是给电子的一端,方向有语义,不能归一;
        // 其余键无向,端点归一成 (小, 大)
        let flip = order_ != BondOrder::Dative && begin > end;
        if flip {
            (end, begin, order_ as u8)
        } else {
            (begin, end, order_ as u8)
        }
    };

    let want: BTreeSet<_> = before
        .bonds()
        .iter()
        .map(|b| key(to_out[b.begin as usize], to_out[b.end as usize], b.order))
        .collect();
    let got: BTreeSet<_> = after
        .bonds()
        .iter()
        .map(|b| key(b.begin, b.end, b.order))
        .collect();

    if want == got {
        return None;
    }
    let missing: Vec<_> = want.difference(&got).take(5).collect();
    let extra: Vec<_> = got.difference(&want).take(5).collect();
    Some(format!(
        "键集合不一致(端点已换算成输出下标,格式 (起, 终, 键级码))\n\
             丢了:{missing:?}\n    多了:{extra:?}"
    ))
}

fn roundtrip_corpus(path: &Path) -> (Stats, Vec<Failure>) {
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("读不到语料 {}: {e}", path.display()));

    let mut stats = Stats::default();
    let mut bad: Vec<Failure> = Vec::new();

    for line in text.lines() {
        let smi = line.split_whitespace().next().unwrap_or("");
        if smi.is_empty() || smi.starts_with('#') {
            continue;
        }
        let Ok(before) = smiles::parse(smi) else {
            continue; // 非法输入本就不该往返,解析器的报错另有测试守着
        };
        stats.parsed += 1;

        if before.num_bonds() >= before.num_atoms() && before.num_atoms() > 0 {
            stats.with_rings += 1;
        }
        if before
            .atoms()
            .iter()
            .any(|a| a.flags.contains(AtomFlags::AROMATIC))
        {
            stats.with_aromatic += 1;
        }
        if before
            .atoms()
            .iter()
            .any(|a| a.flags.contains(AtomFlags::NO_IMPLICIT))
        {
            stats.with_brackets += 1;
        }
        if before.bonds().iter().any(|b| b.order == BondOrder::Dative) {
            stats.with_dative += 1;
        }

        let w = smiles::write(&before);
        if w.smiles.contains('.') {
            stats.with_fragments += 1;
        }
        if w.smiles.contains('%') {
            stats.with_high_ring_labels += 1;
        }

        let mut fail = |why: String| {
            bad.push(Failure {
                smi: smi.to_string(),
                written: w.smiles.clone(),
                why,
            });
        };

        let after = match smiles::parse(&w.smiles) {
            Ok(m) => m,
            Err(e) => {
                fail(format!("写出的 SMILES 无法再解析:{}", e.kind));
                continue;
            }
        };

        if after.num_atoms() != before.num_atoms() {
            fail(format!(
                "原子数 {} → {}",
                before.num_atoms(),
                after.num_atoms()
            ));
            continue;
        }
        if after.num_bonds() != before.num_bonds() {
            fail(format!(
                "键数 {} → {}",
                before.num_bonds(),
                after.num_bonds()
            ));
            continue;
        }
        if w.atom_order.len() != before.num_atoms() {
            fail(format!(
                "输出顺序长度 {} ≠ 原子数 {}",
                w.atom_order.len(),
                before.num_atoms()
            ));
            continue;
        }

        if let Some(why) = compare_atoms(&before, &after, &w.atom_order) {
            fail(why);
            continue;
        }
        if let Some(why) = compare_bonds(&before, &after, &w.atom_order) {
            fail(why);
            continue;
        }
        stats.round_tripped += 1;
    }

    (stats, bad)
}

fn report(bad: &[Failure], limit: usize) -> String {
    let mut out = format!("\n往返失败 {} 条:\n\n", bad.len());
    for f in bad.iter().take(limit) {
        out.push_str(&format!(
            "  原文:{}\n  写出:{}\n  原因:{}\n\n",
            f.smi, f.written, f.why
        ));
    }
    if bad.len() > limit {
        out.push_str(&format!("  ...(另有 {} 条)\n", bad.len() - limit));
    }
    out
}

/// 冒烟语料。含各种语法陷阱:多片段、方括号、配位键、大环标号、非四面体立体。
#[test]
fn roundtrip_smoke() {
    let (stats, bad) = roundtrip_corpus(&corpus("smoke.smi"));
    assert!(bad.is_empty(), "{}", report(&bad, 15));
    assert!(stats.parsed > 100, "语料只解析出 {} 条", stats.parsed);
    assert_eq!(stats.round_tripped, stats.parsed);

    // 判据必须真的走到这些形状,否则"全绿"只说明语料太平凡
    assert!(stats.with_rings > 0, "语料里没有环");
    assert!(stats.with_aromatic > 0, "语料里没有芳香原子");
    assert!(stats.with_fragments > 0, "语料里没有多片段分子");
    assert!(stats.with_brackets > 0, "语料里没有方括号原子");
    assert!(stats.with_dative > 0, "语料里没有配位键");

    println!(
        "写出往返(冒烟):{} 条全部往返成功;含环 {},含芳香 {},多片段 {},\
         方括号 {},配位键 {},%NN 标号 {}",
        stats.round_tripped,
        stats.with_rings,
        stats.with_aromatic,
        stats.with_fragments,
        stats.with_brackets,
        stats.with_dative,
        stats.with_high_ring_labels
    );
}

/// 大语料(~8800 条)。这一档跑得起量,是写出正确性的主力判据。
#[test]
#[ignore = "语料大,用 cargo test -- --ignored 运行"]
fn roundtrip_large() {
    let (stats, bad) = roundtrip_corpus(&corpus("large.smi"));
    assert!(bad.is_empty(), "{}", report(&bad, 20));
    assert!(stats.parsed > 8000, "语料只解析出 {} 条", stats.parsed);
    assert_eq!(stats.round_tripped, stats.parsed);
    // 语料里同时打开的环从不超过 9 个(标号闭合即回收),所以 `%NN` 分支在
    // 这一档是走不到的。它由 [`ring_labels_beyond_nine`] 专门守着。

    println!(
        "写出往返(大语料):{} 条全部往返成功;含环 {},含芳香 {},多片段 {},\
         方括号 {},%NN 标号 {}",
        stats.round_tripped,
        stats.with_rings,
        stats.with_aromatic,
        stats.with_fragments,
        stats.with_brackets,
        stats.with_high_ring_labels
    );
}

/// 优先级决定起点与分支先后 —— 换一组优先级,字面变了但分子不变。
#[test]
fn priority_changes_the_string_but_not_the_molecule() {
    let smi = "OC(=O)c1ccccc1N";
    let m = smiles::parse(smi).unwrap();
    let n = m.num_atoms();

    let forward: Vec<u32> = (0..n as u32).collect();
    let reverse: Vec<u32> = (0..n as u32).rev().collect();
    let a = smiles::write_with_priority(&m, &forward);
    let b = smiles::write_with_priority(&m, &reverse);

    assert_ne!(a.smiles, b.smiles, "两种优先级应写出不同的字面");
    for w in [&a, &b] {
        let back = smiles::parse(&w.smiles)
            .unwrap_or_else(|e| panic!("{} 无法再解析:{}", w.smiles, e.render()));
        assert!(
            compare_atoms(&m, &back, &w.atom_order).is_none(),
            "{}",
            w.smiles
        );
        assert!(
            compare_bonds(&m, &back, &w.atom_order).is_none(),
            "{}",
            w.smiles
        );
    }
}

/// 造一个"轮毂"分子:中心原子与一条长链上的每个原子都成键。
///
/// 写出时中心原子会一次性打开 `spokes` 个环闭合,标号必然越过 9 —— 真实语料
/// 里同时打开的环从不超过 9 个(闭合即回收标号),这个分支只能这样构造。
fn hub_and_chain(spokes: usize) -> MolBuilder {
    use omgkit_core::BondOrder;
    let mut m = MolBuilder::new();
    let hub = m.add_atom(6);
    // 隔一个原子再连回中心,否则链上第一个原子会与中心成两条键
    let mut prev = m.add_atom(6);
    m.add_bond(hub, prev, BondOrder::Single).expect("端点合法");
    for _ in 0..spokes {
        let a = m.add_atom(6);
        m.add_bond(prev, a, BondOrder::Single).expect("端点合法");
        m.add_bond(hub, a, BondOrder::Single).expect("端点合法");
        prev = a;
    }
    m
}

/// 环闭合标号越过 9 与 99 时的两种写法。
#[test]
fn ring_labels_beyond_nine() {
    for (spokes, marker, what) in [(12usize, "%1", "%NN"), (100, "%(", "%(NNN)")] {
        let m = hub_and_chain(spokes);
        let w = smiles::write(&m);
        assert!(
            w.smiles.contains(marker),
            "{spokes} 条辐条应写出 {what} 形式的标号,实际:{}",
            &w.smiles[..w.smiles.len().min(120)]
        );
        let back = smiles::parse(&w.smiles)
            .unwrap_or_else(|e| panic!("{what} 写出后无法再解析:{}", e.render()));
        assert_eq!(back.num_atoms(), m.num_atoms(), "{what}");
        assert_eq!(back.num_bonds(), m.num_bonds(), "{what}");
        assert!(compare_atoms(&m, &back, &w.atom_order).is_none(), "{what}");
        assert!(compare_bonds(&m, &back, &w.atom_order).is_none(), "{what}");
    }
}

/// 标号闭合后要回收。不回收的话,一条长链上的独立小环会把标号一路推高。
#[test]
fn ring_labels_are_recycled() {
    // 20 个互不重叠的三元环串在一起:任一时刻只有一个环是打开的
    let smi = vec!["C1CC1"; 20].join("");
    let m = smiles::parse(&smi).unwrap();
    let w = smiles::write(&m);
    assert!(
        !w.smiles.contains('%'),
        "20 个不重叠的环不该用到两位数标号 —— 标号没有回收。写出:{}",
        w.smiles
    );
}

/// 方括号里的氢数要取**总数**,不能只取显式那一份。
///
/// 从 SMILES 解析出来的分子里这两件事恰好重合 —— 写在方括号里的原子必然
/// 置了 `NO_IMPLICIT`,隐式氢恒为 0。所以整个语料都验不出这条。
///
/// 一旦分子经过净化、或由程序构造,就不再重合:氢记在 `num_implicit_hs` 里,
/// 而原子可能因为**别的**理由(自由基、立体标记、电荷)需要方括号。这时只写
/// `num_explicit_hs` 会凭空丢掉几个氢,而且丢得悄无声息。
#[test]
fn bracketed_atom_writes_total_hydrogen_count() {
    use omgkit_core::AtomFlags;

    // 甲基自由基:氢是推断来的(未置 NO_IMPLICIT),但自由基逼出了方括号
    let mut m = MolBuilder::new();
    let c = m.add_atom(6);
    {
        let a = m.atom_mut(c).expect("原子存在");
        a.num_implicit_hs = 3;
        a.num_radical_electrons = 1;
    }
    assert!(!m.atoms()[0].flags.contains(AtomFlags::NO_IMPLICIT));

    let w = smiles::write(&m);
    assert_eq!(w.smiles, "[CH3]", "自由基要方括号,且氢数取总数");

    let back = smiles::parse(&w.smiles).expect("应能再解析");
    let a = back.atoms()[0];
    assert_eq!(
        u32::from(a.num_explicit_hs) + u32::from(a.num_implicit_hs),
        3,
        "往返后总氢数不能变"
    );
}

/// 空分子写出空串,不 panic。
#[test]
fn empty_molecule_writes_empty_string() {
    let m = MolBuilder::new();
    let w = smiles::write(&m);
    assert_eq!(w.smiles, "");
    assert!(w.atom_order.is_empty());
}

/// 四面体手性要写出,而且要能原样往返 —— 包括那些**必须翻转**标记的形状。
///
/// 输出会重排邻居(环闭合键推到末尾、分支顺序变),标记必须跟着做一次宇称
/// 换算。不换算的话下面这些用例会写出镜像分子,而且**分子式、键集合全都对**,
/// 只有手性是反的 —— 光比拓扑发现不了。
#[test]
fn tetrahedral_chirality_round_trips() {
    for smi in [
        "N[C@@H](C)C(=O)O",
        "N[C@H](C)C(=O)O",
        "[C@](N)(O)(F)Cl",
        "[C@@H](N)(O)F",
        "[C@@H]1CCCCC1O",
        "C[C@H]1CCCCC1O",
        "O[C@H]1CC[C@@H](N)CC1",
        "C[P@H]C",
        "N[C@H](O)F",
        "OC(=O)[C@@H](N)C",
        // 手性中心上挂两个环闭合
        "C1CC[C@H]2CCCC[C@@H]2C1",
    ] {
        let m = smiles::parse(smi).unwrap_or_else(|e| panic!("{smi}: {}", e.render()));
        let w = smiles::write(&m);
        let back = smiles::parse(&w.smiles)
            .unwrap_or_else(|e| panic!("{smi} → {} 无法再解析: {}", w.smiles, e.render()));
        for (i, &orig) in w.atom_order.iter().enumerate() {
            let before = m.atoms()[orig as usize].chiral_tag;
            let after = back.atoms()[i].chiral_tag;
            assert_eq!(
                before, after,
                "{smi} → {}:原子 {orig} 的手性从 {before:?} 变成了 {after:?}",
                w.smiles
            );
        }
        assert!(
            compare_bonds(&m, &back, &w.atom_order).is_none(),
            "{smi} → {}",
            w.smiles
        );
    }
}

/// 双键方向键会写出,而且只写**携带信息**的那些。
#[test]
fn informative_direction_bonds_are_written() {
    for smi in ["F/C=C/F", "F/C=C\\F", "C/C=C/C=C/C", "CC(/F)=C(\\F)C"] {
        let m = smiles::parse(smi).unwrap();
        let w = smiles::write(&m).smiles;
        assert!(
            w.contains('/') || w.contains('\\'),
            "{smi} 写成了 {w},方向键丢了"
        );
    }
}

/// 不携带信息的方向不写。
///
/// 这不只是"输出更干净":颜色细化看不见键方向,一条噪声方向能打破细化分辨
/// 不出的对称性,规范 SMILES 就不再随重排恒定 —— 下面第一条 `C/1CCCCC1`
/// 正是这种分子。
#[test]
fn noise_direction_bonds_are_dropped() {
    for (smi, why) in [
        ("C/1CCCCC1", "根本没有双键"),
        ("F/C=CF", "只有一侧有方向"),
        ("F/C=C(F)F", "双键一端挂着两个相同的取代基"),
        ("C/C=C(/C)C", "同上,两个甲基"),
    ] {
        let m = smiles::parse(smi).unwrap();
        let w = smiles::write(&m).smiles;
        assert!(
            !w.contains('/') && !w.contains('\\'),
            "{smi} 写成了 {w},但这里的方向是噪声:{why}"
        );
    }
}

/// 配位几何(`@SP`/`@TB`/`@OH`)与丙二烯轴手性尚未写出。
///
/// 这条断言把"尚未实现"钉成显式事实 —— 免得将来实现了却没人发现输出里
/// 其实一直没有这些信息。
#[test]
fn coordination_and_axial_stereo_are_not_written_yet() {
    for smi in ["[Pt@SP1](Cl)(Cl)(N)N", "N[C@AL1]=C=C(O)F"] {
        let m = smiles::parse(smi).unwrap();
        let w = smiles::write(&m).smiles;
        assert!(
            !w.contains('@'),
            "{smi} 写成了 {w},但这类立体信息的写出尚未实现"
        );
    }
}

/// 带显式氢的原子必须写在方括号里 —— 简写形式没地方放它们。
///
/// 光靠 [`AtomFlags::NO_IMPLICIT`] 判断是不够的:净化会清掉那个标志,同时把
/// 氢挪进 `num_explicit_hs`。于是"净化之后写出"会把吡咯型氮的 `[nH]` 写成
/// 裸 `n`,氢凭空消失,写出的串连凯库勒化都做不到。
///
/// 这个缺口在两处都看不见 —— 不净化的往返测试里标志还在,L2 的差分比的是
/// 分子对象的字段而不是写出的字符串。实测:8839 条语料净化后写出,633 条
/// 因此坏掉。所以这里直接构造"标志已清、显式氢还在"的状态来守。
#[test]
fn explicit_hydrogens_force_brackets() {
    let mut m = smiles::parse("c1cc[nH]c1").unwrap();
    // 模拟净化后的状态:标志清掉,氢留在 num_explicit_hs 里
    for i in 0..m.num_atoms() as u32 {
        if let Some(a) = m.atom_mut(i) {
            a.flags.remove(omgkit_core::AtomFlags::NO_IMPLICIT);
        }
    }
    let w = smiles::write(&m).smiles;
    assert!(
        w.contains("[nH]"),
        "写成了 {w} —— 显式氢丢了,而这个串外部实现连凯库勒化都做不了"
    );
}
