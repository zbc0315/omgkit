//! 桥环骨架的预存坐标表。
//!
//! # 为什么要有这张表
//!
//! 桥环、笼状体系在平面上没有好解,[`rings::layout_local`](crate::rings) 只能
//! 退化到弹簧松弛。松弛是局部下降,落到哪个极小全看初值 —— 现在的 5 个初值
//! **本身就常常给出自交的解**:实测最常见的 8 个骨架里有 5 个是自交的,包括
//! 双环[2.2.2]辛烷和金刚烷这种最基本的形状。
//!
//! 而这类骨架**极其集中**:全量语料 177 处退化只有 44 种骨架,最常见的 20 种
//! 覆盖 **84.2%**。所以把好解算一次存下来,是划算的。
//!
//! # 坐标从哪来 —— 不是手画的
//!
//! 由 `rings.rs` 里的 `regenerate_templates` 生成:对每个骨架跑**带扰动**
//! 的多起点松弛(基础两万次,仍自交的接着搜到四十万),按**现成的** `Quality`
//! (自交数、最大键长偏差、量化坐标序列)
//! 挑最好的那个。判优的口径与运行时完全一样,只是搜得久得多 ——
//! **这张表就是一次昂贵搜索的缓存**,不是另一套标准。
//!
//! 生成脚本与产物一起进版本库,谁都能重跑一遍核对,与 `harness/gen_elements.py`
//! 生成 `element_data.rs` 是同一个路子。
//!
//! # 命中之后仍然算"退化"
//!
//! 模板给的坐标键长并不严格全等(松弛出来的本就不等)。命中模板只是**换了一个
//! 更好的退化解**,不是把问题解决了,所以 [`Degradation`](crate::rings::Degradation)
//! 照报不误 —— 下游据此拒绝渲染或人工介入的逻辑不受影响。

use std::collections::BTreeMap;

use omgkit_core::{BondOrder, MolBuilder};

use crate::geom::Point2;

/// 环系骨架的指纹:把环系原子抠出来单独成一个分子,取它的规范 SMILES。
///
/// **取代基剥掉、原子一律当碳、键级一律当单键** —— 模板管的是形状,不是化学。
/// 一个甾体骨架无论挂什么取代基、环上是碳还是氮,该摆成同一个形状。
///
/// 原子按**父分子的规范秩**加进去,所以指纹与写法无关。
pub(crate) fn skeleton_of(mol: &MolBuilder, atoms: &[u32], ranks: &[u32]) -> Option<String> {
    let mut order: Vec<u32> = atoms.to_vec();
    order.sort_by_key(|a| (ranks[*a as usize], *a));

    let mut b = MolBuilder::new();
    let mut map: BTreeMap<u32, u32> = BTreeMap::new();
    for a in &order {
        map.insert(*a, b.add_atom(6));
    }
    for bd in mol.bonds() {
        if let (Some(x), Some(y)) = (map.get(&bd.begin), map.get(&bd.end)) {
            b.add_bond(*x, *y, BondOrder::Single).ok()?;
        }
    }
    omgkit_chem::pipeline::sanitize(&mut b).ok()?;
    Some(omgkit_io::canon::canonical_smiles(&b).smiles)
}

/// 查表的结果。**"没命中"要分成两种** —— 它们指向完全不同的动作。
///
/// 先前只有一个 `bool`,把后两种混成一档。而实测全量语料里没命中的那几例
/// **全是 `NoFingerprint`**(骨架全碳化之后度数超 4,`sanitize` 不过)——
/// 报"该补进语料"是条走不通的路,100% 的情形都指错了方向。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Status {
    /// 表里有,坐标就是表给的
    Hit,
    /// 指纹算出来了,表里没有。
    /// **这个骨架该补进 `harness/corpus/bridged.smi` 再重跑生成器。**
    NotInTable,
    /// 指纹根本算不出来 —— 抠出来的骨架自己 sanitize 不过。补语料没用,
    /// 要查的是那个分子本身。
    NoFingerprint,
}

/// 查表:这个环系有没有预存的坐标,以及没有的话是哪一种没有。
///
/// 返回的坐标按**父分子的原子编号**给出,已经对上号了。
///
/// **状态与坐标一起返回**,别让调用方为了知道"命中没有"再调一次 —— 这里面是
/// 建分子 + sanitize + 规范化,每个桥环系统、每种规范、审计里每种写法都要付。
/// 临时顶替表里某一条的坐标。**只给离线的模板生成器用。**
///
/// 生成器要问"把这组坐标装进去之后,真实分子画出来好不好",而 `generate` 会查
/// 这张表 —— 表正是它在生成的东西。这个参数把那层循环拆开。
///
/// **为什么是穿参数,不是全局状态。** 试过 `cfg(test)` + `thread_local`:
/// 它有隐藏状态、panic 之后覆盖会留在线程里、将来 `generate` 内部要是并行就坏。
/// 穿参数多改四个 `pub(crate)` 签名,但没有这些问题,而且看得见。
pub(crate) type Override<'a> = Option<(&'a str, &'a [(f64, f64)])>;

/// 同 [`lookup`],但可以临时顶替表里某一条。见 [`Override`]。
///
/// **状态与坐标一起返回**,别让调用方为了知道"命中没有"再调一次 —— 那里面是
/// 建分子 + sanitize + 规范化,每个桥环系统、每种规范、审计里每种写法都要付。
pub(crate) fn lookup_with(
    mol: &MolBuilder,
    atoms: &[u32],
    ranks: &[u32],
    over: Override<'_>,
) -> (Option<BTreeMap<u32, Point2>>, Status) {
    let Some(skel) = skeleton_of(mol, atoms, ranks) else {
        return (None, Status::NoFingerprint);
    };
    // **给了覆盖时,整张表都被遮住。**
    //
    // 生成器要问"把这组坐标装进去之后,真实分子画出来好不好"。若只遮住匹配的
    // 那一条、其余仍读 `TABLE`,而打分分子里恰好还有**另一个**桥环骨架,那一条
    // 就会读到**正在生成的那张表** —— 生成器就不再是语料的纯函数,
    // 「把 `TABLE` 清空重跑逐字节相同」这条验收当场作废。
    //
    // 遮住整张表,这件事就成了**结构保证**,不再依赖"语料里碰巧没有这种分子"。
    let coords = match over {
        Some((k, v)) if k == skel => v,
        Some(_) => return (None, Status::NotInTable),
        None => match TABLE.iter().find(|(k, _)| *k == skel).map(|(_, v)| *v) {
            Some(v) => v,
            None => return (None, Status::NotInTable),
        },
    };

    // 重建同一个骨架,拿它自己的规范秩去对坐标 —— 存的时候就是按这个存的。
    // 规范秩与存储序无关,所以两边算出来必然一致。
    //
    // 下面任何一步失败都当"表里没有"退回松弛:表在,只是对不上号。
    let mut order: Vec<u32> = atoms.to_vec();
    order.sort_by_key(|a| (ranks[*a as usize], *a));
    let mut b = MolBuilder::new();
    let mut back: Vec<u32> = Vec::with_capacity(order.len());
    for a in &order {
        b.add_atom(6);
        back.push(*a);
    }
    let idx: BTreeMap<u32, u32> = order
        .iter()
        .enumerate()
        .map(|(i, a)| (*a, u32::try_from(i).unwrap_or(0)))
        .collect();
    for bd in mol.bonds() {
        if let (Some(x), Some(y)) = (idx.get(&bd.begin), idx.get(&bd.end)) {
            if b.add_bond(*x, *y, BondOrder::Single).is_err() {
                return (None, Status::NotInTable);
            }
        }
    }
    if omgkit_chem::pipeline::sanitize(&mut b).is_err() {
        return (None, Status::NotInTable);
    }
    let skel_ranks = omgkit_io::canon::canonical_ranks(&b);
    if skel_ranks.len() != coords.len() {
        return (None, Status::NotInTable); // 表坏了 —— 不猜,退回松弛
    }

    let mut out = BTreeMap::new();
    for (i, parent) in back.iter().enumerate() {
        let Some(r) = skel_ranks.get(i) else {
            return (None, Status::NotInTable);
        };
        let Some((x, y)) = coords.get(*r as usize) else {
            return (None, Status::NotInTable);
        };
        out.insert(*parent, Point2::new(*x, *y));
    }
    (Some(out), Status::Hit)
}

/// 查表:这个环系有没有预存的坐标,以及没有的话是哪一种没有。
///
/// 返回的坐标按**父分子的原子编号**给出,已经对上号了。
pub fn lookup(
    mol: &MolBuilder,
    atoms: &[u32],
    ranks: &[u32],
) -> (Option<BTreeMap<u32, Point2>>, Status) {
    lookup_with(mol, atoms, ranks, None)
}

include!("templates_data.rs");

#[cfg(test)]
mod tests {
    use super::*;
    use crate::style::Style;

    fn prep(smi: &str) -> MolBuilder {
        let mut m = omgkit_io::smiles::parse(smi).unwrap();
        omgkit_chem::pipeline::sanitize(&mut m).unwrap();
        m
    }

    /// 表里有、语料里也常见的桥环分子
    const HIT: [&str; 4] = [
        "C1CC2CCC1CC2",                                     // 双环[2.2.2]辛烷
        "C1C2CC3CC1CC(C2)C3",                               // 金刚烷
        "CC1(C)[C@@H]2CC[C@@]1(C)C(=O)C2",                  // 樟脑
        "CN1[C@H]2CC[C@@H]1C[C@@H](C2)OC(=O)C(CO)c1ccccc1", // 阿托品
    ];

    /// 不给覆盖时,`lookup_with` 必须与表里那一行**逐位相同**。
    ///
    /// 这条守的是"加了覆盖这条路,没覆盖时的行为一个字节都没变"。
    /// **不能写成"比较有/无覆盖机制两个版本"** —— 同一个测试二进制里编不出
    /// "没有覆盖机制"的那一版,那种判据根本落不了地。
    #[test]
    fn without_an_override_the_table_is_what_comes_back() {
        let mut checked = 0usize;
        for (skel, coords) in TABLE {
            let mut m = omgkit_io::smiles::parse(skel).expect("表里的骨架该能解析");
            omgkit_chem::pipeline::sanitize(&mut m).expect("该能 sanitize");
            let ranks = omgkit_io::canon::canonical_ranks(&m);
            let atoms: Vec<u32> = (0..u32::try_from(m.num_atoms()).unwrap()).collect();
            let (got, st) = lookup_with(&m, &atoms, &ranks, None);
            let Some(got) = got else { continue };
            assert_eq!(st, Status::Hit);
            // 坐标按骨架自己的规范秩存,取回来逐位比
            for (a, p) in &got {
                let (x, y) = coords[ranks[*a as usize] as usize];
                assert!(
                    p.x.to_bits() == x.to_bits() && p.y.to_bits() == y.to_bits(),
                    "{skel} 的原子 {a}:取回 {p:?},表里是 ({x}, {y})"
                );
            }
            checked += 1;
        }
        assert!(checked >= 40, "只验到 {checked} 条,判据太弱");
    }

    /// 给了覆盖时,拿回来的必须是覆盖的那一份,不是表里的。
    #[test]
    fn an_override_really_replaces_that_one_row() {
        let skel = "C1C2CCC1CC2";
        let mut m = omgkit_io::smiles::parse(skel).expect("该能解析");
        omgkit_chem::pipeline::sanitize(&mut m).expect("该能 sanitize");
        let ranks = omgkit_io::canon::canonical_ranks(&m);
        let atoms: Vec<u32> = (0..u32::try_from(m.num_atoms()).unwrap()).collect();
        let n = m.num_atoms();

        // 一组一眼认得出来的坐标:第 i 个原子放在 (i, -i)
        #[allow(clippy::cast_precision_loss)]
        let fake: Vec<(f64, f64)> = (0..n).map(|i| (i as f64, -(i as f64))).collect();
        let (got, st) = lookup_with(&m, &atoms, &ranks, Some((skel, &fake)));
        assert_eq!(st, Status::Hit);
        let got = got.expect("装了覆盖就该拿得到");
        for (a, p) in &got {
            let (x, y) = fake[ranks[*a as usize] as usize];
            assert!(
                (p.x - x).abs() < 1e-12 && (p.y - y).abs() < 1e-12,
                "{skel} 的原子 {a}:拿回 {p:?},覆盖里是 ({x}, {y})"
            );
        }

        // **覆盖的是另一条骨架时,这一条也读不到表。**
        //
        // 这是有意的:给了覆盖就遮住整张表。生成器打分时,分子里若还有别的
        // 桥环骨架,那一条读到的就会是**正在生成的那张表** —— 纯函数性当场破。
        // 遮住整张表把这件事变成结构保证,不再依赖"语料里碰巧没有这种分子"。
        let (other, st2) = lookup_with(&m, &atoms, &ranks, Some(("C1CC2CCC1CC2", &fake)));
        assert!(other.is_none(), "覆盖了别的骨架,这一条不该还能读到表");
        assert_eq!(st2, Status::NotInTable);
        // 而不给覆盖时照旧
        let (plain, st3) = lookup_with(&m, &atoms, &ranks, None);
        assert!(
            plain.is_some() && st3 == Status::Hit,
            "不给覆盖时该读得到表"
        );
    }

    /// 表里的坐标自己有多少处自交。
    ///
    /// 生成器打出来的 `// 出现 N 次,自交 M` 只是**注释** —— 人改一行坐标它不会
    /// 变红。这条判据把坐标重新算一遍。
    ///
    /// # 这个数**不是**优化目标,只是粗线条的回归闸
    ///
    /// 生成器现在按**整分子**打分挑候选(见 `rings.rs` 的 `score_on_molecules`),
    /// 光骨架的自交数只是平局兜底的一部分。两者会背离,而且这次正是背离的:
    ///
    /// | | 骨架自交总数 | 全量语料的整分子键交叉 |
    /// |---|---:|---:|
    /// | 按骨架挑 | 8 | 62 |
    /// | 按整分子挑 | **9** | **40** |
    ///
    /// 多一处骨架自交,换掉 22 处整分子交叉 —— **这恰恰证明骨架自交是错的
    /// 代理指标**。所以这里的上界跟着新口径走,别把它当成"越小越好"。
    #[test]
    fn the_stored_coordinates_do_not_cross_more_than_they_used_to() {
        let mut total = 0usize;
        let mut worst: Vec<(usize, &str)> = Vec::new();
        for (skel, coords) in TABLE {
            let mut m = omgkit_io::smiles::parse(skel).expect("表里的骨架该能解析");
            omgkit_chem::pipeline::sanitize(&mut m).expect("表里的骨架该能 sanitize");
            let ranks = omgkit_io::canon::canonical_ranks(&m);
            assert_eq!(
                ranks.len(),
                coords.len(),
                "{skel}:表里 {} 组坐标,骨架 {} 个原子",
                coords.len(),
                ranks.len()
            );
            // 坐标是按骨架自己的规范秩存的
            let pos: Vec<Point2> = (0..m.num_atoms())
                .map(|i| {
                    let (x, y) = coords[ranks[i] as usize];
                    Point2::new(x, y)
                })
                .collect();
            let segs: Vec<(Point2, Point2)> = m
                .bonds()
                .iter()
                .map(|b| (pos[b.begin as usize], pos[b.end as usize]))
                .collect();
            let mut cross = 0usize;
            for (k, (u1, v1)) in segs.iter().enumerate() {
                for (u2, v2) in &segs[k + 1..] {
                    if crate::geom::segments_cross(*u1, *v1, *u2, *v2) {
                        cross += 1;
                    }
                }
            }
            if cross > 0 {
                worst.push((cross, skel));
            }
            total += cross;
        }
        worst.sort_unstable();
        // 现值:7 条自交,总数 9。上界是回归闸,不是优化目标 —— 见本判据的文档。
        assert!(
            total <= 9,
            "表里的自交总数涨到了 {total},还剩 {} 条自交:{worst:?}",
            worst.len()
        );
    }

    #[test]
    fn the_table_is_actually_used() {
        // 一条都命中不了的话,下面那些判据全是空过的。
        let mut hits = 0usize;
        for smi in HIT {
            let m = prep(smi);
            let ranks = omgkit_io::canon::canonical_ranks(&m);
            let rs = omgkit_chem::sssr::ring_set(&m);
            for sys in crate::rings::group(&omgkit_chem::rings::fused_ring_systems(&m), &rs) {
                if lookup(&m, &sys.atoms, &ranks).0.is_some() {
                    hits += 1;
                }
            }
        }
        assert!(hits > 0, "这批桥环分子一个都没命中模板 —— 表白建了");
    }

    #[test]
    fn a_template_hit_does_not_depend_on_how_the_molecule_was_written() {
        // 指纹按父分子的规范秩建骨架,坐标按**骨架自己的规范秩**存取。两处都
        // 与存储序无关 —— 沾上一处,同一个分子换种写法就会摆成另一个样子,
        // 而**写法无关是本库的头号契约**。
        for smi in HIT {
            for style in &Style::ALL {
                let m = prep(smi);
                let n = m.num_atoms();
                let want = crate::generate(&m, style);
                let priority: Vec<u32> = (0..n)
                    .map(|i| u32::try_from(n - 1 - i).expect("原子数超出 u32"))
                    .collect();
                let w = omgkit_io::smiles::write_with_priority(&m, &priority);
                let Some(m2) = omgkit_io::smiles::parse(&w.smiles)
                    .ok()
                    .and_then(|mut x| omgkit_chem::pipeline::sanitize(&mut x).ok().map(|()| x))
                else {
                    continue;
                };
                if omgkit_io::canon::canonical_smiles(&m).smiles
                    != omgkit_io::canon::canonical_smiles(&m2).smiles
                {
                    continue;
                }
                let got = crate::generate(&m2, style);
                let q = |c: &[Point2]| {
                    let mut v: Vec<(i64, i64)> = c
                        .iter()
                        .map(|p| ((p.x * 1e4).round() as i64, (p.y * 1e4).round() as i64))
                        .collect();
                    v.sort_unstable();
                    v
                };
                assert_eq!(
                    q(&want.coords),
                    q(&got.coords),
                    "[{}] {smi}:换成 {} 之后画出来不一样了",
                    style.name,
                    w.smiles
                );
            }
        }
    }

    #[test]
    fn a_templated_system_still_reports_itself_degraded() {
        // **命中模板不等于问题解决了。** 模板给的键长并不严格全等,它只是换了
        // 一个更好的退化解。悄悄把 `degraded` 清掉,下游就以为这张图可以放心
        // 用了 —— 那正是这个库最不该做的事。
        for smi in HIT {
            for style in &Style::ALL {
                let m = prep(smi);
                let d = crate::generate(&m, style);
                assert!(
                    !d.degraded.is_empty(),
                    "[{}] {smi}:桥环命中模板之后就不报退化了",
                    style.name
                );
            }
        }
    }
}
