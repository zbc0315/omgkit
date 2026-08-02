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
//! 由 `rings.rs` 里的 `regenerate_templates` 生成:对每个骨架跑两万次**带扰动**
//! 的多起点松弛,按**现成的** `Quality`(自交数、最大键长偏差、量化坐标序列)
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

/// 查表:这个环系有没有预存的坐标。
///
/// 返回的坐标按**父分子的原子编号**给出,已经对上号了。
pub(crate) fn lookup(
    mol: &MolBuilder,
    atoms: &[u32],
    ranks: &[u32],
) -> Option<BTreeMap<u32, Point2>> {
    let skel = skeleton_of(mol, atoms, ranks)?;
    let coords = TABLE.iter().find(|(k, _)| *k == skel).map(|(_, v)| *v)?;

    // 重建同一个骨架,拿它自己的规范秩去对坐标 —— 存的时候就是按这个存的。
    // 规范秩与存储序无关,所以两边算出来必然一致。
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
            b.add_bond(*x, *y, BondOrder::Single).ok()?;
        }
    }
    omgkit_chem::pipeline::sanitize(&mut b).ok()?;
    let skel_ranks = omgkit_io::canon::canonical_ranks(&b);
    if skel_ranks.len() != coords.len() {
        return None; // 表坏了 —— 不猜,退回松弛
    }

    let mut out = BTreeMap::new();
    for (i, parent) in back.iter().enumerate() {
        let r = *skel_ranks.get(i)? as usize;
        let (x, y) = *coords.get(r)?;
        out.insert(*parent, Point2::new(x, y));
    }
    Some(out)
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

    #[test]
    fn the_table_is_actually_used() {
        // 一条都命中不了的话,下面那些判据全是空过的。
        let mut hits = 0usize;
        for smi in HIT {
            let m = prep(smi);
            let ranks = omgkit_io::canon::canonical_ranks(&m);
            let rs = omgkit_chem::sssr::ring_set(&m);
            for sys in crate::rings::group(&omgkit_chem::rings::fused_ring_systems(&m), &rs) {
                if lookup(&m, &sys.atoms, &ranks).is_some() {
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
