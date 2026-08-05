//! 为了画出构型而补出来的显式氢。
//!
//! # 为什么非补不可
//!
//! 楔形的语义是"从窄端那个原子看出去,这根取代基在纸面上方/下方"。由此有两条
//! 约定:窄端必须在立体中心(已有判据守着);**楔形不该画在环键上** —— 环键的
//! 两个原子在读者眼里都躺在环平面里,声明其中一个出平面与读者正在用的环几何
//! 自相矛盾,稠环里更分不清楔形描述的是近端那个原子还是整个环的透视。IUPAC 2006
//! 的立体化学图示建议(Brecher, *Pure Appl. Chem.* 78:1897)就是这么说的:
//! 立体键应当画向**取代基**。
//!
//! 于是数一下这个中心还剩什么键可用:
//!
//! | 中心 | 重原子键 | 结论 |
//! |---|---|---|
//! | 有环外单键 | 至少一根 | 拿它打楔形,**不用画 H** |
//! | 三根键全在环上、**有 H** | 只剩环键 | 唯一合法的楔形是 C–H,**H 必须画出来** |
//! | 四根键全在环上、无 H | 只剩环键 | **无解** —— 补不了也没有合法楔形,如实报出来 |
//!
//! 甾体的稠合碳 C8/C9/C14 正是第二档。全量语料 **149 个中心 / 85 个分子**。
//!
//! # 判据用直接原因,不用拓扑代理
//!
//! RDKit 的口径是"立体中心 **且** 在 ≥2 个环里"
//! (`MolDraw2DUtils::isAtomCandForChiralH`)。那是**代理**:在两个环里通常
//! 意味着三根键都是环键。这里用**直接原因** —— "有没有一根能打楔形的非环单键"。
//!
//! 两者在全量语料上**行为等价**:直接判据从不多补(0 个);RDKit 的代理多圈进
//! 40 个中心,但它们全是 `deg=4, hs=0` 的季碳,`MolOps::addHs` 对没有氢的原子
//! 是空操作。差集全在"补不了"那一档 —— 直接判据把无解情形直接摆到明面上。

use omgkit_core::{AtomData, BondData, BondOrder, ChiralTag, MolBuilder};

/// 补出来的东西。空的话画的就是传入的分子本身。
///
/// **只记增量**,不复制整个分子:原子/键的下标因此与传入的分子逐项对应,
/// 追加的部分排在后面。
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Augmented {
    /// 被改过氢计数的中心:`(原子号, 改后的数据)`。
    ///
    /// 补一个显式氢就要把中心的氢计数减一,否则标签写成 `CH` 而氢又画了一遍,
    /// 价键也对不上。
    pub edited: Vec<(u32, AtomData)>,
    /// 追加的原子 —— 补出来的氢,接在原分子的原子之后。
    pub atoms: Vec<AtomData>,
    /// 追加的键 —— 每个氢一根,接在原分子的键之后。
    pub bonds: Vec<BondData>,
}

impl Augmented {
    /// 把增量贴回去,得到**真正被画的那个分子**。
    ///
    /// 前 `mol.num_atoms()` 个原子、前 `mol.num_bonds()` 根键与传入的分子逐项
    /// 对应,所以按原下标索引仍然是对的。
    #[must_use]
    pub fn apply(&self, mol: &MolBuilder) -> MolBuilder {
        let mut out = mol.clone();
        for (a, data) in &self.edited {
            if let Some(x) = out.atom_mut(*a) {
                *x = *data;
            }
        }
        for a in &self.atoms {
            out.add_atom_data(*a);
        }
        for b in &self.bonds {
            // 增量是自己造出来的,端点必然在范围内
            out.add_bond_data(*b).expect("补出来的键端点越界");
        }
        out
    }
}

/// 这个中心要不要补一个显式氢。
///
/// 见模块文档:没有**能打楔形的非环单键**、而它还有氢,就要补。
fn needs_h(mol: &MolBuilder, a: u32, rings: &[omgkit_chem::sssr::Ring]) -> bool {
    let at = mol.atoms()[a as usize];
    if !matches!(at.chiral_tag, ChiralTag::Cw | ChiralTag::Ccw) {
        return false;
    }
    if u32::from(at.num_explicit_hs) + u32::from(at.num_implicit_hs) == 0 {
        return false; // 没氢可补 —— 那一档无解,不在这里处理
    }
    !mol.neighbors(a).any(|(_, b)| {
        mol.bonds()[b as usize].order == BondOrder::Single
            && !rings.iter().any(|r| r.bonds.contains(&b))
    })
}

/// 给需要的立体中心补上显式氢;一个都不需要就返回 `None`。
///
/// # 原有的原子/键编号一概不变
///
/// [`MolBuilder::add_atom_data`] 是追加,半边链表是尾插 —— 所以补出来的氢排在
/// 原子表末尾、C–H 键排在键表末尾,而且是中心的**最后一个**邻居。原下标全部
/// 对得上,坐标可以直接用回原分子。
///
/// 这与 `as_plain_bonds`(为布局做一份把配位键当单键的副本)是同一个套路。
///
/// # 补的顺序按规范秩
///
/// 补哪些中心是拓扑决定的,与写法无关;但**追加的次序**若跟着存储下标走,
/// 同一个分子换种写法补出来的氢就会拿到不同的原子号,
/// `canonical_ranks` 跟着变,整张图就变了。所以按中心的规范秩排。
#[must_use]
pub fn with_stereo_hs(mol: &MolBuilder) -> Option<Augmented> {
    let rings = omgkit_chem::sssr::ring_set(mol);
    let genuine = omgkit_io::stereo::genuine_tetrahedral(mol);
    let ranks = omgkit_io::canon::canonical_ranks(mol);

    let mut centres: Vec<u32> = (0..u32::try_from(mol.num_atoms()).expect("原子数超出 u32"))
        .filter(|a| genuine[*a as usize] && needs_h(mol, *a, &rings))
        .collect();
    if centres.is_empty() {
        return None;
    }
    centres.sort_by_key(|a| (ranks[*a as usize], *a));

    let mut out = Augmented::default();
    let first = u32::try_from(mol.num_atoms()).expect("原子数超出 u32");
    for (next, a) in (first..).zip(centres) {
        let mut at = mol.atoms()[a as usize];
        // **减非零的那个字段。** 两个字段是互斥的(全仓约定,相加即总氢数),
        // 所以"先减哪个"其实无所谓 —— 真正的坑是**无条件只减某一个**:
        // `[C@H]` 把氢记在 `num_explicit_hs` 上,只减 `num_implicit_hs` 的话
        // 那次 `saturating_sub` 悄悄什么也没做,价键当场超限,`sanitize` 报
        // 「原子 #1(C):显式价超出允许范围,显式价 = 5」。
        if at.num_explicit_hs > 0 {
            at.num_explicit_hs -= 1;
        } else {
            at.num_implicit_hs -= 1;
        }
        out.edited.push((a, at));
        out.atoms.push(AtomData {
            atomic_num: 1,
            ..AtomData::default()
        });
        out.bonds.push(BondData::new(a, next, BondOrder::Single));
    }
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// splitmix64 + Fisher–Yates。与 `examples/audit.rs` 里那个同源 ——
    /// 仿射式的"置换"搅不动东西,见下面判据里的注释。
    fn shuffled(n: usize, seed: u64) -> Vec<u32> {
        let mut s = seed.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut next = || {
            s = s.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = s;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        };
        let mut v: Vec<u32> = (0..u32::try_from(n).unwrap()).collect();
        for i in (1..n).rev() {
            let j = (next() % (i as u64 + 1)) as usize;
            v.swap(i, j);
        }
        v
    }

    fn prep(smi: &str) -> MolBuilder {
        let mut m = omgkit_io::smiles::parse(smi).unwrap();
        omgkit_chem::pipeline::sanitize(&mut m).unwrap();
        m
    }

    /// 甾体:C8/C9/C14 三根键全在环上,只能靠 C–H 打楔形
    const STEROID: &str =
        "CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C";

    #[test]
    fn a_molecule_that_needs_nothing_gets_nothing() {
        // **爆炸半径要钉死。** 不需要补氢的分子必须原样返回 `None` —— 一旦
        // 无条件返回 `Some`,`canonical_ranks` 跟着变,整个语料的布局全都要动。
        for smi in ["CC(=O)Oc1ccccc1C(=O)O", "C[C@H](N)C(=O)O", "c1ccccc1"] {
            assert!(
                with_stereo_hs(&prep(smi)).is_none(),
                "{smi} 不需要补氢,却补了"
            );
        }
    }

    #[test]
    fn a_fused_ring_stereocentre_gets_one() {
        let m = prep(STEROID);
        let aug = with_stereo_hs(&m).expect("甾体该补氢");
        assert_eq!(aug.atoms.len(), 3, "甾体该补三个:C8/C9/C14");
        assert_eq!(aug.bonds.len(), aug.atoms.len());
        for (a, _) in &aug.edited {
            // 补的都该是"三根键全在环上"的中心
            let rings = omgkit_chem::sssr::ring_set(&m);
            assert!(
                m.neighbors(*a)
                    .all(|(_, b)| rings.iter().any(|r| r.bonds.contains(&b))),
                "补氢的中心 {a} 有环外键可用,本不该补"
            );
        }
    }

    #[test]
    fn the_original_numbering_survives() {
        // 坐标是按原下标用回去的,前缀一变就全错位。
        let m = prep(STEROID);
        let aug = with_stereo_hs(&m).expect("甾体该补氢");
        let m2 = aug.apply(&m);
        assert_eq!(m2.num_atoms(), m.num_atoms() + aug.atoms.len());
        assert_eq!(m2.num_bonds(), m.num_bonds() + aug.bonds.len());
        let edited: std::collections::BTreeMap<u32, AtomData> =
            aug.edited.iter().copied().collect();
        for i in 0..m.num_atoms() {
            let a = u32::try_from(i).unwrap();
            let want = edited.get(&a).copied().unwrap_or(m.atoms()[i]);
            assert_eq!(m2.atoms()[i], want, "原子 {i} 的数据变了");
        }
        for i in 0..m.num_bonds() {
            assert_eq!(m2.bonds()[i], m.bonds()[i], "键 {i} 的数据变了");
        }
        // 补出来的氢是中心的**最后一个**邻居 —— 手性标记的槽位论证靠这一条
        for (k, b) in aug.bonds.iter().enumerate() {
            let last = m2.neighbors(b.begin).last().expect("中心有邻居");
            assert_eq!(
                last.0, b.end,
                "第 {k} 个补出来的氢不是中心 {} 的最后一个邻居",
                b.begin
            );
        }
    }

    #[test]
    fn the_augmented_molecule_still_sanitises() {
        // 氢计数减错字段的话价键当场超限 —— 实测报「显式价 = 5」。
        let m = prep(STEROID);
        let mut m2 = with_stereo_hs(&m).expect("甾体该补氢").apply(&m);
        omgkit_chem::pipeline::sanitize(&mut m2).expect("补完氢还得是个合法分子");
        // 总氢数守恒:中心少一个隐式,分子多一个显式原子
        let total = |x: &MolBuilder| -> usize {
            x.atoms()
                .iter()
                .map(|a| {
                    usize::from(a.num_explicit_hs)
                        + usize::from(a.num_implicit_hs)
                        + usize::from(a.atomic_num == 1)
                })
                .sum()
        };
        assert_eq!(total(&m2), total(&m), "补氢前后总氢数不一样");
    }

    #[test]
    fn which_centres_get_hs_does_not_depend_on_how_it_was_written() {
        // 补哪些中心是拓扑决定的,但**追加的次序**若跟着存储下标走,同一分子
        // 换种写法补出来的氢就会拿到不同的原子号 —— `canonical_ranks` 跟着变,
        // 整张图就变了。这是本 crate 的头号契约。
        let m = prep(STEROID);
        let ranks = omgkit_io::canon::canonical_ranks(&m);
        let aug = with_stereo_hs(&m).expect("甾体该补氢");
        // 补氢中心的**规范秩序列**是与写法无关的指纹
        let want: Vec<u32> = aug.edited.iter().map(|(a, _)| ranks[*a as usize]).collect();

        let n = m.num_atoms();
        let mut compared = 0usize;
        for seed in 0..24u64 {
            // **搅拌器要货真价实。** 头一版用 `(i*k+k) % n` 凑优先序 —— 那是
            // 仿射映射,搅出来的置换太规整:把补氢的次序改成按存储下标(本该
            // 打红这条判据的变异)它一次都没抓住。审计里早就记过同一个坑
            // (乘法哈希有 10.85% 的改写原样返回),我又踩了一遍。
            // splitmix64 + Fisher–Yates。
            let priority = shuffled(n, seed);
            let w = omgkit_io::smiles::write_with_priority(&m, &priority);
            let Ok(mut m2) = omgkit_io::smiles::parse(&w.smiles) else {
                continue;
            };
            if omgkit_chem::pipeline::sanitize(&mut m2).is_err() {
                continue;
            }
            if omgkit_io::canon::canonical_smiles(&m2).smiles
                != omgkit_io::canon::canonical_smiles(&m).smiles
            {
                continue; // 改写没保住分子,不算数
            }
            let r2 = omgkit_io::canon::canonical_ranks(&m2);
            let aug2 = with_stereo_hs(&m2).expect("同一个分子,照样该补氢");
            let got: Vec<u32> = aug2.edited.iter().map(|(a, _)| r2[*a as usize]).collect();
            assert_eq!(got, want, "换成 {} 之后补氢的中心变了", w.smiles);
            compared += 1;
        }
        assert!(compared > 0, "一次都没比成 —— 这条判据是空过的");
    }
}
