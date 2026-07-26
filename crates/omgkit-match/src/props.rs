//! 把一个净化过的分子预计算成**逐原子逐键的查询性质**。
//!
//! # 为什么要预计算
//!
//! 匹配是回溯搜索:同一个原子会被反复拿来与不同的查询原子比对,一个 20 原子
//! 的模式配到 40 原子的分子上,单个目标原子被求值几十上百次很常见。
//!
//! 环成员数、最小环大小这些量要遍历环集才能算出来。放在求值里现算,等于把
//! "遍历环集"塞进了搜索的最内层循环 —— 无论搜索本身多快都救不回来。
//!
//! # 环的三个量来自两处
//!
//! `R`(所属环个数)与 `x`(环键数)要数**具体的环**,只能从环集来;
//! `r`(最小环大小)环感知那一步已经算好了,直接取。
//!
//! 这两处必须用**同一个环集**,否则会出现"`R0` 但 `r6`"这种自相矛盾的性质。
//! 所以本模块自己跑一遍环感知和环集,不接受调用方分别传进来。
//!
//! # 氢数要把**图里的氢原子**也算上
//!
//! SMARTS 的 `H<n>` 数的是"这个原子上一共有几个氢",不管氢是记成计数还是
//! 画成图里的一个节点。`[H]N(C)C` 的氮:显式氢 0、隐式氢 0,但 `[N;H1]`
//! 命中它 —— 那个氢是它的**邻居**。
//!
//! 漏掉这一项的失效方式很安静:绝大多数分子不写显式氢,测试全绿,直到遇上
//! 一条 `[H]N(...)` 才莫名其妙地不命中。

use omgkit_chem::{perceive_rings, ring_set};
use omgkit_core::{AtomFlags, BondOrder, MolBuilder};
use omgkit_io::smarts::{AtomProps, BondProps};

/// 一个分子的全部查询性质,按原子/键下标索引。
#[derive(Debug, Clone)]
pub struct MolProps {
    /// 逐原子
    pub atoms: Vec<AtomProps>,
    /// 逐键。配位键的朝向留给匹配时按方向决定,这里 `dative_forward` 恒为真。
    pub bonds: Vec<BondProps>,
}

impl MolProps {
    /// 从一个**已净化**的分子预计算。
    ///
    /// 分子必须已经跑过价键计算与芳香性感知 —— 隐式氢、芳香标志都直接取
    /// 分子上的字段。没净化的分子算出来的性质是错的,而且错得很安静:
    /// 隐式氢全是 0,芳香标志是"作者声称"而非感知结果。
    #[must_use]
    pub fn compute(mol: &MolBuilder) -> Self {
        // 环感知与环集必须来自同一次计算,见模块文档
        let mut scratch = mol.clone();
        let rings = perceive_rings(&mut scratch);
        let cycles = ring_set(&scratch);

        let n = mol.num_atoms();
        let mut ring_count = vec![0u32; n];
        let mut ring_bond_count = vec![0u32; n];
        for ring in &cycles {
            for &a in &ring.atoms {
                ring_count[a as usize] += 1;
            }
        }
        // 环键数要按**键**去重:同一条键属于多个环时,对端点只算一次
        let mut bond_is_ring = vec![false; mol.num_bonds()];
        for ring in &cycles {
            for &b in &ring.bonds {
                bond_is_ring[b as usize] = true;
            }
        }
        for (bi, &in_ring) in bond_is_ring.iter().enumerate() {
            if in_ring {
                let b = mol.bonds()[bi];
                ring_bond_count[b.begin as usize] += 1;
                ring_bond_count[b.end as usize] += 1;
            }
        }

        let atoms = (0..n)
            .map(|i| {
                let a = mol.atoms()[i];
                AtomProps {
                    atomic_num: a.atomic_num,
                    aromatic: a.flags.contains(AtomFlags::AROMATIC),
                    charge: i32::from(a.formal_charge),
                    isotope: a.isotope,
                    degree: mol.degree(i as u32) as u32,
                    total_hs: u32::from(a.num_explicit_hs)
                        + u32::from(a.num_implicit_hs)
                        + neighbour_hydrogens(mol, i as u32),
                    implicit_hs: u32::from(a.num_implicit_hs),
                    valence: valence_of(mol, i as u32),
                    ring_count: ring_count[i],
                    min_ring_size: u32::from(rings.atom_min_ring_size[i]),
                    ring_bonds: ring_bond_count[i],
                    chiral_tag: a.chiral_tag,
                    atom_map: a.atom_map,
                }
            })
            .collect();

        let bonds = (0..mol.num_bonds())
            .map(|i| {
                let b = mol.bonds()[i];
                BondProps {
                    order: b.order,
                    in_ring: rings.bond_in_ring[i],
                    direction: b.direction,
                    dative_forward: true,
                }
            })
            .collect();

        Self { atoms, bonds }
    }
}

/// 总价 = 各键的键级之和(四舍五入)+ 总氢数。
///
/// 芳香键按 1.5 计,所以要先求和再取整 —— 逐键取整会让苯环的碳变成 2+2+1=5。
fn valence_of(mol: &MolBuilder, atom: u32) -> u32 {
    let a = mol.atoms()[atom as usize];
    let sum: f32 = mol
        .neighbors(atom)
        .map(|(_, bi)| mol.bonds()[bi as usize].valence_contribution_to(atom))
        .sum();
    // 芳香环上每个原子的键级和是 x.5,进位与实际价一致
    let bonds = sum.round() as u32;
    bonds + u32::from(a.num_explicit_hs) + u32::from(a.num_implicit_hs)
}

/// 邻居里画成独立节点的氢原子数。
///
/// 只数**真的氢**:同位素(氘、氚)也算,通配原子不算。
fn neighbour_hydrogens(mol: &MolBuilder, atom: u32) -> u32 {
    mol.neighbors(atom)
        .filter(|&(other, _)| mol.atoms()[other as usize].atomic_num == 1)
        .count() as u32
}

/// 该原子上的键在 `BondOrder::Dative` 时,给体是不是 `from` 端。
#[must_use]
pub fn dative_points_from(mol: &MolBuilder, bond: u32, from: u32) -> bool {
    let b = mol.bonds()[bond as usize];
    b.order != BondOrder::Dative || b.begin == from
}
