//! 环感知:判定原子/键是否在环中,以及过每个原子的最短环有多大(净化第 4 步)。
//!
//! 这三个量都是**纯图论量**,与"具体选哪一组环"无关:
//!
//! | 量 | 定义 |
//! |---|---|
//! | 键在环中 | 该边不是桥 |
//! | 原子在环中 | 关联至少一条非桥边 |
//! | 原子的最小环大小 | 过该原子的最短环长度 |
//!
//! 需要**具体是哪些环**的场合(芳香性感知按环判定)由 [`crate::sssr`] 提供。
//!
//! # 配位键不参与成环
//!
//! 配位键是形式上的电子对给予,不构成骨架连接。把它算进环里会在有机金属
//! 分子上给出多余的环 —— 这是可观察的差异,不是学究。

use omgkit_core::{AtomFlags, BondFlags, BondOrder, MolBuilder};

/// 一次环感知的结果。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RingPerception {
    /// 每原子的最小环大小。
    ///
    /// **0 表示"这里没有数",不等于"不在环中"** —— 环大过 [`MAX_RING_SIZE`]
    /// 时同样是 0。要问在不在环里,读 [`atom_in_ring`](Self::atom_in_ring)。
    /// 两者混为一谈过一次:256 元以上的环原子会同时报"在环里"和"不在任何环中"。
    pub atom_min_ring_size: Vec<u8>,
    /// 每原子是否在环中
    pub atom_in_ring: Vec<bool>,
    /// 每键是否在环中(等价于"不是桥")
    pub bond_in_ring: Vec<bool>,
}

/// 环大小的上限:**存这个数的字段有多宽**([`RingPerception::atom_min_ring_size`]
/// 是 `u8`)。比这还大的环记为 0。
///
/// 先前是 20,而那个数字挡住的是真实存在的东西:30 元大环在 `hard.smi` 里就有,
/// 三十元的环内酯与环肽在天然产物里也不罕见。抬到 255 之后压力语料(环原子挂
/// 200 长尾、30 元环每个环原子挂 20 长尾)实测与 20 无法区分 —— BFS 的深度上限
/// 在**第一轮**之后就被 `best` 剪枝接管了,上限只影响第一轮。
///
/// **仍然碰不到的那一档:256 元以上的环。** 它们的 `atom_min_ring_size` 是 0,
/// 而"在不在环里"由 [`RingPerception::atom_in_ring`] 独立回答,不受这个上限影响
/// —— SMARTS 的 `[r0]` / `[R0]` 读的是后者,所以不会把大环原子说成非环原子。
pub const MAX_RING_SIZE: usize = 255;

/// 先按这个上限跑一遍。真实分子的环几乎全在这里面,一遍就出结果。
///
/// **这个数只影响代价,不影响结果** —— 找不到 ≤ 它的环时会用
/// [`MAX_RING_SIZE`] 再跑一遍,答案与只跑一遍大上限相同。判据
/// `a_ring_reports_its_own_size_all_the_way_up_to_the_field_width` 一路扫到
/// 60 元环,正好跨过这条线,所以"两遍给的答案与一遍相同"这件事是被守着的。
const COMMON_RING_MAX: usize = 20;

/// 对分子做环感知,并就地设置 [`AtomFlags::IN_RING`] / [`BondFlags::IN_RING`]。
///
/// 配位键不参与成环判定。
#[must_use]
pub fn perceive_rings(mol: &mut MolBuilder) -> RingPerception {
    let n_atoms = mol.num_atoms();
    let n_bonds = mol.num_bonds();

    // 参与成环的边:排除配位键
    let active: Vec<bool> = mol
        .bonds()
        .iter()
        .map(|b| b.order != BondOrder::Dative)
        .collect();

    let adj = Adjacency::build(mol, &active);
    let is_bridge = find_bridges(&adj);

    let mut bond_in_ring = vec![false; n_bonds];
    for (bi, in_ring) in bond_in_ring.iter_mut().enumerate() {
        *in_ring = active[bi] && !is_bridge[bi];
    }

    let mut atom_in_ring = vec![false; n_atoms];
    for (bi, &in_ring) in bond_in_ring.iter().enumerate() {
        if in_ring {
            let b = mol.bonds()[bi];
            atom_in_ring[b.begin as usize] = true;
            atom_in_ring[b.end as usize] = true;
        }
    }

    let mut atom_min_ring_size = vec![0u8; n_atoms];
    // BFS 的暂存空间只分配一次,所有原子共用 —— 见 `shortest_cycle_through`
    let mut dist = vec![u32::MAX; n_atoms];
    let mut queue: Vec<u32> = Vec::new();
    for a in 0..n_atoms {
        if !atom_in_ring[a] {
            continue;
        }
        // **两遍。** 先用 [`COMMON_RING_MAX`] 跑一遍 —— 真实分子的环几乎全在这个
        // 范围里,一遍就出结果,代价与只有小上限时一模一样。
        //
        // 只有那一遍找不到环的原子才付大上限那笔钱。这一层不是优化,是**必须的**:
        // BFS 的深度上限直接是它的代价上限,而稠合体系里每根键都是环键,剪枝
        // 之前的首轮会一直走到深度上限。上限一律用 255 时,`tests/scaling.rs` 的
        // "单个大稠合体系(线性并苯)"那一档从 2.47 µs/原子涨到 3.16 —— 平方项。
        let small = shortest_cycle_through(
            &adj,
            a as u32,
            &bond_in_ring,
            COMMON_RING_MAX,
            &mut dist,
            &mut queue,
        );
        atom_min_ring_size[a] = if small != 0 {
            small
        } else {
            shortest_cycle_through(
                &adj,
                a as u32,
                &bond_in_ring,
                MAX_RING_SIZE,
                &mut dist,
                &mut queue,
            )
        };
    }
    debug_assert!(
        dist.iter().all(|&d| d == u32::MAX),
        "shortest_cycle_through 没把 dist 复位干净"
    );

    for (a, &in_ring) in atom_in_ring.iter().enumerate() {
        if let Some(atom) = mol.atom_mut(a as u32) {
            atom.flags.set(AtomFlags::IN_RING, in_ring);
        }
    }
    for (bi, &in_ring) in bond_in_ring.iter().enumerate() {
        if let Some(mut bond) = mol.bond_mut(bi as u32) {
            bond.flags_mut().set(BondFlags::IN_RING, in_ring);
        }
    }

    RingPerception {
        atom_min_ring_size,
        atom_in_ring,
        bond_in_ring,
    }
}

// ---------------------------------------------------------------------------

/// CSR 邻接,只含参与成环的边。
pub(crate) struct Adjacency {
    offset: Vec<u32>,
    /// (邻居原子, 键下标)
    nbr: Vec<(u32, u32)>,
    n_atoms: usize,
    /// 分子的**总**键数(含被滤掉的),用来给按键下标寻址的数组定长
    n_bonds: usize,
}

impl Adjacency {
    /// 从 `MolBuilder` 的邻接索引压成 CSR,只保留 `active` 的边。
    ///
    /// 之所以单独建一份而不是直接用 `mol.neighbors()`:配位键要滤掉,
    /// 而且下面的迭代式 Tarjan 需要按下标随机访问邻居切片。
    pub(crate) fn build(mol: &MolBuilder, active: &[bool]) -> Self {
        let n_atoms = mol.num_atoms();
        let mut offset = vec![0u32; n_atoms + 1];
        for a in 0..n_atoms {
            let deg = mol
                .neighbors(a as u32)
                .filter(|&(_, bi)| active[bi as usize])
                .count();
            offset[a + 1] = deg as u32;
        }
        for i in 1..=n_atoms {
            offset[i] += offset[i - 1];
        }

        let mut nbr = Vec::with_capacity(offset[n_atoms] as usize);
        for a in 0..n_atoms {
            nbr.extend(
                mol.neighbors(a as u32)
                    .filter(|&(_, bi)| active[bi as usize]),
            );
        }
        Self {
            offset,
            nbr,
            n_atoms,
            n_bonds: mol.num_bonds(),
        }
    }

    fn neighbors(&self, a: u32) -> &[(u32, u32)] {
        let s = self.offset[a as usize] as usize;
        let e = self.offset[a as usize + 1] as usize;
        &self.nbr[s..e]
    }
}

/// 桥边判定(Tarjan low-link)。
///
/// **迭代实现,不用递归** —— 递归在长链分子(肽、聚合物)上会爆栈,
/// 而这类分子在真实语料里并不罕见。
fn find_bridges(adj: &Adjacency) -> Vec<bool> {
    let n = adj.n_atoms;
    // 按分子的**总**键数定长:被滤掉的键(配位键)也占着下标。
    // 数组短了,越界与否就取决于调用方的短路求值,不该这么依赖
    let mut is_bridge = vec![false; adj.n_bonds];

    let mut disc = vec![u32::MAX; n]; // 发现时间;MAX = 未访问
    let mut low = vec![u32::MAX; n];
    let mut timer = 0u32;

    // 显式栈:(当前原子, 来自哪条边, 已处理到第几个邻居)
    let mut stack: Vec<(u32, u32, usize)> = Vec::new();

    for root in 0..n as u32 {
        if disc[root as usize] != u32::MAX {
            continue;
        }
        disc[root as usize] = timer;
        low[root as usize] = timer;
        timer += 1;
        stack.push((root, u32::MAX, 0));

        while let Some(&mut (v, from_bond, ref mut k)) = stack.last_mut() {
            let nbrs = adj.neighbors(v);
            if *k < nbrs.len() {
                let (u, bond) = nbrs[*k];
                *k += 1;
                if bond == from_bond {
                    continue; // 不沿来路回退(同一条边)
                }
                if disc[u as usize] == u32::MAX {
                    disc[u as usize] = timer;
                    low[u as usize] = timer;
                    timer += 1;
                    stack.push((u, bond, 0));
                } else {
                    low[v as usize] = low[v as usize].min(disc[u as usize]);
                }
            } else {
                stack.pop();
                if let Some(&mut (parent, _, _)) = stack.last_mut() {
                    low[parent as usize] = low[parent as usize].min(low[v as usize]);
                    if low[v as usize] > disc[parent as usize] {
                        is_bridge[from_bond as usize] = true;
                    }
                }
            }
        }
    }
    is_bridge
}

/// 过原子 `v` 的最短环长度;不在环中或超过 [`MAX_RING_SIZE`] 时返回 0。
///
/// **返回 0 有两种意思**,调用方要问"在不在环里"得看 `atom_in_ring`。
///
/// 做法:把 `v` 摘掉,从它的每个邻居做一次 BFS,取任意两个邻居之间的最短距离,
/// 环长 = 该距离 + 2。摘掉 `v` 保证得到的一定是**简单环** —— 这比"从 v 做一次
/// BFS 再撞非树边"的经典写法更容易证明正确。
///
/// # 三处必须做对的开销控制
///
/// 本函数对**每个环原子**都要调一次,所以任何"正比于分子大小"的动作都会
/// 变成 O(原子数²):
///
/// 1. `dist` / `queue` 由调用方提供并复用,不在函数内分配
/// 2. 复位只动**本次访问过**的下标(它们全在 `queue` 里),不整体刷一遍
/// 3. **只走环键。** 过 `v` 的环只可能由环键组成,桥一根也用不上 —— 于是 BFS
///    走不出 `v` 所在的那个环系,深度天然被环系大小卡住。先前靠的是
///    `MAX_RING_SIZE - 2` 这个深度上限,而那让上限的值直接变成性能参数:
///    把它从 20 抬到 255(为了让 30 元大环报得出大小),首轮 BFS 就会溜进无环的
///    尾巴,`tests/scaling.rs` 的"很多个独立小环系"那一档当场从 1.98 µs/原子
///    涨到 2.79。深度上限仍然留着,但它现在只管"装不装得进 `u8`",不再管代价。
///
/// 三条都做之后,代价正比于 `v` 附近的局部邻域,与分子总大小无关。
///
/// 约定:`dist` 进入时必须全为 `u32::MAX`,返回时恢复原样。
/// 这条约定由调用方在**全部原子处理完之后**查一次(见 [`perceive_rings`])——
/// 逐次检查是 O(原子数²);漏清理的项会一直留到最后,末尾查同样抓得到。
fn shortest_cycle_through(
    adj: &Adjacency,
    v: u32,
    bond_in_ring: &[bool],
    max_size: usize,
    dist: &mut [u32],
    queue: &mut Vec<u32>,
) -> u8 {
    let nbrs: Vec<u32> = adj
        .neighbors(v)
        .iter()
        .filter(|&&(_, bi)| bond_in_ring[bi as usize])
        .map(|&(u, _)| u)
        .collect();
    if nbrs.len() < 2 {
        return 0;
    }

    // 环长 = 两邻居间距 + 2,要 ≤ MAX_RING_SIZE 就得间距 ≤ MAX_RING_SIZE - 2
    let max_dist = (max_size - 2) as u32;
    let mut best = usize::MAX;

    for (i, &start) in nbrs.iter().enumerate() {
        // 只需对 i < j 的配对求距离,故从第 i 个邻居出发即可覆盖全部配对
        queue.clear();
        dist[start as usize] = 0;
        queue.push(start);
        let mut head = 0;
        while head < queue.len() {
            let x = queue[head];
            head += 1;
            let dx = dist[x as usize];
            // 深度上限:再远也凑不出足够短的环
            if dx >= max_dist {
                continue;
            }
            // 剪枝:再走下去也超不过当前最优
            if best != usize::MAX && dx as usize + 2 >= best {
                continue;
            }
            for &(y, bi) in adj.neighbors(x) {
                if y == v {
                    continue; // v 已摘除
                }
                if !bond_in_ring[bi as usize] {
                    continue; // 桥凑不出环,走进去只是白走
                }
                if dist[y as usize] == u32::MAX {
                    dist[y as usize] = dx + 1;
                    queue.push(y);
                }
            }
        }
        for &other in &nbrs[i + 1..] {
            let d = dist[other as usize];
            if d != u32::MAX {
                best = best.min(d as usize + 2);
            }
        }
        // 复位:queue 恰好装着本轮访问过的全部节点
        for &x in queue.iter() {
            dist[x as usize] = u32::MAX;
        }
    }

    if best == usize::MAX || best > max_size {
        0
    } else {
        u8::try_from(best).unwrap_or(0)
    }
}

#[cfg(test)]
mod tests {
    use omgkit_io::smiles;

    use super::*;

    fn perceive(smi: &str) -> (MolBuilder, RingPerception) {
        let mut m = smiles::parse(smi).unwrap_or_else(|e| panic!("{}", e.render()));
        let r = perceive_rings(&mut m);
        (m, r)
    }

    #[test]
    fn acyclic_has_no_rings() {
        let (_, r) = perceive("CCO");
        assert!(r.atom_in_ring.iter().all(|&x| !x));
        assert!(r.bond_in_ring.iter().all(|&x| !x));
        assert_eq!(r.atom_min_ring_size, vec![0, 0, 0]);
    }

    #[test]
    fn simple_rings() {
        for (smi, size) in [
            ("C1CC1", 3),
            ("C1CCC1", 4),
            ("C1CCCCC1", 6),
            ("c1ccccc1", 6),
        ] {
            let (_, r) = perceive(smi);
            assert!(r.atom_in_ring.iter().all(|&x| x), "{smi}:全部原子应在环中");
            assert!(r.bond_in_ring.iter().all(|&x| x), "{smi}:全部键应在环中");
            assert!(
                r.atom_min_ring_size.iter().all(|&s| s == size),
                "{smi}:最小环应为 {size},实际 {:?}",
                r.atom_min_ring_size
            );
        }
    }

    #[test]
    fn substituent_is_not_in_ring() {
        // 甲苯:环上 6 个碳在环中,甲基不在
        let (_, r) = perceive("Cc1ccccc1");
        assert!(!r.atom_in_ring[0], "甲基碳不应在环中");
        assert!(r.atom_in_ring[1..].iter().all(|&x| x));
        assert!(!r.bond_in_ring[0], "甲基-环 的键是桥");
        assert_eq!(r.atom_min_ring_size[0], 0);
        assert!(r.atom_min_ring_size[1..].iter().all(|&s| s == 6));
    }

    #[test]
    fn fused_rings_take_the_smallest() {
        // 萘:融合位的两个碳同属两个六元环
        let (_, r) = perceive("c1ccc2ccccc2c1");
        assert!(r.atom_min_ring_size.iter().all(|&s| s == 6));
        assert!(r.atom_in_ring.iter().all(|&x| x));
    }

    #[test]
    fn spiro_and_bridged() {
        // 螺环:共用一个原子的 3+4 环
        let (_, r) = perceive("C1CC12CCC2");
        assert_eq!(r.atom_min_ring_size[2], 3, "螺原子应取较小的环");
        // 桥环双环[2.2.2]辛烷
        let (_, r) = perceive("C1CC2CCC1CC2");
        assert!(r.atom_in_ring.iter().all(|&x| x));
    }

    #[test]
    fn disconnected_fragments() {
        // 一个有环一个无环
        let (_, r) = perceive("C1CC1.CCO");
        assert_eq!(&r.atom_in_ring, &[true, true, true, false, false, false]);
        assert_eq!(r.atom_min_ring_size, vec![3, 3, 3, 0, 0, 0]);
    }

    #[test]
    fn long_chain_does_not_overflow_stack() {
        // 桥判定必须是迭代的:递归在长链上会爆栈
        let smi = "C".repeat(20_000);
        let mut m = smiles::parse(&smi).unwrap();
        let r = perceive_rings(&mut m);
        assert_eq!(m.num_atoms(), 20_000);
        assert!(r.bond_in_ring.iter().all(|&x| !x));
    }

    #[test]
    fn dative_bonds_do_not_form_rings() {
        let mut m = smiles::parse("C1CC1").unwrap();
        // 把其中一条环键改成配位键,环就应当"断开"
        m.bond_mut(0).unwrap().set_order(BondOrder::Dative);
        let r = perceive_rings(&mut m);
        assert!(
            r.bond_in_ring.iter().all(|&x| !x),
            "配位键不参与成环,剩下的边全变成桥"
        );
        assert!(r.atom_in_ring.iter().all(|&x| !x));
    }

    /// **环大小逐个扫过去,不只扫小环。**
    ///
    /// 先前这条不变量只在一张最大到六元环的表上跑,而实现里有一个 20 的上限:
    /// 21 元以上的环 `atom_min_ring_size` 是 0,于是同一个原子同时报"在环里"
    /// (桥判定说的)和"不在任何环中"(0 被这么读)。判据的参照侧那时写的是
    /// `range(3, 21)` —— **与被测常量是同一个数**,两边一起给 0,一次都没红过。
    ///
    /// 这里断的是契约:**环有多大,`atom_min_ring_size` 就报多大**,一直到存这个
    /// 数的字段装不下为止。上限碰不到的那一档(256 元以上)由下一条判据管。
    #[test]
    fn a_ring_reports_its_own_size_all_the_way_up_to_the_field_width() {
        for n in 3..=60usize {
            let smi = format!("C1{}1", "C".repeat(n - 1));
            let (_, r) = perceive(&smi);
            assert!(
                r.atom_in_ring.iter().all(|&x| x),
                "{n} 元环:有原子没被判成环原子"
            );
            for (i, &size) in r.atom_min_ring_size.iter().enumerate() {
                assert_eq!(
                    usize::from(size),
                    n,
                    "{n} 元环的第 {i} 个原子报的最小环大小是 {size}"
                );
            }
        }
    }

    /// **上限之外那一档也要写下来。** 全绿的校准表会被读成"守住了"。
    ///
    /// 256 元以上的环装不进 `u8`,`atom_min_ring_size` 只能是 0。那一档仍然必须
    /// 由 `atom_in_ring` 独立答对 —— 两个字段来自两个算法,大环只让其中一个失效。
    #[test]
    fn a_ring_too_big_for_the_field_still_knows_it_is_a_ring() {
        let n = MAX_RING_SIZE + 40;
        let smi = format!("C1{}1", "C".repeat(n - 1));
        let (_, r) = perceive(&smi);
        assert!(
            r.atom_in_ring.iter().all(|&x| x),
            "{n} 元环的原子没被判成环原子 —— 上限把两个字段一起打翻了"
        );
        assert!(
            r.atom_min_ring_size.iter().all(|&s| s == 0),
            "{n} 元环居然报出了大小,那 MAX_RING_SIZE 的注释就过期了"
        );
    }

    /// 内部一致性:原子在环中 ⟺ 有有限的最小环大小。
    /// 两者由**不同算法**得出(桥判定 vs 最短环 BFS),互为独立校验。
    #[test]
    fn membership_agrees_with_min_ring_size() {
        for smi in [
            "CCO",
            "C1CCCCC1",
            "Cc1ccccc1",
            "c1ccc2ccccc2c1",
            "C1CC2CCC1CC2",
            "C1CC12CCC2",
            "CC(=O)Oc1ccccc1C(=O)O",
            "CN1C=NC2=C1C(=O)N(C)C(=O)N2C",
            "C1CC1.CCO",
        ] {
            let (_, r) = perceive(smi);
            for (i, (&in_ring, &size)) in r
                .atom_in_ring
                .iter()
                .zip(r.atom_min_ring_size.iter())
                .enumerate()
            {
                assert_eq!(
                    in_ring,
                    size > 0,
                    "{smi} 原子 {i}:在环中={in_ring} 但最小环={size}"
                );
            }
        }
    }

    #[test]
    fn flags_are_written_back() {
        let (m, r) = perceive("Cc1ccccc1");
        for (i, a) in m.atoms().iter().enumerate() {
            assert_eq!(a.flags.contains(AtomFlags::IN_RING), r.atom_in_ring[i]);
        }
        for (i, b) in m.bonds().iter().enumerate() {
            assert_eq!(b.flags.contains(BondFlags::IN_RING), r.bond_in_ring[i]);
        }
    }
}

#[cfg(test)]
mod organometallic_root_cause {
    //! 环感知的结果依赖净化第 2 步(有机金属键改配位键)是否已经执行。
    //!
    //! 二茂铁一类的分子里,金属与配体之间的键会被第 2 步改成配位键;
    //! 环感知排除配位键,于是那条键不再成环。
    //!
    //! 本测试把这个依赖关系钉死:同一个分子,跑不跑第 2 步,环感知给出不同的
    //! 结果。它守的是**顺序**而不是环感知本身 —— 谁把第 2 步从管线里挪走或
    //! 挪到环感知之后,这里立刻红。

    use omgkit_io::smiles;

    use super::*;

    const FERROCENE: &str = "CN(C)C[C-]12C3=C4C5=C1[Fe++]23456789[C-]%10C6=C7C8=C9%10";

    #[test]
    fn ring_result_follows_the_bond_becoming_dative() {
        // 不跑第 2 步:键 15 还是单键,判为在环中
        let mut before = smiles::parse(FERROCENE).unwrap();
        crate::clean_up(&mut before);
        let r_before = perceive_rings(&mut before);
        assert_eq!(before.bonds()[15].order, BondOrder::Single);
        assert!(r_before.bond_in_ring[15]);
        assert_eq!(r_before.atom_min_ring_size[4], 3);

        // 跑第 2 步:调**真的**那一步,不手工模拟 —— 手工模拟只能证明
        // "改成配位键会怎样",证明不了"第 2 步确实会把这条键改掉"
        let mut after = smiles::parse(FERROCENE).unwrap();
        crate::clean_up(&mut after);
        assert_eq!(
            crate::cleanup_organometallics(&mut after),
            1,
            "第 2 步应当恰好改动一条键"
        );
        assert_eq!(after.bonds()[15].order, BondOrder::Dative);
        let r_after = perceive_rings(&mut after);

        assert!(!r_after.bond_in_ring[15], "配位键应被排除在环外");
        assert_eq!(r_after.atom_min_ring_size[4], 4, "三元环断开后应退到四元环");
    }
}

/// 融合环系 —— kekulize 的作用单元。
///
/// "融合"的定义是**共享至少一条键**;等价的图论表述是**含环的双连通分量**:
///
/// - 共享一条键的两个环必然落在同一双连通分量里
/// - 螺环只共享一个原子,不算融合;螺原子恰是割点,双连通分解天然把它们分开
/// - 联苯那样由一条桥键连接的两个环也分属不同分量,桥自成一个分量
///
/// 用双连通分解就不必先求出具体的环集。
///
/// 返回的每个分量内原子按下标升序,分量之间按最小原子下标升序 ——
/// 顺序确定,便于复现。
#[must_use]
pub fn fused_ring_systems(mol: &MolBuilder) -> Vec<Vec<u32>> {
    let active: Vec<bool> = mol
        .bonds()
        .iter()
        .map(|b| b.order != BondOrder::Dative)
        .collect();
    let adj = Adjacency::build(mol, &active);
    let mut systems: Vec<Vec<u32>> = biconnected_bond_components(&adj)
        .into_iter()
        .map(|bonds| {
            let mut atoms: Vec<u32> = bonds
                .iter()
                .flat_map(|&bi| {
                    let b = mol.bonds()[bi as usize];
                    [b.begin, b.end]
                })
                .collect();
            atoms.sort_unstable();
            atoms.dedup();
            atoms
        })
        .collect();
    systems.sort_by_key(|s| s.first().copied().unwrap_or(u32::MAX));
    systems
}

/// Tarjan 双连通分量(迭代实现),返回每个分量的**键下标**集合。
///
/// 只返回**含环**的分量:边数 ≥ 2 的分量必含环,单边分量是桥。
///
/// 返回键而不是原子,是因为环由边定义:相关环搜索([`crate::sssr`])需要
/// 逐分量的边集,才能把候选生成限制在分量内部 —— 否则那一步是
/// O(原子数 × 键数),整个分子上就成了平方项。原子集合由调用方从边集导出。
pub(crate) fn biconnected_bond_components(adj: &Adjacency) -> Vec<Vec<u32>> {
    let n = adj.n_atoms;
    let mut disc = vec![u32::MAX; n];
    let mut low = vec![u32::MAX; n];
    let mut timer = 0u32;
    // (起点, 终点, 键下标) —— 弹出时要的是键
    let mut edge_stack: Vec<(u32, u32, u32)> = Vec::new();
    let mut out: Vec<Vec<u32>> = Vec::new();
    let mut stack: Vec<(u32, u32, usize)> = Vec::new();

    for root in 0..n as u32 {
        if disc[root as usize] != u32::MAX {
            continue;
        }
        disc[root as usize] = timer;
        low[root as usize] = timer;
        timer += 1;
        stack.push((root, u32::MAX, 0));

        while let Some(&mut (v, from_bond, ref mut k)) = stack.last_mut() {
            let nbrs = adj.neighbors(v);
            if *k < nbrs.len() {
                let (u, bond) = nbrs[*k];
                *k += 1;
                if bond == from_bond {
                    continue;
                }
                if disc[u as usize] == u32::MAX {
                    edge_stack.push((v, u, bond));
                    disc[u as usize] = timer;
                    low[u as usize] = timer;
                    timer += 1;
                    stack.push((u, bond, 0));
                } else if disc[u as usize] < disc[v as usize] {
                    edge_stack.push((v, u, bond));
                    low[v as usize] = low[v as usize].min(disc[u as usize]);
                }
            } else {
                stack.pop();
                if let Some(&mut (parent, _, _)) = stack.last_mut() {
                    low[parent as usize] = low[parent as usize].min(low[v as usize]);
                    if low[v as usize] >= disc[parent as usize] {
                        // 弹出属于该分量的所有边
                        let mut comp: Vec<u32> = Vec::new();
                        while let Some(&(a, _, bi)) = edge_stack.last() {
                            if disc[a as usize] < disc[v as usize] {
                                break;
                            }
                            edge_stack.pop();
                            comp.push(bi);
                        }
                        if let Some(pos) = edge_stack
                            .iter()
                            .rposition(|&(a, b, _)| a == parent && b == v)
                        {
                            let (_, _, bi) = edge_stack.remove(pos);
                            comp.push(bi);
                        }
                        // 单边分量是桥,不含环
                        if comp.len() >= 2 {
                            comp.sort_unstable();
                            comp.dedup();
                            out.push(comp);
                        }
                    }
                }
            }
        }
    }
    out
}

#[cfg(test)]
mod fused_system_tests {
    use omgkit_io::smiles;

    use super::*;

    fn systems(smi: &str) -> Vec<Vec<u32>> {
        let m = smiles::parse(smi).unwrap_or_else(|e| panic!("{}", e.render()));
        fused_ring_systems(&m)
    }

    #[test]
    fn acyclic_has_no_systems() {
        assert!(systems("CCO").is_empty());
        assert!(systems("CC(C)C(=O)O").is_empty());
    }

    #[test]
    fn single_ring_is_one_system() {
        assert_eq!(systems("c1ccccc1"), vec![vec![0, 1, 2, 3, 4, 5]]);
        // 取代基不进入环系
        assert_eq!(systems("Cc1ccccc1"), vec![vec![1, 2, 3, 4, 5, 6]]);
    }

    #[test]
    fn fused_rings_merge() {
        // 萘:两个环共享一条键 → 一个环系,10 个原子
        let s = systems("c1ccc2ccccc2c1");
        assert_eq!(s.len(), 1);
        assert_eq!(s[0].len(), 10);
    }

    /// 螺环共享一个原子但**不共享键**,不算融合;螺原子恰是割点,
    /// 双连通分解天然把它们分开。
    #[test]
    fn spiro_rings_stay_separate() {
        let s = systems("C1CC12CC2");
        assert_eq!(s.len(), 2, "螺环应是两个独立环系,实际 {s:?}");
        assert!(s.iter().all(|x| x.len() == 3));
        // 螺原子(下标 2)同时属于两个环系
        assert!(s.iter().filter(|x| x.contains(&2)).count() == 2);
    }

    /// 联苯的两个环由一条**桥键**连接 —— 桥自成一个单边分量,被滤掉,
    /// 两个环因此分属不同环系。
    #[test]
    fn biphenyl_rings_stay_separate() {
        let s = systems("c1ccccc1-c1ccccc1");
        assert_eq!(s.len(), 2, "联苯应是两个独立环系,实际 {s:?}");
        assert!(s.iter().all(|x| x.len() == 6));
    }

    #[test]
    fn bridged_bicyclic_is_one_system() {
        // 双环[2.2.2]辛烷:桥环整体是一个双连通分量
        let s = systems("C1CC2CCC1CC2");
        assert_eq!(s.len(), 1);
        assert_eq!(s[0].len(), 8);
    }

    #[test]
    fn disconnected_fragments_give_separate_systems() {
        let s = systems("c1ccccc1.C1CC1");
        assert_eq!(s.len(), 2);
        assert_eq!(s[0].len(), 6);
        assert_eq!(s[1].len(), 3);
    }

    /// 环系必须恰好覆盖"在环中"的原子 —— 与桥判定互为独立校验。
    #[test]
    fn systems_cover_exactly_the_ring_atoms() {
        for smi in [
            "CCO",
            "c1ccccc1",
            "Cc1ccccc1",
            "c1ccc2ccccc2c1",
            "C1CC12CC2",
            "c1ccccc1-c1ccccc1",
            "C1CC2CCC1CC2",
            "CC(=O)Oc1ccccc1C(=O)O",
            "CN1C=NC2=C1C(=O)N(C)C(=O)N2C",
            "c1ccccc1.C1CC1",
        ] {
            let mut m = smiles::parse(smi).unwrap();
            let r = perceive_rings(&mut m);
            let mut from_systems: Vec<u32> = fused_ring_systems(&m).concat();
            from_systems.sort_unstable();
            from_systems.dedup();
            let from_bridges: Vec<u32> = r
                .atom_in_ring
                .iter()
                .enumerate()
                .filter(|(_, &x)| x)
                .map(|(i, _)| i as u32)
                .collect();
            assert_eq!(from_systems, from_bridges, "{smi}: 环系与环成员判定不一致");
        }
    }

    #[test]
    fn long_chain_does_not_overflow_stack() {
        let m = smiles::parse(&"C".repeat(20_000)).unwrap();
        assert!(fused_ring_systems(&m).is_empty());
    }
}
