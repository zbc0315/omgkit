//! 芳香性感知(净化第 7 步)。
//!
//! 输入是[相关环集合](crate::sssr),输出是原子与键上的芳香标志。
//!
//! # 判定流程
//!
//! 1. 给每个原子定**电子给体类型**:它能向 π 体系贡献几个电子
//! 2. 筛出**候选环** —— 环上每个原子都够格参与芳香
//! 3. 把共享**恰好一条键**的候选环归为一组(融合体系)
//! 4. 组内枚举 1..=6 个环的**连通**子集,对每个子集做 Hückel 判定
//! 5. 判定通过的子集,把它的**边界键**标为芳香
//!
//! # 为什么要枚举子集而不是只看单个环
//!
//! 稠合体系里,单个环可能不满足 4n+2,而若干环的并集满足。反过来也有:
//! 并集不满足但单个环满足。两种情况都要覆盖,所以从 1 个环开始逐级放大。
//!
//! 上限取 6 个环:再大的组合在真实分子里没有化学意义,而组合数是 C(n,k),
//! 放开会指数爆炸。
//!
//! # 三个容易写错的地方
//!
//! **1. "融合"要求共享的键**恰好一条**。** 共享两条及以上的(笼状体系里常见)
//! 不算融合,不进同一组 —— 否则会把立体的笼子当成平面 π 体系。
//!
//! **2. 参与 Hückel 计数的原子,只取出现在子集中 1 个或 2 个环里的。**
//! 出现在 3 个及以上环里的是体系的中心原子,它的电子不在环流里。
//!
//! **3. 只有**边界键**被标芳香。** 在当前子集中出现两次的键是融合键,
//! 由更小的子集(单个环)那一轮负责标记。逐级放大的过程保证不会漏。

use omgkit_core::{element, AtomFlags, BondFlags, BondOrder, MolBuilder};

use crate::sssr::{ring_set, Ring};
use crate::valence::explicit_valence_nonstrict;

/// 组合中最多容纳的环数
const MAX_FUSED_RINGS: usize = 6;
/// 参与融合分组的环的最大尺寸(按键数)。更大的环只能单独接受判定。
const MAX_FUSED_RING_SIZE: usize = 24;

/// 原子能向 π 体系贡献的电子数。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Donor {
    /// 不贡献,也没有空轨道 —— 该原子不能参与芳香
    None,
    /// 有空 p 轨道,贡献 0 个电子
    Vacant,
    One,
    Two,
    /// 身份未知(通配原子),1 或 2 个都可以
    Any,
}

impl Donor {
    /// 贡献电子数的取值区间
    fn range(self) -> (i32, i32) {
        match self {
            Self::One => (1, 1),
            Self::Two => (2, 2),
            Self::Any => (1, 2),
            Self::None | Self::Vacant => (0, 0),
        }
    }

    /// 该类型是否允许原子参与芳香
    fn can_be_aromatic(self) -> bool {
        !matches!(self, Self::None)
    }
}

/// 感知芳香性,就地设置原子与键的芳香标志。返回判定为芳香的环数。
///
/// 调用前必须已完成 kekulize 与自由基赋值 —— 电子计数要读具体的键级和
/// 自由基数。
pub fn set_aromaticity(mol: &mut MolBuilder) -> usize {
    let rings = ring_set(mol);
    if rings.is_empty() {
        return 0;
    }

    let n_atoms = mol.num_atoms();
    let donors: Vec<Donor> = (0..n_atoms as u32).map(|i| donor_type(mol, i)).collect();
    let candidate: Vec<bool> = (0..n_atoms as u32)
        .map(|i| donors[i as usize].can_be_aromatic() && is_arom_candidate(mol, i))
        .collect();

    // 候选环:环上每个原子都够格,且不能整环都是通配原子
    let cand_rings: Vec<&Ring> = rings
        .iter()
        .filter(|r| {
            r.atoms.iter().all(|&a| candidate[a as usize])
                && r.atoms
                    .iter()
                    .any(|&a| mol.atoms()[a as usize].atomic_num != 0)
        })
        .collect();
    if cand_rings.is_empty() {
        return 0;
    }

    let neigh = ring_neighbors(&cand_rings, mol.num_bonds());
    let mut n_arom = 0usize;
    let mut done = vec![false; cand_rings.len()];
    // 暂存空间只分配一次,所有环组共用 —— 逐组分配是 O(组数 × 分子规模)
    let mut scratch = Scratch::new(n_atoms, mol.num_bonds());

    for start in 0..cand_rings.len() {
        if done[start] {
            continue;
        }
        let group = collect_group(start, &neigh, &mut done);
        n_arom += huckel_on_group(mol, &cand_rings, &group, &neigh, &donors, &mut scratch);
    }
    n_arom
}

/// 跨环组复用的暂存空间。
///
/// 每一项都按分子规模定长,逐组重新分配就是 O(组数 × 分子规模) —— 一个隐式
/// 的平方项。所以只分配一次,用完**只清理自己碰过的下标**。
struct Scratch {
    /// 每个原子在当前子集中出现于几个环
    ring_count: Vec<u32>,
    /// 键是否已被标为芳香(用于"全部环键都标完就提前收工")
    marked_bonds: Vec<bool>,
    /// 逐组统计环键总数时的去重标记
    bond_seen: Vec<bool>,
}

impl Scratch {
    fn new(n_atoms: usize, n_bonds: usize) -> Self {
        Self {
            ring_count: vec![0; n_atoms],
            marked_bonds: vec![false; n_bonds],
            bond_seen: vec![false; n_bonds],
        }
    }
}

// ---------------------------------------------------------------------------
// 电子计数

/// 原子可向 π 体系贡献的电子数;`None` 表示该原子根本不能参与。
///
/// 共轭判定([`crate::conjugation`])也用它 —— 两处判据共享同一套电子计数。
pub(crate) fn count_pi_electrons(mol: &MolBuilder, idx: u32) -> Option<i32> {
    let atom = mol.atoms()[idx as usize];
    let z = atom.atomic_num;
    let dv = element::by_atomic_num(z)
        .and_then(|e| e.valences.first().copied())
        .map_or(-1, i32::from);
    // 单价元素既不能芳香也不能共轭
    if dv <= 1 {
        return None;
    }

    // 有效配位数:重原子邻居 + 全部氢,再扣掉不贡献价的键(配位键的给体端)
    let mut degree = mol.degree(idx) as i32 + i32::from(total_hs(mol, idx));
    for (_, bi) in mol.neighbors(idx) {
        if mol.bonds()[bi as usize].valence_contribution_to(idx) == 0.0 {
            degree -= 1;
        }
    }
    // 配位数超过 3 就不可能是平面 π 体系的一员
    if degree > 3 {
        return None;
    }

    let n_outer = i32::from(element::by_atomic_num(z).map_or(0, |e| e.outer_electrons));
    // 孤对电子数 = 外层电子 - 默认价,再扣掉形式电荷
    let n_lone_pairs = (n_outer - dv - i32::from(atom.formal_charge)).max(0);
    let n_radicals = i32::from(atom.num_radical_electrons);

    let mut res = (dv - degree) + n_lone_pairs - n_radicals;
    if res > 1 {
        // 有三键及以上时只算 1 个电子。用总不饱和度探测 ——
        // 候选判据已经排除了"多于一处不饱和"的原子,所以这里的 >1 必是高阶键。
        if unsaturations(mol, idx) > 1 {
            res = 1;
        }
    }
    Some(res)
}

/// 总氢数 = 显式 + 隐式
fn total_hs(mol: &MolBuilder, idx: u32) -> u8 {
    let a = mol.atoms()[idx as usize];
    a.num_explicit_hs + a.num_implicit_hs
}

/// 不饱和度 = 显式价 − 重原子度数
fn unsaturations(mol: &MolBuilder, idx: u32) -> i32 {
    explicit_valence_nonstrict(mol, idx) - mol.degree(idx) as i32
}

/// 该原子是否有一条**环外**的多重键;有则返回另一端的原子
fn exocyclic_multiple_bond(mol: &MolBuilder, idx: u32) -> Option<u32> {
    mol.neighbors(idx).find_map(|(other, bi)| {
        let b = mol.bonds()[bi as usize];
        (!b.flags.contains(BondFlags::IN_RING) && b.valence_contribution_to(idx) >= 2.0)
            .then_some(other)
    })
}

/// 该原子是否有一条**环内**的多重键
fn cyclic_multiple_bond(mol: &MolBuilder, idx: u32) -> bool {
    mol.neighbors(idx).any(|(_, bi)| {
        let b = mol.bonds()[bi as usize];
        b.flags.contains(BondFlags::IN_RING) && b.valence_contribution_to(idx) >= 2.0
    })
}

/// 该原子是否有任何多重键
fn has_multiple_bond(mol: &MolBuilder, idx: u32) -> bool {
    let mut deg = mol.degree(idx) as i32 + i32::from(mol.atoms()[idx as usize].num_explicit_hs);
    for (_, bi) in mol.neighbors(idx) {
        if mol.bonds()[bi as usize]
            .valence_contribution_to(idx)
            .round() as i32
            == 0
        {
            deg -= 1;
        }
    }
    explicit_valence_nonstrict(mol, idx) != deg
}

/// 电负性比较:外层电子多者更强;相同则原子序数小者更强。
fn more_electronegative(z1: u8, z2: u8) -> bool {
    let e = |z: u8| element::by_atomic_num(z).map_or(0, |x| x.outer_electrons);
    let (n1, n2) = (e(z1), e(z2));
    n1 > n2 || (n1 == n2 && z1 < z2)
}

fn donor_type(mol: &MolBuilder, idx: u32) -> Donor {
    let atom = mol.atoms()[idx as usize];
    // 通配原子身份未知:环内有多重键时至少能贡献 1 个,否则 1 或 2 都行
    if atom.atomic_num == 0 {
        return if cyclic_multiple_bond(mol, idx) {
            Donor::One
        } else {
            Donor::Any
        };
    }

    let Some(nelec) = count_pi_electrons(mol, idx) else {
        return Donor::None;
    };
    if nelec < 0 {
        return Donor::None;
    }

    if nelec == 0 {
        // 没有电子可给,但可能有空 p 轨道
        if exocyclic_multiple_bond(mol, idx).is_some() {
            Donor::Vacant
        } else if cyclic_multiple_bond(mol, idx) {
            Donor::One
        } else {
            Donor::None
        }
    } else if nelec == 1 {
        if let Some(other) = exocyclic_multiple_bond(mol, idx) {
            // 唯一的那个电子被环外多重键占用。若对端电负性更强,
            // 电子被拉走,只剩空轨道。
            let z_other = mol.atoms()[other as usize].atomic_num;
            if more_electronegative(z_other, atom.atomic_num) {
                Donor::Vacant
            } else {
                Donor::One
            }
        } else if has_multiple_bond(mol, idx) {
            Donor::One
        } else if atom.formal_charge == 1 {
            // 卓鎓离子、环丙烯正离子:空轨道
            Donor::Vacant
        } else {
            Donor::None
        }
    } else {
        let mut nelec = nelec;
        if let Some(other) = exocyclic_multiple_bond(mol, idx) {
            let z_other = mol.atoms()[other as usize].atomic_num;
            if more_electronegative(z_other, atom.atomic_num) {
                nelec -= 1;
            }
        }
        if nelec % 2 == 1 {
            Donor::One
        } else {
            Donor::Two
        }
    }
}

// ---------------------------------------------------------------------------
// 候选判据

/// 该原子是否够格参与芳香(不含电子给体类型的判断)。
fn is_arom_candidate(mol: &MolBuilder, idx: u32) -> bool {
    let atom = mol.atoms()[idx as usize];
    let z = atom.atomic_num;

    // 只允许周期表**前三行**(Z ≤ 18),外加 Se 与 Te。
    //
    // 参照实现 `Aromaticity.cpp:422-424` 的注释写的是"前两行",而它的代码同样
    // 是 `z > 18` —— **那句注释在上游就是错的**,两边的行为一致。照抄行为,
    // 不照抄那句话。
    if z > 18 && z != 34 && z != 52 {
        return false;
    }

    // 价态偏离默认值的原子出局
    let default_valence = |zz: u8| {
        element::by_atomic_num(zz)
            .and_then(|e| e.valences.first().copied())
            .map_or(-1, i32::from)
    };
    let dv = default_valence(z);
    if dv > 0 {
        let eff_z = (i32::from(z) - i32::from(atom.formal_charge)).clamp(0, 118) as u8;
        let total_valence = explicit_valence_nonstrict(mol, idx) + i32::from(total_hs(mol, idx))
            - i32::from(atom.num_explicit_hs);
        if total_valence > default_valence(eff_z) {
            return false;
        }
    }

    // 带自由基的杂原子、带电的碳自由基都出局
    if atom.num_radical_electrons > 0 && (z != 6 || atom.formal_charge != 0) {
        return false;
    }

    // 不允许一个原子带**多于一根**双键或三键 —— 如 `C1=C=NC=N1` 里的累积双键
    if unsaturations(mol, idx) > 1 {
        let n_mult = mol
            .neighbors(idx)
            .filter(|&(_, bi)| {
                matches!(
                    mol.bonds()[bi as usize].order,
                    BondOrder::Double | BondOrder::Triple
                )
            })
            .count();
        if n_mult > 1 {
            return false;
        }
    }

    true
}

// ---------------------------------------------------------------------------
// 环的分组与组合枚举

/// 候选环之间的邻接:共享**恰好一条**键才算相邻。
///
/// 共享两条及以上的出现在笼状体系里,那不是平面融合。尺寸超过
/// [`MAX_FUSED_RING_SIZE`] 的环不参与分组,只能单独接受判定。
///
/// 用"键 → 含它的环"倒排来找候选环对,而不是两两比较 —— 后者是
/// O(环数²),在多环分子上就是平方项。
fn ring_neighbors(rings: &[&Ring], n_bonds: usize) -> Vec<Vec<usize>> {
    let mut bond_rings: Vec<Vec<u32>> = vec![Vec::new(); n_bonds];
    for (i, r) in rings.iter().enumerate() {
        if r.bonds.len() > MAX_FUSED_RING_SIZE {
            continue;
        }
        for &b in &r.bonds {
            bond_rings[b as usize].push(i as u32);
        }
    }

    let mut shared: std::collections::HashMap<(u32, u32), u32> = std::collections::HashMap::new();
    for rs in &bond_rings {
        for a in 0..rs.len() {
            for b in (a + 1)..rs.len() {
                *shared.entry((rs[a], rs[b])).or_default() += 1;
            }
        }
    }

    let mut out = vec![Vec::new(); rings.len()];
    for ((i, j), count) in shared {
        if count == 1 {
            out[i as usize].push(j as usize);
            out[j as usize].push(i as usize);
        }
    }
    // 哈希表遍历顺序不定,排序保证结果确定
    for v in &mut out {
        v.sort_unstable();
    }
    out
}

/// 从 `start` 出发收集一个连通的环组(深度优先,迭代实现)。
fn collect_group(start: usize, neigh: &[Vec<usize>], done: &mut [bool]) -> Vec<usize> {
    let mut group = Vec::new();
    let mut stack = vec![start];
    done[start] = true;
    while let Some(r) = stack.pop() {
        group.push(r);
        for &nb in &neigh[r] {
            if !done[nb] {
                done[nb] = true;
                stack.push(nb);
            }
        }
    }
    group.sort_unstable();
    group
}

/// 逐级生成组内的**连通**子集,每一级内按字典序排列。
///
/// 直接枚举 C(n,k) 再过滤连通性,在长链稠合体系上是组合爆炸:16 个环的并苯
/// 要试 14892 个组合,其中连通的只有 81 个。这里改为从上一级生长出下一级,
/// 只产出连通子集。
///
/// **必须逐级惰性生成。** 判定过程会在"全部环键都已标为芳香"时提前收工;
/// 若一次性把 1..6 级全算出来再遍历,那个收工判断一分钱也省不下 —— 而高级别
/// 的子集数是爆炸的:一张 25×25 的稠合薄片,6 级合计 53 万个子集,可它在
/// 第 1 级就已经收工了。
///
/// 参数 `neigh` 按**组内位置**索引,不是环号。
struct ConnectedSubsets<'a> {
    neigh: &'a [Vec<usize>],
    max_size: usize,
    /// 当前这一级
    cur: Vec<Vec<usize>>,
    size: usize,
}

impl<'a> ConnectedSubsets<'a> {
    fn new(neigh: &'a [Vec<usize>], max_size: usize) -> Self {
        Self {
            neigh,
            max_size,
            cur: Vec::new(),
            size: 0,
        }
    }

    /// 生成下一级;返回 `None` 表示已到顶或再无连通子集。
    fn next_level(&mut self) -> Option<&[Vec<usize>]> {
        if self.size >= self.max_size {
            return None;
        }
        self.size += 1;

        if self.size == 1 {
            self.cur = (0..self.neigh.len()).map(|i| vec![i]).collect();
            return Some(&self.cur);
        }

        let mut next: std::collections::BTreeSet<Vec<usize>> = std::collections::BTreeSet::new();
        for subset in &self.cur {
            for &p in subset {
                for &nb in &self.neigh[p] {
                    if subset.binary_search(&nb).is_ok() {
                        continue;
                    }
                    let mut t = Vec::with_capacity(subset.len() + 1);
                    t.extend_from_slice(subset);
                    t.push(nb);
                    t.sort_unstable();
                    next.insert(t);
                }
            }
        }
        if next.is_empty() {
            return None;
        }
        self.cur = next.into_iter().collect();
        Some(&self.cur)
    }
}

// ---------------------------------------------------------------------------
// Hückel 判定

/// 对一个融合环组做判定,返回其中判为芳香的环数。
fn huckel_on_group(
    mol: &mut MolBuilder,
    rings: &[&Ring],
    group: &[usize],
    neigh: &[Vec<usize>],
    donors: &[Donor],
    scratch: &mut Scratch,
) -> usize {
    // 把邻接换算成**组内位置**索引;group 已排序,可用二分
    let local_neigh: Vec<Vec<usize>> = group
        .iter()
        .map(|&r| {
            let mut v: Vec<usize> = neigh[r]
                .iter()
                .filter_map(|x| group.binary_search(x).ok())
                .collect();
            v.sort_unstable();
            v
        })
        .collect();

    // 本组的环键总数 —— 全部标完就可以提前收工
    let mut n_ring_bonds = 0usize;
    for &r in group {
        for &b in &rings[r].bonds {
            if !scratch.bond_seen[b as usize] {
                scratch.bond_seen[b as usize] = true;
                n_ring_bonds += 1;
            }
        }
    }
    for &r in group {
        for &b in &rings[r].bonds {
            scratch.bond_seen[b as usize] = false;
        }
    }

    let mut aromatic_rings: Vec<usize> = Vec::new();
    let mut n_done_bonds = 0usize;
    let mut touched: Vec<u32> = Vec::new();
    let mut counted: Vec<u32> = Vec::new();

    let mut levels = ConnectedSubsets::new(&local_neigh, MAX_FUSED_RINGS);
    let mut subset: Vec<usize> = Vec::new();
    while n_done_bonds < n_ring_bonds {
        let Some(level) = levels.next_level() else {
            break;
        };
        for positions in level {
            subset.clear();
            subset.extend(positions.iter().map(|&p| group[p]));

            touched.clear();
            for &r in &subset {
                for &a in &rings[r].atoms {
                    if scratch.ring_count[a as usize] == 0 {
                        touched.push(a);
                    }
                    scratch.ring_count[a as usize] += 1;
                }
            }
            // 参与计数的原子:在子集中恰好出现 1 次或 2 次的
            counted.clear();
            counted.extend(
                touched
                    .iter()
                    .copied()
                    .filter(|&a| matches!(scratch.ring_count[a as usize], 1 | 2)),
            );
            let aromatic = huckel(&counted, donors);
            for &a in &touched {
                scratch.ring_count[a as usize] = 0;
            }

            if aromatic {
                mark_aromatic(
                    mol,
                    rings,
                    &subset,
                    &mut scratch.marked_bonds,
                    &mut n_done_bonds,
                );
                aromatic_rings.extend_from_slice(&subset);
            }
        }
    }

    // 已标记集合是**每组一份**的:本组用完就把自己那些键复位,否则下一组会把
    // "上一组已经标过"误当成"本组已经标过",提前收工的判断永远不触发。
    for &r in group {
        for &b in &rings[r].bonds {
            scratch.marked_bonds[b as usize] = false;
        }
    }

    aromatic_rings.sort_unstable();
    aromatic_rings.dedup();
    aromatic_rings.len()
}

/// Hückel 4n+2 判定。
fn huckel(atoms: &[u32], donors: &[Donor]) -> bool {
    let (mut low, mut up) = (0i32, 0i32);
    let mut n_any = 0usize;
    for &a in atoms {
        let d = donors[a as usize];
        if d == Donor::Any {
            n_any += 1;
            // 多于一个身份未知的原子,判定失去意义
            if n_any > 1 {
                return false;
            }
        }
        let (lo, hi) = d.range();
        low += lo;
        up += hi;
    }

    if up >= 6 {
        // 区间里存在满足 4n+2 的取值即可
        (low..=up).any(|e| (e - 2) % 4 == 0)
    } else {
        up == 2
    }
}

/// 把子集的**边界键**及其两端标为芳香。
///
/// 在子集里出现两次的键是融合键,不在此处标记 —— 它由只含其中一个环的
/// 更小子集负责。逐级放大的枚举保证不会漏。
fn mark_aromatic(
    mol: &mut MolBuilder,
    rings: &[&Ring],
    subset: &[usize],
    marked_bonds: &mut [bool],
    n_done_bonds: &mut usize,
) {
    let mut count = std::collections::BTreeMap::<u32, usize>::new();
    for &r in subset {
        for &b in &rings[r].bonds {
            *count.entry(b).or_default() += 1;
        }
    }
    for (bi, c) in count {
        if c != 1 {
            continue;
        }
        let bond = mol.bonds()[bi as usize];
        if !marked_bonds[bi as usize] {
            marked_bonds[bi as usize] = true;
            *n_done_bonds += 1;
        }
        if let Some(mut b) = mol.bond_mut(bi) {
            b.flags_mut().insert(BondFlags::AROMATIC);
            // 只有单键和双键转成芳香键;三键之类保持原样
            if matches!(bond.order, BondOrder::Single | BondOrder::Double) {
                b.set_order(BondOrder::Aromatic);
            } else {
                continue;
            }
        }
        for a in [bond.begin, bond.end] {
            if let Some(at) = mol.atom_mut(a) {
                at.flags.insert(AtomFlags::AROMATIC);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use omgkit_io::smiles;

    use super::*;
    use crate::{assign_radicals, clean_up, kekulize, perceive_rings, update_property_cache};

    /// 跑完整的第 1/3/4/5/6/7 步前缀,返回 (芳香原子, 芳香键)
    fn perceive(smi: &str) -> (Vec<bool>, Vec<bool>) {
        let mut m = smiles::parse(smi).unwrap_or_else(|e| panic!("{}", e.render()));
        clean_up(&mut m);
        update_property_cache(&mut m).expect("价键校验应通过");
        let _ = perceive_rings(&mut m);
        kekulize(&mut m).expect("应能 kekulize");
        assign_radicals(&mut m);
        set_aromaticity(&mut m);
        (
            m.atoms()
                .iter()
                .map(|a| a.flags.contains(AtomFlags::AROMATIC))
                .collect(),
            m.bonds()
                .iter()
                .map(|b| b.flags.contains(BondFlags::AROMATIC))
                .collect(),
        )
    }

    fn n_aromatic_atoms(smi: &str) -> usize {
        perceive(smi).0.iter().filter(|&&x| x).count()
    }

    #[test]
    fn benzene_is_aromatic() {
        let (a, b) = perceive("c1ccccc1");
        assert!(a.iter().all(|&x| x), "苯的全部碳都应芳香");
        assert!(b.iter().all(|&x| x), "苯的全部键都应芳香");
    }

    /// 用凯库勒式写的苯同样应被判为芳香 —— 芳香性是**感知**出来的,
    /// 不是从输入里读出来的。
    #[test]
    fn kekule_input_is_perceived_as_aromatic() {
        assert_eq!(n_aromatic_atoms("C1=CC=CC=C1"), 6);
        assert_eq!(n_aromatic_atoms("C1=CC=NC=C1"), 6, "凯库勒式吡啶");
    }

    #[test]
    fn classic_heteroaromatics() {
        for (smi, n, name) in [
            ("c1ccncc1", 6, "吡啶"),
            ("c1cc[nH]c1", 5, "吡咯"),
            ("c1ccoc1", 5, "呋喃"),
            ("c1ccsc1", 5, "噻吩"),
            ("c1cnc[nH]1", 5, "咪唑"),
            ("c1ccc2ccccc2c1", 10, "萘"),
            ("c1ccc2[nH]ccc2c1", 9, "吲哚"),
        ] {
            assert_eq!(n_aromatic_atoms(smi), n, "{name}");
        }
    }

    #[test]
    fn saturated_rings_are_not_aromatic() {
        for smi in ["C1CCCCC1", "C1CCNCC1", "C1CO1", "C1CCC1"] {
            assert_eq!(n_aromatic_atoms(smi), 0, "{smi}");
        }
    }

    /// 环状但不满足 4n+2 的体系不该被判为芳香。
    #[test]
    fn antiaromatic_and_nonaromatic_rings() {
        assert_eq!(n_aromatic_atoms("C1=CC=C1"), 0, "环丁二烯:4 电子");
        assert_eq!(n_aromatic_atoms("C1=CCCCCCC1"), 0, "环辛烯");
        assert_eq!(n_aromatic_atoms("O=C1C=CC(=O)C=C1"), 0, "对苯醌");
    }

    /// 环丙烯正离子(2 电子)与卓鎓离子(6 电子):都靠空 p 轨道成芳香。
    #[test]
    fn cationic_aromatics() {
        assert_eq!(n_aromatic_atoms("[cH+]1cc1"), 3, "环丙烯正离子");
        assert_eq!(n_aromatic_atoms("[cH+]1cccccc1"), 7, "卓鎓离子");
    }

    /// 环外双键会把电子拉走:环己酮不芳香,而环戊二烯酮同理。
    #[test]
    fn exocyclic_double_bonds_steal_electrons() {
        assert_eq!(n_aromatic_atoms("O=C1CCCCC1"), 0);
        assert_eq!(n_aromatic_atoms("O=C1C=CC=C1"), 0);
    }

    /// 融合体系:边界键与融合键最终都应标为芳香。
    #[test]
    fn fused_systems_mark_all_bonds() {
        for (smi, name) in [
            ("c1ccc2ccccc2c1", "萘"),
            ("c1ccc2c(c1)ccc1ccccc12", "菲"),
            ("c1cc2ccc3cccc4ccc(c1)c2c34", "芘"),
        ] {
            let (a, b) = perceive(smi);
            assert!(a.iter().all(|&x| x), "{name}: 应全部原子芳香");
            assert!(b.iter().all(|&x| x), "{name}: 应全部键芳香,含融合键");
        }
    }

    /// 部分芳香:只有苯环那一半是芳香的。
    #[test]
    fn partially_aromatic_molecules() {
        let (a, _) = perceive("c1ccccc1C1CCCCC1");
        assert_eq!(a.iter().filter(|&&x| x).count(), 6, "只有苯环芳香");
        assert!(a[..6].iter().all(|&x| x));
        assert!(a[6..].iter().all(|&x| !x));
    }

    /// 幂等:跑两遍与跑一遍结果相同。
    #[test]
    fn is_idempotent() {
        for smi in ["c1ccccc1", "c1ccc2ccccc2c1", "c1cc[nH]c1", "C1CCCCC1"] {
            let mut m = smiles::parse(smi).unwrap();
            clean_up(&mut m);
            update_property_cache(&mut m).unwrap();
            let _ = perceive_rings(&mut m);
            kekulize(&mut m).unwrap();
            assign_radicals(&mut m);

            set_aromaticity(&mut m);
            let once: Vec<_> = m.atoms().iter().map(|a| a.flags).collect();
            set_aromaticity(&mut m);
            let twice: Vec<_> = m.atoms().iter().map(|a| a.flags).collect();
            assert_eq!(once, twice, "{smi}: 不幂等");
        }
    }

    /// 无环分子完全不受影响。
    #[test]
    fn acyclic_molecules_are_untouched() {
        for smi in ["CCO", "C=CC=C", "N#CC=O"] {
            assert_eq!(n_aromatic_atoms(smi), 0, "{smi}");
        }
    }
}
