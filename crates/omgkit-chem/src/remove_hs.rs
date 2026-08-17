//! 把图里的显式氢原子并成邻居上的氢计数。
//!
//! # 不属于净化
//!
//! 这一步**改变原子数**,是一次独立的图编辑,不在净化的 12 步之内。净化只改
//! 属性不动图,那条性质是"在解析结果的图上直接做纯图算法"成立的前提,不能
//! 为了合并氢就破掉。
//!
//! # 为什么需要它
//!
//! 真实反应数据库里的 SMILES 常把氢写成独立原子(`[H][C@]12...`)。显式氢会
//! 把邻接原子的度数撑大 —— 一个本该是 `D3` 的碳变成 `D4`,写着 `D3` 的模板
//! 就配不上它。合并之后度数回到化学上该有的值,模板才匹配得上。
//!
//! # 总价不变,所以不必重算隐式氢
//!
//! 邻居少了一根键(价 −1)、多了一个显式氢计数(价 +1),显式价净变化为 0,
//! 隐式氢因此**保持正确**。这也是为什么本函数不去碰 `num_implicit_hs`:
//! 碰了反而要重算一遍价键,而重算会把"净化过与否"的状态搅乱。
//!
//! # 保留哪些氢:判错的代价是不对称的
//!
//! 多留一个氢只是图里多个节点,分子还是对的;删错一个氢会**丢信息** ——
//! 同位素标记、电荷、立体参照都可能一并没了,而且丢得很安静。所以判据偏保守:
//! 拿不准就留着。具体见 [`is_removable`]。

use omgkit_core::{AtomFlags, BondData, BondDirection, BondOrder, ChiralTag, MolBuilder};

/// 把可以合并的显式氢并进邻居的氢计数,返回删掉的氢原子数。
///
/// # 原子下标会变
///
/// 删原子必然重排下标,所以本函数**重建**整个分子。调用方手上任何原子/键
/// 下标在返回之后都失效了 —— 这一点没法在类型上拦住,只能靠这句话。
///
/// 邻居的相对顺序**保持不变**:重建时按原键序重新建键,而邻居顺序等于建键
/// 顺序。手性标记依赖这个顺序,顺序一乱标记就全错了。
pub fn remove_hs(mol: &mut MolBuilder) -> usize {
    let n = mol.num_atoms();
    let doomed: Vec<bool> = (0..n as u32).map(|a| is_removable(mol, a)).collect();
    let n_removed = doomed.iter().filter(|&&d| d).count();
    if n_removed == 0 {
        return 0;
    }

    // 旧下标 → 新下标
    let mut new_idx = vec![u32::MAX; n];
    let mut out = MolBuilder::with_capacity(n - n_removed, mol.num_bonds());

    for a in 0..n as u32 {
        if doomed[a as usize] {
            continue;
        }
        let mut data = mol.atoms()[a as usize];
        // 被删掉的氢落在这个原子的第几个邻居位上
        let positions: Vec<usize> = mol
            .neighbors(a)
            .enumerate()
            .filter(|(_, (other, _))| doomed[*other as usize])
            .map(|(pos, _)| pos)
            .collect();
        if !positions.is_empty() {
            let merged = u8::try_from(positions.len()).unwrap_or(u8::MAX);
            // 氢该进哪个槽,取决于宿主是不是"氢数我说了算"的方括号原子。
            //
            // NO_IMPLICIT 立着时隐式氢恒为 0,总氢数就是 num_explicit_hs,
            // 合并进来的氢只能记在那里(`[nH]` 这类正是如此)。
            //
            // 没立时隐式氢由价推出来,合并进 num_explicit_hs 会有个副作用:
            // 写出器见到非零的显式氢就必须加方括号,于是乙醇被写成
            // `C[CH2][OH]` —— 分子没错,却凭空多了一层括号。记进隐式槽就没这事,
            // 而且价的账也对得上:宿主少一根键(−1)、多一个隐式氢(+1)。
            if data.flags.contains(AtomFlags::NO_IMPLICIT) {
                data.num_explicit_hs = data.num_explicit_hs.saturating_add(merged);
            } else {
                data.num_implicit_hs = data.num_implicit_hs.saturating_add(merged);
            }
            data.chiral_tag = rebased_tag(data.chiral_tag, &positions);
        }
        new_idx[a as usize] = out.add_atom_data(data);
    }

    // 按**原键序**重建,邻居的相对顺序才不会变
    for b in mol.bonds() {
        if doomed[b.begin as usize] || doomed[b.end as usize] {
            continue;
        }
        let mut nb = *b;
        nb.begin = new_idx[b.begin as usize];
        nb.end = new_idx[b.end as usize];
        nb.stereo_atoms = [
            translate(b.stereo_atoms[0], &new_idx),
            translate(b.stereo_atoms[1], &new_idx),
        ];
        let _ = out.add_bond_data(nb);
    }

    if let Some(name) = mol.name() {
        out.set_name(name.to_string());
    }
    *mol = out;
    n_removed
}

/// 参照原子的下标换算。哨兵值不换算 —— 它不是下标。
fn translate(idx: u32, new_idx: &[u32]) -> u32 {
    if idx == BondData::NO_STEREO_ATOM {
        return BondData::NO_STEREO_ATOM;
    }
    new_idx
        .get(idx as usize)
        .copied()
        .filter(|&v| v != u32::MAX)
        .unwrap_or(BondData::NO_STEREO_ATOM)
}

/// 删掉若干邻居位上的氢之后,四面体标记要不要翻。
///
/// # 括号氢在参照系里的位置是固定的
///
/// 标记相对**邻居存储顺序**。氢写成图里的原子时,它就在存储序的某个位置 k;
/// 写成括号里的计数时,按 SMILES 的约定它落在"紧跟着前一个邻居"的位置,
/// 也就是概念上的**第 1 位**(解析器对"手性原子恰好是串首"那种写法做的那次
/// 翻转,正是把第 0 位搬到第 1 位)。
///
/// 所以合并一个氢 = 把它从第 k 位搬到第 1 位,需要 `|k − 1|` 次相邻对换,
/// 奇数次就翻。
///
/// # 一个中心上删掉两个以上的氢时不动标记
///
/// 那时中心至少挂着两个氢,两个相同的取代基交换是自同构 —— 它本就不是手性
/// 中心,标记没有内容可言。硬翻一次只是把一个无意义的值换成另一个。
fn rebased_tag(tag: ChiralTag, removed_positions: &[usize]) -> ChiralTag {
    if !tag.is_tetrahedral() || removed_positions.len() != 1 {
        return tag;
    }
    let k = removed_positions[0];
    if k.abs_diff(1) % 2 == 1 {
        tag.inverted()
    } else {
        tag
    }
}

/// 这个原子是不是"可以并进邻居"的显式氢。
///
/// 判据偏保守,每一条挡的都是一类会**丢信息**的删除:
///
/// | 不删的情形 | 删了会丢什么 |
/// |---|---|
/// | 不是氢 | —— |
/// | 邻居数不是 1 | 0 个:并不进谁;≥2 个:桥氢,并给谁都是猜 |
/// | 到邻居的键不是普通单键 | 配位键、芳香键上的氢不是普通取代氢 |
/// | 键上带方向(`/` `\`) | 那根键正是双键顺反的载体,删了顺反就没了 |
/// | 它是某根双键的立体参照 | 同上,参照没了顺反无从表达 |
/// | 有同位素 | 氘、氚是**另一种核素**,并成氢计数就分不出来了 |
/// | 带电荷 | `[H+]` 是质子,一个独立的物种 |
/// | 有映射号 | 反应模板按号引用它 |
/// | 带自由基电子 | 氢自由基携带信息 |
/// | 邻居也是氢 | 氢分子:两个都删就什么都不剩了 |
/// | 邻居是通配原子 | 把氢计数并进"任意原子"没有意义 |
#[must_use]
pub fn is_removable(mol: &MolBuilder, atom: u32) -> bool {
    let Some(&a) = mol.atoms().get(atom as usize) else {
        return false;
    };
    if a.atomic_num != 1 {
        return false;
    }
    if a.isotope != 0
        || a.formal_charge != 0
        || a.atom_map != 0
        || a.num_radical_electrons != 0
        || a.flags.contains(AtomFlags::AROMATIC)
    {
        return false;
    }
    // 自己还挂着氢计数的氢:那是 `[HH]` 之类,合并不了
    if a.num_explicit_hs != 0 {
        return false;
    }

    let mut it = mol.neighbors(atom);
    let Some((other, bond)) = it.next() else {
        return false; // 孤立的氢,并不进谁
    };
    if it.next().is_some() {
        return false; // 桥氢
    }

    let b = mol.bonds()[bond as usize];
    if b.order != BondOrder::Single || b.direction != BondDirection::None {
        return false;
    }

    let host = mol.atoms()[other as usize];
    if host.atomic_num == 1 || host.atomic_num == 0 {
        return false;
    }

    // 被任何一根键当作立体参照的氢都留着
    !mol.bonds()
        .iter()
        .any(|bb| bb.stereo_atoms[0] == atom || bb.stereo_atoms[1] == atom)
}
