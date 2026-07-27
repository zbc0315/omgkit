//! 判断立体标记是不是**真的**携带信息。
//!
//! # 为什么不能只看"邻居是否两两可区分"
//!
//! 直觉的判据是"四个取代基两两不同才算手性中心"。它对孤立的中心成立,
//! 对**相互依赖**的一对中心不成立:
//!
//! | 分子 | 两个邻居同类? | 是手性中心? |
//! |---|---|---|
//! | `[C@@H]1CCCC1` | 是(环上两条通路等价) | **否** |
//! | `O[C@H]1CC[C@@H](N)CC1` | 是(同样等价) | **是** |
//!
//! 两条分子的中心在纯对称性上完全一样。差别在于:后者环上还有**另一个**
//! 手性中心,两者合起来区分顺式与反式。单独看任一个都不成立,合起来才成立。
//!
//! 拿"两个邻居同类即非真"去删标记,会把顺式与反式塌成同一个分子 ——
//! 这是丢信息,比留一个多余的标记严重得多。
//!
//! # 判据
//!
//! 一个带四面体标记的原子**不是**真手性中心,当且仅当:
//!
//! > 它有两个邻居落在同一个对称等价类里,**且**那两条等价支路里没有别的
//! > 带四面体标记的原子。
//!
//! 后半句正是把相互依赖的那一对放行的条件。
//!
//! # 方向键也是同一类问题
//!
//! [`informative_directions`] 回答的是双键那一侧的对应问题:一条 `/` 单独
//! 说明不了任何事,它描述的是一根双键两侧取代基的相对位置。判据见那里。
//!
//! # 方向是写法,顺反是性质
//!
//! [`perceive_bond_stereo`] 把方向键换成双键**自己的**顺反属性(带参照原子)。
//! 两种表示信息等价,鲁棒性却差得远:方向依附在某根单键上,那根键被图编辑
//! 删掉,信息就跟着没了 —— 哪怕双键本身根本没被碰过。存成双键自己的属性,
//! 只要两个参照原子还在,信息就还在。
//!
//! # 为什么不在净化里
//!
//! 双键立体感知不属于净化的 12 步,由调用方在净化之后显式调用。这也解开一个
//! 分层死结:感知要用对称等价类(本层),而净化在更下层,调不到。
//!
//! # 这只是 L6 的一部分
//!
//! 完整的立体化学感知还要处理轴手性,以及"相互依赖"的更深层情形(一条支路
//! 里的中心本身又依赖于第三个中心)。
//!
//! 四面体那一条**偏保守** —— 拿不准时保留标记。保留一个多余的标记只是输出
//! 啰嗦,删错一个标记会让两个分子变成一个。方向键那一条没有这种不对称性:
//! 多写一条噪声方向会破坏规范性,少写一条真方向会丢顺反,两边都不能松。

use std::collections::{BTreeMap, BTreeSet};

use omgkit_core::{BondData, BondDirection, BondFlags, BondOrder, BondStereo, MolBuilder};

use crate::canon::symmetry_classes;

/// 逐原子判断:该原子上的四面体标记是不是真手性中心。
///
/// 没有四面体标记的原子一律为 `false`。
///
/// # 保守方向
///
/// 判不准时返回 `true`(保留)。理由见模块文档:删错标记会丢分子,
/// 留错标记只是多写几个字符。
#[must_use]
pub fn genuine_tetrahedral(mol: &MolBuilder) -> Vec<bool> {
    let n = mol.num_atoms();
    let mut out = vec![false; n];
    if n == 0 {
        return out;
    }
    let classes = symmetry_classes(mol);
    let tagged: Vec<bool> = mol
        .atoms()
        .iter()
        .map(|a| a.chiral_tag.is_tetrahedral())
        .collect();

    for a in 0..n as u32 {
        if !tagged[a as usize] {
            continue;
        }
        out[a as usize] = match equivalent_neighbour_pair(mol, a, &classes) {
            // 邻居两两可区分 —— 无条件是真手性中心
            None => true,
            // 有一对等价邻居 —— 只有当那两条支路里还有别的手性中心时才算数
            Some((x, y)) => branch_has_other_stereocentre(mol, a, x, y, &tagged),
        };
    }
    out
}

/// 直接从**存储的方向键**算一根双键的顺反,不依赖感知结果。
///
/// 返回 `(顺/反, [begin 侧参照, end 侧参照])`;这根键不是双键、或两侧没有
/// 成对的方向时返回 `None`。
///
/// # 与 [`perceive_bond_stereo`] 的分工
///
/// 那个函数要先过 [`informative_directions`] 的筛(要用对称等价类),因此只能
/// 用在**真分子**上。SMARTS 查询的拓扑里原子是占位符,算不出等价类 ——
/// 而查询里的方向是作者显式写的,本就不需要筛。所以匹配时两边都走这一支。
#[must_use]
pub fn raw_cis_trans(mol: &MolBuilder, bond: u32) -> Option<(BondStereo, [u32; 2])> {
    let db = *mol.bonds().get(bond as usize)?;
    if db.order != BondOrder::Double {
        return None;
    }
    let (ra, da) = raw_outward(mol, db.begin, db.end)?;
    let (rb, dbi) = raw_outward(mol, db.end, db.begin)?;
    let stereo = if da == dbi {
        BondStereo::Cis
    } else {
        BondStereo::Trans
    };
    Some((stereo, [ra, rb]))
}

/// `end` 那一侧带方向的邻居,方向换算到"从 `end` 向外"。
fn raw_outward(mol: &MolBuilder, end: u32, other: u32) -> Option<(u32, BondDirection)> {
    mol.neighbors(end)
        .filter(|&(o, _)| o != other)
        .find(|&(_, bi)| mol.bonds()[bi as usize].direction != BondDirection::None)
        .map(|(o, bi)| {
            let b = mol.bonds()[bi as usize];
            let d = if b.begin == end {
                b.direction
            } else {
                b.direction.flipped()
            };
            (o, d)
        })
}

/// 逐键判断:该单键上的方向(`/` `\`)是否携带信息、值得写出来。
///
/// # 为什么不能"有方向就写"
///
/// 方向键单独没有意义,它描述的是**一根双键两侧取代基的相对位置**。三种情形
/// 里的方向都是噪声:
///
/// | 输入 | 为什么不写 |
/// |---|---|
/// | `C/1CCCCC1` | 根本没有双键 |
/// | `F/C=CF` | 只有一侧有方向,说明不了相对位置 |
/// | `F/C=C(F)F` | 双键一端挂着两个相同的取代基,交换它们是自同构 |
///
/// 照写不误的代价不只是输出啰嗦:颜色细化看不见键方向,一条噪声方向能打破
/// 细化分辨不出的对称性,于是规范 SMILES 不再随重排恒定 —— 上表第一行
/// `C/1CCCCC1` 就是这样一条分子。
#[must_use]
pub fn informative_directions(mol: &MolBuilder) -> Vec<bool> {
    let mut out = vec![false; mol.num_bonds()];
    if !mol
        .bonds()
        .iter()
        .any(|b| b.direction != BondDirection::None)
    {
        return out;
    }
    let classes = symmetry_classes(mol);

    for db in mol.bonds() {
        if db.order != BondOrder::Double || db.flags.contains(BondFlags::AROMATIC) {
            continue;
        }
        // 两端都必须能区分自己的取代基,这根双键才是立体源
        if !end_is_stereogenic(mol, db.begin, db.end, &classes)
            || !end_is_stereogenic(mol, db.end, db.begin, &classes)
        {
            continue;
        }
        let left = directional_bonds_at(mol, db.begin, db.end);
        let right = directional_bonds_at(mol, db.end, db.begin);
        // 只有一侧有方向的话,相对位置无从谈起
        if left.is_empty() || right.is_empty() {
            continue;
        }
        for bi in left.into_iter().chain(right) {
            out[bi as usize] = true;
        }
    }
    out
}

/// 双键的 `end` 这一端能不能区分自己的两个取代基。
///
/// 取代基超过两个就不是普通双键端,一律判否。
fn end_is_stereogenic(mol: &MolBuilder, end: u32, other: u32, classes: &[u32]) -> bool {
    let subs: Vec<u32> = mol
        .neighbors(end)
        .map(|(o, _)| o)
        .filter(|&o| o != other)
        .collect();
    match subs.len() {
        // 剩一个取代基时,另一个位置是氢(或空),天然可区分
        0 | 1 => true,
        2 => classes[subs[0] as usize] != classes[subs[1] as usize],
        _ => false,
    }
}

/// `end` 上除通往 `other` 之外、带方向的键。
///
/// 承载方向的键**不一定是单键**:双键挂在芳香环外时(`[H]/N=c/1\\nc[nH]s1`),
/// 指方向的是那两条芳香环键。只排除双键本身。
fn directional_bonds_at(mol: &MolBuilder, end: u32, other: u32) -> Vec<u32> {
    mol.neighbors(end)
        .filter(|&(o, _)| o != other)
        .filter(|&(_, bi)| {
            let b = mol.bonds()[bi as usize];
            b.direction != BondDirection::None && b.order != BondOrder::Double
        })
        .map(|(_, bi)| bi)
        .collect()
}

/// 找一对落在同一等价类里的邻居。没有就返回 `None`。
fn equivalent_neighbour_pair(mol: &MolBuilder, a: u32, classes: &[u32]) -> Option<(u32, u32)> {
    let nbrs: Vec<u32> = mol.neighbors(a).map(|(other, _)| other).collect();
    for i in 0..nbrs.len() {
        for j in i + 1..nbrs.len() {
            if classes[nbrs[i] as usize] == classes[nbrs[j] as usize] {
                return Some((nbrs[i], nbrs[j]));
            }
        }
    }
    None
}

/// 从 `x` 与 `y` 出发、**不经过 `a`** 能到达的原子里,有没有别的四面体标记。
///
/// 环上的两条支路会汇合,于是这个集合就是整个环 —— 正好是"顺/反异构靠环上
/// 另一个中心成立"所需要的判断。
fn branch_has_other_stereocentre(
    mol: &MolBuilder,
    a: u32,
    x: u32,
    y: u32,
    tagged: &[bool],
) -> bool {
    let mut seen: BTreeSet<u32> = [a].into_iter().collect();
    let mut stack = vec![x, y];
    seen.insert(x);
    seen.insert(y);
    while let Some(cur) = stack.pop() {
        if cur != a && tagged[cur as usize] {
            return true;
        }
        for (other, _) in mol.neighbors(cur) {
            if seen.insert(other) {
                stack.push(other);
            }
        }
    }
    false
}

/// 从方向键感知双键顺反,写入 `stereo` 与 `stereo_atoms`。
///
/// 返回被标注的双键数。
///
/// # 只**新增**标注,已经有的一概不动
///
/// 反应产物里的顺反可能是搬运过来的,那时方向键早已不在,重新感知会把它抹掉。
/// 但"不动"要做到两层:既不擦成 `None`,也**不改成别的值**。
///
/// 后一层是必需的。搬运来的标注用双键**自己的参照原子**表达,与写法无关;而它
/// 旁边的方向键可能来自完全另一处 —— 反应模板写的那一根就是。共轭链上这两者会
/// 撞在一起:模板给 `C=C` 写了一对方向,其中一根**同时**贴着下一根双键,而下一
/// 根的另一侧是从底物搬过来的。拿这一对方向去重算下一根双键,等于用模板的参照系
/// 覆盖一个本来正确的答案 —— 拓扑全对,只有几何被悄悄换掉。
///
/// USPTO-50k 上这一档实测 12 条,全是共轭多烯的酰胺化/酯化。
pub fn perceive_bond_stereo(mol: &mut MolBuilder) -> usize {
    let informative = informative_directions(mol);
    if !informative.iter().any(|&x| x) {
        return 0;
    }

    let mut found: Vec<(u32, BondStereo, [u32; 2])> = Vec::new();
    for (di, db) in mol.bonds().iter().enumerate() {
        if db.order != BondOrder::Double || db.flags.contains(BondFlags::AROMATIC) {
            continue;
        }
        // 已经带着有效标注的不重算,见本函数文档
        if stereo_atoms_are_valid(mol, u32::try_from(di).unwrap_or(u32::MAX)) {
            continue;
        }
        let Some((ref_b, dir_b)) = outward_direction(mol, db.begin, db.end, &informative) else {
            continue;
        };
        let Some((ref_e, dir_e)) = outward_direction(mol, db.end, db.begin, &informative) else {
            continue;
        };
        let stereo = if dir_b == dir_e {
            BondStereo::Cis
        } else {
            BondStereo::Trans
        };
        found.push((di as u32, stereo, [ref_b, ref_e]));
    }

    let n = found.len();
    for (di, stereo, atoms) in found {
        if let Some(mut b) = mol.bond_mut(di) {
            b.set_stereo(stereo);
            b.set_stereo_atoms(atoms);
        }
    }
    n
}

/// `end` 这一端携带方向的那个邻居,以及方向换算到"**从 `end` 向外**"之后的值。
///
/// 换算是关键:存储的方向相对键自己的 `begin → end`,而这里要的是统一
/// 以双键原子为起点。少了它,同一个几何会因为两条键的存储朝向不同而
/// 时而判成顺、时而判成反。
fn outward_direction(
    mol: &MolBuilder,
    end: u32,
    other: u32,
    informative: &[bool],
) -> Option<(u32, BondDirection)> {
    mol.neighbors(end)
        .filter(|&(o, _)| o != other)
        .find(|&(_, bi)| informative[bi as usize])
        .map(|(o, bi)| {
            let b = mol.bonds()[bi as usize];
            let dir = if b.begin == end {
                b.direction
            } else {
                b.direction.flipped()
            };
            (o, dir)
        })
}

/// 参照原子仍然合法吗 —— 两个下标都在范围内,且确实是该双键两端的邻居。
///
/// 图编辑之后参照可能失效(原子被删、键被断)。写出方向键之前要先问这个。
#[must_use]
pub fn stereo_atoms_are_valid(mol: &MolBuilder, bond: u32) -> bool {
    let Some(&b) = mol.bonds().get(bond as usize) else {
        return false;
    };
    if b.stereo == BondStereo::None {
        return false;
    }
    let [ra, rb] = b.stereo_atoms;
    if ra == BondData::NO_STEREO_ATOM || rb == BondData::NO_STEREO_ATOM {
        return false;
    }
    mol.neighbors(b.begin).any(|(o, _)| o == ra) && mol.neighbors(b.end).any(|(o, _)| o == rb)
}

/// 写出时每根键该用什么方向。返回值按键下标索引。
///
/// # 优先用感知过的顺反,退回到写法
///
/// 双键上记着有效的 `stereo` 与 `stereo_atoms` 时,方向由它**重新生成** ——
/// 这样原本承载方向的那根单键即使被图编辑删掉,立体也照样写得出来。
/// 没有感知结果的双键则沿用存储的方向(仍要过
/// [`informative_directions`] 那道筛,噪声方向不写)。
///
/// # 共轭体系要一起定,不能各自为政
///
/// `F/C=C/C=C/F` 中间那根单键同时受**两根**双键约束。逐根双键各自挑方向的话,
/// 第二根会覆盖第一根的选择,写出来的顺反就错一个。
///
/// 所以按共轭片段做广度遍历:片段里第一根键的方向任意选(整体翻转不改变
/// 任何相对关系),其余的由约束推出来。
#[must_use]
pub fn directions_for_writing(mol: &MolBuilder) -> WritingDirections {
    let informative = informative_directions(mol);
    let mut out: Vec<BondDirection> = (0..mol.num_bonds())
        .map(|i| {
            if informative[i] {
                mol.bonds()[i].direction
            } else {
                BondDirection::None
            }
        })
        .collect();

    // 收集所有"有感知结果"的双键:两条参照键,以及各自的**锚点** ——
    // 锚点是该参照键上属于双键的那一端,"向外"就是从它出发。
    // 锚点在这里顺手记下,免得之后再去反查。
    let mut constraints: Vec<(Ref, Ref, BondStereo)> = Vec::new();
    for di in 0..mol.num_bonds() as u32 {
        if !stereo_atoms_are_valid(mol, di) {
            continue;
        }
        let db = mol.bonds()[di as usize];
        let (Some(b1), Some(b2)) = (
            bond_between(mol, db.begin, db.stereo_atoms[0]),
            bond_between(mol, db.end, db.stereo_atoms[1]),
        ) else {
            continue;
        };
        constraints.push((
            Ref {
                bond: b1,
                anchor: db.begin,
            },
            Ref {
                bond: b2,
                anchor: db.end,
            },
            db.stereo,
        ));
        // 这根双键**两端的所有**相邻键都要先清干净,不能只清那两条参照键。
        //
        // 沿用来的方向与生成的方向并存时,同一根双键旁边会冒出三条斜杠,
        // 读串的人只能得到另一个几何。反应产物尤其容易撞上:环上某根键从底物
        // 继承了 `/`,而参照原子是另外两个,三条方向就凑齐了。
        for end in [db.begin, db.end] {
            for (_, bi) in mol.neighbors(end) {
                if bi != di {
                    out[bi as usize] = BondDirection::None;
                }
            }
        }
    }
    // 约束图:每条参照键一个节点,同一根双键的两条参照键相连。
    // 边上带着**两端各自的锚点**与"同向还是反向"。
    //
    // 锚点必须随边走,不能只按键号记一个"向外方向" —— 共轭链里中间那根单键
    // 在相邻两根双键的约束中锚点是**不同的**两个原子(`F/C=C/C=C/F` 里是
    // C2 和 C3)。按键号存一个向外方向的话,两个约束会拿不同的参照系读写同一
    // 个值,第一根双键就被写反。所以下面统一存**存储参照系**的方向,用时换算。
    let mut adj: BTreeMap<u32, Vec<(Ref, Ref, bool)>> = BTreeMap::new();
    for &(r1, r2, stereo) in &constraints {
        let same = stereo == BondStereo::Cis;
        adj.entry(r1.bond).or_default().push((r1, r2, same));
        adj.entry(r2.bond).or_default().push((r2, r1, same));
    }

    // 存储参照系(相对键自己的 begin → end)下的方向
    let mut assigned: BTreeMap<u32, BondDirection> = BTreeMap::new();
    let outward = |mol: &MolBuilder, r: Ref, stored: BondDirection| {
        if mol.bonds()[r.bond as usize].begin == r.anchor {
            stored
        } else {
            stored.flipped()
        }
    };

    // 逐个连通片段广度遍历。片段内第一根键任取 UpRight —— 整体翻转不改变
    // 任何一对的相对关系,所以这个自由度是真的自由。
    //
    // **但它只对分子自由,对字符串不自由。** 种子按键的存储下标取,而存储下标是
    // 输入写法留下的痕迹;规范串跟着它变的话,同一个分子就有了两串。所以这里把
    // 片段编号一并交出去,由写出器按输出顺序把这个自由度定死 —— 见
    // `smiles::write` 的 `WriteStyle::Canonical`。
    let mut component = vec![None; mol.num_bonds()];
    let keys: Vec<u32> = adj.keys().copied().collect();
    let mut next_comp = 0u32;
    for start in keys {
        if assigned.contains_key(&start) {
            continue;
        }
        let comp = next_comp;
        next_comp += 1;
        assigned.insert(start, BondDirection::UpRight);
        component[start as usize] = Some(comp);
        let mut queue = vec![start];
        while let Some(cur) = queue.pop() {
            let cur_stored = assigned[&cur];
            for &(from, to, same) in &adj[&cur] {
                if assigned.contains_key(&to.bond) {
                    // 环状共轭里两条路径可能给出矛盾的约束 —— 那是输入本身
                    // 就自相矛盾。先到先得,不为一个无解的输入去猜。
                    continue;
                }
                // 换到"从 from 的锚点向外",按约束翻或不翻,再换回 to 的存储系
                let out_from = outward(mol, from, cur_stored);
                let out_to = if same { out_from } else { out_from.flipped() };
                assigned.insert(to.bond, outward(mol, to, out_to));
                component[to.bond as usize] = Some(comp);
                queue.push(to.bond);
            }
        }
    }

    for (&bi, &dir) in &assigned {
        out[bi as usize] = dir;
    }

    // 第二遍:**沿用**存储写法的那些方向键也要编上片段号。
    //
    // 它们同样有整体翻转自由度 —— 一根双键两侧的方向键**同时**取反,任何一对
    // 取代基的相对位置都不变。`informative_directions` 已经保证写出来的方向必然
    // 成对地夹着一根立体源双键,所以这里不会把孤零零一根当成一组。
    //
    // 不编号就等于不定死:同一个分子换种写法读进来,存储下标一变,规范串里的
    // `/` 与 `\` 就整体互换。反应产物是这一档的大户 —— `run_reactants` 交出来
    // 的分子还没跑过 `perceive_bond_stereo`,几何全靠沿用来的方向撑着,
    // 一根约束也没有,于是整段自由度悬空。
    let mut stack: Vec<u32> = Vec::new();
    for bi in 0..mol.num_bonds() as u32 {
        if out[bi as usize] == BondDirection::None || component[bi as usize].is_some() {
            continue;
        }
        let comp = next_comp;
        next_comp += 1;
        component[bi as usize] = Some(comp);
        stack.push(bi);
        while let Some(cur) = stack.pop() {
            for other in flanking_directions(mol, cur, &out) {
                if component[other as usize].is_none() {
                    component[other as usize] = Some(comp);
                    stack.push(other);
                }
            }
        }
    }

    WritingDirections {
        dirs: out,
        component,
    }
}

/// 与 `bond` 隔着同一根双键锁在一起的那些方向键。
///
/// 走法:从 `bond` 的任一端出发,找一根双键,再从双键的另一端取所有写了方向的
/// 单键。共轭链上这样一路传下去,整条链就是一个片段。
fn flanking_directions(mol: &MolBuilder, bond: u32, dirs: &[BondDirection]) -> Vec<u32> {
    let b = mol.bonds()[bond as usize];
    let mut out = Vec::new();
    for end in [b.begin, b.end] {
        for (far, di) in mol.neighbors(end) {
            if di == bond || mol.bonds()[di as usize].order != BondOrder::Double {
                continue;
            }
            for (_, oi) in mol.neighbors(far) {
                if oi != di && dirs[oi as usize] != BondDirection::None {
                    out.push(oi);
                }
            }
        }
    }
    out
}

/// 每根键该写什么方向,连同它属于哪个**约束片段**。
///
/// 片段是约束图的连通分量:同一片段内各键的方向互相锁死,整体翻转不改变任何一对
/// 的相对关系。那是个真自由度 —— 对分子而言。对**字符串**而言它必须被定死,否则
/// 同一个分子会有两串规范式。定死这件事要知道输出顺序,所以交给写出器做。
pub struct WritingDirections {
    /// 每根键的方向,相对该键自己的 `begin → end`。
    pub dirs: Vec<BondDirection>,
    /// 该键所属片段的编号;没写方向的键是 `None`。
    ///
    /// 由感知结果重新生成方向的键与**沿用**存储写法的键都编号,两者的整体翻转
    /// 自由度是同一回事:一根双键两侧的方向同时取反,几何不变、字符串却变了。
    pub component: Vec<Option<u32>>,
}

/// 一条参照键,连同它属于双键的那一端。
#[derive(Clone, Copy)]
struct Ref {
    bond: u32,
    /// 该键上属于双键的端点。方向的"向外"就是从这里出发。
    anchor: u32,
}

fn bond_between(mol: &MolBuilder, a: u32, b: u32) -> Option<u32> {
    mol.neighbors(a).find(|&(o, _)| o == b).map(|(_, bi)| bi)
}

/// 把每根双键的**参照原子**改选成规范秩最小的那个邻居,顺反跟着换算。
///
/// # 参照挑在哪一侧,方向符号就落在哪根键上
///
/// 参照原子记的是"相对谁说顺反"。同一根双键换一个参照、把顺反翻一次,说的是
/// 同一件事 —— 几何没变。可 [`directions_for_writing`] 是**按参照键**放方向
/// 符号的,于是换个参照,符号就落到另一根键上,串跟着变。
///
/// 而 [`perceive_bond_stereo`] 挑参照用的是"存储顺序里第一个带方向的邻居",
/// 那是**输入写法留下的痕迹**:同一个分子换种写法读进来,参照就换一个。
///
/// # 后果:规范串不是不动点
///
/// 双键一端挂着两个取代基、而其中一个恰好又是另一根双键的端点时,写出来的
/// 符号可能落在两根不同的键上,那个端点于是带了两根有方向的键。读回去时感知
/// 只认下其中一根做参照,第二次写出就少一个符号 —— 几何一模一样,串不一样。
/// 实测语料 8831 条里有 2 条。
///
/// # 修法:让写出器的输入只取决于(图, 几何)
///
/// 按**规范秩**挑参照。规范秩与输入编号无关,于是写出器看到的 `stereo_atoms`
/// 也与输入编号无关;而图与几何往返无损,`write(parse(write(m))) == write(m)`
/// 就自动成立 —— 不必再去管感知那一步挑了谁。
///
/// 只动 [`BondStereo::Cis`] / [`BondStereo::Trans`]:这两个本就是"相对记录的
/// 参照原子"说的,换参照翻一次天经地义。`Z`/`E` 按 CIP 优先级定义,与参照
/// 原子无关,换参照不该动它们。
///
/// 没有任何一根双键带有效参照时原样返回,不做拷贝 —— 绝大多数分子走这条路。
#[must_use]
pub fn normalized_stereo_refs(mol: &MolBuilder, priority: &[u32]) -> Option<MolBuilder> {
    let targets: Vec<u32> = (0..mol.num_bonds() as u32)
        .filter(|&di| {
            stereo_atoms_are_valid(mol, di)
                && matches!(
                    mol.bonds()[di as usize].stereo,
                    BondStereo::Cis | BondStereo::Trans
                )
        })
        .collect();
    if targets.is_empty() || priority.len() != mol.num_atoms() {
        return None;
    }

    let mut out = mol.clone();
    for di in targets {
        let db = mol.bonds()[di as usize];
        let pick = |end: u32, partner: u32| {
            mol.neighbors(end)
                .map(|(other, _)| other)
                .filter(|&other| other != partner)
                .min_by_key(|&other| priority[other as usize])
        };
        let (Some(x), Some(y)) = (pick(db.begin, db.end), pick(db.end, db.begin)) else {
            continue;
        };
        // 换掉一侧的参照就翻一次;两侧都换等于翻两次,回到原样
        let swaps = u8::from(x != db.stereo_atoms[0]) + u8::from(y != db.stereo_atoms[1]);
        let stereo = if swaps % 2 == 1 {
            match db.stereo {
                BondStereo::Cis => BondStereo::Trans,
                _ => BondStereo::Cis,
            }
        } else {
            db.stereo
        };
        if let Some(mut b) = out.bond_mut(di) {
            b.set_stereo(stereo);
            b.set_stereo_atoms([x, y]);
        }
    }
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles;

    fn genuine(smi: &str) -> Vec<u32> {
        let m = smiles::parse(smi).unwrap_or_else(|e| panic!("{smi}: {}", e.render()));
        genuine_tetrahedral(&m)
            .into_iter()
            .enumerate()
            .filter(|&(_, g)| g)
            .map(|(i, _)| i as u32)
            .collect()
    }

    /// 取代基两两可区分 —— 无条件是真手性中心。
    #[test]
    fn distinguishable_substituents_are_genuine() {
        assert_eq!(genuine("N[C@@H](C)C(=O)O"), vec![1], "丙氨酸");
        assert_eq!(genuine("[C@H](N)(O)F"), vec![0], "三取代 + 一个氢");
        assert_eq!(genuine("C[C@](N)(O)F"), vec![1]);
    }

    /// 有一对等价取代基、且支路里没有别的中心 —— 不是真手性中心。
    #[test]
    fn symmetric_substituents_without_help_are_not_genuine() {
        assert!(genuine("C[C@](C)(N)O").is_empty(), "两个甲基");
        assert!(genuine("[C@@H]1CCCC1").is_empty(), "对称的环戊烷");
        assert!(genuine("[C@@H]1CCCCC1").is_empty(), "对称的环己烷");
    }

    /// **相互依赖**的一对:各自看都有等价支路,合起来却区分顺式与反式。
    ///
    /// 这条是整个模块的存在理由。按"两个邻居同类即非真"处理的话,顺反两个
    /// 分子会塌成同一个。
    #[test]
    fn dependent_pair_is_genuine() {
        assert_eq!(
            genuine("O[C@H]1CC[C@@H](N)CC1"),
            vec![1, 4],
            "1,4-二取代环己烷:两个中心互相成全"
        );
        assert_eq!(
            genuine("C(#C)[C@@]1(CC[C@H](C2=CC=CC=C2)CC1)O").len(),
            2,
            "同一形状的另一个例子"
        );
    }

    /// 没有标记的原子一律为假 —— 本模块只回答"已有的标记算不算数",
    /// 不负责发现未标注的潜在手性中心。
    #[test]
    fn untagged_atoms_are_never_genuine() {
        assert!(genuine("CCO").is_empty());
        assert!(genuine("C1CCCCC1").is_empty());
        assert!(genuine("c1ccccc1").is_empty());
    }

    /// 判据对输入编号不敏感 —— 它建立在对称等价类之上,而等价类本身是
    /// 重排不变的。
    #[test]
    fn verdict_is_renumbering_invariant() {
        // 同一个分子的两种写法
        for (a, b) in [
            ("O[C@H]1CC[C@@H](N)CC1", "N[C@H]1CC[C@@H](O)CC1"),
            ("N[C@@H](C)C(=O)O", "OC(=O)[C@@H](N)C"),
            ("[C@@H]1CCCC1", "C1CC[C@@H]1C"),
        ] {
            let ga = genuine(a).len();
            let gb = genuine(b).len();
            assert_eq!(ga, gb, "{a} 与 {b} 的真手性中心个数应当一致");
        }
    }

    fn perceive(smi: &str) -> Vec<(u32, u32, BondStereo, [u32; 2])> {
        let mut m = smiles::parse(smi).unwrap_or_else(|e| panic!("{smi}: {}", e.render()));
        perceive_bond_stereo(&mut m);
        m.bonds()
            .iter()
            .filter(|b| b.stereo != BondStereo::None)
            .map(|b| (b.begin, b.end, b.stereo, b.stereo_atoms))
            .collect()
    }

    /// 最基本的一对:同一个骨架,只有斜杠方向不同。
    #[test]
    fn cis_and_trans_are_distinguished() {
        assert_eq!(
            perceive("F/C=C/F"),
            vec![(1, 2, BondStereo::Trans, [0, 3])],
            "反式"
        );
        assert_eq!(
            perceive("F/C=C\\F"),
            vec![(1, 2, BondStereo::Cis, [0, 3])],
            "顺式"
        );
    }

    /// 换一种写法写同一个分子,结论必须一样。
    ///
    /// 这条守的是"向外换算"那一步 —— 两种写法里方向键的存储朝向相反,
    /// 少了换算就会一个判顺、一个判反。
    #[test]
    fn verdict_does_not_depend_on_how_it_was_written() {
        let a = perceive("F/C=C/F");
        let b = perceive("C(\\F)=C/F");
        assert_eq!(a.len(), 1);
        assert_eq!(b.len(), 1);
        assert_eq!(a[0].2, b[0].2, "F/C=C/F 与 C(\\F)=C/F 是同一个分子");
    }

    /// 四取代双键:参照原子是**携带方向的**那两个邻居,不是随便挑的。
    #[test]
    fn reference_atoms_are_the_ones_carrying_the_directions() {
        let got = perceive("C/C(F)=C(F)/C");
        assert_eq!(got.len(), 1);
        let (_, _, stereo, refs) = got[0];
        assert_eq!(stereo, BondStereo::Trans);
        assert_eq!(refs, [0, 5], "两个甲基才是参照,不是那两个氟");
    }

    /// 共轭多烯:每根双键各自标注。
    #[test]
    fn each_double_bond_gets_its_own_verdict() {
        let got = perceive("F/C=C/C=C/F");
        assert_eq!(got.len(), 2);
        assert!(got.iter().all(|g| g.2 == BondStereo::Trans));
    }

    /// 方向不携带信息时不标注 —— 与写出侧同一套判据。
    #[test]
    fn non_informative_directions_are_not_perceived() {
        assert!(perceive("C/1CCCCC1").is_empty(), "根本没有双键");
        assert!(perceive("F/C=CF").is_empty(), "只有一侧有方向");
        assert!(perceive("F/C=C(F)F").is_empty(), "一端两个取代基相同");
        assert!(perceive("FC=CF").is_empty(), "没有方向键");
    }

    /// 参照原子的有效性检查。
    #[test]
    fn validity_follows_the_graph() {
        let mut m = smiles::parse("F/C=C/F").unwrap();
        perceive_bond_stereo(&mut m);
        let db = m
            .bonds()
            .iter()
            .position(|b| b.stereo != BondStereo::None)
            .expect("有标注") as u32;
        assert!(stereo_atoms_are_valid(&m, db));

        // 没有标注的键一律无效
        let other = (0..m.num_bonds() as u32)
            .find(|&i| i != db)
            .expect("有别的键");
        assert!(!stereo_atoms_are_valid(&m, other));
    }

    /// 感知之后写出,方向键要能原样再现。
    #[test]
    fn perceived_stereo_regenerates_directions() {
        for smi in ["F/C=C/F", "F/C=C\\F", "C/C=C/C", "C/C(F)=C(F)/C"] {
            let mut m = smiles::parse(smi).unwrap();
            perceive_bond_stereo(&mut m);
            let w = smiles::write(&m).smiles;
            let mut back =
                smiles::parse(&w).unwrap_or_else(|e| panic!("{smi} → {w}: {}", e.render()));
            perceive_bond_stereo(&mut back);
            let a: Vec<_> = m
                .bonds()
                .iter()
                .map(|b| b.stereo)
                .filter(|s| *s != BondStereo::None)
                .collect();
            let b: Vec<_> = back
                .bonds()
                .iter()
                .map(|b| b.stereo)
                .filter(|s| *s != BondStereo::None)
                .collect();
            assert_eq!(a, b, "{smi} → {w}:顺反没有守恒");
        }
    }

    /// **共轭多烯**:中间那根单键同时受两根双键约束。
    ///
    /// 逐根双键各自挑方向的话,第二根会覆盖第一根的选择,写出来就错一个。
    /// 这条专门守传播那一步 —— 单根双键的用例覆盖不到它。
    #[test]
    fn conjugated_chain_directions_stay_consistent() {
        for smi in [
            "F/C=C/C=C/F",
            "F/C=C\\C=C/F",
            "F/C=C/C=C\\F",
            "C/C=C/C=C/C=C/C",
        ] {
            let mut m = smiles::parse(smi).unwrap();
            let n = perceive_bond_stereo(&mut m);
            assert!(n >= 2, "{smi}:应当标注到多根双键,实际 {n}");
            let w = smiles::write(&m).smiles;
            let mut back = smiles::parse(&w).unwrap();
            perceive_bond_stereo(&mut back);
            let a: Vec<_> = m
                .bonds()
                .iter()
                .filter(|b| b.stereo != BondStereo::None)
                .map(|b| (b.stereo, b.stereo_atoms))
                .collect();
            let b: Vec<_> = back
                .bonds()
                .iter()
                .filter(|b| b.stereo != BondStereo::None)
                .map(|b| (b.stereo, b.stereo_atoms))
                .collect();
            assert_eq!(a, b, "{smi} → {w}:共轭链上的顺反没有全部守恒");
        }
    }

    /// 承载方向的那根键被删掉之后,立体仍然写得出来。
    ///
    /// 这是把方向换成"双键自己的属性"的**全部理由**。删掉一根带方向的单键,
    /// 若还按老办法从 direction 写,信息就没了。
    #[test]
    fn stereo_survives_deletion_of_the_direction_bearing_bond() {
        let mut m = smiles::parse("C/C=C/F").unwrap();
        perceive_bond_stereo(&mut m);
        // 把甲基那根方向键的 direction 抹掉,模拟"承载方向的键被图编辑动过"
        for i in 0..m.num_bonds() as u32 {
            if let Some(mut b) = m.bond_mut(i) {
                b.set_direction(BondDirection::None);
            }
        }
        let w = smiles::write(&m).smiles;
        assert!(
            w.contains('/') || w.contains('\\'),
            "写成了 {w} —— direction 抹掉之后立体就丢了,说明写出还在读 direction"
        );
    }

    /// 感知过顺反的双键,两端**所有**相邻键上的陈旧方向都要清掉。
    ///
    /// 只清那两条参照键是不够的:图编辑之后,双键旁边可能还留着从别处继承来的
    /// 方向。它与新生成的方向并存时,同一根双键旁边就有三条斜杠,读串的人
    /// 得到的是另一个几何。反应产物最容易撞上这种情形。
    #[test]
    fn stale_directions_next_to_a_perceived_double_bond_are_cleared() {
        let mut m = smiles::parse("C/C=C/F").unwrap();
        perceive_bond_stereo(&mut m);
        // 在双键的另一个邻居上塞一条与参照无关的陈旧方向
        let db = (0..m.num_bonds() as u32)
            .find(|&i| m.bonds()[i as usize].stereo != BondStereo::None)
            .expect("有标注");
        let ends = [m.bonds()[db as usize].begin, m.bonds()[db as usize].end];
        let extra = ends
            .iter()
            .flat_map(|&e| m.neighbors(e).map(|(_, bi)| bi).collect::<Vec<_>>())
            .find(|&bi| bi != db && m.bonds()[bi as usize].direction == BondDirection::None);
        if let Some(bi) = extra {
            if let Some(mut b) = m.bond_mut(bi) {
                b.set_direction(BondDirection::DownRight);
            }
        }
        let dirs = directions_for_writing(&m).dirs;
        let n = dirs.iter().filter(|d| **d != BondDirection::None).count();
        assert_eq!(n, 2, "一根双键只该有两条方向键,实际 {n} 条");
    }
}
