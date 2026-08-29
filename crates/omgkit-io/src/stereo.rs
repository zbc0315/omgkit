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

use omgkit_core::{
    BondData, BondDirection, BondFlags, BondOrder, BondStereo, ChiralTag, MolBuilder,
};

use crate::wedge::Wedge;

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

/// 顺/反能成立的**最小环**。环再小一点,双键上那两条环内通路就被环锁成顺式,
/// `/` 与 `\` 写什么都不改变几何。
///
/// 八元是化学事实,不是拍脑袋的阈值:反式环辛烯是能分离出来的最小反式环烯,
/// 七元及以下的反式在几何上根本搭不起来。外部实现(RDKit `MinBondRingSize < 8`)
/// 用的是同一条线 —— 实测四到七元一律不标、八元起才标。
const MIN_STEREOGENIC_RING: usize = 8;

/// 这根键落在小于 [`MIN_STEREOGENIC_RING`] 元的环里吗。
///
/// 环的大小 = 这根键 + 绕开它的最短通路,所以只要问"不走这根键,能不能在
/// `MIN_STEREOGENIC_RING - 2` 步之内从一端走到另一端"。搜索深度是常数,
/// 不需要成环感知,也就不必把 `omgkit-chem` 拽进这一层。
fn in_small_ring(mol: &MolBuilder, bond: u32) -> bool {
    let Some(&b) = mol.bonds().get(bond as usize) else {
        return false;
    };
    let mut seen: BTreeSet<u32> = [b.begin].into_iter().collect();
    let mut frontier = vec![b.begin];
    for _ in 0..MIN_STEREOGENIC_RING - 2 {
        let mut next = Vec::new();
        for &cur in &frontier {
            for (other, bi) in mol.neighbors(cur) {
                if bi == bond {
                    continue;
                }
                if other == b.end {
                    return true;
                }
                if seen.insert(other) {
                    next.push(other);
                }
            }
        }
        if next.is_empty() {
            return false;
        }
        frontier = next;
    }
    false
}

/// 直接从**存储的方向键**算一根双键的顺反,不依赖感知结果。
///
/// 返回 `(顺/反, [begin 侧参照, end 侧参照])`;这根键不是双键、**落在八元以下
/// 的环里**(顺反由环定死,见模块内的 `MIN_STEREOGENIC_RING`)、或两侧没有
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
    if db.order != BondOrder::Double || in_small_ring(mol, bond) {
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

/// 这个分子**写明了双键几何,却还没感知过顺反**。
///
/// 也就是:某根双键两侧的方向键(`/` `\`)成对写着,而它自己的
/// [`BondData::stereo`](omgkit_core::BondData::stereo) 还是 `None` ——
/// 调用方漏了 [`perceive_bond_stereo`]。
///
/// # 为什么这一步值得单独有个谓词
///
/// 感知**不在** `omgkit_chem::pipeline::sanitize` 里(它要用对称等价类,
/// 那在净化的上一层,调不到)。漏了不会报错,只会**静默丢几何**:
///
/// - 反应那边:方向键依附在某根单键上,反应一旦删掉那根键,几何就跟着没了 ——
///   产物照样合法、原子数照样对,只有顺反悄悄少了;
/// - 画图那边:每根双键的 `stereo` 都是 `None`,顺反校正整个空转,
///   E/Z 可能画反 —— 而线条本身看着一点毛病没有。
///
/// 两处都是"画错了/算错了"而不是"没做好",所以调用方要在 `debug_assert!` 里
/// 拿它当场拦住。**返回 `false` 不等于分子有顺反** —— 根本没写方向键的分子
/// 也返回 `false`,那是对的:没写就没什么可丢的。
///
/// # 它与 [`perceive_bond_stereo`] **由构造保证一致**
///
/// 两者走同一个 `would_annotate`。这不是洁癖:先前这里判的是
/// "[`raw_cis_trans`] 给得出结果",而那一条**少了 [`informative_directions`] 那道过滤** ——
/// 于是"感知有意不标"的分子被报成"感知没跑过"。
///
/// 实测三条合法 SMILES 会误报,而且**再调一次感知也不会变绿**(它本来就不该标):
///
/// | SMILES | 感知标了几根 | 旧谓词 |
/// |---|---|---|
/// | `F/C=C(\F)F` | 0 | **报错** |
/// | `Cl/C=C(\Cl)Cl` | 0 | **报错** |
/// | `CC(/C=C/C)=C(\C)C` | 1 | **报错** |
///
/// 第三条尤其现实:一根真顺反 + 一根冗余方向,公共库 SMILES 的常见写法。
/// 后果不是"多数一个" —— `omgkit-depict` 的 `generate` 与 `omgkit-match`
/// 的重定基都拿它做 `debug_assert!`,于是 debug 下画这三个分子**直接 panic**。
///
/// 光看"有方向键却没有顺反"本来就不够:一条孤零零的 `/` 说明不了任何事,
/// 一端挂着两个等价取代基时方向同样是噪声(见 [`informative_directions`] 的表)。
/// 判据必须与感知问同一个问题,否则它守的是另一件事。
#[must_use]
pub fn directions_not_perceived(mol: &MolBuilder) -> bool {
    let informative = informative_directions(mol);
    if !informative.iter().any(|&x| x) {
        return false;
    }
    (0..mol.num_bonds()).any(|i| {
        mol.bonds()[i].stereo == BondStereo::None && would_annotate(mol, i, &informative).is_some()
    })
}

/// **感知会不会给这根键留下顺反?** —— [`perceive_bond_stereo`] 与
/// [`directions_not_perceived`] 共用的那一条判断。
///
/// 抽出来是为了让"会不会标"与"标成什么"只有一份实现:两处各写一遍的话,
/// 迟早在某个过滤条件上分岔,而分岔的表现是判据报一件代码没做的事。
fn would_annotate(
    mol: &MolBuilder,
    di: usize,
    informative: &[bool],
) -> Option<(BondStereo, [u32; 2])> {
    let db = *mol.bonds().get(di)?;
    if db.order != BondOrder::Double || db.flags.contains(BondFlags::AROMATIC) {
        return None;
    }
    // 小环里的双键没有顺反可言,见 [`MIN_STEREOGENIC_RING`]
    if in_small_ring(mol, u32::try_from(di).ok()?) {
        return None;
    }
    // 已经带着有效标注的不重算,见 `perceive_bond_stereo` 的文档
    if stereo_atoms_are_valid(mol, u32::try_from(di).ok()?) {
        return None;
    }
    let (ref_b, dir_b) = outward_direction(mol, db.begin, db.end, informative)?;
    let (ref_e, dir_e) = outward_direction(mol, db.end, db.begin, informative)?;
    let stereo = if dir_b == dir_e {
        BondStereo::Cis
    } else {
        BondStereo::Trans
    };
    Some((stereo, [ref_b, ref_e]))
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
/// | `C/1=C\CCCC1` | 双键在八元以下的环里,顺反由环定死(`MIN_STEREOGENIC_RING`)|
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

    for (di, db) in mol.bonds().iter().enumerate() {
        if db.order != BondOrder::Double || db.flags.contains(BondFlags::AROMATIC) {
            continue;
        }
        // 小环里的双键不是立体源,挂在它两侧的方向是噪声 —— 见
        // [`MIN_STEREOGENIC_RING`]。少了这一条,规范串会把一条不携带任何
        // 信息的 `/` 写出去,而外部实现不写。
        if in_small_ring(mol, u32::try_from(di).unwrap_or(u32::MAX)) {
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

/// 双键的 `end` 这一端,哪些邻居可以充当**参照原子**。返回 `(原子, 键)`。
///
/// 参照可以经由**非单键**到达:双键挂在芳香环外时(`[H]/N=c/1\\nc[nH]s1`),
/// 担这个角色的是那两条芳香环键。只排除双键 —— `C=C=C` 那种累积双键的
/// "另一侧"本身又是一根双键,它定的是轴手性([`ChiralTag::Allene`]),
/// 不是顺反。
///
/// # 为什么这条筛选只能有一处
///
/// 从方向键感知(`/` `\\`)与从二维坐标反读([`assign_bond_stereo_2d`])问的是
/// 同一件事:这一端拿谁当参照。两处各写一份的话,一条路认下的参照另一条路
/// 认不下 —— 于是同一个分子读进来、写出去,顺反悄悄换了个参照系。
///
/// # 排除双键这一条**没有判据碰得到**
///
/// 拿掉它,语料级判据与单元测试全绿 —— 两条路各有一层挡在前面:
///
/// * 二维那条路,真实的丙二烯画出来是**直线**,那根双键的另一端正好落在轴上,
///   [`cis_trans_from_points`] 先一步判"读不出来";
/// * SMILES 那条路,`/` `\\` 按语法只写在单键/芳香键上,双键身上根本不会有方向。
///
/// 也就是说这一条挡的是**两条路都到不了**的输入。留着是因为它写明了
/// "累积双键定的是轴手性、不是顺反";但它守住了什么,眼下没人量得出来 ——
/// 照实记,免得全绿被读成"这一条也验过了"。
fn reference_neighbours(mol: &MolBuilder, end: u32, other: u32) -> Vec<(u32, u32)> {
    mol.neighbors(end)
        .filter(|&(o, _)| o != other)
        .filter(|&(_, bi)| mol.bonds()[bi as usize].order != BondOrder::Double)
        .collect()
}

/// `end` 上除通往 `other` 之外、带方向的键。
fn directional_bonds_at(mol: &MolBuilder, end: u32, other: u32) -> Vec<u32> {
    reference_neighbours(mol, end, other)
        .into_iter()
        .filter(|&(_, bi)| mol.bonds()[bi as usize].direction != BondDirection::None)
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

    // **判断哪些键该标,与 [`directions_not_perceived`] 共用 [`would_annotate`]** ——
    // 两处各写一遍的话迟早在某个过滤条件上分岔,而那正是先前发生过的事。
    let mut found: Vec<(u32, BondStereo, [u32; 2])> = Vec::new();
    for di in 0..mol.num_bonds() {
        if let Some((stereo, atoms)) = would_annotate(mol, di, &informative) {
            found.push((u32::try_from(di).unwrap_or(u32::MAX), stereo, atoms));
        }
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

/// 按楔形给分子里每个读得出构型的中心打上手性标记,返回打了几个。
///
/// # 必须在**净化之后**调
///
/// 判一个中心要先知道它有几个氢([`crate::wedge::chirality_from_wedges`] 的 `(3, 1)` 那一支
/// 就是"三根键 + 一个隐式氢"),而隐式氢数是净化算出来的。刚从文件读出来的
/// 分子那一栏还是 0 —— 那时调这个函数,带隐式氢的中心会被整档漏掉,
/// 而且一声不响。
///
/// 顺序因此是:**读文件(L1)→ 净化(L2)→ 回来打标记(L1)**。跨层来回一趟
/// 看着别扭,但把净化搬进 L1 是更坏的选择,而把这一步搬进 L2 会让"楔形怎么读"
/// 有两个住处。
///
/// # 只管二维
///
/// 三维 molblock 的立体在**坐标本身**里,楔形一般是空的 —— 那一档这里一个也
/// 打不出来,得另走一条路。名字里带 `_2d` 就是为了不让人误以为它两种都管。
///
/// # 认出三维就整个不做
///
/// "楔形一般是空的"只是**一般**。三维文件里偶尔留着楔形字段,那时
/// [`crate::wedge::chirality_from_wedges`] 会把 z 一丢、按 xy 投影算体积 —— 算得出一个答案,
/// 而那个答案与分子无关。空答案可以接受,错答案不行,所以任何一个 `z` 不为零
/// 就返回 0。这与顺反那一侧
/// ([`crate::stereo::assign_bond_stereo_2d`])是同一条线。
#[must_use]
pub fn assign_chirality_2d(mol: &mut MolBuilder, coords: &[[f64; 3]], wedges: &[Wedge]) -> usize {
    if coords.iter().any(|p| p[2].abs() > FLAT_TOL) {
        return 0;
    }
    let mut n = 0;
    for a in 0..u32::try_from(mol.num_atoms()).unwrap_or(0) {
        let Some(tag) = crate::wedge::chirality_from_wedges(mol, coords, wedges, a) else {
            continue;
        };
        if let Some(at) = mol.atom_mut(a) {
            at.chiral_tag = tag;
            n += 1;
        }
    }
    n
}

/// 从**三维坐标**给每个四面体中心打上手性标记,返回打了几个。
///
/// 坐标是平的(所有 `z` 都是 0)时一个也不打 —— 那是二维图,走
/// [`assign_chirality_2d`]。四个 `assign_*` 都按这条线分工。
///
/// # 号的约定:**以中心原子为基点**,正是 `@`(逆时针)
///
/// 体积取 `det[l₀−c, l₁−c, l₂−c]`,三个配体按**存储顺序**的前三个。
/// 这与 `omgkit_conf::chiral::center_volume` 是同一个量、同一个号 ——
/// 那边把它标定过:`@` 对应正、`@@` 对应负,真实构象上 127/127 与 120/120。
/// 两处的一致由 `omgkit-conf` 里的
/// `the_center_volume_convention_agrees_with_reading_it_back` 钉住。
///
/// **不用四配体那个行列式**(`det[l₁−l₀, l₂−l₀, l₃−l₀]`,楔形那条路在用)。
/// 它完全不看中心原子在哪:中心被挤到配体四面体外面(伞形翻转)时它一点不变,
/// 而真实构型已经翻了。读真实三维结构必须用中心基点这个 —— RDKit 的
/// `assignChiralTypesFrom3D` 同此。
///
/// # 三配位 + 隐式氢 / 孤对:用同样的三个点
///
/// 那一档的槽位约定是 `[n₀, 看不见的那个, n₁, n₂]`。把看不见的那个从槽位 1
/// 挪到槽位 3 是个 3-轮换(**偶置换,号不变**),于是前三个槽位就是
/// `n₀, n₁, n₂` —— 与四配位那一档同一个式子。看不见的配体在哪不必知道。
///
/// # 三个方向张不出体积就不打,而这把尺量的是**没归一化**的体积
///
/// 单位是 Å³。楔形那条路量的是**单位**方向的三重积(那边的坐标以键长为单位),
/// 这条路不是:三维文件的坐标就是 Å,而外部实现的三维立体感知也按 Å³ 卡同一个
/// `0.1`。两边用同一把尺才谈得上"别人读得回来" —— 换成归一化的话,几乎压平的
/// 中心两边判得不一样。实测语料里正好有一个:亚砜的 S 嵌得几乎共面,
/// 单位三重积 −0.061(我方会拒),原始体积 0.268(外部实现照读)。
///
/// 压平的中心是真实存在的,那时构型无从谈起,猜一个就是编。
///
/// # 只打**真手性中心**,而且分两趟打
///
/// 三维坐标给每个 sp3 原子都算得出一个号 —— 甲基、亚甲基一个不落。那些号不是
/// 立体化学,是噪声,而且**有害**:[`genuine_tetrahedral`] 判"支路里有没有别的
/// 中心"时会把噪声当成中心,于是叔丁基那种明摆着不是中心的原子被判成真中心
/// (实测就撞上了)。
///
/// 可光按"邻居两两分得开"去筛又太狠:环上那种**两条支路组成相同、靠彼此才
/// 成立**的中心会被一起筛掉(4-取代环己烷上那一对)。实测这么筛,语料的分歧
/// 从 17 条涨到 85 条。
///
/// 所以按 [`genuine_tetrahedral`] 那条规矩来,只是**从上往下削**而不是从下往上补:
///
/// 1. 几何算得出号的原子先全当候选;
/// 2. 反复剔掉"有一对等价邻居、而那两条支路里没有别的候选"的;
/// 3. 削到不再变为止。
///
/// 从下往上补是行不通的:环上那种**两个中心互相依赖**的(1,4-二取代环己烷
/// 那一对)谁也当不了种子,补法一个都补不出来 —— 实测那样分歧从 17 涨到 47。
/// 从上往下削则天然处理它:两个互为对方支路里的候选,谁也削不掉。
///
/// 削法对噪声同样有效,而且是逐层的:甲基先掉(两条支路是死路),叔丁基跟着掉
/// (它的等价支路里只剩已经掉了的甲基)。
#[must_use]
pub fn assign_chirality_3d(mol: &mut MolBuilder, coords: &[[f64; 3]]) -> usize {
    if !coords.iter().any(|p| p[2].abs() > FLAT_TOL) {
        return 0;
    }
    let classes = symmetry_classes(mol);
    let na = u32::try_from(mol.num_atoms()).unwrap_or(0);
    // 几何先算好:哪些原子的坐标定得出一个号。是不是真中心是另一回事。
    let from_geometry: Vec<Option<ChiralTag>> = (0..na)
        .map(|a| chirality_from_coords(mol, coords, a))
        .collect();

    let mut tagged: Vec<bool> = from_geometry.iter().map(Option::is_some).collect();
    // 只有"有一对等价邻居"的才可能被削掉,其余无条件成立 —— 先挑出来,
    // 免得每一轮都把全部原子重算一遍等价对。
    let mut shaky: Vec<(u32, u32, u32)> = (0..na)
        .filter(|&a| tagged[a as usize])
        .filter_map(|a| equivalent_neighbour_pair(mol, a, &classes).map(|(x, y)| (a, x, y)))
        .collect();
    loop {
        let before = shaky.len();
        let mut doomed = Vec::new();
        shaky.retain(|&(a, x, y)| {
            if branch_has_other_stereocentre(mol, a, x, y, &tagged) {
                true
            } else {
                doomed.push(a);
                false
            }
        });
        for a in doomed {
            tagged[a as usize] = false;
        }
        if shaky.len() == before {
            break;
        }
    }

    let mut n = 0;
    for a in 0..na {
        if !tagged[a as usize] {
            continue;
        }
        let Some(tag) = from_geometry[a as usize] else {
            continue;
        };
        if let Some(at) = mol.atom_mut(a) {
            at.chiral_tag = tag;
            n += 1;
        }
    }
    n
}

/// 一个中心在三维坐标上的构型。判不出来返回 `None`。
fn chirality_from_coords(mol: &MolBuilder, coords: &[[f64; 3]], a: u32) -> Option<ChiralTag> {
    let nbrs: Vec<u32> = mol.neighbors(a).map(|(n, _)| n).collect();
    let hs = crate::wedge::total_hs(mol, a);
    // 与楔形那条路同一套资格判断:四根键、或三根键 + 一个看不见的第四配体
    // (隐式氢,或亚砜/膦上那对孤对)。
    match (nbrs.len(), hs) {
        (4, 0) => {}
        (3, 1) => {}
        (3, 0) if crate::wedge::has_lone_pair(mol, a) => {}
        _ => return None,
    }

    let c = *coords.get(a as usize)?;
    let mut dirs = [[0.0_f64; 3]; 3];
    for (k, d) in dirs.iter_mut().enumerate() {
        let p = *coords.get(*nbrs.get(k)? as usize)?;
        // **不归一化** —— 见 `assign_chirality_3d` 里那把尺的理由。
        *d = [p[0] - c[0], p[1] - c[1], p[2] - c[2]];
    }
    let (u, v, w) = (dirs[0], dirs[1], dirs[2]);
    let vol = u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0])
        + u[2] * (v[0] * w[1] - v[1] * w[0]);
    if vol.abs() <= crate::wedge::ZERO_VOLUME_TOL {
        return None;
    }
    Some(if vol > 0.0 {
        ChiralTag::Ccw
    } else {
        ChiralTag::Cw
    })
}

/// 参照原子落在双键轴上多近就算"读不出来"。单位是坐标单位的平方(叉积)。
///
/// 卡在纯粹的退化上:二维图里参照原子共线是布局病态,不是化学。留一条窄缝
/// 是为了让浮点噪声不至于把同侧判成异侧。
const AXIS_TOL: f64 = 1e-9;

/// `z` 大到多少就算"这不是一张平面图"。单位是**坐标单位**,与 [`AXIS_TOL`]
/// 那个面积量不是一回事,所以另立一个常数 —— 数值撞在一起是巧合。
/// 与 `molblock` 读取器判 `is_3d` 用的是同一条线。
pub(crate) const FLAT_TOL: f64 = 1e-9;

/// 两个参照原子在双键两侧还是同侧 —— **同侧为顺**。判不出来返回 `None`。
///
/// 只吃四个点,不吃分子:调用方各自挑好参照原子之后,剩下的就是一道平面几何。
///
/// # 这条符号约定只能有一处实现
///
/// 画图那边把画反了的双键掰回来(`omgkit_depict::stereo::read_bond_stereo`),
/// 读文件这边从坐标反读顺反([`assign_bond_stereo_2d`]) —— 两处问的是同一个
/// 几何问题。各写一份的话,一处判顺、另一处判反,**而且两边各自自洽**:
/// 画出来的图和从图里读回来的分子是一对相反的顺反,谁也不报错。
#[must_use]
pub fn cis_trans_from_points(
    begin: [f64; 2],
    end: [f64; 2],
    ref_begin: [f64; 2],
    ref_end: [f64; 2],
) -> Option<BondStereo> {
    let (dx, dy) = (end[0] - begin[0], end[1] - begin[1]);
    let side = |p: [f64; 2]| dx * (p[1] - begin[1]) - dy * (p[0] - begin[0]);
    let (sa, sb) = (side(ref_begin), side(ref_end));
    if sa.abs() < AXIS_TOL || sb.abs() < AXIS_TOL {
        return None;
    }
    Some(if sa * sb > 0.0 {
        BondStereo::Cis
    } else {
        BondStereo::Trans
    })
}

/// 这根双键是不是立体源,是的话两侧各挑一个参照原子。
///
/// **只管资格与参照,不碰几何** —— 二维那条路(投影同侧/异侧)与三维那条路
/// (二面角)接在它后面。资格判断各写一遍的话,同一根键会在一种坐标下算立体源、
/// 另一种下不算,而那种差别不报错。
/// 这根双键**几何上分得出顺反,而分子里还没有顺反信息**。
///
/// 两处在用,而且必须是同一个判断:
///
/// * 从坐标读立体时,它圈出"值得去量一量"的那批键([`stereo_candidate`]);
/// * 往文件里写时,它圈出**必须标成交叉双键**的那批键
///   ([`unspecified_cis_trans`])。
///
/// 两处若各写一套,同一根键就会「读的时候算立体源、写的时候不算」—— 那正好是
/// 凭空造出构型的形状:写出去没标"未知",读回来量出一个确定的值。
fn geometry_could_decide(mol: &MolBuilder, di: u32, classes: &[u32]) -> bool {
    let Some(&db) = mol.bonds().get(di as usize) else {
        return false;
    };
    if db.order != BondOrder::Double || db.flags.contains(BondFlags::AROMATIC) {
        return false;
    }
    // 什么算立体源,与方向键那条路是同一套:小环、已有标注两条来自
    // `would_annotate`,两端能否区分取代基那条来自 `informative_directions`
    // (`would_annotate` 通过它的筛选结果间接用上)。各写一遍的话,同一根键
    // 从 SMILES 进来标、从文件进来不标。
    if in_small_ring(mol, di) {
        return false;
    }
    if stereo_atoms_are_valid(mol, di) {
        return false;
    }
    end_is_stereogenic(mol, db.begin, db.end, classes)
        && end_is_stereogenic(mol, db.end, db.begin, classes)
}

/// 逐键:**写进文件时该不该标成交叉双键**(V2000 键块第四列的 `3`,含义是
/// 「顺反未知」)。
///
/// # 不标的后果是凭空造出一个构型
///
/// 二维图也好、三维构象也好,**画出来的每根双键都必然有一个确定的几何** ——
/// 布局算法总得把取代基摆在某一侧。作者没写顺反的那些键,摆完之后从图上就量
/// 得出一个值了。不标交叉的话,读的一方(别人的工具箱,以及我们自己)会把那个
/// 值当成化学信息读走。
///
/// 实测:大语料 8831 个分子里 **551 个(6.2%)** 写出去再读回来会多出顺反标记,
/// 而原串里一个都没有。方向是**造信息**,比丢信息更难发现 —— 拿到的文件看不出
/// 任何毛病,只是多了一句作者从没说过的话。
///
/// # 只标"分得出顺反"的那些
///
/// 苯环里的双键、两端取代基相同的双键、小环内的双键都不标:它们本来就没有
/// 顺反可言,标上去等于说"这里有个未知的构型",同样是假话。判断走
/// `geometry_could_decide`,与从坐标读立体那一侧用的是同一个 —— 各写一套的话,
/// 同一根键会「读的时候算立体源、写的时候不算」,那正好是凭空造构型的形状。
///
/// 已经有顺反信息的键也不标 —— 那种键写出来的是真实构型,不是未知。
#[must_use]
pub fn unspecified_cis_trans(mol: &MolBuilder) -> Vec<bool> {
    let classes = crate::canon::symmetry_classes(mol);
    (0..mol.num_bonds() as u32)
        .map(|di| geometry_could_decide(mol, di, &classes))
        .collect()
}

fn stereo_candidate(
    mol: &MolBuilder,
    unknown: &[bool],
    di: u32,
    classes: &[u32],
) -> Option<[u32; 2]> {
    if !geometry_could_decide(mol, di, classes) {
        return None;
    }
    let db = *mol.bonds().get(di as usize)?;
    // 文件明说这根键的立体未知(交叉双键)—— 坐标照样画得出一个确定的样子,
    // 照读就等于把"作者说不知道"改写成"作者说是顺式"。
    //
    // 切片短了一截时按"未知"处理:那一档不标,总数当场掉下来;反过来按
    // "已知"兜底的话,少传的信息会变成静默多标出来的立体。
    let unsure = |bi: u32| unknown.get(bi as usize).copied().unwrap_or(true);
    if unsure(di) {
        return None;
    }
    // 参照按存储序挑第一个合格邻居。挑哪个不改变分子 —— 换一个参照只是换一套
    // 坐标系,顺反的值跟着挑中的那个一起变。
    let first = |end: u32, other: u32| {
        reference_neighbours(mol, end, other)
            .into_iter()
            .find(|&(_, bi)| !unsure(bi))
            .map(|(o, _)| o)
    };
    Some([first(db.begin, db.end)?, first(db.end, db.begin)?])
}

/// 走一遍所有双键,把算得出来的顺反写进 `stereo` 与 `stereo_atoms`。
///
/// `geometry` 拿到 `(双键起点, 双键终点, begin 侧参照, end 侧参照)` 四个原子号,
/// 给出顺/反或者"读不出来"。
fn assign_bond_stereo(
    mol: &mut MolBuilder,
    unknown: &[bool],
    geometry: impl Fn(u32, u32, u32, u32) -> Option<BondStereo>,
) -> usize {
    let classes = symmetry_classes(mol);
    let found: Vec<(u32, BondStereo, [u32; 2])> = (0..mol.num_bonds())
        .filter_map(|di| {
            let di = u32::try_from(di).ok()?;
            let [ra, rb] = stereo_candidate(mol, unknown, di, &classes)?;
            let db = *mol.bonds().get(di as usize)?;
            let stereo = geometry(db.begin, db.end, ra, rb)?;
            Some((di, stereo, [ra, rb]))
        })
        .collect();
    let n = found.len();
    for (di, stereo, atoms) in found {
        if let Some(mut b) = mol.bond_mut(di) {
            b.set_stereo(stereo);
            b.set_stereo_atoms(atoms);
        }
    }
    n
}

/// 从**二维坐标**反读双键顺反,写入 `stereo` 与 `stereo_atoms`。返回标注的根数。
///
/// `unknown` 按键下标索引,标出输入里**明说立体未知**的键(molblock 键块第四列
/// 的交叉双键 `3`、波浪单键 `4`)。没有这类信息时传一个全 `false` 的切片。
///
/// # 与 [`perceive_bond_stereo`] 的分工:同一个问题,两种输入
///
/// SMILES 把顺反写在方向键(`/` `\`)上,molblock 把它画在坐标里 —— **图里
/// 没有方向键可读**。两条路只在"参照原子从哪来"上不同:那边是携带方向的那根
/// 键,这边是按存储序挑的第一个合格邻居。什么算立体源(小环、芳香、两端能否
/// 区分取代基)由 `reference_neighbours` / [`stereo_atoms_are_valid`] 与
/// `end_is_stereogenic` 共用,两条路一模一样。
///
/// # 前置条件:先净化
///
/// 要用对称等价类判"这一端的两个取代基分不分得开",而等价类要芳香性定下来
/// 才算得准。刚从文件读出来的分子还没有,那时调这个函数,芳香环外的双键会
/// 被当成普通双键。顺序与 [`assign_chirality_2d`] 一样:
/// **读文件(L1)→ 净化(L2)→ 回来打标记(L1)**。
///
/// # 只管二维,而且**当场认出三维就整个不做**
///
/// 三维坐标投到 xy 平面上照样算得出"同侧/异侧",而那个答案与分子无关 ——
/// 一根真正的反式双键投影下来完全可能落成同侧。这与手性那一侧不对称:三维
/// molblock 的楔形一般是空的,[`assign_chirality_2d`] 自然什么
/// 也标不出来;顺反这一侧不拦的话给出的是**错的答案**,不是空答案。
///
/// 所以任何一个 `z` 不为零就整个不做,返回 0。三维的顺反在扭转角里,是另一
/// 条路。
///
/// # 波浪单键这一档只有单元测试守着
///
/// 交叉双键在语料里有 625 根,判据踩得实。波浪单键**一根都没有** ——
/// 外部实现给这批分子写出的立体码只有 `0`/`1`/`3`/`6`。也就是说这一档只有
/// 单元测试,没有语料级判据。照实记下来,免得"全绿"被读成两档都验过了。
pub fn assign_bond_stereo_2d(mol: &mut MolBuilder, coords: &[[f64; 3]], unknown: &[bool]) -> usize {
    if coords.iter().any(|p| p[2].abs() > FLAT_TOL) {
        return 0;
    }
    let xy = |a: u32| coords.get(a as usize).map(|p| [p[0], p[1]]);
    assign_bond_stereo(mol, unknown, |b, e, ra, rb| {
        cis_trans_from_points(xy(b)?, xy(e)?, xy(ra)?, xy(rb)?)
    })
}

/// 两个参照原子绕双键轴**转到了同一边还是对面** —— 同一边为顺。
///
/// 判不出来(某个参照几乎与双键轴共线)返回 `None`。
///
/// 只吃四个点,不吃分子 —— 与 [`cis_trans_from_points`] 分工相同,那个是二维
/// 投影,这个是二面角。
///
/// # 没有"接近 90° 就不判"的死区
///
/// 二面角正好 90° 的双键在化学上没有顺反可言,可**文件里不会有** ——
/// 双键是平面的。留一条死区反而会在几何本来就无意义的地方与外部实现分歧
/// (它按 90° 一刀切,没有死区)。所以这里也只看符号,另外挡住真正的退化:
/// 参照与轴共线时垂直分量为零,那时"哪一边"根本不成立。
#[must_use]
pub fn cis_trans_from_torsion(
    begin: [f64; 3],
    end: [f64; 3],
    ref_begin: [f64; 3],
    ref_end: [f64; 3],
) -> Option<BondStereo> {
    let sub = |p: [f64; 3], q: [f64; 3]| [p[0] - q[0], p[1] - q[1], p[2] - q[2]];
    let dot = |p: [f64; 3], q: [f64; 3]| p[0] * q[0] + p[1] * q[1] + p[2] * q[2];

    let axis = sub(end, begin);
    let n2 = dot(axis, axis);
    if n2 < FLAT_TOL {
        return None; // 双键两端重合
    }
    // 各自去掉沿轴的分量,剩下的就是"绕轴指向哪边"
    let perp = |v: [f64; 3]| {
        let k = dot(v, axis) / n2;
        [v[0] - k * axis[0], v[1] - k * axis[1], v[2] - k * axis[2]]
    };
    let (pa, pb) = (perp(sub(ref_begin, begin)), perp(sub(ref_end, end)));
    let (la, lb) = (dot(pa, pa).sqrt(), dot(pb, pb).sqrt());
    if la < FLAT_TOL || lb < FLAT_TOL {
        return None; // 参照与轴共线(累积双键那一档就是这样)
    }
    let cos = dot(pa, pb) / (la * lb);
    if cos.abs() < FLAT_TOL {
        return None;
    }
    Some(if cos > 0.0 {
        BondStereo::Cis
    } else {
        BondStereo::Trans
    })
}

/// 从**三维坐标**反读双键顺反,写入 `stereo` 与 `stereo_atoms`。返回标注的根数。
///
/// 与 [`assign_bond_stereo_2d`] 的分工只在几何那一步:那边看两个参照在双键
/// 两侧还是同侧(平面投影),这边看它们绕轴的**二面角**。什么算立体源、
/// 参照原子怎么挑,两边共用 `stereo_candidate` —— 各写一遍的话,同一根键在
/// 二维文件里算立体源、在三维文件里不算,而那种差别不报错。
///
/// 坐标是平的(所有 `z` 都是 0)时一根也不标:那是二维图,投影下来两个参照的
/// 二面角只会是 0° 或 180°,算得出答案但那是二维那条路的答案,该由它来给。
///
/// `unknown` 的含义见 [`assign_bond_stereo_2d`]。
pub fn assign_bond_stereo_3d(mol: &mut MolBuilder, coords: &[[f64; 3]], unknown: &[bool]) -> usize {
    if !coords.iter().any(|p| p[2].abs() > FLAT_TOL) {
        return 0;
    }
    let at = |a: u32| coords.get(a as usize).copied();
    assign_bond_stereo(mol, unknown, |b, e, ra, rb| {
        cis_trans_from_torsion(at(b)?, at(e)?, at(ra)?, at(rb)?)
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

    // **同一个原子上的两根方向键必须反号。**
    //
    // 一个双键端只有两个取代基,它们分处双键两侧 —— 两根方向键都朝外指同一
    // 个方向的话,读串的人只能得到"两个取代基在同一侧",几何上不成立。
    //
    // # 这一条不是多余的,而且它管的是共轭
    //
    // 一根双键的参照原子挑哪个不影响分子,可**挑法一变,方向符号就落到另一根
    // 键上**。共轭链里中间那根单键同时是两根双键的候选参照:两边都挑它,一根
    // 键就够了(先前一直是这个样子);其中一根改挑别的取代基,那个端点上就
    // 冒出**两根**方向键 —— 而 BFS 只按"双键"连边,压根不知道这两根有关系,
    // 于是各自算各自的,同号了也没人拦。
    //
    // 实测:三维语料里补过显式氢的分子改变了规范秩,参照挑法跟着变,当场
    // 撞上这一档 —— `[N-]/C=N/C(C#N)=C(\N)C#N` 那一族的顺反被外部实现整个丢掉
    // (它读到自相矛盾的一对方向,只能放弃)。隐式氢那条路上碰不到,因为那时
    // 秩恰好让两边都挑中共用的那根键。
    // **按双键的端点收,不按参照键的锚点收。** 一根参照键锚在哪一端是它作为
    // 某根双键的参照时的角色;而它同时**碰着**另一根双键的端点,在那里它就是
    // 一根"这一端的方向键",要与那一端的另一根方向键反号。共轭链里冲突正是
    // 这么出来的:共用单键锚在前一根双键那侧,冲突却在后一根双键的端点上。
    let refbonds: BTreeSet<u32> = constraints
        .iter()
        .flat_map(|&(r1, r2, _)| [r1.bond, r2.bond])
        .collect();
    for di in 0..mol.num_bonds() as u32 {
        if !stereo_atoms_are_valid(mol, di) {
            continue;
        }
        let db = mol.bonds()[di as usize];
        for end in [db.begin, db.end] {
            let here: Vec<Ref> = mol
                .neighbors(end)
                .filter(|&(_, bi)| bi != di && refbonds.contains(&bi))
                .map(|(_, bi)| Ref {
                    bond: bi,
                    anchor: end,
                })
                .collect();
            for i in 0..here.len() {
                for j in i + 1..here.len() {
                    let (a, b) = (here[i], here[j]);
                    adj.entry(a.bond).or_default().push((a, b, false));
                    adj.entry(b.bond).or_default().push((b, a, false));
                }
            }
        }
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
        // **先挑非氢的,再按规范秩。** 方向符号落在哪根键上,由参照挑在哪一侧
        // 决定;落在一个显式氢上的话,那个符号就挂在一个**随时会被删掉**的原子上。
        //
        // 实测:补过显式氢的共轭多烯,方向全落在氢上之后,外部实现按默认设置
        // (解析时去显式氢)读回来会翻掉一根双键 —— `C/C=C/C=C\C` 读成
        // `C/C=C/C=C/C`。同一个串让它**保留氢**去读又是对的,所以我方的编码
        // 并不违规,只是**易碎**:它把立体挂在了别人会丢掉的原子上。外部实现
        // 自己写显式氢时把方向放在重原子键上,扛得住去氢。
        //
        // 二维/隐式氢的分子里根本没有氢原子,这一条不改变任何东西 ——
        // 语料 8825 条逐条一致的结果照旧。
        let pick = |end: u32, partner: u32| {
            let key = |other: &u32| {
                let h = u8::from(mol.atoms()[*other as usize].atomic_num == 1);
                (h, priority[*other as usize])
            };
            mol.neighbors(end)
                .map(|(other, _)| other)
                .filter(|&other| other != partner)
                .min_by_key(key)
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

    /// 走产品那条路读一张二维图:文本 → `read_v2000` → 净化 → 反读顺反。
    ///
    /// 自己拼 `MolBuilder` 会绕开读取器与净化,量的就是另一件事了。
    fn ez_from_block(block: &str) -> (usize, String) {
        let got = crate::molblock::read_v2000(block).expect("读 molblock");
        let mut m = got.mol;
        omgkit_chem::pipeline::sanitize(&mut m).expect("净化");
        let n = assign_bond_stereo_2d(&mut m, &got.coords, &got.unknown_stereo);
        (n, crate::canon::canonical_smiles(&m).smiles)
    }

    /// 同一个分子**从 SMILES 那条路**进来的规范串。
    ///
    /// 拿它当靶子,而不是把某个具体的规范串写死在测试里:后者会在规范化算法
    /// 稍一改动时变红,而那时代码并没有错。这里要钉的是"文件那条路与 SMILES
    /// 那条路读出同一个分子"。跨实现那一侧由 `harness/check_molblock_read.py`
    /// 守着。
    fn ez_from_smiles(smi: &str) -> String {
        let mut m = smiles::parse(smi).unwrap_or_else(|e| panic!("{smi}: {}", e.render()));
        omgkit_chem::pipeline::sanitize(&mut m).expect("净化");
        perceive_bond_stereo(&mut m);
        crate::canon::canonical_smiles(&m).smiles
    }

    const TRANS_DIFLUOROETHENE: &str = "\
F/C=C/F
     RDKit          2D

  4  3  0  0  0  0  0  0  0  0999 V2000
   -1.9796   -0.1365    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5994    0.4508    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5994   -0.4508    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9796    0.1365    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  3  4  1  0
M  END
";

    /// 图上画着反式,读出来就得是反式 —— 顺反的**值**读对了,不只是"读到了"。
    #[test]
    fn a_trans_double_bond_is_read_as_trans() {
        let (n, smi) = ez_from_block(TRANS_DIFLUOROETHENE);
        assert_eq!(n, 1, "该标一根");
        assert_eq!(smi, ez_from_smiles("F/C=C/F"));
    }

    /// **交叉双键(键块第四列 `3`)是"作者说不知道"**,不是"作者没写"。
    ///
    /// 坐标照样画得出一个确定的样子 —— 照读就等于替作者把话说死。语料里
    /// 外部实现写出 625 根这样的键,一根都不该被读成构型。
    #[test]
    fn a_crossed_double_bond_is_not_read_as_a_configuration() {
        let crossed = TRANS_DIFLUOROETHENE.replace("  2  3  2  0", "  2  3  2  3");
        assert_ne!(crossed, TRANS_DIFLUOROETHENE, "键块那一行没改到");
        let (n, smi) = ez_from_block(&crossed);
        assert_eq!(n, 0, "不该标");
        assert_eq!(smi, ez_from_smiles("FC=CF"));
    }

    /// **波浪单键(键块第四列 `4`)**同样是"作者说不知道",挨着它的那根双键
    /// 读不出构型。
    ///
    /// 语料里一根波浪键都没有 —— 这一档只有这个测试守着。
    #[test]
    fn a_wavy_single_bond_blocks_its_neighbouring_double_bond() {
        let wavy = TRANS_DIFLUOROETHENE.replace("  1  2  1  0", "  1  2  1  4");
        assert_ne!(wavy, TRANS_DIFLUOROETHENE, "键块那一行没改到");
        let (n, smi) = ez_from_block(&wavy);
        assert_eq!(n, 0, "不该标");
        assert_eq!(smi, ez_from_smiles("FC=CF"));
    }

    /// 三维坐标整个不做 —— 投影到 xy 平面上算出来的"同侧"与分子无关。
    #[test]
    fn three_dimensional_coordinates_are_refused_outright() {
        let lifted = TRANS_DIFLUOROETHENE.replace(
            "    1.9796    0.1365    0.0000 F",
            "    1.9796    0.1365    1.2000 F",
        );
        assert_ne!(lifted, TRANS_DIFLUOROETHENE, "原子块那一行没改到");
        let (n, _) = ez_from_block(&lifted);
        assert_eq!(n, 0, "三维不该走这条路");
    }

    /// 补上显式氢之后,**参照原子不许挑到氢头上**。
    ///
    /// 挑到氢,方向符号就落在一个"别人随手会删掉"的原子上。实测:外部实现按
    /// 默认设置(解析时去显式氢)读回来会翻掉一根双键 —— `C/C=C/C=C\\C` 读成
    /// `C/C=C/C=C/C`。同一个串让它保留氢去读又是对的,所以我方的编码并不违规,
    /// 只是把立体挂在了会被丢掉的原子上。
    #[test]
    fn a_reference_atom_is_never_a_hydrogen_when_a_heavy_one_exists() {
        let mut m = smiles::parse("C/C=C/C=C\\C").expect("解析");
        omgkit_chem::pipeline::sanitize(&mut m).expect("净化");
        perceive_bond_stereo(&mut m);
        let ranks = crate::canon::classed_ranks(&m);
        omgkit_chem::add_explicit_hs(&mut m, &ranks);
        let pr = crate::canon::canonical_ranks(&m);
        let norm = normalized_stereo_refs(&m, &pr).expect("有双键要规范化");
        let mut checked = 0;
        for b in norm.bonds() {
            if !matches!(b.stereo, BondStereo::Cis | BondStereo::Trans) {
                continue;
            }
            for r in b.stereo_atoms {
                assert_ne!(
                    norm.atoms()[r as usize].atomic_num,
                    1,
                    "参照挑到了氢:{:?}",
                    b.stereo_atoms
                );
                checked += 1;
            }
        }
        assert_eq!(checked, 4, "该查两根双键、四个参照");
    }

    /// **同一个双键端上的两根方向键必须反号。**
    ///
    /// 一个双键端只有两个取代基,它们分处双键两侧;两根方向键朝外指同一个方向
    /// 的话,读串的人只能得到"两个取代基在同一侧",几何上不成立。
    ///
    /// 共轭链里才出得来:中间那根单键同时是两根双键的候选参照,一根双键改挑了
    /// 别的取代基,那个端点上就冒出两根方向键 —— 而约束图先前只按双键连边,
    /// 压根不知道这两根有关系。外部实现读到自相矛盾的一对,只能把整根双键的
    /// 立体丢掉。
    #[test]
    fn two_direction_bonds_at_one_end_never_point_the_same_way() {
        for smi in [
            "N#CCc1ccc([N-]/C=N/C(C#N)=C(\\N)C#N)cc1",
            "C/C=C/C=C\\C",
            "F/C=C/C=C/F",
            "CCOC(=O)/C([O-])=C(C#N)\\C=N\\c1ccccc1",
        ] {
            let mut m = smiles::parse(smi).unwrap_or_else(|e| panic!("{smi}: {}", e.render()));
            omgkit_chem::pipeline::sanitize(&mut m).expect("净化");
            perceive_bond_stereo(&mut m);
            let ranks = crate::canon::classed_ranks(&m);
            omgkit_chem::add_explicit_hs(&mut m, &ranks);
            let pr = crate::canon::canonical_ranks(&m);
            let norm = normalized_stereo_refs(&m, &pr).unwrap_or(m);
            let d = directions_for_writing(&norm);
            for a in 0..u32::try_from(norm.num_atoms()).expect("原子数") {
                let out: Vec<(u32, BondDirection)> = norm
                    .neighbors(a)
                    .filter(|&(_, bi)| d.dirs[bi as usize] != BondDirection::None)
                    .filter(|&(_, bi)| norm.bonds()[bi as usize].order != BondOrder::Double)
                    .map(|(_, bi)| {
                        let b = norm.bonds()[bi as usize];
                        let dir = if b.begin == a {
                            d.dirs[bi as usize]
                        } else {
                            d.dirs[bi as usize].flipped()
                        };
                        (bi, dir)
                    })
                    .collect();
                for i in 0..out.len() {
                    for j in i + 1..out.len() {
                        assert_ne!(
                            out[i].1, out[j].1,
                            "{smi}:原子 {a} 上的键 {} 与 {} 朝外同号",
                            out[i].0, out[j].0
                        );
                    }
                }
            }
        }
    }

    /// 一端挂着两个相同的取代基 —— 交换它们是自同构,这根键没有顺反可言。
    ///
    /// **这一档语料级判据够不着**:写出侧的 `informative_directions` 会把这根
    /// 噪声方向再滤掉一次,于是标错了也看不出来。所以钉在这里。
    #[test]
    fn an_end_with_two_identical_substituents_is_not_a_stereo_source() {
        let (n, _) = ez_from_block(
            "\
CC=C(Cl)Cl
     RDKit          2D

  5  4  0  0  0  0  0  0  0  0999 V2000
   -2.0785    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7794    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5196   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8187    0.7500    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    0.5196   -1.5000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  3  4  1  0
  3  5  1  0
M  END
",
        );
        assert_eq!(n, 0, "1,1-二氯丙烯没有顺反");
    }

    /// 累积双键不记顺反。
    ///
    /// **钉的是结论,不是理由。** 真正拦住它的是几何:丙二烯画出来是直线,
    /// 那根双键的另一端正好落在轴上,[`cis_trans_from_points`] 判"读不出来"。
    /// `reference_neighbours` 里排除双键那一条即使拿掉,这个测试照样绿 ——
    /// 见那个函数的文档。
    #[test]
    fn cumulated_double_bonds_are_not_cis_trans() {
        let (n, _) = ez_from_block(
            "\
FC=C=CF
     RDKit          2D

  5  4  0  0  0  0  0  0  0  0999 V2000
   -2.2500    0.7794    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000   -0.5196    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000   -0.5196    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000   -0.5196    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2500    0.7794    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  3  4  2  0
  4  5  1  0
M  END
",
        );
        assert_eq!(n, 0, "丙二烯型的两根双键都不是顺反");
    }

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

    /// 净化过的分子。这条判据要真双键(`order == Double`),而芳香写法的
    /// 键级要靠凯库勒化才定下来 —— 所以这里比同文件别处多跑一步净化。
    fn sanitized(smi: &str) -> MolBuilder {
        let mut m = smiles::parse(smi).unwrap_or_else(|e| panic!("{smi}: {}", e.render()));
        omgkit_chem::sanitize(&mut m).expect("净化失败");
        m
    }

    /// 漏了顺反感知要能被发现 —— 那是会**静默丢几何**的一步。
    ///
    /// 反面同样要守:一条孤零零的 `/` 说明不了任何事,感知跑过也不会留下顺反,
    /// 那种分子不能被报成"漏了"。误报比不报更糟:一条老是错响的断言,下一个人
    /// 就把它关掉了。
    #[test]
    fn missing_stereo_perception_is_detected_without_false_alarms() {
        for smi in ["C/C=C/C", "C/C=C\\C", "F/C(Cl)=C(/Br)C"] {
            let m = sanitized(smi);
            assert!(
                directions_not_perceived(&m),
                "{smi}:两端方向成对却没顺反,该报出来"
            );
            let mut perceived = sanitized(smi);
            perceive_bond_stereo(&mut perceived);
            assert!(
                !directions_not_perceived(&perceived),
                "{smi}:感知跑过了还在报"
            );
        }

        for smi in [
            "CC=CC",       // 压根没有方向键
            "C/C=CC",      // 只有一侧有方向,几何定不下来
            "CCO",         // 连双键都没有
            "C/C=C/C=C/C", // 多根共轭双键,方向成对
        ] {
            let mut m = sanitized(smi);
            perceive_bond_stereo(&mut m);
            assert!(!directions_not_perceived(&m), "{smi}:感知跑过了不该报");
        }

        let lone = sanitized("C/C=CC");
        assert!(
            !directions_not_perceived(&lone),
            "只有一侧写了方向,几何本就定不下来 —— 没感知过也不该报"
        );
    }

    #[test]
    fn 谓词与感知问的是同一个问题() {
        // **每一条都断言"感知标了几根"与"谓词说没说漏"是一致的。**
        //
        // 先前 `directions_not_perceived` 走的是 `raw_cis_trans`,少了
        // `informative_directions` 那道过滤,于是下面第 3~5 条(方向是噪声、
        // 感知有意不标)被报成"漏了调感知"—— 而再调一次也不会变绿。
        // `omgkit-depict` 与 `omgkit-match` 拿这个谓词做 `debug_assert!`,
        // 那三条合法 SMILES 在 debug 下会直接 panic。
        for smi in [
            "F/C=C/F",            // 正常:两侧都有信息
            "F/C=C\\F",           // 同上,另一侧
            "F/C=C(\\F)F",        // 一端两个相同取代基 —— 方向是噪声
            "Cl/C=C(\\Cl)Cl",     // 同上
            "CC(/C=C/C)=C(\\C)C", // 一根真顺反 + 一根冗余方向(公共库常见写法)
            "F/C=CF",             // 只有一侧有方向,说明不了相对位置
            "C/1CCCCC1",          // 根本没有双键
            "CCO",                // 一根方向键都没有
        ] {
            let mut m = crate::smiles::parse(smi).expect("解析");
            omgkit_chem::pipeline::sanitize(&mut m).expect("净化");
            // 感知之**前**:谓词说"漏了"当且仅当感知真的会标出东西
            let will_annotate = {
                let informative = informative_directions(&m);
                (0..m.num_bonds())
                    .filter(|&i| m.bonds()[i].stereo == BondStereo::None)
                    .filter(|&i| would_annotate(&m, i, &informative).is_some())
                    .count()
            };
            assert_eq!(
                directions_not_perceived(&m),
                will_annotate > 0,
                "{smi}:谓词与感知不一致(感知会标 {will_annotate} 根)"
            );
            let n = perceive_bond_stereo(&mut m);
            assert_eq!(n, will_annotate, "{smi}:实际标的根数与预判不符");
            // 感知之**后**:谓词必须闭嘴,否则调用方无论如何都修不绿
            assert!(
                !directions_not_perceived(&m),
                "{smi}:感知跑过了谓词还在报 —— 这种红没有任何修法"
            );
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

    /// **小环里的双键没有顺反**,阈值是八元 —— 而且三条通路必须同时认这条线:
    /// 感知([`perceive_bond_stereo`])、写出([`informative_directions`])、
    /// 匹配([`raw_cis_trans`])。
    ///
    /// 三条各自问一遍是必需的:它们是三个入口,先前只在感知那一路加过滤,
    /// 匹配那一路照旧,于是差分测试上冒出 6 条"只有本实现命中"。
    ///
    /// 四元到十元逐个走一遍,阈值挪一格就红。
    #[test]
    fn 小环里的双键不给顺反() {
        // n 元环,环内一根 C=C,两侧的方向都落在环上
        let ring = |n: usize| format!("C/1=C\\{}1", "C".repeat(n - 2));
        let double_bond = |m: &MolBuilder| {
            m.bonds()
                .iter()
                .position(|b| b.order == BondOrder::Double)
                .expect("有双键") as u32
        };
        for n in 4..=10 {
            let smi = ring(n);
            let m = smiles::parse(&smi).unwrap_or_else(|e| panic!("{smi}: {}", e.render()));
            assert_eq!(m.num_atoms(), n, "{smi}:环大小没写对,这条用例就测错了东西");
            // **这里的 8 是写死的字面量,不许换成 MIN_STEREOGENIC_RING。**
            // 引用被测常量的话,把阈值改成 7 或 9 判据两侧一起动,照样全绿 ——
            // 那是条自证的断言。八元是化学事实(反式环辛烯是最小的可分离
            // 反式环烯),写死才拦得住"顺手挪一格"。
            let stereogenic = n >= 8;
            assert_eq!(
                perceive(&smi).len(),
                usize::from(stereogenic),
                "{n} 元环 {smi}:感知"
            );
            assert_eq!(
                informative_directions(&m).iter().any(|&x| x),
                stereogenic,
                "{n} 元环 {smi}:写出"
            );
            assert_eq!(
                raw_cis_trans(&m, double_bond(&m)).is_some(),
                stereogenic,
                "{n} 元环 {smi}:匹配"
            );
        }

        // **反过来:挂在小环上的环外双键照给。** 少了这一组,"凡是沾着小环的
        // 双键一律不给"也能让上面全绿,而那会把一大批真顺反抹掉。
        let with_raw = |m: &MolBuilder| {
            (0..m.num_bonds() as u32)
                .filter(|&b| raw_cis_trans(m, b).is_some())
                .count()
        };
        for smi in [r"F/C=C1\CCCC(F)C1", r"O=C1CCC/C1=C/F"] {
            let m = smiles::parse(smi).unwrap();
            assert_eq!(perceive(smi).len(), 1, "{smi}:环外双键,顺反是真的");
            assert_eq!(with_raw(&m), 1, "{smi}:匹配那一路也要给");
        }

        // 语料里真实撞上这条规则的那个分子(`large.smi` 第 5707 条,文件第 5731
        // 行):两个稠合
        // 五元环,融合处的 C=C 两侧写着 `\` 与 `/`,合起来要求"反式" —— 五元环
        // 里根本搭不出来。外部实现同样判它没有顺反(实测 RDKit 2025.09.2)。
        let smi = r"CN1CCC\\2=C1/C(=N\\O)/S/C2=N\\c3ccc(cc3)F";
        let got = perceive(smi);
        assert_eq!(got.len(), 2, "{smi}:两根环外 C=N 有顺反,环内那根 C=C 没有");
        let m = smiles::parse(smi).unwrap();
        let ring_db = m
            .bonds()
            .iter()
            .position(|b| {
                b.order == BondOrder::Double
                    && m.atoms()[b.begin as usize].atomic_num == 6
                    && m.atoms()[b.end as usize].atomic_num == 6
            })
            .expect("有那根 C=C") as u32;
        assert!(raw_cis_trans(&m, ring_db).is_none(), "五元环里的 C=C");
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
