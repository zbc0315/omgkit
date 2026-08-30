//! 立体保真:楔形键指派,以及从坐标反读立体化学。
//!
//! # 这一层针对的是"画错",不是"画得难看"
//!
//! Mayfield 在 RDKit UGM 2016 上的原话是,所有布局算法都会产出拥挤甚至
//! **misrepresent(画错)** 的结构。画错里最要命的两种:
//!
//! - 双键在图里是 Z,画出来成了 E
//! - 手性中心画上楔形,反读回来是对映体
//!
//! 两种都不改拓扑、不改原子数、不改键长 —— 任何"数得出来"的判据都发现不了。
//! 本模块把它们变成**可判定**的:从坐标反读一遍,与图里记的比对。
//!
//! # 参照系:手性标记相对**存储序**,隐式氢占槽位 1
//!
//! [`omgkit_io::smiles`] 的模块文档把这条约定连同实测用例写清楚了:标记相对
//! 邻居的**存储顺序**,而解析时隐式氢**不参与置换**,由一个 `degree == 3` 的
//! 特判补偿。
//!
//! 从那份文档给的两个用例可以反推出存储侧的语义:`N[C@H](O)F` 不翻(存 `Ccw`)、
//! `[C@H](N)(O)F` 要翻(存 `Cw`),而两者的书写序只差一次对换 —— 所以隐式氢在
//! 存储序里**占槽位 1**。这与 `omgkit-match` 反应重定基里的
//! `aligned.insert(1, IMPLICIT_H)` 是同一条约定。
//!
//! 记错这个位置不会让任何一步报错,只会让一半的手性中心画成对映体。
//!
//! # 前置条件:顺反要先感知过
//!
//! [`read_bond_stereo`] 读的是键自己记的 [`BondData::stereo`],而**那一步不在
//! `omgkit_chem::pipeline::sanitize` 里** —— 它由 [`omgkit_io::stereo::perceive_bond_stereo`]
//! 单独完成。只跑了净化的分子,每根双键的 `stereo` 都是 `None`,于是顺反这一
//! 档会**静默地什么都不检查**。
//!
//! (Python 绑定的 `Mol.sanitize()` 把两步合在一起,所以从 Python 看不出这个
//! 区别 —— 直接用 Rust 接口时要自己补上那一步。)

use omgkit_core::{BondData, BondFlags, BondOrder, BondStereo, ChiralTag, MolBuilder};

/// 楔形。**定义搬到了 [`omgkit_io::wedge`]** —— 楔形是 molblock 键块第四列的
/// 字段,读文件和画图共用同一套语义,住在 L1 才够得着。这里 `pub use` 过来,
/// 本 crate 的用法一字不改。
pub use omgkit_io::wedge::Wedge;

/// 从坐标与楔形反读一个中心的手性。判不出来返回 `None`。
///
/// **实现在 [`omgkit_io::wedge::chirality_from_wedges`]**,这里只把本 crate 的
/// [`Point2`] 翻成它要的 `[f64; 3]`。两处各写一份的话,
/// "读文件"与"画图"迟早对同一张图给出不同的构型。
#[must_use]
pub fn read_chirality(
    mol: &MolBuilder,
    coords: &[Point2],
    wedges: &[Wedge],
    a: u32,
) -> Option<ChiralTag> {
    let xyz: Vec<[f64; 3]> = coords.iter().map(|p| [p.x, p.y, 0.0]).collect();
    omgkit_io::wedge::chirality_from_wedges(mol, &xyz, wedges, a)
}

use crate::geom::Point2;

/// 楔形指派的结果。
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Wedges {
    /// 逐键,下标与 [`MolBuilder`] 的键下标一致
    pub bonds: Vec<Wedge>,
    /// **没能画出构型的立体中心**。如实报出来,不假装画好了。
    pub unwedged: Vec<u32>,
}

/// 给所有四面体立体中心指派楔形键。
///
/// # 做法
///
/// 每个中心挑一根键,试 [`Wedge::Up`] 与 [`Wedge::Down`],取**反读回来等于
/// 图里记的那个标记**的那一个。这是"构造即正确",简单且不会错。
///
/// 代价是这样一来 `assign` 与 [`read_chirality`] 就有了共谋:拿"反读一致"去
/// 检验它是**空过的**。真正的判据是**区分力**:一个分子和它的对映体必须得到
/// 相反的楔形。反读函数若与几何无关,两者会得到同一个楔形 —— 那一条测得出来。
/// 绝对意义上对不对由外部判官定,见 `harness/check_wedge_readback.py`。
///
/// # 候选键的两条忌讳互相冲突,所以两种序都跑
///
/// 楔形该躲开**环键**(IUPAC 的图示建议:立体键应当画向取代基;环键的两个原子
/// 在读者眼里都躺在环平面里),也该躲开**与另一个立体中心共用的键**(读者会问
/// 这个楔形在说谁)。
///
/// 两条常常冲不到一起,但**相邻的两个中心只有一根共用的非环键时就冲了**:谁先
/// 拿走谁受益,另一个被饿死。实测把"环键最差"提到最前面,某个分子里修好了 1 个
/// 却让另一个中心**根本画不出来** —— 净 −3 环上楔形、+1 `unwedged`。
///
/// **哪一条该优先没有普遍答案**,所以不猜:两种序各指派一遍,按结果挑 ——
/// 先比画不出来的中心数(少的赢,信息不能丢),再比落在环键上的楔形数。
/// 两个都跑过一遍才挑,所以结果**不会比任何一种单独跑更差**。
#[must_use]
pub fn assign_wedges(mol: &MolBuilder, coords: &[Point2], ranks: &[u32]) -> Wedges {
    let ring_first = assign_with(mol, coords, ranks, Taboo::RingWorst);
    let shared_first = assign_with(mol, coords, ranks, Taboo::SharedWorst);
    // 落在环键上的楔形有几个 —— 越少越好
    let on_ring = |w: &Wedges| {
        w.bonds
            .iter()
            .enumerate()
            .filter(|(bi, x)| {
                x.narrow().is_some() && mol.bonds()[*bi].flags.contains(BondFlags::IN_RING)
            })
            .count()
    };
    let key = |w: &Wedges| (w.unwedged.len(), on_ring(w));
    // 打平取 `RingWorst` —— 要一个定死的方向,不能让它随机
    if key(&shared_first) < key(&ring_first) {
        shared_first
    } else {
        ring_first
    }
}

/// 候选键的忌讳排序里,哪一条排在最前。见 [`assign_wedges`]。
#[derive(Clone, Copy, PartialEq, Eq)]
enum Taboo {
    /// 环键最忌讳(IUPAC 的口径)
    RingWorst,
    /// 与另一个立体中心共用的键最忌讳(先前的口径)
    SharedWorst,
}

fn assign_with(mol: &MolBuilder, coords: &[Point2], ranks: &[u32], taboo: Taboo) -> Wedges {
    let mut out = Wedges {
        bonds: vec![Wedge::None; mol.num_bonds()],
        unwedged: Vec::new(),
    };
    let mut used: Vec<bool> = vec![false; mol.num_bonds()];
    // 已经指派好的中心,连同它当时想要的构型 —— 全部指派完要回头复核一遍
    let mut picked: Vec<(u32, ChiralTag)> = Vec::new();

    // **只画真手性中心。** `C[C@H](C)O` 这类标记两个取代基一模一样,谈不上
    // 构型;规范 SMILES 会把它抹掉,而这里若照画不误,画出来的楔形方向取决于
    // 邻居的存储顺序 —— 同一个分子换个写法,楔形就换一边。
    //
    // `genuine_tetrahedral` 拿不准时返回 `true`(保留),所以这一道只会滤掉
    // 确定没意义的,不会把真中心滤没。
    let genuine = omgkit_io::stereo::genuine_tetrahedral(mol);
    let mut centres: Vec<u32> = (0..mol.num_atoms())
        .map(|i| u32::try_from(i).expect("原子数超出 u32"))
        .filter(|a| {
            genuine[*a as usize]
                && matches!(
                    mol.atoms()[*a as usize].chiral_tag,
                    ChiralTag::Cw | ChiralTag::Ccw
                )
        })
        .collect();
    // 按规范秩处理 —— 谁先挑走一根键会影响后面的选择
    centres.sort_by_key(|a| (ranks[*a as usize], *a));

    // **配位几何画不出来,那就说出来。**
    //
    // `@SP` / `@TB` / `@OH` 的构型没法用一根楔形表达(平面四方形在纸面上根本
    // 不是"出平面/入平面"的事,三角双锥与八面体要两根以上),这一版不画。
    // 但先前它们**连报都不报**:上面那道 filter 只收 `Cw | Ccw`,于是这些中心
    // 既不指派楔形,也不进 `unwedged`,`is_clean()` 照样为真 ——
    // `[Pt@SP1](Cl)(Cl)(N)N` 画出来读回去就是 `[Pt](Cl)(Cl)(N)N`,构型整个
    // 消失而诊断全绿。README 说"画不好会说出来",这里没说。
    out.unwedged.extend(
        (0..u32::try_from(mol.num_atoms()).expect("原子数超出 u32")).filter(|a| {
            let tag = mol.atoms()[*a as usize].chiral_tag;
            tag != ChiralTag::Unspecified && !tag.is_tetrahedral()
        }),
    );

    for a in centres {
        let want = mol.atoms()[a as usize].chiral_tag;
        let mut done = false;
        // **一根不成就换下一根。** 先前只挑最合适的那一根试,Up、Down 都反读
        // 不对就放弃报 `unwedged` —— 而换一根往往就成了。实测:抗坏血酸两个
        // 立体中心里的那个侧链碳因此没画出来,图上看不出任何异常。
        for bond in candidate_bonds(mol, a, &used, ranks, taboo) {
            for w in [Wedge::Up { narrow: a }, Wedge::Down { narrow: a }] {
                out.bonds[bond as usize] = w;
                if read_chirality(mol, coords, &out.bonds, a) == Some(want) {
                    used[bond as usize] = true;
                    done = true;
                    break;
                }
            }
            if done {
                break;
            }
            out.bonds[bond as usize] = Wedge::None; // 复原,再试下一根
        }
        if !done {
            out.unwedged.push(a);
        } else {
            picked.push((a, want));
        }
    }

    // **回头复核。** 一个中心读得对不对,取决于它周围**所有**楔形 —— 而后面的
    // 中心可能会占走与它共用的那根键,把它读成别的构型甚至读不出来。当时是对的,
    // 全部指派完就未必了。
    //
    // 撤掉一根楔形只会让别处的 z 更少,不会凭空造出新的错,所以这个循环一定收敛。
    while let Some(k) = picked
        .iter()
        .position(|(a, want)| read_chirality(mol, coords, &out.bonds, *a) != Some(*want))
    {
        let (a, _) = picked.remove(k);
        for (_, bi) in mol.neighbors(a) {
            if out.bonds[bi as usize].narrow() == Some(a) {
                out.bonds[bi as usize] = Wedge::None;
            }
        }
        out.unwedged.push(a);
    }
    out.unwedged.sort_unstable();
    out
}

/// 立体中心 `a` 可以打楔形的键,**按优先级从好到差排好**。
///
/// 两条忌讳,哪条排前面由 `taboo` 定 —— 它们会冲突,见 [`assign_wedges`]:
///
/// - **环上的键**:IUPAC 的图示建议说立体键该画向取代基。环键的两个原子在读者
///   眼里都躺在环平面里,声明其中一个出平面与读者正在用的环几何自相矛盾。
/// - **对端也是立体中心的键**:读者会问这个楔形在说谁的构型。几何上它是明确的
///   ([`Wedge`] 带着窄端),但能躲就躲。
///
/// 后面两档不变:**对端是端基原子**(打在端基上最清楚)> 规范秩小。
///
/// 返回的是整个候选序列而不是最好的那一根:最好的那根未必读得回正确的构型,
/// 那时要接着试下一根。平局一律按规范秩打破,拿存储下标打破会引入写法依赖。
fn candidate_bonds(
    mol: &MolBuilder,
    a: u32,
    used: &[bool],
    ranks: &[u32],
    taboo: Taboo,
) -> Vec<u32> {
    let mut cands: Vec<(u8, u8, u8, u32, u32, u32)> = mol
        .neighbors(a)
        .filter(|(_, bi)| !used[*bi as usize])
        // **楔形只画得到单键上。** 双键、三键、芳香键在 `render` 里走的是另一条
        // 分支,那里根本不看楔形 —— 指派到那种键上,楔形就无声无息地没了,而
        // `unwedged` 还是空的,诊断全绿。四配位的 P(V)、S(VI) 中心会碰到。
        .filter(|(_, bi)| mol.bonds()[*bi as usize].order == BondOrder::Single)
        .map(|(n, bi)| {
            let in_ring = mol.bonds()[bi as usize].flags.contains(BondFlags::IN_RING);
            let other_is_centre = matches!(
                mol.atoms()[n as usize].chiral_tag,
                ChiralTag::Cw | ChiralTag::Ccw
            );
            let (first, second) = match taboo {
                Taboo::RingWorst => (in_ring, other_is_centre),
                Taboo::SharedWorst => (other_is_centre, in_ring),
            };
            (
                u8::from(first),
                u8::from(second),
                u8::from(mol.degree(n) > 1),
                ranks[n as usize],
                n,
                bi,
            )
        })
        .collect();
    cands.sort_unstable();
    cands.into_iter().map(|c| c.5).collect()
}

/// 从坐标反读一根双键的顺反。判不出来返回 `None`。
///
/// 用的是键自己记的**参照原子**([`BondData::stereo_atoms`]),不涉及 CIP ——
/// 与 [`BondStereo::Cis`] / [`BondStereo::Trans`] 的定义一致。
///
/// "同侧算顺"这条符号约定在
/// [`omgkit_io::stereo::cis_trans_from_points`],这里只把点递过去 ——
/// 读 molblock 那条路问的是同一个几何问题,约定各写一份的话两条路会给出
/// 一对相反的顺反,而且各自自洽。
#[must_use]
pub fn read_bond_stereo(mol: &MolBuilder, coords: &[Point2], bond: u32) -> Option<BondStereo> {
    let b = &mol.bonds()[bond as usize];
    let [ra, rb] = b.stereo_atoms;
    if ra == BondData::NO_STEREO_ATOM || rb == BondData::NO_STEREO_ATOM {
        return None;
    }
    let p = |a: u32| {
        let q = coords[a as usize];
        [q.x, q.y]
    };
    omgkit_io::stereo::cis_trans_from_points(p(b.begin), p(b.end), p(ra), p(rb))
}

/// 把几何画反了的双键掰回来:把一侧的子树沿双键轴整体镜像。
///
/// 返回被掰过的键。镜像是**等距变换**,键长键角一点不变 —— 这是选它而不是
/// 挪原子的理由。
///
/// # 为什么要在消冲突**之前**做
///
/// 消冲突靠翻转可旋转键,而双键旁边的单键一翻,那一侧的参照原子就换了边 ——
/// 顺反跟着反。所以顺序必须是:先把顺反摆对,再在**不许破坏它**的前提下消冲突
/// (见 `refine` 里的立体守卫)。反过来做的话,冲突消完顺反又被弄反了。
pub(crate) fn fix_cis_trans(mol: &MolBuilder, coords: &mut [Point2], ranks: &[u32]) -> Vec<u32> {
    let mut fixed = Vec::new();
    // **次序与"镜像哪一侧"都不能跟着存储序走。**
    //
    // 掰一根顺反键的做法是把一侧的子树整体镜像。多根顺反键时,**先掰哪根、
    // 掰哪一侧,决定最终的几何** —— 两次镜像不对易。先前一是按键的存储下标
    // 遍历,二是一律镜像 `b.end` 那一侧,而 `begin`/`end` 只是书写痕迹。
    // 于是同一个分子换种写法就摆成另一个样子。
    //
    // 实测:全量语料 78 个写法相关的分子里,**40 个的首次分岔就在这一步**,
    // 是最大的一块。
    //
    // 改成两条:按两端的规范秩排着掰;每根键镜像**最小规范秩更大**的那一侧 ——
    // 让"更靠前"的那半边钉住不动。
    //
    // **两条的分量不一样,语料级变异量过:**
    //
    // | 变异 | 写法无关违例 |
    // |---|---:|
    // | 去掉排序(退回按存储下标掰) | **77 —— 没影响** |
    // | 一律镜像 `end` 那一侧 | **155 —— 全部效果在这** |
    //
    // 所以 155 → 77 全部来自"挑哪一侧",排序那一条在本语料上**一次都没触发**。
    // 留着不是因为量到了收益,是因为按存储序遍历本身就是头号契约的隐患 ——
    // 这一点照实记,不假装它有贡献。
    let mut todo: Vec<u32> = (0..u32::try_from(mol.num_bonds()).expect("键数超出 u32"))
        .filter(|bi| {
            matches!(
                mol.bonds()[*bi as usize].stereo,
                BondStereo::Cis | BondStereo::Trans
            )
        })
        .collect();
    todo.sort_by_key(|bi| {
        let b = &mol.bonds()[*bi as usize];
        let (x, y) = (ranks[b.begin as usize], ranks[b.end as usize]);
        (x.min(y), x.max(y), *bi)
    });
    for bi in todo {
        let want = mol.bonds()[bi as usize].stereo;
        if read_bond_stereo(mol, coords, bi) == Some(want) {
            continue;
        }
        let b = &mol.bonds()[bi as usize];
        // 两侧都算出来,镜像最小规范秩更大的那一侧
        let low = |s: &Vec<u32>| {
            s.iter()
                .map(|a| ranks[*a as usize])
                .min()
                .unwrap_or(u32::MAX)
        };
        let side = match (subtree(mol, bi, b.begin), subtree(mol, bi, b.end)) {
            (Some(x), Some(y)) => {
                if low(&x) > low(&y) {
                    x
                } else {
                    y
                }
            }
            (Some(x), None) => x,
            (None, Some(y)) => y,
            (None, None) => continue,
        };
        let (pa, pb) = (coords[b.begin as usize], coords[b.end as usize]);
        for a in &side {
            coords[*a as usize] = coords[*a as usize].mirrored(pa, pb - pa);
        }
        if read_bond_stereo(mol, coords, bi) == Some(want) {
            fixed.push(bi);
        } else {
            // 掰不过来(多半是这根键落在环上,一侧就是另一侧)—— 原样退回,
            // 让判据把它报出来,不留一个"动过手脚但仍然不对"的中间态
            for a in &side {
                coords[*a as usize] = coords[*a as usize].mirrored(pa, pb - pa);
            }
        }
    }
    fixed
}

/// 画出来的几何与记录的顺反**不符**的那些双键。
///
/// 照**最终坐标**量,所以掰不动的(环上的键,两侧是同一片原子)与掰对了又被消
/// 冲突挪反的都在里面。记着顺反、而几何读不出确定值的也算 —— 那同样是"读图的人
/// 拿不到作者写的那句话"。
pub(crate) fn stereo_mismatches(mol: &MolBuilder, coords: &[Point2]) -> Vec<u32> {
    (0..u32::try_from(mol.num_bonds()).expect("键数超出 u32"))
        .filter(|&bi| {
            let want = mol.bonds()[bi as usize].stereo;
            matches!(want, BondStereo::Cis | BondStereo::Trans)
                && read_bond_stereo(mol, coords, bi) != Some(want)
        })
        .collect()
}

/// 断开键 `bond` 之后,`start` 那一侧的原子。绕回去(环上的键)返回 `None`。
fn subtree(mol: &MolBuilder, bond: u32, start: u32) -> Option<Vec<u32>> {
    let b = &mol.bonds()[bond as usize];
    let blocked = if start == b.end { b.begin } else { b.end };
    let mut seen = std::collections::BTreeSet::from([start]);
    let mut stack = vec![start];
    while let Some(a) = stack.pop() {
        for (n, bi) in mol.neighbors(a) {
            if bi == bond {
                continue;
            }
            if n == blocked {
                return None;
            }
            if seen.insert(n) {
                stack.push(n);
            }
        }
    }
    let mut out: Vec<u32> = seen.into_iter().collect();
    out.sort_unstable();
    Some(out)
}

/// 所有记了顺反的双键,几何是不是都与记录一致。
///
/// 消冲突每翻一次键都要过这一关 —— 翻转会把参照原子换到另一侧。
#[must_use]
pub fn cis_trans_intact(mol: &MolBuilder, coords: &[Point2]) -> bool {
    (0..u32::try_from(mol.num_bonds()).expect("键数超出 u32")).all(|bi| {
        let want = mol.bonds()[bi as usize].stereo;
        !matches!(want, BondStereo::Cis | BondStereo::Trans)
            || read_bond_stereo(mol, coords, bi) == Some(want)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{generate, style::Style};

    fn prep(smi: &str) -> MolBuilder {
        let mut m = omgkit_io::smiles::parse(smi).unwrap();
        omgkit_chem::pipeline::sanitize(&mut m).unwrap();
        m
    }

    #[test]
    fn a_wedge_pointing_the_other_way_reads_as_the_opposite_configuration() {
        // 楔形是**有方向**的:窄端那头在纸面上,宽端那头翘起来。同一根实楔形,
        // 窄端在中心自己这头读作"邻居在上面",窄端在邻居那头就该读作"邻居在
        // 下面" —— 两种读法必须给出相反的构型。
        //
        // 不看窄端的话,相邻的两个立体中心里后画的那个会读成对映体。指派本身
        // 现在会尽量躲开这种键,所以这条只能在这一层单独测 —— 走完整流水线
        // 是碰不到它的。
        let m = prep("N[C@@H](C)O");
        let d = generate(&m, &Style::ACS_1996);
        let a = 1u32;
        assert!(
            matches!(
                m.atoms()[a as usize].chiral_tag,
                ChiralTag::Cw | ChiralTag::Ccw
            ),
            "原子 {a} 该是立体中心"
        );
        let (n, bi) = m.neighbors(a).next().expect("立体中心总有邻居");

        let mut ws = vec![Wedge::None; m.num_bonds()];
        ws[bi as usize] = Wedge::Up { narrow: a };
        let from_centre = read_chirality(&m, &d.coords, &ws, a);
        ws[bi as usize] = Wedge::Up { narrow: n };
        let from_other = read_chirality(&m, &d.coords, &ws, a);

        assert!(from_centre.is_some(), "窄端在中心这头该读得出构型");
        assert!(from_other.is_some(), "窄端在邻居那头也该读得出构型");
        assert_ne!(
            from_centre, from_other,
            "同一根实楔形,窄端换到另一头,读出来的构型却没变"
        );
    }

    #[test]
    fn a_centre_left_wedged_still_reads_back_after_all_the_others_are_placed() {
        // 一个中心读得对不对,取决于它周围**所有**楔形。后面的中心可能占走与它
        // 共用的那根键,把它读成别的构型、甚至读不出来 —— 指派当时是对的,全部
        // 指派完就未必了。
        //
        // 于是会出现:图上那个中心画着楔形,读出来却是另一个构型,而 `unwedged`
        // 是空的。这比不画还糟 —— 不画只是缺信息,画错是给了假信息。
        for smi in [
            "C[C@@H]1[C@H]2[C@H]3C[C@@H](O1)O[C@@H]2OC=C3C(=O)OC",
            "OC[C@H](O)[C@H]1OC(=O)C(O)=C1O",
            "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
            "CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C",
        ] {
            for style in &Style::ALL {
                let m = prep(smi);
                let d = generate(&m, style);
                // `coords`/`wedges` 的下标相对**被画的那个分子** —— 为画出构型补
                // 出来的氢也在里面。`scene` 自己会补,所以它拿的是原分子;这里影子化
                // 之后,下面按下标索引才不会错位甚至越界。
                let m = d.drawn(&m);
                let mut checked = 0;
                for (i, a) in m.atoms().iter().enumerate() {
                    let at = u32::try_from(i).expect("原子数超出 u32");
                    if !matches!(a.chiral_tag, ChiralTag::Cw | ChiralTag::Ccw)
                        || d.unwedged.contains(&at)
                    {
                        continue;
                    }
                    assert_eq!(
                        read_chirality(&m, &d.coords, &d.wedges, at),
                        Some(a.chiral_tag),
                        "[{}] {smi}:中心 {at} 画出来了,反读却不是它该有的构型",
                        style.name
                    );
                    checked += 1;
                }
                assert!(
                    checked >= 2,
                    "[{}] {smi}:只查到 {checked} 个中心",
                    style.name
                );
            }
        }
    }

    #[test]
    fn a_centre_with_four_drawn_neighbours_reads_even_when_the_wedges_cancel() {
        // 四个邻居都画出来的中心不需要摆隐式氢,两个方向相反的楔形照样定得下
        // 构型。把"z 之和为零就读不出"这条无差别地用上去,这种中心会被误判成
        // 读不出 —— 而它的几何一点不含糊。
        let m = prep("O[C@](N)(F)Cl");
        let d = generate(&m, &Style::ACS_1996);
        let a = 1u32;
        let bs: Vec<u32> = m.neighbors(a).map(|(_, bi)| bi).collect();
        assert_eq!(bs.len(), 4, "这个中心该有四个画得出来的邻居");

        let mut ws = vec![Wedge::None; m.num_bonds()];
        ws[bs[0] as usize] = Wedge::Up { narrow: a };
        let one = read_chirality(&m, &d.coords, &ws, a);
        assert!(one.is_some(), "一个楔形就该读得出来");

        ws[bs[1] as usize] = Wedge::Down { narrow: a };
        assert!(
            read_chirality(&m, &d.coords, &ws, a).is_some(),
            "一实一虚两个楔形,z 之和为零,却被判成读不出来"
        );
    }

    #[test]
    fn every_stereocentre_is_either_drawn_or_reported() {
        // **不画出来还不报,是最坏的一种结果**:读者拿到的是一个构型未定的
        // 分子,而线条本身看着一点毛病没有,退化、冲突、交叉全是 0。
        //
        // 这一条数的是:立体中心的个数 = 画了楔形的键数 + `unwedged` 里的个数。
        // 两边对不上就说明有中心被悄悄漏掉了。
        for smi in [
            "OC[C@H](O)[C@H]1OC(=O)C(O)=C1O", // 抗坏血酸:2 个,一个在环上
            "C[C@H](N)C(=O)O",                // 丙氨酸:1 个
            "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O", // 葡萄糖:5 个,全在环上
            "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O", // 青霉素 G:3 个
            "CN1CCC[C@H]1c1cccnc1",           // 尼古丁:1 个,在环上
        ] {
            for style in &Style::ALL {
                let m = prep(smi);
                let d = generate(&m, style);
                // `coords`/`wedges` 的下标相对**被画的那个分子** —— 为画出构型补
                // 出来的氢也在里面。`scene` 自己会补,所以它拿的是原分子;这里影子化
                // 之后,下面按下标索引才不会错位甚至越界。
                let m = d.drawn(&m);
                let centres = m
                    .atoms()
                    .iter()
                    .filter(|a| matches!(a.chiral_tag, ChiralTag::Cw | ChiralTag::Ccw))
                    .count();
                let drawn = d.wedges.iter().filter(|w| **w != Wedge::None).count();
                assert_eq!(
                    centres,
                    drawn + d.unwedged.len(),
                    "[{}] {smi}:{centres} 个立体中心,画了 {drawn} 个,报了 {} 个画不了",
                    style.name,
                    d.unwedged.len()
                );
                // 这几个分子的每个中心都还有键可打 —— 一个都不该落下。
                // 上面那条只保证"漏了会报",这条保证"根本不漏"。
                assert!(
                    d.unwedged.is_empty(),
                    "[{}] {smi}:立体中心 {:?} 没画出构型",
                    style.name,
                    d.unwedged
                );
            }
        }
    }

    #[test]
    fn the_reference_tetrahedron_pins_the_sign() {
        // 有向体积的符号对应哪个手性,靠推是推不放心的 —— 这里用一个**手算得
        // 出朝向**的构型钉死它。
        //
        // 0 号配体在 +z(朝观察者),其余三个在下方平面上依次位于
        // 右 / 左上 / 左下。从 0 号看向中心时,1→2→3 是**逆时针**,按 SMILES
        // 的约定那就是 `@`,即 Ccw。
        let pts = [
            [0.0, 0.0, 1.0],
            [1.0, 0.0, -0.33],
            [-0.5, 0.866, -0.33],
            [-0.5, -0.866, -0.33],
        ];
        let d = |i: usize, j: usize| pts[i][j] - pts[0][j];
        let det = d(1, 0) * (d(2, 1) * d(3, 2) - d(2, 2) * d(3, 1))
            - d(1, 1) * (d(2, 0) * d(3, 2) - d(2, 2) * d(3, 0))
            + d(1, 2) * (d(2, 0) * d(3, 1) - d(2, 1) * d(3, 0));
        assert!(det < 0.0, "参照构型的有向体积应当为负,实得 {det}");
    }

    /// 手搭一个中心:三个邻居按给定角度摆,第一根键带一个实楔形。
    fn cramped(angles_deg: [f64; 3]) -> (MolBuilder, Vec<Point2>, Vec<Wedge>) {
        let mut m = MolBuilder::new();
        let c = m.add_atom(6);
        let mut coords = vec![Point2::ORIGIN];
        let mut wedges = Vec::new();
        for (i, a) in angles_deg.iter().enumerate() {
            // 三个邻居各不相同,免得被 `genuine_tetrahedral` 当成假中心
            let n = m.add_atom(if i == 0 {
                8
            } else if i == 1 {
                7
            } else {
                9
            });
            m.add_bond(c, n, BondOrder::Single).unwrap();
            let r = a.to_radians();
            coords.push(Point2::new(r.cos(), r.sin()));
            wedges.push(if i == 0 {
                Wedge::Up { narrow: c }
            } else {
                Wedge::None
            });
        }
        if let Some(at) = m.atom_mut(c) {
            at.num_implicit_hs = 1;
            at.chiral_tag = ChiralTag::Cw;
        }
        (m, coords, wedges)
    }

    #[test]
    fn a_cramped_stereocentre_is_not_read_as_if_the_h_sat_on_the_centre() {
        // 三个邻居全挤在中心的一侧时,**隐式氢不可能投影到中心上** —— 四面体
        // 的四个键方向之和为零,所以它们的 2D 投影之和也为零,三个都在半平面
        // 里,第四个必然在对面的空扇区。
        //
        // 判据的形式:同一个几何,把氢摆在中心 vs 摆进空扇区,**读出来必须不同**。
        // 不同才说明这个位置是**吃劲**的 —— 先前摆在中心,全量语料上 21 个中心
        // 因此画成了对映体(外部判官 `harness/check_wedge_readback.py` 量的)。
        //
        // **哪一个才对不由这条判据说了算**,它只钉住"位置吃劲"。对错由外部判官
        // 定:改成摆进空扇区之后,那 21 处不一致全部消失(459 → 480 一致,
        // 0 处不一致)。
        // **角度是搜出来的,不是随手写的。** 头一版写了 `[0, 60, 120]`(空隙
        // 240°,确实挤在一侧),两个模型读出来**一样** —— 那时判据是空过的。
        // 穷举 10° 栅格上全部"挤在一侧且三重积够大"的组合,3484 组里能区分的
        // 才有效;这里取空隙 240° 的那一组,两边 det 是 −0.366 / +2.000。
        let angles = [0.0, 30.0, 240.0];
        let (m, coords, wedges) = cramped(angles);
        let got = read_chirality(&m, &coords, &wedges, 0).expect("这个几何读得出来");

        // 把氢摆回中心,手算同一个行列式
        let z = 1.0; // 只有第一根键带楔形,窄端在中心
        let p: Vec<[f64; 3]> = (1..4)
            .map(|i| {
                let q = coords[i];
                [q.x, q.y, if i == 1 { z } else { 0.0 }]
            })
            .collect();
        // 参照序:氢占槽位 1
        let old = [p[0], [0.0, 0.0, -z], p[1], p[2]];
        let d = |i: usize, j: usize| old[i][j] - old[0][j];
        let det = d(1, 0) * (d(2, 1) * d(3, 2) - d(2, 2) * d(3, 1))
            - d(1, 1) * (d(2, 0) * d(3, 2) - d(2, 2) * d(3, 0))
            + d(1, 2) * (d(2, 0) * d(3, 1) - d(2, 1) * d(3, 0));
        let old_tag = if det < 0.0 {
            ChiralTag::Ccw
        } else {
            ChiralTag::Cw
        };
        assert_ne!(
            got, old_tag,
            "把氢摆在中心与摆进空扇区读出了同一个构型 —— 这个几何没有区分力,\
             判据是空过的,换一组角度"
        );
    }

    #[test]
    fn two_nearly_collinear_bonds_cannot_pin_a_configuration() {
        // 三根键里有两根几乎**共线**时,它们张不出面积,楔形再怎么画也定不出
        // 手性 —— 与氢摆在哪无关。这时必须**如实返回 `None`**,不许猜。
        //
        // 实测踩到过:一个笼状分子的中心,两根键相差 184.4°,三重积只有
        // 0.0905。先前照读不误,而外部实现(RDKit,同一个口径、同一个常数
        // `ZERO_VOLUME_TOL = 0.1`)读出来是"未指定" —— 我们说画出来了,别人
        // 读不到,那就是在骗人。
        //
        // 拒掉之后 `assign_wedges` 会接着试下一根候选键,实测那个中心因此**换了
        // 一根键**,外部实现读出了正确的构型 —— 不是被消音。
        let (m, coords, wedges) = cramped([90.0, 0.1, 180.0]); // 后两根差 179.9°
        assert_eq!(
            read_chirality(&m, &coords, &wedges, 0),
            None,
            "两根键几乎共线,这张图定不出手性,该返回 None"
        );

        // 同样的三根键张开一点就读得出来 —— 证明上面那条不是因为别的原因返回 None
        let (m2, c2, w2) = cramped([90.0, 20.0, 180.0]);
        assert!(
            read_chirality(&m2, &c2, &w2, 0).is_some(),
            "张开之后应当读得出来"
        );
    }

    #[test]
    fn a_lone_pair_counts_as_the_fourth_ligand() {
        // 亚砜的硫是三配位 + 一对**孤对电子**,构型照样确定,画法也和碳上的
        // 隐式氢一样 —— 楔形打在三根键之一,孤对在它的反面。
        //
        // 先前 `read_chirality` 只认 `(4 邻居, 0 氢)` 与 `(3, 1)`,`(3, 0)` 一律
        // 返回 `None`,于是这些中心全进 `unwedged`:全量语料 18 个画不出构型的
        // 中心里 **14 个是这一档**,补上之后只剩 4 个。
        //
        // **孤对占哪个槽位由外部判官定,不由这条判据定。** 实测:挪到槽位 3 是
        // 偶置换,读出来一模一样(所以那个变异是空的);挪到槽位 0 是奇置换,
        // 14 个中心**全部翻成对映体**,判官当场红。
        for smi in ["C[S@@](=O)CC", "C[S@](=O)CC"] {
            let m = prep(smi);
            let ranks = omgkit_io::canon::canonical_ranks(&m);
            let d = generate(&m, &Style::ACS_1996);
            let w = assign_wedges(&m, &d.coords, &ranks);
            assert!(
                w.unwedged.is_empty(),
                "{smi}:亚砜的硫没画出构型 {:?}",
                w.unwedged
            );
        }

        // **区分力**:一个分子和它的对映体必须得到相反的楔形。反读函数若与
        // 几何无关,两者会拿到同一个 —— 那时上面那条照样绿。
        let wedge_of = |smi: &str| {
            let m = prep(smi);
            let ranks = omgkit_io::canon::canonical_ranks(&m);
            let d = generate(&m, &Style::ACS_1996);
            let w = assign_wedges(&m, &d.coords, &ranks);
            w.bonds
                .iter()
                .find_map(|x| match x {
                    Wedge::Up { .. } => Some(true),
                    Wedge::Down { .. } => Some(false),
                    Wedge::None => None,
                })
                .expect("该有一个楔形")
        };
        assert_ne!(
            wedge_of("C[S@@](=O)CC"),
            wedge_of("C[S@](=O)CC"),
            "亚砜的一对对映体拿到了同一个楔形 —— 反读没有区分力"
        );
    }

    #[test]
    fn a_three_coordinate_nitrogen_is_not_given_a_wedge() {
        // 三配位氮的孤对**翻转极快**(氨的翻转垒只有 24 kJ/mol),常温下两个
        // 构型互变 —— 画出楔形是在断言一个不存在的构型。
        //
        // 少数被环卡住的(氮杂环丙烷、Tröger 碱)确实稳定,但那要看环张力,
        // 不是看元素。这里宁可**漏报**也不乱画:漏的如实进 `unwedged`,下游
        // 拿到的是"未指定",不会当真。
        //
        // **三个取代基必须各不相同。** 头一版写的 `C[N@](C)CC` 有两个甲基,
        // `genuine_tetrahedral` 先把它当假中心滤掉了 —— 那时把氮加进孤对表
        // 这条判据照样绿,是空过的。
        let m = prep("C[N@](CC)CCC");
        let ranks = omgkit_io::canon::canonical_ranks(&m);
        let d = generate(&m, &Style::ACS_1996);
        let w = assign_wedges(&m, &d.coords, &ranks);
        assert!(
            w.bonds.iter().all(|x| x.narrow().is_none()),
            "给三配位氮画了楔形 —— 那个构型常温下不存在"
        );
    }

    #[test]
    fn when_the_two_taboos_conflict_the_better_of_both_orders_wins() {
        // 楔形该躲开**环键**,也该躲开**与另一个立体中心共用的键**。相邻的两个
        // 中心只有一根共用的非环键时这两条就冲了:谁先拿走谁受益,另一个被饿死。
        //
        // 下面两个分子各站一边 —— 单跑任何一种序都会在其中一个上吃亏:
        //
        // | | 只 SharedWorst | 只 RingWorst | 两种都跑 |
        // |---|---|---|---|
        // | 缩醛那个 | 2 根环上楔形 | **1 根** | 1 根 |
        // | 噻嗪那个 | 2 根环上楔形 | 1 根,但**丢一个中心** | 2 根 |
        //
        // 全量语料上:只 SharedWorst 159 根环上楔形 / 18 个 unwedged;
        // 只 RingWorst 156 / **19**;两种都跑 158 / 18 —— **信息不能丢**,
        // 所以先比 `unwedged` 再比环上楔形。
        for (smi, want_ring) in [
            ("CC1(OC[C@@H](O1)[C@@H]2CC23SCCCS3)C", 1usize),
            ("COCCCNC1=C(C(=O)N[C@@H](S1)[C@@H]2CCC=CC2)C#N", 2),
        ] {
            let m = prep(smi);
            let ranks = omgkit_io::canon::canonical_ranks(&m);
            let d = generate(&m, &Style::ACS_1996);
            let w = assign_wedges(&m, &d.coords, &ranks);
            assert!(
                w.unwedged.is_empty(),
                "{smi}:有中心没画出构型 {:?} —— 丢信息换环上楔形是亏的",
                w.unwedged
            );
            let on_ring = w
                .bonds
                .iter()
                .enumerate()
                .filter(|(bi, x)| {
                    x.narrow().is_some() && m.bonds()[*bi].flags.contains(BondFlags::IN_RING)
                })
                .count();
            assert_eq!(on_ring, want_ring, "{smi}:落在环键上的楔形数不对");
        }
    }

    #[test]
    fn a_stereocentre_gets_a_wedge() {
        for smi in [
            "N[C@@H](C)O",
            "N[C@H](C)O",
            "C[C@H](N)C(=O)O",
            "CN1CCC[C@H]1c1cccnc1",
        ] {
            let m = prep(smi);
            let ranks = omgkit_io::canon::canonical_ranks(&m);
            let d = generate(&m, &Style::ACS_1996);
            let w = assign_wedges(&m, &d.coords, &ranks);
            assert!(
                w.unwedged.is_empty(),
                "{smi} 有立体中心没画出构型:{:?}",
                w.unwedged
            );
            assert!(
                w.bonds.iter().any(|x| *x != Wedge::None),
                "{smi} 一根楔形键都没画"
            );
        }
    }

    #[test]
    fn a_molecule_and_its_mirror_image_get_opposite_wedges() {
        // **这是这一层唯一不空过的判据。**
        //
        // `assign_wedges` 是靠 `read_chirality` 挑方向的,所以"指派完再反读一致"
        // 必然成立 —— 拿它当判据是走过场。区分力才是真的:两个对映体的图完全
        // 全等(镜像),唯一的差别就在楔形的朝向。反读函数若与几何无关,两者
        // 会拿到同一个楔形,这一条立刻红。
        for (l, r) in [
            ("N[C@@H](C)O", "N[C@H](C)O"),
            ("C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O"),
            ("O[C@@H](Cl)Br", "O[C@H](Cl)Br"),
        ] {
            let (ml, mr) = (prep(l), prep(r));
            let rl = omgkit_io::canon::canonical_ranks(&ml);
            let rr = omgkit_io::canon::canonical_ranks(&mr);
            let dl = generate(&ml, &Style::ACS_1996);
            let dr = generate(&mr, &Style::ACS_1996);
            let wl = assign_wedges(&ml, &dl.coords, &rl);
            let wr = assign_wedges(&mr, &dr.coords, &rr);

            // 两者的坐标必须全等(对映体的 2D 图只差楔形)
            assert_eq!(dl.coords.len(), dr.coords.len());
            for (a, b) in dl.coords.iter().zip(&dr.coords) {
                assert!(a.dist(*b) < 1e-9, "{l} 与 {r} 的坐标不一致");
            }
            assert_ne!(
                wl.bonds, wr.bonds,
                "{l} 与 {r} 拿到了同一套楔形 —— 图上分不出对映体"
            );
        }
    }

    #[test]
    fn reading_back_gives_the_recorded_configuration() {
        // 这一条与 `assign` 有共谋(见 `assign_wedges` 的文档),但仍值得留着:
        // 它守的是"指派之后没有被别的步骤改坏",不是"指派本身对"。
        for smi in ["N[C@@H](C)O", "C[C@H](N)C(=O)O", "CN1CCC[C@H]1c1cccnc1"] {
            let m = prep(smi);
            let ranks = omgkit_io::canon::canonical_ranks(&m);
            let d = generate(&m, &Style::ACS_1996);
            let w = assign_wedges(&m, &d.coords, &ranks);
            for a in 0..u32::try_from(m.num_atoms()).unwrap() {
                let tag = m.atoms()[a as usize].chiral_tag;
                if matches!(tag, ChiralTag::Cw | ChiralTag::Ccw) {
                    assert_eq!(
                        read_chirality(&m, &d.coords, &w.bonds, a),
                        Some(tag),
                        "{smi} 的原子 {a} 反读回来不是原构型"
                    );
                }
            }
        }
    }

    #[test]
    fn a_centre_with_no_wedge_reads_as_unknown_rather_than_guessing() {
        // 一根楔形都没有时,图上根本读不出构型。返回 None 才诚实 —— 猜一个
        // 出来会让"反读一致"这条判据变成掷硬币。
        let m = prep("N[C@@H](C)O");
        let d = generate(&m, &Style::ACS_1996);
        let none = vec![Wedge::None; m.num_bonds()];
        let centre = (0..u32::try_from(m.num_atoms()).unwrap())
            .find(|a| {
                matches!(
                    m.atoms()[*a as usize].chiral_tag,
                    ChiralTag::Cw | ChiralTag::Ccw
                )
            })
            .expect("有立体中心");
        assert_eq!(read_chirality(&m, &d.coords, &none, centre), None);
    }

    #[test]
    fn cis_and_trans_are_read_back_from_the_geometry() {
        // 双键顺反画错是"画错"里最要命的一种:拓扑、原子数、键长全对,
        // 只有几何反了。这条把它变成可判定的。
        for smi in ["C/C=C/C", r"C/C=C\C", "C/C=C/CC", r"CC/C=C\C"] {
            // **顺反感知不在 `sanitize` 里**,要单独跑一步 —— 见模块文档。
            // 少了它,`stereo` 全是 `None`,这一档就成了空过(第一版正是如此,
            // 靠下面那个 `checked > 0` 守卫才发现)。
            let mut m = prep(smi);
            omgkit_io::stereo::perceive_bond_stereo(&mut m);
            let m = m;
            let d = generate(&m, &Style::ACS_1996);
            let mut checked = 0;
            for b in 0..u32::try_from(m.num_bonds()).unwrap() {
                let want = m.bonds()[b as usize].stereo;
                if !matches!(want, BondStereo::Cis | BondStereo::Trans) {
                    continue;
                }
                checked += 1;
                let got = read_bond_stereo(&m, &d.coords, b);
                assert_eq!(got, Some(want), "{smi} 的键 {b} 几何与记录不符");
            }
            assert!(
                checked > 0,
                "{smi} 一根带顺反的双键都没查到 —— 这一档在空过"
            );
            assert!(
                d.misdrawn_stereo.is_empty(),
                "{smi} 明明画对了,却被报成画错"
            );
        }
    }

    /// **画不出来的要说出来。** 环内的顺反就是画不出来的那一档。
    ///
    /// 掰顺反靠"把一侧整个镜像过去",而环上的键两侧是同一片原子,镜像动不了它;
    /// 环按凸多边形画,环内双键一律成顺式。八元以上的环里记着反式的,画出来必反。
    ///
    /// 先前这一档谁都不报 —— `fix_cis_trans` 的返回值被丢掉,四个诊断字段一个也
    /// 装不下,`is_clean()` 照样为真。断的是**契约**("画不好就说出来"),不是
    /// "环内顺反画不出来"这个当时的能力上限:哪天环画法能容下反式了,这条判据
    /// 会因为 `checked == 0` 而红,逼人来改,而不是悄悄空过。
    #[test]
    fn a_ring_double_bond_it_cannot_draw_is_reported_not_hidden() {
        let mut checked = 0;
        for n in 8..=14usize {
            // n 元环里一根记着反式的双键
            let smi = format!("C/1{}\\C=C1", "C".repeat(n - 3));
            let Ok(mut m) = omgkit_io::smiles::parse(&smi) else {
                continue;
            };
            if omgkit_chem::pipeline::sanitize(&mut m).is_err() {
                continue;
            }
            omgkit_io::stereo::perceive_bond_stereo(&mut m);
            let want: Vec<u32> = (0..u32::try_from(m.num_bonds()).unwrap())
                .filter(|&b| {
                    matches!(
                        m.bonds()[b as usize].stereo,
                        BondStereo::Cis | BondStereo::Trans
                    )
                })
                .collect();
            if want.is_empty() {
                continue;
            }
            let d = generate(&m, &Style::ACS_1996);
            for b in want {
                checked += 1;
                if read_bond_stereo(&m, &d.coords, b) != Some(m.bonds()[b as usize].stereo) {
                    assert!(
                        d.misdrawn_stereo.contains(&b),
                        "{smi}: 键 {b} 画反了,而诊断里没有它"
                    );
                    assert!(!d.is_clean(), "{smi}: 画反了还报 is_clean");
                }
            }
        }
        assert!(checked > 0, "一根环内的顺反键都没查到 —— 这一档在空过");
    }

    /// **配位几何画不出来,要报进 `unwedged`,不许一声不响地丢掉。**
    ///
    /// 先前挑立体中心那道 filter 只收 `Cw | Ccw`,于是 `@SP`/`@TB`/`@OH` 既不
    /// 指派楔形也不进诊断:`[Pt@SP1](Cl)(Cl)(N)N` 画出来读回去构型整个消失,
    /// 而 `is_clean()` 为真。
    ///
    /// 断的是契约,不是"这一版画不出配位几何"这个能力上限 —— 哪天画得出来了,
    /// 这条会因为中心不在 `unwedged` 里而红,逼人来改判据,而不是悄悄空过。
    #[test]
    fn a_coordination_geometry_it_cannot_draw_goes_into_unwedged() {
        let mut checked = 0;
        for smi in [
            "F[Pt@SP1](Cl)(Br)I",
            "S[As@TB1](F)(Cl)(Br)N",
            "O[Co@OH1](Cl)(C)(N)(F)P",
        ] {
            let Ok(mut m) = omgkit_io::smiles::parse(smi) else {
                continue;
            };
            if omgkit_chem::pipeline::sanitize(&mut m).is_err() {
                continue;
            }
            let centres: Vec<u32> = (0..u32::try_from(m.num_atoms()).unwrap())
                .filter(|a| {
                    let t = m.atoms()[*a as usize].chiral_tag;
                    t != ChiralTag::Unspecified && !t.is_tetrahedral()
                })
                .collect();
            if centres.is_empty() {
                continue;
            }
            let d = generate(&m, &Style::ACS_1996);
            for a in centres {
                checked += 1;
                assert!(
                    d.unwedged.contains(&a),
                    "{smi}: 原子 {a} 的配位构型没画出来,诊断里也没有它"
                );
            }
            assert!(!d.is_clean(), "{smi}: 构型丢了还报 is_clean");
        }
        assert!(checked > 0, "一个配位中心都没查到 —— 这一档在空过");
    }
}
