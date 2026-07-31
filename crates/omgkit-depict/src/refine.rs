//! 消冲突:先用离散算子,不轻易动几何。
//!
//! # 算子的优先级
//!
//! 沿用 Shelly (1983) / Helson (1999) 的算子集,顺序取 CoordGen 的经验 ——
//! **先试简单且更美观的**:
//!
//! | 算子 | 代价 | 本模块 |
//! |---|---|---|
//! | **翻转**(可旋转键的一侧整体镜像) | 键长、键角**一点不变** | 实现 |
//! | 开角 | 键角偏离理想值 | 未做 |
//! | 伸缩键长 | 键长不再全等 | 未做 |
//!
//! 只做翻转是有意的:它是唯一**不损失任何几何性质**的算子,而后两个一旦用上,
//! "键长全等""键角标准"这两条判据就守不住了。翻转解决不了的情形如实报出来
//! ([`Report::unresolved`]),不靠拉扯几何把数字做好看。
//!
//! # 碰撞半径来自标签,而标签尺寸来自规范
//!
//! 两个原子撞没撞上,取决于它们的**标签**占多大 —— 这就是
//! [`Style`] 必须参与布局的地方。同一张图在 ACS 规范下
//! (标签占 0.69 个键长)会撞,在 ChemDraw 默认规范下(0.33)可能不撞。

use std::collections::{BTreeMap, BTreeSet};

use omgkit_core::{BondFlags, MolBuilder};

use crate::geom::{segments_cross, Point2};
use crate::label::{label_for, HSide};
use crate::style::Style;

/// 没有标签的骨架碳,占位半径。
///
/// 两个裸碳中心靠得比一个键长的一半还近,图上就读不清了 —— 取 0.25 使阈值恰好
/// 落在半个键长。纯靠标签留白算出来的半径(ACS 下约 0.11)太小,骨架上的碰撞
/// 会漏判。
const BARE_RADIUS: f64 = 0.25;

/// 一次消冲突的结果。
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Report {
    /// 翻转过的键
    pub flipped: Vec<u32>,
    /// **仍然撞着的原子对**。翻转解决不了的情形如实留在这里,不假装消掉了。
    pub unresolved: Vec<(u32, u32)>,
    /// 仍然交叉的键对
    pub crossings: Vec<(u32, u32)>,
}

/// 反复翻转可旋转键,直到没有改善。
///
/// `pos` 就地修改。返回仍未解决的部分。
pub(crate) fn relieve(
    mol: &MolBuilder,
    pos: &mut BTreeMap<u32, Point2>,
    ranks: &[u32],
    style: &Style,
) -> Report {
    let radii = radii(mol, style);
    let mut best = score(mol, pos, &radii);
    let mut report = Report::default();
    if best.0 == 0.0 && best.1 == 0 {
        return report;
    }

    // 候选键按规范秩排 —— 拿存储下标排,同一分子的不同写法会翻不同的键
    let mut cands: Vec<u32> = rotatable(mol, pos);
    cands.sort_by_key(|&b| {
        let bd = &mol.bonds()[b as usize];
        let (x, y) = (ranks[bd.begin as usize], ranks[bd.end as usize]);
        (x.min(y), x.max(y), b)
    });

    // 上限防止在数值抖动上来回翻。每轮至多改善一次,轮数够覆盖所有键即可。
    let max_rounds = cands.len().max(1) * 2;
    for _ in 0..max_rounds {
        let mut improved = false;
        for &b in &cands {
            let Some(side) = far_side(mol, pos, b) else {
                continue;
            };
            let bd = &mol.bonds()[b as usize];
            let (u, v) = (pos[&bd.begin], pos[&bd.end]);

            let saved: Vec<(u32, Point2)> = side.iter().map(|a| (*a, pos[a])).collect();
            for a in &side {
                let p = pos[a].mirrored(u, v - u);
                pos.insert(*a, p);
            }
            // **立体守卫**:翻转会把双键旁的参照原子换到另一侧,顺反跟着反。
            // 消掉一处碰撞、同时把 Z 画成 E,是拿"看着好一点"换"画错了"。
            let now = score(mol, pos, &radii);
            let keeps_stereo = {
                let mut flat = vec![Point2::ORIGIN; mol.num_atoms()];
                for (a, q) in pos.iter() {
                    flat[*a as usize] = *q;
                }
                crate::stereo::cis_trans_intact(mol, &flat)
            };
            if keeps_stereo && better(now, best) {
                best = now;
                report.flipped.push(b);
                improved = true;
            } else {
                for (a, p) in saved {
                    pos.insert(a, p);
                }
            }
        }
        if !improved {
            break;
        }
    }

    let (pairs, crossings) = remaining(mol, pos, &radii);
    report.unresolved = pairs;
    report.crossings = crossings;
    report
}

/// 每个原子的碰撞半径,单位是键长。
pub(crate) fn radii(mol: &MolBuilder, style: &Style) -> Vec<f64> {
    (0..mol.num_atoms())
        .map(|i| {
            let a = u32::try_from(i).expect("原子数超出 u32");
            // 氢挂哪一侧此刻还定不下来(它要看最终坐标),取两侧中更宽的那个 ——
            // 半径宁可**偏大**:偏大只是把原子推得开一点,偏小会漏判碰撞
            [HSide::Right, HSide::Left]
                .iter()
                .filter_map(|s| label_for(mol, a, style, *s))
                .map(|l| l.half_w.hypot(l.half_h))
                .fold(BARE_RADIUS, f64::max)
        })
        .collect()
}

/// `now` 是不是**确实**比 `best` 好。
///
/// # 为什么不能直接写 `now < best`
///
/// 碰撞深度是一串浮点求和,而求和次序取决于原子编号 —— 编号又取决于 SMILES
/// 怎么写。一次对深度毫无影响的翻转,两种写法算出来会是
/// `0.2499999999999985` 与 `0.2499999999999984`,一边判"更好"接受、另一边
/// 判"更差"拒绝。同一个分子于是画成两张图。
///
/// 实测:阿司匹林的两种写法就差在这一位上,布局、候选序、翻转的那一侧全都
/// 一模一样,只有第 16 位有效数字不同。
///
/// 所以要带容差,而且**分不出高下时不动** —— 无意义的翻转不做,保守方向。
fn better(now: (f64, usize), best: (f64, usize)) -> bool {
    const EPS: f64 = 1e-9;
    if now.1 != best.1 {
        return now.1 < best.1;
    }
    now.0 < best.0 - EPS
}

/// 打分。**越小越好**,按(碰撞深度平方和, 交叉键对数)的字典序比较。
///
/// 交叉排在后面而不是加权求和:两者量纲不同,加权就要引入一个说不清的系数,
/// 而系数一变结论就变。字典序不需要系数。
fn score(mol: &MolBuilder, pos: &BTreeMap<u32, Point2>, radii: &[f64]) -> (f64, usize) {
    let (pairs, crossings) = remaining(mol, pos, radii);
    let depth: f64 = pairs
        .iter()
        .map(|(i, j)| {
            let want = radii[*i as usize] + radii[*j as usize];
            let d = pos[i].dist(pos[j]);
            (want - d).max(0.0).powi(2)
        })
        .sum();
    (depth, crossings.len())
}

/// 仍在碰撞的原子对,与仍然交叉的键对。
type Trouble = (Vec<(u32, u32)>, Vec<(u32, u32)>);

/// 仍在碰撞的原子对与交叉的键对。
fn remaining(mol: &MolBuilder, pos: &BTreeMap<u32, Point2>, radii: &[f64]) -> Trouble {
    let bonded: BTreeSet<(u32, u32)> = mol
        .bonds()
        .iter()
        .map(|b| (b.begin.min(b.end), b.begin.max(b.end)))
        .collect();
    let mut atoms: Vec<u32> = pos.keys().copied().collect();
    atoms.sort_unstable();

    let mut pairs = Vec::new();
    for (k, &i) in atoms.iter().enumerate() {
        for &j in &atoms[k + 1..] {
            if bonded.contains(&(i, j)) {
                continue;
            }
            let want = radii[i as usize] + radii[j as usize];
            if pos[&i].dist(pos[&j]) < want {
                pairs.push((i, j));
            }
        }
    }

    // 键交叉:只看两端都已放置的键
    let live: Vec<u32> = (0..mol.num_bonds())
        .map(|i| u32::try_from(i).expect("键数超出 u32"))
        .filter(|b| {
            let bd = &mol.bonds()[*b as usize];
            pos.contains_key(&bd.begin) && pos.contains_key(&bd.end)
        })
        .collect();
    let mut crossings = Vec::new();
    for (k, &b1) in live.iter().enumerate() {
        for &b2 in &live[k + 1..] {
            let (x, y) = (&mol.bonds()[b1 as usize], &mol.bonds()[b2 as usize]);
            if segments_cross(pos[&x.begin], pos[&x.end], pos[&y.begin], pos[&y.end]) {
                crossings.push((b1, b2));
            }
        }
    }
    (pairs, crossings)
}

/// 可翻转的键:不在环里,且两端度数都大于 1。
///
/// 环上的键翻不动(翻了会把环撕开);端点键翻了等于什么都没做 —— 那一侧只有
/// 它自己,镜像回原位。
fn rotatable(mol: &MolBuilder, pos: &BTreeMap<u32, Point2>) -> Vec<u32> {
    (0..mol.num_bonds())
        .map(|i| u32::try_from(i).expect("键数超出 u32"))
        .filter(|b| {
            let bd = &mol.bonds()[*b as usize];
            !bd.flags.contains(BondFlags::IN_RING)
                && mol.degree(bd.begin) > 1
                && mol.degree(bd.end) > 1
                && pos.contains_key(&bd.begin)
                && pos.contains_key(&bd.end)
        })
        .collect()
}

/// 断开键 `b` 之后,**原子更少**的那一侧。
///
/// 返回 `None` 表示这根键其实在环上(断开后两端仍连通),翻不得。环感知标记
/// 之外再判一次是**故意的**:标记来自净化,而调用方未必净化过。
///
/// # 为什么必须取更少的那一侧,而不是 `end` 那一侧
///
/// 哪一端是 `end` 依**写法**而定。取 `end` 那一侧的话,同一根化学键在一种写法
/// 里镜像的是一个甲基、在另一种写法里镜像的是整个苯环 —— 两者相差一次全局
/// 反射,而接受的翻转**次数**也会跟着不同,最后坐标就对不上。
///
/// 实测:阿司匹林的两种写法,布局阶段已经完全一致了,却在这里分岔 ——
/// 一种翻两次、另一种翻一次。
///
/// 取更少的那一侧也更自然:不该为了挪一个甲基把整个分子翻过来。
fn far_side(mol: &MolBuilder, pos: &BTreeMap<u32, Point2>, b: u32) -> Option<Vec<u32>> {
    let bd = &mol.bonds()[b as usize];
    let (start, blocked) = (bd.end, bd.begin);
    let mut seen: BTreeSet<u32> = BTreeSet::from([start]);
    let mut stack = vec![start];
    while let Some(a) = stack.pop() {
        for (n, bi) in mol.neighbors(a) {
            if bi == b || !pos.contains_key(&n) {
                continue;
            }
            if n == blocked {
                return None; // 绕回去了 —— 是环上的键
            }
            if seen.insert(n) {
                stack.push(n);
            }
        }
    }
    let mut out: Vec<u32> = seen.into_iter().collect();
    out.sort_unstable();

    // 取更少的那一侧。平局(两侧一样多)时按**这一侧最小的规范秩**定 ——
    // 拿存储下标定就又把写法依赖引回来了。
    let placed = pos.len();
    if out.len() * 2 > placed {
        let other: Vec<u32> = {
            let mine: BTreeSet<u32> = out.iter().copied().collect();
            let mut o: Vec<u32> = pos.keys().copied().filter(|a| !mine.contains(a)).collect();
            o.sort_unstable();
            o
        };
        if !other.is_empty() {
            return Some(other);
        }
    }
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::layout;

    fn prep(smi: &str) -> MolBuilder {
        let mut m = omgkit_io::smiles::parse(smi).unwrap();
        omgkit_chem::pipeline::sanitize(&mut m).unwrap();
        m
    }

    fn laid(smi: &str, style: &Style) -> (MolBuilder, BTreeMap<u32, Point2>, Report) {
        let m = prep(smi);
        let ranks = omgkit_io::canon::canonical_ranks(&m);
        let mut pos: BTreeMap<u32, Point2> = BTreeMap::new();
        for p in layout::layout_all(&m, &ranks, style) {
            pos.extend(p.pos);
        }
        let rep = relieve(&m, &mut pos, &ranks, style);
        (m, pos, rep)
    }

    #[test]
    fn the_aspirin_overlap_is_gone() {
        // 这是判据先抓到、再修的那一个:`OC(=O)c1ccccc1OC(C)=O` 的两个羰基氧
        // 原本落在同一个点 (2.50, -0.87) 上,距离 0.0000。
        let (m, pos, rep) = laid("OC(=O)c1ccccc1OC(C)=O", &Style::ACS_1996);
        let n = m.num_atoms();
        for i in 0..n {
            for j in (i + 1)..n {
                let (a, b) = (i as u32, j as u32);
                assert!(pos[&a].dist(pos[&b]) > 0.3, "原子 {a} 与 {b} 仍然几乎重合");
            }
        }
        assert!(!rep.flipped.is_empty(), "应当至少翻转过一根键");
    }

    #[test]
    fn flipping_keeps_every_bond_exactly_one_unit() {
        // 翻转是唯一**不损失几何性质**的算子。若不小心写成了缩放或平移,
        // 冲突照样能消掉,但键长会悄悄变 —— 那正是选它而不选伸缩的理由。
        for smi in [
            "OC(=O)c1ccccc1OC(C)=O",
            "CC(C)(C)c1ccccc1C(C)(C)C",
            "CCCCCCCC",
        ] {
            let (m, pos, _) = laid(smi, &Style::ACS_1996);
            for b in m.bonds() {
                let d = pos[&b.begin].dist(pos[&b.end]);
                assert!(
                    (d - 1.0).abs() < 1e-9,
                    "{smi} 键 {}–{} 长 {d}",
                    b.begin,
                    b.end
                );
            }
        }
    }

    #[test]
    fn a_ring_bond_is_never_flipped() {
        // 翻环上的键会把环撕开。`far_side` 在环感知标记之外再判一次连通性,
        // 这里守的是那一层。
        let m = prep("c1ccccc1");
        let ranks = omgkit_io::canon::canonical_ranks(&m);
        let mut pos: BTreeMap<u32, Point2> = BTreeMap::new();
        for p in layout::layout_all(&m, &ranks, &Style::ACS_1996) {
            pos.extend(p.pos);
        }
        for b in 0..u32::try_from(m.num_bonds()).unwrap() {
            assert!(
                far_side(&m, &pos, b).is_none(),
                "环上的键 {b} 不该给出可翻的一侧"
            );
        }
    }

    #[test]
    fn what_cannot_be_fixed_is_reported_not_hidden() {
        // 翻转解决不了的必须留在 `unresolved` 里。悄悄清空它,图上还是挤的,
        // 而调用方以为一切正常 —— 那比报出来糟得多。
        //
        // 六个叔丁基围着一个苯环,平面上无论如何都排不开。
        let (_, _, rep) = laid(
            "CC(C)(C)c1c(C(C)(C)C)c(C(C)(C)C)c(C(C)(C)C)c(C(C)(C)C)c1C(C)(C)C",
            &Style::ACS_1996,
        );
        assert!(
            !rep.unresolved.is_empty(),
            "六个叔丁基挤在一个苯环上不可能全排开,却报告说没有冲突"
        );
    }

    #[test]
    fn the_two_styles_can_disagree_about_whether_it_clashes() {
        // **这是 Style 参与布局的落点。** 同一张图,ACS 的标签占 0.69 个键长、
        // ChemDraw 默认占 0.33 —— 判定必须不同,否则把 Style 传进来就是白传。
        let m = prep("OC(=O)c1ccccc1OC(C)=O");
        let ranks = omgkit_io::canon::canonical_ranks(&m);
        let mut pos: BTreeMap<u32, Point2> = BTreeMap::new();
        for p in layout::layout_all(&m, &ranks, &Style::ACS_1996) {
            pos.extend(p.pos);
        }
        let acs = remaining(&m, &pos, &radii(&m, &Style::ACS_1996)).0.len();
        let cd = remaining(&m, &pos, &radii(&m, &Style::CHEMDRAW_DEFAULT))
            .0
            .len();
        assert!(
            acs >= cd,
            "ACS 的标签更大,判出的碰撞不该少于 ChemDraw 默认:{acs} vs {cd}"
        );
        assert!(
            radii(&m, &Style::ACS_1996)[0] > radii(&m, &Style::CHEMDRAW_DEFAULT)[0],
            "ACS 的碰撞半径应当更大"
        );
    }
}
