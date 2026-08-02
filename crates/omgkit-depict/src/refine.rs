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
//! | **螺环翻转**(绕螺原子把一个环整体镜像) | 同上 | 实现 |
//! | **度 4 置换**(把两个分支对调) | 同上 | 实现 |
//! | 开角 | 键角偏离理想值 | 未做 |
//! | 伸缩键长 | 键长不再全等 | 未做 |
//!
//! 只做前三个是有意的:它们都是**不损失任何几何性质**的等距变换(轴过某个不动的
//! 原子,反射保持到它的距离),而后两个一旦用上,"键长全等""键角标准"这两条
//! 判据就守不住了。解决不了的情形如实报出来([`Report::unresolved`]),不靠
//! 拉扯几何把数字做好看。
//!
//! 后两个各补一个**能力空洞**,收益都要如实说(全量语料 17662 个分子×规范):
//!
//! | 算子 | 补的空洞 | 转干净 | 键交叉 |
//! |---|---|---:|---:|
//! | 螺环翻转 | 可翻转的键排除环上的键,螺环两侧的相对朝向动不了 | **+6** | 0 |
//! | 度 4 置换 | 垂直于翻转轴的那一对邻居,怎么翻都换不过来 | **+25** | −2 |
//!
//! 干净率 91.4% → **91.6%**。做它们是因为便宜且不碰契约,不是因为它们大。
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
    /// **翻转过的螺原子**。绕螺原子把一个环整体镜像过去,见 `SpiroFlip`。
    pub spiro_flipped: Vec<u32>,
    /// **做过分支对调的度 4 原子**,见 `PairSwap`。
    pub swapped: Vec<u32>,
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

    // 螺环翻转与度 4 置换的候选。**枚举一次就够** —— 只依赖拓扑,不依赖坐标。
    let spiros = spiro_flips(mol, pos, ranks);
    let swaps = pair_swaps(mol, pos, ranks);

    // 立体守卫抽出来:两个算子用同一套判定,免得改了一处漏了另一处
    let keeps_stereo = |pos: &BTreeMap<u32, Point2>| {
        let mut flat = vec![Point2::ORIGIN; mol.num_atoms()];
        for (a, q) in pos.iter() {
            flat[*a as usize] = *q;
        }
        crate::stereo::cis_trans_intact(mol, &flat)
    };

    // 上限防止在数值抖动上来回翻。每轮至多改善一次,轮数够覆盖所有键即可。
    let max_rounds = (cands.len() + spiros.len() + swaps.len()).max(1) * 2;
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
            if keeps_stereo(pos) && better(now, best) {
                best = now;
                report.flipped.push(b);
                improved = true;
            } else {
                for (a, p) in saved {
                    pos.insert(a, p);
                }
            }
        }

        // 螺环翻转。放在键翻转之后:键翻转更常见也更便宜,先让它把能消的消掉。
        for sf in &spiros {
            let (Some(c), Some(p1), Some(p2)) = (
                pos.get(&sf.centre).copied(),
                pos.get(&sf.ring_nbrs.0).copied(),
                pos.get(&sf.ring_nbrs.1).copied(),
            ) else {
                continue;
            };
            // 轴每轮重算 —— 前面的键翻转可能已经把这两个邻居挪过了
            let axis = (p1 + p2) * 0.5 - c;
            if axis.norm() < 1e-9 {
                continue; // 两个邻居正好对称在螺原子两侧,轴退化,镜像没有意义
            }
            // **两根环键一样长,这次反射才真的是"把 N1 与 N2 对调"。**
            //
            // 先说清楚一件容易想岔的事:**键长与轴怎么选无关**。绕任何一条过 S
            // 的直线反射,都保持每个点到 S 的距离,而跨越边界的键只有 S–N1 与
            // S–N2 两根 —— 所以这个算子无论轴选哪条都是保键长的,被反射那一侧
            // 内部的键长键角更是原样(等距变换)。
            //
            // 中点这条轴买到的是别的东西:`|S−N1| = |S−N2|` 时它恰好是 ∠N1SN2 的
            // 角平分线,反射把 N1 与 N2 精确对调 —— 于是这个环落在自己的镜像位置
            // 上,是一次**对称操作**,而不是把环转到一个任意朝向。长度不等时
            // (退化布局里就会)中点偏向短的那一边,对调不成立,环会歪到一个没有
            // 道理的角度上;`better()` 仍然会拦下变差的结果,但那已经不是这个
            // 算子想做的事了。
            //
            // **这条守卫没有被全量语料证伪过也没有被证实过**:去掉它,语料上
            // 键长全等仍是 0 违例、干净率一个不差。留着是因为上面那个论证要它,
            // 不是因为量到了收益。
            if (c.dist(p1) - c.dist(p2)).abs() > 1e-9 {
                continue;
            }
            let saved: Vec<(u32, Point2)> = sf.side.iter().map(|a| (*a, pos[a])).collect();
            for a in &sf.side {
                let p = pos[a].mirrored(c, axis);
                pos.insert(*a, p);
            }
            let now = score(mol, pos, &radii);
            if keeps_stereo(pos) && better(now, best) {
                best = now;
                report.spiro_flipped.push(sf.centre);
                improved = true;
            } else {
                for (a, p) in saved {
                    pos.insert(a, p);
                }
            }
        }

        // 度 4 置换。排最后:它最贵(每个结点六对),而前两个算子能消的先消掉。
        for sw in &swaps {
            let (Some(c), Some(pa), Some(pb)) = (
                pos.get(&sw.centre).copied(),
                pos.get(&sw.ends.0).copied(),
                pos.get(&sw.ends.1).copied(),
            ) else {
                continue;
            };
            let axis = (pa + pb) * 0.5 - c;
            if axis.norm() < 1e-9 {
                continue; // 两个邻居正对着,轴退化
            }
            // 两根键一样长,反射才真的是"把 A 与 B 对调" —— 与螺环那条同一个道理
            if (c.dist(pa) - c.dist(pb)).abs() > 1e-9 {
                continue;
            }
            let saved: Vec<(u32, Point2)> = sw
                .sides
                .0
                .iter()
                .chain(sw.sides.1.iter())
                .map(|a| (*a, pos[a]))
                .collect();
            for a in sw.sides.0.iter().chain(sw.sides.1.iter()) {
                let p = pos[a].mirrored(c, axis);
                pos.insert(*a, p);
            }
            let now = score(mol, pos, &radii);
            if keeps_stereo(pos) && better(now, best) {
                best = now;
                report.swapped.push(sw.centre);
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

/// 打分。**越小越好**,按(交叉键对数, 碰撞深度平方和)的字典序比较 ——
/// 注意**交叉在前**,见 [`better`]。
///
/// 用字典序而不是加权求和:两者量纲不同,加权就要引入一个说不清的系数,而
/// 系数一变结论就变。字典序不需要系数。
///
/// 交叉排在深度前面是量出来的,不是拍的:把两者对调之后,全语料上有交叉的
/// 图从 381 涨到 415,而未解冲突一个没少。
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

/// 一次螺环翻转:绕螺原子把某一个环连同挂在它上面的一切镜像过去。
///
/// # 为什么这个算子不损失任何几何
///
/// 螺原子 `S` 在目标环里有两个邻居 `N1`、`N2`,两根都是环上的键,所以
/// `|S−N1| = |S−N2|` —— **`S` 到两个邻居中点的连线恰好就是 ∠N1SN2 的角平分线**。
/// 绕它反射,`N1` 与 `N2` 精确对调。
///
/// 被翻的那一侧(整个环 + 稠上去的环 + 挂着的全部取代基)受的是**同一个反射**,
/// 而反射是等距变换 —— 那一侧内部所有键长、键角逐个不变;跨界的键只有 `S–N1`
/// 与 `S–N2` 两根,`S` 在轴上不动,两根长度都保持,`S` 处的键角多重集也保持。
///
/// 所以它与"绕可旋转键翻转"同级:**不损失任何几何性质**,符合本模块只用这类
/// 算子的规矩。
///
/// # 为什么现在的算子够不着这些图
///
/// [`rotatable`] 明确排除环上的键(翻了会把环撕开),而螺原子两侧全是环上的键
/// —— **螺环两侧的相对朝向根本动不了**。实测:含螺原子的图在基线里占 1.2%,
/// 在有键交叉的图里占 13.7%(富集 11.7 倍),在有未解冲突的图里占 8.6%。
///
/// # 两处平局都不许照抄 RDKit
///
/// RDKit 的 `flipAboutSpiroCenter` 固定翻 `rings[0]`,而 `rings` 来自
/// `getRingInfo()->atomRings()` 是**存储序**;要翻哪一侧又是从
/// `atomNeighbors` 的第一个取的,同样是存储序。照抄这两处,同一个分子换种
/// 写法就会翻不同的环、得到不同的图。这里两处都按规范秩定。
struct SpiroFlip {
    /// 螺原子
    centre: u32,
    /// 要镜像的那一侧
    side: Vec<u32>,
    /// 螺原子在目标环里的两个邻居。**存编号不存坐标** —— 前面的键翻转会改坐标,
    /// 枚举时算好的中点到这一步就过期了。
    ring_nbrs: (u32, u32),
}

/// 枚举所有可做的螺环翻转,**按规范秩定序**。
fn spiro_flips(mol: &MolBuilder, pos: &BTreeMap<u32, Point2>, ranks: &[u32]) -> Vec<SpiroFlip> {
    let rings = omgkit_chem::sssr::ring_set(mol);
    let mut out: Vec<(Vec<u32>, SpiroFlip)> = Vec::new();

    for a in 0..u32::try_from(mol.num_atoms()).expect("原子数超出 u32") {
        if !pos.contains_key(&a) {
            continue;
        }
        // 含这个原子的环,按规范秩定序 —— 不看 SSSR 给出来的次序
        let mut mine: Vec<&omgkit_chem::sssr::Ring> =
            rings.iter().filter(|r| r.atoms.contains(&a)).collect();
        if mine.len() < 2 {
            continue;
        }
        mine.sort_by_key(|r| ring_key(r, ranks));

        for r in &mine {
            // 这个环必须**只**和其它含 a 的环共用 a 自己 —— 否则不是螺,是稠合
            let set: BTreeSet<u32> = r.atoms.iter().copied().collect();
            let spiro_here = mine.iter().any(|o| {
                !std::ptr::eq(*o, *r) && o.atoms.iter().filter(|x| set.contains(x)).count() == 1
            });
            if !spiro_here {
                continue;
            }
            // a 在这个环里的两个邻居
            let nbrs: Vec<u32> = r
                .atoms
                .iter()
                .copied()
                .filter(|x| mol.neighbors(a).any(|(n, _)| n == *x))
                .collect();
            if nbrs.len() != 2 {
                continue;
            }
            if !pos.contains_key(&nbrs[0]) || !pos.contains_key(&nbrs[1]) {
                continue;
            }
            // 从其中一个邻居出发、把 a 挡住,收集这一侧
            let Some(side) = one_side(mol, pos, nbrs[0], a) else {
                continue;
            };
            // **另一个环内邻居必须在这一侧**(顺着环绕回来)。不在的话说明这两
            // 个环之间还有别的通路,a 不是割点,翻了会把结构撕开。
            if !side.contains(&nbrs[1]) {
                continue;
            }
            // 这一侧不许覆盖到别的环内邻居之外的东西 —— 上一条已经保证连通性,
            // 这里只需再确认 a 自己没被卷进去
            if side.contains(&a) {
                continue;
            }
            let key: Vec<u32> = {
                let mut k = vec![ranks[a as usize]];
                k.extend(ring_key(r, ranks));
                k
            };
            out.push((
                key,
                SpiroFlip {
                    centre: a,
                    side,
                    ring_nbrs: (nbrs[0], nbrs[1]),
                },
            ));
        }
    }
    out.sort_by(|x, y| x.0.cmp(&y.0));
    out.into_iter().map(|x| x.1).collect()
}

/// 环的确定性排序键:环上规范秩的有序多重集。
fn ring_key(r: &omgkit_chem::sssr::Ring, ranks: &[u32]) -> Vec<u32> {
    let mut k: Vec<u32> = r.atoms.iter().map(|a| ranks[*a as usize]).collect();
    k.sort_unstable();
    k
}

/// 从 `start` 出发、把 `blocked` 挡住能走到的那些原子(不含 `blocked`)。
///
/// 走不通(图里根本没放置这些原子)时返回 `None`。
fn one_side(
    mol: &MolBuilder,
    pos: &BTreeMap<u32, Point2>,
    start: u32,
    blocked: u32,
) -> Option<Vec<u32>> {
    if !pos.contains_key(&start) {
        return None;
    }
    let mut seen: BTreeSet<u32> = BTreeSet::from([start]);
    let mut stack = vec![start];
    while let Some(x) = stack.pop() {
        for (n, _) in mol.neighbors(x) {
            if n == blocked || !pos.contains_key(&n) {
                continue;
            }
            if seen.insert(n) {
                stack.push(n);
            }
        }
    }
    Some(seen.into_iter().collect())
}

/// 一次度 4 置换:把中心原子的两个分支对调。
///
/// # 为什么它是翻转够不着的
///
/// 绕键 `C–D` 翻转,镜像的是远侧那一整块,轴是 `C–D` 这条线 —— 它对调的是
/// **关于这条轴对称的那一对**邻居。垂直于轴的那一对怎么翻都换不过来。
/// RDKit 的注释把这一点说得很明白,并且据此只把这个算子用在度 4 的结点上:
/// 度 3 时绕第三根键翻转正好就能对调另外两根,用不着单独的算子。
///
/// # 它同样不损失几何
///
/// 轴过中心 `C`,反射保持每个点到 `C` 的距离,所以四根键的长度都不变;
/// 被反射的两块内部是等距变换,块与块之间也是同一个反射 —— 距离全保。
/// `|CA| = |CB|` 时轴恰好是 ∠ACB 的角平分线,`A` 与 `B` 精确对调,
/// `C` 处的键角多重集也保持,这才是"置换"该有的样子。
///
/// # 判据不照抄 RDKit
///
/// `findBondsPairsToPermuteDeg4` 用 `fabs(dp) < 1e-3` 判两根键是否垂直,
/// **硬假设度 4 结点画成 90° 十字**。omgkit 的 `chains::allocate` 在已占方向
/// ≥2 时是劈最大空隙,`free_direction` 还会按 ±30° 躲让 —— 十字没有保证。
/// 照抄那个点积判据会静默掉进 `else` 分支,返回的恰好是**本来就能被翻转达到**
/// 的那一对,算子等于白做。
///
/// 这里改成**六对全枚举**,让 [`better`] 自己判优:冗余的那几对至多是空转,
/// 不会给出错的结果。
struct PairSwap {
    /// 中心原子
    centre: u32,
    /// 要对调的两个邻居
    ends: (u32, u32),
    /// 两个分支各自的原子
    sides: (Vec<u32>, Vec<u32>),
}

/// 枚举所有可做的度 4 置换,**按规范秩定序**。
fn pair_swaps(mol: &MolBuilder, pos: &BTreeMap<u32, Point2>, ranks: &[u32]) -> Vec<PairSwap> {
    let mut out: Vec<(Vec<u32>, PairSwap)> = Vec::new();
    for c in 0..u32::try_from(mol.num_atoms()).expect("原子数超出 u32") {
        if mol.degree(c) != 4 || !pos.contains_key(&c) {
            continue;
        }
        // 中心原子完全不在环上 —— 与 RDKit 的守卫同一条,但它在这里只是个**便宜
        // 的前置过滤**,不是安全底线。真正兜底的是下面那句"两块必须不相交":
        // 去掉这一句,判据 `a_deg4_swap_never_tears_the_structure_apart` 照样绿,
        // 因为重叠的那些对会被显式查出来跳掉。
        //
        // 留着它是因为便宜(一次遍历换掉一堆无用的 `one_side` 搜索),而不是
        // 因为少了它会出错 —— 这一点变异验证说得很清楚。
        if mol
            .neighbors(c)
            .any(|(_, bi)| mol.bonds()[bi as usize].flags.contains(BondFlags::IN_RING))
        {
            continue;
        }
        let mut nbrs: Vec<u32> = mol.neighbors(c).map(|(x, _)| x).collect();
        nbrs.sort_by_key(|x| (ranks[*x as usize], *x));
        nbrs.dedup();
        if nbrs.len() != 4 {
            continue;
        }
        for i in 0..4 {
            for j in (i + 1)..4 {
                let (a, b) = (nbrs[i], nbrs[j]);
                let (Some(sa), Some(sb)) = (one_side(mol, pos, a, c), one_side(mol, pos, b, c))
                else {
                    continue;
                };
                // 两块必须**不相交**,而且都不含中心
                let set_a: BTreeSet<u32> = sa.iter().copied().collect();
                if sa.contains(&c) || sb.contains(&c) || sb.iter().any(|x| set_a.contains(x)) {
                    continue;
                }
                out.push((
                    vec![
                        ranks[c as usize],
                        ranks[a as usize].min(ranks[b as usize]),
                        ranks[a as usize].max(ranks[b as usize]),
                    ],
                    PairSwap {
                        centre: c,
                        ends: (a, b),
                        sides: (sa, sb),
                    },
                ));
            }
        }
    }
    out.sort_by(|x, y| x.0.cmp(&y.0));
    out.into_iter().map(|x| x.1).collect()
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

    // 另一侧:从 `blocked` 出发、同样不跨这根键能走到的那些原子。
    //
    // **不能拿 `pos` 取补集。** `relieve` 跑在所有分量合并之后的 `pos` 上
    // (分量之间也可能撞上),补集因此会**跨过分量边界** —— 翻一根键会把另一个
    // 不相干的片段整个镜像走。实测这种翻转从没被 `better()` 接受过(它必然
    // 制造碰撞),所以图上看不出问题;但 `far_side` 返回的东西本身就是错的,
    // 判据 `a_flip_never_reaches_into_another_fragment` 一上来就把它照出来了。
    let mut other_seen: BTreeSet<u32> = BTreeSet::from([blocked]);
    let mut stack = vec![blocked];
    while let Some(a) = stack.pop() {
        for (n, bi) in mol.neighbors(a) {
            if bi == b || !pos.contains_key(&n) {
                continue;
            }
            if other_seen.insert(n) {
                stack.push(n);
            }
        }
    }
    let mut other: Vec<u32> = other_seen.into_iter().collect();
    other.sort_unstable();

    // 取更少的那一侧。平局(两侧一样多)时按**这一侧最小的规范秩**定 ——
    // 拿存储下标定就又把写法依赖引回来了。
    if out.len() > other.len() && !other.is_empty() {
        return Some(other);
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

    /// 螺环分子。前几个取自语料,后几个是教科书上的螺环。
    const SPIRO: [&str; 8] = [
        // 前六个是从全量语料里**扫出来确实会触发螺环翻转**的。判据要求至少
        // 触发一次,所以这批不能随手换成"看着像螺环"的分子 —— 第一版列了八个
        // 教科书螺环,一个都没触发,那时"键长不变"那条是恒真的。
        "c1ccc2c(c1)C3C[C@@]4(C2c5c3cccc5)C=CS4(=O)=O",
        "c1ccc2c(c1)CN(C(=[NH2+])C23CCOCC3)N",
        "C1COC2(O1)[C@@]3(C[C@@]3(C(=[NH+]2)N)C#N)C#N",
        "Cc1cccc(c1)NC2=C(C(=O)NC3(S2)CCCC3)C#N",
        "Cc1ccccc1[C@@H]2[C@@H]([C@@]23C(=NN=C3O)N)C#N",
        "N1=C(SC2(CCCCC2)C3=C1CCCC3)N",
        // 这两个不触发,留着让"等距"那条也覆盖到不触发的路径
        "C1CC2(CC1)CCCC2",
        "C[C@@]12C(=C[C@@H](O1)C(=O)C23CC3)C(=O)OC",
    ];

    /// 度 4 置换会触发的分子。**从全量语料扫出来的**,不是挑的。
    const SWAP: [&str; 6] = [
        "[O-][N+](=O)O[Cu](O[N+]([O-])=O)([N+]1=C2C=CC=CC2=CC=C1)[N+]3=C4C=CC=CC4=CC=C3",
        "[O-]S([O-])(=O)=O.C1CN[Cr+3]23(N1)(NCCN2)NCCN3",
        "BrC1=CC=C(NC(=N)NC(=N)NC2=CC=C(C=C2)S(=O)(=O)NC3=CN=CC=N3)C=C1",
        "C([C@@H]1[C@H]([C@@H]([C@@H]([C@@H]([NH2+]1)S(=O)(=O)[O-])O)O)O)O",
        "C[C](O)([CH](C(O)=O)C1=CC=CC=C1)C2=CC=CC=C2",
        "C[C]1(CS(O)(=O)=O)[CH]2CC[C]1(C)C(=O)[CH]2Br",
    ];

    #[test]
    fn a_deg4_swap_is_exactly_isometric_and_actually_fires() {
        // 度 4 置换与前两个算子同级:轴过中心原子,反射保持每个点到中心的距离,
        // 四根键的长度都不变;被反射的两块内部是等距变换,块与块之间是同一个
        // 反射 —— 距离全保。
        //
        // 与螺环那条同样的分工:这里守"整个消冲突不改键长的多重集"与"它确实
        // 触发过",守不到"轴选错"(绕任何一条过中心的直线反射都保键长)。
        let lengths = |m: &MolBuilder, pos: &BTreeMap<u32, Point2>| -> Vec<i64> {
            let mut v: Vec<i64> = m
                .bonds()
                .iter()
                .filter_map(|b| {
                    let (u, v) = (pos.get(&b.begin)?, pos.get(&b.end)?);
                    Some((u.dist(*v) * 1e9).round() as i64)
                })
                .collect();
            v.sort_unstable();
            v
        };
        let mut fired = 0usize;
        for smi in SWAP {
            for style in &Style::ALL {
                let m = prep(smi);
                let ranks = omgkit_io::canon::canonical_ranks(&m);
                let mut pos: BTreeMap<u32, Point2> = BTreeMap::new();
                for p in layout::layout_all(&m, &ranks, style) {
                    pos.extend(p.pos);
                }
                let before = lengths(&m, &pos);
                let rep = relieve(&m, &mut pos, &ranks, style);
                let after = lengths(&m, &pos);
                if !rep.swapped.is_empty() {
                    fired += 1;
                }
                assert_eq!(
                    before, after,
                    "[{}] {smi}:消冲突改变了键长 —— 用的算子不是等距的",
                    style.name
                );
            }
        }
        assert!(
            fired > 0,
            "这批分子上度 4 置换一次都没触发 —— 那它是白加的,上面那条键长判据也是恒真的"
        );
    }

    #[test]
    fn a_deg4_swap_never_tears_the_structure_apart() {
        // **两个分支必须不相交。** 中心原子只要沾上环,从某个邻居出发的搜索就会
        // 顺着环绕回来,把另一个分支整个吞进去 —— 两块重叠,交集里的原子被反射
        // 两次(等于没动),结构当场撕开,而键长判据未必看得出来。
        //
        // RDKit 的守卫是"度 4 且完全不在环上",这里守的是同一件事,但直接查
        // 枚举出来的两块**确实不相交**,而不是相信那个代理条件。
        for smi in SWAP
            .iter()
            .chain(["C1CC2(CC1)CCCC2", "CC(C)(C)C(C)(C)C", "C1CCC2(CC1)CCCCC2"].iter())
        {
            for style in &Style::ALL {
                let m = prep(smi);
                let ranks = omgkit_io::canon::canonical_ranks(&m);
                let mut pos: BTreeMap<u32, Point2> = BTreeMap::new();
                for p in layout::layout_all(&m, &ranks, style) {
                    pos.extend(p.pos);
                }
                for sw in pair_swaps(&m, &pos, &ranks) {
                    let a: BTreeSet<u32> = sw.sides.0.iter().copied().collect();
                    let b: BTreeSet<u32> = sw.sides.1.iter().copied().collect();
                    let both: Vec<u32> = a.intersection(&b).copied().collect();
                    assert!(
                        both.is_empty(),
                        "[{}] {smi}:中心 {} 的两个分支重叠在 {both:?} 上",
                        style.name,
                        sw.centre
                    );
                    assert!(
                        !a.contains(&sw.centre) && !b.contains(&sw.centre),
                        "[{}] {smi}:分支里含了中心原子 {}",
                        style.name,
                        sw.centre
                    );
                }
            }
        }
    }

    #[test]
    fn a_spiro_flip_is_exactly_isometric_and_actually_fires() {
        // 螺环翻转与"绕可旋转键翻转"同级:**不损失任何几何性质**。轴过螺原子
        // 与它在目标环里两个邻居的中点,而两根都是环上的键、长度相同,所以那
        // 条轴恰好是角平分线,反射把两个邻居精确对调。
        //
        // 这条判据守两件事:
        //
        // - **消冲突前后键长的多重集一模一样**。守的是"只用不损失几何的算子"
        //   这条规矩 —— 哪天有人往 `relieve` 里塞一个开角或缩键的算子,这里
        //   立刻红。
        // - **螺环翻转确实翻过**。一次都不翻的话上一条恒真,判据就是空过的;
        //   第一版列了八个教科书螺环,一个都没触发。
        //
        // **它抓不到"镜像轴选错"** —— 绕任何一条过螺原子的直线反射都保持
        // 到螺原子的距离,键长怎么都不会变。轴选错由
        // `a_spiro_flip_does_not_depend_on_how_the_molecule_was_written` 抓到
        // (变异验证过:把轴改成指向 N1,那条立刻红)。
        // **判据比的是"消冲突前后键长的多重集"**,不是"每根键都等于 1"。
        // 后者在退化布局(桥环走弹簧松弛)上本来就不成立,拿它当判据会把
        // "松弛给的键长本来就不等"误报成"翻转弄坏了几何"。多重集这个说法对
        // 任何起始几何都成立,因为**等距变换保长度**。
        let lengths = |m: &MolBuilder, pos: &BTreeMap<u32, Point2>| -> Vec<i64> {
            let mut v: Vec<i64> = m
                .bonds()
                .iter()
                .filter_map(|b| {
                    let (u, v) = (pos.get(&b.begin)?, pos.get(&b.end)?);
                    Some((u.dist(*v) * 1e9).round() as i64)
                })
                .collect();
            v.sort_unstable();
            v
        };

        let mut fired = 0usize;
        for smi in SPIRO {
            for style in &Style::ALL {
                let m = prep(smi);
                let ranks = omgkit_io::canon::canonical_ranks(&m);
                let mut pos: BTreeMap<u32, Point2> = BTreeMap::new();
                for p in layout::layout_all(&m, &ranks, style) {
                    pos.extend(p.pos);
                }
                let before = lengths(&m, &pos);
                let rep = relieve(&m, &mut pos, &ranks, style);
                let after = lengths(&m, &pos);
                if !rep.spiro_flipped.is_empty() {
                    fired += 1;
                }
                assert_eq!(
                    before, after,
                    "[{}] {smi}:消冲突改变了键长 —— 用的算子不是等距的",
                    style.name
                );
            }
        }
        assert!(
            fired > 0,
            "这批螺环分子上一次都没翻过 —— 那这个算子是白加的,上面那条键长判据也就是恒真的"
        );
    }

    #[test]
    fn a_spiro_flip_does_not_depend_on_how_the_molecule_was_written() {
        // 翻哪个环、翻哪一侧,RDKit 那边两处都取存储序(`rings[0]` 与
        // `atomNeighbors` 的第一个)。照抄的话,同一个分子换种写法就会翻不同的
        // 环、得到不同的图 —— 而**写法无关是本库的头号契约**。
        //
        // 这里两处都按规范秩定。判据直接比最终坐标:同一个分子的两种写法,
        // 逐点相同。
        for smi in SPIRO {
            for style in &Style::ALL {
                let m = prep(smi);
                let n = m.num_atoms();
                let want = crate::generate(&m, style);
                // 换一种写法:把优先序整个倒过来,足以让存储序变样
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
                    continue; // 改写出来不是同一个分子,比不了
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
    fn a_flip_never_reaches_into_another_fragment() {
        // `relieve` 跑在**所有分量合并之后**的 `pos` 上(分量之间也可能撞上),
        // 而 `far_side` 取的是"原子更少的那一侧" —— 更少的那一侧算不出来时它
        // 返回 `pos` 的**补集**,那个补集会跨过分量边界。
        //
        // 真翻下去就是:动一根键,把**另一个不相干的片段**整个镜像走。
        //
        // 实测这件事在全量语料上**没有发生过**(280 个多分量的图,用分离轴定理
        // 判,没有一个分量互相穿插):这种翻转会制造碰撞,`better()` 把它挡掉了。
        // 所以它是个**潜在陷阱**,不是活着的缺陷 —— 但机制在那里,钉住它。
        for smi in [
            "[Na+].[Cl-]",
            "CCCCCCCC.c1ccccc1",
            "CC(=O)Oc1ccccc1C(=O)O.CC(C)(C)c1ccccc1",
            "[O-]S([O-])(=O)=O.CCCCCCCCCC",
            "c1ccccc1.c1ccccc1.CCCC",
        ] {
            for style in &Style::ALL {
                let m = prep(smi);
                let ranks = omgkit_io::canon::canonical_ranks(&m);
                let mut pos: BTreeMap<u32, Point2> = BTreeMap::new();
                for p in layout::layout_all(&m, &ranks, style) {
                    pos.extend(p.pos);
                }
                // 原子 → 分量号
                let n = u32::try_from(m.num_atoms()).unwrap();
                let mut comp = vec![usize::MAX; n as usize];
                let mut c = 0usize;
                for s in 0..n {
                    if comp[s as usize] != usize::MAX {
                        continue;
                    }
                    let mut st = vec![s];
                    comp[s as usize] = c;
                    while let Some(x) = st.pop() {
                        for (y, _) in m.neighbors(x) {
                            if comp[y as usize] == usize::MAX {
                                comp[y as usize] = c;
                                st.push(y);
                            }
                        }
                    }
                    c += 1;
                }
                for b in 0..u32::try_from(m.num_bonds()).unwrap() {
                    let Some(side) = far_side(&m, &pos, b) else {
                        continue;
                    };
                    let want = comp[m.bonds()[b as usize].begin as usize];
                    let spans: BTreeSet<usize> = side.iter().map(|a| comp[*a as usize]).collect();
                    assert!(
                        spans.len() == 1 && spans.contains(&want),
                        "[{}] {smi}:翻键 {b} 会动到别的片段(涉及分量 {spans:?},键在 {want})",
                        style.name
                    );
                }
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
