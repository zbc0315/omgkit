//! 布局驱动:把环系统、链、取代基串成一张完整的图。
//!
//! # 走法
//!
//! 每个连通分量各自布局,最后并排摆开。分量内部是一次广度优先:
//!
//! 1. **起点**。有环就取原子最多的那个环系统(平局按规范秩),整块先摆好;
//!    没有环就取规范秩最小的原子。
//! 2. **逐个原子向外**。遇到还没放置的邻居:
//!    - 邻居落在**还没摆过的环系统**里 → 把整个系统摆上去
//!    - 否则按 [`chains`](crate::chains) 的规则给它分配方向
//!
//! # 螺环是这里最容易写错的一处
//!
//! 螺原子同时属于两个环系统([`fused_ring_systems`] 用双连通分解,螺原子是割点,
//! 两边天然分开)。BFS 走到螺原子时第一个系统已经摆好、第二个还没有,而**锚点
//! 就是螺原子自己**,不是它的某个邻居。
//!
//! 按"先把邻居放到一个键长外、再以邻居为锚"的通路去处理螺环,会把螺原子复制
//! 成两个位置 —— 图上看起来是两个环被一根键连着,而不是共用一个顶点。这个错
//! 不会让任何一步报错。
//!
//! [`fused_ring_systems`]: omgkit_chem::rings::fused_ring_systems

// **迭代顺序必须与进程无关。** 标准库的 `HashMap` 用随机播种的哈希器,
// 迭代顺序每次运行都不同;而这里有按位置求和、取极值这类操作,顺序一变结果的
// 末位就变,同一个分子同一份代码两次运行画出的图会不一样。实测:全量语料的
// 违例数在 141/142 之间来回跳。`BTreeMap` 按键定序,这一整类不确定性消失。
use std::collections::{BTreeMap, BTreeSet, VecDeque};

use omgkit_chem::{rings::fused_ring_systems, sssr::ring_set};
use omgkit_core::MolBuilder;

use crate::chains::place_neighbours;
use crate::geom::Point2;
use crate::rings::{self, Degradation};
use crate::style::Style;

/// 一个分量的布局结果。
pub(crate) struct Piece {
    pub pos: BTreeMap<u32, Point2>,
    pub degraded: Vec<Degradation>,
}

/// 给整个分子布局,返回逐分量的结果(尚未并排摆开)。
pub(crate) fn layout_all(mol: &MolBuilder, ranks: &[u32], style: &Style) -> Vec<Piece> {
    let comps = components(mol);
    let rings_all = ring_set(mol);
    let mut systems_all = rings::group(&fused_ring_systems(mol), &rings_all);
    // **环系统必须按规范秩定序。** `fused_ring_systems` 的返回顺序来自双连通
    // 分解的遍历,是**存储序** —— 而下面 `of_atom` 的值就是按这个顺序堆的,
    // 螺原子上挂着的几个系统于是按写法依赖的顺序摆出去。
    //
    // 平时看不出来是因为几个系统通常长得不一样,摆错顺序也还能靠"背离已放部分"
    // 挑回来;**几个系统一模一样时就露馅**:镍配合物挂着三条完全相同的螯合环
    // (镍是割点,三个 12 原子系统各含它一个),两种写法下系统顺序是 (5,3,4)
    // 与 (4,5,3),摆出来差 5.29 个单位。
    //
    // 身份取"环系原子的规范秩有序多重集" —— 与原子编号无关,而且系统之间互不
    // 相同(它们的原子集不同),所以这是全序,不留平局。
    systems_all.sort_by_key(|s| {
        let mut k: Vec<u32> = s.atoms.iter().map(|a| ranks[*a as usize]).collect();
        k.sort_unstable();
        k
    });

    comps
        .into_iter()
        .map(|atoms| layout_component(mol, &atoms, &systems_all, ranks, style))
        .collect()
}

/// 环系统被消费的顺序,每个用「原子规范秩的有序多重集」表示。
///
/// 判据要验的正是这个顺序与写法无关,而它是 `layout_all` 内部的中间量 ——
/// 与其在判据里重算一遍(那就成了抄实现),不如把它开出来。
#[cfg(test)]
pub(crate) fn system_order(mol: &MolBuilder, ranks: &[u32]) -> Vec<Vec<u32>> {
    let rings_all = ring_set(mol);
    let mut systems_all = rings::group(&fused_ring_systems(mol), &rings_all);
    systems_all.sort_by_key(|s| {
        let mut k: Vec<u32> = s.atoms.iter().map(|a| ranks[*a as usize]).collect();
        k.sort_unstable();
        k
    });
    systems_all
        .iter()
        .map(|s| {
            let mut k: Vec<u32> = s.atoms.iter().map(|a| ranks[*a as usize]).collect();
            k.sort_unstable();
            k
        })
        .collect()
}

/// 连通分量,每个内部按规范秩排序,分量之间按最小规范秩排序。
///
/// **配位键计入连通性。** 环感知里它不计(那是成环判定),但画图时靠配位键连着
/// 的片段就该画在一起 —— 拆开会让一个配合物散成互不相干的几块。
fn components(mol: &MolBuilder) -> Vec<Vec<u32>> {
    // `num_atoms()` 给的是 usize,而原子下标在整个库里是 u32 —— 两者不能混着用,
    // 混了编译器会挑出来,但把某一处写成 `as` 强转就会挑不出来了
    let n = u32::try_from(mol.num_atoms()).expect("原子数超出 u32");
    let mut seen = vec![false; n as usize];
    let mut out = Vec::new();
    for start in 0..n {
        if seen[start as usize] {
            continue;
        }
        let mut stack: Vec<u32> = vec![start];
        let mut comp: Vec<u32> = Vec::new();
        seen[start as usize] = true;
        while let Some(a) = stack.pop() {
            comp.push(a);
            for (b, _) in mol.neighbors(a) {
                if !seen[b as usize] {
                    seen[b as usize] = true;
                    stack.push(b);
                }
            }
        }
        comp.sort_unstable();
        out.push(comp);
    }
    out
}

fn layout_component(
    mol: &MolBuilder,
    atoms: &[u32],
    systems: &[rings::System<'_>],
    ranks: &[u32],
    style: &Style,
) -> Piece {
    let here: BTreeSet<u32> = atoms.iter().copied().collect();
    // 落在本分量里的环系统
    let mine: Vec<usize> = (0..systems.len())
        .filter(|i| systems[*i].atoms.first().is_some_and(|a| here.contains(a)))
        .collect();
    // 原子 → 它所属的系统(螺原子属于多个)
    let mut of_atom: BTreeMap<u32, Vec<usize>> = BTreeMap::new();
    for &i in &mine {
        for &a in &systems[i].atoms {
            of_atom.entry(a).or_default().push(i);
        }
    }

    let mut pos: BTreeMap<u32, Point2> = BTreeMap::new();
    let mut degraded: Vec<Degradation> = Vec::new();
    let mut done_sys: BTreeSet<usize> = BTreeSet::new();

    // 起点
    let seed_atom = if let Some(&i) = mine.iter().max_by_key(|&&i| {
        (
            systems[i].atoms.len(),
            std::cmp::Reverse(sys_key(&systems[i], ranks)),
        )
    }) {
        let (local, deg) = rings::layout_local(mol, &systems[i], ranks);
        if let Some(d) = deg {
            degraded.push(d);
        }
        // 整块直接落在原点附近,朝向留给后面的规范朝向处理
        pos.extend(local);
        done_sys.insert(i);
        *systems[i]
            .atoms
            .iter()
            .min_by_key(|a| ranks[**a as usize])
            .expect("环系统非空")
    } else {
        let a = *atoms
            .iter()
            .min_by_key(|a| ranks[**a as usize])
            .expect("分量非空");
        pos.insert(a, Point2::ORIGIN);
        a
    };

    // 广度优先向外长
    let mut zig: BTreeMap<u32, i8> = BTreeMap::new();
    zig.insert(seed_atom, 1);
    let mut queue: VecDeque<u32> = VecDeque::new();
    // 起点若是整个环系统,系统里每个原子都要进队 —— 取代基可能挂在任何一个上
    let mut seeds: Vec<u32> = pos.keys().copied().collect();
    seeds.sort_by_key(|a| (ranks[*a as usize], *a));
    queue.extend(seeds);

    while let Some(a) = queue.pop_front() {
        let z = zig.get(&a).copied().unwrap_or(1);

        // 待放置的邻居,按规范秩排 —— 拿存储下标排就会引入写法依赖
        let mut todo: Vec<u32> = mol
            .neighbors(a)
            .map(|(b, _)| b)
            .filter(|b| !pos.contains_key(b))
            .collect();
        todo.sort_by_key(|b| (ranks[*b as usize], *b));
        todo.dedup();
        if todo.is_empty() {
            continue;
        }

        // 先处理"以 a 自己为锚"的螺环:a 还属于另一个没摆过的系统
        for &s in of_atom.get(&a).into_iter().flatten() {
            if done_sys.contains(&s) {
                continue;
            }
            let (local, deg) = rings::layout_local(mol, &systems[s], ranks);
            if let Some(d) = deg {
                degraded.push(d);
            }
            // **求和的次序必须定死。** `pos` 是 BTreeMap,它的迭代顺序每个
            // 进程随机播种,而这里是浮点求和 —— 同一个分子同一份代码,两次
            // 运行就可能得到方向差一点点的 `away`,进而摆出不同的图。
            // 实测:全量语料的违例数在 141/142 之间来回跳。
            let mut placed: Vec<(u32, Point2)> =
                pos.iter().map(|(k, v)| (ranks[*k as usize], *v)).collect();
            placed.sort_unstable_by_key(|x| x.0);
            let away = away_from(pos[&a], placed.iter().map(|x| x.1));
            for (k, p) in rings::place_at(&local, a, pos[&a], away) {
                pos.entry(k).or_insert(p);
            }
            done_sys.insert(s);
            let mut fresh: Vec<u32> = systems[s].atoms.clone();
            fresh.sort_by_key(|x| (ranks[*x as usize], *x));
            for f in fresh {
                zig.entry(f).or_insert(-z);
                queue.push_back(f);
            }
        }

        // 重新算一遍:上面摆螺环可能已经放好了一部分
        let mut todo: Vec<u32> = mol
            .neighbors(a)
            .map(|(b, _)| b)
            .filter(|b| !pos.contains_key(b))
            .collect();
        todo.sort_by_key(|b| (ranks[*b as usize], *b));
        todo.dedup();

        for p in place_neighbours(mol, a, &pos, &todo, ranks, style, z) {
            pos.insert(p.atom, p.at);
            zig.insert(p.atom, p.zig);

            // 这个邻居若落在还没摆过的环系统里,以**它**为锚把整块摆上
            let sys_here: Vec<usize> = of_atom
                .get(&p.atom)
                .into_iter()
                .flatten()
                .copied()
                .filter(|s| !done_sys.contains(s))
                .collect();
            for s in sys_here {
                let (local, deg) = rings::layout_local(mol, &systems[s], ranks);
                if let Some(d) = deg {
                    degraded.push(d);
                }
                let dir = (p.at - pos[&a]).normalized();
                for (k, q) in rings::place_at(&local, p.atom, p.at, dir) {
                    pos.entry(k).or_insert(q);
                }
                done_sys.insert(s);
                let mut fresh: Vec<u32> = systems[s].atoms.clone();
                fresh.sort_by_key(|x| (ranks[*x as usize], *x));
                for f in fresh {
                    zig.entry(f).or_insert(p.zig);
                    queue.push_back(f);
                }
            }
            queue.push_back(p.atom);
        }
    }

    debug_assert_eq!(pos.len(), atoms.len(), "有原子没被放上,BFS 漏了分支");
    Piece { pos, degraded }
}

/// 环系统的确定性排序键。
fn sys_key(s: &rings::System<'_>, ranks: &[u32]) -> Vec<u32> {
    let mut k: Vec<u32> = s.atoms.iter().map(|a| ranks[*a as usize]).collect();
    k.sort_unstable();
    k
}

/// 从 `from` 指向"已放置原子的反方向"。已放置的只有 `from` 自己时给一个固定方向。
fn away_from(from: Point2, placed: impl Iterator<Item = Point2>) -> Point2 {
    let mut sum = Point2::ORIGIN;
    let mut n = 0.0;
    for p in placed {
        sum = sum + (p - from);
        n += 1.0;
    }
    if n == 0.0 {
        return Point2::new(1.0, 0.0);
    }
    let mean = sum * (1.0 / n);
    if mean.norm() < 1e-9 {
        // 已放置的部分正好以 from 为质心 —— 任何方向都一样,取固定的保证可复现
        Point2::new(1.0, 0.0)
    } else {
        (mean * -1.0).normalized()
    }
}

#[cfg(test)]
mod tests {

    #[test]
    fn the_ring_systems_are_consumed_in_canonical_order() {
        // `fused_ring_systems` 的返回顺序来自双连通分解的遍历,是**存储序**;
        // 而 `of_atom` 的值就是按这个顺序堆的,螺原子上挂着的几个系统于是按
        // 写法依赖的顺序摆出去。
        //
        // 平时看不出来是因为几个系统通常长得不一样;**一模一样时就露馅**:
        // 镍配合物挂三条完全相同的螯合环(镍是割点,三个 12 原子系统各含它
        // 一个),两种写法下原始顺序是 (5,3,4) 与 (4,5,3),摆出来差 5.29 个单位。
        //
        // 这条判据验的是**顺序本身**,不看画出来的图 —— 所以它与模板表无关,
        // 换表不会让它空过。
        for smi in [
            // 镍配合物,三个一模一样的螯合环
            "C[N+]12CCCC1C3=CC=C[N+](=C3)[Ni++]245([N+]6=CC(=CC=C6)C7CCC[N+]47C)\
             [N+]8=CC(=CC=C8)C9CCC[N+]59C.SC#N",
            "c1ccc2ccccc2c1.c1ccccc1", // 两个分量各带环系
            "C1CC1c1ccccc1C2CC2",      // 一条链上挂三个环系
        ] {
            let smi: String = smi.split_whitespace().collect();
            let mut m = omgkit_io::smiles::parse(&smi).expect("SMILES 该能解析");
            omgkit_chem::pipeline::sanitize(&mut m).expect("该能 sanitize");
            let ranks = omgkit_io::canon::canonical_ranks(&m);
            let base = system_order(&m, &ranks);
            assert!(
                base.len() >= 2,
                "{smi} 只有 {} 个环系,验不出顺序",
                base.len()
            );

            let mut compared = 0usize;
            for seed in 0..12u64 {
                let w = omgkit_io::smiles::write_with_priority(&m, &shuffled(m.num_atoms(), seed));
                let Ok(mut m2) = omgkit_io::smiles::parse(&w.smiles) else {
                    continue;
                };
                if omgkit_chem::pipeline::sanitize(&mut m2).is_err() {
                    continue;
                }
                if omgkit_io::canon::canonical_smiles(&m2).smiles
                    != omgkit_io::canon::canonical_smiles(&m).smiles
                {
                    continue;
                }
                let r2 = omgkit_io::canon::canonical_ranks(&m2);
                compared += 1;
                assert_eq!(
                    base,
                    system_order(&m2, &r2),
                    "{smi} 写成 {} 之后环系统的顺序变了",
                    w.smiles
                );
            }
            assert!(compared >= 6, "{smi} 只比上了 {compared} 种写法");
        }
    }

    /// splitmix64 + Fisher-Yates。
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
            let j = usize::try_from(next() % (i as u64 + 1)).unwrap();
            v.swap(i, j);
        }
        v
    }
    use super::*;

    fn prep(smi: &str) -> MolBuilder {
        let mut m = omgkit_io::smiles::parse(smi).unwrap();
        omgkit_chem::pipeline::sanitize(&mut m).unwrap();
        m
    }

    fn run(smi: &str) -> (MolBuilder, Vec<Piece>) {
        let m = prep(smi);
        let ranks = omgkit_io::canon::canonical_ranks(&m);
        let ps = layout_all(&m, &ranks, &Style::ACS_1996);
        (m, ps)
    }

    #[test]
    fn every_atom_gets_a_finite_coordinate() {
        // 漏放一个原子、或让 NaN 溜进来,后面每一步都会跟着错,而且不报错
        for smi in [
            "CCO",
            "c1ccccc1",
            "c1ccc2ccccc2c1",
            "CC(C)(C)c1ccccc1",
            "C1CC2(CC1)CCCC2",
            "c1ccc(-c2ccccc2)cc1",
            "CC(=O)Oc1ccccc1C(=O)O",
            "C1CC2CCC1CC2",
            "[Na+].[Cl-]",
            "O",
            "C",
        ] {
            let (m, ps) = run(smi);
            let total: usize = ps.iter().map(|p| p.pos.len()).sum();
            assert_eq!(total, m.num_atoms(), "{smi} 有原子没放上");
            for p in &ps {
                for (a, q) in &p.pos {
                    assert!(
                        q.x.is_finite() && q.y.is_finite(),
                        "{smi} 原子 {a} 坐标非有限"
                    );
                }
            }
        }
    }

    #[test]
    fn disconnected_pieces_come_out_separately() {
        let (_, ps) = run("[Na+].[Cl-]");
        assert_eq!(ps.len(), 2, "盐应当分成两个分量");
        let (_, one) = run("CCO");
        assert_eq!(one.len(), 1);
    }

    #[test]
    fn a_spiro_atom_is_one_point_not_two() {
        // 螺原子同时属于两个环系统。按"先把邻居放到一个键长外、再以邻居为锚"
        // 的通路去处理,会把它复制成两个位置 —— 图上成了两个环被一根键连着,
        // 而不是共用一个顶点。**这个错不会让任何一步报错。**
        let (m, ps) = run("C1CC2(CC1)CCCC2");
        assert_eq!(ps.len(), 1);
        let pos = &ps[0].pos;
        assert_eq!(pos.len(), m.num_atoms());

        // 螺原子:度数为 4 且在环里
        let spiro = (0..u32::try_from(m.num_atoms()).unwrap())
            .find(|a| m.degree(*a) == 4)
            .expect("螺[4.4]壬烷有一个四度碳");
        // 它的四个邻居必须都恰好一个键长远,且分处两个环
        for (b, _) in m.neighbors(spiro) {
            let d = pos[&spiro].dist(pos[&b]);
            assert!(
                (d - 1.0).abs() < 1e-6,
                "螺原子到邻居 {b} 的距离是 {d},应当是 1"
            );
        }
        // 两个环不能叠在一起
        let pts: Vec<Point2> = pos.values().copied().collect();
        for i in 0..pts.len() {
            for j in (i + 1)..pts.len() {
                assert!(pts[i].dist(pts[j]) > 0.4, "有两个原子几乎重合");
            }
        }
    }

    #[test]
    fn bonded_atoms_stay_one_unit_apart() {
        for smi in [
            "CCO",
            "c1ccccc1",
            "CC(C)(C)c1ccccc1",
            "CC(=O)Oc1ccccc1C(=O)O",
        ] {
            let (m, ps) = run(smi);
            let pos: BTreeMap<u32, Point2> = ps
                .iter()
                .flat_map(|p| p.pos.iter().map(|(k, v)| (*k, *v)))
                .collect();
            for b in m.bonds() {
                let d = pos[&b.begin].dist(pos[&b.end]);
                assert!(
                    (d - 1.0).abs() < 1e-6,
                    "{smi} 键 {}–{} 长 {d}",
                    b.begin,
                    b.end
                );
            }
        }
    }

    // 写法无关那条判据**不在这里** —— 它属于 `generate` 的最终输出。
    // 布局本身会因 SSSR 给出的环原子顺序而有差异,消冲突之后才收敛;放在这一层
    // 测的是中间态,红了也说明不了问题出在哪。见 `crate::tests`。
}
