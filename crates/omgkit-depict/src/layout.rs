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
pub(crate) fn layout_all(
    mol: &MolBuilder,
    ranks: &[u32],
    style: &Style,
    over: crate::templates::Override<'_>,
) -> Vec<Piece> {
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
        .map(|atoms| layout_component(mol, &atoms, &systems_all, ranks, style, over))
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

/// 记一笔退化,**同时**把那一块的原子记进「离网」。
///
/// 两件事必须一起做,所以合成一个函数 —— 分开写过一版,四个记录点里漏掉任何
/// 一处,现有判据**一条都不会红**(实测)。合起来之后"漏掉一处"这个状态在结构
/// 上就表示不出来了。
///
/// 「离网」的含义与用处见 [`crate::chains::Env::off_grid`]。
fn note_degraded(
    deg: Option<Degradation>,
    atoms: impl IntoIterator<Item = u32>,
    degraded: &mut Vec<Degradation>,
    off_grid: &mut BTreeSet<u32>,
) {
    if let Some(d) = deg {
        degraded.push(d);
        off_grid.extend(atoms);
    }
}

fn layout_component(
    mol: &MolBuilder,
    atoms: &[u32],
    systems: &[rings::System<'_>],
    ranks: &[u32],
    style: &Style,
    over: crate::templates::Override<'_>,
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
    // **坐标不是从理想 30° 栅格算出来的那些原子。** 松弛出来的环系统做种子,
    // 再沿摆放传下去:从一个离网的枢纽挑方向,挑出来的位置同样离网。
    // 用处见 `chains::Env::off_grid`。
    let mut off_grid: BTreeSet<u32> = BTreeSet::new();
    let mut done_sys: BTreeSet<usize> = BTreeSet::new();

    // 摆环系统时要能判"撞没撞上",用的是消冲突那一套半径与"成键的不算"。
    // 同一套口径,免得布局觉得没撞、消冲突觉得撞了。
    let radii = crate::refine::radii(mol, style);
    let bonded: BTreeSet<(u32, u32)> = mol
        .bonds()
        .iter()
        .map(|b| (b.begin.min(b.end), b.begin.max(b.end)))
        .collect();

    // 起点
    let seed_atom = if let Some(&i) = mine.iter().max_by_key(|&&i| {
        (
            systems[i].atoms.len(),
            std::cmp::Reverse(sys_key(&systems[i], ranks)),
        )
    }) {
        let (local, deg) = rings::layout_local(mol, &systems[i], ranks, over);
        note_degraded(deg, local.keys().copied(), &mut degraded, &mut off_grid);
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
            let (local, deg) = rings::layout_local(mol, &systems[s], ranks, over);
            note_degraded(deg, local.keys().copied(), &mut degraded, &mut off_grid);
            // **求和的次序必须定死 —— 按秩,不按存储下标。**
            //
            // `pos` 是 `BTreeMap`,迭代顺序按键定序,是确定的(会随机播种的是
            // `HashMap`,见本文件模块注释);这里要挡的不是进程间的不确定,是
            // **写法依赖**:同一个分子换个 SMILES 写法,原子的存储下标就变,
            // 浮点求和的次序跟着变,`away` 的方向差一点点,图就摆得不一样。
            // 实测:全量语料的违例数在 141/142 之间来回跳。
            let mut placed: Vec<(u32, Point2)> =
                pos.iter().map(|(k, v)| (ranks[*k as usize], *v)).collect();
            placed.sort_unstable_by_key(|x| x.0);
            let away = away_from(pos[&a], placed.iter().map(|x| x.1));
            // **锚点自己就是接口,方向是自由的** —— 所以整圈都能试。
            // `away` 排第一,但这**不是**"从前的行为原样保留":参照物同时
            // 从质心改成了平分线,两笔改动一起动了这条通路。
            let dirs = ring_of_dirs(away);
            let around = Around {
                pos: &pos,
                radii: &radii,
                bonded: &bonded,
            };
            let put = place_clear(mol, &local, a, pos[&a], &dirs, &around);
            for (k, p) in put {
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

        // **挑方向之前先把待摆的那几块算出来。**
        //
        // `place_neighbours` 从前只看"这一个原子落在这儿撞不撞",而这些邻居
        // 里有的是一整块环系统的接口。要让它把整块算进去,就得先有那一块的
        // 局部坐标 —— 而 `layout_local` 很贵,所以算一次、既给前瞻用、也给
        // 真正摆放用(见 `Placed::block`),**不是多算一遍**。
        //
        // 一个原子理论上可以同时是几个未摆系统的接口(螺原子挂在链上)。这里
        // 只规划**规范序里的第一个**,其余仍走下面的老路 —— 那种情形语料里
        // 没出现过,不为它加复杂度。
        /// 一个待摆的环系统:`(系统下标, 局部坐标, 要记的退化)`。
        type Plan = (usize, crate::chains::Block, Option<Degradation>);
        let mut plans: BTreeMap<u32, Plan> = BTreeMap::new();
        // **同一个系统只规划一次。** 两个待放邻居可能落在**同一个**未摆系统里
        // —— 配位键被环感知排除在外(`omgkit_chem::rings`),而它算进连通性,
        // 所以一个金属可以用两根配位键咬住同一个稠环系的两个给体原子。
        //
        // 不去重的话,第二个给体也会拿到一份"块",而那一块**根本不会被画**
        // (第一个已经把整个系统摆上了)。它却会被 `place_neighbours` 累积进
        // "已占",后面的兄弟就是对着一个不存在的环挑方向;退化也会被记两遍。
        //
        // 语料里没出现(全量 8831 个分子只有 1 条配位键),但构造得出来,所以
        // 这里挡住。
        let mut planned: BTreeSet<usize> = BTreeSet::new();
        for &b in &todo {
            let Some(&s) = of_atom
                .get(&b)
                .into_iter()
                .flatten()
                .find(|s| !done_sys.contains(s) && !planned.contains(s))
            else {
                continue;
            };
            let (local, deg) = rings::layout_local(mol, &systems[s], ranks, over);
            planned.insert(s);
            plans.insert(b, (s, local, deg));
        }
        let blocks: BTreeMap<u32, crate::chains::Block> = plans
            .iter()
            .map(|(b, (_, local, _))| (*b, local.clone()))
            .collect();
        let env = crate::chains::Env {
            mol,
            ranks,
            style,
            radii: &radii,
            bonded: &bonded,
            off_grid: &off_grid,
            blocks: &blocks,
        };

        for p in place_neighbours(&env, a, &pos, &todo, z) {
            // **兄弟那一块可能已经把它摆上了。** 上面去重之后,落在同一个系统里
            // 的第二个邻居没有自己的块,但第一个邻居的块会把整个系统(含它)摆
            // 好。这时它的坐标由环几何定,链上分给它的那个方向作废 —— 直接盖
            // 上去会把它从环里拽出来。它也已经被那一块的 `fresh` 进过队了。
            if pos.contains_key(&p.atom) {
                continue;
            }
            pos.insert(p.atom, p.at);
            zig.insert(p.atom, p.zig);

            // 前瞻已经把这一块摆好了 —— **原样用**。重算会重新挑镜像,而那时
            // `pos` 已经变了,挑出来的可能是另一块,与前瞻累积进去的对不上。
            if let (Some(put), Some((s, _, deg))) = (p.block, plans.remove(&p.atom)) {
                // 退化只在**真正摆下去**这一刻记一笔 —— 规划阶段不记,否则同一
                // 个系统会被记两遍
                note_degraded(deg, put.keys().copied(), &mut degraded, &mut off_grid);
                for (k, q) in put {
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

            // 同一个原子上还挂着别的未摆系统时走老路(规划阶段只规划了第一个)
            let sys_here: Vec<usize> = of_atom
                .get(&p.atom)
                .into_iter()
                .flatten()
                .copied()
                .filter(|s| !done_sys.contains(s))
                .collect();
            for s in sys_here {
                let (local, deg) = rings::layout_local(mol, &systems[s], ranks, over);
                note_degraded(deg, local.keys().copied(), &mut degraded, &mut off_grid);
                // **这一路方向不自由**:环外键已经画在那儿了,平分线必须落在它
                // 的延长线上,否则两根环键就不对称了(那正是 79.1° 那一族的
                // 来历,见 `rings::place_candidates`)。能挑的只有镜像。
                let dir = (p.at - pos[&a]).normalized();
                let around = Around {
                    pos: &pos,
                    radii: &radii,
                    bonded: &bonded,
                };
                let put = place_clear(mol, &local, p.atom, p.at, &[dir], &around);
                for (k, q) in put {
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
///
/// 它只是 [`place_clear`] 的**头号候选**,不是最终答案 —— 下面那条退化分支
/// 从前是致命的,现在只是"第一发打偏了"。
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
        // **已放置的部分以 `from` 为质心时,这个方向是瞎给的。**
        //
        // 注释从前写着"任何方向都一样" —— 那是错的。质心为零只说明已放的东西
        // 关于 `from` 大致对称,不说明它们均匀铺满;固定取 (1,0) 会正对着某一块
        // 已经画好的东西。三条一模一样的螯合环就是这么塌的:第一条摆好,第二条
        // 摆到对面(质心归零),第三条拿到 (1,0),**和第一条逐位重合**。
        //
        // 全量语料 79 张图有原子画在同一点上,重合的 291 对里 **274 对是两个不同
        // 环系统**,这一处是主犯。现在它只是候选之一,撞了会被换掉。
        Point2::new(1.0, 0.0)
    } else {
        (mean * -1.0).normalized()
    }
}

/// 以 `first` 为首、按 30° 一档往两边铺开的一整圈方向。
///
/// 30° 是全库统一的栅格(取代基避让也走它,见 `chains::free_direction`),环系统
/// 没有理由用别的。
///
/// `first` 排在最前,但**这不等于"它够用时结果不变"**:[`place_clear`] 只在
/// 「不撞**且**不交叉」时才提前收工,`first` 不撞却带一处交叉的话,后面的方向
/// 仍可能凭交叉更少胜出。
fn ring_of_dirs(first: Point2) -> Vec<Point2> {
    let step = std::f64::consts::FRAC_PI_6;
    let mut out = vec![first];
    for k in 1..=6i32 {
        for sign in [1.0, -1.0] {
            if k == 6 && sign < 0.0 {
                continue; // +180° 与 −180° 是同一个方向
            }
            out.push(first.rotated(step * f64::from(k) * sign));
        }
    }
    out
}

/// 把一个环系统摆到锚点上,在「候选方向 × 两个镜像」里挑**最不撞**的那个。
///
/// # 为什么非挑不可
///
/// 先前是算一个方向就直接摆,**摆完不看撞没撞**。语料上的后果有两档:
///
/// | | |
/// |---|---:|
/// | 原子画在同一点上(重合对) | 291,其中**两个不同环系统 274**(94%) |
/// | 「键角不过窄」剩下的违例 | 25,其中 79.1° 那一族 **16** |
///
/// 后者是参照物选错(质心而非外角平分线),已在 [`rings::place_candidates`]
/// 里修掉;前者是**根本没找过别的位置**,这个函数补的就是这一步。
///
/// 那 291 对重合的距离**全是 0.0000**,不是"挤了一点"—— 键长定死为 1、键角走
/// 30° 栅格,所有原子都落在同一张格点上,两块折回来就是逐位重合。(这说的是
/// 改动**之前**;改完之后剩下的那些里已经有距离不为零的了。)
///
/// # 这一步只补得上一半:第二个调用点几乎没得挑
///
/// 锚点自由的那一路(螺环 / 割点)有整整一圈方向可试,而**从一根键挂上去**的
/// 那一路只剩镜像,单环上那个镜像还是空的([`rings::place_candidates`] 里说了
/// 理由:正多边形关于平分线对称)。所以"两个吡啶配体逐位重合"这类图,这个
/// 函数**够不着** —— 要挪只能回到上游,让 `chains` 给锚点换一个方向,而那要
/// 先知道后面会挂一整个环。记在这儿,没做。
///
/// # 打分:深度在前,交叉破平局 —— 与 `refine` 相反,而且是量出来的
///
/// 五种排法各跑了一遍全量语料(17662 个分子×规范):
///
/// | 打分 | 原子不重合 | 有键交叉 | 取代基挤压 | 未解冲突 |
/// |---|---:|---:|---:|---:|
/// | 基线:质心参照,**不选位** | 79 | **48** | 28 | 1130 |
/// | 平分线参照,不选位 | 79 | 52 | 28 | 1135 |
/// | **深度在前,交叉破平局**(取它) | **57** | 72 | **12** | 1120 |
/// | 只看深度 | **57** | 73 | **12** | 1119 |
/// | 交叉在前 | 77 | 54 | 22 | 1120 |
/// | `(撞上的对数 + 交叉数, 深度)` | 79 | 58 | 28 | 1120 |
/// | `(撞上的对数, 交叉数, 深度)` | 72 | 63 | 33 | 1120 |
///
/// **重合与交叉几乎是 1:1 的取舍**,五种排法只是在同一条前沿上滑动,没有一种
/// 能把两者一起压下去。选位买到的每一处重合,大约要还回一处交叉。
///
/// 取深度在前,理由是两者**不同级**:「原子不重合」是判据里的**硬性质**,
/// 两个原子叠在一点会让图上多出一个分子里没有的环,而读者没有任何办法看出
/// 那个环是假的;「键交叉」只是质量分档 —— 交叉是**看得见**的,读者知道那里
/// 有两根键。宁可多一处看得见的丑,不要多一处看不见的假。
///
/// `refine::score` 把交叉排在深度**前面**,那也是量出来的,不矛盾:
/// 消冲突动的是已经摆好的图,那时深度差一点只是"挤";这边是**在决定摆哪儿**,
/// 深度大到极限就是重合。
///
/// 交叉只做破平局 —— 一大堆候选深度都是 0,谁在前面就选谁,那个次序本来是
/// 任意的,让交叉来定至少省下一处(73 → 72)。
///
/// # 平局怎么破
///
/// 逐个跟"当前最好"比:深度小 1e-9 以上就换,差在容差内就比交叉数,都平就
/// 留着前面那个 —— 而候选序是写法无关的(`dirs` 由 `away` 或键方向生成,
/// 两者都已经定死)。深度本身也**先排序再求和**,免得同一份几何在两种写法下
/// 差最后一位。
///
/// **这不完全等于按 `(深度, 交叉数)` 排序**:带容差的"相等"不传递,深度
/// `(0, 0.9e-9, 1.7e-9)` 配交叉 `(5, 0, 0)` 时,正着扫和倒着扫会选中不同的
/// 候选。它**不破坏写法无关**(候选序是定死的),只是结果不等于一个真正的
/// 字典序。现实里深度要么恰好 0、要么远大于 1e-9,碰不到;记在这儿,免得
/// 以后照字典序去推它的行为。
fn place_clear(
    mol: &MolBuilder,
    local: &BTreeMap<u32, Point2>,
    anchor: u32,
    at: Point2,
    dirs: &[Point2],
    around: &Around<'_>,
) -> BTreeMap<u32, Point2> {
    const EPS: f64 = 1e-9;
    // 已经画出来的键 —— 每个候选都要跟它们比,只算一次
    let drawn: Vec<(Point2, Point2)> = mol
        .bonds()
        .iter()
        .filter_map(|b| Some((*around.pos.get(&b.begin)?, *around.pos.get(&b.end)?)))
        .collect();
    let mut best: Option<((f64, usize), BTreeMap<u32, Point2>)> = None;
    for d in dirs {
        for cand in rings::place_candidates(mol, local, anchor, at, *d) {
            let s = (
                clash(&cand, around),
                new_crossings(mol, &cand, around.pos, &drawn),
            );
            let better = match &best {
                None => true,
                Some((old, _)) => {
                    if s.0 < old.0 - EPS {
                        true
                    } else if s.0 > old.0 + EPS {
                        false
                    } else {
                        s.1 < old.1
                    }
                }
            };
            if better {
                let done = s.0 == 0.0 && s.1 == 0;
                best = Some((s, cand));
                if done {
                    return best.expect("刚放进去").1; // 不撞不交叉,不用再找了
                }
            }
        }
    }
    best.expect("`dirs` 非空,候选至少有两个").1
}

/// 已经画在纸上的那部分,连同判"撞没撞上"要用的两样东西。
struct Around<'a> {
    pos: &'a BTreeMap<u32, Point2>,
    /// 每个原子的碰撞半径,来自 [`crate::refine::radii`]
    radii: &'a [f64],
    /// 成键的原子对 —— 它们本来就靠在一起,不算撞
    bonded: &'a BTreeSet<(u32, u32)>,
}

/// 一块待放置的坐标与**已放置部分**的碰撞深度平方和。
///
/// 口径与 [`crate::refine`] 一致:半径来自标签,成键的一对不算(相邻原子本来
/// 就靠在一起)。已经在 `pos` 里的原子(锚点自己)不参与 —— 它不是新放的。
fn clash(cand: &BTreeMap<u32, Point2>, around: &Around<'_>) -> f64 {
    let (pos, radii, bonded) = (around.pos, around.radii, around.bonded);
    let mut parts: Vec<f64> = Vec::new();
    for (i, p) in cand {
        if pos.contains_key(i) {
            continue;
        }
        for (j, q) in pos {
            if bonded.contains(&((*i).min(*j), (*i).max(*j))) {
                continue;
            }
            let want = radii[*i as usize] + radii[*j as usize];
            let d = p.dist(*q);
            if d < want {
                parts.push((want - d).powi(2));
            }
        }
    }
    // **先排序再求和。** 两种写法给的是同一个几何、同一个多重集,但迭代序是
    // 存储序 —— 不排的话和会差最后一位,而上面的平局判定就靠这一位。
    parts.sort_by(f64::total_cmp);
    parts.iter().sum()
}

/// 这块坐标摆上去会新添几处键交叉。
///
/// 只数「新系统的键 × 已画出的键」这一类。系统**内部**的交叉与摆在哪儿无关
/// (那是 [`rings::layout_local`] 的事),数了也是每个候选一样多。
fn new_crossings(
    mol: &MolBuilder,
    cand: &BTreeMap<u32, Point2>,
    pos: &BTreeMap<u32, Point2>,
    drawn: &[(Point2, Point2)],
) -> usize {
    mol.bonds()
        .iter()
        .filter(|b| !(pos.contains_key(&b.begin) && pos.contains_key(&b.end)))
        .filter_map(|b| Some((*cand.get(&b.begin)?, *cand.get(&b.end)?)))
        .map(|(u, v)| {
            drawn
                .iter()
                .filter(|(x, y)| crate::geom::segments_cross(u, v, *x, *y))
                .count()
        })
        .sum()
}

#[cfg(test)]
mod tests {

    /// 同一个原子上挂着的几个环系统,不许有两个被摆到同一个位置上。
    ///
    /// 分子取自真语料(`harness/corpus/large.smi` 第 2536 行):三(乙二胺)合钴。
    /// 钴是割点,三条螯合环是**三个各含它一个**的环系统 —— 它们走
    /// `layout_component` 里"以锚点自己为锚"那一路。
    ///
    /// 先前那一路只算一个 `away_from` 方向就直接摆下去,**摆完不看撞没撞**:
    /// 第一条摆好,第二条摆到对面,此时已放置部分关于钴的质心**归零**,第三条
    /// 于是拿到退化兜底的固定方向 `(1, 0)` —— 与第一条**逐位重合**。
    ///
    /// 距离正好 0.0000 不是巧合:键长定死为 1、键角走 30° 栅格,所有原子都落在
    /// 同一张格点上,折回来就是精确重合。图上于是多出一个分子里没有的环,而
    /// 读者没有任何办法看出那个环是假的。
    ///
    /// 变异:把 `place_clear` 的候选截成第一个(`dirs.iter().take(1)` 且
    /// `place_candidates(..).into_iter().take(1)`)→ 这条当场红。
    #[test]
    fn three_chelate_rings_on_one_metal_do_not_land_on_top_of_each_other() {
        use omgkit_chem::rings::fused_ring_systems;

        let smi = "C1CN[Co]23(N1)(NCCN2)NCCN3";
        let mut m = omgkit_io::smiles::parse(smi).expect("SMILES 该能解析");
        omgkit_chem::pipeline::sanitize(&mut m).expect("该能 sanitize");

        // **前提要自己成立**:必须真有三个环系统共用同一个原子,否则这条判据
        // 查的根本不是它想查的那一路。
        let systems = fused_ring_systems(&m);
        let shared = (0..u32::try_from(m.num_atoms()).expect("原子数超出 u32"))
            .map(|a| systems.iter().filter(|s| s.contains(&a)).count())
            .max()
            .unwrap_or(0);
        assert!(
            systems.len() >= 3 && shared >= 3,
            "{smi} 只有 {} 个环系、最多共用 {shared} 个 —— 走不到「以锚点自己为锚」那一路",
            systems.len()
        );

        for style in &crate::style::Style::ALL {
            let d = crate::generate(&m, style);
            for i in 0..d.coords.len() {
                for j in (i + 1)..d.coords.len() {
                    let dist = d.coords[i].dist(d.coords[j]);
                    assert!(
                        dist >= 0.05,
                        "{}:原子 {i} 与 {j} 相距 {dist:.4} 个键长",
                        style.name
                    );
                }
            }
        }
    }

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
        let ps = layout_all(&m, &ranks, &Style::ACS_1996, None);
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
