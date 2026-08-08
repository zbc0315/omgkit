//! 环系统布局 —— 整个算法里最难的一段。
//!
//! # 分三档,而且**第三档如实承认自己是退化解**
//!
//! | 形状 | 做法 |
//! |---|---|
//! | 单环 | 正多边形 |
//! | 邻稠(逐个环只与已放置部分共用**一根键**) | 沿那根键把新多边形拼到外侧 |
//! | 桥环 / 笼状 | 规则给不出好解 —— 退化到弹簧松弛,并记进 [`Degradation`] |
//!
//! 第三档是所有工具箱共同的软肋(见 Mayfield, RDKit UGM 2016:桥环与拥挤小环
//! 是 11 类障碍中反复失手的两类)。**这里不假装它成功了**:退化的地方明确
//! 报出来,下游可以选择拒绝渲染或人工介入。悄悄给一张看着还行、其实构型
//! 读不出来的图,比明说"这一块我画不好"糟得多。
//!
//! # 平局一律按规范秩打破
//!
//! 从哪个环起手、共用键取哪一根 —— 这些选择只要沾上原子的**存储下标**,同一个
//! 分子换一种 SMILES 写法就会得到另一张图。全部改用规范秩,写法就影响不到结果。

use std::collections::{BTreeMap, BTreeSet};

use omgkit_chem::sssr::Ring;
use omgkit_core::MolBuilder;

use crate::geom::{regular_polygon, Point2, BOND_LEN};

/// 布局中不得不退化的地方。
///
/// 标了 `non_exhaustive`:以后加新的退化种类不该是下游的破坏性变更。
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum Degradation {
    /// 桥环或笼状体系:没有规则能给出平面上的好解,坐标由弹簧松弛得到。
    ///
    /// 环内键角、键长都不再保证,重叠也可能消不掉。
    BridgedRingSystem {
        /// 涉及的原子
        atoms: Vec<u32>,
        /// 坐标是不是查表来的,以及没命中时是哪一种没命中。
        ///
        /// # 为什么要分开报
        ///
        /// 命中模板的解是一次昂贵搜索的结果(两万次带扰动的多起点),通常没有
        /// 自交;没命中就只有运行时那 5 个初值的松弛,**实测常常是自交的**。
        /// 两者都叫"退化",但坏的程度差着量级 —— 下游要拒绝渲染还是人工介入,
        /// 分不出来就没法定。
        ///
        /// 没命中还意味着一件可操作的事:**这个骨架该补进
        /// `harness/corpus/bridged.smi` 再重跑生成器。**
        ///
        /// # 没命中之后还分两种,这里**分不出来**
        ///
        /// 弧法接进运行时之后,`NotInTable` / `NoFingerprint` 可能是两种完全
        /// 不同的坐标:
        ///
        /// | | 键长 | 132 个体系里自交的 |
        /// |---|---|---:|
        /// | **弧法**(等张角圆弧) | **精确 1** | **3** |
        /// | 松弛(弧法摆不了时) | 偏差 20%~60% | 75 |
        ///
        /// 两者被压进同一个状态,下游拿它决定"拒绝渲染还是人工介入"时分不出
        /// 好坏。**要分,得再加一档状态** —— 记在这儿,没做。
        template: crate::templates::Status,
    },
    /// η<sup>n</sup> 配位:金属与一整个环之间有 n 根键(二茂铁的 Fe 有 10 根)。
    ///
    /// # 为什么这必然是退化的
    ///
    /// 那 n 根键在三维里等长(金属在环平面的正上方),**在平面上做不到** ——
    /// 把金属摆在环外,离各个环原子的距离必然不同;摆进环里又会压在环上。
    /// 所以这类图**键长一定不全等**,如实报出来。
    ///
    /// 画法本身是好的:布局只留一根代表键,于是 Cp 成了普通五元环、金属成了
    /// 两个五边形之间的连接原子,也就是夹心式。见 `crate::hapto_extras`。
    /// 实测二茂铁由此从「退化 2 / 未解冲突 2 / 交叉 8、两个环叠在一起」变成
    /// 「两个干净的五边形 + 居中的 Fe,交叉 4」。
    HaptoCoordination {
        /// 金属原子
        metal: u32,
        /// 与它 η 配位的那个环上的原子
        ring: Vec<u32>,
    },
}

/// 一个环系统连同落在它里面的 SSSR 环。
pub(crate) struct System<'a> {
    pub atoms: Vec<u32>,
    pub rings: Vec<&'a Ring>,
}

/// 把 SSSR 环按所属的稠环系统归类。
///
/// `fused_ring_systems` 用的是双连通分解,所以**螺环与单键相连的环各成一个
/// 系统**(实测:螺[4.4]壬烷给出两个系统,共用那个螺原子;联苯给出两个系统,
/// 中间那根键不自成系统)。这对布局是好事 —— 各自摆好再接起来即可。
pub(crate) fn group<'a>(systems: &[Vec<u32>], rings: &'a [Ring]) -> Vec<System<'a>> {
    systems
        .iter()
        .map(|atoms| {
            let set: BTreeSet<u32> = atoms.iter().copied().collect();
            System {
                atoms: atoms.clone(),
                rings: rings
                    .iter()
                    .filter(|r| r.atoms.iter().all(|a| set.contains(a)))
                    .collect(),
            }
        })
        .collect()
}

/// 在**局部坐标系**里给一个环系统布局。返回逐原子坐标与(可能的)退化记录。
///
/// 调用方拿到之后再整体平移旋转到该去的位置,见 [`place_at`]。
pub(crate) fn layout_local(
    mol: &MolBuilder,
    sys: &System<'_>,
    ranks: &[u32],
    over: crate::templates::Override<'_>,
) -> (BTreeMap<u32, Point2>, Option<Degradation>) {
    let mut pos: BTreeMap<u32, Point2> = BTreeMap::new();

    if sys.rings.is_empty() {
        // 环感知说这里有环、SSSR 却一个都没给出来。不猜,直接走退化。
        let (pos, st) = relax(mol, &sys.atoms, ranks, &sys.rings, over);
        return (pos, Some(bridged(&sys.atoms, st)));
    }

    // 起手环:先按大小(大环更能定住整体形状),再按规范秩 —— 不看存储下标
    let mut order: Vec<&Ring> = sys.rings.clone();
    order.sort_by_key(|r| (std::cmp::Reverse(r.atoms.len()), ring_key(r, ranks)));
    // **起手环怎么摆,必须完全由规范秩决定。**
    //
    // SSSR 给出的环原子顺序依赖存储序:同一个环,不同写法可能从不同的原子起、
    // 甚至朝相反方向绕。直接拿它去对多边形顶点,两种写法就落在**不同的构型**上
    // (不只是旋转或镜像)—— 后面的取代基方向、消冲突翻哪根键全跟着分岔。
    //
    // 实测:阿司匹林的两种写法,一种消冲突一次没翻,另一种翻了两次,最后坐标
    // 对不上。而"两两距离的多重集"那种指纹**看不出来**,因为点集确实全等。
    let first = canonical_cycle(&order[0].atoms, ranks);
    for (a, p) in first.iter().zip(regular_polygon(first.len(), 0.0)) {
        pos.insert(*a, p);
    }

    let mut placed: BTreeSet<usize> = BTreeSet::from([0]);
    let mut degraded = None;

    // 反复找"只与已放置部分共用一根键"的环,拼上去
    while placed.len() < order.len() {
        let mut best: Option<(usize, u32, u32)> = None;
        for (i, r) in order.iter().enumerate() {
            if placed.contains(&i) {
                continue;
            }
            let shared: Vec<u32> = r
                .atoms
                .iter()
                .copied()
                .filter(|a| pos.contains_key(a))
                .collect();
            if shared.len() != 2 {
                continue; // 共用 1 个是螺(不会同系统)、>2 个是桥
            }
            // **共用的那两个原子要按规范秩定序。** `shared` 保留的是
            // `r.atoms` 的顺序 —— 那是 SSSR 的输出序,依赖存储序。
            //
            // 交换 `u`、`v` 在**非平局**时不改结果(两个候选环心跟着互换,
            // "取远离质心的那个"仍然选中同一个点,转向也自适应)。可一旦
            // **平局** —— 已放置的质心落在这根键的中垂线上 —— `fuse_on_bond`
            // 走 `else` 取 `c2`,而 `c1`/`c2` 正是被 `u`/`v` 交换换掉的那一对,
            // 于是两种写法把新环拼到了相反的一侧。
            //
            // 实测语料第 7879 行(一个 Ni 的四齿配合物,五个环共用同一个金属,
            // 对称度高、质心常落在中垂线上):四次拼环的 `(u,v)` **每一次都
            // 正好反过来**,最后 40/40 个图元全不同。
            let (u, v) = if (ranks[shared[0] as usize], shared[0])
                <= (ranks[shared[1] as usize], shared[1])
            {
                (shared[0], shared[1])
            } else {
                (shared[1], shared[0])
            };
            if !adjacent_in_ring(r, u, v) {
                continue; // 共用两个原子却不相邻 —— 那也是桥
            }
            let key = (ring_key(r, ranks), ranks[u as usize].min(ranks[v as usize]));
            // 不用 `Option::is_none_or` —— 它到 1.82 才稳定,而工作区 MSRV 是 1.75
            let better = match best {
                None => true,
                Some((bi, bu, bv)) => {
                    key < (
                        ring_key(order[bi], ranks),
                        ranks[bu as usize].min(ranks[bv as usize]),
                    )
                }
            };
            if better {
                best = Some((i, u, v));
            }
        }

        let Some((i, u, v)) = best else {
            // 剩下的环都不是邻稠 —— 桥环。整个系统重排,**丢掉部分结果**
            // (见 `relax` 的注释)。三条路按可信度排:
            //
            // 1. **查表**:一次昂贵的离线搜索,而且是按**整分子**打分挑出来的
            //    —— 语料里的桥环骨架 173/177 命中这一档。
            // 2. **弧法**:表没命中时的正解。它有几何保证(键长精确 1),而
            //    松弛是随机重启的抽奖。
            // 3. **松弛**:兜底。弧法自己摆不了时报 `None`,退到这里。
            //
            // 第 2 档的收益全在**语料之外** —— 语料内 173/177 走第 1 档,把
            // 弧法接上前后审计输出逐字节相同。而把表遮住(模拟一个语料里没有
            // 的新骨架)量,弧法能摆的那 132 个体系里:
            //
            // | | 弧法 | 松弛 |
            // |---|---:|---:|
            // | 内部有自交的 | **3** | **75** |
            // | 逐体系比谁交叉少 | **73 胜** | 1 胜 · 58 平 |
            //
            // **鸽笼不必先解。** 立项时以为要先解决"三条等长桥挤两个镜像位"
            // 才能接线;实测那是个不必要的前提 —— 摆不了就退回松弛,而摆得了
            // 的那 74.6% 已经把松弛打得很惨。
            let (pos, st) = relax(mol, &sys.atoms, ranks, &sys.rings, over);
            degraded = Some(bridged(&sys.atoms, st));
            return (pos, degraded);
        };

        fuse_on_bond(order[i], u, v, ranks, &mut pos);
        placed.insert(i);
    }

    (pos, degraded)
}

/// 把一个环的原子序列旋转/翻转到**规范起点与规范方向**。
///
/// 起点取规范秩最小的原子;方向取"沿着走一圈得到的秩序列字典序更小"的那一边。
/// 两个自由度都被定死,同一个环无论怎么写都得到同一个序列。
fn canonical_cycle(atoms: &[u32], ranks: &[u32]) -> Vec<u32> {
    let n = atoms.len();
    let start = (0..n)
        .min_by_key(|i| (ranks[atoms[*i] as usize], atoms[*i]))
        .expect("环非空");
    let fwd: Vec<u32> = (0..n).map(|k| atoms[(start + k) % n]).collect();
    let bwd: Vec<u32> = (0..n).map(|k| atoms[(start + n - k) % n]).collect();
    let key = |v: &[u32]| -> Vec<u32> { v.iter().map(|a| ranks[*a as usize]).collect() };
    if key(&bwd) < key(&fwd) {
        bwd
    } else {
        fwd
    }
}

fn bridged(atoms: &[u32], template: crate::templates::Status) -> Degradation {
    let mut a = atoms.to_vec();
    a.sort_unstable();
    Degradation::BridgedRingSystem { atoms: a, template }
}

/// 环的确定性排序键:环上规范秩的**有序**多重集。
///
/// 用规范秩而不是原子下标,同一分子的不同写法才会选出同一个起手环。
fn ring_key(r: &Ring, ranks: &[u32]) -> Vec<u32> {
    let mut k: Vec<u32> = r.atoms.iter().map(|a| ranks[*a as usize]).collect();
    k.sort_unstable();
    k
}

fn adjacent_in_ring(r: &Ring, u: u32, v: u32) -> bool {
    let n = r.atoms.len();
    (0..n).any(|i| {
        let (a, b) = (r.atoms[i], r.atoms[(i + 1) % n]);
        (a == u && b == v) || (a == v && b == u)
    })
}

/// 沿已放置的键 `u–v` 把环 `r` 拼到外侧。
fn fuse_on_bond(r: &Ring, u: u32, v: u32, ranks: &[u32], pos: &mut BTreeMap<u32, Point2>) {
    let n = r.atoms.len();
    // 把环的原子序列转到以 u 开头、v 紧随其后
    let start = r.atoms.iter().position(|a| *a == u).expect("u 在环上");
    let forward = r.atoms[(start + 1) % n] == v;
    let seq: Vec<u32> = (0..n)
        .map(|k| {
            let i = if forward { start + k } else { start + n - k };
            r.atoms[i % n]
        })
        .collect();
    debug_assert_eq!(seq[0], u);
    debug_assert_eq!(seq[1], v);

    let (pu, pv) = (pos[&u], pos[&v]);
    let mid = (pu + pv) * 0.5;
    let along = (pv - pu).normalized();
    let normal = Point2::new(-along.y, along.x);
    // 边心距:边长 s 的正 n 边形,中心到边的距离是 s / (2 tan(π/n))
    let apothem = BOND_LEN / (2.0 * (std::f64::consts::PI / n as f64).tan());

    // 两个候选中心,取**远离已放置质心**的那个 —— 新环要长在外侧。
    //
    // # 平局要显式判,不能交给浮点比较
    //
    // `c1`/`c2` 沿**法线**对称,所以到两者等距的点集是**过 `mid` 沿键方向的那条
    // 直线**(键所在的直线及其延长线),不是键的中垂线。质心落在它上面时,
    // `c1.dist(anchor)` 与 `c2.dist(anchor)` 数学上相等,**胜负全由舍入决定**。
    //
    // 而 `anchor` 是一串浮点累加,累加次序一变末位就变 —— 那正是写法依赖。
    // 所以两件事都做:
    //
    // 1. 质心**按规范秩累加**(与 `relax` 里那条同一个道理);
    // 2. 用**有符号投影**显式判平局,平局时由规范量选边,不进浮点比较。
    //
    // **第 2 条有判据**(`a_tie_between_the_two_ring_centres_is_decided_by_a_rule_not_by_rounding`,
    // 变异回浮点比较当场红)。**第 1 条没有** —— 全量语料上把它换回按存储序
    // 累加,147 条判据全绿、审计一处不动。它是构造性的防御,不是量到的收益:
    // 语料里没有"累加次序真的翻了符号"的例子。如实记着,不编一个来凑。
    //
    // 实测语料第 7879 行(Ni 四齿配合物,五个环共用同一个金属):四次拼环的
    // 投影绝对值依次是 1.3764 / **0** / 0.9346 / 0.3753 —— **只有第二次是
    // 平局**,另外两次拼歪是因为拼在一个已经分岔了的局部布局上。拿 400 种
    // 真置换改写量:那一次里 **109 种(27%)** 的 gap 已经不是精确 0 了,
    // 只是符号碰巧没翻 —— 语料过得去靠的是运气,不是构造。
    //
    // 非平局那一侧不在刀锋上:全量 5398 次拼环里最小的非零间隔是
    // **0.167**,比浮点噪声高十五个量级,中间一次都没落进 1e-6 以内。
    let mut order: Vec<u32> = pos.keys().copied().collect();
    order.sort_by_key(|a| (ranks[*a as usize], *a));
    let anchor = centroid(order.iter().map(|a| pos[a]));
    let c1 = mid + normal * apothem;
    let c2 = mid - normal * apothem;
    // `normal` 由 `(u, v)` 定,而调用处已按规范秩把它们定序 —— 所以 `c1`
    // 与写法无关,平局时无条件取它就是个规范的选择。
    const TIE: f64 = 1e-9;
    let s = (anchor - mid).dot(normal);
    let center = if s < -TIE {
        c1
    } else if s > TIE {
        c2
    } else {
        c1
    };

    // 转向的正负:取能把 pu 转到 pv 的那一个。这一步顺带验证了几何 ——
    // 边心距或法线写错的话,两个方向都转不到 pv,debug 下会直接断言失败。
    let step = std::f64::consts::TAU / n as f64;
    let sign = if pu.rotated_about(center, step).dist(pv) < 1e-6 {
        1.0
    } else {
        debug_assert!(
            pu.rotated_about(center, -step).dist(pv) < 1e-6,
            "拼环的几何不自洽:两个转向都到不了对面那个原子"
        );
        -1.0
    };

    for (k, a) in seq.iter().enumerate().skip(2) {
        pos.insert(*a, pu.rotated_about(center, sign * step * k as f64));
    }
}

fn centroid(pts: impl Iterator<Item = Point2>) -> Point2 {
    let mut sum = Point2::ORIGIN;
    let mut n = 0.0;
    for p in pts {
        sum = sum + p;
        n += 1.0;
    }
    if n == 0.0 {
        Point2::ORIGIN
    } else {
        sum * (1.0 / n)
    }
}

/// 桥环的兜底:弹簧松弛。
///
/// 键长拉向 [`BOND_LEN`],非键原子互斥。给不出标准键角,也不保证消得掉重叠 ——
/// 这正是调用方要把它记进 [`Degradation`] 的原因。
pub(crate) fn relax(
    mol: &MolBuilder,
    atoms: &[u32],
    ranks: &[u32],
    rings: &[&Ring],
    over: crate::templates::Override<'_>,
) -> (BTreeMap<u32, Point2>, crate::templates::Status) {
    // **先查表。** 松弛是局部下降,5 个初值本身就常常给出自交的解 —— 实测最常见
    // 的 8 个骨架里 5 个自交,双环[2.2.2]辛烷和金刚烷都在内。表里存的是同一个
    // `quality` 口径下搜得久得多的结果,见 [`crate::templates`]。
    // **查表的状态一并带出去。** 先前调用方为了知道"命中没有"又调了一次
    // `lookup`,而那里面是建分子 + sanitize + 规范化 —— 每个桥环系统、每种
    // 规范、审计里每种写法都白付一遍。
    let (hit, status) = crate::templates::lookup_with(mol, atoms, ranks, over);
    if let Some(p) = hit {
        return (p, status);
    }
    // **表没命中就先试弧法。** 它有几何保证(键长精确 1),而下面的松弛是随机
    // 重启的抽奖。摆不了它自己报 `None`,照样落到松弛。
    //
    // 排在**这里**而不是调用处,是为了短路掉松弛:实测遮表跑全量,弧法总耗时
    // 1.56 ms、松弛 71.9 ms(**便宜 46 倍**),而其中 **44.4 ms(61.7%)** 是
    // "弧法赢了、松弛白跑"的。
    //
    // 代价是**弧法赢了就不再与松弛比**。逐体系比是「73 胜 1 负 58 平」,那 1 负
    // 是白输的。不比是有意的:现成的 `quality` 只数**光骨架**的自交,而这一路
    // 反复撞的就是"光骨架的分数不等于整分子的好坏"(见「模板生成器改成按整分子
    // 打分」那一节)。真要挑得对,得按整分子打分 —— 那正是模板表离线在做的事,
    // 也正是查表排第一档的原因。
    if let Some(p) = crate::arcs::place(rings, ranks) {
        return (p, status);
    }
    // **原子按规范秩排序,不按存储下标。** 初值、乃至浮点求和的次序都因此固定,
    // 于是同一个分子的任何写法得到同一张图。
    //
    // 这里刻意**不**接受"贪心走到一半"的部分结果当种子。那个部分结果依赖遍历
    // 顺序,拿它当初值会把写法依赖直接带进来 —— 实测:苊的两种写法就是这样
    // 给出了两个不同形状,而萘、菲、蒽因为太对称,根本触发不到,看着一切正常。
    let mut sorted: Vec<u32> = atoms.to_vec();
    sorted.sort_by_key(|a| (ranks[*a as usize], *a));

    // **多起点。** 松弛是局部下降,落到哪个局部极小全看初值。单一初值下
    // 实测 177 个桥环系统里 172 个(97%)自身有键交叉 —— 那不是消冲突没做好,
    // 消冲突根本够不着:环系统是 2-连通的刚性块,翻转动不了它内部的相对位置。
    //
    // 换几个初值再挑最好的,算法一个字不用改。每个初值都由规范秩派生,
    // 挑选时的平局用量化坐标序列打破,写法无关这条不受影响。
    let mut best: Option<(Quality, BTreeMap<u32, Point2>)> = None;
    for seed in 0..SEEDS {
        let out = relax_from(mol, &sorted, seed, rings, ranks);
        let key = quality(mol, &out, ranks);
        let take = match &best {
            None => true,
            Some((b, _)) => key < *b,
        };
        if take {
            best = Some((key, out));
        }
    }
    (best.expect("SEEDS 至少为 1").1, status)
}

/// 试几个初值。**每多一个都要有个说法**,而且要拿全量语料量过。
///
/// 试过第 6 个(多边形起手 + 新原子放进最大空隙):交叉多消 6 处,写法无关
/// 却多 6 处违例。**写法无关是本库的头号契约,不拿它换**,所以没要。
const SEEDS: usize = 5;

/// 一个松弛解好不好:(系统内自交的键对数, 最大键长偏差, 量化坐标序列)。
///
/// **越小越好。**
type Quality = (usize, i64, Vec<(i64, i64)>);

/// 给一个松弛解打分,口径见 [`Quality`]。
///
/// 键长偏差排第二是因为这条路径上"键长全等"本来就不成立(实测全部 177 个
/// 桥环系统 relax 之后偏差都 ≥20%),但同样糟的两个解里该挑偏差小的。
/// 第三项是平局兜底:不留任意性,同一个分子的任何写法挑到同一个解。
///
/// # 试过把"共线的二度原子数"插进第二档,亏了
///
/// 共线的顶点在图上看不见(渲染那边补一个元素符号才读得出来),所以直觉上
/// 该躲。但**模板的自交是在光骨架上数的,而真正要紧的交叉是带取代基的整个
/// 分子上的** —— 骨架上同样 0 自交的几个解,挂上取代基之后交叉数并不一样,
/// 第二档挑谁对整分子的交叉是不可预测的。全量语料实测(TOP=24):
///
/// | 第二档 | 有键交叉 | 骨架 180° |
/// |---|---:|---:|
/// | 键长偏差(保留) | **116** | 224 |
/// | 共线数(否掉) | 150 | **72** |
///
/// +34 处交叉换 −152 处共线。按本库一贯的次序 —— **交叉可能让人读错,而共线
/// 补了符号之后并不误导** —— 这笔买卖是亏的,所以没要。
fn quality(mol: &MolBuilder, pos: &BTreeMap<u32, Point2>, ranks: &[u32]) -> Quality {
    let live: Vec<(u32, Point2, Point2)> = mol
        .bonds()
        .iter()
        .enumerate()
        .filter_map(|(i, b)| {
            Some((
                u32::try_from(i).ok()?,
                *pos.get(&b.begin)?,
                *pos.get(&b.end)?,
            ))
        })
        .collect();
    let mut cross = 0usize;
    for (k, (_, u1, v1)) in live.iter().enumerate() {
        for (_, u2, v2) in &live[k + 1..] {
            if crate::geom::segments_cross(*u1, *v1, *u2, *v2) {
                cross += 1;
            }
        }
    }
    #[allow(clippy::cast_possible_truncation)]
    let dev = live
        .iter()
        .map(|(_, u, v)| ((u.dist(*v) - BOND_LEN).abs() * 1e6).round() as i64)
        .max()
        .unwrap_or(0);
    // **按规范秩排,不按原子下标。** `BTreeMap` 的迭代序是下标序,拿它当平局
    // 兜底就把写法依赖带了进来 —— 两种写法会在同分的几个解里挑到不同的那个。
    let mut by_rank: Vec<(u32, Point2)> =
        pos.iter().map(|(a, p)| (ranks[*a as usize], *p)).collect();
    by_rank.sort_by_key(|x| x.0);
    #[allow(clippy::cast_possible_truncation)]
    let seq: Vec<(i64, i64)> = by_rank
        .iter()
        .map(|(_, p)| ((p.x * 1e6).round() as i64, (p.y * 1e6).round() as i64))
        .collect();
    (cross, dev, seq)
}

/// 从第 `seed` 个初值出发做一遍松弛。`atoms` 已按规范秩排好。
fn relax_from(
    mol: &MolBuilder,
    atoms: &[u32],
    seed: usize,
    rings: &[&Ring],
    ranks: &[u32],
) -> BTreeMap<u32, Point2> {
    let idx: BTreeMap<u32, usize> = atoms.iter().enumerate().map(|(i, a)| (*a, i)).collect();
    let n = atoms.len();

    // **按规范秩下标定序,不按键的存储序。** `settle` 是照这个次序累加力的,
    // 而浮点加法不满足结合律 —— 存储序随写法变,同一分子的两种写法算出的坐标
    // 就会差最后几位。平时看不出来,**坐标恰好落在四舍五入边界上时就会翻**;
    // 模板换成几何求解之后正是这样炸出来的(镍配合物,差 10.6 个单位)。
    //
    // `idx` 是按规范秩排好的原子在 `atoms` 里的下标,所以按 `(小, 大)` 排序
    // 就与写法无关了。
    //
    // **这一处也没有判据守着。** 与 `place_at` 同理:`settle` 跑 400 步,末位
    // 差别在迭代里既可能放大也可能被吃掉,造不出稳定会红的样本。留着是因为
    // "顺序必须与写法无关"这条本身成立,不是因为量到了收益。
    let mut bonded: Vec<(usize, usize)> = mol
        .bonds()
        .iter()
        .filter_map(|b| Some((*idx.get(&b.begin)?, *idx.get(&b.end)?)))
        .map(|(u, v)| if u <= v { (u, v) } else { (v, u) })
        .collect();
    bonded.sort_unstable();

    // 初值 4:**最大的那个环先摆成正多边形**,其余原子沿着已放好的邻居向外
    // 铺开。前四个初值都是"所有原子摆在一个圆上",拓扑上太像,弹簧下降往往
    // 收敛到同一批坏极小;这个起手的形状不一样,实测它才是降幅的主要来源。
    if seed >= 4 {
        if let Some(p) = polygon_seed(mol, atoms, rings, ranks, &idx) {
            return settle(p, n, &bonded, atoms);
        }
    }

    // 其余初值全部由规范秩派生,不看存储下标:
    //   0 圆环、规范秩序      1 圆环、逆序
    //   2 圆环、BFS 序        3 圆环、隔一个取一个(把成键的原子在圆上分开)
    let order: Vec<usize> = match seed {
        1 => (0..n).rev().collect(),
        2 => bfs_order(n, &bonded),
        3 => (0..n).step_by(2).chain((1..n).step_by(2)).collect(),
        _ => (0..n).collect(),
    };
    let r = BOND_LEN * n as f64 / std::f64::consts::TAU.max(1.0);
    let mut p: Vec<Point2> = vec![Point2::ORIGIN; n];
    for (slot, &i) in order.iter().enumerate() {
        p[i] = Point2::new(r, 0.0).rotated(std::f64::consts::TAU * slot as f64 / n as f64);
    }

    settle(p, n, &bonded, atoms)
}

/// 弹簧松弛本体:键长拉到 1,靠得太近的推开。400 步。
fn settle(
    mut p: Vec<Point2>,
    n: usize,
    bonded: &[(usize, usize)],
    atoms: &[u32],
) -> BTreeMap<u32, Point2> {
    for _ in 0..400 {
        let mut force = vec![Point2::ORIGIN; n];
        for &(i, j) in bonded {
            let d = p[j] - p[i];
            let len = d.norm().max(1e-6);
            let f = d.normalized() * ((len - BOND_LEN) * 0.35);
            force[i] = force[i] + f;
            force[j] = force[j] - f;
        }
        for i in 0..n {
            for j in (i + 1)..n {
                let d = p[j] - p[i];
                let len = d.norm().max(1e-6);
                if len < BOND_LEN * 1.2 {
                    let f = d.normalized() * ((BOND_LEN * 1.2 - len) * 0.25);
                    force[i] = force[i] - f;
                    force[j] = force[j] + f;
                }
            }
        }
        for i in 0..n {
            p[i] = p[i] + force[i];
        }
    }

    atoms.iter().copied().zip(p).collect()
}

/// 初值:系统里最大的那个环摆成正多边形,其余原子沿已放好的邻居向外铺开。
///
/// 放不出来(系统里一个环都没有)时返回 `None`,退回圆环初值。
fn polygon_seed(
    mol: &MolBuilder,
    atoms: &[u32],
    rings: &[&Ring],
    ranks: &[u32],
    idx: &BTreeMap<u32, usize>,
) -> Option<Vec<Point2>> {
    // 起手环:先按大小,再按规范秩 —— 平局不许看存储下标
    let first = rings
        .iter()
        .filter(|r| r.atoms.iter().all(|a| idx.contains_key(a)))
        .min_by_key(|r| (std::cmp::Reverse(r.atoms.len()), ring_key(r, ranks)))?;
    let cyc = canonical_cycle(&first.atoms, ranks);

    let n = atoms.len();
    let mut p = vec![Point2::ORIGIN; n];
    let mut placed = vec![false; n];
    for (a, q) in cyc.iter().zip(regular_polygon(cyc.len(), 0.0)) {
        let i = *idx.get(a)?;
        p[i] = q;
        placed[i] = true;
    }

    // 其余的沿已放好的邻居向外铺:方向取"背离已放好那堆的质心"
    loop {
        let next = atoms.iter().enumerate().find(|(i, a)| {
            !placed[*i]
                && mol
                    .neighbors(**a)
                    .any(|(nb, _)| idx.get(&nb).is_some_and(|j| placed[*j]))
        });
        let Some((i, a)) = next else { break };
        // **锚点按规范秩挑,不看 `neighbors` 的存储序。** 有两个已放好的邻居
        // 可选时,拿存储序挑就把写法依赖直接带了进来 —— 实测全量语料的写法
        // 无关违例会从 129 涨到 349。
        let anchor = mol
            .neighbors(*a)
            .filter_map(|(nb, _)| Some((ranks[nb as usize], nb, *idx.get(&nb)?)))
            .filter(|(_, _, j)| placed[*j])
            .min()
            .map(|(_, _, j)| j)?;
        // 背离已放好那堆的质心
        let dir = {
            let (mut c, mut k) = (Point2::ORIGIN, 0.0_f64);
            for (j, on) in placed.iter().enumerate() {
                if *on {
                    c = c + p[j];
                    k += 1.0;
                }
            }
            let away = (p[anchor] - c * (1.0 / k.max(1.0))).normalized();
            if away.norm() < 1e-9 {
                0.0
            } else {
                away.angle()
            }
        };
        p[i] = p[anchor] + Point2::new(BOND_LEN, 0.0).rotated(dir);
        placed[i] = true;
    }
    // 还有没连上的(理论上不会 —— 环系统是连通的),摊在圆上兜底
    for (i, on) in placed.iter().enumerate() {
        if !on {
            p[i] = Point2::new(BOND_LEN * n as f64, 0.0)
                .rotated(std::f64::consts::TAU * i as f64 / n as f64);
        }
    }
    Some(p)
}

/// 从 0 号原子出发的 BFS 序。邻接表按下标升序,而下标已经是规范秩序,
/// 所以这个序也与写法无关。
fn bfs_order(n: usize, bonded: &[(usize, usize)]) -> Vec<usize> {
    let mut adj = vec![Vec::new(); n];
    for &(i, j) in bonded {
        adj[i].push(j);
        adj[j].push(i);
    }
    for a in &mut adj {
        a.sort_unstable();
    }
    let mut seen = vec![false; n];
    let mut out = Vec::with_capacity(n);
    for start in 0..n {
        if seen[start] {
            continue;
        }
        seen[start] = true;
        let mut q = std::collections::VecDeque::from([start]);
        while let Some(x) = q.pop_front() {
            out.push(x);
            for &y in &adj[x] {
                if !seen[y] {
                    seen[y] = true;
                    q.push_back(y);
                }
            }
        }
    }
    out
}

/// 把一个局部布局整体搬到位:让 `anchor` 落在 `at`,并让整体质心朝 `dir`。
///
/// 两个自由度(平移 + 旋转)刚好被这两个条件定死,不留任意性。
pub(crate) fn place_at(
    local: &BTreeMap<u32, Point2>,
    anchor: u32,
    at: Point2,
    dir: Point2,
) -> BTreeMap<u32, Point2> {
    let a = local[&anchor];
    // **求和的次序必须与写法无关。** `local` 是按原子下标建的 `BTreeMap`,
    // 迭代序就是存储序;浮点加法不满足结合律,同一分子的两种写法算出的质心
    // 会差最后一位(~1e-16),`theta` 跟着差那么一点,**整个环系被转了 1e-16**。
    //
    // 平时看不出来,但下游会把它放大:`orient` 在 24 个候选姿态里挑最小的键,
    // 1e-16 的差别足以让另一个姿态胜出,最终差出 10 个单位。实测就是这么炸的
    // (镍配合物,三条一模一样的配体)。
    //
    // 点集本身与写法无关(几何一样,只是标号不同),所以**按坐标排序再求和**
    // 就定死了 —— 不需要把 `ranks` 传进来。
    //
    // **这一处没有判据守着,如实说。** 试过写一条"同一组点换个键序,`place_at`
    // 输出逐位相同"的判据:造了几组点,质心末位确实差,但那点差被后面
    // `(c - a).normalized()` 的归一化吸收了,输出逐位相同 —— 判据在打不打这个
    // 补丁下都是绿的,**空过的判据不留**。`orient::canonicalise` 那处同类修复
    // 有判据(`shuffling_the_storage_order_does_not_move_the_picture`),因为它
    // 的质心直接进坐标,没有归一化这一步。
    let mut pts: Vec<Point2> = local.values().copied().collect();
    pts.sort_by(|u, v| u.x.total_cmp(&v.x).then(u.y.total_cmp(&v.y)));
    let c = centroid(pts.into_iter());
    let from = (c - a).normalized();
    let to = dir.normalized();
    // from 是零向量只可能出现在"质心恰好落在锚点上"的对称情形,那时转多少都一样
    let theta = if from.norm() < 1e-9 {
        0.0
    } else {
        to.angle() - from.angle()
    };
    local
        .iter()
        .map(|(k, p)| (*k, (*p - a).rotated(theta) + at))
        .collect()
}

#[cfg(test)]
mod tests {
    #[test]
    fn a_tie_between_the_two_ring_centres_is_decided_by_a_rule_not_by_rounding() {
        // **两个候选环心到已放置质心等距时,不能交给浮点比较拍板。**
        //
        // `c1`/`c2` 沿法线对称,所以等距点集是**过 `mid` 沿键方向的那条直线**
        // (不是键的中垂线 —— 中垂线上只有 `mid` 那一个点等距)。质心落在它上面
        // 时 `c1.dist(anchor) > c2.dist(anchor)` 两边数学上相等,谁赢全看舍入,
        // 而质心是一串随写法变次序的浮点累加。
        //
        // 这里把平局摆得干干净净:已放置的只有键的两端,质心**精确等于** `mid`。
        // 断言取 `c1`(= `mid + normal × 边心距`)—— `normal` 由已按规范秩定序的
        // `(u, v)` 决定,所以这是个规范的选择。
        //
        // 变异:把那三支换回 `if c1.dist(anchor) > c2.dist(anchor) { c1 } else
        // { c2 }` → 精确平局走 `else` 取 `c2`,这条当场红。
        use super::*;
        let r = Ring {
            atoms: vec![0, 1, 2, 3],
            bonds: vec![0, 1, 2, 3],
        };
        let mut pos: BTreeMap<u32, Point2> = BTreeMap::new();
        pos.insert(0, Point2::new(0.0, 0.0));
        pos.insert(1, Point2::new(1.0, 0.0));
        let ranks = [0u32, 1, 2, 3];
        // 前提要自己成立:质心必须**精确**落在 mid 上,否则这不是平局
        let anchor = centroid(pos.values().copied());
        assert!(
            anchor.x == 0.5 && anchor.y == 0.0,
            "这个摆法下质心该精确等于 mid,实得 ({}, {})",
            anchor.x,
            anchor.y
        );
        fuse_on_bond(&r, 0, 1, &ranks, &mut pos);
        // normal = (-along.y, along.x) = (0, 1),所以 c1 在上方
        let c2 = pos[&2];
        assert!(
            c2.y > 0.0,
            "平局时该取 `c1`(法线正向那个),实得原子 2 落在 y={:.4}",
            c2.y
        );
    }

    #[test]
    fn the_ring_layout_does_not_care_how_sssr_wrote_the_cycles() {
        // **SSSR 给出的环原子序列只有两个自由度随写法变**:从哪个原子起、
        // 朝哪边绕。这条判据把它们**穷举**掉 —— 每个环转 k 步、按需反向,
        // 断言 `layout_local` 的输出逐点相同。
        //
        // 拿语料第 7879 行(Ni 四齿配合物,五个环共用同一个金属)。它是这条
        // 判据唯一在全量语料上暴露过的分子:拼环时"共用的那两个原子"取自
        // SSSR 输出序,而选环心的那个浮点比较在**平局**时由舍入拍板 ——
        // 两种写法把新环拼到了相反的一侧,40/40 个图元全不同。
        //
        // **不经过 refine/orient/render**,所以断言直接指着根因。
        use super::*;
        let smi = "O=C1C[N+]23CC[N+]45CC(=O)O[Ni]24(O1)(OC(=O)C3)OC(=O)C5";
        let mut m = omgkit_io::smiles::parse(smi).unwrap();
        omgkit_chem::pipeline::sanitize(&mut m).unwrap();
        let ranks = crate::ranks_of(&m);
        let rings_all = omgkit_chem::sssr::ring_set(&m);
        let systems = group(&omgkit_chem::rings::fused_ring_systems(&m), &rings_all);
        let sys = systems
            .iter()
            .max_by_key(|s| s.rings.len())
            .expect("该有一个环系统");
        assert!(sys.rings.len() >= 5, "这个分子该有五个环共用一个金属");

        let base = layout_local(&m, sys, &ranks, None).0;
        assert!(!base.is_empty(), "布局该给出坐标");

        // 把每个环的序列转 k 步、按需反向,重跑
        for k in 0..6usize {
            for rev in [false, true] {
                let rotated: Vec<Ring> = sys
                    .rings
                    .iter()
                    .map(|r| {
                        let n = r.atoms.len();
                        let idx: Vec<usize> = (0..n)
                            .map(|i| if rev { (k + n - i) % n } else { (k + i) % n })
                            .collect();
                        Ring {
                            atoms: idx.iter().map(|i| r.atoms[*i]).collect(),
                            // 键跟着走:`bonds[i]` 连 `atoms[i]` 与 `atoms[i+1]`
                            bonds: (0..n)
                                .map(|i| {
                                    let (x, y) = (idx[i], idx[(i + 1) % n]);
                                    r.bonds[if rev { y } else { x }]
                                })
                                .collect(),
                        }
                    })
                    .collect();
                let shuffled = System {
                    atoms: sys.atoms.clone(),
                    rings: rotated.iter().collect(),
                };
                let got = layout_local(&m, &shuffled, &ranks, None).0;
                assert_eq!(
                    base.len(),
                    got.len(),
                    "环的序列转 {k} 步 rev={rev} 之后原子数变了"
                );
                for (a, p) in &base {
                    let q = got[a];
                    assert!(
                        (p.x - q.x).abs() < 1e-9 && (p.y - q.y).abs() < 1e-9,
                        "环的原子序列转 {k} 步 rev={rev} 之后布局就变了:\
                         原子 {a} 从 ({:.4},{:.4}) 挪到 ({:.4},{:.4})",
                        p.x,
                        p.y,
                        q.x,
                        q.y
                    );
                }
            }
        }
    }

    #[test]
    fn a_bridged_system_is_relaxed_from_several_starts() {
        // 松弛是局部下降,落到哪个局部极小全看初值。单一初值下实测 177 个桥环
        // 系统里 172 个自身有键交叉 —— 而**消冲突根本够不着**:环系统是 2-连通
        // 的刚性块,翻转只动挂在外面的子树,动不了它内部的相对位置。
        //
        // 这条要求多起点确实起作用:把候选砍到一个,下面这些分子的桥环系统
        // 内部就会出现自交。
        let mut won = 0;
        for smi in [
            "CC1(C)[C@@H]2CC[C@@]1(C)C(=O)C2",                          // 樟脑
            "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5", // 吗啡
            "CN1[C@H]2CC[C@@H]1C[C@@H](C2)OC(=O)C(CO)c1ccccc1",         // 阿托品
        ] {
            let mut m = omgkit_io::smiles::parse(smi).unwrap();
            omgkit_chem::pipeline::sanitize(&mut m).unwrap();
            let ranks = omgkit_io::canon::canonical_ranks(&m);
            let systems = omgkit_chem::rings::fused_ring_systems(&m);
            let rs = omgkit_chem::sssr::ring_set(&m);
            let mut checked = 0;
            for sys in group(&systems, &rs) {
                if sys.rings.is_empty() || sys.atoms.len() < 6 {
                    continue;
                }
                let (pos, deg) = layout_local(&m, &sys, &ranks, None);
                // **只算真正走了松弛那条路的系统。** 邻稠系统走的是正多边形
                // 拼接,拿它去和强行松弛比,当然赢 —— 那样这条判据就是空过的。
                if deg.is_none() {
                    continue;
                }
                let single = relax_from(
                    &m,
                    &{
                        let mut a = sys.atoms.clone();
                        a.sort_by_key(|x| (ranks[*x as usize], *x));
                        a
                    },
                    0,
                    &sys.rings,
                    &ranks,
                );
                let (best, _, _) = quality(&m, &pos, &ranks);
                let (one, _, _) = quality(&m, &single, &ranks);
                assert!(
                    best <= one,
                    "{smi}:多起点挑出来的解({best} 处自交)还不如单起点({one} 处)"
                );
                // `best <= one` 是恒真的(初值 0 本来就在候选里),单靠它这条
                // 判据是空过的。真正要守的是**多起点确实赢过单起点**。
                if best < one {
                    won += 1;
                }
                checked += 1;
            }
            assert!(checked >= 1, "{smi}:一个环系统都没查到");
        }
        assert!(
            won >= 1,
            "多起点在这三个桥环分子上一次都没赢过单起点 —— 那它就是白跑的"
        );
    }

    use super::*;
    use omgkit_chem::{rings::fused_ring_systems, sssr::ring_set};

    fn prep(smi: &str) -> MolBuilder {
        let mut m = omgkit_io::smiles::parse(smi).unwrap();
        omgkit_chem::pipeline::sanitize(&mut m).unwrap();
        m
    }

    fn layout(smi: &str) -> (BTreeMap<u32, Point2>, Option<Degradation>) {
        let m = prep(smi);
        let ranks = omgkit_io::canon::canonical_ranks(&m);
        let rings = ring_set(&m);
        let sys = group(&fused_ring_systems(&m), &rings);
        layout_local(&m, &sys[0], &ranks, None)
    }

    fn bond_lengths(m: &MolBuilder, pos: &BTreeMap<u32, Point2>) -> Vec<f64> {
        m.bonds()
            .iter()
            .filter_map(|b| Some(pos.get(&b.begin)?.dist(*pos.get(&b.end)?)))
            .collect()
    }

    #[test]
    fn a_single_ring_is_a_regular_polygon() {
        for smi in ["C1CC1", "C1CCC1", "c1ccccc1", "C1CCCCCC1"] {
            let m = prep(smi);
            let (pos, deg) = layout(smi);
            assert_eq!(deg, None, "{smi} 不该退化");
            assert_eq!(pos.len(), m.num_atoms(), "{smi} 有原子没放上");
            for d in bond_lengths(&m, &pos) {
                assert!((d - BOND_LEN).abs() < 1e-9, "{smi} 键长 {d}");
            }
        }
    }

    #[test]
    fn ortho_fused_rings_share_exactly_one_bond_and_keep_unit_bonds() {
        // 萘、吲哚、芴 —— 逐个环只与已放置部分共用一根键的典型
        for smi in [
            "c1ccc2ccccc2c1",
            "c1ccc2[nH]ccc2c1",
            "c1ccc2c(c1)Cc1ccccc1-2",
        ] {
            let m = prep(smi);
            let (pos, deg) = layout(smi);
            assert_eq!(deg, None, "{smi} 不该退化");
            let ring_atoms: BTreeSet<u32> = pos.keys().copied().collect();
            for b in m.bonds() {
                if ring_atoms.contains(&b.begin) && ring_atoms.contains(&b.end) {
                    let d = pos[&b.begin].dist(pos[&b.end]);
                    assert!((d - BOND_LEN).abs() < 1e-9, "{smi} 环内键长 {d}");
                }
            }
            // 稠环的原子必须两两分开 —— 拼错方向会让新环叠回旧环上,而键长
            // 全都还是 1.0,只看键长发现不了
            let pts: Vec<Point2> = pos.values().copied().collect();
            for i in 0..pts.len() {
                for j in (i + 1)..pts.len() {
                    assert!(pts[i].dist(pts[j]) > 0.5, "{smi} 有两个原子挤在一起");
                }
            }
        }
    }

    #[test]
    fn not_in_the_table_and_no_fingerprint_at_all_are_reported_apart() {
        // 两种"没命中"指向完全不同的动作:一个是"补进 `bridged.smi` 重跑生成器",
        // 另一个是"去查那个分子本身"。混成一档的话审计给的指路是错的 ——
        // 实测全量语料里没命中的那 8 例**全是** `NoFingerprint`。
        use crate::templates::Status;

        // 指纹算得出来、表里没有:一个编出来的大笼
        let m = prep("C1CC2CCC3CCC4CCC5CCC1C1C2C3C4C51");
        let ranks = omgkit_io::canon::canonical_ranks(&m);
        let rs = omgkit_chem::sssr::ring_set(&m);
        let syss = group(&omgkit_chem::rings::fused_ring_systems(&m), &rs);
        let sys = syss.iter().max_by_key(|s| s.atoms.len()).expect("有环系统");
        assert_eq!(
            crate::templates::lookup(&m, &sys.atoms, &ranks).1,
            Status::NotInTable,
            "指纹算得出来、表里没有,该报 NotInTable"
        );

        // 指纹根本算不出来:**二茂铁**。铁与环戊二烯基的五个碳全成键,骨架
        // 全碳化之后那个原子度数 9,`sanitize` 过不去。补语料对它没有用 ——
        // 要改的是骨架抽取本身怎么对待这类 η5 配位。
        //
        // **分子取自 `harness/corpus/large.smi` 第 2135、5558 行**,不是编的:
        // 全量审计报的 8 处 `NoFingerprint` 全部出自这两个分子。
        let mut found = 0usize;
        for smi in [
            "C12C3=C4C5=C1[Fe]23456789C%10C6=C7C8=C9%10",
            "CN(C)C[C-]12C3=C4C5=C1[Fe++]23456789[C-]%10C6=C7C8=C9%10",
        ] {
            let m = prep(smi);
            let ranks = omgkit_io::canon::canonical_ranks(&m);
            let rs = omgkit_chem::sssr::ring_set(&m);
            for sys in group(&omgkit_chem::rings::fused_ring_systems(&m), &rs) {
                if crate::templates::lookup(&m, &sys.atoms, &ranks).1 == Status::NoFingerprint {
                    found += 1;
                }
            }
        }
        assert!(
            found > 0,
            "这两个分子该有环系统报 NoFingerprint,实际一个都没有 —— 判据验不到东西了"
        );
    }

    #[test]
    fn a_bridged_skeleton_says_whether_its_coordinates_came_from_the_table() {
        // 命中模板与只能松弛,坏的程度差着量级:前者是两万次带扰动多起点搜出来
        // 的,通常不自交;后者只有运行时那 5 个初值。两者都叫"退化",下游要
        // 拒绝渲染还是人工介入,分不出来就没法定。
        //
        // 前三个取自 `harness/corpus/bridged.smi`(表里有),最后一个刻意不在表里。
        for (smi, want) in [
            ("C1CC2CCC1CC2", true),
            ("C1C2CC3CC1CC(C2)C3", true),
            (
                "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5",
                true,
            ),
            ("C1CC2CCC3CCC4CCC5CCC1C1C2C3C4C51", false),
        ] {
            let m = prep(smi);
            let ranks = omgkit_io::canon::canonical_ranks(&m);
            let rs = omgkit_chem::sssr::ring_set(&m);
            let syss = group(&omgkit_chem::rings::fused_ring_systems(&m), &rs);
            let sys = syss.iter().max_by_key(|s| s.atoms.len()).expect("有环系统");
            let (_, deg) = layout_local(&m, sys, &ranks, None);
            let Some(Degradation::BridgedRingSystem { atoms, template }) = deg else {
                panic!("{smi} 该报桥环退化,得到 {deg:?}");
            };
            assert_eq!(
                template == crate::templates::Status::Hit,
                want,
                "{smi} 的查表状态报成了 {template:?}"
            );
            // **光报"表里有"是不够的** —— 注释说的是"坐标是不是查表来的"。
            // 实测把 `relax` 里那个查表短路关掉(坐标改由 5 起点松弛给出、
            // 标志位仍报 Hit),先前这条判据是绿的。所以还要验坐标真的等于
            // 表里那一组。
            if want {
                let (pos, _) = layout_local(&m, sys, &ranks, None);
                let (tpl, _) = crate::templates::lookup(&m, &atoms, &ranks);
                let tpl = tpl.expect("报了 Hit 就该查得到");
                for (a, p) in &tpl {
                    let got = pos.get(a).expect("每个原子都该有坐标");
                    assert!(
                        got.dist(*p) < 1e-9,
                        "{smi} 原子 {a}:画出来的 {got:?} 不是表里的 {p:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn a_bridged_system_says_so_instead_of_pretending() {
        // 双环[2.2.2]辛烷:三个六元环两两共用不止一根键,平面上没有好解。
        // 判据不是"画得好看",是"**如实说自己画不好**"。
        let (pos, deg) = layout("C1CC2CCC1CC2");
        assert!(
            matches!(deg, Some(Degradation::BridgedRingSystem { .. })),
            "桥环必须报退化,得到 {deg:?}"
        );
        assert_eq!(pos.len(), 8, "退化也要把每个原子都放上");
        for p in pos.values() {
            assert!(p.x.is_finite() && p.y.is_finite(), "退化解不能给出 NaN");
        }
    }

    #[test]
    fn the_layout_does_not_depend_on_how_the_ring_was_written() {
        // **这条测试挑分子要挑对。** 萘、菲、蒽换写法都给出同一形状 —— 但那
        // 不是因为算法写法无关,而是因为它们太对称,起手环已经被"按大小降序"
        // 定死,平局判据压根没被触发。拿它们当判据是走过场(实测:把 ring_key
        // 改回用存储下标,这三个仍然全绿)。
        //
        // 苊是桥式系统,会落到 relax() 兜底,而兜底以"贪心走到哪一步"为初值
        // —— 写法依赖真正藏在那里。用它才守得住。
        let shapes: Vec<Vec<i64>> = ["C1Cc2cccc3cccc1c23", "c1cc2CCc3cccc(c1)c23"]
            .iter()
            .map(|smi| shape_key(smi))
            .collect();
        assert_eq!(shapes[0], shapes[1], "同一分子的两种写法给出了不同形状");
    }

    /// 形状指纹:两两距离排序后的多重集。与原子编号、平移旋转都无关。
    fn shape_key(smi: &str) -> Vec<i64> {
        let m = prep(smi);
        let ranks = omgkit_io::canon::canonical_ranks(&m);
        let rings = ring_set(&m);
        let sys = group(&fused_ring_systems(&m), &rings);
        let s = sys.iter().max_by_key(|s| s.atoms.len()).expect("有环系统");
        let (pos, _) = layout_local(&m, s, &ranks, None);
        let pts: Vec<Point2> = pos.values().copied().collect();
        let mut ds: Vec<i64> = (0..pts.len())
            .flat_map(|i| ((i + 1)..pts.len()).map(move |j| (i, j)))
            .map(|(i, j)| (pts[i].dist(pts[j]) * 1e4).round() as i64)
            .collect();
        ds.sort_unstable();
        ds
    }
}

/// 桥环骨架坐标表的**生成器**。平时不跑。
///
/// ```shell
/// cargo test -p omgkit-depict --release --lib -- --ignored regenerate_templates --nocapture
/// ```
///
/// 把输出贴进 [`crate::templates`]。与 `harness/gen_elements.py` 生成
/// `element_data.rs` 是同一个路子:**生成脚本进版本库,产物也进版本库**,
/// 谁都能重跑一遍核对。
#[cfg(test)]
mod generator {
    use super::*;
    use crate::geom::Point2;

    const TOP: usize = 50;
    /// 短名单:每个骨架留几个候选交给整分子打分。
    ///
    /// # 这个数是量出来的
    ///
    /// 光骨架的 `Quality` 常常分不出高下 —— 双环[2.2.2]辛烷的 8 个候选**前两档
    /// 完全相同**(自交 0、偏差 0.203),而它们造成的整分子交叉是 2 到 20。
    /// 名单太短就把好解筛掉了:取 8 时有 3 条骨架选中最后一名(说明边界还在
    /// 起约束作用);取 16 之后选中名次一路用到 #8/#9/#11/#13/#14。
    ///
    /// **但 16 不是"被证明够用",是碰巧落得好 —— 别随手调大。** 实测调到 48:
    /// 目标侧继续变好(逐骨架交叉总和 34→30,金刚烷 8→4),而**全量审计的键
    /// 交叉反而从 40 涨到 44**。要动这个数,先确认打分的口径与审计报的量是
    /// 同一个 —— 见 `score_on_molecules` 里那段"数有没有,不是数几处"。
    const SHORTLIST: usize = 16;
    /// 每个骨架的基础搜索预算。
    const TRIES: usize = 20_000;
    /// 基础预算跑完仍自交时,最多再搜到这个数,**一到 0 交叉就停**。
    ///
    /// # 这个数是量出来的
    ///
    /// 吗啡骨架先前定格在自交 1,当时的说法是"这是几何下限"—— **那是错的**。
    /// 它的骨架图是**平面图**(18 原子 22 键,`networkx.check_planarity` 为真),
    /// 按 Fáry 定理必有无交叉的直线画法。加大预算实测:
    ///
    /// | 预算 | 最好的自交数 | 键长偏差 |
    /// |---:|---:|---:|
    /// | 2 万(原) | 1 | 0.230 |
    /// | 22.4 万 | **0** | 0.620 |
    /// | 50 万 | 0 | 0.443 |
    ///
    /// 0 交叉的解键长更不齐,但 [`Quality`] 的次序本来就是交叉优先 ——
    /// 交叉会让人读错,键长不齐只是难看。
    ///
    /// **只有还在自交的骨架才付这笔钱**,而且一到 0 就停,所以总耗时涨得有限。
    const ESCALATE: usize = 400_000;

    fn splitmix(state: &mut u64) -> u64 {
        *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = *state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// 一个候选的整分子分数。**越小越好。**
    ///
    /// # 为什么必须按整分子算
    ///
    /// `Quality` 全在**光骨架**上算,而骨架好看不代表挂上取代基好看 —— 取代基
    /// 从环上哪个方向伸出去取决于环的形状,`Quality` 对此一无所知。实测:
    /// 双环[2.2.2]辛烷的 8 个候选骨架质量**完全等价**(自交 0、偏差 0.203),
    /// 造成的整分子交叉却是 2 到 20,而表里存的正是最差的那个。
    ///
    /// 更根本的是:其中三个候选连**两两距离多重集都相同**(同一个形状,差别只在
    /// 哪个骨架原子落在哪个位置,也就是自同构的选取)。**任何骨架级的几何量对
    /// 自同构都是不变的,结构上就看不见这个差别** —— 只能拿真实分子去问。
    ///
    /// 试过的骨架级代理指标(第二档换成它们,在 187 个分子上量整分子交叉/未解):
    /// 键长偏差 92/329、最近原子距离 92/339、靠太近的原子对数 92/337、
    /// 共线原子数 92/333、回转半径 80/300 —— **整分子打分 60/226**。
    type MolScore = (usize, usize, usize, usize, usize, usize);

    /// 拿真实分子给一个候选打分:`(交叉, 未解冲突, 原子重合, 标签塞不下, 骨架 180°)`。
    ///
    /// 后两档是审核实测补上的:只按前三档挑,**骨架 180° 会从 240 涨到 250、
    /// 标签塞不下从 49 涨到 61**;补上之后前三档一处不掉,这两档反而好过现状
    /// (198 / 44)。
    fn score_on_molecules(mols: &[String], skel: &str, coords: &[(f64, f64)]) -> MolScore {
        let mut out: MolScore = (0, 0, 0, 0, 0, 0);
        for smi in mols {
            let Ok(mut m) = omgkit_io::smiles::parse(smi) else {
                continue;
            };
            if omgkit_chem::pipeline::sanitize(&mut m).is_err() {
                continue;
            }
            omgkit_io::stereo::perceive_bond_stereo(&mut m);
            for style in &crate::style::Style::ALL {
                let d = crate::generate_with(&m, style, Some((skel, coords)));
                let grown = d.drawn(&m);
                let mol = &*grown;
                // **数「这张图有没有」,不是「有几处」。** 审计报的是 17662 张
                // 图里的**发生率**,而打分先前求的是 185 个分子上的**总数** ——
                // 两者不是一回事:把总数从 8 压到 4,完全可能是把 2 张各 4 处
                // 交叉的图变成 4 张各 1 处,总数减半而发生率翻倍。
                //
                // 实测过:按总数打分时把短名单从 16 放到 48,目标侧交叉总和
                // 34→30(变好),而全量审计的键交叉 40→**44**(变坏)。
                // **优化的量必须与报告的量是同一个。**
                out.0 += usize::from(!d.crossings.is_empty());
                // **取代基挤到另一根键上**,与交叉同属"读错结构"那一类:
                // 两根键只差几度,画出来就是一根,读者会整个漏掉一个取代基。
                //
                // 这一档是**看图看出来的** —— 樟脑的偕二甲基塌成一根,而当时
                // 所有指标都说没事(`未解冲突` 按原子间距离判,挤成 6° 的两个
                // 原子相距 0.10 个键长,够不着阈值)。实测把它接进来之前,全量
                // 语料从 15 个分子 62 处涨到了 38 个分子 156 处,没有一条指标
                // 报得出来。
                const CRAMPED: f64 = 15.0;
                let mut cramped = 0usize;
                for at in 0..u32::try_from(mol.num_atoms()).expect("原子数超出 u32") {
                    let here = d.coords[at as usize];
                    let mut angs: Vec<f64> = mol
                        .neighbors(at)
                        .map(|(nb, _)| {
                            (d.coords[nb as usize] - here)
                                .angle()
                                .to_degrees()
                                .rem_euclid(360.0)
                        })
                        .collect();
                    if angs.len() < 3 {
                        continue;
                    }
                    angs.sort_by(|x, y| x.partial_cmp(y).expect("角度不会是 NaN"));
                    for k in 0..angs.len() {
                        if (angs[(k + 1) % angs.len()] - angs[k]).rem_euclid(360.0) < CRAMPED {
                            cramped += 1;
                        }
                    }
                }
                out.5 += usize::from(cramped > 0);
                out.1 += usize::from(!d.unresolved.is_empty());
                // **阈值与 `audit.rs::no_atom_sits_on_another` 同口径(0.05 个
                // 键长)。** 先前写的是 1e-6,严了五万倍 —— 实测 864 个候选里
                // 只有 1 个非零,这一档从没决定过任何一次选择,而审计报的
                // 79 处「原子不重合」违例(距离在 1e-6 与 0.05 之间)它一个
                // 都看不见。**优化的量必须与报告的量是同一个。**
                const OVERLAP: f64 = 0.05;
                for i in 0..d.coords.len() {
                    for j in (i + 1)..d.coords.len() {
                        if d.coords[i].dist(d.coords[j]) < OVERLAP {
                            out.2 += 1;
                        }
                    }
                }
                let labels: Vec<Option<crate::label::Label>> = (0..mol.num_atoms())
                    .map(|a| {
                        crate::render::label_at(
                            mol,
                            u32::try_from(a).expect("原子数超出 u32"),
                            style,
                            &d.coords,
                        )
                    })
                    .collect();
                for b in mol.bonds() {
                    if crate::render::is_squeezed(
                        d.coords[b.begin as usize],
                        d.coords[b.end as usize],
                        labels[b.begin as usize].as_ref(),
                        labels[b.end as usize].as_ref(),
                        style,
                    ) {
                        out.3 += 1;
                    }
                }
                // **要过滤掉 sp 原子。** 审计报的「骨架原子被摆成 180°」用的是
                // `accidental_collinear`,它多一道 sp 过滤(有三键、或两根双键
                // 的本来就该 180°)。裸 `is_collinear` 实测命中 228 次,其中
                // **74 次是真正的 sp** —— 三成是在惩罚正确的丙二烯几何。
                // 这一档实测决定了 54 条骨架里 8 条的选择,不是可有可无的。
                for a in 0..u32::try_from(mol.num_atoms()).expect("原子数超出 u32") {
                    if !crate::render::is_collinear(mol, a, &d.coords) {
                        continue;
                    }
                    let mut doubles = 0usize;
                    let mut triple = false;
                    for (_, bi) in mol.neighbors(a) {
                        match mol.bonds()[bi as usize].order {
                            omgkit_core::BondOrder::Triple => triple = true,
                            omgkit_core::BondOrder::Double => doubles += 1,
                            _ => {}
                        }
                    }
                    if !(triple || doubles >= 2) {
                        out.4 += 1;
                    }
                }
            }
        }
        out
    }

    /// 把一个候选放进短名单:按量化坐标序列去重,按 `Quality` 排序,只留前
    /// [`SHORTLIST`] 个。
    fn offer(
        pool: &mut Vec<(Quality, BTreeMap<u32, Point2>)>,
        q: Quality,
        p: BTreeMap<u32, Point2>,
    ) {
        if pool.iter().any(|(b, _)| b.2 == q.2) {
            return; // 同一个解,不重复入选
        }
        pool.push((q, p));
        pool.sort_by(|x, y| x.0.cmp(&y.0));
        pool.truncate(SHORTLIST);
    }

    /// 把候选按**规范秩**摊平成表里那种坐标数组。
    fn flatten(p: &BTreeMap<u32, Point2>, ranks: &[u32]) -> Vec<(f64, f64)> {
        let mut v: Vec<(u32, Point2)> = p.iter().map(|(a, q)| (ranks[*a as usize], *q)).collect();
        v.sort_by_key(|x| x.0);
        v.iter().map(|(_, q)| (q.x, q.y)).collect()
    }

    /// 从短名单里挑一个:**整分子打分定胜负,骨架 `Quality` 只当平局兜底。**
    ///
    /// 抽成函数是为了判据能调它 —— 判据自己再写一遍选择逻辑的话,改坏了选择
    /// 它照样绿(实测:把整分子分数从排序键里去掉,自己算分数的那版判据不红)。
    fn pick_best(
        pool: Vec<(Quality, BTreeMap<u32, Point2>)>,
        mols: &[String],
        skel: &str,
        ranks: &[u32],
    ) -> Option<(Quality, BTreeMap<u32, Point2>)> {
        pool.into_iter()
            .map(|(q, p)| {
                let flat = flatten(&p, ranks);
                (score_on_molecules(mols, skel, &flat), q, p)
            })
            // # 次序:挤压 → 交叉 → 冲突 → 重合 → 塞不下 → 共线 → `Quality`
            //
            // **挤压排第一,连交叉都排在它后面。** 本库明写过的轻重是"共线可以
            // 补符号(渲染会给它补一个元素符号,读者还看得见),交叉会让人读错"
            // —— 而取代基挤成一根是**读错那一类,且没有任何补救**:那个原子在
            // 图上彻底消失,没有符号可补。
            //
            // 三种次序在全量语料上实测:
            //
            // | 排序 | 交叉 | 挤压 | 共线 180° |
            // |---|---:|---:|---:|
            // | 旧表(按骨架挑) | 62 | ~30 | 184 |
            // | 交叉优先 | **38** | 74 | **72** |
            // | **挤压优先(采用)** | 50 | **28** | 102 |
            //
            // 挤压优先在三项上**全面好过旧表**;相对"交叉优先"是用 +12 处交叉、
            // +30 处共线换掉 46 处挤压 —— 按上面的轻重,这笔买卖是赚的。
            .min_by(|x, y| {
                let key = |v: &(MolScore, Quality, BTreeMap<u32, Point2>)| {
                    (
                        (v.0 .5, v.0 .0, v.0 .1, v.0 .2, v.0 .3, v.0 .4),
                        v.1.clone(),
                    )
                };
                key(x).cmp(&key(y))
            })
            .map(|(_, q, p)| (q, p))
    }

    /// 搜一条骨架,返回它在表里那一行。搜不出来(解析/sanitize 失败等)返回 `None`。
    ///
    /// 抽成函数是为了能并行 —— 各条骨架彼此独立,`std::thread::scope` 一 spawn
    /// 就行,每条自己一条种子流,确定性不受影响。
    fn one_skeleton(skel: &str, n: usize, mols: &[String]) -> Option<String> {
        let mut m = omgkit_io::smiles::parse(skel).ok()?;
        omgkit_chem::pipeline::sanitize(&mut m).ok()?;
        let ranks = omgkit_io::canon::canonical_ranks(&m);
        let atoms: Vec<u32> = (0..u32::try_from(m.num_atoms()).ok()?).collect();
        let mut sorted = atoms.clone();
        sorted.sort_by_key(|a| (ranks[*a as usize], *a));
        let cnt = sorted.len();
        let idx: BTreeMap<u32, usize> = sorted.iter().enumerate().map(|(i, a)| (*a, i)).collect();
        let bonded: Vec<(usize, usize)> = m
            .bonds()
            .iter()
            .filter_map(|b| Some((*idx.get(&b.begin)?, *idx.get(&b.end)?)))
            .collect();

        let rs = omgkit_chem::sssr::ring_set(&m);
        let sys = group(&omgkit_chem::rings::fused_ring_systems(&m), &rs);
        let s = sys.iter().max_by_key(|s| s.atoms.len())?;

        // **不能调 `relax`** —— 它先查表,于是 `best` 的初值来自它正在生成的
        // 那张表,生成器就不是语料的纯函数了。更要命的是升不升级由 `best` 定:
        // 表里已经是 0 交叉的骨架,一进来就判"不用升级",**永远搜不到更好的解**。
        // 实测就是这么栽的:morphine 的偏差卡在 0.620,而跑满能到 0.443。
        //
        // 所以这里把 `relax` 的多起点部分照抄一遍,**只是不查表**。
        // **留一份短名单,不是只留一个。** 光骨架的 `Quality` 常常分不出高下,
        // 真正的差别要拿真实分子才问得出来 —— 见 `score_on_molecules`。
        let mut pool: Vec<(Quality, BTreeMap<u32, Point2>)> = Vec::new();
        for seed in 0..SEEDS {
            let out = relax_from(&m, &sorted, seed, &s.rings, &ranks);
            let q = quality(&m, &out, &ranks);
            offer(&mut pool, q, out);
        }

        let mut st = 0x51ED_270B_D5AB_C0DEu64 ^ (cnt as u64);
        let r = BOND_LEN * cnt as f64 / std::f64::consts::TAU;
        // 基础预算 `TRIES`;跑完还自交才接着搜到 `ESCALATE`。
        //
        // **升不升级在 `k == TRIES` 处一次性决定,决定了就跑满。** 先前写的是
        // "一到 0 交叉就 break",而 [`Quality`] 是三档的 —— 自交归零之后第二档
        // "键长偏差"就不再优化了,拿到的是**第一个** 0 交叉解而不是**最好的**
        // 那个。实测差得很多:morphine 那条偏差 0.620,跑满是 0.443;34 原子
        // 那条 0.384 → 0.293。全量语料的交叉/退化/冲突一处不动,纯赚。
        for k in 0..ESCALATE.max(TRIES) {
            if k == TRIES && pool.first().is_some_and(|(q, _)| q.0 == 0) {
                break;
            }
            let mut p = vec![Point2::ORIGIN; cnt];
            for (i, q) in p.iter_mut().enumerate() {
                let j = (splitmix(&mut st) % 1000) as f64 / 1000.0 - 0.5;
                let t = std::f64::consts::TAU * (i as f64 + j * 3.0) / cnt as f64;
                let rad = r * (1.0 + ((splitmix(&mut st) % 1000) as f64 / 1000.0 - 0.5) * 0.6);
                *q = Point2::new(rad, 0.0).rotated(t);
            }
            let out = settle(p, cnt, &bonded, &sorted);
            let q = quality(&m, &out, &ranks);
            offer(&mut pool, q, out);
        }
        // 弧法的解也进名单(它摆不出来时自己报 `None`)
        if let Some(p) = crate::arcs::place(&s.rings, &ranks) {
            let q = quality(&m, &p, &ranks);
            offer(&mut pool, q, p);
        }

        let best = pick_best(pool, mols, skel, &ranks)?;

        // 按**骨架自己的规范秩**存坐标,查表时才对得上
        let mut by_rank: Vec<(u32, Point2)> = best
            .1
            .iter()
            .map(|(a, p)| (ranks[*a as usize], *p))
            .collect();
        by_rank.sort_by_key(|x| x.0);
        let mut line = format!("    (\"{skel}\", &[");
        for (_, p) in &by_rank {
            line.push_str(&format!("({:.6}, {:.6}), ", p.x, p.y));
        }
        // **偏差也打出来。** 先前只打自交数,而"提前退出"丢的恰恰是偏差 ——
        // 不打出来,那个错就不会出现在 diff 里。
        #[allow(clippy::cast_precision_loss)]
        let dev = best.0 .1 as f64 / 1e6;
        line.push_str(&format!(
            "]),   // 出现 {n} 次,自交 {},偏差 {dev:.3}",
            best.0 .0
        ));
        Some(line)
    }

    /// 整分子打分**确实在定胜负**,而不是摆设。
    ///
    /// # 为什么拿双环[2.2.2]辛烷
    ///
    /// 它的候选里有一批**骨架 `Quality` 前两档完全相同**的(自交 0、偏差 0.203),
    /// 而它们造成的整分子交叉是 2 到 20 —— 骨架级的任何量都分不出高下,只能拿
    /// 真实分子去问。其中三个候选连**两两距离多重集都相同**(同一个形状,差别
    /// 只在哪个骨架原子落在哪个位置,也就是自同构的选取),**任何骨架级几何量
    /// 对自同构都是不变的,结构上就看不见**。
    ///
    /// 这条比"去找骨架自交更少但整分子更差的对立"结实得多 —— 那种严格对立
    /// 全语料只有 1 条骨架撑着,搜索预算一改就可能整条消失。
    #[test]
    fn the_whole_molecule_score_is_what_decides() {
        let skel = "C1CC2CCC1CC2";
        let mut m = omgkit_io::smiles::parse(skel).expect("该能解析");
        omgkit_chem::pipeline::sanitize(&mut m).expect("该能 sanitize");
        let ranks = omgkit_io::canon::canonical_ranks(&m);
        let rs = omgkit_chem::sssr::ring_set(&m);
        let sys = group(&omgkit_chem::rings::fused_ring_systems(&m), &rs);
        let s = sys.iter().max_by_key(|s| s.atoms.len()).expect("有环系");

        // **取自 `harness/corpus/large.smi`,不是编的。** 这一点要紧:头一版
        // 我自己写了四个简单取代的双环[2.2.2]辛烷,结果按 `Quality` 排第一的
        // 候选就已经是 0 交叉 —— 判据自己的守卫当场拦下,说"验不出东西"。
        // 真实分子挂着大取代基,才问得出候选之间的差别。
        let mols: Vec<String> = [
            r"C1CC2CCN1C(=C/c1cnccc1)\C2=O",
            "C1C[S+]2CC[S+]1CC2",
            "CCCCC12CCC(CC1=O)(CC2)O",
            "CCOC(C1C(C2CCC1CC2)C(=O)OCC)=O",
            "COC(C1=C[C@H]2[C@@H](C[C@@H]1OC2=O)C#N)=O",
            "COc1c(c(ccc1/C=C1/C(C2CCN1CC2)=O)OC)OC",
            "COc1cc2c(ccnc2cc1)[C@H](C1CC2CCN1CC2CC)O",
            "COc1ccccc1/C=C1/C(C2CCN1CC2)=O",
        ]
        .iter()
        .map(|s| (*s).to_string())
        .collect();

        // 攒一批候选:5 个初值 + 一批带扰动的多起点
        let mut sorted: Vec<u32> = (0..u32::try_from(m.num_atoms()).unwrap()).collect();
        sorted.sort_by_key(|a| (ranks[*a as usize], *a));
        let cnt = sorted.len();
        let idx: BTreeMap<u32, usize> = sorted.iter().enumerate().map(|(i, a)| (*a, i)).collect();
        let mut bonded: Vec<(usize, usize)> = m
            .bonds()
            .iter()
            .filter_map(|b| Some((*idx.get(&b.begin)?, *idx.get(&b.end)?)))
            .map(|(u, v)| if u <= v { (u, v) } else { (v, u) })
            .collect();
        bonded.sort_unstable();

        let mut pool: Vec<(Quality, BTreeMap<u32, Point2>)> = Vec::new();
        for seed in 0..SEEDS {
            let out = relax_from(&m, &sorted, seed, &s.rings, &ranks);
            let q = quality(&m, &out, &ranks);
            offer(&mut pool, q, out);
        }
        let mut st = 0x51ED_270B_D5AB_C0DEu64 ^ (cnt as u64);
        let r = BOND_LEN * cnt as f64 / std::f64::consts::TAU;
        for _ in 0..4000 {
            let mut p = vec![Point2::ORIGIN; cnt];
            for (i, q) in p.iter_mut().enumerate() {
                let j = (splitmix(&mut st) % 1000) as f64 / 1000.0 - 0.5;
                let t = std::f64::consts::TAU * (i as f64 + j * 3.0) / cnt as f64;
                let rad = r * (1.0 + ((splitmix(&mut st) % 1000) as f64 / 1000.0 - 0.5) * 0.6);
                *q = Point2::new(rad, 0.0).rotated(t);
            }
            let out = settle(p, cnt, &bonded, &sorted);
            let q = quality(&m, &out, &ranks);
            offer(&mut pool, q, out);
        }
        assert!(pool.len() >= 4, "只攒到 {} 个候选,验不出东西", pool.len());

        let flat = |p: &BTreeMap<u32, Point2>| flatten(p, &ranks);
        let scores: Vec<MolScore> = pool
            .iter()
            .map(|(_, p)| score_on_molecules(&mols, skel, &flat(p)))
            .collect();

        // 候选之间的整分子交叉**确实不同** —— 否则这条判据是空过的
        let (lo, hi) = (
            scores.iter().map(|s| s.0).min().expect("非空"),
            scores.iter().map(|s| s.0).max().expect("非空"),
        );
        assert!(hi > lo, "候选的整分子交叉全是 {lo},分不出高下,判据空过");

        // 按 `Quality` 排第一的那个,**不是**整分子最好的那个
        let by_quality = scores[0].0;
        assert!(
            by_quality > lo,
            "`Quality` 排第一的整分子交叉是 {by_quality},已经是最好的 {lo} —— \
             这条判据在这个骨架上验不出东西了,该换骨架"
        );

        // **实现真的挑了整分子最好的那个。** 这一句必须调 `pick_best`,不能
        // 自己再写一遍选择逻辑 —— 自己写的话,把整分子分数从排序键里去掉,
        // 判据照样是绿的(实测过)。
        let picked = pick_best(pool, &mols, skel, &ranks).expect("名单非空");
        let got = score_on_molecules(&mols, skel, &flat(&picked.1)).0;
        assert_eq!(
            got, lo,
            "实现挑出来的整分子交叉是 {got},名单里最好的是 {lo}"
        );
    }

    #[test]
    #[ignore]
    fn regenerate_templates() {
        // 一、扫语料,按出现次数排出最常见的桥环骨架
        //
        // **两份语料。** `large.smi` 是通用语料,按频次取前 `TOP`;`bridged.smi`
        // 是专门收"经典难画"的桥环骨架(生物碱、萜类、教科书上的笼),里面每一
        // 个骨架**无条件全收** —— 它们在通用语料里各出现零到一次,按频次排永远
        // 挤不进前 `TOP`,而正是它们让吗啡那类分子画成一团乱麻。
        let text = std::fs::read_to_string("../../harness/corpus/large.smi").unwrap();
        let extra = std::fs::read_to_string("../../harness/corpus/bridged.smi").unwrap();
        let mut freq: BTreeMap<String, usize> = BTreeMap::new();
        let mut must: BTreeSet<String> = BTreeSet::new();
        // 骨架 → 用它的真实分子(规范 SMILES)。打分就拿这些分子跑。
        let mut by_skel: BTreeMap<String, Vec<String>> = BTreeMap::new();
        for (line, is_extra) in text
            .lines()
            .map(|l| (l, false))
            .chain(extra.lines().map(|l| (l, true)))
        {
            let smi = line.split_whitespace().next().unwrap_or("");
            if smi.is_empty() || smi.starts_with('#') {
                continue;
            }
            let Ok(mut m) = omgkit_io::smiles::parse(smi) else {
                continue;
            };
            if omgkit_chem::pipeline::sanitize(&mut m).is_err() {
                continue;
            }
            if m.num_atoms() < 2 {
                continue;
            }
            let ranks = omgkit_io::canon::canonical_ranks(&m);
            let rs = omgkit_chem::sssr::ring_set(&m);
            // 这个分子里退化的环系有几个 —— **打分只收恰好一个的**,见下
            let mut degraded_here: Vec<String> = Vec::new();
            for sys in group(&omgkit_chem::rings::fused_ring_systems(&m), &rs) {
                let (_, deg) = layout_local(&m, &sys, &ranks, None);
                if deg.is_none() {
                    continue;
                }
                if let Some(k) = crate::templates::skeleton_of(&m, &sys.atoms, &ranks) {
                    degraded_here.push(k.clone());
                    // **频次只由 `large.smi` 定。** 额外语料若也计数,它贡献的
                    // 那一两次会把排名搅动 —— 实测原本第 47~50 名的 4 个骨架
                    // 被挤出了前 `TOP`,凭空丢了模板。额外语料只管覆盖面。
                    if is_extra {
                        must.insert(k);
                    } else {
                        *freq.entry(k).or_default() += 1;
                    }
                }
            }
            // **打分分子只收"恰好含一个退化环系"的。**
            //
            // 含两个的话,打分时另一个会去查**正在生成的那张表** —— 生成器就
            // 不再是语料的纯函数,「把 `TABLE` 清空重跑逐字节相同」这条验收
            // 当场作废。实测目前一个都没有,但那是语料的偶然性质,不是结构
            // 保证:往 `bridged.smi` 里加一个双桥环分子就破。
            if degraded_here.len() == 1 {
                by_skel
                    .entry(degraded_here.remove(0))
                    .or_default()
                    .push(omgkit_io::canon::canonical_smiles(&m).smiles);
            }
        }
        let mut v: Vec<(String, usize)> = freq.into_iter().collect();
        // 频次降序、同频按骨架字典序 —— 与写法无关,重跑逐字节可复现
        v.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
        // 前 `TOP` 名,加上 `bridged.smi` 里的全部
        let mut keep: Vec<(String, usize)> = v
            .iter()
            .enumerate()
            .filter(|(i, (k, _))| *i < TOP || must.contains(k))
            .map(|(_, kv)| kv.clone())
            .collect();
        // `must` 里可能有 `large.smi` 一次都没出现过的骨架 —— 它们不在 `v` 里
        for k in &must {
            if !keep.iter().any(|(s, _)| s == k) {
                keep.push((k.clone(), 0));
            }
        }
        keep.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));

        // **打分分子的清单必须可复现。** 排序 + 去重,与文件行序、与写法都无关。
        for v in by_skel.values_mut() {
            v.sort_unstable();
            v.dedup();
        }

        // 二、对每个骨架跑带扰动的多起点,按现成的 Quality 挑最好的
        //
        // **各条骨架彼此完全独立,所以并行跑。** 每条自己一条 splitmix 种子流、
        // 结果按下标收回再统一打印 —— 确定性一点不掉。去掉提前退出之后串行要
        // 15 分钟,并行之后由最长的那条决定。
        let lines: Vec<Option<String>> = std::thread::scope(|sc| {
            let handles: Vec<_> = keep
                .iter()
                .map(|(skel, n)| {
                    let ms = by_skel.get(skel).cloned().unwrap_or_default();
                    sc.spawn(move || one_skeleton(skel, *n, &ms))
                })
                .collect();
            handles.into_iter().map(|h| h.join().unwrap()).collect()
        });

        println!("// 本表由 rings.rs 的 `regenerate_templates` 生成,勿手改。");
        println!("pub(crate) const TABLE: &[(&str, &[(f64, f64)])] = &[");
        let mut kept = 0usize;
        for l in lines.into_iter().flatten() {
            println!("{l}");
            kept += 1;
        }
        println!("];");
        println!("// 共 {kept} 条");
    }
}
