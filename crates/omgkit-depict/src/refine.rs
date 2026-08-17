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
//! | **撑开**(把支点一侧转开一档) | 键长仍**一点不变**,只有**支点处**的键角变 | 实现 |
//! | 伸缩键长 | 键长不再全等 | 未做 |
//!
//! 前三个是**等距变换**(轴过某个不动的原子,反射保持到它的距离),不损失任何
//! 几何性质,所以先试。
//!
//! **撑开是唯一动几何的算子,排在最后,而且带闸门。** 它付的账只有一处,而且
//! 可以证明只有这一处:键长**全部**保住,变的只是支点那一个角(证明写在
//! `Splay` 的文档里)。所以上面那两条判据里,"键长全等"照旧守得住,松动的只有
//! "键角标准" —— 而本库在布局阶段早就允许键角按 ±30° 一档铺开,撑开只是把同一
//! 把尺子搬过来用,且**只走一档**。
//!
//! 仍然解决不了的情形如实报出来([`Report::unresolved`]),不靠拉扯几何把数字
//! 做好看。
//!
//! 后三个各补一个**能力空洞**,收益都要如实说(全量语料 17662 个分子×规范):
//!
//! | 算子 | 补的空洞 | 转干净 | 键交叉 |
//! |---|---|---:|---:|
//! | 螺环翻转 | 可翻转的键排除环上的键,螺环两侧的相对朝向动不了 | **+6** | 0 |
//! | 度 4 置换 | 垂直于翻转轴的那一对邻居,怎么翻都换不过来 | **+25** | −2 |
//! | 撑开 | 三个等距算子都改不动"两个原子精确重合" | **+2** | **+4** |
//!
//! **三行不是同一次实测**:前两行是加那两个算子时量的,撑开那一行是加撑开时量的
//! (`原子不重合` 8 → **2** 是同一次的另一面)。后来「挑方向」那一步把基线整体
//! 挪过一次,**当前值一律见 `harness/README.md`**,不要拿这三行相加。
//!
//! 做它们是因为便宜且不碰契约,不是因为它们大。
//!
//! # 碰撞半径来自标签,而标签尺寸来自规范
//!
//! 两个原子撞没撞上,取决于它们的**标签**占多大 —— 这就是
//! [`Style`] 必须参与布局的地方。同一张图在 ACS 规范下
//! (标签占 0.69 个键长)会撞,在 ChemDraw 默认规范下(0.33)可能不撞。

use std::collections::{BTreeMap, BTreeSet};

use omgkit_core::{BondFlags, BondOrder, MolBuilder};

use crate::geom::{segments_cross, Point2};
use crate::label::{label_for, HSide, LabelPlace};
use crate::style::Style;

/// 没有标签的骨架碳,占位半径。
///
/// 两个裸碳中心靠得比一个键长的一半还近,图上就读不清了 —— 取 0.25 使阈值恰好
/// 落在半个键长。纯靠标签留白算出来的半径(ACS 下约 0.11)太小,骨架上的碰撞
/// 会漏判。
const BARE_RADIUS: f64 = 0.25;

/// 一次消冲突的结果。
///
/// `#[non_exhaustive]`:这里记的是"用过哪几个算子",而算子还会加(撑开就是后加
/// 的)。外部按字段全列构造它没有意义 —— 只有本模块的消冲突造得出有内容的实例。
#[derive(Debug, Clone, PartialEq, Eq, Default)]
#[non_exhaustive]
pub struct Report {
    /// 翻转过的键
    pub flipped: Vec<u32>,
    /// **翻转过的螺原子**。绕螺原子把一个环整体镜像过去,见 `SpiroFlip`。
    pub spiro_flipped: Vec<u32>,
    /// **做过分支对调的度 4 原子**,见 `PairSwap`。
    pub swapped: Vec<u32>,
    /// **撑开过的支点原子**,见 `Splay`。与上面三个同级 —— 都是"这一步做了什么"
    /// 的记录,都不进 [`crate::Depiction`]。
    pub splayed: Vec<u32>,
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
    if best == (0, 0, 0.0) {
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
    //
    // **撑开的候选也要算进来。** 它自己不会失控(每次接受 `same` 至少 −1,是有界
    // 的非负整数),但它会给前三个算子多制造几轮机会,轮数不够就会在半路停下。
    //
    // 这里按 `2 * cands` 估(每根键两个转向),而撑开的键集**比 `cands` 大** ——
    // 它还收端点键,见 `splays`。所以这只是个够用的下界估计,不是精确计数;
    // 上限本来就是防抖动的兜底,估宽估窄都不影响正确性,只影响什么时候放弃。
    //
    // 用尽上限时结果是"当前状态",不是错的:每一步都只在严格改善时才接受,
    // 而候选序全按规范秩定,所以停在哪儿都仍然与写法无关。
    let max_rounds = (cands.len() * 3 + spiros.len() + swaps.len()).max(1) * 2;
    for _ in 0..max_rounds {
        let mut improved = false;
        for &b in &cands {
            let Some(side) = far_side(mol, pos, b, ranks) else {
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

        // **撑开排最后,而且只在前三个这一轮都没改善时才动。**
        //
        // 前三个不花几何代价,谁先谁后无所谓;这一个花代价,所以先让免费的算子
        // 把能消的消掉 —— 否则一处本可由一次翻转解决的重合,会先被撑开付掉一个
        // 键角。
        //
        // `best.0 > 0` 是闸门:**只许用来消掉原子精确重合**,那是唯一值得拿键角
        // 去换的一档。口径提醒:`score.0` 数的是 `remaining()` 给的对,而那一步
        // **跳过成键的对**,硬判据 `原子不重合` 不跳 —— 两边阈值都是 0.05,但
        // 只在**非键对**上口径一致。两个成键原子叠在一起时这个闸门是瞎的。
        if !improved && best.0 > 0 {
            let mut ranked: Vec<(SplayKey, Splay)> = Vec::new();
            for sp in splays(mol, pos, ranks) {
                let before = narrowest_at(mol, pos, sp.pivot);
                let saved: Vec<(u32, Point2)> = sp.moved.iter().map(|a| (*a, pos[a])).collect();
                let c = pos[&sp.pivot];
                for a in &sp.moved {
                    let p = pos[a].rotated_about(c, sp.by);
                    pos.insert(*a, p);
                }
                let now = score(mol, pos, &radii);
                let after = narrowest_at(mol, pos, sp.pivot);
                let ok = now.0 < best.0 && keeps_stereo(pos) && angle_survives(before, after);
                // 闸门蕴含 `better`:字典序第一位就赢。钉住这条兼容性 —— 它是
                // "四个算子混在一起一定收敛"的依据(`best` 严格单调下降)。
                debug_assert!(!ok || better(now, best), "撑开被接受了,却不算改善");
                for (a, p) in saved {
                    pos.insert(a, p);
                }
                if ok {
                    ranked.push((
                        (
                            now.0,
                            now.1,
                            -q6(after),
                            ranks[sp.pivot as usize],
                            sp.moved
                                .iter()
                                .map(|a| ranks[*a as usize])
                                .min()
                                .unwrap_or(u32::MAX),
                            sp.toward,
                        ),
                        sp,
                    ));
                }
            }
            // **选最好的一个,不是选第一个能用的。** 前三个算子不花代价,先到先得
            // 无所谓;这一个花代价,必须挑代价最小的 —— 同样消掉一处重合,可能
            // 白白把一个角压到 90°,而另一个候选一个角都不压窄。
            //
            // 定序键里**故意没有 `depth`**:那是一串按原子下标求和的浮点,同一个
            // 分子换种写法末位会不同。实测 1068 有 4 个候选并列在
            // `(same=0, 交叉=3)`、最窄角全是 90.0° —— 真让 `depth` 进来,胜负就
            // 完全压在一个量化边界上,而头号契约现在是 0 违例,没有余量可赌。
            // 平局改由规范秩接住,那是单射的,一定分得出。
            if let Some((_, sp)) = ranked.into_iter().min_by(|x, y| x.0.cmp(&y.0)) {
                let c = pos[&sp.pivot];
                for a in &sp.moved {
                    let p = pos[a].rotated_about(c, sp.by);
                    pos.insert(*a, p);
                }
                best = score(mol, pos, &radii);
                report.splayed.push(sp.pivot);
                improved = true;
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
            // 氢挂哪一侧此刻还定不下来(它要看最终坐标)。**左右两侧其实一样宽**
            // —— 同一个多重集,只在求和次序造成的最后一位上差一点(实测
            // `[13CH4]`、`C[SiH3]` 会差)。`max` 留着是为了确定性,不是为了取宽的。
            // 半径宁可**偏大**:偏大只是把原子推得开一点,偏小会漏判碰撞
            // **这里仍然用"盒心在原子上"的近似,而盒心其实偏开了 `dx`**
            // (整串挪开好让元素符号落在原子位置上)。所以长的那一侧被低估约
            // `dx`(ACS 下 `OH` 约 0.25 个键长),短的那一侧被高估同样多。
            //
            // 试过改成"原子到盒最远那个角"—— 各向同性地取最大,长边诚实了、
            // 短边跟着过估,**全量语料干净率 91.6% → 89.3%,未解冲突 +401**。
            // 亏的,所以没要。
            //
            // **"改成偏心的盒对盒"这条路走不通,不是没空做,是流水线不允许。**
            //
            // 标签盒在**画布**里是轴对齐的。可 `refine` 跑完之后还要过
            // [`crate::orient::canonicalise`],它在 **30° 的倍数**里挑姿态 ——
            // 也就是说 `refine` 这个坐标系里的"水平"并不是最终的水平。实测
            // 语料前 300 个分子 × 2 规范共 600 次:最终转角落在 90° 倍数上的
            // 只有 **38%**,**62% 是斜的**。在这里摆一个轴对齐的矩形,六成
            // 情形下方向就是错的。
            //
            // 圆是"对旋转角一无所知"时唯一诚实的形状 —— 它正是盒在所有朝向上
            // 的**上界**。所以这里的各向同性不是偷懒,是这个位置能给出的正解。
            //
            // 那个前置条件 —— **把 `canonicalise` 挪到 `refine` 之前** —— 试过
            // 了,**代价是头号契约**:全量语料 `写法无关` 从 9 处涨到 **21 处**
            // (多出来的 12 处全是"形状相同、摆位不同"),而干净率、未解冲突、
            // 交叉一档没动。摆正一旦不在最后一步,消冲突之后就没有东西再把姿态
            // 归位了 —— 同一个分子的两种写法会停在两个姿态上。
            //
            // 所以这条路是:**先付 12 处头号契约,才换得到试盒模型的资格**。
            // 除非有办法既在最终坐标系里消冲突、又保证最后的姿态规范,否则
            // 圆就是这里的正解。
            //
            // **竖排的那一档不进来。** 消冲突跑的时候坐标还在动,横竖之分那时
            // 判不了(`render::label_dir` 在这个阶段没有定义,与 `label_at` 同
            // 一个理由)。而且把竖排也算进来取最大,就是上面那条已经量过并
            // 否掉的"各向同性取最大"。实测 ACS 下 `NH` 的外接圆:横排 0.806 em、
            // 竖排 **0.836 em** —— 竖排几乎不改外接圆(+3.7%),它做的是把各向
            // 异性**转了 90°**。所以对竖排原子沿用横排值,最坏方向上只差 0.05
            // 个键长,比"取最大"那 2.3 个百分点便宜得多。
            [
                LabelPlace::Horizontal(HSide::Right),
                LabelPlace::Horizontal(HSide::Left),
            ]
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
fn better(now: Score, best: Score) -> bool {
    const EPS: f64 = 1e-9;
    if (now.0, now.1) != (best.0, best.1) {
        return (now.0, now.1) < (best.0, best.1);
    }
    now.2 < best.2 - EPS
}

/// 一次布局有多坏:`(画在同一点上的对数, 交叉的键对数, 碰撞深度平方和)`。
///
/// 越小越好,按字典序比。
type Score = (usize, usize, f64);

/// 打分。**越小越好**,按 `(画在同一点上的对数, 交叉的键对数, 碰撞深度)` 的
/// 字典序比较,见 [`better`]。
///
/// # 「重合」这一位是后补的,补它是因为**消冲突自己会造出重合**
///
/// 先前只有 `(交叉, 深度)`,交叉排第一。于是**一次消掉交叉的翻转会被接受,
/// 哪怕它把两个苯环叠成一个** —— 深度排在交叉后面,拦不住。
///
/// 实测坐实:全量语料最后剩的 8 处「原子不重合」违例里,1068/1069/5000 那
/// 6 处**布局阶段一对重合都没有**,是走完消冲突才冒出来的(整整六对,两个
/// 对位取代的苯环逐位叠上)。
///
/// 补上这一位之后次序与全库一致:**假环 > 看得见的丑 > 挤**。两个原子叠在
/// 一点会让图上凭空多一个环、读者没有任何办法看出它是假的;交叉只是难读。
///
/// **它在本语料上一次都没触发** —— 插桩数过,"只因这一位被拦下的翻转"是 0 次
/// (那 6 处重合是 [`crate::stereo::fix_cis_trans`] 造的,不是翻转造的)。留着
/// 不是因为量到了收益,是因为**次序不一致本身就是隐患**:布局一变,"拿假环换
/// 一处交叉"随时可能变成真事。判据 `a_phantom_ring_is_worse_than_any_number_of_crossings`
/// 直接钉这个次序,不靠语料。
///
/// 用字典序而不是加权求和:两者量纲不同,加权就要引入一个说不清的系数,而
/// 系数一变结论就变。字典序不需要系数。
///
/// 交叉排在深度前面是量出来的,不是拍的:把两者对调之后,全语料上有交叉的
/// 图从 381 涨到 415,而未解冲突一个没少。
fn score(mol: &MolBuilder, pos: &BTreeMap<u32, Point2>, radii: &[f64]) -> Score {
    /// 多近算"画在同一点上" —— 与硬判据 `原子不重合` 同一个阈值。
    const SAME: f64 = 0.05;
    let (pairs, crossings) = remaining(mol, pos, radii);
    let same = pairs
        .iter()
        .filter(|(i, j)| pos[i].dist(pos[j]) < SAME)
        .count();
    let depth: f64 = pairs
        .iter()
        .map(|(i, j)| {
            let want = radii[*i as usize] + radii[*j as usize];
            let d = pos[i].dist(pos[j]);
            (want - d).max(0.0).powi(2)
        })
        .sum();
    (same, crossings.len(), depth)
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

    (pairs, crossings(mol, pos))
}

/// 交叉的键对。只看两端都已放置的键。
///
/// # 为什么要单独拿出来
///
/// 布局有时跑在一份**改过拓扑的副本**上(η5 配位的多余键被摘掉了,见
/// `crate::hapto_extras`)。那时消冲突报的交叉少算了被摘掉的那些键 ——
/// 而它们照样会画出来。**画不好就要说出来**,所以 `generate` 在最终坐标上
/// 拿**原分子**再算一遍。
pub(crate) fn crossings(mol: &MolBuilder, pos: &BTreeMap<u32, Point2>) -> Vec<(u32, u32)> {
    let live: Vec<u32> = (0..mol.num_bonds())
        .map(|i| u32::try_from(i).expect("键数超出 u32"))
        .filter(|b| {
            let bd = &mol.bonds()[*b as usize];
            pos.contains_key(&bd.begin) && pos.contains_key(&bd.end)
        })
        .collect();
    let mut out = Vec::new();
    for (k, &b1) in live.iter().enumerate() {
        for &b2 in &live[k + 1..] {
            let (x, y) = (&mol.bonds()[b1 as usize], &mol.bonds()[b2 as usize]);
            if segments_cross(pos[&x.begin], pos[&x.end], pos[&y.begin], pos[&y.end]) {
                out.push((b1, b2));
            }
        }
    }
    out
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

/// 撑开候选的定序键。
///
/// `(结果重合数, 结果交叉数, −量化后的支点最窄角, 支点秩, 动的那一侧最小秩, 转向)`
/// —— 前三项是**等距不变的几何量**,后三项是规范量。取最小者。见 `relieve` 里
/// 用它的那一段:`depth` 是**故意不在里面**的。
type SplayKey = (usize, usize, i64, u32, u32, bool);

/// 一次撑开:把 `moved` 整体绕 `pivot` 转 `by`。
///
/// # 它与前三个算子不同:**会改变一个键角**
///
/// 键翻转、螺环翻转、度 4 置换都是等距变换,几何一点不损。这一个不是:
/// 支点处含被转那根键的几个角各变了 `by`。用它的理由是严重度次序 ——
/// 两个原子叠在一点,图上凭空多一个环,读者没有任何办法看出它是假的;
/// 一个键角从 120° 变成 90° 或 150° 只是不理想。**假环 > 读错结构 > 看得见的丑**。
///
/// # 但这不是"破例"
///
/// 本库**早就允许键角偏离理想值**:`chains::free_direction` 躲让时按
/// ±30° 一档铺开,最多铺到 5 档;硬判据 `键角不过窄` 的 89° 地板正是为它设的。
/// 撑开只是把同一把尺子搬到消冲突阶段用,而且**只走一档**,比布局阶段温和得多。
///
/// 放在这里而不是放进布局,还有一个量的理由:这里有 `best.0 > 0`(还有原子精确
/// 重合)这个闸门,触发面是全语料 4 个分子;放进布局就没有闸门,触发面是全语料。
///
/// # 它要付的账,如实记
///
/// 撑开会**把原本不交叉的键推成交叉**。加它那一次的全量语料实测:`有键交叉`
/// 86 → **90** —— 第 1068、1069 行的四张图从 0 处交叉变成有交叉。按严重度次序
/// 这笔账该付(假环是读者无从察觉的错,交叉只是难读),但不能不说。
///
/// **这一节的数是加撑开那一次的实测,不是当前值**(下面「同一次实测的另一面」
/// 那几个同样如此)。 后来「挑方向」那一步把交叉
/// 从硬门槛降成排序键第三位,基线整体挪过一次;当前值见 `harness/README.md`。
///
/// 同一次实测的另一面:`原子不重合` **8 → 2**(只剩第 7880 行那两张 —— 它能
/// 构造出 **8** 个候选,全部因"没消掉重合"被拒),`干净`
/// 17124 → **17126**,`有未解冲突` 187 → **185**,`键角不过窄` 违例仍是 8。
/// 改动半径 **6 张 / 17662**,正是那 3 个分子 × 2 套规范,别处一张没动。
///
/// # 几何后果只有一处,可以证明
///
/// 设 `R` 是绕支点 `p` 转 `θ`,只作用在 `moved` 上。`far_side` 返回的那一侧
/// 恒不含支点,所以 `p ∉ moved`;跨越 `moved` 与补集的键只有 `p–c` 一根
/// (否则 `far_side` 返 `None`,根本不出候选)。于是:
///
/// - **键长全保**:`moved` 内部同受 `R`;补集不动;`|p − R(c)| = |p − c|`
///   (`p` 是旋转中心)。
/// - **只有支点处的角变**:`moved` 内部是刚体变换;补集不动;**`c` 处的角也不变**
///   —— `c→x` 变成 `Rot_θ(x−c)`,而 `c→p = p − R(c) = Rot_θ(p−c)`,两个向量同转,
///   夹角不动。
struct Splay {
    pivot: u32,
    moved: Vec<u32>,
    by: f64,
    /// 转向:`false` = 背着参照邻居,`true` = 朝着它。见 [`splays`] ——
    /// 这个标法是**整体反射下等变**的,所以拿它当定序键与坐标系无关。
    toward: bool,
}

/// 撑开一档 30°,与布局用的方向网格同一个刻度。
const SPLAY: f64 = std::f64::consts::FRAC_PI_6;

/// 键角地板,与硬判据 `no_angle_is_pinched` 的 `FLOOR` 是同一个数。
const ANGLE_FLOOR_DEG: f64 = 89.0;

/// 浮点量化。阈值比较与定序都过它 —— 与 [`crate::chains`] 里 `pick` 同一个写法,
/// 理由也一样:裸比浮点会让"够不够 89°"取决于算到这一步的运算次序。
fn q6(x: f64) -> i64 {
    (x * 1e6).round() as i64
}

/// 撑开之后支点的最窄角**够不够格**。
///
/// 两条:不许跌破硬判据 `no_angle_is_pinched` 的地板(89°),**也不许比原来更窄**。
///
/// 后半句是为已经退化的布局留的 —— 那里的角本来就可能不到 89°,一刀切会把这个
/// 算子从那些图上整个挡掉,而它们恰恰最需要它。
///
/// 比之前先量化,与 [`crate::chains`] 里 `pick` 同一个写法:裸比浮点会让"够不够
/// 89°"取决于算到这一步的运算次序,而同一个分子的两种写法那个次序不同。
///
/// # 度数 ≥ 4 的支点会被它整个挡掉,这是想要的
///
/// 那些中心(金属)的角本来就密,任何一档 30° 都会跌破地板。而硬判据
/// `no_angle_is_pinched` **根本不查度数 ≥ 4 的原子** —— 也就是说在那里撑开
/// **拿不到判据上的任何收益,却要付真实的几何代价**。不做。
///
/// (别拿"金属的理想角本来就是 90°"去论证:`PairSwap` 的文档里写着
/// `chains::allocate` 在已占方向 ≥ 2 时是劈最大空隙,**十字没有保证** ——
/// 那会是一个没量过的说法。)
///
/// # 它在本语料上一次都没拦下过
///
/// 插桩数过(口径 `audit … large.smi 1`):全量语料共评估 **280** 个撑开候选,
/// 其中 **192** 个因"没消掉重合"被拒,余下 **88** 个全部合格。
/// `keeps_stereo` 判否 **0** 次;这条守卫判否 **4** 次,而那 4 个**同时也没消掉
/// 重合**(全在第 7880 行)—— 也就是说**它一次都没有单独拦下过任何候选**。
///
/// 留着不是因为量到了收益,是因为少了它这个算子就能把一个角压到 60°,那是
/// 读错结构。判据 `a_splay_may_not_pinch_the_pivot_below_the_floor` 直接钉这条,
/// 不靠语料 —— 拿分子去验会空过。
fn angle_survives(before: f64, after: f64) -> bool {
    // 写成"够得着地板 **或** 没比原来更窄"这个析取,而不是
    // `q6(after) >= q6(FLOOR.min(before))` —— 两者只在一处不同:地板那一支要
    // **比硬判据严一格**。判据是 `deg < 89.0` 即违例,而量化到 1e-6 会把
    // 88.9999995 归到 89.0 上,写成严格大于把这一格让出去,不然会放行一个判据当场要报的角。
    // 严一格的代价是最多多挡掉宽 1e-6 度的一条缝。
    //
    // **"不比原来更窄"这一支有个后果,记在这里:被转的那根键可以从支点的另一个
    // 邻居身上扫过去。** 原来 ∠(c,p,n) = a < 15° 时,朝 n 转 30° 之后是 30−a,
    // 只要 30−a ≥ a 就放行 —— 而这时 c 与 n 的先后已经反了。构型不会画错(楔形
    // 是走完全程之后按最终坐标重新指派的,外部判官 496 一致 / 0 读反没动),
    // 但图上会出现两根键几乎叠在一起。语料上不触发:这条守卫拦下 0 次,而所有
    // 合格候选的 `before` 全部 ≥ 90°。
    q6(after) > q6(ANGLE_FLOOR_DEG) || q6(after) >= q6(before)
}

/// `p` 处最窄的键角(度)。没有两个已放置的邻居时给 180°(无角可窄)。
///
/// **跳过三元环内角** —— 三根键等长时它必然恰好 60°,与 `键长全等` 不可兼得,
/// 硬判据 `no_angle_is_pinched` 那边同样放行,两边口径要一致。
fn narrowest_at(mol: &MolBuilder, pos: &BTreeMap<u32, Point2>, p: u32) -> f64 {
    let nb: Vec<u32> = mol
        .neighbors(p)
        .map(|(n, _)| n)
        .filter(|n| pos.contains_key(n))
        .collect();
    let Some(c) = pos.get(&p) else {
        return 180.0;
    };
    let mut worst = 180.0f64;
    for i in 0..nb.len() {
        for j in (i + 1)..nb.len() {
            if mol.neighbors(nb[i]).any(|(n, _)| n == nb[j]) {
                continue; // 三元环内角
            }
            let u = (pos[&nb[i]] - *c).normalized();
            let v = (pos[&nb[j]] - *c).normalized();
            worst = worst.min(u.dot(v).clamp(-1.0, 1.0).acos().to_degrees());
        }
    }
    worst
}

/// 这个原子**本来就该画成直的**:有三键,或者有两根双键(sp 中心)。
///
/// 撑开不许拿它当支点 —— 把炔或累积双键画成弯的是**读错结构**,而 89° 那道
/// 地板拦不住(150° 远在地板之上)。语料上不触发,留着是因为这一类陷阱要钉住。
fn is_linear_centre(mol: &MolBuilder, a: u32) -> bool {
    // **度数必须是 2。** 少了这一条,砜/硫酸酯的硫会被算进来:它有两根 S=O
    // 而度数是 4,是四面体,不是 sp。实测 `CS(C)(=O)=O`、`O=S(=O)(O)O`
    // 的硫都会被误判。方向上那是保守的(只多挡候选,而度 4 的支点本来也会被
    // 角地板挡掉),但按错的理由排除,哪天出现度 3 的双双键中心就会静默漏判。
    if mol.degree(a) != 2 {
        return false;
    }
    let mut doubles = 0;
    let mut triple = false;
    for (_, bi) in mol.neighbors(a) {
        match mol.bonds()[bi as usize].order {
            BondOrder::Triple => triple = true,
            BondOrder::Double => doubles += 1,
            _ => {}
        }
    }
    triple || doubles >= 2
}

/// 枚举所有可做的撑开。
///
/// 要动的那一侧由 [`far_side`] 给(更少的一边,平局按最小规范秩),支点是
/// **不在那一侧**的那个端点。
///
/// # 键集**不能**直接照搬 [`rotatable`]
///
/// 那个函数排除端点键,理由写在它自己的文档里:「端点键翻了等于什么都没做 ——
/// 那一侧只有它自己,镜像回原位」。**那条理由只对镜像成立。** 绕邻居转 30°
/// 不是对合变换,那个端基真的移动了(弦长 2·sin15° ≈ 0.518 个键长)。
///
/// 而且这不是纸上谈兵:`键角不过窄` 剩的那 11 处窄角,支点**全都挂着端点键**
/// —— 照搬 `rotatable` 的话它们**一处都消不掉**(实测 0/11);放开之后能消
/// 4 处(另有 2 处只能靠制造一处重合来消,那是拿假环换窄角,该拒)。
///
/// 所以这里自己筛:不在环里、两端都已放置、支点度数 ≥ 2。环那一条仍然要 ——
/// 环上的键转了会把环撕开,而 `far_side` 本来就对它返 `None`。
///
/// **度数那一条是冗余的早退,如实记着**:度 1 的支点只有 `root` 一个邻居,
/// 下面 `others` 必空、`sign` 的 `find` 返 `None`,照样 `continue`。变异验过 ——
/// 删掉它判据全绿、语料改动半径 0。留着只因为它便宜,不是因为少了它会出错。
///
/// **放开端点键在本语料上一张图都没改**(实测改动半径 0 / 17662):闸门只在
/// 4 个分子上开,那几个的赢家没换。改它不是因为量到了收益,是因为照搬那条
/// 排除等于**按一个不成立的理由**把一整类候选丢掉。
///
/// # 两个转向怎么定名,以及为什么必须这么定
///
/// "顺时针 30°"是在**全局坐标系**里说的话。而消冲突之前的那个坐标系,
/// **今天确实是逐写法一致的,靠的是构造**:起手环落在
/// `regular_polygon(n, 0.0)` 上(起始角是**常数**,环上原子序由规范秩定),
/// 无环分量的种子按规范秩取、落在原点。
///
/// **但没有任何判据守着这一条。** 端到端那条 `写法无关` 排在
/// [`crate::orient::canonicalise`] **之后**,而 `canonicalise` 会把姿态(含镜像)
/// 归一 —— 前置坐标系哪天被改成镜像,现有判据**一条都不会红**,而带符号的
/// ±30° 会当场分岔。依赖一个没有判据守着的前提,与判据空过是同一件事。
///
/// 所以两个转向按**内蕴**的方式命名:取支点除 `root` 之外规范秩最小、且与
/// `root` 不共线的那个邻居 `n`,把两个转向标成「朝 n」与「背 n」。整体反射会
/// 同时翻转 `cross` 的符号和"朝 n"对应的转向 —— 这个标法在反射下**等变**,
/// 于是与坐标系无关。
///
/// **这一条只有一个判据钉得住,而且必须是那一个。** 实测:把它改回固定的
/// `+30° 在先`,全量语料写法无关仍是 0 违例、端到端那条写法判据(含 1068 的
/// 三种写法)照样绿 —— 正因为前提今天成立。红的只有
/// `the_two_splay_directions_are_named_without_looking_at_the_canvas`:
/// 它把坐标整体反射一次,直接验等变性,不经过 `canonicalise`。
///
/// 全共线时(找不到这样的 `n`)两个转向分不出规范的先后,**整个候选跳过**:
/// 宁可不做,也不留一处依赖写法的分岔。
fn splays(mol: &MolBuilder, pos: &BTreeMap<u32, Point2>, ranks: &[u32]) -> Vec<Splay> {
    let mut out = Vec::new();
    for b in 0..u32::try_from(mol.num_bonds()).expect("键数超出 u32") {
        let bd = &mol.bonds()[b as usize];
        if bd.flags.contains(BondFlags::IN_RING)
            || !pos.contains_key(&bd.begin)
            || !pos.contains_key(&bd.end)
        {
            continue;
        }
        let Some(moved) = far_side(mol, pos, b, ranks) else {
            continue;
        };
        let (pivot, root) = if moved.contains(&bd.end) {
            (bd.begin, bd.end)
        } else {
            (bd.end, bd.begin)
        };
        if mol.degree(pivot) < 2 || is_linear_centre(mol, pivot) {
            continue;
        }
        let c = pos[&pivot];
        let u = pos[&root] - c;
        let mut others: Vec<u32> = mol
            .neighbors(pivot)
            .map(|(n, _)| n)
            .filter(|n| *n != root && pos.contains_key(n))
            .collect();
        others.sort_by_key(|n| (ranks[*n as usize], *n));
        let Some(sign) = others
            .iter()
            .map(|n| u.cross(pos[n] - c))
            .find(|x| x.abs() > 1e-9)
            .map(f64::signum)
        else {
            continue; // 支点周围全与 `root` 共线 —— 两个转向定不出规范的先后
        };
        for toward in [false, true] {
            let by = if toward { sign * SPLAY } else { -sign * SPLAY };
            out.push(Splay {
                pivot,
                moved: moved.clone(),
                by,
                toward,
            });
        }
    }
    out
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

/// 断开键 `b` 之后该翻的那一侧:**原子更少**的那一边,一样多时按**最小规范秩**。
///
/// 返回 `None` 表示这根键其实在环上(断开后两端仍连通),翻不得。环感知标记
/// 之外再判一次是**故意的**:标记来自净化,而调用方未必净化过。
///
/// # 为什么必须取更少的那一侧(以及平局时为什么不能看 `end`)
///
/// 哪一端是 `end` 依**写法**而定。取 `end` 那一侧的话,同一根化学键在一种写法
/// 里镜像的是一个甲基、在另一种写法里镜像的是整个苯环 —— 两者相差一次全局
/// 反射,而接受的翻转**次数**也会跟着不同,最后坐标就对不上。
///
/// 实测:阿司匹林的两种写法,布局阶段已经完全一致了,却在这里分岔 ——
/// 一种翻两次、另一种翻一次。
///
/// 取更少的那一侧也更自然:不该为了挪一个甲基把整个分子翻过来。
fn far_side(
    mol: &MolBuilder,
    pos: &BTreeMap<u32, Point2>,
    b: u32,
    ranks: &[u32],
) -> Option<Vec<u32>> {
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

    // 取更少的那一侧。**平局时按这一侧最小的规范秩定。**
    //
    // 先前这里只有"更少的那一侧"这一条,平局就直接返回 `out` —— 而 `out` 是从
    // `bd.end` 走出来的,`begin`/`end` 谁在前正是书写痕迹。注释当时已经写着
    // "平局按最小规范秩定",**可代码里没有那一步**。
    //
    // # 一次翻错侧就够了,`orient` 兜不住
    //
    // 绕同一根轴翻这侧还是那侧,两套坐标差**一次整体反射**。直觉上
    // [`crate::orient::canonicalise`] 该能归一它 —— 它的 24 个候选姿态里有镜像。
    // **但那是 D₁₂,镜面只落在 15° 的整数倍上**:绕角度 φ 的轴做整体反射,归一
    // 之后的残差是 `2φ mod 30°`。轴不在 15° 网格上,残差就消不掉。
    //
    // 实测(临时判据,拿阿司匹林量):
    //
    // | 轴 | 0° | 15° | 30° | 45° | 7° | 22.5° | 84° | 96° |
    // |---|---|---|---|---|---|---|---|---|
    // | 残差 | 0 | 0 | 0 | 0 | **16°** | **15°** | **12°** | **18°** |
    //
    // 语料第 6458 行正是这一支:两种写法**各只翻了 1 根键**,布局阶段完全一致
    // (都是 −14.6099°),消冲突后一个 −91.5135°、一个 −76.4865° —— 关于
    // **−84.0° 对称**,是镜像不是旋转;2×(−84) ≡ 12 (mod 30),最终图上差的
    // 正是 12°。而稠环的布局姿态本来就不在 30° 网格上,所以这不是罕见情形。
    //
    // # 前提:`ranks` 是单射
    //
    // 两侧不相交,所以秩单射时两边的最小值必不相等,`<` 一定分得出胜负。
    // 秩若有重复,`<` 为假会**静默退回 `out`** —— 又回到写法依赖,且不报信号。
    // 生产上的 [`crate::ranks_of`] 造的是 `0..n` 的双射;实测全量语料
    // 136439 次调用里重复 0 次、两侧最小秩相等 0 次。
    debug_assert!(
        {
            let mut seen: BTreeSet<u32> = BTreeSet::new();
            out.iter()
                .chain(other.iter())
                .all(|a| seen.insert(ranks[*a as usize]))
        },
        "规范秩有重复,平局会静默退回存储序"
    );
    let smallest = |v: &[u32]| v.iter().map(|a| ranks[*a as usize]).min();
    // `other` 恒非空(它从 `blocked` 起头),`Ordering::Equal` 那一支因此总取得到
    // 两个最小值。实测 136439 次调用里 `other` 为空 0 次。
    let take_other = match out.len().cmp(&other.len()) {
        std::cmp::Ordering::Greater => true,
        std::cmp::Ordering::Equal => smallest(&other) < smallest(&out),
        std::cmp::Ordering::Less => false,
    };
    Some(if take_other { other } else { out })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::layout;

    /// 解析 + 净化 + **顺反感知**。第三步不在净化的 12 步里,漏了的话每根双键的
    /// `stereo` 都是 `None`,顺反校正整个空转 —— 顺式反式画成同一张图而判据照样
    /// 绿。[`generate`] 的 `debug_assert!` 现在会当场拦住这种输入。
    fn prep(smi: &str) -> MolBuilder {
        let mut m = omgkit_io::smiles::parse(smi).unwrap();
        omgkit_chem::pipeline::sanitize(&mut m).unwrap();
        omgkit_io::stereo::perceive_bond_stereo(&mut m);
        m
    }

    fn laid(smi: &str, style: &Style) -> (MolBuilder, BTreeMap<u32, Point2>, Report) {
        let m = prep(smi);
        let ranks = omgkit_io::canon::canonical_ranks(&m);
        let mut pos: BTreeMap<u32, Point2> = BTreeMap::new();
        for p in layout::layout_all(&m, &ranks, style, None) {
            pos.extend(p.pos);
        }
        let rep = relieve(&m, &mut pos, &ranks, style);
        (m, pos, rep)
    }

    /// 布局 → **顺反校正** → 消冲突。[`laid`] 少了顺反那一步,而撑开正是为它
    /// 服务的 —— 那 6 处重合是 [`crate::stereo::fix_cis_trans`] 造出来的,
    /// 布局摆开时一对都没有。秩用 `crate::ranks_of` 而不是 `canonical_ranks`,
    /// 与生产同源。
    ///
    /// # 它**不是**完整的流水线,少的这几步限定了能拿什么分子来验
    ///
    /// 生产(`generate_with`)还有:补立体氢、配位键摊平、η 摘细,以及**把各个
    /// 连通分量按包围盒并排摆开**。这里是 `pos.extend`,所有分量原样堆在原点。
    ///
    /// 所以**只能拿单分量、无配位键、无 η、不需要补氢的分子来验**。往下面几条
    /// 判据里加盐或配合物会踩坑:堆在原点会凭空造出重合,于是凭空造出撑开,
    /// 而判据看着还是绿的。现用的那些分子都满足这个条件(已逐个核对),而且
    /// 语料级的改动半径反证了生产确实对第 1068、1069、5000 行撑开了。
    fn posed_with_cis_trans(
        smi: &str,
        style: &Style,
    ) -> (MolBuilder, BTreeMap<u32, Point2>, Vec<u32>) {
        let m = prep(smi);
        let ranks = crate::ranks_of(&m);
        let mut pos: BTreeMap<u32, Point2> = BTreeMap::new();
        for p in layout::layout_all(&m, &ranks, style, None) {
            pos.extend(p.pos);
        }
        let mut flat = vec![Point2::ORIGIN; m.num_atoms()];
        for (a, q) in &pos {
            flat[*a as usize] = *q;
        }
        crate::stereo::fix_cis_trans(&m, &mut flat, &ranks);
        for (a, q) in pos.iter_mut() {
            *q = flat[*a as usize];
        }
        (m, pos, ranks)
    }

    fn laid_with_cis_trans(
        smi: &str,
        style: &Style,
    ) -> (MolBuilder, BTreeMap<u32, Point2>, Report) {
        let (m, mut pos, ranks) = posed_with_cis_trans(smi, style);
        let rep = relieve(&m, &mut pos, &ranks, style);
        (m, pos, rep)
    }

    /// 顺反规格把两个环钉到同一个位置上时,撑开要把它们拉开。
    ///
    /// 这三个分子是全量语料最后剩下的顺反类重合(第 1068、1069、5000 行)。
    /// 它们的共同构型是**两个大取代基挂在相邻的环原子上**:外向键只岔开 60°,
    /// 而顺反规格把两个 ipso 位置钉死在相距 1.0 的地方 —— 正是同一个六边形的
    /// 相邻顶点,于是两个对位取代的环逐位叠上。
    ///
    /// 理想几何下这**画不出来**,三条出路都堵死:换镜像的另一侧(两侧互为整体
    /// 反射,等距,距离逐一相同)、翻转那个环(单点连接的六元环关于该轴自对称,
    /// 翻了等于没翻)、改布局次序(位置由几何加规格联合钉死)。所以只剩键角。
    #[test]
    fn two_rings_that_cis_trans_forces_onto_each_other_get_pulled_apart() {
        let mut splayed_any = false;
        for smi in [
            "C1(/C(NC2=C(N1)C=CC=C2)=N\\C3=CC=C(C(=O)OCC)C=C3)=N/C4=CC=C(C(=O)OCC)C=C4",
            "C1(/C(NC2=C(N1)C=CC=C2)=N\\C3C=CC(=CC=3)OCC)=N/C4=CC=C(C=C4)OCC",
            "CCn\\1ccnc(/c1=N\\c2ccc(cc2)OC)N3CCOCC3",
        ] {
            for style in &Style::ALL {
                let (m, pos, rep) = laid_with_cis_trans(smi, style);
                splayed_any |= !rep.splayed.is_empty();
                let atoms: Vec<u32> = pos.keys().copied().collect();
                for (k, &i) in atoms.iter().enumerate() {
                    for &j in &atoms[k + 1..] {
                        assert!(
                            pos[&i].dist(pos[&j]) >= 0.05,
                            "[{}] {smi}:原子 {i} 与 {j} 还画在同一点上",
                            style.name
                        );
                    }
                }
                // 撑开不许拿顺反去换:那是把"画不好"换成"画错了"
                let mut flat = vec![Point2::ORIGIN; m.num_atoms()];
                for (a, q) in &pos {
                    flat[*a as usize] = *q;
                }
                assert!(
                    crate::stereo::cis_trans_intact(&m, &flat),
                    "[{}] {smi}:撑开之后顺反不对了",
                    style.name
                );
                // 也不许把角压到硬判据的地板以下
                for a in pos.keys() {
                    let deg = narrowest_at(&m, &pos, *a);
                    assert!(
                        q6(deg) >= q6(ANGLE_FLOOR_DEG),
                        "[{}] {smi}:原子 {a} 处的键角被压到 {deg:.1}°",
                        style.name
                    );
                }
            }
        }
        // **防空过。** 这三个分子里没有一个真撑开过的话,上面那些断言全是恒真的
        // —— 它们描述的本来就是"没毛病"。
        assert!(splayed_any, "一次撑开都没发生,这条判据在空过");
    }

    /// 撑开**只改支点那一处的角**,别处一个都不动;键长一根不变。
    ///
    /// 这是这个算子唯一的几何代价,必须钉死它就这么多 —— 论证见 [`Splay`],
    /// 这里逐原子把角的多重集比一遍,不靠论证。
    #[test]
    fn a_splay_only_changes_the_angle_at_its_pivot() {
        let angles = |m: &MolBuilder, pos: &BTreeMap<u32, Point2>, a: u32| -> Vec<i64> {
            let nb: Vec<u32> = m.neighbors(a).map(|(n, _)| n).collect();
            let c = pos[&a];
            let mut v = Vec::new();
            for i in 0..nb.len() {
                for j in (i + 1)..nb.len() {
                    let u = (pos[&nb[i]] - c).normalized();
                    let w = (pos[&nb[j]] - c).normalized();
                    v.push(q6(u.dot(w).clamp(-1.0, 1.0).acos()));
                }
            }
            v.sort_unstable();
            v
        };
        let lengths = |m: &MolBuilder, pos: &BTreeMap<u32, Point2>| -> Vec<i64> {
            let mut v: Vec<i64> = m
                .bonds()
                .iter()
                .map(|b| q6(pos[&b.begin].dist(pos[&b.end])))
                .collect();
            v.sort_unstable();
            v
        };
        let mut tried = 0usize;
        for smi in ["CCc1ccccc1", "CC(C)COC(C)C", "OC(=O)c1ccccc1OC(C)=O"] {
            let (m, pos, _) = laid(smi, &Style::ACS_1996);
            let ranks = crate::ranks_of(&m);
            for sp in splays(&m, &pos, &ranks) {
                tried += 1;
                let mut after = pos.clone();
                let c = pos[&sp.pivot];
                for a in &sp.moved {
                    let p = pos[a].rotated_about(c, sp.by);
                    after.insert(*a, p);
                }
                assert_eq!(
                    lengths(&m, &pos),
                    lengths(&m, &after),
                    "{smi}:绕 {} 撑开改了键长",
                    sp.pivot
                );
                for a in pos.keys() {
                    if *a == sp.pivot {
                        continue;
                    }
                    assert_eq!(
                        angles(&m, &pos, *a),
                        angles(&m, &after, *a),
                        "{smi}:绕 {} 撑开,却把原子 {a} 处的键角也改了",
                        sp.pivot
                    );
                }
            }
        }
        assert!(tried > 0, "一个撑开候选都没枚举出来,这条判据在空过");
    }

    /// 端点键**要**进撑开的候选 —— 转一个甲基不是空操作。
    ///
    /// [`rotatable`] 把端点键排除在外,理由是"翻了等于没翻";那对镜像成立,
    /// 对旋转不成立。照搬它的话,支点挂着端基的那些位置一个候选都构造不出来,
    /// 而实测语料里剩下的窄角**恰恰都在那种位置上**。
    ///
    /// 这条钉的是"候选集里确实有一个动的只是端基"。变异:把枚举换回
    /// `for b in rotatable(mol, pos)`,当场红。
    #[test]
    fn a_terminal_bond_is_a_real_splay_candidate() {
        let (m, pos, _) = laid("CCc1ccccc1", &Style::ACS_1996);
        let ranks = crate::ranks_of(&m);
        let cands = splays(&m, &pos, &ranks);
        assert!(!cands.is_empty(), "一个候选都没有,这条判据在空过");
        // 甲基那根键:动的一侧只有它自己,支点是它的邻居
        let terminal = cands.iter().find(|sp| {
            sp.moved.len() == 1 && m.degree(sp.moved[0]) == 1 && m.degree(sp.pivot) >= 2
        });
        assert!(
            terminal.is_some(),
            "候选里没有一个是「只动一个端基」—— 端点键被整类丢掉了"
        );
        // 而且它真的动得了:转过去之后那个端基不在原处
        let sp = terminal.expect("上面已经断言过有");
        let c = pos[&sp.pivot];
        let moved = pos[&sp.moved[0]].rotated_about(c, sp.by);
        assert!(
            moved.dist(pos[&sp.moved[0]]) > 0.1,
            "端基转了 30° 却几乎没动 —— 那才真是空操作"
        );
    }

    /// 撑开不许拿**本来就该画成直的**原子当支点。
    ///
    /// 炔、累积双键的中心画弯了是**读错结构**,而 89° 那道地板拦不住它
    /// (撑开之后是 150°,远在地板之上)。
    ///
    /// # 为什么要手摆坐标
    ///
    /// 拿布局出来的坐标验这条是**空过的**:布局把炔画直了,于是 sp 中心的另一个
    /// 邻居与 `root` 共线,`splays` 里那条"共线就跳过"先一步把它挡掉了 ——
    /// 去掉 `is_linear_centre` 判据照样绿(实测,变异验过)。
    ///
    /// 可 `pos` 是**外部传进来的**,守卫不能依赖布局的好意。所以这里手摆一份把
    /// 炔中心画弯的坐标:那时共线规则拦不住,只剩这条守卫。下面先证明"共线拦不住"
    /// 确实成立,再验守卫 —— 否则这条判据又空过了。
    ///
    /// (只挪一个原子,别处的键长会不成样子 —— 无所谓:`splays` 的筛选是拓扑加
    /// 一次叉乘,不看键长。)
    #[test]
    fn a_splay_never_pivots_on_an_atom_that_should_be_drawn_straight() {
        let mut exposed = 0usize;
        for smi in ["CC#CCC", "CCC#CC#CCC", "CC=C=CCC"] {
            let (m, pos0, _) = laid(smi, &Style::ACS_1996);
            let ranks = crate::ranks_of(&m);
            let mut pos = pos0.clone();
            // 把每个 sp 中心的一个邻居绕它转 40°,让那里不再是直的
            for a in 0..u32::try_from(m.num_atoms()).unwrap() {
                if !is_linear_centre(&m, a) {
                    continue;
                }
                if let Some((n, _)) = m.neighbors(a).next() {
                    let c = pos[&a];
                    let q = pos[&n].rotated_about(c, 40f64.to_radians());
                    pos.insert(n, q);
                }
            }
            // **先证明共线那条已经拦不住了** —— 不然下面的断言是恒真的
            for b in rotatable(&m, &pos) {
                let Some(moved) = far_side(&m, &pos, b, &ranks) else {
                    continue;
                };
                let bd = &m.bonds()[b as usize];
                let (pivot, root) = if moved.contains(&bd.end) {
                    (bd.begin, bd.end)
                } else {
                    (bd.end, bd.begin)
                };
                if !is_linear_centre(&m, pivot) {
                    continue;
                }
                let c = pos[&pivot];
                let u = pos[&root] - c;
                if m.neighbors(pivot)
                    .map(|(n, _)| n)
                    .filter(|n| *n != root && pos.contains_key(n))
                    .any(|n| u.cross(pos[&n] - c).abs() > 1e-9)
                {
                    exposed += 1;
                }
            }
            for sp in splays(&m, &pos, &ranks) {
                assert!(
                    !is_linear_centre(&m, sp.pivot),
                    "{smi}:撑开拿 sp 原子 {} 当了支点,会把它画弯",
                    sp.pivot
                );
            }
        }
        assert!(
            exposed > 0,
            "没有一个 sp 支点是共线规则拦不住的 —— 这条判据在空过"
        );
    }

    /// 撑开不许把支点的角压到硬判据的地板以下,也不许比原来更窄。
    ///
    /// **拿分子验这条会空过** —— 语料上这道守卫一次都没拦下过任何候选
    /// (最窄角全部 ≥ 90°)。所以直接钉这个谓词,像
    /// [`a_phantom_ring_is_worse_than_any_number_of_crossings`] 那样。
    #[test]
    fn a_splay_may_not_pinch_the_pivot_below_the_floor() {
        assert!(!angle_survives(120.0, 60.0), "把 120° 压成 60° 被放行了");
        assert!(!angle_survives(120.0, 88.9), "压到地板下面一点也不许");
        // **量化的那半格**:硬判据是 `deg < 89.0` 即违例,而 q6 会把 88.9999996
        // 归到 89.0 上。守卫不比判据严一格的话,这个角会被放行,而判据当场要报。
        assert!(
            !angle_survives(120.0, 88.999_999_6),
            "88.9999996° 被量化抹平成 89° 放行了 —— 硬判据那边是违例"
        );
        assert!(
            angle_survives(120.0, 90.0),
            "120° 变 90° 该放行 —— 地板是 89°"
        );
        assert!(angle_survives(120.0, 150.0), "撑宽了反倒被挡下");
        // 已经退化的布局:角本来就不到地板,只要不更窄就该放行 —— 一刀切会把
        // 这个算子从最需要它的那些图上整个挡掉
        assert!(angle_survives(70.0, 70.0), "本来就 70°、没更窄,却被挡下");
        assert!(angle_survives(70.0, 100.0), "从 70° 撑到 100° 反倒被挡下");
        assert!(!angle_survives(70.0, 69.0), "本来 70° 又压窄了,该挡下");
    }

    /// 两个转向的命名**不看画布**:整套坐标做一次整体反射,候选跟着等变。
    ///
    /// # 为什么这条必须有
    ///
    /// "顺时针 30°"是在全局坐标系里说的话,而 [`crate::orient::canonicalise`]
    /// 排在消冲突**之后** —— 消冲突之前的姿态一旦哪天变成镜像,带符号的 ±30°
    /// 当场分岔,而端到端那条 `写法无关` 判据**一条都不会红**(它排在
    /// `canonicalise` 之后,姿态已经归一)。依赖一个没有判据守着的前提,
    /// 与判据空过是同一件事。
    ///
    /// 所以这里直接验等变性:把 `pos` 整体反射一次(那是同一张图的另一种坐标
    /// 表示),`toward` 相同的候选转出来的几何,必须与反射前那个候选的反射像
    /// 逐点相同。
    #[test]
    fn the_two_splay_directions_are_named_without_looking_at_the_canvas() {
        let flip = |p: Point2| Point2::new(p.x, -p.y);
        let mut pairs = 0usize;
        for smi in ["CCc1ccccc1", "CC(C)COC(C)C", "OC(=O)c1ccccc1OC(C)=O"] {
            let (m, pos, _) = laid(smi, &Style::ACS_1996);
            let ranks = crate::ranks_of(&m);
            let mirrored: BTreeMap<u32, Point2> = pos.iter().map(|(a, p)| (*a, flip(*p))).collect();
            let here = splays(&m, &pos, &ranks);
            let there = splays(&m, &mirrored, &ranks);
            assert!(!here.is_empty(), "{smi} 一个候选都没有,这一档验不了");
            assert_eq!(here.len(), there.len(), "{smi}:反射之后候选数变了");
            for (a, b) in here.iter().zip(&there) {
                assert_eq!(
                    (a.pivot, &a.moved, a.toward),
                    (b.pivot, &b.moved, b.toward),
                    "{smi}:反射之后候选的次序或标法变了"
                );
                pairs += 1;
                let (ca, cb) = (pos[&a.pivot], mirrored[&b.pivot]);
                for x in &a.moved {
                    let got = mirrored[x].rotated_about(cb, b.by);
                    let want = flip(pos[x].rotated_about(ca, a.by));
                    assert!(
                        got.dist(want) < 1e-9,
                        "{smi}:绕 {} 撑开,反射之后落点对不上 —— 转向的命名看了画布",
                        a.pivot
                    );
                }
            }
        }
        assert!(pairs > 0, "一对都没比,这条判据在空过");
    }

    /// 撑开**只花在原子精确重合上**,不许拿它去消交叉或消挤压。
    ///
    /// 这是全模块唯一要付几何代价的算子,闸门就是它的全部纪律:交叉只是难读,
    /// 假环是读者无从察觉的错。闸门松掉的话,一堆本来只是挤一点的图会白白被
    /// 改掉一个键角。
    ///
    /// # 前提要自己成立,而且这一条被判据审核抓过一次
    ///
    /// 先前挑的四个分子(阿司匹林等)进 `relieve` 时打分**就是 `(0, 0, 0)`**,
    /// 于是函数在入口那句 `if best == (0, 0, 0.0) { return }` 直接返回,
    /// **算子循环一次都没跑** —— `splayed.is_empty()` 是恒真的,把闸门变异成
    /// "有交叉就开"这条判据也不红。典型的空过。
    ///
    /// 所以现在两条前提都当场断言:
    ///
    /// 1. 进消冲突时**没有任何一对原子叠在一点**(`entry.0 == 0`)—— 闸门不该开;
    /// 2. 打分**不是** `(0, 0, 0)` —— 否则消冲突在入口就返回了,下面白断言。
    ///
    /// 分子取自真实语料(挑的是"有挤压、没重合"那一档),不是自己造的。
    #[test]
    fn a_splay_is_only_spent_on_atoms_drawn_on_top_of_each_other() {
        let mut checked = 0usize;
        // 三条都取自真实语料,而且**两套规范下入口打分都非平凡**(这一点是
        // 断言出来的,不是挑的时候拍的)。第三条尤其要紧:它入口就带着一处
        // **键交叉**,闸门照样不许开 —— 交叉不是撑开该管的事。
        for smi in [
            "C1(C(N2C=CC=CC(=NC=1O)2)=O)SC(=S)N(C)C",
            "c1ccc(c(c1)/C=N\\[C@@H]2CONC2=O)O",
            "C1CN[Ni]23(N1)(NCCN2)NCCN3",
            "[O-][N+](=O)C1=CC(=C(NC2=C(C=C(C=C2[N+]([O-])=O)[N+]([O-])=O)\
             [N+]([O-])=O)C(=C1)[N+]([O-])=O)[N+]([O-])=O",
        ] {
            for style in &Style::ALL {
                let (m, mut pos, ranks) = posed_with_cis_trans(smi, style);
                let radii = radii(&m, style);
                let entry = score(&m, &pos, &radii);
                assert_eq!(
                    entry.0, 0,
                    "[{}] {smi} 进消冲突时就有原子叠着 —— 闸门本来就该开,选错例子了",
                    style.name
                );
                assert!(
                    entry != (0, 0, 0.0),
                    "[{}] {smi} 进消冲突时打分是 (0,0,0),消冲突会在入口返回 —— \
                     算子循环一次都不跑,这条判据是空过的",
                    style.name
                );
                checked += 1;
                let rep = relieve(&m, &mut pos, &ranks, style);
                assert!(
                    rep.splayed.is_empty(),
                    "[{}] {smi}:没有一对原子叠着,却撑开了 {:?}",
                    style.name,
                    rep.splayed
                );
            }
        }
        assert!(checked > 0, "一个例子都没查,这条判据在空过");
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
        // **先前这里断言"至少翻转过一根键"** —— 那是防判据空过的:重合本来就不
        // 发生的话,上面那句是恒真的。
        //
        // 现在不该再这么断言了,因为**修法搬到了上游**:`chains::place_neighbours`
        // 在"只有一个已占方向"时会比较两侧的拥挤度,乙酰基那条臂不再朝苯环卷
        // 回去,两个羰基氧从一开始就不落在同一点上,`relieve` 无事可做。
        //
        // 防空过换一个更强的说法:**布局阶段(消冲突之前)就已经没有重合**。
        // 它比"翻过一根键"更贴近现在的事实,而且一旦有人把上游那个改动去掉,
        // 这一条会立刻红。
        let mut before: BTreeMap<u32, Point2> = BTreeMap::new();
        for p in layout::layout_all(
            &m,
            &omgkit_io::canon::canonical_ranks(&m),
            &Style::ACS_1996,
            None,
        ) {
            before.extend(p.pos);
        }
        for i in 0..n {
            for j in (i + 1)..n {
                let (a, b) = (i as u32, j as u32);
                assert!(
                    before[&a].dist(before[&b]) > 0.3,
                    "布局阶段原子 {a} 与 {b} 就重合了 —— 上游那个「挑空的一侧」没起作用"
                );
            }
        }
        let _ = &rep;
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
                for p in layout::layout_all(&m, &ranks, style, None) {
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
                for p in layout::layout_all(&m, &ranks, style, None) {
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
                for p in layout::layout_all(&m, &ranks, style, None) {
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
        let mut checked = 0usize;
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
                checked += 1;
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
        // 改写之后解析/sanitize 失败、或者规范 SMILES 对不上,都会整段跳过 ——
        // 一次都没比成时这条判据是空过的
        assert!(checked > 0, "一次都没比成,判据空过了");
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
    fn which_side_gets_flipped_does_not_depend_on_how_it_was_written() {
        // **`far_side` 平局那一档的直接判据。** 两侧一样多时翻哪一侧,先前由
        // `begin`/`end` 说了算 —— 而谁在前正是书写痕迹。函数注释当时已经写着
        // "平局按最小规范秩定",**可代码里没有那一步**。
        //
        // 把每根键映成「两端的规范秩 → 返回那一侧的秩集合」。这张表与原子编号
        // 无关,所以同一个分子的不同写法必须给出同一张表。
        //
        // **秩要用生产上那个 `ranks_of`**,不是 `canonical_ranks` —— 后者的深层
        // 平局是任取的(见 `crate::ranks_of`),拿它当判据会把两种毛病混在一起。
        //
        // 分子从两处来:正丁烷(最小的两侧等大的例子)与**语料里现成的四个**。
        // 变异(平局退回 `out`)在每一个上都红。
        type Table = std::collections::BTreeMap<(u32, u32), BTreeSet<u32>>;
        let table = |smi: &str| -> Table {
            let m = prep(smi);
            let ranks = crate::ranks_of(&m);
            let mut pos: BTreeMap<u32, Point2> = BTreeMap::new();
            for p in layout::layout_all(&m, &ranks, &Style::ACS_1996, None) {
                pos.extend(p.pos);
            }
            let mut t = Table::new();
            for b in 0..u32::try_from(m.num_bonds()).unwrap() {
                let Some(side) = far_side(&m, &pos, b, &ranks) else {
                    continue;
                };
                let bd = &m.bonds()[b as usize];
                let (x, y) = (ranks[bd.begin as usize], ranks[bd.end as usize]);
                t.insert(
                    (x.min(y), x.max(y)),
                    side.iter().map(|a| ranks[*a as usize]).collect(),
                );
            }
            t
        };
        for ws in [
            // 正丁烷:中间那根键两边各两个碳
            vec!["CCCC", "C(CC)C"],
            // 语料里现成的(4881 / 4824 / 2631 / 4719 行)
            vec!["CCCOS(O)(=O)=O", "OS(=O)(=O)OCCC"],
            vec!["CCCCS(O)(=O)=O", "OS(=O)(=O)CCCC"],
            vec!["CC(=C)CS(O)(=O)=O", "OS(=O)(=O)CC(=C)C"],
            vec!["CCCCCOS(O)(=O)=O", "OS(=O)(=O)OCCCCC"],
        ] {
            // 前提要自己成立:两种写法真的换了存储序
            let seq = |s: &str| -> Vec<(u32, u32)> {
                prep(s).bonds().iter().map(|b| (b.begin, b.end)).collect()
            };
            assert_ne!(seq(ws[0]), seq(ws[1]), "{} 与 {} 存储序一样", ws[0], ws[1]);
            let t0 = table(ws[0]);
            assert!(!t0.is_empty(), "{} 一根可翻的键都没有,验不了东西", ws[0]);
            assert_eq!(
                t0,
                table(ws[1]),
                "{} 与 {}:同一根键翻的不是同一侧",
                ws[0],
                ws[1]
            );
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
        let mut checked = 0usize;
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
                for p in layout::layout_all(&m, &ranks, style, None) {
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
                    let Some(side) = far_side(&m, &pos, b, &ranks) else {
                        continue;
                    };
                    checked += 1;
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
        // 一根可翻的键都没有时,下面那圈断言一次都不跑
        assert!(checked > 0, "一根可翻的键都没有,判据空过了");
    }

    #[test]
    fn a_ring_bond_is_never_flipped() {
        // 翻环上的键会把环撕开。`far_side` 在环感知标记之外再判一次连通性,
        // 这里守的是那一层。
        let m = prep("c1ccccc1");
        let ranks = omgkit_io::canon::canonical_ranks(&m);
        let mut pos: BTreeMap<u32, Point2> = BTreeMap::new();
        for p in layout::layout_all(&m, &ranks, &Style::ACS_1996, None) {
            pos.extend(p.pos);
        }
        for b in 0..u32::try_from(m.num_bonds()).unwrap() {
            assert!(
                far_side(&m, &pos, b, &ranks).is_none(),
                "环上的键 {b} 不该给出可翻的一侧"
            );
        }
    }

    /// **一个假环比任意多处交叉都糟。**
    ///
    /// [`score`] 的排序键是 `(画在同一点上的对数, 交叉数, 碰撞深度)`。这条判据
    /// 直接钉这个次序 —— 不靠语料,因为那一位在本语料上一次都没触发(插桩
    /// 实测 0 次),拿分子去验会是空过的。
    ///
    /// 理由见 `score` 的文档:两个原子叠在一点会让图上凭空多出一个环,而读者
    /// **没有任何办法看出那个环是假的**;交叉只是难读,读者知道那儿有两根键。
    ///
    /// 变异:`better` 里去掉第一位(退回 `(交叉, 深度)`)→ 这条当场红。
    #[test]
    fn a_phantom_ring_is_worse_than_any_number_of_crossings() {
        // 一处重合、零交叉、零深度  vs  零重合、五处交叉、深度 9.9
        assert!(
            !better((1, 0, 0.0), (0, 5, 9.9)),
            "拿一个假环换掉五处交叉被判成了改善"
        );
        assert!(
            better((0, 5, 9.9), (1, 0, 0.0)),
            "从假环换成五处交叉反倒没被判成改善"
        );
        // 重合一样多时,交叉才说话;交叉也一样时才比深度
        assert!(better((1, 0, 0.0), (1, 1, 0.0)), "重合相同时交叉少的该赢");
        assert!(
            better((1, 1, 0.0), (1, 1, 1.0)),
            "重合与交叉都相同时深度小的该赢"
        );
    }

    #[test]
    fn what_cannot_be_fixed_is_reported_not_hidden() {
        // 翻转解决不了的必须留在 `unresolved` 里。悄悄清空它,图上还是挤的,
        // 而调用方以为一切正常 —— 那比报出来糟得多。
        //
        // 六个硝酸根配位在一个铈上(真语料第 449 行)。金属只有六个方向,每个
        // 方向上挂的是 `O`—`N⁺`(`O⁻`)`=O` —— 一个键长的间距塞不下这些字,**这与
        // 布局挑得好不好无关**,是键长定死为 1 之后的算术。
        //
        // **这里原本用的是六叔丁基苯,已经换掉。** 那个分子在化学上确实极度
        // 拥挤,但**图上排得开** —— 六个叔丁基的甲基落在半径 3 那一圈,周长
        // 18.8 个单位、十八个甲基各要半个单位,绰绰有余。挑方向那一步一学会
        // 看标签(见 `chains::free_direction`),它就画干净了,这条判据的前提
        // 跟着失效。前提是化学直觉,而判据要的是几何事实,不是一回事。
        let (_, _, rep) = laid(
            "[O-][N+](=O)O[Ce](O[N+]([O-])=O)(O[N+]([O-])=O)(O[N+]([O-])=O)\
             (O[N+]([O-])=O)O[N+]([O-])=O",
            &Style::ACS_1996,
        );
        assert!(
            !rep.unresolved.is_empty(),
            "六个硝酸根挤在一个铈上不可能全排开,却报告说没有冲突"
        );
    }

    #[test]
    fn the_two_styles_can_disagree_about_whether_it_clashes() {
        // **这是 Style 参与布局的落点。** 同一张图,ACS 的标签占 0.69 个键长、
        // ChemDraw 默认占 0.33 —— 判定必须不同,否则把 Style 传进来就是白传。
        let m = prep("OC(=O)c1ccccc1OC(C)=O");
        let ranks = omgkit_io::canon::canonical_ranks(&m);
        let mut pos: BTreeMap<u32, Point2> = BTreeMap::new();
        for p in layout::layout_all(&m, &ranks, &Style::ACS_1996, None) {
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
