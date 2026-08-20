//! **界矩阵** —— 把化学事实翻译成"每一对原子之间必须离多远"。
//!
//! 这是全部化学进入算法的**唯一入口**。后面的光滑化、嵌入、精修都只认距离,
//! 不认元素、不认环、不认金属。所以:
//!
//! > **新来一类分子,要改的是这张表,不是代码分支。**
//!
//! 这条是可判定的 —— 嵌入器里出现 `is_macrocycle` / `is_metal` 这种谓词就算违规。
//!
//! # 四类约束
//!
//! | 拓扑距离 | 从哪来 |
//! |---|---|
//! | 1-2(成键) | 实测键长表的 `p05`/`p95` |
//! | 1-3(隔一个原子) | 实测键角表 + 余弦定理 |
//! | 1-4(隔两个原子) | 扭转角从顺式转到反式,取距离的两端 |
//! | ≥1-5 | 下限 = vdW 半径之和 ×按拓扑距离分档的 vdW 系数;上限交给三角光滑化去收 |
//!
//! **环不需要特殊处理**:环上两个原子之间有两条路径,光滑化取最短那条,
//! 上限自然被压到"绕不过去"的值。稠环、螺环、桥环、大环同理 —— 一条代码都不用加。
//!
//! # 与 RDKit 的两处实质差别
//!
//! 1. **配位数 > 4 的中心。** RDKit(`BoundsMatrixBuilder.cpp:541-552`)把 1-3 界写成
//!    `[1.0, 1.2×(b₁+b₂)]`,注释原话是 "an arbitrary min angle" —— 那个 `1.0` 没有
//!    任何物理依据,而上限比 180° 还宽。等于**六个配体之间毫无角度约束**,
//!    金属必然落在质心,然后被嵌入那一步的检查判死。实测 RDKit 剩下的失败**全是**
//!    金属配合物。这里改用按配位数给的角度包络([`coord_angle_envelope`]),
//!    同一段代码服务所有配位数。
//! 2. **查不到参数的原子。** RDKit 在 UFF 认不出元素时把键长界放成 `[0.5×d, 1.5×d]`
//!    (`BoundsMatrixBuilder.cpp:244-252`,实测 SF₆ 的 S–F 是 `[0.82, 2.47]`)。
//!    这里退到共价半径之和 ±10%,并且 [`Source`] 会**说出自己走的哪一级**。

use crate::params::{self, Source};
use crate::smooth::Bounds;
use omgkit_core::MolBuilder;

/// **非键下限的系数,按拓扑距离分档**(与 RDKit `BoundsMatrixBuilder.cpp:262-288` 同法)。
///
/// | 拓扑距离 | 系数 | 理由 |
/// |---|---|---|
/// | 4(1-5) | 0.70 | 被链拽在一起,挤得动 |
/// | 5(1-6) | 0.85 | 少挤一点 |
/// | ≥6 | 1.00 | 没理由比 vdW 接触还近 |
///
/// 头一版一律用 0.75,理由写的是"真实分子里非键接触本来就压得比 vdW 之和近"。
/// **那句话只对拓扑上近的成立** —— 远处的原子对没有链把它们拽在一起,
/// 而把它们的下限压到 0.75 就等于白白放宽 6.7 万对(语料里占总数六成)的区间。
#[must_use]
pub fn vdw_frac(topo_dist: usize) -> f64 {
    match topo_dist {
        0..=4 => 0.70,
        5 => 0.85,
        _ => 1.00,
    }
}

/// 没有任何约束的原子对,上限先给这个数(Å),交给三角光滑化去收。
///
/// 与 RDKit 的 `MAX_UPPER` 同量级。取多大都行 —— 光滑化会把它压到
/// "沿着键网络走过去的最短路",所以这个数只是个占位。
pub const MAX_UPPER: f64 = 1000.0;

/// 1-2 距离的绝对容差(Å):区间取 `中位 ± 这个数`。
///
/// **与 RDKit 的 `DIST12_DELTA`(`BoundsMatrixBuilder.cpp:27`)取同一个值**,
/// 它的 1-2 界宽中位 0.020 就是 `2 × 0.01`,逐位对得上。
///
/// # 为什么不用统计分位
///
/// 头一版用中位 ±1.2%(相对),得到 0.032 —— 比 RDKit 松 1.6 倍。那个宽度里
/// 装的是**查表键(元素+键级+环尺寸)分辨不了的桶内化学差异**,而界矩阵要的是
/// "这一根键该多长",不是"这一类键能有多不一样"。
///
/// RDKit 的自信来自 UFF 给的键长就是一个确定值。我们的中位同样是一个确定值,
/// 差别只在它偏多少 —— 而这个目标是**给力场的起点**,偏一点由精修拉回中点。
/// 偏太多会在判据一("真实构象必须落在界内")上现形,所以这个数是有闸看着的。
pub const DIST12_TOL: f64 = 0.01;

/// 1-3 距离的绝对容差(Å)。同 [`DIST12_TOL`],对应 RDKit 的 `DIST13_TOL`。
///
/// **1-3 由中位键长 + 中位键角按余弦定理算出一个值,再 ± 这个容差** ——
/// 不是把两端键长区间的端点乘进去。头一版那么做,键长的松会**复利**进 1-3。
pub const DIST13_TOL: f64 = 0.04;

/// 1-4 距离的绝对容差(Å)。对应 RDKit 的 `GEN_DIST_TOL`(`BoundsMatrixBuilder.cpp:33`),
/// 它把解掉析取的 1-4 钉成宽度 `2 × 0.06 = 0.12`,与实测中位逐位吻合。
pub const DIST14_TOL: f64 = 0.06;

/// 1-5 链式约束的绝对容差(Å)。对应 RDKit 的 `DIST15_TOL`(`BoundsMatrixBuilder.cpp:33`)。
pub const DIST15_TOL: f64 = 0.08;

/// 键长区间的相对容差 —— 只在**查不到表、只能用共价半径模型**时用。
///
/// 模型的不确定度本来就大,给绝对容差不诚实。
///
/// **不用 p05/p95。** 那个跨度装的是查表键(元素+键级+环尺寸)**分辨不了**的
/// 真实化学差异 —— 拿它当界,等于把"这一类键的全部变化"都允许给每一根键。
/// 实测我们的 1-2 界宽中位 0.081 Å,而 RDKit 是 0.020 Å,松了 4 倍。
///
/// 界要紧,精修才有东西可依;真不该这么紧的会在判据一("真实构象落在界内")
/// 上现形 —— 两条闸互相顶着,所以这个数是量出来的,不是拍的。
pub const BOND_REL: f64 = 0.012;

/// 查得到表时,键角区间取 `中位 ± 这个数`(度)。理由同 [`BOND_REL`]。
pub const ANGLE_TOL: f64 = 2.5;

/// **环内 1-4 的扭转角区间**(度),按共处环的尺寸。
///
/// 这是光滑化**推不出来**的信息:它只知道两条路径多长,不知道环是平的。
/// 实测我们的 1-4 界宽中位 0.758 Å 而 RDKit 是 0.120 Å,松了 6 倍,
/// 根子就在这里 —— 我给的是"顺式到反式"的全程,而环上的扭转是被锁死的。
///
/// | 共处环 | 扭转 | 依据 |
/// |---|---|---|
/// | 芳环 | 0° | 平面 |
/// | 3、4 元 | 0° | 平面(三元必然,四元近似) |
/// | 5 元 | 0–40° | 信封/半椅 |
/// | 6 元 | 0–60° | 平面(芳/共轭)到椅式 |
/// | 7 元及以上 | 0–90° | 柔性 |
/// | 不共处一环 | 0–180° | 自由旋转 |
/// **环的内扭转角范围**(度) —— 从语料实测,不是拍的。
///
/// 表在 `data/mmff.ringtorsion.tsv`,由 `harness/measure_ring_torsion.py` 生成:
/// 8526 个 ETKDGv3 + MMFF94 **收敛**的分子(与键长键角表同一口径),
/// 逐环沿圈取连续四元组记 `|扭转角|`,按 **环尺寸 + 是否芳环** 分桶,取 `p05`/`p95`。
///
/// # 拍出来的数错在哪
///
/// 头一版这几个数是我按化学常识写的,实测下来小环还行、**大环错得离谱**:
///
/// | 环 | 我拍的上界 | 实测 p95 |
/// |---|---|---|
/// | 4 元 | 25 | 33.3 |
/// | 5 元 | 40 | 45.4 |
/// | 6 元(饱和) | 65 | 61.2 |
/// | 8 元 | 90 | **116.9** |
/// | 9 元 | 90 | **138.0** |
/// | 16 元 | 90 | **178** |
///
/// 大环卡在 90° 就是**把真实几何排除在界外** —— 判据一迟早要红。
/// 芳环也不是精确的 0(五元 p95 = 3.0、六元 2.2)。
///
/// 表里没有的尺寸(10~15、17 以上)按"越大越柔"退到全程 —— 大环本来就柔,
/// 硬给一个内插值是在编数。
#[must_use]
pub fn ring_internal_torsion(size: usize, aromatic: bool, sp3: bool) -> (f64, f64) {
    /// 一行:`(中位, p05, p95)`,单位度。
    type Row = (f64, f64, f64);
    static T: std::sync::OnceLock<std::collections::HashMap<(usize, bool, bool), Row>> =
        std::sync::OnceLock::new();
    let t = T.get_or_init(|| {
        let mut m = std::collections::HashMap::new();
        for line in include_str!("../data/mmff.ringtorsion.tsv").lines() {
            if line.starts_with('#') {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 7 {
                continue;
            }
            let (Ok(sz), Ok(ar), Ok(sp), Ok(p05), Ok(p95)) = (
                f[0].parse::<usize>(),
                f[1].parse::<u8>(),
                f[2].parse::<u8>(),
                f[5].parse::<f64>(),
                f[6].parse::<f64>(),
            ) else {
                continue;
            };
            let Ok(med) = f[4].parse::<f64>() else {
                continue;
            };
            m.insert((sz, ar == 1, sp == 1), (med, p05, p95));
        }
        m
    });
    // **取中位那一个值,不取 p05/p95 区间** —— 与 1-2、1-3 同一条道理:
    // 分位跨度装的是"这一类环能有多不一样",而界矩阵要的是"这一个环该多扭"。
    // 偏了由精修拉回,偏太多在判据一上现形。
    //
    // 饱和六元环的中位是 20.9° 而不是椅式的 55°,因为这一档里混着共轭近平面的环 ——
    // 我们只生成**一个**构型,挑中位是诚实的选择,力场会把它带到该去的极小。
    //
    // 表里没有的尺寸(10~15、17 以上)退到全程 —— 大环本来就柔,内插是在编数。
    t.get(&(size, aromatic, sp3))
        .map_or((0.0, 180.0), |&(med, p05, p95)| {
            if aromatic {
                // 芳环真的是**一个值**:六元芳环 p05=0.0、p95=2.2,整桶挤在 0 附近。
                // 这一档钉成中位是诚实的。
                (med, med)
            } else {
                // **饱和环不能钉成一个值** —— 它有多个构象。六元全 sp³ 那一桶
                // p05=19.0、p95=64.4(椅 ~55°、扭船 ~30°),钉在中位 53.8° 上
                // 就把扭船排除在界外了:实测判据一从 0.135% 涨到 0.431%,
                // 而 1-4 界宽只从 7.68 收到 7.54 —— **代价大、收益几乎没有**。
                //
                // 改用实测的 [p05, p95]:比"顺式到反式全程"紧得多,又包得住真实构象。
                (p05, p95)
            }
        })
}

/// **中心键是有立体标记的双键时,一条 1-4 路径的扭转角是确定值。**
///
/// 这是 RDKit `_setChain14Bounds`(`BoundsMatrixBuilder.cpp:1020-1038`)那一支:
/// 它调 `_getAtomStereo(bnd2, aid1, aid4)` **逐对**问"这两个原子是顺还是反",
/// 然后 `dl = du`。
///
/// **注意别把这条路的功劳记大了。** 先前这里写着"它靠这类解析把 59.1% 的 1-4 钉住" ——
/// 那是**记错了**:实测语料里 9298 根键上只有 **51 根**带立体标记(0.55%),
/// 这条路撑不起 59% 的钉住率。RDKit 那 58.7% 的来源是 **sp²–sp² 键**
/// (芳环、共轭、酰胺 —— 由共轭定平面),见 `_setInRing14Bounds` 里的 `preferCis`。
/// 我们目前只钉了芳**环**,没钉一般的 sp²–sp² 键 —— 那才是 1-4 界宽还差 7.5 倍的原因。
///
/// # 顺反离开参照没有意义
///
/// 一根双键两端各有两个取代基,说"同侧"总得回答"谁和谁同侧" ——
/// 所以 [`BondData::stereo`](omgkit_core::BondData::stereo) 必须配
/// [`stereo_atoms`](omgkit_core::BondData::stereo_atoms)(两侧各一个参照)。
///
/// 任意一对 `(i, j)` 的顺反,由立体值经**两次参照翻转**得到:`i` 不是这一侧的
/// 参照就翻一次,`j` 不是就再翻一次,翻偶数次等于没翻。
///
/// 返回 `None` 表示这根键没有立体标记(或参照缺失),由调用方退回环/自由那两条。
fn stereo_path_torsion(mol: &MolBuilder, bidx: usize, i: u32, k: u32, j: u32) -> Option<f64> {
    use omgkit_core::BondStereo;
    let bd = mol.bonds().get(bidx)?;
    let same_side = match bd.stereo {
        BondStereo::Z | BondStereo::Cis => true,
        BondStereo::E | BondStereo::Trans => false,
        BondStereo::None => return None,
    };
    let (ra, rb) = (bd.stereo_atoms[0], bd.stereo_atoms[1]);
    if ra == omgkit_core::BondData::NO_STEREO_ATOM || rb == omgkit_core::BondData::NO_STEREO_ATOM {
        return None;
    }
    // 参照是按 (begin 侧, end 侧) 存的;`k` 是路径上 `i` 那一头的原子
    let (ref_i, ref_j) = if k == bd.begin { (ra, rb) } else { (rb, ra) };
    let flips = usize::from(i != ref_i) + usize::from(j != ref_j);
    let cis = if flips % 2 == 0 {
        same_side
    } else {
        !same_side
    };
    Some(if cis { 0.0 } else { 180.0 })
}

/// **中心键在环上时,一条 1-4 路径的扭转角范围。**
///
/// 这是从 RDKit 读来的一课(`BoundsMatrixBuilder.cpp:1005-1038`):
/// 它**不取凸包,而是用化学把析取解掉** —— 双键上问立体描述符"这一对是顺是反",
/// 然后 `dl = du`,宽度是 `2×GEN_DIST_TOL = 0.12 Å`。
/// (它总共钉住 58.7% 的 1-4,但**主要不是靠立体标记** —— 见
/// [`stereo_path_torsion`] 里那段更正。)
///
/// 环上不用立体描述符就能解:环有一个内扭转 `τ`([`ring_internal_torsion`]),
/// 而环外的取代基相对环内的路径是**反过来的**。于是
///
/// | `i`、`j` 的归属 | 扭转 |
/// |---|---|
/// | 都在环里 | `τ` |
/// | 都在环外(邻位取代基,都朝外) | `τ` |
/// | 一个在环外 | `180 − τ` |
///
/// 芳环是 `τ = 0` 的特例,正好给出 0 / 0 / 180 —— **不需要单独一条分支**。
///
/// 返回 `None` 表示中心键不在任何环上,由调用方给自由旋转的全程。
struct Rings<'a> {
    sets: &'a [Vec<u32>],
    aromatic: &'a [bool],
    sp3: &'a [bool],
}

fn ring_path_torsion(
    mol: &MolBuilder,
    rings: &Rings<'_>,
    i: u32,
    k: u32,
    l: u32,
    j: u32,
) -> Option<(f64, f64)> {
    let (ring_sets, ring_aromatic, ring_sp3) = (rings.sets, rings.aromatic, rings.sp3);
    // 中心键所在的**最小**环说了算:环越小越硬,给出的约束越紧
    let mut best: Option<(usize, usize)> = None; // (尺寸, 下标)
    for (r, set) in ring_sets.iter().enumerate() {
        let has = |a: u32| set.binary_search(&a).is_ok();
        if has(k) && has(l) && best.map_or(true, |(sz, _)| set.len() < sz) {
            best = Some((set.len(), r));
        }
    }
    let (size, r) = best?;
    // **只有近平面的环才把扭转钉成一个值。**
    //
    // # 为什么饱和环不能钉
    //
    // 表里存的是 `|τ|` 的**中位**,而它与同一张表的键角**几何上不相容**:
    // 六元饱和环的中位扭转是 20.9°,可是拿同表的键长 1.527、键角 111.6° 解闭环,
    // 六个内扭转必然是 `|τ| = 54.4°`;反过来要让扭转都等于 20.9°,键角得是 118.9°。
    // 两个中位取自不同的子总体(前者混进了共轭近平面环),拼在一起是**摆不出来的构型**。
    //
    // 后果不是"界松一点",是**约束集自相矛盾**,而三角光滑化查不出来
    // (它只保证两两自洽,不保证三维摆得出来)。实测端到端:
    //
    // | | 能量 | 1-2 越界 | 1-3 | 1-4 |
    // |---|---|---|---|---|
    // | 无环分子(五个) | **0.0000** | 0% | 0% | 0% |
    // | 苯 | **0.0000** | 0% | 0% | 0% |
    // | 环己烷 | 1.0073 | 33.3% | 83.3% | **94.1%** |
    // | 环戊烷 | 1.0941 | 66.7% | 83.3% | **100%** |
    //
    // 无环与芳环精修到**恰好零**,饱和环全线崩 —— 优化器没问题,是约束在打架。
    //
    // RDKit 在这里也只对 sp²/sp² 钉住(`_setTwoInSameRing14Bounds`,
    // `BoundsMatrixBuilder.cpp:708-790`),其余"here we will assume anything is
    // possible",给顺式到反式的全程。这里照同一条规矩:**芳环才钉**,
    // 饱和环退回全程,宁可界松也不要摆不出来。
    // **只有总体齐整的那两桶才钉:芳环、全 sp³ 环。**
    // 混合桶(非芳、非全 sp³)里既有共轭近平面环又有半椅,中位描述的是两者都不是
    // 的东西 —— 那正是先前把六元环钉成 20.9° 的来源。
    // **中心键两端都是 sp² 且环不大于 8 元 → 共轭把这一段定成平面,扭转钉 0。**
    //
    // 这是 RDKit 的规矩(`_setInRing14Bounds`:`ringSize <= 8 && ahyb2 == SP2 &&
    // ahyb3 == SP2` 就 `preferCis`),而它 58.7% 的 1-4 钉住率**主要就来自这一条** ——
    // 不是来自立体标记(语料 9298 根键只有 51 根带立体标记,0.55%),
    // 也不是来自饱和环(那一档 RDKit 同样不钉)。
    //
    // 先前这里只认**芳环**,于是非芳的共轭环(环己烯酮、马来酰亚胺、内酰胺……)
    // 全漏了 —— 那正是我们 1-4 界宽还差 7.5 倍的来源。
    let sp2 = |a: u32| mol.atoms()[a as usize].hybridization == omgkit_core::Hybridization::Sp2;
    let conj_planar = size <= 8 && sp2(k) && sp2(l);
    if !ring_aromatic[r] && !ring_sp3[r] && !conj_planar {
        return None;
    }
    let set = &ring_sets[r];
    let has = |a: u32| set.binary_search(&a).is_ok();
    let (t_lo, t_hi) = if conj_planar && !ring_aromatic[r] {
        // 共轭平面:钉 0(顺式)。与 RDKit 的 `compute14DistCis` 同一个意思。
        (0.0, 0.0)
    } else {
        ring_internal_torsion(size, ring_aromatic[r], ring_sp3[r])
    };
    Some(if has(i) == has(j) {
        (t_lo, t_hi)
    } else {
        (180.0 - t_hi, 180.0 - t_lo)
    })
}

/// 解不掉时按共处环尺寸给的扭转包络(度)。
#[must_use]
pub fn torsion_envelope(shared_ring: usize, aromatic: bool) -> (f64, f64) {
    // **四个原子不在同一个环里就没有环的约束。**
    //
    // 头一版只看"中心键是不是芳香的"就把扭转钉成 0,结果**芳环上的取代基**
    // (第一个原子在环外)被摁到了错的一侧 —— 它的扭转是 0 **或** 180,
    // 是个析取,塌成 0 就把上限压死了。实测真实距离 13.018 Å 越界 2.256 Å,
    // 判据一从 0.4% 炸到 15.3%。
    if shared_ring == 0 {
        return (0.0, 180.0);
    }
    if aromatic {
        return (0.0, 0.0);
    }
    match shared_ring {
        3 | 4 => (0.0, 0.0),
        5 => (0.0, 40.0),
        6 => (0.0, 60.0),
        _ => (0.0, 90.0),
    }
}

/// 一次界矩阵构建的分级计数。**每一条提前返回都要在这儿留痕。**
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Stats {
    /// 原子数。
    pub atoms: usize,
    /// 写下的 1-2 约束条数。
    pub n12: usize,
    /// 写下的 1-3 约束条数。
    pub n13: usize,
    /// 写下的 1-4 约束条数。
    pub n14: usize,
    /// **两条 1-3 估计交空、只好退回并集的次数。**
    ///
    /// 小环里同一对原子会同时是两个不同中心的 1-3(氧杂环丁烷 `C1COC1` 的
    /// C1···C3 既经 C0 也经 O2)。两条估计各带 ±[`DIST13_TOL`],
    /// 两张表行差过 `2×DIST13_TOL` 就交不上。
    ///
    /// **这个数量的是参数表的自相矛盾,不是分子的毛病。** 先前硬交集会把它
    /// 变成一个空区间,于是查表精度问题伪装成"分子不可行"计进覆盖率损失 ——
    /// 全语料 8831 个分子里有 21 个是这么死的。现在退并集保住可行性,
    /// 但**必须留下这个计数**,否则表越写越矛盾也没人看得见。
    pub n13_conflict: usize,
    /// **几何退化、这条 1-4 路径没写下来的次数。**
    ///
    /// `Stats` 的规矩是"每一条提前返回都要在这儿留痕"(见结构体文档),
    /// 而这一处先前是个**不留痕的 `continue`** —— 丢掉约束只会让界更松,
    /// 判据一反而更容易绿,于是覆盖率的损失对任何判据都不可见。
    /// 实测全语料曾因此丢掉 1697 条 1-4(0.31%),`CC#CC` 这类分子是**全丢**。
    ///
    /// **改用中点之后这个数是 0**,而且是结构性的 0:`d_il` 本来就是由 `d_ik`、
    /// `d_kl` 与夹角按余弦定理算出来的,三角形按构造自洽。留着计数与
    /// `torsion_span` 里那道 `.max(0.0)` 是防浮点噪声,不是防已知情形 ——
    /// 这个数一旦不再是 0,说明有别的东西变了。
    pub n14_degenerate: usize,
    /// **由 1-5 链式约束写下的条数。**
    ///
    /// 一条 5 原子路径上两个扭转都被钉住时,整段几何就定死了 ——
    /// 这是三角光滑化推不出来的(它只知道两两的界,不知道两个扭转是联动的)。
    pub n15: usize,
    /// 键长逐项命中表 / 退到"不在环里" / 只能用共价半径模型。
    pub bond_table: usize,
    /// 见 [`Stats::bond_table`]。
    pub bond_relaxed: usize,
    /// 见 [`Stats::bond_table`]。
    pub bond_model: usize,
    /// 键角逐项命中 / 放宽 / 兜底。
    pub angle_table: usize,
    /// 见 [`Stats::angle_table`]。
    pub angle_relaxed: usize,
    /// 见 [`Stats::angle_table`]。
    pub angle_model: usize,
    /// **走了配位数包络**([`coord_angle_envelope`])的 1-3 约束条数。
    ///
    /// 这一档是给配位数 ≥ 5 的中心用的 —— 它们的角不是一个值而是一**组**
    /// (八面体是 90° 与 180°),单个区间只能给出析取约束的松弛。
    /// 计数是为了知道有多少约束落在这条较松的路上。
    pub angle_envelope: usize,
}

/// **按配位数给出配体之间夹角的包络**(度)。
///
/// 配位数 ≥ 5 的中心,两个配体之间的夹角不是一个值而是一**组**:
/// 八面体是 90° 或 180°,三角双锥是 90°/120°/180°。单个区间表达不了"或",
/// 所以这里给的是**包络** `[最小可能, 最大可能]` —— 比真值松,
/// 但比 RDKit 那个 `[1.0 Å, 1.2×(b₁+b₂)]` 紧得多,而且有物理依据。
///
/// | 配位数 | 理想多面体 | 配体夹角 | 包络 |
/// |---|---|---|---|
/// | 2 | 直线 | 180 | 180–180 |
/// | 3 | 三角平面 | 120 | 120–120 |
/// | 4 | 四面体 | 109.47 | 109.47–109.47 |
/// | 5 | 三角双锥 | 90 / 120 / 180 | 90–180 |
/// | 6 | 八面体 | 90 / 180 | 90–180 |
/// | 7 | 五角双锥 | 72 / 90 / 144 / 180 | 72–180 |
/// | ≥8 | 四方反棱柱等 | — | 70–180 |
///
/// **这不是"金属分支"** —— 它只认配位数,不认元素。碳要是有六个邻居,
/// 走的也是同一行。
#[must_use]
pub fn coord_angle_envelope(degree: usize) -> (f64, f64) {
    match degree {
        0 | 1 => (180.0, 180.0),
        2 => (180.0, 180.0),
        3 => (120.0, 120.0),
        4 => (109.4712, 109.4712),
        5 | 6 => (90.0, 180.0),
        7 => (72.0, 180.0),
        _ => (70.0, 180.0),
    }
}

/// 由两条边长与夹角算第三边(余弦定理)。
fn third_side(a: f64, b: f64, theta_deg: f64) -> f64 {
    let t = theta_deg.to_radians();
    (a * a + b * b - 2.0 * a * b * t.cos()).max(0.0).sqrt()
}

/// 一根键的键长区间。
fn bond_range(mol: &MolBuilder, i: u32, j: u32, min_ring: usize, st: &mut Stats) -> (f64, f64) {
    let ord = mol
        .neighbors(i)
        .find(|(y, _)| *y == j)
        .map_or(omgkit_core::BondOrder::Single, |(_, bi)| {
            mol.bonds()[bi as usize].order
        });
    let p = params::bond_length(
        mol.atoms()[i as usize].atomic_num,
        mol.atoms()[j as usize].atomic_num,
        ord,
        min_ring,
    );
    match p.source {
        Source::Table => st.bond_table += 1,
        Source::RingRelaxed => st.bond_relaxed += 1,
        Source::Model => st.bond_model += 1,
    }
    // 查得到表就收紧到中位 ± 相对容差;只能用模型时保留它自己那个较宽的区间
    if p.source == Source::Model {
        (p.value * (1.0 - BOND_REL), p.value * (1.0 + BOND_REL))
    } else {
        (p.value - DIST12_TOL, p.value + DIST12_TOL)
    }
}

/// 一个中心的键角区间(度)。配位数 ≥ 5 走包络。
fn angle_range(
    mol: &MolBuilder,
    c: u32,
    ring_self: usize,
    ring_shared: usize,
    st: &mut Stats,
) -> (f64, f64) {
    let deg = mol.neighbors(c).count();
    // **配位数 ≥ 5 一律走包络。** 实测表是在 MMFF 优化过的结构上量的,
    // 那批结构里没有超配位中心,所以表对它们本来就没有发言权。
    if deg >= 5 {
        st.angle_envelope += 1;
        return coord_angle_envelope(deg);
    }
    let p = params::angle(
        mol.atoms()[c as usize].atomic_num,
        deg,
        mol.atoms()[c as usize]
            .flags
            .contains(omgkit_core::AtomFlags::AROMATIC),
        ring_self,
        ring_shared,
    );
    match p.source {
        Source::Table => st.angle_table += 1,
        Source::RingRelaxed => st.angle_relaxed += 1,
        Source::Model => st.angle_model += 1,
    }
    if p.source == Source::Model {
        (p.lo.to_degrees(), p.hi.to_degrees())
    } else {
        let v = p.value.to_degrees();
        ((v - ANGLE_TOL).max(1.0), (v + ANGLE_TOL).min(180.0))
    }
}

/// 一对原子拿到的是哪一档约束。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kind {
    /// 还没有任何拓扑约束 —— 最后铺 vdW 地板。
    None,
    /// 成键。
    B12,
    /// 隔一个原子。
    B13,
    /// 隔两个原子。
    B14,
}

/// 把 `[lo, hi]` **交**进已有区间(取更紧的那一侧)。
fn tighten(b: &mut Bounds, i: usize, j: usize, lo: f64, hi: f64) -> bool {
    let (clo, chi) = (b.lower(i, j), b.upper(i, j));
    let (nlo, nhi) = (clo.max(lo), chi.min(hi));
    if nlo <= nhi {
        b.set_lower(i, j, nlo);
        b.set_upper(i, j, nhi);
        return false;
    }
    // **交空了不等于这个分子摆不出来。**
    //
    // 两条约束各自带 ±DIST13_TOL = 0.04 Å,只要两张表行差过 0.08 Å 就交不上。
    // 那说的是**参数表自相矛盾**,不是几何不可能 —— 而硬交集会把它伪装成
    // "分子不可行",于是一个查表精度问题被记成了覆盖率损失。
    //
    // 交空时退回**并集**(RDKit 的 `_checkAndSetBounds` 一直是这么做的,
    // 它旁边的注释写着 "conservative bound setting")。这样比 RDKit 严:
    // 交得上时我们取交集拿到更紧的界,交不上才退到与它相同的并集。
    //
    // **退了必须记账** —— 见 `Stats::n13_conflict`。
    b.set_lower(i, j, clo.min(lo));
    b.set_upper(i, j, chi.max(hi));
    true
}

/// 每个原子所在的最小环尺寸(不在环里记 0),以及环的原子集合。
fn ring_info(mol: &MolBuilder) -> (Vec<usize>, Vec<Vec<u32>>) {
    let n = mol.num_atoms();
    let rings = omgkit_chem::sssr::ring_set(mol);
    let mut min_ring = vec![0usize; n];
    let mut sets: Vec<Vec<u32>> = Vec::with_capacity(rings.len());
    for r in &rings {
        let sz = r.atoms.len();
        for a in &r.atoms {
            let a = *a as usize;
            if a < n && (min_ring[a] == 0 || sz < min_ring[a]) {
                min_ring[a] = sz;
            }
        }
        let mut v = r.atoms.clone();
        v.sort_unstable();
        sets.push(v);
    }
    (min_ring, sets)
}

/// 含这几个原子的最小环尺寸;不共处一环记 0。
fn shared_ring(sets: &[Vec<u32>], atoms: &[u32]) -> usize {
    let mut best = 0usize;
    for s in sets {
        if atoms.iter().all(|a| s.binary_search(a).is_ok()) && (best == 0 || s.len() < best) {
            best = s.len();
        }
    }
    best
}

/// **建界矩阵。**
///
/// 返回的矩阵**尚未光滑** —— 光滑化是下一步([`crate::smooth::triangle_smooth`]),
/// 分开是为了让判据能拿 RDKit 的未光滑矩阵直接比。
///
/// # 环尺寸是查表键的一列,不是一条分支
///
/// 头一版这里传 `min_ring = 0`,理由写的是"环张力由光滑化处理"。**那是错的:
/// 光滑化只会收紧,它压不下一个本来就抬太高的下限。** 三元环的 C–C–C 角真值是
/// 60°,查表按无环给约 109°,推出 1-3 下限 2.45 Å —— 而那两个原子在三元环里
/// **本身就成键**(1.5 Å)。实测螺环 `C1CCC12CCC2` 当场出现"下限 2.79 > 上限 1.57"。
///
/// 所以环尺寸必须进查表键。零期量表时就分了这两维(中心自己的最小环、
/// 三个原子共处的最小环),正是为了这件事。这仍然是**表驱动** ——
/// 环尺寸只是键的一列,与元素、配位数同级,不是 `if (in_ring)`。
#[must_use]
pub fn build(mol: &MolBuilder) -> (Bounds, Stats) {
    let n = mol.num_atoms();
    let mut st = Stats {
        atoms: n,
        ..Stats::default()
    };
    let mut b = Bounds::new(n, 0.0, MAX_UPPER);
    if n == 0 {
        return (b, st);
    }
    let (min_ring, ring_sets) = ring_info(mol);
    // 每个环是不是芳环 —— 环上的原子全带芳香标志才算
    let aromatic_ring: Vec<bool> = ring_sets
        .iter()
        .map(|set| {
            set.iter().all(|a| {
                mol.atoms()[*a as usize]
                    .flags
                    .contains(omgkit_core::AtomFlags::AROMATIC)
            })
        })
        .collect();
    // 每个环是不是**全 sp³** —— 这一维把"共轭近平面"与"椅/船"两个总体分开。
    // 不分开的话中位落在两者之间,与同一张表的键角**几何上不相容**(见
    // `ring_path_torsion` 里那段账)。
    let sp3_ring: Vec<bool> = ring_sets
        .iter()
        .map(|set| {
            set.iter()
                .all(|a| mol.atoms()[*a as usize].hybridization == omgkit_core::Hybridization::Sp3)
        })
        .collect();
    let mut kind = vec![Kind::None; n * n];

    // ---- 1-2 ----
    for bd in mol.bonds() {
        let (i, j) = (bd.begin as usize, bd.end as usize);
        let r = shared_ring(&ring_sets, &[bd.begin, bd.end]);
        let (lo, hi) = bond_range(mol, bd.begin, bd.end, r, &mut st);
        tighten(&mut b, i, j, lo, hi);
        kind[i * n + j] = Kind::B12;
        kind[j * n + i] = Kind::B12;
        st.n12 += 1;
    }

    // ---- 1-3:中心 c 的两个邻居 ----
    for (c, &self_ring) in min_ring.iter().enumerate().take(n) {
        let Ok(cu) = u32::try_from(c) else { continue };
        let nb: Vec<u32> = mol.neighbors(cu).map(|(y, _)| y).collect();
        if nb.len() < 2 {
            continue;
        }
        for x in 0..nb.len() {
            for y in (x + 1)..nb.len() {
                let (i, j) = (nb[x] as usize, nb[y] as usize);
                if i == j {
                    continue;
                }
                let shared = shared_ring(&ring_sets, &[nb[x], cu, nb[y]]);
                let (ang_lo, ang_hi) = angle_range(mol, cu, self_ring, shared, &mut st);
                // **由中点算一个值,再 ± 绝对容差。**
                // 头一版拿两端键长区间的端点配最小/最大角,键长的松会**复利**进 1-3。
                let b1 = ((b.lower(i, c)) + (b.upper(i, c))) / 2.0;
                let b2 = ((b.lower(c, j)) + (b.upper(c, j))) / 2.0;
                // **角也取中点,不取区间端点。** `DIST13_TOL` 本来就是用来吸收
                // 角度不确定度的,再叠一个角度区间等于把同一件事算两遍 ——
                // 实测那样 1-3 界宽 0.141,而只用中点是 0.080(与 RDKit 相同)。
                //
                // 配位数 ≥5 的中心是例外:它的"角"是个**包络**(八面体 90 或 180),
                // 那是真的析取,不能塌成中点。
                let (lo, hi) = if ang_hi - ang_lo > 2.0 * ANGLE_TOL + 1e-9 {
                    (
                        third_side(b1, b2, ang_lo) - DIST13_TOL,
                        third_side(b1, b2, ang_hi) + DIST13_TOL,
                    )
                } else {
                    let d = third_side(b1, b2, ((ang_lo) + (ang_hi)) / 2.0);
                    (d - DIST13_TOL, d + DIST13_TOL)
                };
                // **这一对已经是键了就别写 1-3。** 三元/四元环里同一对可以既是
                // 1-2 又是 1-3(环硫乙烷 `C1CS1` 的两个碳既成键、又同时连着硫)。
                // 直接量到的键长比"两根键 + 夹角推出来的距离"可靠得多,
                // 让角推的那个去挤键,是拿差的信息覆盖好的信息。
                //
                // 头一版是硬交集,于是 `C1CS1` 的 C–C 被推成 [1.4980 > 1.3829] ——
                // **连成键那一对都被 1-3 干掉了**,分子当场不可行。
                if kind[i * n + j] == Kind::B12 {
                    continue;
                }
                // 取交集,不是覆盖;交空了退并集并记账(见 `tighten`)
                if tighten(&mut b, i, j, lo, hi) {
                    st.n13_conflict += 1;
                }
                if kind[i * n + j] == Kind::None {
                    kind[i * n + j] = Kind::B13;
                    kind[j * n + i] = Kind::B13;
                }
                st.n13 += 1;
            }
        }
    }

    // ---- 1-4:路径 i–k–l–j,扭转角从顺式(0°)到反式(180°) ----
    for (bidx, bd) in mol.bonds().iter().enumerate() {
        let (k, l) = (bd.begin, bd.end);
        let nk: Vec<u32> = mol
            .neighbors(k)
            .map(|(y, _)| y)
            .filter(|y| *y != l)
            .collect();
        let nl: Vec<u32> = mol
            .neighbors(l)
            .map(|(y, _)| y)
            .filter(|y| *y != k)
            .collect();
        for &i in &nk {
            for &j in &nl {
                if i == j {
                    continue; // 三元环:这一对已经是 1-3
                }
                let (iu, ju) = (i as usize, j as usize);
                let (ku, lu) = (k as usize, l as usize);
                // 已经是 1-2 或 1-3 的对(四元环等),它们的约束更紧,别插手
                if kind[iu * n + ju] != Kind::None {
                    continue;
                }
                // **五个距离全取中点,不取上限。**
                //
                // 这与上面 1-3 那段是同一条规矩(那里写着"由中点算一个值,
                // 再 ± 绝对容差"),而这里原先全用 `b.upper` —— 同一个文件里
                // 两套做法。后果有两个:
                //
                // 1. **整条 1-4 系统性外移。** 苯的 para C···C 用上限算出
                //    [2.8372, 2.9572],而同一张表摆出的正六边形是 2×1.397 = 2.794 ——
                //    **界把自己表里的几何排除在外**。用中点得 [2.7382, 2.8582],正好包住。
                // 2. **经过近线性中心时整条 1-4 被静默丢弃。** sp 碳的角是 178.9°,
                //    于是 `d_il` 的上限 ≈ 两键中点之和 + 0.04,而
                //    `d_ik上限 + d_kl上限` = 两键中点之和 + 0.02 —— 三角形不成立,
                //    `torsion_span` 返回 `None`。实测 `CC#CC`、`CC#N`、`C(#N)C#N`
                //    的 1-4 **一条都没写下来**,全语料丢 1697 条(0.31%)。
                //
                // 中点没有这两个毛病:`d_il` 本来就是由 `d_ik`、`d_kl` 与夹角
                // 按余弦定理算出来的,所以这个三角形**按构造自洽**。
                let mid = |x: usize, y: usize| (b.lower(x, y) + b.upper(x, y)) / 2.0;
                let Some((cis, trans)) = torsion_span(
                    mid(iu, ku),
                    mid(ku, lu),
                    mid(lu, ju),
                    mid(iu, lu),
                    mid(ku, ju),
                ) else {
                    st.n14_degenerate += 1;
                    continue;
                };
                // **环上的扭转是锁死的,不是顺式到反式的全程。**
                // 光滑化推不出这件事(它只知道路径长度,不知道环是平的)。
                let arom = mol.bonds()[bidx].order == omgkit_core::BondOrder::Aromatic;
                // 芳环上的扭转是**确定值**(见 `planar_ring_torsion`),不是区间;
                // 解不掉的才退回按环尺寸的包络
                // 顺序:立体标记 > 环 > 自由旋转。**有立体标记的双键最硬**,
                // 它给的是确定值;环给区间;都没有才退回顺式到反式的全程。
                let (t_lo, t_hi) = stereo_path_torsion(mol, bidx, i, k, j)
                    .map(|t| (t, t))
                    .or_else(|| {
                        ring_path_torsion(
                            mol,
                            &Rings {
                                sets: &ring_sets,
                                aromatic: &aromatic_ring,
                                sp3: &sp3_ring,
                            },
                            i,
                            k,
                            l,
                            j,
                        )
                    })
                    .unwrap_or_else(|| {
                        torsion_envelope(shared_ring(&ring_sets, &[i, k, l, j]), arom)
                    });
                // d² 随 cos(扭转) 单调 —— 端点取遍区间即可
                let f = |t: f64| {
                    let c = t.to_radians().cos();
                    let (a2, b2) = (cis * cis, trans * trans);
                    (a2 + (b2 - a2) * (1.0 - c) / 2.0).max(0.0).sqrt()
                };
                // 扭转被定死时,区间宽度就来自这个容差(与 RDKit 的 GEN_DIST_TOL 同值)
                tighten(&mut b, iu, ju, f(t_lo) - DIST14_TOL, f(t_hi) + DIST14_TOL);
                kind[iu * n + ju] = Kind::B14;
                kind[ju * n + iu] = Kind::B14;
                st.n14 += 1;
            }
        }
    }

    // 拓扑距离(封顶 6:分档只用到 4/5/≥6)
    let mut topo = vec![6u8; n * n];
    for start in 0..n {
        let mut d = vec![u8::MAX; n];
        d[start] = 0;
        let mut q = std::collections::VecDeque::from([start]);
        while let Some(x) = q.pop_front() {
            if d[x] >= 6 {
                continue;
            }
            let Ok(xu) = u32::try_from(x) else { continue };
            for (y, _) in mol.neighbors(xu) {
                let y = y as usize;
                if y < n && d[y] == u8::MAX {
                    d[y] = d[x] + 1;
                    q.push_back(y);
                }
            }
        }
        for j in 0..n {
            topo[start * n + j] = d[j].min(6);
        }
    }

    // ---- 1-5 链式约束:**试过,撤了** ----
    //
    // 想法(与 RDKit `BoundsMatrixBuilder.cpp:1997-2045` 同):一条 5 原子路径上
    // 两个扭转都被钉住时,整段几何就定死了,`a···e` 可以直接摆出来量,再 ± 0.08。
    // 实测确实把总体宽度比从 1.020 收到 1.005。
    //
    // **但它在环上是错的,所以撤掉了。** 椅式环己烷的内扭转是 +55、−55、+55
    // 交替的,而我们的实测表存的是 `|τ|` 的中位 —— **符号被扔了**。
    // 链式把两个 `+55` 接起来,摆出来的不是椅式,算出的 1-5 距离是错的,
    // 于是与别处的约束打架:实测界不可行的分子从 3 个涨到 **11 个**,
    // 也就是 8 个分子**完全嵌不出来**。
    //
    // 覆盖率是这个项目的头号指标(整件事就是为了赢 RDKit 那 0.52% 的失败率),
    // 拿 2% 的分子换 1.5% 的界宽是笔坏买卖。要做对得先有**带符号**的环构象
    // 分析(椅/船/信封各自的扭转序列),那是另一块活,不是调容差能补的。
    //
    // 判官那边留下了两样东西:界不可行的分子数现在**有闸**(先前只报不闸,
    // 于是界可以越写越不自洽而只表现为参与统计的对数在掉);而那次对数
    // 从 108981 崩到 9317 时,剩下一成样本上的比值反而更好看 ——
    // **样本被腰斩,数就不作数**。

    // ---- vdW 地板:**只铺给没有拓扑约束的对** ----
    // 1-2/1-3/1-4 是被键角扭转钉死的,它们本来就可以比 vdW 之和近
    // (sp 中心的 1-3 只有 2.3 Å,而 C···C 的 0.75×vdW 是 2.55)——
    // 给它们铺地板会把界矩阵直接铺成矛盾的。
    for i in 0..n {
        for j in (i + 1)..n {
            if kind[i * n + j] != Kind::None {
                continue;
            }
            let d0 = vdw_frac(topo[i * n + j] as usize)
                * (params::vdw_radius(mol.atoms()[i].atomic_num)
                    + params::vdw_radius(mol.atoms()[j].atomic_num));
            if d0 > b.lower(i, j) {
                b.set_lower(i, j, d0);
            }
        }
    }

    (b, st)
}

/// 由 1-2 与 1-3 距离算 1-4 距离在**顺式**与**反式**两端的值。
///
/// 把 `i–k–l–j` 摆进以 `k–l` 为轴的柱坐标:`i` 与 `j` 各自到轴的垂距 `ρ`
/// 与沿轴坐标 `z` 由三角形定死,扭转角只改变它们的方位差。于是
///
/// ```text
/// d² = Δz² + ρᵢ² + ρⱼ² − 2ρᵢρⱼ·cos(扭转角)
/// ```
///
/// 顺式(0°)取 `Δz² + (ρᵢ−ρⱼ)²`,反式(180°)取 `Δz² + (ρᵢ+ρⱼ)²`。
/// 这与消撞里那个"绕轴转动的距离下界"是同一个式子。
fn torsion_span(d_ik: f64, d_kl: f64, d_lj: f64, d_il: f64, d_kj: f64) -> Option<(f64, f64)> {
    if d_kl <= 1e-9 {
        return None;
    }
    // i 在以 k 为原点、k→l 为 +z 的柱坐标里
    let zi = (d_ik * d_ik + d_kl * d_kl - d_il * d_il) / (2.0 * d_kl);
    let ri2 = d_ik * d_ik - zi * zi;
    // j 同理,但以 l 为原点量,再换算到 k 的坐标
    let zj = (d_lj * d_lj + d_kl * d_kl - d_kj * d_kj) / (2.0 * d_kl);
    let rj2 = d_lj * d_lj - zj * zj;
    // **几何退化不等于没信息。** `ri2 < 0` 只可能来自浮点噪声或近共线,
    // 而共线的正确答案是半径 0(原子落在 k→l 轴上),不是"这条路径作废"。
    // 原先这里 `return None`,于是 sp 中心上整条 1-4 被静默丢掉 ——
    // 而丢掉只会让界更松,判据一反而更容易绿,没有任何闸看得见。
    // 现在用 `.max(0.0)` 夹住;真正非有限的情况仍由下面那道 `is_finite` 拦。
    let (ri, rj) = (ri2.max(0.0).sqrt(), rj2.max(0.0).sqrt());
    let dz = d_kl - zi - zj;
    let cis = (dz * dz + (ri - rj) * (ri - rj)).max(0.0).sqrt();
    let trans = (dz * dz + (ri + rj) * (ri + rj)).max(0.0).sqrt();
    (cis.is_finite() && trans.is_finite()).then_some((cis, trans))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn prep(smi: &str) -> MolBuilder {
        let mut m = omgkit_io::smiles::parse(smi).expect("SMILES 该解析得了");
        omgkit_chem::pipeline::sanitize(&mut m).expect("该 sanitize 得了");
        // **立体感知要在建界之前跑。** 界矩阵消费 `bond.stereo`,
        // 没感知的话双键那一支永远走不到,而判据看不出区别 —— 只是界更松。
        omgkit_io::stereo::perceive_bond_stereo(&mut m);
        let r = omgkit_io::canon::classed_ranks(&m);
        omgkit_chem::add_explicit_hs(&mut m, &r);
        m
    }

    /// **每一对原子的下限都不许超过上限。** 界矩阵自己就矛盾的话,后面全白搭。
    #[test]
    fn every_pair_has_a_non_empty_interval() {
        let mut checked = 0;
        let mut n_mol = 0;
        let smis = [
            "CCO",                // 链
            "c1ccccc1",           // 芳环
            "C1CCCCC1",           // 饱和环
            "C1CC2CCC1CC2",       // 桥环
            "C1CCC12CCC2",        // 螺环
            "c1ccc2ccccc2c1",     // 稠环
            "CC(=O)N(C)C",        // 酰胺
            "C#N",                // 三键(共线)
            "CC(Cl)=C=C(C)Cl",    // 累积双键
            "FS(F)(F)(F)(F)F",    // 超配位
            "N[Co](N)(N)(N)(N)N", // 金属配合物
            "C1CCCCCCCCCCC1",     // 大环
            // ---- 小环与小杂环:先前这张表里**一个三元环、一个小杂环都没有** ----
            //
            // 恰好表里那几个(螺环、桥环)全是**全碳**且查表都命中,所以全绿,
            // 而全语料上有 21 个分子建完界就是空区间。下面这几个是实测会死的:
            "C1CC1",         // 环丙烷
            "C1CS1",         // 环硫乙烷 —— 两个碳既成键、又同为 S 的 1-3
            "C1CO1",         // 环氧乙烷
            "C1NN1",         // 三元双氮
            "C1COC1",        // 氧杂环丁烷 —— C1···C3 同时是 C0 与 O2 的 1-3
            "C1C(CS1)O",     // 3-羟基硫杂环丁烷(取自 large.smi)
            "C1CC2C1C2",     // 双环[1.1.0]丁烷 —— 角的 ring_self=3 而 ring_shared=4
            "C1CCC2(C1)CC2", // 螺[2.4] —— 三元环与五元环共用一个碳
        ];
        for smi in smis {
            let m = prep(smi);
            let (b, _) = build(&m);
            n_mol += 1;
            for i in 0..b.len() {
                for j in (i + 1)..b.len() {
                    assert!(
                        b.lower(i, j) <= b.upper(i, j),
                        "{smi}:第 {i}/{j} 对的下限 {} 超过了上限 {}",
                        b.lower(i, j),
                        b.upper(i, j)
                    );
                    checked += 1;
                }
            }
        }
        // 防空转:分子一个都不许漏,对数也得像样。
        // **按分子数断言,不按对数** —— 对数一改分子列表就得跟着调,那是会陈的判据。
        assert_eq!(n_mol, smis.len(), "有分子没验到");
        assert!(checked > 300, "只验了 {checked} 对");
    }

    /// **1-2 的区间就是表里那一行。** 乙烷的 C–C 该落在实测的 p05~p95 之间。
    #[test]
    fn a_bond_gets_the_measured_interval() {
        let m = prep("CC");
        let (b, st) = build(&m);
        assert!(st.n12 >= 7, "乙烷该有 7 根键,写了 {}", st.n12);
        let bd = &m.bonds()[0];
        let want = params::bond_length(
            m.atoms()[bd.begin as usize].atomic_num,
            m.atoms()[bd.end as usize].atomic_num,
            bd.order,
            0,
        );
        let (i, j) = (bd.begin as usize, bd.end as usize);
        // **断契约:区间必须包住表里的中位,而且要紧。**
        // 头一版断的是"区间恰好等于 p05/p95" —— 那是把当时的实现抄了一遍,
        // 一收紧就红,而收紧正是要做的事。
        assert!(
            b.lower(i, j) <= want.value && want.value <= b.upper(i, j),
            "区间 [{:.4}, {:.4}] 没包住表里的中位 {:.4}",
            b.lower(i, j),
            b.upper(i, j),
            want.value
        );
        let w = b.upper(i, j) - b.lower(i, j);
        assert!(w < 0.06, "键长区间宽 {w:.4} —— 太松,精修就没东西可依");
    }

    /// **1-3 要由余弦定理算出来。** 水的 H···H:两条 O–H 各约 0.97,夹角约 104.5°。
    #[test]
    fn a_one_three_distance_comes_from_the_law_of_cosines() {
        let m = prep("O");
        let (b, _) = build(&m);
        let hs: Vec<usize> = (0..m.num_atoms())
            .filter(|i| m.atoms()[*i].atomic_num == 1)
            .collect();
        assert_eq!(hs.len(), 2, "水该有两个氢");
        let (lo, hi) = (b.lower(hs[0], hs[1]), b.upper(hs[0], hs[1]));
        // **参照取表里的中位,不另找来源。**
        //
        // 头一版写死了实验水的 0.97 Å / 104.5°(给 1.534),而区间是 [1.545, 1.633] ——
        // 差 0.011 Å 被判红。查下来不是代码错:参数表是从 **MMFF 优化过的结构**
        // 量出来的,**MMFF 的水不是实验的水**。拿另一个来源的数当参照,
        // 量的就不是"余弦定理有没有用对",而是"两个来源合不合"。
        //
        // 这里断的是契约:区间必须以"表里的中位键长 + 中位键角按余弦定理算出的值"为心。
        let bl = params::bond_length(8, 1, omgkit_core::BondOrder::Single, 0).value;
        let ang = params::angle(8, 2, false, 0, 0).value.to_degrees();
        let d = third_side(bl, bl, ang);
        assert!(
            lo <= d && d <= hi,
            "H···H 区间 [{lo:.3}, {hi:.3}] 没包住表算的 {d:.3}(键长 {bl:.3}、键角 {ang:.1}°)"
        );
        assert!(hi - lo < 0.2, "区间 [{lo:.3}, {hi:.3}] 太宽了");
    }

    /// **配位数 ≥5 走包络,而且比 RDKit 那个 `[1.0, 1.2×(b₁+b₂)]` 紧。**
    ///
    /// SF₆:六个氟两两之间,八面体的实际值是 90° 或 180°,包络给 [90, 180]。
    /// 对应距离约 [2.2, 3.1] Å —— 而 RDKit 给的是 [1.0, 3.7]。
    #[test]
    fn a_hypervalent_centre_gets_a_polyhedral_envelope_not_an_arbitrary_floor() {
        let m = prep("FS(F)(F)(F)(F)F");
        let (b, st) = build(&m);
        assert!(st.angle_envelope > 0, "六配位的硫该走包络");
        let fs: Vec<usize> = (0..m.num_atoms())
            .filter(|i| m.atoms()[*i].atomic_num == 9)
            .collect();
        assert_eq!(fs.len(), 6);
        let (lo, hi) = (b.lower(fs[0], fs[1]), b.upper(fs[0], fs[1]));
        assert!(
            lo > 1.5,
            "F···F 下限 {lo:.3} —— RDKit 那个无依据的 1.0 不该出现在这里"
        );
        assert!(hi < 3.6, "F···F 上限 {hi:.3} 该被 180° 压住");
        assert!(lo < hi);
    }

    /// **环不需要特殊处理:光滑化会把上限压到"绕不过去"的值。**
    ///
    /// 苯的对位碳,沿环走是三根键(约 4.2 Å),但真实距离约 2.8 Å ——
    /// 建界时给的是路径上界,光滑化不会把它压到 2.8(那要环闭合的信息),
    /// 但**两条路径**取短的这件事必须发生:苯的 1-4 有两条等长路径,
    /// 所以上限该是"三根键"而不是 `MAX_UPPER`。
    #[test]
    fn a_ring_needs_no_special_case_the_smoothing_finds_the_short_way_round() {
        let m = prep("c1ccccc1");
        let (mut b, _) = build(&m);
        crate::smooth::triangle_smooth(&mut b).expect("苯的界该可行");
        let cs: Vec<usize> = (0..m.num_atoms())
            .filter(|i| m.atoms()[*i].atomic_num == 6)
            .collect();
        assert_eq!(cs.len(), 6);
        // 找一对对位碳:环上距离最远的那对
        let mut worst = 0.0f64;
        for x in 0..6 {
            for y in (x + 1)..6 {
                worst = worst.max(b.upper(cs[x], cs[y]));
            }
        }
        assert!(
            worst < 6.0,
            "苯里最远的一对碳上限 {worst:.3} —— 光滑化没把 MAX_UPPER 压下来"
        );
    }

    /// **芳环上的 1-4 扭转是确定值,所以那些原子对的区间宽度该近乎为零。**
    ///
    /// 苯:环内的 `C–C–C–C` 扭转 0、环上氢的 `H–C–C–C` 扭转 180 ——
    /// 两者都是确定的,不是"顺式到反式"的区间。这条盯住 `planar_ring_torsion`
    /// 真的接上了:它不生效的话宽度会是 0.7 Å 上下。
    #[test]
    fn 三元环里成键那一对不许被_1_3_挤宽() {
        // 三元环上任意两个重原子**既成键、又同为第三个原子的 1-3**。
        // 直接量到的键长比"两根键 + 夹角推出来的距离"可靠得多,所以 1-3
        // 不许往已经成键的那一对上写 —— 让角推的去挤键,是拿差的信息盖好的。
        //
        // **这条规则先前没有任何判据看着。** 语料中位一位都不动(带杂原子的
        // 三元环太少),而实测撤掉它之后环硫乙烷的三根键宽从 0.020 涨到
        // 0.243 / 0.132 / 0.163 —— 最狠的一根**宽了 12 倍**。
        // 这已经是这一轮里第三次"中位藏住子群"了,所以这里直接钉具体分子。
        for smi in ["C1CS1", "C1CO1", "C1CC1", "C1NN1"] {
            let m = prep(smi);
            let (b, _) = build(&m);
            let z: Vec<u8> = m.atoms().iter().map(|a| a.atomic_num).collect();
            for bd in m.bonds() {
                let (i, j) = (bd.begin as usize, bd.end as usize);
                if z[i] == 1 || z[j] == 1 {
                    continue;
                }
                let w = b.upper(i, j) - b.lower(i, j);
                // 键宽应当就是 2×DIST12_TOL。写成常数的式子,不写死 0.02 ——
                // 容差改了这条判据要跟着改,而不是变成一句过期的断言。
                assert!(
                    w <= 2.0 * DIST12_TOL + 1e-9,
                    "{smi} 的键 {i}-{j} 宽 {w:.4},超过 2×DIST12_TOL —— 1-3 挤到键上了"
                );
                assert!(
                    b.lower(i, j) <= b.upper(i, j),
                    "{smi} 的键 {i}-{j} 区间空了"
                );
            }
        }
    }

    #[test]
    fn a_one_four_across_an_aromatic_ring_has_a_pinned_torsion() {
        let m = prep("c1ccccc1");
        let (b, _) = build(&m);
        let n = m.num_atoms();
        // 拓扑距离恰好 3 的原子对
        let mut widths = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                if topo_dist(&m, i, j) != 3 {
                    continue;
                }
                widths.push(b.upper(i, j) - b.lower(i, j));
            }
        }
        assert!(
            widths.len() >= 6,
            "苯该有不少 1-4 对,只找到 {}",
            widths.len()
        );
        let worst = widths.iter().fold(0.0f64, |a, x| a.max(*x));
        // **用常数表达,不写死数字。** 头一版这里写 `< 0.05`,而"钉住"的宽度
        // 后来变成 `2 × DIST14_TOL = 0.12` —— 断言当场红,可代码是**对的**
        // (那正是与 RDKit 相同的钉死宽度)。判据里的阈值一旦与被判的量脱钩,
        // 就会在下一次改进时冒充失败。
        assert!(
            worst <= 2.0 * DIST14_TOL + 1e-9,
            "苯的 1-4 最宽 {worst:.4} Å,超过钉死宽度 {:.4} —— 芳环扭转没被钉住",
            2.0 * DIST14_TOL
        );
    }

    /// 两点之间的拓扑距离(封顶 4)。
    fn topo_dist(m: &MolBuilder, a: usize, b: usize) -> u8 {
        let n = m.num_atoms();
        let mut d = vec![u8::MAX; n];
        d[a] = 0;
        let mut q = std::collections::VecDeque::from([a]);
        while let Some(x) = q.pop_front() {
            if d[x] >= 4 {
                continue;
            }
            for (y, _) in m.neighbors(u32::try_from(x).expect("下标")) {
                let y = y as usize;
                if d[y] == u8::MAX {
                    d[y] = d[x] + 1;
                    q.push_back(y);
                }
            }
        }
        d[b].min(4)
    }

    /// **顺式与反式必须给出不同的界。**
    ///
    /// 这条是防"立体那一支接了等于没接":如果 `stereo_path_torsion` 从不生效,
    /// 两个写法会拿到**同一个**顺式到反式的区间,判据照样绿 —— 只是界更松。
    ///
    /// 2-丁烯:顺式的两个甲基碳约 3.0 Å,反式约 3.7 Å。
    #[test]
    fn cis_and_trans_get_different_bounds() {
        let cis = prep(r"C/C=C\C");
        let trans = prep(r"C/C=C/C");
        let mut got = Vec::new();
        for m in [&cis, &trans] {
            let (b, _) = build(m);
            // 两个甲基碳:与双键碳成键、且自己只连一个重原子
            let cs: Vec<usize> = (0..m.num_atoms())
                .filter(|i| {
                    m.atoms()[*i].atomic_num == 6
                        && m.neighbors(u32::try_from(*i).expect("下标"))
                            .filter(|(y, _)| m.atoms()[*y as usize].atomic_num > 1)
                            .count()
                            == 1
                })
                .collect();
            assert_eq!(cs.len(), 2, "2-丁烯该有两个端甲基");
            let (lo, hi) = (b.lower(cs[0], cs[1]), b.upper(cs[0], cs[1]));
            assert!(hi - lo < 0.30, "区间 [{lo:.3}, {hi:.3}] 没被立体标记钉住");
            got.push(((lo) + (hi)) / 2.0);
        }
        assert!(
            got[1] - got[0] > 0.5,
            "顺式 {:.3}、反式 {:.3} —— 立体那一支没生效,两个写法拿到了同一个区间",
            got[0],
            got[1]
        );
    }

    /// 空分子、单原子:不许 panic。
    #[test]
    fn tiny_molecules_are_fine() {
        for smi in ["C", "[H][H]", "O"] {
            let m = prep(smi);
            let (b, st) = build(&m);
            assert_eq!(b.len(), m.num_atoms());
            assert_eq!(st.atoms, m.num_atoms());
        }
    }
}
