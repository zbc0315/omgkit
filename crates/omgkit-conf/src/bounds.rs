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
//! | ≥1-5 | 下限 = vdW 半径之和 ×[`VDW_FRAC`];上限交给三角光滑化去收 |
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

/// 非键原子对的下限 = 两个 vdW 半径之和乘这个系数。
///
/// 不取 1.0:真实分子里非键接触**本来就压得比 vdW 之和近**(氢键、堆积),
/// 取 1.0 会让界矩阵与现实矛盾。0.75 是本仓一贯用的值。
pub const VDW_FRAC: f64 = 0.75;

/// 没有任何约束的原子对,上限先给这个数(Å),交给三角光滑化去收。
///
/// 与 RDKit 的 `MAX_UPPER` 同量级。取多大都行 —— 光滑化会把它压到
/// "沿着键网络走过去的最短路",所以这个数只是个占位。
pub const MAX_UPPER: f64 = 1000.0;

/// 查得到表时,键长区间取 `中位 × (1 ± 这个数)`。
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
/// **平面环上,一条 1-4 路径的扭转角是确定值,不是区间。**
///
/// 这是从 RDKit 读来的关键一课(`BoundsMatrixBuilder.cpp:1005-1038`):
/// 它**不取凸包,而是用化学把析取解掉** —— 双键上问立体描述符"这一对是顺是反",
/// 然后 `dl = du`,宽度为 0。
///
/// 平面环上这件事更简单,连立体描述符都不用:环是平的,于是每个原子有个确定的
/// **侧** —— 环内原子在圆心那一侧,环外取代基朝外。中心键 `k–l` 在环上时:
///
/// | `i`、`j` 的归属 | 扭转 | 例子(苯) |
/// |---|---|---|
/// | 都在环里 | **0°** | `C6–C1–C2–C3` |
/// | 一个在环外 | **180°** | 取代基 `X–C1–C2–C3` |
/// | 都在环外 | **0°** | 邻位两个取代基 `X–C1–C2–Y`,都朝外,同侧 |
///
/// 返回 `None` 表示这条路径的中心键不在**平面**环上,由调用方退回包络。
///
/// # 只认芳环
///
/// 非芳香环的平面性没有保证(环己烷是椅式),所以这里只对芳环下确定值,
/// 其余交给 [`torsion_envelope`] 按环尺寸给区间。
fn planar_ring_torsion(
    ring_sets: &[Vec<u32>],
    aromatic_ring: &[bool],
    i: u32,
    k: u32,
    l: u32,
    j: u32,
) -> Option<f64> {
    for (r, set) in ring_sets.iter().enumerate() {
        if !aromatic_ring[r] {
            continue;
        }
        let has = |a: u32| set.binary_search(&a).is_ok();
        if !(has(k) && has(l)) {
            continue;
        }
        return Some(if has(i) == has(j) { 0.0 } else { 180.0 });
    }
    None
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
        (p.lo, p.hi)
    } else {
        (p.value * (1.0 - BOND_REL), p.value * (1.0 + BOND_REL))
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
fn tighten(b: &mut Bounds, i: usize, j: usize, lo: f64, hi: f64) {
    if lo > b.lower(i, j) {
        b.set_lower(i, j, lo);
    }
    if hi < b.upper(i, j) {
        b.set_upper(i, j, hi);
    }
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
                let (b1lo, b1hi) = (b.lower(i, c), b.upper(i, c));
                let (b2lo, b2hi) = (b.lower(c, j), b.upper(c, j));
                let lo = third_side(b1lo, b2lo, ang_lo);
                let hi = third_side(b1hi, b2hi, ang_hi);
                // **取交集,不是覆盖。** 三元/四元环里同一对可以既是 1-2 又是 1-3,
                // 头一版直接覆盖,把成键的那一对写成了 1-3 的距离。
                tighten(&mut b, i, j, lo, hi);
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
                let Some((cis, trans)) = torsion_span(
                    b.upper(iu, ku),
                    b.upper(ku, lu),
                    b.upper(lu, ju),
                    b.upper(iu, lu),
                    b.upper(ku, ju),
                ) else {
                    continue;
                };
                // **环上的扭转是锁死的,不是顺式到反式的全程。**
                // 光滑化推不出这件事(它只知道路径长度,不知道环是平的)。
                let arom = mol.bonds()[bidx].order == omgkit_core::BondOrder::Aromatic;
                // 芳环上的扭转是**确定值**(见 `planar_ring_torsion`),不是区间;
                // 解不掉的才退回按环尺寸的包络
                let (t_lo, t_hi) = planar_ring_torsion(&ring_sets, &aromatic_ring, i, k, l, j)
                    .map_or_else(
                        || torsion_envelope(shared_ring(&ring_sets, &[i, k, l, j]), arom),
                        |t| (t, t),
                    );
                // d² 随 cos(扭转) 单调 —— 端点取遍区间即可
                let f = |t: f64| {
                    let c = t.to_radians().cos();
                    let (a2, b2) = (cis * cis, trans * trans);
                    (a2 + (b2 - a2) * (1.0 - c) / 2.0).max(0.0).sqrt()
                };
                tighten(&mut b, iu, ju, f(t_lo), f(t_hi));
                kind[iu * n + ju] = Kind::B14;
                kind[ju * n + iu] = Kind::B14;
                st.n14 += 1;
            }
        }
    }

    // ---- vdW 地板:**只铺给没有拓扑约束的对** ----
    // 1-2/1-3/1-4 是被键角扭转钉死的,它们本来就可以比 vdW 之和近
    // (sp 中心的 1-3 只有 2.3 Å,而 C···C 的 0.75×vdW 是 2.55)——
    // 给它们铺地板会把界矩阵直接铺成矛盾的。
    for i in 0..n {
        for j in (i + 1)..n {
            if kind[i * n + j] != Kind::None {
                continue;
            }
            let d0 = VDW_FRAC
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
    if ri2 < 0.0 || rj2 < 0.0 {
        return None; // 三角形不成立,这条路径给不出信息
    }
    let (ri, rj) = (ri2.sqrt(), rj2.sqrt());
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
        assert!(
            worst < 0.05,
            "苯的 1-4 最宽 {worst:.4} Å —— 芳环扭转没被钉住"
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
