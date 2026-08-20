//! **自穿检测** —— 链有没有从环里穿过去、两根键有没有交叉。
//!
//! # 为什么非有这一条不可
//!
//! 精修阶段有个决定:力场里放**全部** `N²` 对原子,而不是像 RDKit 那样只放
//! `u − l ≤ 5.0` 的对。理由是被滤掉的那些(实测占 16.8%,柔性大分子上最高 51%)
//! 正是拓扑上远、几何上可能撞在一起的对 —— 也就是自穿会发生的地方。
//!
//! **但那个改动的目的一直没有对应的判据。** 收益证明不了,也防不住后续改动把它
//! 弄丢。判据一(越界)看不见自穿:链穿过环时,每一对原子的距离都可以完全合法 ——
//! 穿过去的那根键与环上的键**不共享原子**,它们的距离约束只有一条很松的 vdW 下限,
//! 而"穿过"与"贴着"在两两距离上几乎没有区别。
//!
//! 所以要直接量**几何**:线段与线段、线段与环面。
//!
//! # 两个量
//!
//! - **键-键穿插**:两根**互相够远**的键(不共享原子,且拓扑距离 > [`RIGID_TOPO`]),
//!   算两条线段的最短距离。贴得太近就是穿了 —— 阈值不是拍的,见 [`CROSS_TOL`]。
//! - **环穿刺**:每个环按质心扇形三角化,数有多少条键的线段穿过任一三角形。
//!
//! # 环穿刺**会**假阳性 —— 这里先前写着"零假阳性",是错的
//!
//! 质心扇形三角化配 `.any()` 数的是"有没有交点",不是**交点的奇偶**。
//! 环非凸时链可以从凹口外穿过 z 平面、与两个扇形面各交一次 —— 相交 2 次,
//! mod-2 环绕数为 0(拓扑上没穿过去),`.any()` 照样报真。
//! 已复现:18 元月牙环(外弧 r=4、内弧 r=2.6,θ∈[−150°,150°],质心半径 0.238
//! 落在环外),链在半径 0.9 处竖直穿过 → 报 `pierces = 1`。
//!
//! 语料上没触发,因为里面的环都小到基本是凸的。**柔性大环、以及精修之前的
//! 嵌入坐标**最容易踩。修法是数交点奇偶而不是 `.any()`,尚未做 ——
//! 所以 `feasibility` 那边**故意没给环穿刺上闸**。
//!
//! 两个都是确定的,复杂度 `O(键²)`,56 原子的分子上几千次运算,可忽略。

use omgkit_core::MolBuilder;

/// 两根不共享原子的键靠到多近就算穿插(Å)。
///
/// **这个数是量出来的,不是拍的。** 拿语料里 MMFF 优化过的真实构象量:
/// 真实分子里不共享原子的键对,最短距离的**最小值**是多少?低于那个值的构型
/// 现实中不存在,所以阈值取在它下面就不会有假阳性。
///
/// 取 **1.1**。先前取 1.2 —— 那只比当时的实测最小值低 0.017 Å,换一批分子就可能误报,
/// 而**校准一旦误报,后面量什么都不作数**。物理上也说得通:两根键贴到 1.1 Å 以内,
/// 中间连一个氢都塞不下(H 的 vdW 半径 1.2 Å)。
///
/// # 头一版的实测下界 1.217 Å 是**假的** —— 它来自根本不可能互穿的键对
///
/// 设下界的那一档是 `[Na]N=[N+]=[N-]`:最近的键对是 `Na–N1` 与 `N2–N3`,
/// 中间隔着 `N1=N2` **一根键**。这一档的最短距离**恒等于中间那根键的键长**
/// (实测 1.20 / 1.14 / 1.10 / 1.13 / 1.15),与"穿没穿"毫无关系。
/// 四面体烷 `C12C3C1C23` 更极端:三对**对棱**几何上就相距 `a/√2 = 1.066 Å`
/// (`a` = C–C 1.508),已经低于这个阈值,而它们被环系锁死、无法互穿。
///
/// 所以修法不是挪这个常数,是**别把那一档算进来**(见 [`RIGID_TOPO`])。
/// 排除前后,400 个分子真实构象上够格的键对最近距离:
///
/// | | 最小 | p05 | 中位 |
/// |---|---|---|---|
/// | 排除前 | 1.217 | 1.286 | 1.352 |
/// | **排除后** | **1.290** | **1.981** | **2.289** |
///
/// 中位翻了近一倍 —— 旧口径里绝大多数"键对"根本就是刚性那一族。
/// 阈值 1.1 的余量因此从 10% 变成 **17%**,而检测器没有变瞎:
/// 同一批分子上我们自己嵌出来的坐标仍报 643 处交叉。
pub const CROSS_TOL: f64 = 1.1;

/// 两根键的原子集之间,拓扑距离小到多少就**不再算作可能互穿**。
///
/// `0` 是共享原子(那是键角);`1` 是被**一根键**连起来的两根键
/// (`A–B` 与 `C–D`,`B–C` 成键)。后者穿不过去 —— 要互穿,两条线段必须相交,
/// 而它们被中间那根键钉在一起,相交只能发生在退化几何上,不是拓扑上的穿插。
/// 而这一档恰恰是把 [`CROSS_TOL`] 的实测下界压到 1.217 Å 的那一档。
///
/// 取 `1`(即排除拓扑距离 ≤ 1 的键对)。`2` 就过头了:五、六元环上隔着两根键的
/// 两条边靠得也近,但那时"穿插"已经开始有意义了。
pub const RIGID_TOPO: u8 = 1;

/// 一个构型的自穿账。
///
/// # **只数同一个连通片里的**
///
/// 盐、共晶这类分子有多个互不相连的片段,而**片段之间摆在哪儿是另一个问题**,
/// 与"链穿过自己"没有关系,修法也不同(那是片段摆放,不是几何精修)。
///
/// 这不是理论上的顾虑:实测 RDKit 自己给多片段盐的构象里,
/// `[O-][N+]([O-])=O.F[Co+]12(F)(NCCN1)NCCN2` 的硝酸根与钴配合物**直接叠在一起**,
/// 两个原子只差 0.3 Å。把这种情形算进"自穿",会让检测器在**真实构象**上
/// 报出 48 次交叉 —— 然后校准这一步就废了,后面量什么都不作数。
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Threading {
    /// 够远的键对(见 [`RIGID_TOPO`])里,最短距离低于 [`CROSS_TOL`] 的对数。
    pub crossings: usize,
    /// **(键, 环) 对**的个数,不是键数 —— 一根键竖直穿过立方烷笼会记 2
    /// (上下两个面各一个环),稠合/桥环体系里同一次穿插会被多个环各记一次。
    pub pierces: usize,
    /// 所有够远的键对中,最短的那个距离(Å)。
    ///
    /// [`detect`] 在没有这样的键对时留 `f64::MAX`,但 [`Threading::default()`]
    /// 给的是 **0.0** —— 拿默认值当"没量到"会读成"两根键完全重合"。
    pub min_gap: f64,
    /// 检查了多少对键。分母 —— **没有它,`crossings = 0` 可能只是因为没在看**。
    pub pairs: usize,
}

/// 两条线段之间的最短距离。
///
/// 经典的夹紧法:先求两条**直线**上的最近点参数,夹到 `[0, 1]`,再互相回代一次。
/// 平行、退化成点这些情形都由分母的下限兜住。
#[must_use]
pub fn segment_distance(p1: [f64; 3], q1: [f64; 3], p2: [f64; 3], q2: [f64; 3]) -> f64 {
    let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let d1 = sub(q1, p1);
    let d2 = sub(q2, p2);
    let r = sub(p1, p2);
    let (a, e, f) = (dot(d1, d1), dot(d2, d2), dot(d2, r));
    const EPS: f64 = 1e-12;

    // 两条都退化成点
    if a <= EPS && e <= EPS {
        return dot(r, r).sqrt();
    }
    let (s, t);
    if a <= EPS {
        // 第一条退化成点
        s = 0.0;
        t = (f / e).clamp(0.0, 1.0);
    } else {
        let c = dot(d1, r);
        if e <= EPS {
            // 第二条退化成点
            t = 0.0;
            s = (-c / a).clamp(0.0, 1.0);
        } else {
            let b = dot(d1, d2);
            let denom = a * e - b * b;
            // `denom == 0` 就是两条平行 —— 这时 s 取 0,靠下面的回代把 t 摆对
            let s0 = if denom > EPS {
                ((b * f - c * e) / denom).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let t0 = (b * s0 + f) / e;
            // t 夹出界的话,把 s 按夹住的 t 重算一遍(这一步不能省,
            // 省了在"两条线段错开"的情形下会给出偏大的距离)
            if t0 < 0.0 {
                t = 0.0;
                s = (-c / a).clamp(0.0, 1.0);
            } else if t0 > 1.0 {
                t = 1.0;
                s = ((b - c) / a).clamp(0.0, 1.0);
            } else {
                t = t0;
                s = s0;
            }
        }
    }
    let c1 = [p1[0] + d1[0] * s, p1[1] + d1[1] * s, p1[2] + d1[2] * s];
    let c2 = [p2[0] + d2[0] * t, p2[1] + d2[1] * t, p2[2] + d2[2] * t];
    let d = sub(c1, c2);
    dot(d, d).sqrt()
}

/// 线段有没有穿过三角形(Möller–Trumbore,带线段参数的夹紧)。
#[must_use]
pub fn segment_hits_triangle(
    p: [f64; 3],
    q: [f64; 3],
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
) -> bool {
    let sub = |x: [f64; 3], y: [f64; 3]| [x[0] - y[0], x[1] - y[1], x[2] - y[2]];
    let dot = |x: [f64; 3], y: [f64; 3]| x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
    let cross = |x: [f64; 3], y: [f64; 3]| {
        [
            x[1] * y[2] - x[2] * y[1],
            x[2] * y[0] - x[0] * y[2],
            x[0] * y[1] - x[1] * y[0],
        ]
    };
    let dir = sub(q, p);
    let (e1, e2) = (sub(b, a), sub(c, a));
    let h = cross(dir, e2);
    let det = dot(e1, h);
    if det.abs() < 1e-12 {
        return false; // 与三角形所在平面平行
    }
    let inv = 1.0 / det;
    let s = sub(p, a);
    let u = inv * dot(s, h);
    if !(0.0..=1.0).contains(&u) {
        return false;
    }
    let qv = cross(s, e1);
    let v = inv * dot(dir, qv);
    if v < 0.0 || u + v > 1.0 {
        return false;
    }
    let t = inv * dot(e2, qv);
    // 交点要落在**线段**上,不是整条直线上
    (0.0..=1.0).contains(&t)
}

/// 连通片划分:同一个片段里的原子拿到同一个编号。
fn components(mol: &MolBuilder) -> Vec<usize> {
    let n = mol.num_atoms();
    let mut comp = vec![usize::MAX; n];
    let mut next = 0;
    for start in 0..n {
        if comp[start] != usize::MAX {
            continue;
        }
        let mut q = std::collections::VecDeque::from([start]);
        comp[start] = next;
        while let Some(x) = q.pop_front() {
            let Ok(xu) = u32::try_from(x) else { continue };
            for (y, _) in mol.neighbors(xu) {
                let y = y as usize;
                if y < n && comp[y] == usize::MAX {
                    comp[y] = next;
                    q.push_back(y);
                }
            }
        }
        next += 1;
    }
    comp
}

/// 数一个构型里的自穿。
///
/// # Panics
///
/// 坐标数与原子数对不上时 panic。
#[must_use]
pub fn detect(mol: &MolBuilder, coords: &[[f64; 3]]) -> Threading {
    assert_eq!(coords.len(), mol.num_atoms(), "坐标数与原子数对不上");
    let bonds: Vec<(usize, usize)> = mol
        .bonds()
        .iter()
        .map(|b| (b.begin as usize, b.end as usize))
        .collect();
    let mut t = Threading {
        min_gap: f64::MAX,
        ..Threading::default()
    };
    let comp = components(mol);
    let n = mol.num_atoms();
    // 拓扑距离,封顶 `RIGID_TOPO + 1` —— 只用来判"够不够远",不需要真值
    let cap = RIGID_TOPO + 1;
    let mut topo = vec![cap; n * n];
    for start in 0..n {
        let mut d = vec![u8::MAX; n];
        d[start] = 0;
        let mut q = std::collections::VecDeque::from([start]);
        while let Some(x) = q.pop_front() {
            if d[x] >= cap {
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
            topo[start * n + j] = d[j].min(cap);
        }
    }

    // ---- 键-键 ----
    for (x, &(i, j)) in bonds.iter().enumerate() {
        for &(k, l) in &bonds[(x + 1)..] {
            // **共享原子的键对不算** —— 它们本来就该挨着,那是键角不是穿插
            if i == k || i == l || j == k || j == l {
                continue;
            }
            // **跨片段的不算** —— 那是片段摆放的问题,见 `Threading` 的文档
            if comp[i] != comp[k] {
                continue;
            }
            // **被一根键连起来的两根键不算** —— 它们穿不过去,而这一档正是把
            // `CROSS_TOL` 的实测下界压到 1.217 Å 的那一档。见 `RIGID_TOPO`。
            let near = [(i, k), (i, l), (j, k), (j, l)]
                .iter()
                .map(|&(a, b)| topo[a * n + b])
                .min()
                .unwrap_or(u8::MAX);
            if near <= RIGID_TOPO {
                continue;
            }
            let d = segment_distance(coords[i], coords[j], coords[k], coords[l]);
            t.pairs += 1;
            t.min_gap = t.min_gap.min(d);
            if d < CROSS_TOL {
                t.crossings += 1;
            }
        }
    }

    // ---- 环穿刺 ----
    for ring in omgkit_chem::sssr::ring_set(mol) {
        let atoms: Vec<usize> = ring.atoms.iter().map(|a| *a as usize).collect();
        if atoms.len() < 3 {
            continue;
        }
        // 质心扇形三角化。环不是平面也没关系 —— 扇形三角面片合起来仍然把环口封住,
        // 而"穿过环口"正是要数的事。
        let n = atoms.len();
        #[allow(clippy::cast_precision_loss)]
        let nf = n as f64;
        let mut cen = [0.0; 3];
        for &a in &atoms {
            for k in 0..3 {
                cen[k] += coords[a][k] / nf;
            }
        }
        for (x, &(i, j)) in bonds.iter().enumerate() {
            let _ = x;
            // 环上的键、以及与环共享原子的键,不算穿刺
            if atoms.contains(&i) || atoms.contains(&j) {
                continue;
            }
            // 跨片段的不算,同上
            if comp[i] != comp[atoms[0]] {
                continue;
            }
            let hit = (0..n).any(|k| {
                segment_hits_triangle(
                    coords[i],
                    coords[j],
                    cen,
                    coords[atoms[k]],
                    coords[atoms[(k + 1) % n]],
                )
            });
            if hit {
                t.pierces += 1;
            }
        }
    }
    t
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn 线段距离的解析解() {
        // 两条正交且错开的线段:x 轴上的 [0,1] 与 z=1 处 y 轴上的 [0,1],
        // 最近点是 (0,0,0) 与 (0,0,1),距离 1
        let d = segment_distance(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
        );
        assert!((d - 1.0).abs() < 1e-12, "{d}");
        // 平行
        let d = segment_distance(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [1.0, 2.0, 0.0],
        );
        assert!((d - 2.0).abs() < 1e-12, "平行线段 {d}");
        // 共线但不重叠:[0,1] 与 [3,4],距离 2
        let d = segment_distance(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
        );
        assert!((d - 2.0).abs() < 1e-12, "共线 {d}");
        // 真的相交
        let d = segment_distance(
            [-1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 1.0, 0.0],
        );
        assert!(d < 1e-12, "相交的两条线段距离应当是 0,实得 {d}");
        // 退化成点
        let d = segment_distance([0.0; 3], [0.0; 3], [3.0, 4.0, 0.0], [3.0, 4.0, 0.0]);
        assert!((d - 5.0).abs() < 1e-12, "两个点 {d}");
    }

    /// 造一个分子:`n` 个原子,给定键表,给定坐标。
    fn mk(n: usize, bonds: &[(u32, u32)], xyz: &[[f64; 3]]) -> (MolBuilder, Vec<[f64; 3]>) {
        let mut m = MolBuilder::new();
        for _ in 0..n {
            m.add_atom_data(omgkit_core::AtomData::new(6));
        }
        for &(i, j) in bonds {
            m.add_bond(i, j, omgkit_core::BondOrder::Single).unwrap();
        }
        (m, xyz.to_vec())
    }

    #[test]
    fn 被一根键连起来的两根键不算交叉() {
        // `A–B–C–D` 四原子链,把它折到 `A–B` 与 `C–D` 只差 0.3 Å ——
        // 远低于 CROSS_TOL,但它们穿不过去,不该记成交叉。
        // 摆法:B、C 在 x 轴上相距 1;A 在 B 上方 0.3、C 上方对称放 D。
        let (m, xyz) = mk(
            4,
            &[(0, 1), (1, 2), (2, 3)],
            &[
                [-1.0, 0.3, 0.0], // A
                [0.0, 0.0, 0.0],  // B
                [1.0, 0.0, 0.0],  // C
                [2.0, 0.3, 0.0],  // D
            ],
        );
        let t = detect(&m, &xyz);
        // 先确认这一对**几何上**确实贴得够近 —— 否则这条测试测了个寂寞
        let d = segment_distance(xyz[0], xyz[1], xyz[2], xyz[3]);
        assert!(d < CROSS_TOL, "构型没摆够近({d}),这条测试白测");
        assert_eq!(t.crossings, 0, "被一根键连起来的两根键不该记成交叉");
        assert_eq!(t.pairs, 0, "这一对应当连查都不查");
    }

    #[test]
    fn 真正够远的两根键照样查得出来() {
        // 两条互不相连的链,各两个原子,摆成十字且只差 0.2 Å —— 必须报交叉。
        // **但它们要在同一连通片里**(跨片段不算),所以用一条长链把两端接起来:
        // 0–1 与 4–5 是要查的两根键,中间靠 1–2–3–4 连着(拓扑距离 ≥ 2)。
        let (m, xyz) = mk(
            6,
            &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)],
            &[
                [-1.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [3.0, 2.0, 0.0],
                [3.0, 6.0, 0.0],
                [0.0, -1.0, 0.2],
                [0.0, 1.0, 0.2],
            ],
        );
        let t = detect(&m, &xyz);
        assert!(t.pairs > 0, "一对都没查,那个 0 只说明没在看");
        assert!(
            t.crossings >= 1,
            "0–1 与 4–5 只差 0.2 Å 且拓扑上隔着 3 根键,必须报交叉;实得 {t:?}"
        );
    }

    #[test]
    fn 四面体烷的对棱不算交叉() {
        // 正四面体的三对**对棱**互相垂直、相距 `a/√2`。取 a = 1.508(C–C),
        // 对棱距离 1.066 Å < CROSS_TOL —— 而它们被环系锁死,无法互穿。
        let a = 1.508;
        let s = a / 2.0_f64.sqrt() / 2.0;
        let (m, xyz) = mk(
            4,
            &[(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)],
            &[
                [-a / 2.0, 0.0, -s],
                [a / 2.0, 0.0, -s],
                [0.0, -a / 2.0, s],
                [0.0, a / 2.0, s],
            ],
        );
        let d = segment_distance(xyz[0], xyz[1], xyz[2], xyz[3]);
        assert!(
            (d - a / 2.0_f64.sqrt()).abs() < 1e-9,
            "对棱距离该是 a/√2 = {},实得 {d}",
            a / 2.0_f64.sqrt()
        );
        assert!(d < CROSS_TOL, "对棱距离 {d} 该低于阈值,否则这条测试白测");
        assert_eq!(detect(&m, &xyz).crossings, 0, "四面体烷的对棱不该记成交叉");
    }

    #[test]
    fn 错开的线段不能给出偏大的距离() {
        // 这一组专治"t 夹出界之后不回代 s"那个错:两条线段在参数域上错开,
        // 不回代的话会拿端点硬算,给出偏大的值。
        // x 轴 [0,1] 与 从 (2,0,0) 到 (2,0,1) 的竖直段 —— 最近点是 (1,0,0) 与 (2,0,0),距离 1
        let d = segment_distance(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 0.0, 1.0],
        );
        assert!((d - 1.0).abs() < 1e-12, "{d}");
    }

    #[test]
    fn 线段穿三角形() {
        let (a, b, c) = ([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]);
        // 从上往下穿过三角形内部
        assert!(segment_hits_triangle(
            [0.5, 0.5, 1.0],
            [0.5, 0.5, -1.0],
            a,
            b,
            c
        ));
        // 从三角形外面过
        assert!(!segment_hits_triangle(
            [5.0, 5.0, 1.0],
            [5.0, 5.0, -1.0],
            a,
            b,
            c
        ));
        // 方向对但线段太短,够不着平面
        assert!(!segment_hits_triangle(
            [0.5, 0.5, 1.0],
            [0.5, 0.5, 0.5],
            a,
            b,
            c
        ));
        // 与平面平行
        assert!(!segment_hits_triangle(
            [0.5, 0.5, 1.0],
            [1.5, 0.5, 1.0],
            a,
            b,
            c
        ));
    }
}
