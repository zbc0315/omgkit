//! **手性:有符号四点体积,以及嵌入之后那一次全局定向。**
//!
//! # 为什么全局定向必须在这里定,不能交给精修
//!
//! 嵌入给出的坐标系定向是**任意的** —— 一半的概率拿到镜像(见
//! [`crate::linalg::symmetric_eigen`] 里符号约定那一段)。而翻转整体手性是一次
//! **反射**(`det = −1`),不在 `SO(3)` 的连通分支里:要用连续下降从一个手性走到
//! 它的镜像,必须经过"所有手性体积同时为 0"的构型,也就是把整个分子压平。
//! **下降法不会付这个势垒**,手性罚项的权重调多大都没用。
//!
//! 所以定向只能**离散地**定一次:两个镜像各算一遍总罚,取小的那个。
//! 这一步确定、便宜(只与手性中心数有关),而且必须做。
//!
//! (RDKit 一有手性中心就把分子嵌到四维,`Embedder.cpp:1632`,就是为了绕开这件事:
//! 四维里在 `(x₃, x₄)` 平面转 π 就把 `x₃` 送到 `−x₃`,而**四维两两距离精确不变** ——
//! 三维里的反射在四维里是一次免费的连续旋转。)
//!
//! # 符号约定
//!
//! 有符号体积取 `det[p₁−p₀, p₂−p₀, p₃−p₀]`,四个配体按**槽位顺序**给。
//! `det < 0` 对应 `@`([`ChiralTag::Ccw`](omgkit_core::ChiralTag::Ccw))、
//! `det > 0` 对应 `@@`([`Cw`](omgkit_core::ChiralTag::Cw))。
//!
//! 这与 `omgkit-depict` 的 `read_chirality` 是**同一个约定**,不是另起一套 ——
//! 那边由 `the_reference_tetrahedron_pins_the_sign` 钉住符号,并且有外部判官
//! (`harness/check_wedge_readback.py`,拿 RDKit 从导出的 molblock 反读)验过。
//! 两套约定并存迟早会在某个交界处翻号,所以这里明确复用它。
//!
//! # 抽中心那一半:**槽位约定是实测出来的,不是推的**

//!
//! "四个配体按什么槽位排"这个约定,排错的话**错法整批一致** —— 于是"符号正确率"
//! 要么 0% 要么 100%,两个数看起来都像"约定定死了"。这种错推不出来,只能拿真值比。
//!
//! 真值取自 `harness/dump_chirality.py`:每个中心的有符号体积**在真实三维构象上
//! 的实际符号**(不是标记推出来的号 —— 那正是待验的东西)。实测:
//!
//! | 标记 | 有符号体积 | 样本 |
//! |---|---|---|
//! | `@`(Ccw) | **负** | 22 / 22 |
//! | `@@`(Cw) | **正** | 17 / 17 |
//!
//! 与 `omgkit-depict` 的约定一致,而且那一致性写成了机器可验的断言。
//! 全量判官见 `examples/chiral_oracle.rs`:247 个中心,符号错 0、漏抽 0。
//! 变异验证过它抓得住整体翻号(那一下是 475/475 全错,不是"差不多对")。
//!
//! # 调用方要负责的那一条
//!
//! [`centers`] 要求**标记与当前键序一致**。`omgkit_chem::add_explicit_hs` 把补出来
//! 的氢追加到原子表末尾、并且**明确不碰 `chiral_tag`**,而 SMILES 里隐式氢占第 1 槽 ——
//! 所以"解析 SMILES → 补氢 → 直接调 [`centers`]"中间**缺一步槽位重排**。
//! 按一张四配位齐全的连接表建分子则天然满足(判官走的正是这条路)。

/// 一个四面体手性中心:四个配体的槽位顺序,以及有符号体积**该是什么号**。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Center {
    /// 中心原子。
    pub atom: u32,
    /// 四个配体原子,**按槽位顺序**(见模块文档)。
    pub ligands: [u32; 4],
    /// 目标符号:`-1.0` 对应 `@`、`+1.0` 对应 `@@`。
    pub sign: f64,
}

/// 从分子里抽出四面体手性中心。
///
/// # 前置条件 —— 调用方负责,判官验的就是这两条
///
/// 1. **氢必须已经是显式原子。** 中心得凑够四个邻居才谈得上四面体。
/// 2. **`chiral_tag` 必须与当前的键序一致**,也就是"按邻居迭代顺序取四个配体"
///    正好是标记所指的槽位顺序。
///
/// 第 2 条**不是自动成立的**。`omgkit_chem::add_explicit_hs` 把补出来的氢追加到
/// 原子表末尾,并且**明确不碰 `chiral_tag`**(它的文档写着"谁要在补氢之后用手性,
/// 必须自己把这一层想清楚")—— 而 SMILES 里隐式氢占的是第 1 槽。
/// 所以"解析 SMILES → 补氢 → 直接调这里"是**错**的,中间缺一步槽位重排。
///
/// 反过来,按一张**四配位齐全的连接表**建出来的分子天然满足这两条:
/// 没有隐式氢,键序就是槽位序。判官走的正是这条路
/// (`examples/chiral_oracle.rs`,真值取自真实构象的有符号体积)。
///
/// 抽不出四个邻居、或者标记不是四面体的原子,直接跳过 —— 不猜。
#[must_use]
pub fn centers(mol: &omgkit_core::MolBuilder) -> Vec<Center> {
    use omgkit_core::ChiralTag;
    let mut out = Vec::new();
    for (idx, a) in mol.atoms().iter().enumerate() {
        let sign = match a.chiral_tag {
            ChiralTag::Ccw => -1.0,
            ChiralTag::Cw => 1.0,
            _ => continue,
        };
        let Ok(id) = u32::try_from(idx) else { continue };
        let nb: Vec<u32> = mol.neighbors(id).map(|(y, _)| y).collect();
        let Ok(ligands) = <[u32; 4]>::try_from(nb.as_slice()) else {
            continue; // 不是四配位:氢没补、或者根本不是四面体中心
        };
        out.push(Center {
            atom: id,
            ligands,
            sign,
        });
    }
    out
}

/// 四个点围出的**有符号**体积(实为体积的 6 倍,只用它的号与相对大小)。
///
/// `det[p₁−p₀, p₂−p₀, p₃−p₀]`。只吃前三维 —— 即便将来在四维里精修,
/// 手性也只由前三维定(RDKit 的 `calcChiralVolume` 同此)。
#[must_use]
pub fn signed_volume(p0: [f64; 3], p1: [f64; 3], p2: [f64; 3], p3: [f64; 3]) -> f64 {
    let d = |p: [f64; 3]| [p[0] - p0[0], p[1] - p0[1], p[2] - p0[2]];
    let (a, b, c) = (d(p1), d(p2), d(p3));
    a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
}

/// 一个中心在给定坐标下的有符号体积。
///
/// # Panics
///
/// 配体下标越界时 panic。
#[must_use]
pub fn center_volume(coords: &[[f64; 3]], c: &Center) -> f64 {
    let p = |k: usize| coords[c.ligands[k] as usize];
    signed_volume(p(0), p(1), p(2), p(3))
}

/// 一组中心里有几个的号是对的。
#[must_use]
pub fn correct_count(coords: &[[f64; 3]], centers: &[Center]) -> usize {
    centers
        .iter()
        .filter(|c| {
            let v = center_volume(coords, c);
            v != 0.0 && v.signum() == c.sign
        })
        .count()
}

/// **嵌入之后定一次全局定向:要不要把整个结构镜像过来。**
///
/// 判据是"镜像之后号对的中心更多"。镜像把每个中心的有符号体积整体变号,
/// 所以这一步只是数一遍,不需要真的把坐标翻过去再算。
///
/// 平局时返回 `false`(不翻)——**必须有个确定的平局规则**,否则同一个分子
/// 两次跑可能给出互为镜像的答案,而这个 crate 承诺确定性。
#[must_use]
pub fn needs_reflection(coords: &[[f64; 3]], centers: &[Center]) -> bool {
    let ok = correct_count(coords, centers);
    // 镜像后每个体积变号,所以"号对的"恰好换成"号错且非零的"
    let nonzero = centers
        .iter()
        .filter(|c| center_volume(coords, c) != 0.0)
        .count();
    nonzero - ok > ok
}

/// 把整个结构镜像过来(翻第一根坐标轴)。
///
/// 翻哪一根都一样 —— 差别只是一次旋转,而旋转不改变分子。
pub fn reflect(coords: &mut [[f64; 3]]) {
    for p in coords.iter_mut() {
        p[0] = -p[0];
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 一个标准四面体:`p0` 在 +z,另外三个在下方平面上,按 0°/120°/240° 排。
    ///
    /// 从 `p0` 往中心看(即沿 −z 方向看),`p1 → p2 → p3` 是**逆时针**,
    /// 所以这就是 `@`([`ChiralTag::Ccw`](omgkit_core::ChiralTag::Ccw))。
    fn reference_tetrahedron() -> [[f64; 3]; 4] {
        let r = (8.0_f64).sqrt() / 3.0;
        let z = -1.0 / 3.0;
        [
            [0.0, 0.0, 1.0],
            [r, 0.0, z],
            [
                r * (2.0 * std::f64::consts::PI / 3.0).cos(),
                r * (2.0 * std::f64::consts::PI / 3.0).sin(),
                z,
            ],
            [
                r * (4.0 * std::f64::consts::PI / 3.0).cos(),
                r * (4.0 * std::f64::consts::PI / 3.0).sin(),
                z,
            ],
        ]
    }

    #[test]
    fn 参照四面体把符号钉死() {
        // **这条测试是整个模块的符号约定本身。** 它必须与 omgkit-depict 的
        // `the_reference_tetrahedron_pins_the_sign` 给出同一个号 —— 两套约定
        // 并存的话,迟早在某个交界处翻出对映体。
        let t = reference_tetrahedron();
        let v = signed_volume(t[0], t[1], t[2], t[3]);
        assert!(v < 0.0, "逆时针(@ / Ccw)的有符号体积应当是负的,实际 {v}");
        // 换两个配体 = 一次对换 = 手性翻转
        let v2 = signed_volume(t[0], t[2], t[1], t[3]);
        assert!(v2 > 0.0, "对换两个配体之后应当变号,实际 {v2}");
        assert!((v + v2).abs() < 1e-12, "一次对换只该变号,大小不变");

        // **把与 omgkit-depict 的一致性变成机器可验的,而不是一句注释。**
        // 下面这四个点是 depict 的 `the_reference_tetrahedron_pins_the_sign`
        // 里逐字写着的那一组(它断言 `det < 0` ⇒ Ccw)。两个 crate 各有一套
        // 约定却没人比过,迟早在交界处翻出对映体 —— 所以这里直接比。
        let depict_pts = [
            [0.0, 0.0, 1.0],
            [1.0, 0.0, -0.33],
            [-0.5, 0.866, -0.33],
            [-0.5, -0.866, -0.33],
        ];
        let dv = signed_volume(depict_pts[0], depict_pts[1], depict_pts[2], depict_pts[3]);
        assert!(
            dv < 0.0,
            "与 omgkit-depict 的符号约定对不上了:同一组点它判 Ccw(det<0),这里得 {dv}"
        );
    }

    #[test]
    fn 镜像把每个体积都变号() {
        let t = reference_tetrahedron();
        let mut m = t;
        reflect(&mut m);
        let (a, b) = (
            signed_volume(t[0], t[1], t[2], t[3]),
            signed_volume(m[0], m[1], m[2], m[3]),
        );
        assert!((a + b).abs() < 1e-12, "镜像后应当恰好变号:{a} vs {b}");
    }

    #[test]
    fn 旋转不改变符号() {
        // 有符号体积在**旋转**下不变、在**反射**下变号 —— 这正是它能当手性判据的理由。
        let t = reference_tetrahedron();
        let (c, s) = (0.6_f64, 0.8_f64);
        let rot = |p: [f64; 3]| [c * p[0] - s * p[1], s * p[0] + c * p[1], p[2]];
        let r = [rot(t[0]), rot(t[1]), rot(t[2]), rot(t[3])];
        let (a, b) = (
            signed_volume(t[0], t[1], t[2], t[3]),
            signed_volume(r[0], r[1], r[2], r[3]),
        );
        assert!((a - b).abs() < 1e-12, "旋转不该改变有符号体积:{a} vs {b}");
    }

    fn ctr(sign: f64) -> Center {
        Center {
            atom: 0,
            ligands: [0, 1, 2, 3],
            sign,
        }
    }

    #[test]
    fn 全局定向按多数决() {
        let t = reference_tetrahedron();
        let coords: Vec<[f64; 3]> = t.to_vec();
        // 参照四面体的号是负的
        assert!(
            !needs_reflection(&coords, &[ctr(-1.0)]),
            "号已经对了,不该翻"
        );
        assert!(needs_reflection(&coords, &[ctr(1.0)]), "号反了,应当翻");
        // 两个中心一对一错 —— 平局,规则是**不翻**(必须确定,不能随实现摆动)
        assert!(
            !needs_reflection(&coords, &[ctr(-1.0), ctr(1.0)]),
            "平局时的规则是不翻"
        );
        // 二比一
        assert!(
            needs_reflection(&coords, &[ctr(1.0), ctr(1.0), ctr(-1.0)]),
            "二比一应当翻"
        );
    }

    #[test]
    fn 压平的中心不算数() {
        // 四点共面 → 体积为 0 → 这个中心的号无从谈起,两边都不该把它计进去。
        // **这一条正是"判据只查符号会被压平刷绿"的根** —— 号为 0 时
        // `signum()` 给 0,与 ±1 都不等,所以它既不算对也不算错。
        let flat = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        ];
        assert_eq!(correct_count(&flat, &[ctr(1.0)]), 0);
        assert_eq!(correct_count(&flat, &[ctr(-1.0)]), 0);
        assert!(
            !needs_reflection(&flat, &[ctr(1.0)]),
            "全是零体积,翻了也没用"
        );
    }

    #[test]
    fn 没有手性中心时不翻() {
        let coords = vec![[0.0; 3]; 4];
        assert!(!needs_reflection(&coords, &[]));
    }
}
