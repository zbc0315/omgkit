//! Gasteiger–Marsili 部分电荷(PEOE,逐轨道电负性均衡)。
//!
//! 原始文献:J. Gasteiger, M. Marsili, *Tetrahedron* **36** (1980) 3219,
//! "Iterative Equalization of Orbital Electronegativity — A Rapid Access to
//! Atomic Charges"。
//!
//! # 算法一句话
//!
//! 每个原子的轨道电负性写成自身电荷的二次式 `E = a + b·q + c·q²`。相邻原子
//! 电负性不等时电荷沿键流动,流量按较**正**一端的电离能标度;每轮流量乘一个
//! 阻尼因子(首轮 0.5,此后逐轮减半),迭代 12 轮收敛。
//!
//! # 隐式氢不在图里,但参与迭代
//!
//! 一个重原子带的 n 个隐式氢没有各自的节点。它们彼此不相邻、只连着同一个重
//! 原子,所以可以**当作一个整体**在同一轮里更新:`h_charges[i]` 记的是第 i
//! 个重原子所带全部隐式氢的电荷之和。这不是近似,是把 n 个完全等价的节点合
//! 并同类项。
//!
//! # 为什么会算出 NaN,以及为什么不许把它藏起来
//!
//! 参数表只覆盖有机合成常见的十几种元素。表外元素(多数金属)拿到的是零参数,
//! 于是它的电离能标度 `ionx` 为 0 —— 若一根键的两端都是表外元素,流量的分母
//! 就是 0,结果为 `inf` 或 `NaN`,并沿图扩散。
//!
//! RDKit 的行为与此逐位相同(`throwOnParamFailure=false` 时落到 `X *` 那行零
//! 参数)。**这里照抄,不做兜底** —— 把 NaN 换成 0 会让"这个原子算不出电荷"
//! 和"这个原子电荷恰好是 0"变成同一件事,而下游的特征化恰恰要靠这个区别决定
//! 该不该屏蔽这一维。调用方用 [`f64::is_finite`] 自己判。
//!
//! # 参数表溯源
//!
//! 转录自 RDKit `Code/GraphMol/PartialCharges/GasteigerParams.cpp`
//! (BSD-3-Clause),含该文件的 `defaultParamData` 与 `additionalParamData`
//! 两段;后者由该文件注明是用 PyBabel 的算法补算的。表里没有的
//! (元素, 模式) 组合落到 `X *` 那一行,即零参数。
//!
//! 转录对不对不靠人眼核:`harness/check_descriptors.py` 拿全语料
//! 8831 个分子逐原子与 RDKit 的 `ComputeGasteigerCharges` 比,任何一个数字
//! 抄错都会在那里变红。

use omgkit_core::{element, Hybridization, MolBuilder};

/// 氢的电离能标度。其余元素由三个参数之和现算,唯独氢是个常数。
const IONX_H: f64 = 20.02;
/// 阻尼因子的初值
const DAMP: f64 = 0.5;
/// 每轮之后阻尼因子乘上它
const DAMP_SCALE: f64 = 0.5;

/// 迭代轮数。文献与 RDKit 的默认值都是 12。
pub const DEFAULT_ITERATIONS: usize = 12;

/// `(元素符号, 模式, a, b, c)`。模式即杂化,硫另有两个氧化态档。
///
/// **不含 `X *` 那一行**(全零)—— 它在 RDKit 里是查表失败的兜底,写在表里
/// 会让查表自己引用自己。兜底直接由 [`params_for`] 返回零。
static PARAMS: &[(&str, &str, f64, f64, f64)] = &[
    // defaultParamData
    ("H", "*", 7.17, 6.24, -0.56),
    ("C", "sp3", 7.98, 9.18, 1.88),
    ("C", "sp2", 8.79, 9.32, 1.51),
    ("C", "sp", 10.39, 9.45, 0.73),
    ("N", "sp3", 11.54, 10.82, 1.36),
    ("N", "sp2", 12.87, 11.15, 0.85),
    ("N", "sp", 15.68, 11.7, -0.27),
    ("O", "sp3", 14.18, 12.92, 1.39),
    ("O", "sp2", 17.07, 13.79, 0.47),
    ("F", "sp3", 14.66, 13.85, 2.31),
    ("Cl", "sp3", 11.00, 9.69, 1.35),
    ("Br", "sp3", 10.08, 8.47, 1.16),
    ("I", "sp3", 9.9, 7.96, 0.96),
    ("S", "sp3", 10.14, 9.13, 1.38),
    ("S", "so", 10.14, 9.13, 1.38),
    ("S", "so2", 12.00, 10.81, 1.20),
    ("S", "sp2", 10.88, 9.49, 1.33),
    ("P", "sp3", 8.90, 8.24, 0.96),
    // additionalParamData
    ("P", "sp2", 9.665, 8.530, 0.735),
    ("Si", "sp3", 7.300, 6.567, 0.657),
    ("Si", "sp2", 7.905, 6.748, 0.443),
    ("Si", "sp", 9.065, 7.027, -0.002),
    ("B", "sp3", 5.980, 6.820, 1.605),
    ("B", "sp2", 6.420, 6.807, 1.322),
    ("Be", "sp3", 3.845, 6.755, 3.165),
    ("Be", "sp2", 4.005, 6.725, 3.035),
    ("Mg", "sp2", 3.565, 5.572, 2.197),
    ("Mg", "sp3", 3.300, 5.587, 2.447),
    ("Mg", "sp", 4.040, 5.472, 1.823),
    ("Al", "sp3", 5.375, 4.953, 0.867),
    ("Al", "sp2", 5.795, 5.020, 0.695),
];

/// 查参数;查不到给零参数(即 RDKit 的 `X *` 兜底行)。
fn params_for(symbol: &str, mode: &str) -> [f64; 3] {
    for &(s, m, a, b, c) in PARAMS {
        if s == symbol && m == mode {
            return [a, b, c];
        }
    }
    [0.0, 0.0, 0.0]
}

/// 原子的查表模式。
///
/// 杂化是 sp/sp²/sp³ 时直接用它。**其余情况只有两个元素有说法**:氢用 `*`;
/// 硫按邻接的氧个数落到 `so2` / `so` / `sp3`。别的元素落到空串,查不到 ——
/// 那正是零参数的入口,不是漏写。
fn mode_of(mol: &MolBuilder, idx: u32) -> &'static str {
    let atom = &mol.atoms()[idx as usize];
    match atom.hybridization {
        Hybridization::Sp3 => "sp3",
        Hybridization::Sp2 => "sp2",
        Hybridization::Sp => "sp",
        _ => {
            if atom.atomic_num == 1 {
                "*"
            } else if atom.atomic_num == 16 {
                let n_oxygen = mol
                    .neighbors(idx)
                    .filter(|&(nbr, _)| mol.atoms()[nbr as usize].atomic_num == 8)
                    .count();
                match n_oxygen {
                    2 => "so2",
                    1 => "so",
                    _ => "sp3",
                }
            } else {
                ""
            }
        }
    }
}

/// 迭代开始前,把共轭体系里的形式电荷摊到同种原子上。
///
/// 例:苯脒的两个氮各带 0.5 起步,而不是一个 1、一个 0。少了这一步,同一个
/// 分子写成两个共振式会给出两组不同的电荷 —— 而它们是同一个分子。
///
/// "同种"只看原子序数,"同一体系"只看**隔着两根共轭键**能不能走到,不是全局
/// 的离域计算。这是 RDKit `splitChargeConjugated` 的口径,照抄。
fn split_charge_conjugated(mol: &MolBuilder, charges: &mut [f64]) {
    let n = mol.num_atoms();
    let mut marker: Vec<usize> = Vec::new();
    for aix in 0..n {
        let atom = &mol.atoms()[aix];
        let mut formal = f64::from(atom.formal_charge);
        if formal.abs() <= f64::EPSILON || charges[aix].abs() >= f64::EPSILON {
            continue;
        }
        marker.clear();
        marker.push(aix);
        for (aax, b1) in mol.neighbors(u32::try_from(aix).unwrap_or(0)) {
            if !mol.bonds()[b1 as usize]
                .flags
                .contains(omgkit_core::BondFlags::CONJUGATED)
            {
                continue;
            }
            for (yax, b2) in mol.neighbors(aax) {
                if b2 == b1 {
                    continue;
                }
                if !mol.bonds()[b2 as usize]
                    .flags
                    .contains(omgkit_core::BondFlags::CONJUGATED)
                {
                    continue;
                }
                if mol.atoms()[yax as usize].atomic_num == atom.atomic_num {
                    formal += f64::from(mol.atoms()[yax as usize].formal_charge);
                    marker.push(yax as usize);
                }
            }
        }
        #[allow(clippy::cast_precision_loss)]
        let share = formal / marker.len() as f64;
        for &m in &marker {
            charges[m] = share;
        }
    }
}

/// 算出每个原子的 Gasteiger 部分电荷。
///
/// 返回的向量与 `mol.atoms()` 等长、同序。**隐式氢的电荷不计在所属重原子头上**
/// —— 与 RDKit 的 `_GasteigerCharge` 同一口径(它把隐式氢那份单独记在
/// `_GasteigerHCharge` 里)。
///
/// 表外元素上会返回非有限值,见模块文档。
///
/// # 前置
///
/// 分子必须**净化过**:要读杂化、共轭标志和隐式氢数,三样都只有净化才会填。
/// 拿一个没净化的分子进来不会报错,只会让每个数悄悄偏移。
#[must_use]
pub fn gasteiger_charges(mol: &MolBuilder, n_iter: usize) -> Vec<f64> {
    let n = mol.num_atoms();
    let mut charges = vec![0.0f64; n];
    if n == 0 {
        return charges;
    }
    split_charge_conjugated(mol, &mut charges);

    // 每个原子的三参数与电离能标度
    let mut atm_ps = Vec::with_capacity(n);
    let mut ionx = Vec::with_capacity(n);
    for (idx, atom) in mol.atoms().iter().enumerate() {
        let z = atom.atomic_num;
        let symbol = element::by_atomic_num(z).map_or("*", |e| e.symbol);
        let p = params_for(symbol, mode_of(mol, u32::try_from(idx).unwrap_or(0)));
        ionx.push(if z == 1 { IONX_H } else { p[0] + p[1] + p[2] });
        atm_ps.push(p);
    }

    // 隐式氢:每个重原子一格,存它带的全部隐式氢的电荷之和
    let mut h_charges = vec![0.0f64; n];
    let h_params = params_for("H", "*");
    let mut energy = vec![0.0f64; n];
    let mut damp = DAMP;

    for _ in 0..n_iter {
        for idx in 0..n {
            let p = atm_ps[idx];
            energy[idx] = p[0] + charges[idx] * (p[1] + p[2] * charges[idx]);
        }
        for idx in 0..n {
            let mut dq = 0.0;
            for (nbr, _) in mol.neighbors(u32::try_from(idx).unwrap_or(0)) {
                let nbr = nbr as usize;
                let dx = energy[nbr] - energy[idx];
                // 流量按较**正**一端的电离能标度。写法逐字照抄 RDKit:
                // `(sgn·(ionx[i] − ionx[j])) + ionx[j]`。数学上 sgn=1 时它就
                // 是 `ionx[i]`,但浮点里 `(a − b) + b ≠ a` —— 化简会让结果
                // 与参照差在末几位,而判据比的正是同一批数。
                let sgn = f64::from(u8::from(dx >= 0.0));
                dq += dx / (sgn * (ionx[idx] - ionx[nbr]) + ionx[nbr]);
            }
            let n_hs = u32::from(mol.atoms()[idx].num_explicit_hs)
                + u32::from(mol.atoms()[idx].num_implicit_hs);
            if n_hs > 0 {
                let n_hs_f = f64::from(n_hs);
                let q_hs = h_charges[idx] / n_hs_f;
                let e_h = h_params[0] + q_hs * (h_params[1] + h_params[2] * q_hs);
                let dx = e_h - energy[idx];
                let sgn = f64::from(u8::from(dx >= 0.0));
                let dq_h = dx / (sgn * (ionx[idx] - IONX_H) + IONX_H);
                dq += n_hs_f * dq_h;
                // 这些氢彼此不相邻,可以与重原子在同一轮里一起更新
                h_charges[idx] -= n_hs_f * dq_h * damp;
            }
            charges[idx] += damp * dq;
        }
        damp *= DAMP_SCALE;
    }
    charges
}

#[cfg(test)]
mod tests {
    use super::*;
    use omgkit_io::smiles;

    fn sanitized(smi: &str) -> MolBuilder {
        let mut m = smiles::parse(smi).unwrap_or_else(|e| panic!("{smi}: {}", e.render()));
        crate::pipeline::sanitize(&mut m).unwrap_or_else(|e| panic!("{smi}: {e}"));
        m
    }

    #[test]
    fn methane_carbon_is_negative_and_sums_to_zero() {
        let mol = sanitized("C");
        let q = gasteiger_charges(&mol, DEFAULT_ITERATIONS);
        assert_eq!(q.len(), 1);
        // 碳比氢电负,拿到负电荷
        assert!(q[0] < 0.0, "甲烷的碳应带负电,实得 {}", q[0]);
    }

    /// 与 RDKit 2025.09.2 的 `ComputeGasteigerCharges` 逐位对上。
    ///
    /// 全语料判据在 `harness/check_descriptors.py`(要 RDKit),这里只钉几条
    /// 走不同分支的:烷烃、羧酸(sp²/sp³ 氧)、砜(硫的 `so2` 档)、
    /// 胍鎓(共轭体系摊电荷 —— 三个氮必须相等)。**不装 RDKit 也能跑**,
    /// 所以 CI 的主 job 也守得住。
    #[test]
    fn charges_match_the_external_reference() {
        let cases: &[(&str, &[f64])] = &[
            ("C", &[-0.077_558]),
            ("CCO", &[-0.041_838, 0.040_221, -0.396_664]),
            ("CC(=O)O", &[0.033_768, 0.299_685, -0.252_820, -0.481_433]),
            (
                "CS(=O)(=O)C",
                &[0.038_503, 0.144_104, -0.229_414, -0.229_414, 0.038_503],
            ),
            (
                "NC(N)=[NH2+]",
                &[-0.291_178, 0.335_948, -0.291_178, -0.291_178],
            ),
        ];
        for &(smi, want) in cases {
            let got = gasteiger_charges(&sanitized(smi), DEFAULT_ITERATIONS);
            assert_eq!(got.len(), want.len(), "{smi}:原子数对不上");
            for (i, (&g, &w)) in got.iter().zip(want).enumerate() {
                assert!(
                    (g - w).abs() < 5e-7,
                    "{smi} 第 {i} 个原子:本实现 {g},参照 {w}"
                );
            }
        }
    }

    #[test]
    fn an_element_outside_the_table_yields_a_non_finite_charge() {
        // 两端都是表外元素 ⇒ 分母为 0。这条钉的是"不做兜底"这个决定本身:
        // 谁要是加一行 `if !q.is_finite() { q = 0.0 }`,这里立刻红。
        let mol = sanitized("[Na][Na]");
        let q = gasteiger_charges(&mol, DEFAULT_ITERATIONS);
        assert!(
            q.iter().any(|v| !v.is_finite()),
            "表外元素应给出非有限值,实得 {q:?}"
        );
    }
}
