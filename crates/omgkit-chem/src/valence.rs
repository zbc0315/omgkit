//! 价键计算与隐式氢推断(净化第 3 步)。
//!
//! 严格模式:超价直接判为失败,整条净化终止。
//!
//! # 三个容易写错的地方
//!
//! **1. `dv` 与 `valens` 取自不同的有效原子序数。**
//! `dv`(默认价)用**调整前**的有效原子序数取,超价元素调整之后再用**调整后**
//! 的取 `valens`(允许的价态表)。两步顺序颠倒,结果就不同。
//!
//! **2. 无价约束的元素不进芳香分支。**
//! 这类元素的默认价是 `-1`,判据里的 `dv >= 0 &&` 就是把它们挡在外面。
//! 少了这个条件,过渡金属之类会误入只对有机元素成立的芳香价态修正。
//!
//! **3. 配位键的价贡献不对称。**
//! 对起点(给体)算 0,对终点(受体)算 1。见
//! [`BondData::valence_contribution_to`](omgkit_core::BondData::valence_contribution_to)。
//!
//! # 自由基电子数
//!
//! 隐式氢的推断要读自由基电子数,而它由第 6 步 [`crate::radicals`] 填充 ——
//! 排在本步**之后**。所以净化管线中本步看到的必然是 0。
//!
//! 即便如此,这里也**读字段而不是写死 0**:本函数是公开 API,在净化完成之后
//! 再调一次是正常用法,那时拿到的才应是正确结果。写死 0 能在管线内蒙对,
//! 却会在管线外静默给出错误的氢数。

use omgkit_core::{element, AtomFlags, MolBuilder};

/// 第 3 步的产出。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ValenceResult {
    /// 每原子的显式价
    pub explicit_valence: Vec<i32>,
    /// 每原子的隐式氢数
    pub implicit_hs: Vec<u8>,
}

/// 价键校验失败。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ValenceError {
    /// 出问题的原子下标
    pub atom: u32,
    /// 元素符号
    pub symbol: &'static str,
    /// 算出的显式价
    pub valence: i32,
    /// 具体原因
    pub kind: ValenceErrorKind,
}

/// 价键校验失败的类别。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum ValenceErrorKind {
    /// 显式价超出该元素允许的最大值
    ExplicitValenceTooHigh,
    /// 芳香原子的显式价不等于任何允许的价态
    AromaticValenceNotAllowed,
    /// 形式电荷不合理(仅 H 的特殊分支)
    UnreasonableFormalCharge,
}

impl core::fmt::Display for ValenceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let what = match self.kind {
            ValenceErrorKind::ExplicitValenceTooHigh => "显式价超出允许范围",
            ValenceErrorKind::AromaticValenceNotAllowed => "芳香原子的显式价不在允许的价态中",
            ValenceErrorKind::UnreasonableFormalCharge => "形式电荷不合理",
        };
        write!(
            f,
            "原子 #{}({}):{},显式价 = {}",
            self.atom, self.symbol, what, self.valence
        )
    }
}

impl std::error::Error for ValenceError {}

/// 计算全部原子的显式价与隐式氢数,并就地写回 `num_implicit_hs`。
///
/// # Errors
/// 任一原子超价时返回 [`ValenceError`],此时整条净化应当判为失败。
pub fn update_property_cache(mol: &mut MolBuilder) -> Result<ValenceResult, ValenceError> {
    let n = mol.num_atoms();
    let mut explicit_valence = vec![0i32; n];
    let mut implicit_hs = vec![0u8; n];

    for i in 0..n as u32 {
        let ev = explicit_valence_of(mol, i, true)?;
        explicit_valence[i as usize] = ev;
        let ih = implicit_hs_of(mol, i, ev, true)?;
        implicit_hs[i as usize] = ih;
    }

    for (i, &h) in implicit_hs.iter().enumerate() {
        if let Some(a) = mol.atom_mut(i as u32) {
            a.num_implicit_hs = h;
        }
    }

    Ok(ValenceResult {
        explicit_valence,
        implicit_hs,
    })
}

// ---------------------------------------------------------------------------

/// 原子自身带芳香标志,**或**任一关联键是芳香键。
fn is_aromatic_atom(mol: &MolBuilder, idx: u32) -> bool {
    if mol.atoms()[idx as usize]
        .flags
        .contains(AtomFlags::AROMATIC)
    {
        return true;
    }
    mol.neighbors(idx).any(|(_, bi)| {
        let b = mol.bonds()[bi as usize];
        b.flags.contains(omgkit_core::BondFlags::AROMATIC)
            || b.order == omgkit_core::BondOrder::Aromatic
    })
}

/// 有效原子序数:`clamp(Z - 形式电荷, 0, 最大原子序数)`。
/// 形式电荷从 `from` 变到 `to` 时,这个元素**能用的价**变了多少。
///
/// # 为什么不能写成"负电荷减一、正电荷加一"
///
/// 价由**有效原子序数**定(带电原子按 Z∓q 那个元素的价表算),所以电荷对价的
/// 作用完全取决于元素:
///
/// | | 中性 | 带电 | 变化 |
/// |---|---|---|---|
/// | N⁺ | 3 | 4 | **+1** |
/// | O⁻ | 2 | 1 | −1 |
/// | C⁻ | 4 | 3 | −1 |
/// | **C⁺** | 4 | 3 | **−1**(与 N⁺ 反向) |
/// | **B⁻** | 3 | 4 | **+1**(与 O⁻ 反向) |
///
/// 按"负减正加"写死的话,碳正离子与硼酸根这两类会算反。返回 0 表示这个元素
/// 没有价约束(或电荷离谱到查不出价),调用方应当把它当作"说不准"。
#[must_use]
pub fn valence_shift(z: u8, from: i8, to: i8) -> i32 {
    let ovalens = valences_of(z);
    if ovalens.len() == 1 && ovalens[0] == -1 {
        return 0; // 无价约束的元素
    }
    let at = |q: i8| -> Option<i32> {
        let eff = effective_atomic_num(z, q);
        if eff == 0 {
            return None;
        }
        let v = default_valence(eff);
        if v < 0 {
            None
        } else {
            Some(v)
        }
    };
    match (at(from), at(to)) {
        (Some(a), Some(b)) => b - a,
        _ => 0,
    }
}

fn effective_atomic_num(z: u8, charge: i8) -> u8 {
    let max = (element::count() - 1) as i32;
    (i32::from(z) - i32::from(charge)).clamp(0, max) as u8
}

/// 带负电的 P/S/As/Se 可以保留超价形式,尽管它们与不支持超价的
/// Cl/Ar、Br/Kr 等电子。
fn can_be_hypervalent(z: u8, eff_z: u8) -> bool {
    (eff_z > 16 && (z == 15 || z == 16)) || (eff_z > 34 && (z == 33 || z == 34))
}

fn valences_of(z: u8) -> &'static [i8] {
    element::by_atomic_num(z).map_or(&[-1][..], |e| e.valences)
}

/// 价列表的首项,即默认价。
fn default_valence(z: u8) -> i32 {
    valences_of(z).first().map_or(-1, |&v| i32::from(v))
}

/// 价列表的末项,即最大允许价。
fn last_valence(z: u8) -> i32 {
    valences_of(z).last().map_or(-1, |&v| i32::from(v))
}

/// 非严格模式的显式价:跳过超价校验,永不失败。
///
/// 净化第 1 步全程用它 —— 那一步跑在价键校验之前,分子里本就存在
/// `N(=O)=O` 这类超价写法,而它的职责恰恰是把这些写法修成合法形式。
pub(crate) fn explicit_valence_nonstrict(mol: &MolBuilder, idx: u32) -> i32 {
    explicit_valence_of(mol, idx, false).unwrap_or(0)
}

/// 计算显式价。`strict` 为假时跳过超价校验,永不返回 `Err`。
fn explicit_valence_of(mol: &MolBuilder, idx: u32, strict: bool) -> Result<i32, ValenceError> {
    let atom = mol.atoms()[idx as usize];
    let z = atom.atomic_num;

    let mut accum: f32 = mol
        .neighbors(idx)
        .map(|(_, bi)| mol.bonds()[bi as usize].valence_contribution_to(idx))
        .sum();
    accum += f32::from(atom.num_explicit_hs);

    let ovalens = valences_of(z);
    // 只有当该元素本身有价约束时才启用"有效原子序数"(即扣除形式电荷)
    let eff_z = if ovalens.len() > 1 || ovalens[0] != -1 {
        effective_atomic_num(z, atom.formal_charge)
    } else {
        z
    };
    let dv = default_valence(eff_z);
    let valens = valences_of(eff_z);

    // `dv >= 0` 把无价约束的元素挡在外面,见模块文档第 2 点
    if dv >= 0 && accum > dv as f32 && is_aromatic_atom(mol, idx) {
        let mut pval = dv;
        for &val in valens {
            let val = i32::from(val);
            if val == -1 || val as f32 > accum {
                break;
            }
            pval = val;
        }
        // 差在 1.5 以内就取该价态 —— 针对 c1cccn1C 这类"按芳香键算像 4 价、
        // kekulize 之后其实是 3 价"的芳香原子
        if accum - pval as f32 <= 1.5 {
            accum = pval as f32;
        }
    }

    // x.5 要向上取整(1.5 → 2):加 0.1 再四舍五入
    accum += 0.1;
    let res = accum.round() as i32;

    // -- 严格校验(strict = false 时整段跳过,与 C++ 的 `if (strict || checkIt)` 一致)--
    if !strict {
        return Ok(res);
    }
    let mut max_valence = last_valence(eff_z);
    let mut offset = 0i32;
    if can_be_hypervalent(z, eff_z) {
        max_valence = last_valence(z);
        offset -= i32::from(atom.formal_charge);
    }
    // 历史遗留:双配位的 [H-] 一直被接受
    if z == 1 && atom.formal_charge == -1 {
        max_valence = 2;
    }
    // max_valence == -1 表示高端不设限
    if max_valence >= 0 && last_valence(z) >= 0 && (res + offset) > max_valence {
        return Err(ValenceError {
            atom: idx,
            symbol: element::by_atomic_num(z).map_or("?", |e| e.symbol),
            valence: res,
            kind: ValenceErrorKind::ExplicitValenceTooHigh,
        });
    }

    Ok(res)
}

/// 非严格模式的隐式氢数:跳过校验,永不失败。
pub(crate) fn implicit_hs_nonstrict(mol: &MolBuilder, idx: u32, ev: i32) -> u8 {
    implicit_hs_of(mol, idx, ev, false).unwrap_or(0)
}

/// 非严格模式的总价 = 显式价 + 隐式氢数。
///
/// Kekulize 用它做前后快照比对:kekulize 不应改变任何原子的总价。
pub(crate) fn total_valence_nonstrict(mol: &MolBuilder, idx: u32) -> i32 {
    let ev = explicit_valence_nonstrict(mol, idx);
    ev + i32::from(implicit_hs_nonstrict(mol, idx, ev))
}

/// 计算隐式氢数。`strict` 为假时不返回错误。
fn implicit_hs_of(mol: &MolBuilder, idx: u32, ev: i32, strict: bool) -> Result<u8, ValenceError> {
    let atom = mol.atoms()[idx as usize];
    if atom.flags.contains(AtomFlags::NO_IMPLICIT) {
        return Ok(0);
    }
    let z = atom.atomic_num;
    if z == 0 {
        return Ok(0); // 通配原子 `*`
    }

    // 自由基电子数由第 6 步 [`assign_radicals`](crate::radicals::assign_radicals)
    // 填充。净化管线里第 6 步排在第 3 步之后,所以这里读到的必然是 0;
    // 但读字段而不是写死 0 —— 用户在净化之后再调一次本函数时,
    // 拿到的才是正确结果。
    let n_radicals = i32::from(atom.num_radical_electrons);

    // 氢的特殊分支
    if ev == 0 && n_radicals == 0 && z == 1 {
        return match atom.formal_charge {
            1 | -1 => Ok(0),
            0 => Ok(1),
            _ if strict => Err(ValenceError {
                atom: idx,
                symbol: "H",
                valence: ev,
                kind: ValenceErrorKind::UnreasonableFormalCharge,
            }),
            _ => Ok(0),
        };
    }

    let mut explicit_plus_rad = ev + n_radicals;

    let ovalens = valences_of(z);
    let mut eff_z = if ovalens.len() > 1 || ovalens[0] != -1 {
        effective_atomic_num(z, atom.formal_charge)
    } else {
        z
    };
    if eff_z == 0 {
        return Ok(0);
    }

    // 注意:`dv` 取自**调整前**的 eff_z,`valens` 取自**调整后**的 —— 见模块文档第 1 点
    let dv = default_valence(eff_z);
    if dv == -1 {
        return Ok(0); // d 区 / f 区元素无默认价
    }
    if can_be_hypervalent(z, eff_z) {
        eff_z = z;
        explicit_plus_rad -= i32::from(atom.formal_charge);
    }
    let valens = valences_of(eff_z);

    let res: i32 = if is_aromatic_atom(mol, idx) {
        if explicit_plus_rad <= dv {
            dv - explicit_plus_rad
        } else {
            // 芳香原子被假定已处于某个允许的价态,不再补氢
            let satisfied = valens
                .iter()
                .map(|&v| i32::from(v))
                .take_while(|&v| v > 0)
                .any(|v| explicit_plus_rad == v);
            if !satisfied && strict {
                return Err(ValenceError {
                    atom: idx,
                    symbol: element::by_atomic_num(z).map_or("?", |e| e.symbol),
                    valence: ev,
                    kind: ValenceErrorKind::AromaticValenceNotAllowed,
                });
            }
            0
        }
    } else {
        // 非芳香:允许非默认价,取下一个不小于当前价的允许价态
        let found = valens
            .iter()
            .map(|&v| i32::from(v))
            .take_while(|&v| v >= 0)
            .find(|&v| explicit_plus_rad <= v);
        match found {
            Some(v) => v - explicit_plus_rad,
            None => {
                if strict && last_valence(eff_z) != -1 && last_valence(z) > 0 {
                    return Err(ValenceError {
                        atom: idx,
                        symbol: element::by_atomic_num(z).map_or("?", |e| e.symbol),
                        valence: ev,
                        kind: ValenceErrorKind::ExplicitValenceTooHigh,
                    });
                }
                0
            }
        }
    };

    Ok(u8::try_from(res.max(0)).unwrap_or(0))
}

#[cfg(test)]
mod tests {
    use omgkit_io::smiles;

    use super::*;

    fn calc(smi: &str) -> (Vec<i32>, Vec<u8>) {
        let mut m = smiles::parse(smi).unwrap_or_else(|e| panic!("{}", e.render()));
        let r = update_property_cache(&mut m).unwrap_or_else(|e| panic!("{smi}: {e}"));
        (r.explicit_valence, r.implicit_hs)
    }

    #[test]
    fn simple_organics() {
        // (显式价, 隐式氢)
        assert_eq!(calc("CCO"), (vec![1, 2, 1], vec![3, 2, 1]));
        assert_eq!(calc("C"), (vec![0], vec![4]));
        assert_eq!(calc("CC(=O)O"), (vec![1, 4, 2, 1], vec![3, 0, 0, 1]));
    }

    #[test]
    fn bracket_atoms_get_no_implicit_hs() {
        // [CH4] 的氢是作者显式给的,不再推断
        let (ev, ih) = calc("[CH4]");
        assert_eq!(ev, vec![4]);
        assert_eq!(ih, vec![0], "NO_IMPLICIT 置位时隐式氢恒为 0");
    }

    #[test]
    fn aromatic_ring_atoms() {
        // 苯:芳香碳的显式价 3(1.5+1.5),补 1 个氢
        let (ev, ih) = calc("c1ccccc1");
        assert!(ev.iter().all(|&v| v == 3), "实际 {ev:?}");
        assert!(ih.iter().all(|&h| h == 1), "实际 {ih:?}");
    }

    #[test]
    fn kekulized_form_gives_same_counts() {
        // 凯库勒式与芳香式的隐式氢数应一致
        let (_, ih_arom) = calc("c1ccccc1");
        let (_, ih_kek) = calc("C1=CC=CC=C1");
        assert_eq!(ih_arom, ih_kek);
    }

    #[test]
    fn charged_atoms_use_effective_atomic_number() {
        // [NH4+] 的有效原子序数是 7-1=6(碳),故允许 4 价
        let (ev, ih) = calc("[NH4+]");
        assert_eq!(ev, vec![4]);
        assert_eq!(ih, vec![0]);
        // 铵的非方括号写法不存在,用 [O-] 验证负电
        let (ev, _) = calc("CC(=O)[O-]");
        assert_eq!(ev[3], 1);
    }

    #[test]
    fn hypervalent_sulfur() {
        // 硫酸:S 是 6 价
        let (ev, ih) = calc("O=S(=O)(O)O");
        assert_eq!(ev[1], 6, "S 应为 6 价");
        assert_eq!(ih[1], 0);
    }

    #[test]
    fn wildcard_gets_no_implicit_hs() {
        let (_, ih) = calc("*CC");
        assert_eq!(ih[0], 0, "通配原子不补氢");
    }

    #[test]
    fn overvalent_carbon_is_rejected() {
        // 五键碳:严格模式下必须判为失败
        let mut m = smiles::parse("C(C)(C)(C)(C)C").unwrap();
        let err = update_property_cache(&mut m).expect_err("五键碳应当被拒绝");
        assert_eq!(err.atom, 0);
        assert_eq!(err.kind, ValenceErrorKind::ExplicitValenceTooHigh);
        assert_eq!(err.valence, 5);
    }

    #[test]
    fn implicit_hs_written_back_to_atoms() {
        let mut m = smiles::parse("CCO").unwrap();
        let r = update_property_cache(&mut m).unwrap();
        for (i, a) in m.atoms().iter().enumerate() {
            assert_eq!(a.num_implicit_hs, r.implicit_hs[i]);
        }
    }

    #[test]
    fn dative_bond_asymmetry_flows_through() {
        // 配位键:给体不计价,受体计 1
        let mut m = smiles::parse("CC").unwrap();
        m.bond_mut(0)
            .unwrap()
            .set_order(omgkit_core::BondOrder::Dative);
        let r = update_property_cache(&mut m).unwrap();
        assert_eq!(r.explicit_valence[0], 0, "给体(起点)显式价为 0");
        assert_eq!(r.explicit_valence[1], 1, "受体(终点)显式价为 1");
        assert_eq!(r.implicit_hs, vec![4, 3]);
    }
}
