//! 净化管线的唯一入口。
//!
//! # 为什么要有这个模块
//!
//! 步骤之间有顺序依赖,而且**不是所有依赖都写在类型里**。手抄一遍调用序列
//! 很容易漏掉某一步或调错顺序,而漏掉的后果往往不是报错,是某个字段悄悄
//! 停在中间状态 —— 例如少了第 12 步,吡咯氮的总氢数会从 1 变成 0,而比隐式氢
//! 的差分测试照样全绿(那一列恰好两边都是 0)。
//!
//! 所以调用序列只写一遍,放在这里。需要"只跑到某一步"的测试自己拼,
//! 但那是刻意为之的例外,不是默认做法。
//!
//! # 顺序里两处不能动的地方
//!
//! **芳香化之后要重算一次价键。** 芳香化把键改回芳香键,键级和随之改变,
//! 隐式氢要跟着更新。
//!
//! **第 12 步要用重算之前的隐式氢。** 它的判据是"重算后少了多少",所以那份
//! 旧值必须在重算前抓下来 —— 见 [`mod@crate::adjust_hs`] 的模块文档。

use omgkit_core::MolBuilder;

use crate::{
    adjust_hs, assign_radicals, clean_up, cleanup_chirality, cleanup_organometallics, kekulize,
    perceive_rings, set_aromaticity, set_conjugation, set_hybridization, update_property_cache,
    KekulizeError, ValenceError,
};

/// 净化失败的原因。
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SanitizeError {
    /// 价键校验不通过(超价等)
    Valence(ValenceError),
    /// 芳香环无法写成交替单双键
    Kekulize(KekulizeError),
}

impl std::fmt::Display for SanitizeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Valence(e) => write!(f, "价键计算失败:{e}"),
            Self::Kekulize(e) => write!(f, "kekulize 失败:{e}"),
        }
    }
}

impl std::error::Error for SanitizeError {}

/// 跑完整条净化管线。
///
/// 当前实现的步骤:1、2、3、4、5、6、7、8、9、11、12。只有第 10 步
/// (阻转异构)没做 —— 实测它在 8839 条语料上改动 0 条分子,没有任何用例
/// 能守着它,所以在补到能触发的语料之前不实现。
///
/// 注意第 11 步只剔除**几何上不成立**的立体标记;"剔除非真手性中心"要跑
/// 对称性判定,由 `omgkit_io::stereo::genuine_tetrahedral` 负责,不在净化里。
///
/// # Errors
/// 价键校验或 kekulize 失败时返回 [`SanitizeError`]。失败时分子可能已被
/// 部分修改 —— 需要"要么全成功要么不动"的调用方应当自己先克隆。
pub fn sanitize(mol: &mut MolBuilder) -> Result<(), SanitizeError> {
    clean_up(mol);
    // 第 2 步必须排在价键计算之前 —— 它的作用正是让那些超价原子不再超价,
    // 放在后面的话价键计算已经先一步拒绝了整个分子
    cleanup_organometallics(mol);
    update_property_cache(mol).map_err(SanitizeError::Valence)?;
    let _ = perceive_rings(mol);
    kekulize(mol).map_err(SanitizeError::Kekulize)?;
    assign_radicals(mol);
    set_aromaticity(mol);

    // 第 12 步要用**芳香化之前**算出的隐式氢,所以在重算之前抓一份快照。
    // 顺序反过来的话第 12 步会变成拿新值比新值,恒等,静默失效。
    let implicit_before: Vec<u8> = mol.atoms().iter().map(|a| a.num_implicit_hs).collect();

    update_property_cache(mol).map_err(SanitizeError::Valence)?;
    set_conjugation(mol);
    set_hybridization(mol);
    // 第 11 步用到杂化,必须排在第 9 步之后
    cleanup_chirality(mol);
    adjust_hs(mol, &implicit_before);
    Ok(())
}

#[cfg(test)]
mod tests {
    use omgkit_io::smiles;

    use super::*;

    fn total_hs(m: &MolBuilder, i: usize) -> u32 {
        let a = m.atoms()[i];
        u32::from(a.num_explicit_hs) + u32::from(a.num_implicit_hs)
    }

    /// 整条管线跑完,吡咯氮的氢要还在。
    ///
    /// 这正是第 12 步缺失时丢掉的那个氢。
    #[test]
    fn pyrrole_keeps_its_nh() {
        let mut m = smiles::parse("c1cc[nH]c1").unwrap();
        sanitize(&mut m).expect("应能净化");
        assert_eq!(total_hs(&m, 3), 1, "吡咯氮上的氢");
        assert_eq!(m.atoms()[3].num_explicit_hs, 1, "且记在显式那一侧");
    }

    /// 常见分子跑得通,且总氢数对。
    #[test]
    fn common_molecules() {
        for (smi, atom, hs) in [
            ("CCO", 0, 3),
            ("CCO", 2, 1),
            ("c1ccccc1", 0, 1),
            ("c1ccc2[nH]ccc2c1", 4, 1),
            ("CC(=O)O", 3, 1),
        ] {
            let mut m = smiles::parse(smi).unwrap();
            sanitize(&mut m).unwrap_or_else(|e| panic!("{smi}: {e}"));
            assert_eq!(total_hs(&m, atom), hs, "{smi} 的原子 {atom}");
        }
    }

    /// 超价的分子要报错,而不是给出一个悄悄错着的结果。
    #[test]
    fn hypervalent_is_rejected() {
        let mut m = smiles::parse("C(C)(C)(C)(C)C").unwrap();
        assert!(matches!(sanitize(&mut m), Err(SanitizeError::Valence(_))));
    }

    /// 幂等:整条管线跑两遍,结果不变。
    ///
    /// 第 12 步会改写自己下次要读的字段,是这里最容易出问题的一步。
    #[test]
    fn is_idempotent() {
        for smi in [
            "c1cc[nH]c1",
            "c1ccccc1",
            "CCO",
            "c1ccc2[nH]ccc2c1",
            "CC(=O)O",
        ] {
            let mut m = smiles::parse(smi).unwrap();
            sanitize(&mut m).unwrap();
            let once: Vec<_> = m.atoms().to_vec();
            let once_b: Vec<_> = m.bonds().to_vec();
            sanitize(&mut m).unwrap();
            assert_eq!(m.atoms(), &once[..], "{smi}:原子不幂等");
            assert_eq!(m.bonds(), &once_b[..], "{smi}:键不幂等");
        }
    }
}
