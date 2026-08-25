//! 从 JSONL 基准的连接表重建分子 —— **四个判官共用这一份**。
//!
//! # 为什么必须只有一份
//!
//! 先前四个判官各抄一遍,于是分岔了,而且是静默分岔:
//!
//! | 判官 | 读几列键 | 跑 `sanitize` |
//! |---|---|---|
//! | `bounds_oracle` | 6(用了 `stereo`/`stereo_atoms`) | 是 |
//! | `threading_oracle` | **3**(基准明明有 6 列) | **否** |
//! | `chiral_oracle` | 3 | **否** |
//! | `conformer_oracle` | 3 | **否** |
//!
//! 两处都咬人:
//!
//! **一、漏掉 `sanitize`。** `bounds::build` 读 `AtomFlags::AROMATIC` 与
//! `hybridization`,这两样只有净化才会填 —— 不跑净化,建出来的界是另一张表,
//! 判官量的就不是产品那条路(`feasibility` 走的是
//! `sanitize → perceive_bond_stereo → add_explicit_hs → conformer`)。
//! 实测 `conformer_oracle` 补上净化之后:嵌入阶段的键交叉 **405 → 280**、
//! 精修步数 168 → 179、手性自洽口径 83.8% → 85.0%。
//!
//! **二、丢掉顺反。** `smoke.bounds.jsonl` 的键元组是 6 列,第 4–6 列是
//! `stereo` 与两个参照原子(`dump_bounds.py` 里那段注释专门解释了为什么要导:
//! 界矩阵靠它把 1-4 的顺反析取解掉)。`threading_oracle` 读同一个文件却只取前
//! 3 列。
//!
//! # 已知限制:`smoke.chirality.jsonl` / `smoke.lonepair.jsonl` 只有 3 列
//!
//! `dump_chirality.py` 压根不导 `stereo`。实测(RDKit 2025.09.2)那 150 个分子
//! 里 **23 个带顺反、共 28 根双键**;lonepair 那 15 个一根都没有。
//! 也就是说端到端那两条判官在 15% 的分子上是**看不见顺反的**。
//!
//! 补它要重导基准,而那份文件是 **2022.09.5** 导的(`dump_chirality.py` 的
//! 文件头写着为什么故意与 2025.09.2 的界基准分开),重导会连坐标带真值一起换,
//! 得单独一块来做。这里先把它**报出来**:读到 3 列的分子会计数并打印。

// 四个判官各自 #[path] 引一份,某个判官用不到全部函数是正常的
#![allow(dead_code)]

use omgkit_core::{BondOrder, BondStereo, MolBuilder};

/// 基准里的一根键。`stereo` 为 `None` 表示这份基准没导这一列。
pub struct BondRec {
    pub i: u32,
    pub j: u32,
    pub order: u8,
    /// (RDKit 的 `BondStereo` 号, 参照原子 0, 参照原子 1)
    pub stereo: Option<(i64, i64, i64)>,
}

/// 重建的结果。**失败要分类计数**,不能一律 `continue` —— 那会让分母悄悄变小。
pub enum BuildFail {
    /// 连接表本身建不起来(下标越界、价键冲突)
    Topology,
    /// 净化失败
    Sanitize,
}

/// 读 `v["bonds"]`。第 4–6 列缺席时 `stereo` 是 `None`。
pub fn parse_bonds(v: &serde_json::Value) -> Vec<BondRec> {
    v["bonds"]
        .as_array()
        .into_iter()
        .flatten()
        .filter_map(|e| {
            let t = e.as_array()?;
            let i = u32::try_from(t.first()?.as_u64()?).ok()?;
            let j = u32::try_from(t.get(1)?.as_u64()?).ok()?;
            let order = u8::try_from(t.get(2)?.as_u64()?).ok()?;
            let stereo = match (t.get(3), t.get(4), t.get(5)) {
                (Some(a), Some(b), Some(c)) => Some((
                    a.as_i64()?,
                    b.as_i64().unwrap_or(-1),
                    c.as_i64().unwrap_or(-1),
                )),
                _ => None,
            };
            Some(BondRec {
                i,
                j,
                order,
                stereo,
            })
        })
        .collect()
}

/// 这份连接表带顺反列吗。
pub fn has_stereo_column(bonds: &[BondRec]) -> bool {
    bonds.iter().any(|b| b.stereo.is_some())
}

/// 重建出来的分子,外加**真的写回了几根顺反**。
///
/// 那个计数不是诊断,是闸:实测把"写回顺反"整段删掉,四个判官**一个都不红**
/// —— 而 `dump_bounds.py` 导这一列的唯一理由就是让界矩阵解掉 1-4 的顺反析取。
/// 界那一层自己有判据(`bounds.rs` 的 `cis_and_trans_get_different_bounds`),
/// 缺的正是"判官到底有没有把这一列喂进去",所以数出来、上闸。
pub struct Built {
    pub mol: MolBuilder,
    pub stereo_applied: usize,
}

/// 按**产品那条路**重建:建原子与键 → `sanitize` → 写回顺反。
///
/// 顺反必须写在净化**之后**(净化可能重排键),而且要在建界**之前**。
pub fn build(z: &[u8], chg: &[i8], rad: &[u8], bonds: &[BondRec]) -> Result<Built, BuildFail> {
    let mut m = MolBuilder::new();
    for (k, &a) in z.iter().enumerate() {
        let mut ad = omgkit_core::AtomData::new(a);
        ad.formal_charge = chg.get(k).copied().unwrap_or(0);
        ad.num_radical_electrons = rad.get(k).copied().unwrap_or(0);
        m.add_atom_data(ad);
    }
    for b in bonds {
        let ord = match b.order {
            2 => BondOrder::Double,
            3 => BondOrder::Triple,
            4 => BondOrder::Aromatic,
            _ => BondOrder::Single,
        };
        m.add_bond(b.i, b.j, ord).map_err(|_| BuildFail::Topology)?;
    }
    omgkit_chem::pipeline::sanitize(&mut m).map_err(|_| BuildFail::Sanitize)?;
    let mut stereo_applied = 0usize;
    for (bi, b) in bonds.iter().enumerate() {
        let Some((st, sa0, sa1)) = b.stereo else {
            continue;
        };
        if sa0 < 0 || sa1 < 0 {
            continue;
        }
        // RDKit 的 `Bond::BondStereo`:0 无 2 Z 3 E 4 cis 5 trans
        let s = match st {
            2 => BondStereo::Z,
            3 => BondStereo::E,
            4 => BondStereo::Cis,
            5 => BondStereo::Trans,
            _ => continue,
        };
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        if let Some(mut bond) = m.bond_mut(bi as u32) {
            bond.set_stereo(s);
            bond.set_stereo_atoms([sa0 as u32, sa1 as u32]);
            stereo_applied += 1;
        }
    }
    Ok(Built {
        mol: m,
        stereo_applied,
    })
}
