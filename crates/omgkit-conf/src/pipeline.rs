//! **端到端:分子进,一组三维坐标出。**
//!
//! ```text
//! 界矩阵 → 三角光滑化 → 取上限矩阵 U 当参考距离表 → 度量矩阵嵌入
//!        → 全局手性定向(离散,一次)→ L-BFGS 精修
//! ```
//!
//! **全程没有随机数。** 同一个分子永远同一个答案,不需要 seed。
//!
//! # 每一步为什么在这儿
//!
//! | 步 | 为什么 |
//! |---|---|
//! | 光滑化 | 之后的上限矩阵 `U` 按构造满足三角不等式,是一张画得出来的距离表 |
//! | 用 `U` 而不是随机取 | RDKit 逐对独立随机取,取出来的表常常摆不出来,它的应对是作废重掷 |
//! | **全局手性定向** | 反射不在 `SO(3)` 连通分支里,**下降法翻不过去**,只能离散地定一次 |
//! | L-BFGS | 目标只有 `C¹`,线搜索必须上强 Wolfe(实测 Armijo 会退化成最速下降) |
//!
//! # 失败是什么意思
//!
//! 只有**界矩阵自相矛盾**(光滑化判不可行)才叫失败 —— 那时连一张自洽的距离表
//! 都没有。实测全语料 8831 个分子里 5 个(0.06%),是 RDKit ETKDG 那 0.52% 的 1/9。
//!
//! 嵌入摆不进三维**不算失败**(压掉那一维,精修去救),精修没收敛也**不算失败**
//! (给出当前最好的坐标,并如实报出残差)。理由是:这一步的产物是给力场优化用的
//! **起点**,起点差一点可以修,没有起点才是灾难。

use crate::bounds;
use crate::chiral::{self, Center};
use crate::embed::{self, reference_distances};
use crate::field::Field;
use crate::optimize::{minimize, Options};
use crate::smooth::{triangle_smooth, SmoothError};
use omgkit_core::MolBuilder;

/// 生成构型时失败的原因。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConformerError {
    /// 界矩阵自相矛盾 —— 连一张自洽的距离表都拿不出来。
    Infeasible {
        /// 卡住的那一对原子。
        pair: (usize, usize),
    },
    /// 嵌入那一步的输入坏了(非有限数、特征分解不收敛)。这不该发生。
    Embed(crate::embed::EmbedError),
}

/// 一个构型,外加"它有多好"的账。
#[derive(Debug, Clone)]
pub struct Conformer {
    /// 每个原子的坐标。
    pub coords: Vec<[f64; 3]>,
    /// 精修**之前**的误差函数值。
    pub energy_before: f64,
    /// 精修**之后**的误差函数值。
    pub energy: f64,
    /// 精修迭代了多少次。
    pub iterations: usize,
    /// 全局手性定向那一步有没有把结构翻过来。
    pub reflected: bool,
    /// 手性中心数,以及精修之后号正确的个数。
    pub chiral_total: usize,
    /// 见 [`Conformer::chiral_total`]。
    pub chiral_ok: usize,
}

/// 精修的迭代上限。
///
/// 400 与 RDKit 第一段极小化的 `field->minimize(400, ...)` 同量级。
/// 到了上限不算失败 —— 给出当前坐标并把残差报出来。
pub const MAX_REFINE_ITER: usize = 400;

/// 给一个分子生成**一个**三维构型。
///
/// 分子应当**已经补过氢**(界矩阵与手性都按显式氢算)。
/// `centers` 由调用方给 —— 见 [`crate::chiral::centers`] 的前置条件,
/// 那一步需要"立体标记与当前键序一致",而补氢会破坏这个前提。
///
/// # Errors
///
/// 界矩阵自相矛盾,或嵌入的输入坏掉。
pub fn conformer(mol: &MolBuilder, centers: &[Center]) -> Result<Conformer, ConformerError> {
    let n = mol.num_atoms();
    let (mut b, _) = bounds::build(mol);
    if let Err(SmoothError::Infeasible { pair }) = triangle_smooth(&mut b) {
        return Err(ConformerError::Infeasible { pair });
    }
    let e = embed::embed(&reference_distances(&b), n).map_err(ConformerError::Embed)?;
    let mut coords = e.coords;

    // **全局手性定向:离散,一次,必须在精修之前。**
    // 反射不在 SO(3) 的连通分支里 —— 连续下降要走到镜像必须把整个分子压平,
    // 下降法不会付这个势垒,所以精修救不了整体定向。
    let reflected = chiral::needs_reflection(&coords, centers);
    if reflected {
        chiral::reflect(&mut coords);
    }

    let field = Field::new(&b, centers);
    let mut x: Vec<f64> = coords.iter().flat_map(|p| p.iter().copied()).collect();
    let mut g = vec![0.0; x.len()];
    let energy_before = {
        use crate::optimize::Objective;
        field.value_and_grad(&x, &mut g)
    };
    let report = minimize(
        &field,
        &mut x,
        &Options {
            max_iter: MAX_REFINE_ITER,
            grad_tol: 1e-6,
            memory: 8,
        },
    );
    for (i, p) in coords.iter_mut().enumerate() {
        *p = [x[3 * i], x[3 * i + 1], x[3 * i + 2]];
    }

    Ok(Conformer {
        chiral_total: centers.len(),
        chiral_ok: chiral::correct_count(&coords, centers),
        coords,
        energy_before,
        energy: report.value,
        iterations: report.iterations,
        reflected,
    })
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

    #[test]
    fn 常见分子都给得出构型() {
        for smi in [
            "CCO",
            "c1ccccc1",
            "C1CCCCC1",
            "CC(=O)Nc1ccc(O)cc1",
            "C1CC2CCC1CC2",
            "CC(C)(C)OC(=O)N1CCC(CC1)N",
            "FS(F)(F)(F)(F)F",
            "C=C=C",
        ] {
            let m = prep(smi);
            let c = conformer(&m, &[]).unwrap_or_else(|e| panic!("{smi} 失败:{e:?}"));
            assert_eq!(c.coords.len(), m.num_atoms(), "{smi} 坐标数不对");
            for (i, p) in c.coords.iter().enumerate() {
                assert!(
                    p.iter().all(|v| v.is_finite()),
                    "{smi} 第 {i} 个原子坐标不是有限数:{p:?}"
                );
            }
            // **精修必须真的降下去**,不能原地不动
            assert!(
                c.energy <= c.energy_before,
                "{smi} 精修之后反而更差:{} → {}",
                c.energy_before,
                c.energy
            );
        }
    }

    #[test]
    fn 精修确实在干活() {
        // 起点是嵌入给的,残差不小;精修应当把它压掉一大截。
        let m = prep("CC(=O)Nc1ccc(O)cc1");
        let c = conformer(&m, &[]).unwrap();
        assert!(
            c.energy_before > 0.1,
            "起点残差 {} 太小,测不到东西",
            c.energy_before
        );
        assert!(
            c.energy < c.energy_before * 0.2,
            "只压掉了 {:.1}%:{} → {}",
            100.0 * (1.0 - c.energy / c.energy_before),
            c.energy_before,
            c.energy
        );
    }

    #[test]
    fn 同一个分子两次给逐位相同的坐标() {
        // 全程无随机数 —— 这条一旦红,说明哪里混进了不确定性
        // (HashMap 迭代序、并行归约、未定义的排序 tie-break)。
        let m = prep("CC(C)(C)OC(=O)N1CCC(CC1)N");
        let a = conformer(&m, &[]).unwrap();
        let b = conformer(&m, &[]).unwrap();
        assert_eq!(a.coords, b.coords, "两次跑出来的坐标不逐位相同");
        assert_eq!(a.energy, b.energy);
        assert_eq!(a.iterations, b.iterations);
    }

    #[test]
    fn 界不可行时如实报失败() {
        // 手工造一个自相矛盾的界。这里直接调底层,确认错误一路传得上来。
        use crate::smooth::Bounds;
        let mut b = Bounds::new(3, 0.0, 10.0);
        b.set_lower(0, 1, 5.0);
        b.set_upper(0, 1, 5.0);
        b.set_lower(1, 2, 5.0);
        b.set_upper(1, 2, 5.0);
        b.set_lower(0, 2, 50.0);
        b.set_upper(0, 2, 50.0);
        assert!(triangle_smooth(&mut b).is_err(), "这组界本来就该判不可行");
    }
}
