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
    /// 破对称动了几个原子(嵌入给出重合坐标的那些)。
    ///
    /// **这个数不该常年是 0** —— 对称分子本来就会撞上简并。它一直是 0
    /// 反倒说明这一步没接上,或者判据的样本里没有对称分子。
    pub spread: usize,
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
/// # 调用方要先做的两件事
///
/// 1. **补氢**(界矩阵与手性都按显式氢算)。
/// 2. **把 SMILES 的 `/` `\` 折算成双键自己的 `BondStereo`** ——
///    `omgkit_io::stereo::perceive_bond_stereo`。顺反记在**相邻单键**的
///    `direction` 上,不折算的话 `bounds::stereo_path_torsion` 一次都不发力,
///    双键的 1-4 扭转整档退回"顺式到反式的全程",交付的几何会有一半站错边。
///
/// 第 2 条整条流水线先前**压根没做**:实测全语料 405 条双键受影响,
/// 外部判据(RDKit 从坐标读回)上 10 个分子交付的是错的几何。
/// 这不是能靠文档解决的 —— 所以它同时被
/// `omgkit_io::stereo::directions_not_perceived` 这个谓词看着
/// (那个谓词与感知**由构造保证一致**),`examples/feasibility.rs` 拿它当闸。
///
/// `centers` 由调用方给 —— 见 [`crate::chiral::centers`]。
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

    // **破对称必须在优化器之前。** 对称分子的 Gram 矩阵有重特征值,等价原子会拿到
    // 逐位相同的坐标 —— 而完全重合的两个原子**梯度恰好为零**(方向向量是零向量),
    // 优化器永远分不开它们。实测语料里 0.50% 的分子这样,全语料 44 个,
    // 而且是静默的:坐标照样返回,只是废的。见 `crate::spread`。
    let spread = crate::spread::break_coincidence(&mut coords);

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
        spread,
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
        // **这一步先前没有,整条流水线都没有。** 见 `conformer` 的前置条件那一节。
        omgkit_io::stereo::perceive_bond_stereo(&mut m);
        let r = omgkit_io::canon::classed_ranks(&m);
        omgkit_chem::add_explicit_hs(&mut m, &r);
        m
    }

    /// 一条 1-2-3-4 路径在给定坐标下的扭转角(度)。
    fn torsion(p: &[[f64; 3]], a: usize, b: usize, c: usize, d: usize) -> f64 {
        let sub = |u: [f64; 3], v: [f64; 3]| [u[0] - v[0], u[1] - v[1], u[2] - v[2]];
        let dot = |u: [f64; 3], v: [f64; 3]| u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
        let cross = |u: [f64; 3], v: [f64; 3]| {
            [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ]
        };
        let (b0, b1, b2) = (sub(p[a], p[b]), sub(p[c], p[b]), sub(p[d], p[c]));
        let n = dot(b1, b1).sqrt();
        let u = [b1[0] / n, b1[1] / n, b1[2] / n];
        let proj = |v: [f64; 3]| {
            let t = dot(v, u);
            [v[0] - t * u[0], v[1] - t * u[1], v[2] - t * u[2]]
        };
        let (v, w) = (proj(b0), proj(b2));
        // `y.atan2(x)`,**别把两个参数写反** —— 反了得到的是 90° − τ,
        // 而那种错在"离 0 还是离 180"上极难看出来(诊断时踩过)。
        dot(cross(u, v), w).atan2(dot(v, w)).to_degrees()
    }

    #[test]
    fn 双键顺反必须落到正确的一侧() {
        // **期望值写死在这里,不从 `bd.stereo` 读。**
        // 从 `bd.stereo` 读的话,"整批顺反反了"这种变异会让期望和几何一起翻,
        // 永远自洽 —— 实测把 `perceive_bond_stereo` 的 Cis/Trans 对调,
        // 那种写法照样全绿。成对给出顺式与反式也救不了,原因同上。
        //
        // 第二列是**几何事实**:`/F ... /F` 两个 F 在双键两侧(反式,|τ| > 90°)。
        for (smi, cis) in [
            ("F/C=C/F", false),
            ("F/C=C\\F", true),
            ("C/C=C/C", false),
            ("C/C=C\\C", true),
            ("Cl/C(C)=C(C)/Br", false),
            ("Cl/C(C)=C(C)\\Br", true),
            // 环上的那一档 —— 先前交付的几何就是在这里站错边的
            ("[H]/N=C/1\\N[C@]2(CSC(=[NH+]2)N)CS1", true),
            ("[H]/N=C/1\\N=C([C@H](S1)CC(=O)[O-])O", true),
            ("CCOC(=O)[C@@H]1C(=N/C(=N/CC=C)/S1)C", false),
        ] {
            let m = prep(smi);
            // 先确认标记真的被折算出来了 —— 否则这条测试测了个寂寞
            assert!(
                !omgkit_io::stereo::directions_not_perceived(&m),
                "{smi}:有方向键没折算,`prep` 漏了 perceive_bond_stereo"
            );
            let marked: Vec<_> = m
                .bonds()
                .iter()
                .filter(|b| b.stereo != omgkit_core::BondStereo::None)
                .copied()
                .collect();
            assert!(!marked.is_empty(), "{smi}:一根带立体标记的双键都没有");
            let centers = chiral::centers(&m);
            let c = conformer(&m, &centers).unwrap_or_else(|e| panic!("{smi} 失败:{e:?}"));
            // 每个分子的第一根带标记的双键 —— 就是 SMILES 里写明的那一根,
            // 拿**写死的**期望去比
            let bd = marked[0];
            let (i, j) = (bd.stereo_atoms[0] as usize, bd.stereo_atoms[1] as usize);
            let t = torsion(&c.coords, i, bd.begin as usize, bd.end as usize, j);
            assert_eq!(
                t.abs() < 90.0,
                cis,
                "{smi} 键 {}={}({:?}):参照 {i}/{j} 的扭转 {t:.1}°,应当在 {} 一侧",
                bd.begin,
                bd.end,
                bd.stereo,
                if cis {
                    "顺式(|τ|<90°)"
                } else {
                    "反式(|τ|>90°)"
                }
            );
        }
    }

    #[test]
    fn 漏了顺反折算会被谓词看见() {
        // `directions_not_perceived` 是把"前置条件"变成机器可查的那一条。
        // 它一旦恒为 false,`feasibility` 那道闸就什么都没守住。
        let mut m = omgkit_io::smiles::parse("F/C=C/F").expect("解析");
        omgkit_chem::pipeline::sanitize(&mut m).expect("净化");
        // **故意不调** perceive_bond_stereo
        let r = omgkit_io::canon::classed_ranks(&m);
        omgkit_chem::add_explicit_hs(&mut m, &r);
        assert!(
            omgkit_io::stereo::directions_not_perceived(&m),
            "漏了折算却没被谓词看见 —— 那道前置条件闸是瞎的"
        );
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
