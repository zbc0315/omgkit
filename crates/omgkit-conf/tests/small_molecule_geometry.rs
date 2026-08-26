//! **小分子的几何**:被闸盯着的两份语料全是药物样大分子,这一档从来没人看。
//!
//! 是这么发现的:`feasibility` 的四档越界改成逐语料棘轮之后,顺手拿粗闸跑了
//! 一遍解析冒烟语料 `harness/corpus/smoke.smi`(里面才有小分子),当场红 ——
//! 1-2 键越界 0.932%,而 `large` 上只有 0.026%。追下去是**甲醇**:
//!
//! ```text
//! CO 的 C–O 出来 1.211 Å,而界是 [1.374, 1.394] —— 短了 0.163 Å
//! ```
//!
//! # 根因不是界,也不是线搜索
//!
//! 界是对的(手算能摆出满足全部 15 对的构型),`CO.CO`(同一个分子凑成两个
//! 片段)也完全正常。精修**收敛了**(梯度 4.4e-7),从产物原地再起一次走 0 步
//! —— 它停在一个**真的局部极小**上,残差 4.04e-2。
//!
//! 换参考距离表(在上下界之间插值)没用:键的界宽只有 0.02 Å,插值改不动形状,
//! 六种插值落到同一个极小、三位有效数字一致。而把**起点扰动**一下,
//! 12 次里 11 次收到 ~0。
//!
//! 修法是 `pipeline` 里的确定性重试阶梯。这个文件是它的判据:**小分子的四档
//! 越界必须是 0**。少了它,下一次动优化器就没人拦得住这一档退回去。

use omgkit_conf::{bounds, pipeline, smooth};

/// 拿这些分子跑一遍,每一根键、每一对都必须落在界内(容差 0.1 Å)。
///
/// 名单里前四个是当初报出问题的那几个(甲醇、甲硫醇、叠氮甲烷、高氯酸),
/// 后面几个是同量级的常见小分子 —— 一个都不许越界。
#[test]
fn small_molecules_land_inside_their_bounds() {
    const TOL: f64 = 0.1;
    let mut bad = Vec::new();
    for smi in [
        "CO",
        "CS",
        "CN=N#N",
        "OCl(=O)(=O)=O",
        "OBr(=O)(=O)=O",
        "CN",
        "CF",
        "CC",
        "CCO",
        "C=O",
        "OO",
        "C#N",
        "CCl",
        "CBr",
        "C=C",
        "C#C",
        "NO",
        "SO",
        "OS(=O)(=O)O",
        "CC(=O)O",
    ] {
        let mut mol = omgkit_io::smiles::parse(smi).unwrap_or_else(|e| panic!("{smi}: {e:?}"));
        let conf = pipeline::conformer_for(&mut mol)
            .unwrap_or_else(|e| panic!("{smi}:生成不出构型 {e:?}"));

        let (mut b, _) = bounds::build(&mol);
        let _ = smooth::triangle_smooth(&mut b);
        let n = mol.num_atoms();
        for i in 0..n {
            for j in (i + 1)..n {
                let d = (0..3)
                    .map(|k| (conf.coords[i][k] - conf.coords[j][k]).powi(2))
                    .sum::<f64>()
                    .sqrt();
                let over = (b.lower(i, j) - d).max(d - b.upper(i, j));
                if over > TOL {
                    bad.push(format!(
                        "{smi}:{i}-{j} 实测 {d:.3},界 [{:.3},{:.3}],超 {over:.3} Å",
                        b.lower(i, j),
                        b.upper(i, j)
                    ));
                }
            }
        }
    }
    assert!(
        bad.is_empty(),
        "小分子越界 {} 处:\n  {}",
        bad.len(),
        bad.join("\n  ")
    );
}
