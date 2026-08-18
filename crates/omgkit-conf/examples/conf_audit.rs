//! **语料判据**:全语料跑一遍构造法,量的是"能不能用",不是"像不像 MMFF"。
//!
//! # 判什么
//!
//! 按方案 §6 的粗线,只守**硬失败**:
//!
//! | 判据 | 为什么是它 |
//! |---|---|
//! | 覆盖率 | 摆不出来,下游拿到的是个洞 |
//! | 永不 panic | 摆不了要返回"摆不了",不能崩 |
//! | 非键最小距离 | 两个原子叠在一起,`1/r¹²` 直接爆掉,优化器起不来 |
//! | 键长 | 构造法唯一"按定义精确"的量,它不准说明接错了 |
//! | 键角 | 粗线(15°),**不比 MMFF 的分位数 band** |
//!
//! **不判**:是否接近 MMFF 优化后的几何、环的构象族、扭转角的能量 ——
//! 那些是后续力场的事。
//!
//! # 这一版的范围
//!
//! 一期只摆**无环**分子。有环的如实计数、不判 —— 是范围不是失败。

use omgkit_conf::{build, params};

/// 非键最小距离的绝对下限(Å)。
///
/// 方案初稿写的是 0.5 Å,**送审用实测否掉了**:消撞之后最坏的纯有机例子是
/// 0.71 Å,稳稳过闸;而 0.71 Å 的 H···H 在 LJ 里是约 4×10⁶ 倍 ε,
/// 优化器根本起不来 —— 而"起得来"正是这条闸自称的用途。
/// RDKit 在同一批分子上的地板是 1.89 Å、中位 2.07 Å。
const MIN_NONBONDED: f64 = 1.6;

/// 非键距离相对 vdW 之和的下限。与 [`MIN_NONBONDED`] 两条都要过。
const MIN_VDW_FRAC: f64 = 0.75;

/// 键长相对误差上限。构造法这一条**按定义**应当是机器精度。
const MAX_BOND_REL: f64 = 1e-9;

/// 键角**超出解析边界**多少才算红(度)。
///
/// 头一版这里是"与表值偏差 < 15°"一刀切,**口径是错的**:
/// 四配位中心有 6 个夹角、只有 5 个自由度,构造法让"父–子"那几个精确等于表值,
/// **兄弟之间的是推出来的**(`cos φ = cos²θ + sin²θ·cos120°`)。
/// 实测磷那个中心表值 98.20°、实得 118.00°,而公式给的正是 118.00° ——
/// 拿它去比表值等于在判一件本来就不成立的事。
///
/// 所以改成:逐中心算出**解析上应当偏多少**([`omgkit_conf::vsepr::sibling_skew`]),
/// 实测偏差超出它这么多才算红。这样既挡得住真的几何错,又不冤枉超定带来的偏差。
const MAX_ANGLE_EXCESS: f64 = 0.5;

/// 兄弟角被推歪超过这个值的中心,单独计数(度)。**这才是真正的质量信号** ——
/// θ 离 109.47° 越远,兄弟角错得越狠(120° 时差 22.8°)。
const STRAIN_WARN: f64 = 5.0;

/// 因超配位中心(配位数 ≥ 5)而摆不到的原子数上限。贴着现值。
///
/// 这一档是方案 §4.5 写明的**范围外** —— 不进覆盖率的分母,但要有闸:
/// 它一涨就说明有别的东西掉进了这条兜底路,而覆盖率**看不出来**。
const MAX_HYPERVALENT: usize = 285;

fn main() {
    let corpus = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "harness/corpus/large.smi".into());
    let text = std::fs::read_to_string(&corpus).unwrap_or_else(|e| {
        eprintln!("读不了语料 {corpus}:{e}");
        std::process::exit(1);
    });

    let (mut n_line, mut n_parse_bad, mut n_sanitize_bad) = (0u64, 0u64, 0u64);
    let (mut n_cyclic, mut n_acyclic, mut n_full) = (0u64, 0u64, 0u64);
    let (mut worst_bond, mut worst_angle) = (0.0f64, 0.0f64);
    let mut agg = build::Stats::default();
    // 非键最小距离的分布
    let mut mind_all: Vec<f64> = Vec::new();
    let (mut n_below_abs, mut n_below_frac) = (0u64, 0u64);
    let mut worst_case: Option<(f64, String)> = None;
    let mut worst_ang_case: Option<(f64, String, String)> = None;

    for line in text.lines() {
        let smi = line.split('\t').next().unwrap_or("").trim();
        if smi.is_empty() {
            continue;
        }
        n_line += 1;
        let Ok(mut mol) = omgkit_io::smiles::parse(smi) else {
            n_parse_bad += 1;
            continue;
        };
        if omgkit_chem::pipeline::sanitize(&mut mol).is_err() {
            n_sanitize_bad += 1;
            continue;
        }
        let r = omgkit_io::canon::classed_ranks(&mol);
        omgkit_chem::add_explicit_hs(&mut mol, &r);
        let r = omgkit_io::canon::classed_ranks(&mol);

        if !omgkit_chem::sssr::ring_set(&mol).is_empty() {
            n_cyclic += 1;
            continue; // 一期范围之外
        }
        n_acyclic += 1;

        let out = build::place(&mol, &r);
        agg.atoms += out.stats.atoms;
        agg.placed += out.stats.placed;
        agg.degenerate += out.stats.degenerate;
        agg.disconnected += out.stats.disconnected;
        agg.bond_table += out.stats.bond_table;
        agg.bond_relaxed += out.stats.bond_relaxed;
        agg.bond_model += out.stats.bond_model;
        agg.angle_table += out.stats.angle_table;
        agg.angle_relaxed += out.stats.angle_relaxed;
        agg.angle_model += out.stats.angle_model;
        agg.degree_ge5 += out.stats.degree_ge5;
        agg.skipped_hypervalent += out.stats.skipped_hypervalent;
        agg.skipped_ring += out.stats.skipped_ring;
        agg.angle_strained += out.stats.angle_strained;
        if out.complete() {
            n_full += 1;
        }

        // 键长:逐根比表值
        for b in mol.bonds() {
            if !out.placed[b.begin as usize] || !out.placed[b.end as usize] {
                continue;
            }
            let want = params::bond_length(
                mol.atoms()[b.begin as usize].atomic_num,
                mol.atoms()[b.end as usize].atomic_num,
                b.order,
                0,
            )
            .value;
            let got = out.coords[b.begin as usize].dist(out.coords[b.end as usize]);
            worst_bond = worst_bond.max(((got - want) / want).abs());
        }

        // 键角:逐个比表值(只比父–子那种,兄弟角是推出来的,见 Stats::angle_strained)
        for k in 0..mol.num_atoms() {
            let Ok(ku) = u32::try_from(k) else { continue };
            let nb: Vec<u32> = mol.neighbors(ku).map(|(y, _)| y).collect();
            if nb.len() < 2 || !out.placed[k] {
                continue;
            }
            let want = params::angle(
                mol.atoms()[k].atomic_num,
                nb.len(),
                mol.atoms()[k]
                    .flags
                    .contains(omgkit_core::AtomFlags::AROMATIC),
                0,
                0,
            )
            .value
            .to_degrees();
            for i in 0..nb.len() {
                for j in (i + 1)..nb.len() {
                    if !out.placed[nb[i] as usize] || !out.placed[nb[j] as usize] {
                        continue;
                    }
                    if let Some(a) = omgkit_conf::geom::angle_at(
                        out.coords[nb[i] as usize],
                        out.coords[k],
                        out.coords[nb[j] as usize],
                    ) {
                        // 解析上这个中心的兄弟角**应当**偏多少
                        let arr =
                            omgkit_conf::vsepr::arrangement(mol.atoms()[k].hybridization, nb.len());
                        let bound =
                            omgkit_conf::vsepr::expected_sibling_skew(arr, want.to_radians())
                                .to_degrees();
                        let dev = ((a.to_degrees() - want).abs() - bound).max(0.0);
                        if dev > worst_angle {
                            worst_angle = dev;
                            worst_ang_case = Some((
                                dev,
                                smi.to_string(),
                                format!(
                                    "中心元素 {} 配位 {} 表值 {want:.2}° 实得 {:.2}° 解析边界 {bound:.2}°(邻居 {} / {})",
                                    mol.atoms()[k].atomic_num,
                                    nb.len(),
                                    a.to_degrees(),
                                    nb[i],
                                    nb[j]
                                ),
                            ));
                        }
                    }
                }
            }
        }

        // 非键最小距离(拓扑距离 ≥ 3 才算,1-2 与 1-3 是键长键角管的)
        let n = mol.num_atoms();
        let mut far = vec![vec![false; n]; n];
        for (i, row) in far.iter_mut().enumerate() {
            let Ok(iu) = u32::try_from(i) else { continue };
            for (y, _) in mol.neighbors(iu) {
                row[y as usize] = true;
                for (z, _) in mol.neighbors(y) {
                    row[z as usize] = true;
                }
            }
            row[i] = true;
        }
        let mut mind = f64::INFINITY;
        for (i, row) in far.iter().enumerate() {
            for (j, &near) in row.iter().enumerate().skip(i + 1) {
                if near || !out.placed[i] || !out.placed[j] {
                    continue;
                }
                let d = out.coords[i].dist(out.coords[j]);
                let want = MIN_VDW_FRAC
                    * (params::vdw_radius(mol.atoms()[i].atomic_num)
                        + params::vdw_radius(mol.atoms()[j].atomic_num));
                if d < want {
                    n_below_frac += 1;
                }
                mind = mind.min(d);
            }
        }
        if mind.is_finite() {
            mind_all.push(mind);
            if mind < MIN_NONBONDED {
                n_below_abs += 1;
                if worst_case.as_ref().map_or(true, |(w, _)| mind < *w) {
                    worst_case = Some((mind, smi.to_string()));
                }
            }
        }
    }

    println!("语料 {n_line} 行:parse 失败 {n_parse_bad}、sanitize 失败 {n_sanitize_bad}");
    println!("  有环(一期范围外,不判):{n_cyclic}");
    println!("  无环:{n_acyclic},其中全部摆好 {n_full}");
    println!(
        "  原子 {} 个,摆好 {}(退化 {}、连不上 {})",
        agg.atoms, agg.placed, agg.degenerate, agg.disconnected
    );
    println!(
        "  参数来源 —— 键长:表 {} / 放宽 {} / 模型 {};键角:表 {} / 放宽 {} / 模型 {}",
        agg.bond_table,
        agg.bond_relaxed,
        agg.bond_model,
        agg.angle_table,
        agg.angle_relaxed,
        agg.angle_model
    );
    println!(
        "  配位数 ≥5 的中心 {};兄弟角被推歪 >{STRAIN_WARN}° 的中心 {}(这是超定的必然,不是 bug)",
        agg.degree_ge5, agg.angle_strained
    );
    println!("  键长最大相对误差 {worst_bond:.3e}(上限 {MAX_BOND_REL:.0e})");
    println!("  键角**超出解析边界**最多 {worst_angle:.3}°(上限 {MAX_ANGLE_EXCESS}°)");

    mind_all.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    if !mind_all.is_empty() {
        let q = |f: f64| {
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let i = ((mind_all.len() as f64 - 1.0) * f).round() as usize;
            mind_all[i]
        };
        println!(
            "  非键最小距离:地板 {:.3} / p05 {:.3} / 中位 {:.3} Å",
            mind_all[0],
            q(0.05),
            q(0.5)
        );
        println!(
            "    低于 {MIN_NONBONDED} Å 的分子 {n_below_abs} 个({:.1}%);低于 {MIN_VDW_FRAC}×vdW 的原子对 {n_below_frac} 个",
            100.0 * n_below_abs as f64 / mind_all.len() as f64
        );
    }
    if let Some((d, smi, what)) = &worst_ang_case {
        println!("    角最狠的:{d:.2}°  {what}  {smi}");
    }
    if let Some((d, smi)) = &worst_case {
        println!("    最狠的:{d:.3} Å  {smi}");
    }

    let mut fatal = false;
    // **范围内**的原子必须 100% 摆好。超配位那一档是方案 §4.5 写明的范围外
    // (`vsepr` 对它只有"均分了事"的 Spread,实测钴的六配位差 53°),
    // 所以它不进分母 —— 但它**单独有闸**,不许悄悄变多。
    let in_scope = agg
        .atoms
        .saturating_sub(agg.skipped_ring + agg.skipped_hypervalent);
    if agg.placed != in_scope {
        eprintln!(
            "\n范围内的原子没摆全:{} / {in_scope} —— 一期这一条必须 100%",
            agg.placed
        );
        fatal = true;
    }
    if agg.degenerate != 0 {
        eprintln!(
            "\n有 {} 个原子因参考点共线摆不出来 —— 垂直参考点那条兜底该接住它们",
            agg.degenerate
        );
        fatal = true;
    }
    if agg.skipped_hypervalent > MAX_HYPERVALENT {
        eprintln!(
            "\n因超配位中心摆不到的原子 {} 个,超过上限 {MAX_HYPERVALENT} —— 这一档是范围外,但不许悄悄变多",
            agg.skipped_hypervalent
        );
        fatal = true;
    }
    // 守恒:摆好的 + 环上的 + 超配位挡住的 + 退化的 = 总数
    if agg.placed + agg.skipped_ring + agg.skipped_hypervalent + agg.degenerate != agg.atoms {
        eprintln!("\n原子的账对不上 —— 有没计数的分支");
        fatal = true;
    }
    if worst_bond > MAX_BOND_REL {
        eprintln!(
            "\n键长最大相对误差 {worst_bond:.3e} 超过 {MAX_BOND_REL:.0e} —— 构造法这条该是机器精度"
        );
        fatal = true;
    }
    if worst_angle > MAX_ANGLE_EXCESS {
        eprintln!(
            "\n键角超出解析边界 {worst_angle:.3}°,过了上限 {MAX_ANGLE_EXCESS}° —— \
             这不是超定带来的偏差,是真的摆错了"
        );
        fatal = true;
    }
    if agg.disconnected != 0 {
        eprintln!(
            "\n有 {} 个原子既没摆好、又说不出原因 —— 分量是预先算好的,这条该恒为 0",
            agg.disconnected
        );
        fatal = true;
    }
    if fatal {
        std::process::exit(1);
    }
    println!("\n(非键距离这一条本期只报不闸 —— 消撞还没做,见方案 §4.6)");
}
