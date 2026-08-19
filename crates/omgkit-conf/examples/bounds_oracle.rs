//! **外部判官:界矩阵对不对。两条闸,缺一不可。**
//!
//! | 闸 | 问什么 | 不设它会怎样 |
//! |---|---|---|
//! | **正确性** | 真实构象(MMFF 优化过的)是不是落在界内 | — |
//! | **信息量** | 界的宽度是不是不比 RDKit 更松 | 只设第一条的话,**把界写成 `[0, ∞)` 就能满分** |
//!
//! 第二条是必须的:第一条是**单向**的 —— 界越宽越容易过。单向的闸必须配一道
//! 反向的上限,不然它奖励的正是"什么都不约束"。
//!
//! # 原子顺序为什么要靠导出的连接表
//!
//! RDKit 的 `AddHs` 顺序与 omgkit 补氢的顺序不保证一致,两边各自解析 SMILES
//! 会错位 —— 而错位之后判据仍然"跑得通",只是量的是别的分子。所以这里
//! **按 `harness/dump_bounds.py` 导出的 `z` + `bonds` 直接建分子**,下标天生对齐。
//!
//! ```shell
//! python3 harness/dump_bounds.py harness/corpus/large.smi harness/baseline/rdkit_bounds.jsonl 400
//! cargo run -p omgkit-conf --release --example bounds_oracle -- harness/baseline/rdkit_bounds.jsonl
//! ```

use omgkit_conf::bounds;
use omgkit_conf::embed::{embed, metric_matrix, reference_distances};
use omgkit_conf::smooth::{triangle_smooth, Bounds};
use omgkit_core::{BondOrder, MolBuilder};

/// 真实构象允许越界多少(Å)。**不是容差,是给参数表分位数留的余量** ——
/// 表用的是 p05/p95,按定义就有约 10% 的真实值落在外面,越界量应当很小。
const MAX_VIOLATION: f64 = 0.35;

/// 越界超过 [`MAX_VIOLATION`] 的原子对占比上限。
const MAX_VIOLATION_FRAC: f64 = 0.02;

/// 按连接表建不出来的分子数上限。**样本被腰斩不许静悄悄。**
const MAX_BUILD_FAIL: u64 = 8;

/// 界矩阵自相矛盾(光滑化判不可行)的分子数上限。
///
/// **这一项先前只报不闸,是个洞**:界可以越写越不自洽,而判据只会看到
/// 参与统计的对数在掉 —— 而对数掉了之后,剩下样本上的比值反而更好看。
/// 实测加 1-5 链式约束时它从 3 涨到 11,正是这么暴露的。
const MAX_INFEASIBLE: u64 = 5;

/// 我们的界宽相对 RDKit 的中位比值上限。
///
/// # 这是**棘轮**,不是"达标"
///
/// 立这条闸时写的是 1.0(与 RDKit 持平),那是目标不是现状。一路收下来:
/// 1.599 → 1.313 → 1.269 → 1.162 → 1.043 → **1.020**,四档里三档已经逐位相同:
///
/// | 档 | 我们 | RDKit | 比 |
/// |---|---|---|---|
/// | 1-2 | 0.020 | 0.020 | 1.00 |
/// | 1-3 | 0.080 | 0.080 | 1.00 |
/// | 1-4 | 0.120 | 0.120 | 1.00 |
/// | ≥1-5 | 3.595 | 3.400 | **1.06** |
///
/// 剩下这 2% 的来源**已经定位**:RDKit 有一层 1-5 的**链式约束**
/// (`BoundsMatrixBuilder.cpp:1997-2045`)—— 一条 5 原子路径上两个扭转都被钉住时
/// (cis/cis、cis/trans、trans/trans),1-5 距离可以直接算出来再 `± DIST15_TOL`。
/// 我们还没有这一层,那一档的上界全靠三角光滑化推。
///
/// 闸设在 1.05 是**贴着现值的棘轮**:它拦得住回退,而不是把目标改成现状 ——
/// 一条永远红的闸,所有人都会学会忽略它。**1-5 链式约束落地后这个数要跟着降。**
const MAX_WIDTH_RATIO: f64 = 1.05;

/// 判据三:`U` 的前三个特征值占正谱质量的比例,中位下限。
///
/// # 这条闸在守什么
///
/// 整个算法的立论是"**光滑化之后的上限矩阵 `U` 本身就是一张画得出来的距离表**,
/// 所以不需要 RDKit 那一步逐对独立随机采样"。这句话是个**可以量的断言**,
/// 不是修辞 —— 这里就量它,并且与 RDKit 真正在做的事(在区间里随机取)
/// 在**同一张界矩阵**上对比。
///
/// 不设这条闸的话,界矩阵可以一边把宽度收得很漂亮(判据二全绿),
/// 一边把 `U` 变成一张摆不出来的表 —— 而那正好抽掉了立论的地基。
const MIN_FIT3: f64 = 0.85;

/// 判据三:`U` 的负特征值占谱绝对质量的比例,中位上限。
const MAX_NEG_SHARE: f64 = 0.08;

/// 判据三:`U` 必须比"区间内随机取"好多少倍(负份额之比)。
///
/// 只卡绝对值不够:万一哪天界矩阵整体变紧,两种取法会一起变好,
/// 而"换掉随机采样"这件事本身值不值就看不出来了。
const MIN_NEG_SHARE_GAIN: f64 = 2.0;

/// 判据三:嵌出来的坐标里,越界超过 0.1 Å 的原子对占比上限。
///
/// **这是给精修阶段的起点质量,也是判据三里唯一直接量几何的一条。**
/// `fit3` 与负份额都是谱上的代理量 —— 谱好看不等于坐标好用。
/// 这条闸先立在现值附近当棘轮,精修落地后应当能大幅收紧。
const MAX_EMBED_VIOL_FRAC: f64 = 30.0;

/// 按导出的连接表建分子。
///
/// **形式电荷必须一起带。** 头一版只带了原子序数,于是 `[NH3+]` 变成一个
/// 带四根键的中性氮,价键检查当场判死 —— 400 个分子里 **201 个**建不出来,
/// 而判据仍然"跑得通",只是在剩下那一半上量。**判据的样本被腰斩却不报警,
/// 比判据本身写错更危险。**
fn build_mol(
    z: &[u8],
    chg: &[i8],
    rad: &[u8],
    bonds: &[(u32, u32, u8, i64, i64, i64)],
) -> Option<MolBuilder> {
    let mut m = MolBuilder::new();
    for (k, &a) in z.iter().enumerate() {
        let mut ad = omgkit_core::AtomData::new(a);
        ad.formal_charge = chg.get(k).copied().unwrap_or(0);
        ad.num_radical_electrons = rad.get(k).copied().unwrap_or(0);
        m.add_atom_data(ad);
    }
    for &(i, j, o, _, _, _) in bonds {
        let ord = match o {
            2 => BondOrder::Double,
            3 => BondOrder::Triple,
            4 => BondOrder::Aromatic,
            _ => BondOrder::Single,
        };
        m.add_bond(i, j, ord).ok()?;
    }
    omgkit_chem::pipeline::sanitize(&mut m).ok()?;
    // 立体标记要在 sanitize **之后**写(它可能重排键),而且要在建界**之前**
    for (bi, &(_, _, _, st, sa0, sa1)) in bonds.iter().enumerate() {
        if sa0 < 0 || sa1 < 0 {
            continue;
        }
        // RDKit 的 Bond::BondStereo:0 无 2 Z 3 E 4 cis 5 trans
        let s = match st {
            2 => omgkit_core::BondStereo::Z,
            3 => omgkit_core::BondStereo::E,
            4 => omgkit_core::BondStereo::Cis,
            5 => omgkit_core::BondStereo::Trans,
            _ => continue,
        };
        #[allow(clippy::cast_possible_truncation)]
        if let Some(mut b) = m.bond_mut(bi as u32) {
            b.set_stereo(s);
            #[allow(clippy::cast_sign_loss, clippy::cast_possible_truncation)]
            b.set_stereo_atoms([sa0 as u32, sa1 as u32]);
        }
    }
    Some(m)
}

/// 模拟 RDKit 的 `pickRandomDistMat`:在每一对的区间里**各自独立**取一个值。
///
/// 这是判据三的对照组,不是产品路径 —— 所以用一个确定的 LCG 顶替均匀采样,
/// 判官才能次次跑出同一个数。取法本身与 RDKit 一致:`d = lb + r·(ub − lb)`。
fn pick_random_dist(b: &Bounds, seed: u64) -> Vec<f64> {
    let n = b.len();
    let mut d = vec![0.0; n * n];
    let mut st = seed;
    for i in 0..n {
        for j in (i + 1)..n {
            st = st.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            #[allow(clippy::cast_precision_loss)]
            let r = ((st >> 11) as f64) / ((1u64 << 53) as f64);
            let (lo, hi) = (b.lower(i, j), b.upper(i, j));
            let v = r.mul_add(hi - lo, lo);
            d[i * n + j] = v;
            d[j * n + i] = v;
        }
    }
    d
}

/// RDKit 的作废条件:`sqD0i[i] < EIGVAL_TOL` 且 `N > 3` 就把**整次尝试**判死
/// (`DistGeomUtils.cpp:110-115`)。`sqD0i` 就是 Gram 矩阵的对角。
///
/// 这里只是**记账**,看这条规则如果照搬会打掉多少分子 —— 本算法不用它。
fn rdkit_would_abort(t: &[f64], n: usize) -> bool {
    /// RDKit 的 `EIGVAL_TOL`。注意它是**绝对**阈值,单位是 Å²,
    /// 而同一个常数在那边还兼任特征值的零判据 —— 两个量纲不同的东西共用一个数。
    const EIGVAL_TOL: f64 = 0.001;
    n > 3 && (0..n).any(|i| t[i * n + i] < EIGVAL_TOL)
}

/// 嵌出来的坐标离界矩阵有多远:返回 `(越界的对数, 总对数, 最狠一处越界 Å)`。
///
/// # 为什么这个数比 `fit3` 重要
///
/// `fit3` 是**谱**上的指标 —— 它说的是"这张距离表离三维有多近",
/// 而下一阶段(精修)真正要面对的是**几何**:坐标违反了多少条约束、违反得多厉害。
/// 两者不是一回事:一张谱很漂亮的表照样可以嵌出一堆原子撞在一起的坐标。
fn violations(coords: &[[f64; 3]], b: &Bounds) -> (u64, u64, f64) {
    let n = b.len();
    let (mut bad, mut total, mut worst) = (0u64, 0u64, 0.0f64);
    for i in 0..n {
        for j in (i + 1)..n {
            let d = ((coords[i][0] - coords[j][0]).powi(2)
                + (coords[i][1] - coords[j][1]).powi(2)
                + (coords[i][2] - coords[j][2]).powi(2))
            .sqrt();
            let over = (b.lower(i, j) - d).max(d - b.upper(i, j)).max(0.0);
            total += 1;
            if over > 0.1 {
                bad += 1;
            }
            worst = worst.max(over);
        }
    }
    (bad, total, worst)
}

fn quantile(v: &[f64], f: f64) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let i = ((v.len() as f64 - 1.0) * f).round() as usize;
    v[i]
}

fn main() {
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "harness/baseline/rdkit_bounds.jsonl".into());
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| {
        eprintln!("读不了判官基准 {path}:{e}");
        eprintln!("先跑:python3 harness/dump_bounds.py harness/corpus/large.smi {path} 400");
        std::process::exit(1);
    });

    let (mut n, mut n_build_fail, mut n_infeasible) = (0u64, 0u64, 0u64);
    // 判据一:真实构象越界
    let (mut n_pairs, mut n_viol) = (0u64, 0u64);
    let mut worst_viol = 0.0f64;
    let mut worst_viol_case = String::new();
    // 判据二:界宽(按拓扑档分)
    let mut ratios: Vec<f64> = Vec::new();
    let (mut w_ours, mut w_rdkit) = (Vec::new(), Vec::new());
    // **按拓扑距离拆开** —— "整体 1.6 倍"没法指导修改,得知道松在哪一档
    let mut by_class: [(Vec<f64>, Vec<f64>); 5] = Default::default();
    // 判据三:U 摆不摆得进三维,以及它比"区间内随机取"好多少
    let (mut fit3_u, mut neg_u) = (Vec::new(), Vec::new());
    let (mut fit3_r, mut neg_r) = (Vec::new(), Vec::new());
    let (mut n_degenerate, mut n_neg_centroid, mut n_atoms) = (0u64, 0u64, 0u64);
    // 照搬 RDKit 那条作废条件的话会打掉多少分子(两张表各记一笔)
    let (mut n_abort_u, mut n_abort_r) = (0u64, 0u64);
    // 嵌出来的坐标违反了多少条界 —— 这是给精修阶段的起点质量
    let (mut v_bad_u, mut v_tot_u, mut v_bad_r, mut v_tot_r) = (0u64, 0u64, 0u64, 0u64);
    let (mut v_worst_u, mut v_worst_r): (Vec<f64>, Vec<f64>) = (Vec::new(), Vec::new());

    for line in text.lines() {
        let Ok(v) = serde_json::from_str::<serde_json::Value>(line) else {
            continue;
        };
        let smi = v["smiles"].as_str().unwrap_or("").to_string();
        let nat = v["n"].as_u64().unwrap_or(0) as usize;
        #[allow(clippy::cast_possible_truncation)]
        let z: Vec<u8> = v["z"]
            .as_array()
            .map(|a| {
                a.iter()
                    .filter_map(|x| x.as_u64().map(|y| y as u8))
                    .collect()
            })
            .unwrap_or_default();
        #[allow(clippy::cast_possible_truncation)]
        let bl: Vec<(u32, u32, u8, i64, i64, i64)> = v["bonds"]
            .as_array()
            .map(|a| {
                a.iter()
                    .filter_map(|e| {
                        let t = e.as_array()?;
                        Some((
                            t.first()?.as_u64()? as u32,
                            t.get(1)?.as_u64()? as u32,
                            t.get(2)?.as_u64()? as u8,
                            t.get(3).and_then(serde_json::Value::as_i64).unwrap_or(0),
                            t.get(4).and_then(serde_json::Value::as_i64).unwrap_or(-1),
                            t.get(5).and_then(serde_json::Value::as_i64).unwrap_or(-1),
                        ))
                    })
                    .collect()
            })
            .unwrap_or_default();
        let rd: Vec<f64> = v["smoothed"]
            .as_array()
            .map(|a| a.iter().filter_map(serde_json::Value::as_f64).collect())
            .unwrap_or_default();
        if z.len() != nat || rd.len() != nat * nat {
            continue;
        }
        n += 1;
        #[allow(clippy::cast_possible_truncation)]
        let chg: Vec<i8> = v["charge"]
            .as_array()
            .map(|a| {
                a.iter()
                    .filter_map(|x| x.as_i64().map(|y| y as i8))
                    .collect()
            })
            .unwrap_or_default();
        #[allow(clippy::cast_possible_truncation)]
        let rad: Vec<u8> = v["radical"]
            .as_array()
            .map(|a| {
                a.iter()
                    .filter_map(|x| x.as_u64().map(|y| y as u8))
                    .collect()
            })
            .unwrap_or_default();
        let Some(mol) = build_mol(&z, &chg, &rad, &bl) else {
            n_build_fail += 1;
            continue;
        };
        if mol.num_atoms() != nat {
            n_build_fail += 1;
            continue;
        }
        let (mut b, _) = bounds::build(&mol);
        if triangle_smooth(&mut b).is_err() {
            n_infeasible += 1;
            continue;
        }
        let Some(rdb) = Bounds::from_row_major(nat, rd) else {
            continue;
        };

        // 拓扑距离(封顶 4:0=自己 1=键 2=1-3 3=1-4 4=更远)
        let mut topo = vec![4u8; nat * nat];
        for start in 0..nat {
            let mut dist = vec![u8::MAX; nat];
            dist[start] = 0;
            let mut q = std::collections::VecDeque::from([start]);
            while let Some(x) = q.pop_front() {
                if dist[x] >= 3 {
                    continue;
                }
                let Ok(xu) = u32::try_from(x) else { continue };
                for (y, _) in mol.neighbors(xu) {
                    let y = y as usize;
                    if y < nat && dist[y] == u8::MAX {
                        dist[y] = dist[x] + 1;
                        q.push_back(y);
                    }
                }
            }
            for j in 0..nat {
                topo[start * nat + j] = dist[j].min(4);
            }
        }

        // ---- 判据三:U 能不能摆进三维,与"区间内随机取"同表对照 ----
        let du = reference_distances(&b);
        if let Ok(e) = embed(&du, nat) {
            fit3_u.push(e.fit3);
            neg_u.push(e.negative_share);
            if e.degenerate_axes > 0 {
                n_degenerate += 1;
            }
            n_neg_centroid += e.negative_centroid_sq as u64;
            n_atoms += nat as u64;
            let (bad, tot, worst) = violations(&e.coords, &b);
            v_bad_u += bad;
            v_tot_u += tot;
            v_worst_u.push(worst);
        }
        if rdkit_would_abort(&metric_matrix(&du, nat).0, nat) {
            n_abort_u += 1;
        }
        let dr = pick_random_dist(&b, 0xf00d);
        if let Ok(e) = embed(&dr, nat) {
            fit3_r.push(e.fit3);
            neg_r.push(e.negative_share);
            let (bad, tot, worst) = violations(&e.coords, &b);
            v_bad_r += bad;
            v_tot_r += tot;
            v_worst_r.push(worst);
        }
        if rdkit_would_abort(&metric_matrix(&dr, nat).0, nat) {
            n_abort_r += 1;
        }

        // ---- 判据二:界宽比 ----
        for i in 0..nat {
            for j in (i + 1)..nat {
                let wo = b.upper(i, j) - b.lower(i, j);
                let wr = rdb.upper(i, j) - rdb.lower(i, j);
                w_ours.push(wo);
                w_rdkit.push(wr);
                if wr > 1e-9 {
                    ratios.push(wo / wr);
                }
                let c = topo[i * nat + j] as usize;
                by_class[c].0.push(wo);
                by_class[c].1.push(wr);
            }
        }

        // ---- 判据一:真实构象必须落在界内 ----
        let coords: Vec<[f64; 3]> = v["coords"]
            .as_array()
            .map(|a| {
                a.iter()
                    .filter_map(|e| {
                        let t = e.as_array()?;
                        Some([
                            t.first()?.as_f64()?,
                            t.get(1)?.as_f64()?,
                            t.get(2)?.as_f64()?,
                        ])
                    })
                    .collect()
            })
            .unwrap_or_default();
        if coords.len() != nat {
            continue; // 这个分子没嵌出构象,判据一跳过(判据二已经记了)
        }
        for i in 0..nat {
            for j in (i + 1)..nat {
                let d = ((coords[i][0] - coords[j][0]).powi(2)
                    + (coords[i][1] - coords[j][1]).powi(2)
                    + (coords[i][2] - coords[j][2]).powi(2))
                .sqrt();
                let over = (b.lower(i, j) - d).max(d - b.upper(i, j)).max(0.0);
                n_pairs += 1;
                if over > MAX_VIOLATION {
                    n_viol += 1;
                }
                if over > worst_viol {
                    worst_viol = over;
                    worst_viol_case = format!("{smi}  第 {i}/{j} 对  实距 {d:.3}");
                }
            }
        }
    }

    ratios.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    w_ours.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    w_rdkit.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    #[allow(clippy::cast_precision_loss)]
    let viol_frac = if n_pairs == 0 {
        0.0
    } else {
        n_viol as f64 / n_pairs as f64
    };

    println!("判官:界矩阵,分子 {n} 个(建不出来 {n_build_fail}、界不可行 {n_infeasible})");
    println!("  ── 判据一:真实构象落在界内 ──");
    println!(
        "    原子对 {n_pairs},越界 >{MAX_VIOLATION} Å 的 {n_viol}({:.3}%,上限 {:.1}%)",
        100.0 * viol_frac,
        100.0 * MAX_VIOLATION_FRAC
    );
    println!("    最狠一处越界 {worst_viol:.3} Å  {worst_viol_case}");
    println!("  ── 判据二:界宽不许比 RDKit 松 ──");
    println!(
        "    我们的界宽 中位 {:.3} / p90 {:.3};RDKit 中位 {:.3} / p90 {:.3}",
        quantile(&w_ours, 0.5),
        quantile(&w_ours, 0.9),
        quantile(&w_rdkit, 0.5),
        quantile(&w_rdkit, 0.9)
    );
    println!(
        "    逐对宽度比(我们/RDKit)中位 {:.3}(上限 {MAX_WIDTH_RATIO});p90 {:.3}",
        quantile(&ratios, 0.5),
        quantile(&ratios, 0.9)
    );
    println!("    ── 按拓扑距离拆开(中位界宽) ──");
    for (c, name) in [
        (1usize, "1-2 键"),
        (2, "1-3 角"),
        (3, "1-4 扭转"),
        (4, "≥1-5"),
    ] {
        let (o, r) = &mut by_class[c];
        o.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        r.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        // **中位会把子群里的改进整个藏起来。** 实测:把芳环上的 1-4 扭转从
        // "顺式到反式的区间"改成确定值之后,苯的 1-4 宽度掉到 0.05 以下,
        // 而语料的 1-4 中位**一位没动**(0.762)—— 因为芳环上的 1-4 只占少数,
        // 中位落在 sp³ 链上的氢那一群里。所以这里同时报"被钉住的比例"。
        // 阈值取 0.15:RDKit 把"解掉了析取"的对钉成宽度 **2×GEN_DIST_TOL = 0.12 Å**
        // (`BoundsMatrixBuilder.cpp:32` 与 `:1011-1012`)。头一版这里取 0.05,
        // 低于它的钉死宽度,于是把 RDKit 显示成"钉住 0.0%" —— **阈值定错了,
        // 结论就反了**:它其实把过半的 1-4 都钉住了,那正是它中位 0.120 的来源。
        let pinned = o.iter().filter(|w| **w < 0.15).count();
        let pinned_r = r.iter().filter(|w| **w < 0.15).count();
        #[allow(clippy::cast_precision_loss)]
        let (pf, pfr) = (
            100.0 * pinned as f64 / o.len().max(1) as f64,
            100.0 * pinned_r as f64 / r.len().max(1) as f64,
        );
        println!(
            "      {name:9} 对数 {:7}  中位 我们 {:6.3} / RDKit {:6.3}(比 {:5.2})  钉住 我们 {pf:5.1}% / RDKit {pfr:5.1}%",
            o.len(),
            quantile(o, 0.5),
            quantile(r, 0.5),
            quantile(o, 0.5) / quantile(r, 0.5).max(1e-9)
        );
    }

    for v in [&mut fit3_u, &mut neg_u, &mut fit3_r, &mut neg_r] {
        v.sort_by(f64::total_cmp);
    }
    let (m_fit3, m_neg) = (quantile(&fit3_u, 0.5), quantile(&neg_u, 0.5));
    let (r_fit3, r_neg) = (quantile(&fit3_r, 0.5), quantile(&neg_r, 0.5));
    println!("  ── 判据三:U 摆不摆得进三维(对照组 = RDKit 那样在区间里随机取)──");
    println!(
        "    前三特征值占正谱质量 中位  我们 {m_fit3:.3} / 随机取 {r_fit3:.3}(下限 {MIN_FIT3})"
    );
    println!(
        "    负特征值占谱绝对质量 中位  我们 {m_neg:.3} / 随机取 {r_neg:.3}(上限 {MAX_NEG_SHARE},且须优于对照 {MIN_NEG_SHARE_GAIN} 倍)"
    );
    #[allow(clippy::cast_precision_loss)]
    let neg_frac = 100.0 * n_neg_centroid as f64 / n_atoms.max(1) as f64;
    println!("    有退化轴的分子 {n_degenerate} 个;到质心平方距离出负的原子 {n_neg_centroid}/{n_atoms}({neg_frac:.2}%)");
    // **这一行是记账,不是判据。** RDKit 只要发现一个原子的 sqD0i 低于阈值就把
    // 整次尝试作废(`DistGeomUtils.cpp:110-115`),而那是个"全或无"的规则:
    // 一个原子的毛病判死整个分子。本算法不作废,只把差距如实报出来。
    println!(
        "    照搬 RDKit 那条作废条件的话会打掉:U 表 {n_abort_u}/{n} 个分子、随机取 {n_abort_r}/{n} 个"
    );
    v_worst_u.sort_by(f64::total_cmp);
    v_worst_r.sort_by(f64::total_cmp);
    #[allow(clippy::cast_precision_loss)]
    let (fu, fr) = (
        100.0 * v_bad_u as f64 / v_tot_u.max(1) as f64,
        100.0 * v_bad_r as f64 / v_tot_r.max(1) as f64,
    );
    // **这一行才是给精修阶段的起点质量**,`fit3` 只是谱上的代理量。
    println!("    嵌出来的坐标越界 >0.1 Å 的对:U 表 {fu:.1}% / 随机取 {fr:.1}%");
    println!(
        "      每分子最狠越界 中位  U 表 {:.2} Å / 随机取 {:.2} Å(上限 {MAX_EMBED_VIOL_FRAC:.0}%)",
        quantile(&v_worst_u, 0.5),
        quantile(&v_worst_r, 0.5)
    );

    let mut fatal = false;
    if n == 0 {
        eprintln!("\n一个分子都没读到 —— 基准文件是空的?");
        fatal = true;
    }
    // **样本被腰斩不许静悄悄。** 建不出来的分子会被跳过,而判据照样打印一个
    // 好看的百分比 —— 那个百分比是在剩下的分子上量的。头一版这里没有闸,
    // 实测 400 个里 201 个建不出来(漏了形式电荷),判据一照样报 0.4%。
    if n_infeasible > MAX_INFEASIBLE {
        eprintln!(
            "\n有 {n_infeasible} 个分子的界自相矛盾(光滑化判不可行),超过上限 {MAX_INFEASIBLE} —— 约束之间在打架"
        );
        fatal = true;
    }
    if n_build_fail > MAX_BUILD_FAIL {
        eprintln!(
            "\n有 {n_build_fail} 个分子按连接表建不出来,超过上限 {MAX_BUILD_FAIL} —— 判据是在剩下的分子上量的,那个数不作数"
        );
        fatal = true;
    }
    if viol_frac > MAX_VIOLATION_FRAC {
        eprintln!(
            "\n真实构象越界的原子对占 {:.3}%,超过上限 —— 界矩阵把真实几何排除在外了",
            100.0 * viol_frac
        );
        fatal = true;
    }
    // **这一条是防"把界写宽了蒙混过关"的。** 判据一是单向的:界越宽越容易过。
    if quantile(&ratios, 0.5) > MAX_WIDTH_RATIO {
        eprintln!(
            "\n界宽中位比 {:.3} > {MAX_WIDTH_RATIO} —— 比 RDKit 还松,判据一那条绿是靠放宽换来的",
            quantile(&ratios, 0.5)
        );
        fatal = true;
    }
    // **判据三守的是算法的立论。** 判据一二只管界矩阵自己好不好,管不到
    // "拿 U 当参考距离表"这个决定还成不成立 —— 界可以一边收紧一边变得摆不出来。
    if m_fit3 < MIN_FIT3 {
        eprintln!("\nU 的三维贴合度中位 {m_fit3:.3} < {MIN_FIT3} —— 拿 U 当参考距离表这个前提塌了");
        fatal = true;
    }
    if m_neg > MAX_NEG_SHARE {
        eprintln!(
            "\nU 的负特征值份额中位 {m_neg:.3} > {MAX_NEG_SHARE} —— U 正在离能摆出来越来越远"
        );
        fatal = true;
    }
    if fu > MAX_EMBED_VIOL_FRAC {
        eprintln!("\n嵌出来的坐标越界的对占 {fu:.1}% > {MAX_EMBED_VIOL_FRAC:.0}% —— 起点太差,精修要从很深的坑里爬");
        fatal = true;
    }
    if m_neg * MIN_NEG_SHARE_GAIN > r_neg {
        eprintln!(
            "\nU 的负份额 {m_neg:.3} 相对随机取 {r_neg:.3} 不足 {MIN_NEG_SHARE_GAIN} 倍优势 —— 换掉随机采样这一步不再值得"
        );
        fatal = true;
    }
    if fatal {
        std::process::exit(1);
    }
    println!("\n三条都过。");
}
