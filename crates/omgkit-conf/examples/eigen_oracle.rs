//! **外部判官:特征分解与嵌入对不对。两条闸,问的是两件事。**
//!
//! | 闸 | 问什么 | 判官是谁 |
//! |---|---|---|
//! | **判据一** | 特征值算得对不对 | numpy `linalg.eigvalsh`(LAPACK,分治/QR) |
//! | **判据二** | 嵌入公式对不对 | 几何本身:真实构象的精确距离表必须能精确嵌回去 |
//!
//! # 为什么两条都要
//!
//! 判据一比的是**特征分解**。Gram 矩阵的公式两边写的是同一个,所以它管不住公式 ——
//! 公式写错了两边一起错,判据一照样全绿。
//!
//! 判据二补的正是这一块,而且**它不依赖任何外部实现**:一组真实的三维坐标算出
//! 精确距离表,再嵌回去,还原出来的两两距离必须逐对相同。这条把
//! `d_0i²` 与 `T_ij` 两个公式完全钉死 —— 任何一项写错都会当场红。
//!
//! 反过来判据二也管不住判据一:一个精度很差的特征分解在"距离表本来就精确是三维"
//! 的输入上照样能还原得不错,因为那时候第 4 个及以后的特征值真的是 0,
//! 数值噪声无处藏身。真实流水线喂进去的 `U` 不是这样 —— 它带着一整条负尾巴。
//!
//! ```shell
//! python3 harness/dump_gram.py harness/baseline/rdkit_bounds.jsonl harness/baseline/gram_eigs.jsonl
//! cargo run -p omgkit-conf --release --example eigen_oracle
//! ```

use omgkit_conf::embed::{embed, metric_matrix};
use omgkit_conf::linalg::symmetric_eigen;

/// 特征值与 LAPACK 的相对偏差上限(以谱半径归一)。
///
/// Jacobi 与 LAPACK 是完全不同的算法,不可能逐位相同;两边的绝对误差都是
/// `eps·‖A‖` 量级,归一之后约 1e-14。闸设在 1e-9 留了五个数量级的余量 ——
/// 它拦的是"算错了",不是"最后几位不一样"。
const MAX_EIG_DEV: f64 = 1e-9;

/// 精确回嵌之后,两两距离允许差多少(Å)。
///
/// 这一条是**几何恒等式**,不是统计:输入距离表本来就是三维的,嵌入必须精确还原。
/// 实际偏差在 1e-12 Å 量级,闸设在 1e-6 Å。
const MAX_ROUNDTRIP_DEV: f64 = 1e-6;

/// 允许有几个分子对不上。**样本被腰斩不许静悄悄** —— 这是先前栽过两次的坑:
/// 一次是形式电荷没导致 400 个分子里 201 个建不出来,一次是 1-5 链式约束
/// 把参与统计的原子对从 108981 打到 9317,而剩下那一成的比值反而更好看。
const MAX_MISSING: u64 = 0;

fn quantile(v: &[f64], f: f64) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let i = ((v.len() as f64 - 1.0) * f).round() as usize;
    v[i]
}

fn floats(v: &serde_json::Value) -> Vec<f64> {
    v.as_array()
        .map(|a| a.iter().filter_map(serde_json::Value::as_f64).collect())
        .unwrap_or_default()
}

/// 两张距离表的最大逐对偏差。
fn max_dist_dev(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .fold(0.0_f64, |w, (x, y)| w.max((x - y).abs()))
}

/// 由坐标算精确距离表。
fn distances(x: &[[f64; 3]]) -> Vec<f64> {
    let n = x.len();
    let mut d = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            d[i * n + j] = (0..3)
                .map(|k| (x[i][k] - x[j][k]).powi(2))
                .sum::<f64>()
                .sqrt();
        }
    }
    d
}

/// 比一组特征值,返回以谱半径归一的最大偏差。
fn eig_dev(mine: &[f64], theirs: &[f64]) -> Option<f64> {
    if mine.len() != theirs.len() || mine.is_empty() {
        return None;
    }
    let scale = theirs
        .iter()
        .fold(0.0_f64, |w, v| w.max(v.abs()))
        .max(1e-12);
    Some(
        mine.iter()
            .zip(theirs)
            .fold(0.0_f64, |w, (a, b)| w.max((a - b).abs() / scale)),
    )
}

fn main() {
    let bounds_path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "harness/baseline/rdkit_bounds.jsonl".into());
    let eig_path = std::env::args()
        .nth(2)
        .unwrap_or_else(|| "harness/baseline/gram_eigs.jsonl".into());

    let read = |p: &str| {
        std::fs::read_to_string(p).unwrap_or_else(|e| {
            eprintln!("读不了判官基准 {p}:{e}");
            eprintln!("先跑:python3 harness/dump_gram.py harness/baseline/rdkit_bounds.jsonl harness/baseline/gram_eigs.jsonl");
            std::process::exit(1);
        })
    };
    let bounds_text = read(&bounds_path);
    let eig_text = read(&eig_path);

    let mut n = 0u64;
    let mut n_missing = 0u64;
    // 判据一
    let mut dev_u: Vec<f64> = Vec::new();
    let mut dev_x: Vec<f64> = Vec::new();
    let mut worst_eig = (0.0f64, String::new());
    // 判据二
    let mut roundtrip: Vec<f64> = Vec::new();
    let mut worst_rt = (0.0f64, String::new());

    for (bl, el) in bounds_text.lines().zip(eig_text.lines()) {
        let (Ok(bv), Ok(ev)) = (
            serde_json::from_str::<serde_json::Value>(bl),
            serde_json::from_str::<serde_json::Value>(el),
        ) else {
            n_missing += 1;
            continue;
        };
        // 两个文件必须是同一批分子、同一个顺序 —— 错位之后判据照样"跑得通",
        // 只是在比别的分子
        if bv["smiles"] != ev["smiles"] {
            n_missing += 1;
            continue;
        }
        let smi = bv["smiles"].as_str().unwrap_or("").to_string();
        let Some(nat) = bv["n"].as_u64().map(|x| x as usize) else {
            n_missing += 1;
            continue;
        };
        let sm = floats(&bv["smoothed"]);
        if sm.len() != nat * nat {
            n_missing += 1;
            continue;
        }
        n += 1;

        // —— 判据一(a):光滑化上限矩阵 U 建的 Gram
        let mut u = vec![0.0; nat * nat];
        for i in 0..nat {
            for j in (i + 1)..nat {
                let v = sm[i * nat + j];
                u[i * nat + j] = v;
                u[j * nat + i] = v;
            }
        }
        let (t, _) = metric_matrix(&u, nat);
        let eig = symmetric_eigen(&t, nat).expect("光滑化后的界矩阵不该含非有限数");
        if let Some(d) = eig_dev(&eig.values, &floats(&ev["eig_u"])) {
            dev_u.push(d);
            if d > worst_eig.0 {
                worst_eig = (d, format!("{smi}(U)"));
            }
        } else {
            n_missing += 1;
        }

        // —— 判据一(b)与判据二:真实构象
        let coords: Vec<[f64; 3]> = bv["coords"]
            .as_array()
            .map(|a| {
                a.iter()
                    .filter_map(|p| {
                        let q = p.as_array()?;
                        Some([
                            q.first()?.as_f64()?,
                            q.get(1)?.as_f64()?,
                            q.get(2)?.as_f64()?,
                        ])
                    })
                    .collect()
            })
            .unwrap_or_default();
        if coords.len() != nat {
            continue;
        }
        let dx = distances(&coords);
        let (tx, _) = metric_matrix(&dx, nat);
        let eigx = symmetric_eigen(&tx, nat).expect("真实坐标不该含非有限数");
        if let Some(d) = eig_dev(&eigx.values, &floats(&ev["eig_x"])) {
            dev_x.push(d);
            if d > worst_eig.0 {
                worst_eig = (d, format!("{smi}(X)"));
            }
        }
        let back = embed(&dx, nat).expect("真实构象必须嵌得出来");
        let dev = max_dist_dev(&dx, &distances(&back.coords));
        roundtrip.push(dev);
        if dev > worst_rt.0 {
            worst_rt = (dev, smi);
        }
    }

    dev_u.sort_by(f64::total_cmp);
    dev_x.sort_by(f64::total_cmp);
    roundtrip.sort_by(f64::total_cmp);

    println!("判官:特征分解与嵌入,分子 {n} 个(对不上 {n_missing} 个)");
    println!("  判据一:特征值 vs numpy/LAPACK,以谱半径归一");
    println!(
        "    U 建的 Gram   中位 {:.2e}  最大 {:.2e}  (样本 {})",
        quantile(&dev_u, 0.5),
        dev_u.last().copied().unwrap_or(0.0),
        dev_u.len()
    );
    println!(
        "    真实构象 Gram 中位 {:.2e}  最大 {:.2e}  (样本 {})",
        quantile(&dev_x, 0.5),
        dev_x.last().copied().unwrap_or(0.0),
        dev_x.len()
    );
    println!(
        "    最差:{} = {:.2e}(上限 {MAX_EIG_DEV:.0e})",
        worst_eig.1, worst_eig.0
    );
    println!("  判据二:真实构象精确回嵌,两两距离偏差(Å)");
    println!(
        "    中位 {:.2e}  最大 {:.2e}  最差分子 {}(上限 {MAX_ROUNDTRIP_DEV:.0e})",
        quantile(&roundtrip, 0.5),
        roundtrip.last().copied().unwrap_or(0.0),
        worst_rt.1
    );

    let mut bad = false;
    if n_missing > MAX_MISSING {
        println!("\n【判据红】对不上的分子 {n_missing} 个,上限 {MAX_MISSING}");
        bad = true;
    }
    if n == 0 {
        println!("\n【判据红】一个分子都没量到");
        bad = true;
    }
    if worst_eig.0 > MAX_EIG_DEV {
        println!(
            "\n【判据一红】特征值偏差 {:.3e} 超过上限 {MAX_EIG_DEV:.0e}({})",
            worst_eig.0, worst_eig.1
        );
        bad = true;
    }
    if worst_rt.0 > MAX_ROUNDTRIP_DEV {
        println!(
            "\n【判据二红】回嵌距离偏差 {:.3e} Å 超过上限 {MAX_ROUNDTRIP_DEV:.0e}({})",
            worst_rt.0, worst_rt.1
        );
        bad = true;
    }
    if bad {
        std::process::exit(1);
    }
    println!("\n两条判据都过。");
}
