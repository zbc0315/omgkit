//! **外部判官:三角光滑化对不对,由 RDKit 说了算。**
//!
//! 自己写的光滑化拿自己的性质去验("再光滑一次不变"、"满足三角不等式"),
//! 两边同错就同绿 —— 那两条性质对一个**什么都不做**的实现同样成立。
//! 所以要一个不知道我们怎么写的判官。
//!
//! 做法:`harness/dump_bounds.py` 把 RDKit 的界矩阵导两份 ——
//! `raw`(未光滑)与 `smoothed`(RDKit 自己光滑过的)。这里把 `raw` 喂给
//! [`omgkit_conf::smooth::triangle_smooth`],结果与 `smoothed` **逐元素**比。
//!
//! ```shell
//! python3 harness/dump_bounds.py harness/corpus/large.smi harness/baseline/rdkit_bounds.jsonl 400
//! cargo run -p omgkit-conf --release --example smooth_oracle -- harness/baseline/rdkit_bounds.jsonl
//! ```

use omgkit_conf::smooth::{triangle_smooth, Bounds};

/// 逐元素允许的绝对偏差(Å)。**不是容差,是浮点噪声的余量** ——
/// 两边做的是同一串加减比大小,理应到机器精度。
const TOL: f64 = 1e-9;

/// 允许有多少个分子对不上。**贴着现值**:光滑化是纯算术,没有"大致相同"这回事。
const MAX_MISMATCH: usize = 0;

/// 基准里允许有几行**没真正比到**。
///
/// # 分母也要有闸
///
/// 先前 `n += 1` 排在 `Bounds::from_row_major` 那个 `continue` **之前**,于是
/// 读到的分子照样计进 `n`,而 `triangle_smooth` 一次都没被调用。变异实测
/// (让 `from_row_major` 恒返回 `None`):
///
/// ```text
/// 判官:RDKit 的界矩阵(未光滑 → 光滑),分子 27 个
///   逐元素最大偏差 0.000e0(上限 1e-9)
///   对不上的分子 0 个(上限 0)
/// 与 RDKit 逐元素一致。            ← 退出码 0
/// ```
///
/// 那个"完美"的 `0.000e0` 没有任何闸在读。
///
/// 现在 `n` 数的是**真正比到的分子数**(计数点挪到全部 `continue` 之后),
/// 于是两档各有各的闸:
///
/// - **全都没比到** → `n == 0`,由原有那条挡住(先前 `n` 恒等于行数,它挡不住);
/// - **只有一部分没比到** → 这条 `MAX_UNCOMPARED`。变异实测(让
///   `from_row_major` 只对超过 25 个原子的分子返回 `None`):
///   "真正比到 21 个""6 行没真正比到",退 1。
///
/// 实测冒烟档 27 行比到 27 个,不多不少,所以钉死 0。
///
/// **基准本身被截短**是另一半 —— 判官无从知道自己本该收到 27 行。
/// 那一半由 `crates/omgkit-conf/tests/baseline_sizes.rs` 把行数当契约钉住。
const MAX_UNCOMPARED: u64 = 0;

fn main() {
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "harness/baseline/rdkit_bounds.jsonl".into());
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| {
        eprintln!("读不了判官基准 {path}:{e}");
        eprintln!("先跑:python3 harness/dump_bounds.py harness/corpus/large.smi {path} 400");
        std::process::exit(1);
    });

    // `n_lines` 是基准里的分子数(分母),`n` 是**真正比到**的(分子)。
    let (mut n_lines, mut n, mut n_mismatch, mut n_infeasible) = (0u64, 0u64, 0usize, 0u64);
    let mut worst = 0.0f64;
    let mut worst_case: Option<(f64, String, usize, usize)> = None;
    let mut worst_mismatch_smi = String::new();

    for line in text.lines() {
        // **行计数排在解析之前。** 排在之后的话,一行非法 JSON 既不进 `n_lines`
        // 也不进任何一档 —— 判官与 `tests/baseline_sizes.rs` 的行数契约会**同时**
        // 失明(那条契约数的是非空行,非法行照样算一行)。放在这里,非法行会
        // 落进 `unread`,被分母闸抓住。
        if line.trim().is_empty() {
            continue;
        }
        n_lines += 1;
        let Ok(v) = serde_json::from_str::<serde_json::Value>(line) else {
            continue;
        };
        let smi = v["smiles"].as_str().unwrap_or("").to_string();
        let nat = v["n"].as_u64().unwrap_or(0) as usize;
        let raw: Vec<f64> = v["raw"]
            .as_array()
            .map(|a| a.iter().filter_map(serde_json::Value::as_f64).collect())
            .unwrap_or_default();
        let want: Vec<f64> = v["smoothed"]
            .as_array()
            .map(|a| a.iter().filter_map(serde_json::Value::as_f64).collect())
            .unwrap_or_default();
        if raw.len() != nat * nat || want.len() != nat * nat {
            continue;
        }
        let Some(mut b) = Bounds::from_row_major(nat, raw) else {
            continue;
        };
        // **`n` 必须在所有 `continue` 之后加** —— 它是"真正比到的分子数",
        // 不是"读到的行数"。见 `MAX_UNCOMPARED`。
        n += 1;
        if triangle_smooth(&mut b).is_err() {
            // RDKit 导得出 smoothed 就说明它认为可行 —— 我们判不可行是**分歧**
            n_infeasible += 1;
            n_mismatch += 1;
            continue;
        }
        let got = b.as_row_major();
        let mut bad = false;
        for k in 0..(nat * nat) {
            let d = (got[k] - want[k]).abs();
            if d > worst {
                worst = d;
                worst_case = Some((d, smi.clone(), k / nat, k % nat));
            }
            if d > TOL {
                bad = true;
            }
        }
        if bad {
            n_mismatch += 1;
            if worst_mismatch_smi.is_empty() {
                worst_mismatch_smi = smi.clone();
            }
        }
    }

    let uncompared = n_lines.saturating_sub(n);
    println!("判官:RDKit 的界矩阵(未光滑 → 光滑),基准 {n_lines} 行,真正比到 {n} 个");
    println!("  逐元素最大偏差 {worst:.3e}(上限 {TOL:.0e})");
    println!("  对不上的分子 {n_mismatch} 个(上限 {MAX_MISMATCH});我们判不可行而 RDKit 判可行的 {n_infeasible} 个");
    if let Some((d, smi, i, j)) = &worst_case {
        println!("    偏差最大处:{d:.3e}  第 {i}/{j} 项  {smi}");
    }
    if !worst_mismatch_smi.is_empty() {
        println!("    第一个对不上的分子:{worst_mismatch_smi}");
    }

    if n == 0 {
        eprintln!("\n一个分子都没读到 —— 基准文件是空的?");
        std::process::exit(1);
    }
    if uncompared > MAX_UNCOMPARED {
        eprintln!(
            "\n基准里有 {uncompared} 行没真正比到(上限 {MAX_UNCOMPARED})—— \
             对不上的个数是分子,覆盖面是分母。少比几个分子,上面那些数都会变好看"
        );
        std::process::exit(1);
    }
    if n_mismatch > MAX_MISMATCH {
        eprintln!("\n有 {n_mismatch} 个分子与 RDKit 的光滑结果对不上,超过上限 {MAX_MISMATCH}");
        std::process::exit(1);
    }
    println!("\n与 RDKit 逐元素一致。");
}
