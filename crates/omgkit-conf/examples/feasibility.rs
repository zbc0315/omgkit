//! **全语料的界可行率** —— 一个分子在进嵌入之前就死掉的比例。
//!
//! # 为什么单独立一个
//!
//! `bounds_oracle` 只跑 400 个分子(受 RDKit 导出基准的规模所限),而这个数
//! 是整个项目的**头号指标**:要赢的正是 RDKit 那 0.52% 的失败率。
//! 界矩阵自相矛盾的分子连嵌入都进不去,它直接计入失败。
//!
//! 这个判据不需要任何外部基准,只要一份 SMILES 语料:
//!
//! ```shell
//! cargo run -p omgkit-conf --release --example feasibility -- harness/corpus/large.smi
//! ```

use omgkit_conf::{bounds, smooth::triangle_smooth};

/// 建界之后就有区间是空的(下限 > 上限)——**上限 0,一个都不许有**。
///
/// 空区间意味着**参数表自相矛盾**,与几何无关:两条约束各自查表、各带自己的容差,
/// 凑不到一起。这类问题必须当场暴露,不能让它伪装成"这个分子摆不出来"。
///
/// 实测从 21 个降到 0 个,靠的是两条:1-3 不再往已经成键的那一对上写
/// (直接量的键长比角推的距离可靠),以及两条 1-3 估计交空时退回并集并记账。
const MAX_EMPTY: u64 = 0;

/// 光滑化判不可行的分子占比上限。**这是整个项目的头号指标。**
///
/// # 闸为什么必须比目标严
///
/// 要赢的是 RDKit ETKDG 的 **0.52%** 失败率;界不可行的分子连嵌入都进不去,
/// 直接计入失败。所以这一项只要碰到 0.52%,算法在还没开始的地方就已经输了。
///
/// 先前这条闸在 `bounds_oracle` 里定的是 **2%** —— 比目标松了近四倍,
/// 于是判据一路是绿的,而实测 0.65% 早已越过目标。**闸松于目标等于没有闸。**
///
/// 现在:0.65% → 0.54%(修了角的回退顺序)→ **0.34%**(修了 1-3 的写入规则)。
/// 闸设在 0.40%,是贴着现值的棘轮,离 0.52% 的目标还留着余量。
///
/// 另外这条闸放在**全语料 8831 个分子**上,不放在 400 个分子的判官里:
/// 400 个样本上真实率 0.34% 只对应 1.4 个分子,泊松噪声足以让闸随机红绿。
const MAX_INFEASIBLE_FRAC: f64 = 0.0040;

fn main() {
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "harness/corpus/large.smi".into());
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| {
        eprintln!("读不了语料 {path}:{e}");
        std::process::exit(1);
    });

    let (mut n, mut n_parse_fail, mut n_empty, mut n_infeasible) = (0u64, 0u64, 0u64, 0u64);
    let mut empty_cases: Vec<String> = Vec::new();
    let mut infeasible_cases: Vec<String> = Vec::new();

    for line in text.lines() {
        let smi = line.split('\t').next().unwrap_or("").trim();
        if smi.is_empty() {
            continue;
        }
        let Ok(mut mol) = omgkit_io::smiles::parse(smi) else {
            n_parse_fail += 1;
            continue;
        };
        if omgkit_chem::pipeline::sanitize(&mut mol).is_err() {
            n_parse_fail += 1;
            continue;
        }
        // 补氢要给一个与写法无关的秩;这里只关心界可不可行,用恒等秩即可
        let order: Vec<u32> = (0..mol.num_atoms() as u32).collect();
        omgkit_chem::explicit_hs::add_explicit_hs(&mut mol, &order);
        n += 1;
        let (mut b, _) = bounds::build(&mol);
        // 建完界先看有没有区间当场就空 —— 那是**表自相矛盾**,与几何无关
        let nat = b.len();
        let mut empty = false;
        for i in 0..nat {
            for j in (i + 1)..nat {
                if b.lower(i, j) > b.upper(i, j) {
                    empty = true;
                }
            }
        }
        if empty {
            n_empty += 1;
            if empty_cases.len() < 8 {
                empty_cases.push(smi.to_string());
            }
        }
        if triangle_smooth(&mut b).is_err() {
            n_infeasible += 1;
            if infeasible_cases.len() < 8 {
                infeasible_cases.push(smi.to_string());
            }
        }
    }

    #[allow(clippy::cast_precision_loss)]
    let (pe, pi) = (
        100.0 * n_empty as f64 / n.max(1) as f64,
        100.0 * n_infeasible as f64 / n.max(1) as f64,
    );
    println!("语料 {path}:建界 {n} 个分子(解析/净化失败 {n_parse_fail})");
    println!("  建界即空区间   {n_empty}({pe:.2}%)");
    println!("  光滑化判不可行 {n_infeasible}({pi:.2}%)  ← 这些分子连嵌入都进不去");
    println!("  RDKit ETKDG 的失败率是 0.52% —— 这一行必须低于它,否则算法还没开始就输了");
    if !empty_cases.is_empty() {
        println!("  空区间的例子:{}", empty_cases.join("  "));
    }
    if !infeasible_cases.is_empty() {
        println!("  不可行的例子:{}", infeasible_cases.join("  "));
    }

    let mut fatal = false;
    if n == 0 {
        eprintln!("\n一个分子都没读到 —— 语料是空的?");
        fatal = true;
    }
    if n_empty > MAX_EMPTY {
        eprintln!(
            "\n有 {n_empty} 个分子建完界就有空区间(上限 {MAX_EMPTY})—— 参数表自相矛盾,\
             别让它伪装成分子不可行"
        );
        fatal = true;
    }
    if pi / 100.0 > MAX_INFEASIBLE_FRAC {
        eprintln!(
            "\n界不可行 {pi:.2}% > {:.2}% —— 这是头号指标,要赢的是 RDKit 的 0.52%",
            100.0 * MAX_INFEASIBLE_FRAC
        );
        fatal = true;
    }
    if fatal {
        std::process::exit(1);
    }
    println!("\n两条都过。");
}
