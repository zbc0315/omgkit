//! 规范化的增长曲线守卫:每原子耗时不随分子规模上升。
//!
//! 与 `omgkit-chem` 的 `scaling.rs` 同一套判据和同样的两个坑(规模要够大、
//! 测量要交错),那边的文件头有完整说明。这里只记本模块特有的压力形状。
//!
//! # 两种形状压的是两处不同的平方项
//!
//! | 形状 | 压的是什么 |
//! |---|---|
//! | 很多个小环系串成长链 | **细化的轮数**。朴素的"每轮重算全部颜色"要跑 O(直径) 轮,长链上就是平方 |
//! | 单个大环 | **打破对称的分支数**。第一格含全部 n 个原子,逐个试就是平方 |
//! | 带方向键的共轭多烯 | **方向判定里那趟细化**。前两种形状里它被提前返回跳过了 |
//!
//! 前两个平方项都真实存在过,而且**互相看不见**:
//!
//! - 长链上第一格很小,分支数的问题显不出来
//! - 大环上直径虽大但细化一轮就稳定(所有原子等价),轮数的问题也显不出来
//!
//! 第三种形状是量出来补的:`informative_directions` 在分子没有方向键时直接
//! 返回,而前两种形状**一条方向键都没有** —— 守卫从没走进过那段代码。
//! 实测带方向键时它占写出耗时的 64~74%,却完全没有被守着。
//!
//! 这是"判据空过"的又一个样子:测试跑得好好的,只是没碰到要守的那条路。
//!
//! 实测过:没有自同构剪枝时,800 元大环要 199 ms;有了剪枝是 0.81 ms。
//! 长链那一侧则靠分裂器工作表 + Hopcroft 技巧压到 O((n+m) log n)。

use std::time::{Duration, Instant};

use omgkit_io::{canon, smiles};

const ROUNDS: usize = 8;
/// 每原子耗时相对**理论增长**的最大允许倍数。
///
/// # 分母为什么不是 1
///
/// 规范化是 `O((n+m) log n)`(文件头自己写着),所以**每原子耗时本来就随 `log n` 涨**。
/// 原先这里的判据是"每原子耗时应基本持平",与实现的复杂度对不上 ——
/// 长链那一档 699 → 2799 原子,光 `log` 因子就是 **1.21**,而闸设在 1.35,
/// 只剩 12% 的余量留给缓存与测量噪声。实测它三次里红两次,
/// 而逐档数据显示第一次翻倍只涨 6~10%、第二次翻倍涨 32% —— 那是过 L2 的坎,
/// 不是算法退化。
///
/// **先核判据算的与它说的是不是一回事。** 现在把 `log` 因子除掉,
/// 判据量的才是"除掉理论增长之后还多涨了多少"。退化成平方仍然会被抓住:
/// 那时候原始涨幅是 4 倍、除掉 log 之后 3.3 倍,离这条闸很远。
const MAX_GROWTH: f64 = 1.35;

/// N 元单环:第一个等价格就是全部原子,专压打破对称的分支数。
fn macrocycle(n: usize) -> String {
    let mut s = String::from("C1");
    for _ in 1..n {
        s.push('C');
    }
    s.push('1');
    s
}

/// N 个苯环用亚甲基串起来:直径正比于 N,专压细化的轮数。
fn many_rings(n: usize) -> String {
    vec!["c1ccccc1"; n].join("C")
}

/// 带方向键的共轭多烯:专压方向判定里那趟颜色细化。
///
/// 双键两端各有两个碳取代基、局部不可区分,所以那趟细化非跑不可 ——
/// 换成 `F/C=C/F` 这种一端只有一个取代基的,判定会走捷径,又白测了。
fn directional_polyene(n: usize) -> String {
    let mut s = String::from("F");
    for _ in 0..n / 2 {
        s.push_str("/C=C");
    }
    s.push_str("/F");
    s
}

/// 让调用方指定量的是**哪一段**。
///
/// 守卫要压的那段代码若与别的开销混在一起,曲线上看到的涨幅就说不清是谁的。
/// 例如方向判定那一档:跑整个 `canonical_smiles` 会把它与打破对称的枚举揉在
/// 一起,量出来的涨幅多半来自后者的噪声,与被守的那段无关。
fn run_with(smi: &str, f: &dyn Fn(&omgkit_core::MolBuilder)) -> usize {
    let m = smiles::parse(smi).expect("语料应能解析");
    f(&m);
    m.num_atoms()
}

/// 交错测量:外层轮次、内层规模,让各档经历同样的 CPU 状态轨迹。
/// 顺序测的话最小档会系统性偏慢,足以把曲线整个翻过来。
fn measure(inputs: &[String]) -> Vec<(usize, f64)> {
    measure_with(inputs, &|m| {
        std::hint::black_box(canon::canonical_smiles(m));
    })
}

fn measure_with(inputs: &[String], f: &dyn Fn(&omgkit_core::MolBuilder)) -> Vec<(usize, f64)> {
    for _ in 0..3 {
        for smi in inputs {
            run_with(smi, f);
        }
    }
    let mut best = vec![Duration::MAX; inputs.len()];
    let mut n_atoms = vec![0usize; inputs.len()];
    for _ in 0..ROUNDS {
        for (i, smi) in inputs.iter().enumerate() {
            let t = Instant::now();
            n_atoms[i] = run_with(smi, f);
            best[i] = best[i].min(t.elapsed());
        }
    }
    best.iter()
        .enumerate()
        .map(|(i, &t)| {
            let per_atom = t.as_secs_f64() * 1e6 / n_atoms[i] as f64;
            println!("{:>6} 原子  {t:>12?}  {per_atom:>6.2} µs/原子", n_atoms[i]);
            (n_atoms[i], per_atom)
        })
        .collect()
}

fn assert_flat(rows: &[(usize, f64)], what: &str) {
    let (n_min, c_min) = rows[0];
    assert!(
        c_min * n_min as f64 > 50.0,
        "{what}:最小规模只跑了 {:.0} µs,比值无意义 —— 请调大规模",
        c_min * n_min as f64
    );
    for &(n, c) in &rows[1..] {
        // 理论上每原子耗时正比于 log n,先把这一份除掉
        let expected = (n as f64).ln() / (n_min as f64).ln();
        let growth = c / c_min / expected;
        assert!(
            growth < MAX_GROWTH,
            "{what}:每原子耗时从 {c_min:.2} µs({n_min} 原子)涨到 {c:.2} µs({n} 原子),\
             除掉 log n 的理论增长({expected:.2} 倍)之后仍多涨 {growth:.2} 倍。\
             逐档数据:{rows:?}"
        );
    }
}

#[test]
fn per_atom_cost_does_not_grow_in_a_large_macrocycle() {
    println!("形状:单个大环(压打破对称的分支数)");
    let inputs: Vec<String> = [300usize, 600, 1200]
        .iter()
        .map(|&n| macrocycle(n))
        .collect();
    assert_flat(&measure(&inputs), "大环");
}

#[test]
fn per_atom_cost_does_not_grow_in_a_long_chain() {
    println!("形状:很多个小环系串成长链(压细化的轮数)");
    let inputs: Vec<String> = [100usize, 200, 400]
        .iter()
        .map(|&n| many_rings(n))
        .collect();
    assert_flat(&measure(&inputs), "长链");
}

/// 方向判定那趟细化的增长曲线。
///
/// 少了这条,`informative_directions` 退化成平方也不会有任何测试变红 ——
/// 另外两个形状里它被提前返回跳过了。
#[test]
fn per_atom_cost_does_not_grow_with_directional_bonds() {
    let inputs: Vec<String> = [400usize, 1600, 6400]
        .iter()
        .map(|&n| directional_polyene(n))
        .collect();
    // 判据得先确认压力真的施加上了:一条方向键都没有的话这个测试什么也没测
    for smi in &inputs {
        let m = smiles::parse(smi).expect("应能解析");
        let dirs = omgkit_io::stereo::informative_directions(&m);
        assert!(
            dirs.iter().filter(|d| **d).count() > 10,
            "构造的多烯没有携带信息的方向键,压力没施加上"
        );
    }
    // 只量方向判定这一段。混进 canonical_smiles 的话,打破对称枚举的测量噪声
    // 会盖过要守的东西,曲线上的涨幅说不清是谁的。
    assert_flat(
        &measure_with(&inputs, &|m| {
            std::hint::black_box(omgkit_io::stereo::directions_for_writing(m));
        }),
        "带方向键的共轭多烯(只量方向判定)",
    );
}
