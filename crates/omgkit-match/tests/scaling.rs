//! 子结构匹配的增长曲线守卫。
//!
//! 回溯搜索是这个项目里**最容易在复杂度上出事**的一块:它天然是指数的,
//! 全靠剪枝把它压住。剪枝一旦失效,结果照样全对,只是慢 —— 而且在小分子上
//! 完全看不出来,差分测试也抓不到。
//!
//! # 三种压力形状,压的是三处不同的东西
//!
//! | 形状 | 压的是什么 | 剪枝失效时的样子 |
//! |---|---|---|
//! | 分子变大、模式不变 | **候选生成**。每层的候选必须来自已映射原子的邻居 | 每层都扫全分子 → 平方 |
//! | 模式变长、分子不变 | **定序**。新加的查询原子必须与已映射的相连 | 断开的层要全分子回溯 → 指数 |
//! | 对称分子 + 早停 | **早停**。`max_matches` 要真的截断搜索 | 先枚举完再截断 → 阶乘 |
//!
//! 前两种互相看不见:分子大而模式只有两个原子时,定序问题显不出来;
//! 模式长而分子只有十几个原子时,候选生成的问题也显不出来。
//!
//! # 判据与阈值
//!
//! 判据是**每原子耗时不随规模上升**,与 `omgkit-chem/tests/scaling.rs` 同一套。
//! 那边的文件头写了两个把这类测试做废的坑(规模不够、没交错测量),这里同样
//! 适用。阈值拿真实缺陷标定过:把候选生成改成全分子扫描,大分子那一档立刻
//! 报出 3 倍以上的涨幅。

use std::time::{Duration, Instant};

use omgkit_chem::sanitize;
use omgkit_io::{smarts, smiles};
use omgkit_match::{substructure_matches, MatchOptions, MolProps};

const ROUNDS: usize = 5;
/// 每原子耗时的最大允许涨幅
const MAX_GROWTH: f64 = 1.6;

/// N 个苯环用亚甲基串起来 —— 分子规模这一维。
fn many_rings(n: usize) -> String {
    vec!["c1ccccc1"; n].join("C")
}

/// N 元碳链 —— 模式规模这一维用它当靶子,匹配位置多但每处都很浅。
fn chain(n: usize) -> String {
    "C".repeat(n)
}

/// 一次测量:(规模, 每单位耗时微秒)
fn measure<F>(sizes: &[usize], mut run: F) -> Vec<(usize, f64)>
where
    F: FnMut(usize) -> usize,
{
    // 预热:线程会先落在能效核上再迁到性能核,频率也要爬
    for _ in 0..2 {
        for &n in sizes {
            run(n);
        }
    }
    let mut best = vec![Duration::MAX; sizes.len()];
    let mut units = vec![0usize; sizes.len()];
    // 交错测量:外层轮次、内层规模,让各档经历同样的 CPU 状态轨迹
    for _ in 0..ROUNDS {
        for (i, &n) in sizes.iter().enumerate() {
            let t = Instant::now();
            units[i] = run(n);
            best[i] = best[i].min(t.elapsed());
        }
    }
    best.iter()
        .enumerate()
        .map(|(i, &t)| {
            let per = t.as_secs_f64() * 1e6 / units[i] as f64;
            println!("{:>6} 单位  {t:>12?}  {per:>7.3} µs/单位", units[i]);
            (units[i], per)
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
        let growth = c / c_min;
        assert!(
            growth < MAX_GROWTH,
            "{what}:每单位耗时从 {c_min:.3} µs({n_min} 单位)涨到 {c:.3} µs({n} 单位),\
             涨了 {growth:.2} 倍。线性实现应基本持平;涨上去说明剪枝失效了。逐档:{rows:?}"
        );
    }
}

/// 分子这一维:模式固定,分子越来越大。
///
/// 压的是候选生成 —— 每层的候选必须来自已映射原子的邻居,而不是全分子。
#[test]
fn cost_per_atom_does_not_grow_with_molecule_size() {
    println!("形状:模式固定,分子变大");
    let q = smarts::parse("c1ccccc1").expect("模式应能解析");
    let opts = MatchOptions::default();
    let rows = measure(&[40usize, 80, 160], |n| {
        let mut m = smiles::parse(&many_rings(n)).expect("语料应能解析");
        sanitize(&mut m).expect("应能净化");
        let props = MolProps::compute(&m);
        let hits = substructure_matches(&q, &m, &props, opts);
        std::hint::black_box(&hits);
        m.num_atoms()
    });
    assert_flat(&rows, "分子规模");
}

/// 模式这一维:分子固定,模式越来越长。
///
/// 压的是定序 —— 新加的查询原子必须与已映射的相连,否则每层都要全分子回溯。
#[test]
fn cost_per_query_atom_does_not_grow_with_pattern_size() {
    println!("形状:分子固定,模式变长");
    let mut m = smiles::parse(&chain(300)).expect("语料应能解析");
    sanitize(&mut m).expect("应能净化");
    let props = MolProps::compute(&m);
    let opts = MatchOptions {
        max_matches: 200,
        uniquify: true,
        use_chirality: true,
    };
    let rows = measure(&[8usize, 16, 32], |k| {
        let q = smarts::parse(&chain(k)).expect("模式应能解析");
        let hits = substructure_matches(&q, &m, &props, opts);
        std::hint::black_box(&hits);
        k
    });
    assert_flat(&rows, "模式规模");
}

/// 早停必须真的截断搜索,而不是先枚举完再截断。
///
/// # 判据要用**枚举不完**的用例
///
/// 用 `cc` 在 640 个苯环上测是不够的:不去重也只有约 7680 个匹配,全枚举
/// 只要几微秒,早停失效根本压不出来 —— 去掉早停测试照样绿。
///
/// 换成一条 12 个通配原子的链:每一步都能走向任意邻居,匹配数按分支因子的
/// 12 次方增长,全枚举根本跑不完。早停生效时**微秒级**返回,失效时会挂住。
/// 所以这里断言的是**绝对耗时上界**,不是相对增长 —— 相对比较在"跑不完"
/// 面前没有意义。
#[test]
fn max_matches_truncates_the_search() {
    // 通配原子 + 任意键:分支因子最大化
    let q = smarts::parse("*~*~*~*~*~*~*~*~*~*~*~*").expect("模式应能解析");
    let opts = MatchOptions {
        max_matches: 1,
        uniquify: false,
        use_chirality: true,
    };
    let mut m = smiles::parse(&many_rings(60)).expect("语料应能解析");
    sanitize(&mut m).expect("应能净化");
    let props = MolProps::compute(&m);

    let t = Instant::now();
    let hits = substructure_matches(&q, &m, &props, opts);
    let el = t.elapsed();

    assert_eq!(hits.len(), 1, "早停应当只返回一个匹配");
    assert!(
        el < Duration::from_millis(50),
        "早停没有截断搜索:{} 原子上找第一个匹配花了 {el:?}。\n\
         这条模式的全枚举是指数的,耗时一旦上到毫秒级就说明它在往下枚举。",
        m.num_atoms()
    );
    println!("{:>6} 原子  找到首个匹配 {el:?}", m.num_atoms());
}

/// 定序要按**候选数**挑起点,不是按度数。
///
/// 模式末尾挂一个稀有元素时,从它起头候选只有 1 个,从碳起头候选是几百个。
/// 两者的差别在链状模式上完全体现不出来(度数都一样),必须用这种形状测。
#[test]
fn rare_atoms_are_anchored_first() {
    // 500 个碳里只有一个溴
    let mut m = smiles::parse(&format!("{}Br", chain(500))).expect("语料应能解析");
    sanitize(&mut m).expect("应能净化");
    let props = MolProps::compute(&m);
    let q = smarts::parse("CCCCCCCCBr").expect("模式应能解析");

    let t = Instant::now();
    let hits = substructure_matches(&q, &m, &props, MatchOptions::default());
    let el = t.elapsed();

    assert_eq!(hits.len(), 1);
    assert!(
        el < Duration::from_micros(150),
        "从稀有原子起头的话这一步几乎不用回溯,实际花了 {el:?} —— \n\
         说明起点挑的是碳,于是几百个碳都要各试一遍。\n\
         阈值卡在两者之间:定序完全不看候选数时实测 1.96ms,看候选数时实测 31µs"
    );
    println!("稀有起点:{el:?}");
}
