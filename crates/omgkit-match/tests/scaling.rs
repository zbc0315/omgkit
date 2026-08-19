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
        el < Duration::from_micros(600),
        "从稀有原子起头的话这一步几乎不用回溯,实际花了 {el:?} —— \n\
         说明起点挑的是碳,于是几百个碳都要各试一遍。\n\
         阈值卡在两者之间:定序完全不看候选数时实测 1.96ms,看候选数时实测 31µs"
    );
    // 阈值放在 600µs 而不是 150µs:健康值 31µs、退化值 1.96ms,**两者差 63 倍**,
    // 阈值在这中间有很大余地。原先的 150µs 只比健康值高 4.8 倍 —— 那在机器有负载时
    // 不够:实测跟着别的闸门一起跑时量到 226µs,单独跑 5 次全是绿的。
    // 600µs 距健康值 19 倍、距退化值仍有 3.3 倍,两头都留得住。
    //
    // (与 omgkit-io 的 canon_scaling 是同一类问题:墙钟阈值必须放在**实测的
    //  健康值与退化值之间**,而不是紧贴健康值。)
    println!("稀有起点:{el:?}");
}

// ---------------------------------------------------------------------------
// 副产物收口的增长曲线
// ---------------------------------------------------------------------------

/// 收口必须**线性于片段规模**。
///
/// # 为什么这里特别容易出事
///
/// 收口的一批操作天然是"按位点"的:摘哪个原子的氢、哪两处空价配对、电荷落在
/// 哪个原子。而位点表是**每个片段原子一条**,不是每处空价一条 —— 于是"遍历
/// 位点"这件事的代价正比于整个片段,而不是正比于真正欠着价的那几个原子。
/// 在按键或按位点的循环里再做一次这样的遍历,就是本仓库反复警告的那个形状。
///
/// 差分测试抓不到这类问题:结果全对,只是慢,而且真实的离去片段都很小,
/// 在语料上完全看不出来。
///
/// # 阈值拿真实缺陷标定过
///
/// 判据是**每原子耗时不随规模上升**,与本文件其余部分同一套。实测(release,
/// 交错测量、取最小值):
///
/// | | n=800 | 1600 | 3200 | 6400 | 涨幅 |
/// |---|---|---|---|---|---|
/// | 现在 | 230 | 219 | 208 | 204 ns/原子 | **1.13×** |
/// | 塞回一个提不出去的平方项 | 238 | 337 | 530 | 931 | **3.92×** |
///
/// 阈值卡在两者之间。**标定时踩过一个坑值得记下**:第一次把平方项塞回去时曲线
/// 一点没动 —— 因为塞的那个形状(`if sites[i].opens == 0 { continue }`)是循环
/// 不变量,被编译器提到内层循环外了。用"塞回去看它红不红"标定阈值时,得确认
/// 塞的东西**真的**在跑,否则标出来的阈值毫无意义。
#[test]
fn closing_byproducts_is_linear_in_fragment_size() {
    use omgkit_match::{byproduct, run_reactants};

    // 阈值用模块级那个 `MAX_GROWTH`(同样是 1.6)—— 改用 `assert_flat` 之后
    // 本函数不再自己判,局部那份就成了死代码。
    let sizes = [800usize, 1600, 3200, 6400];

    // 长链酯:模板删掉 OCH2C,整条链失去落脚点成为一个很大的片段,
    // 而且要摘一个氢、成一根键 —— 摘氢与配对两条路都走得到。
    let prep = |n: usize| {
        let rxn = smarts::parse_reaction("[C:1](=[O:2])[O:3][CH2]C>>[C:1](=[O:2])[OH:3]")
            .expect("模板应能解析");
        let mut m = smiles::parse(&format!("CC(=O)O{}", "C".repeat(n))).expect("底物应能解析");
        sanitize(&mut m).expect("应能净化");
        let props = MolProps::compute(&m);
        let outs = run_reactants(&rxn, &[(m.clone(), props)], 1, false);
        assert!(!outs.is_empty(), "n={n} 应当出产物");
        (m, outs)
    };
    let cases: Vec<_> = sizes.iter().map(|&n| prep(n)).collect();

    // 预热:线程会先落在能效核上再迁到性能核,频率也要爬
    for _ in 0..2 {
        for (m, o) in &cases {
            let _ = byproduct::reconstruct(std::slice::from_ref(m), &o[0]);
        }
    }
    // 交错测量:外层轮次、内层规模,让各档经历同样的 CPU 状态轨迹
    let mut best = vec![Duration::MAX; sizes.len()];
    for _ in 0..ROUNDS {
        for (i, (m, o)) in cases.iter().enumerate() {
            let t = Instant::now();
            let by = byproduct::reconstruct(std::slice::from_ref(m), &o[0]);
            best[i] = best[i].min(t.elapsed());
            assert!(by.verdict.is_closed(), "n={} 应当收得了口", sizes[i]);
        }
    }

    // **比的是"涨没涨",不是"散不散"。**
    //
    // 头一版写的是 `hi / lo`(全档的最大比最小),那测的是**离散度**不是**增长**。
    // CI 上真栽过一次:实测 800→6400 原子的每原子耗时是 449 / 435 / 365 / 276 ns,
    // **单调往下走**(小规模那一档被固定开销摊薄了),明明比线性还好,
    // 而 `hi/lo = 1.62` 照样把断言顶红,报的还是"有一处正比于整个片段的操作"——
    // 数据说的正相反。同一个文件里的 `assert_flat` 一直是对的:拿每一档比**最小
    // 那一档**,只看增长方向。这里改成用它,顺带白拿那条"最小规模够不够大"的闸。
    #[allow(clippy::cast_precision_loss)]
    let rows: Vec<(usize, f64)> = sizes
        .iter()
        .zip(&best)
        .map(|(&n, d)| (n, d.as_nanos() as f64 / n as f64 / 1000.0))
        .collect();
    for (n, c) in &rows {
        println!("{n:>6} 原子  每原子 {:.0} ns", c * 1000.0);
    }
    assert_flat(&rows, "收口");
}
