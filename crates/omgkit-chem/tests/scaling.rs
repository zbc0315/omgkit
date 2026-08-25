//! 增长曲线守卫:整条管线必须**线性于分子规模**。
//!
//! 差分测试抓不到复杂度问题 —— 结果全对,只是慢,而且在小分子上完全看不出来。
//! 所以要有一条专门盯增长曲线的测试。
//!
//! # 判据:每原子耗时不随规模上升
//!
//! 比"两个规模的耗时之比"更直接:线性实现的每原子耗时是常数,失败信息也一眼
//! 能看懂。阈值 [`MAX_GROWTH`] 是拿真实缺陷标定出来的 —— 把已知的几种平方项
//! 逐个塞回代码测过,正常约 1.0,缺陷版 1.4~3.1。改动本测试时必须重做这个标定,
//! 否则它随时可能变成一条永远绿着、什么也没守住的测试。
//!
//! # 两个把这类测试做废的坑
//!
//! **1. 规模不够。** 平方项要到几千原子才占主导;规模太小时缺陷版与正常版的
//! 差别淹没在噪声里。
//!
//! **2. 没交错测量。** 顺序测(小的先、大的后)时,最小档会系统性偏慢 ——
//! 线程起初落在能效核上、频率还没爬满 —— 足以把曲线整个翻过来,反而像是
//! "每原子耗时随规模下降"。所以外层轮次、内层规模地交错跑,让各档经历同样的
//! CPU 状态轨迹。
//!
//! # 规模和形状,两个维度都要覆盖
//!
//! | 形状 | 语料 | 压的是什么 |
//! |---|---|---|
//! | 很多个独立小环系 | N 个苯环用亚甲基串起来 | 解析、逐环系的暂存复用、逐原子的扫描 |
//! | 单个大稠合体系 | N 环线性并苯 | 分量内的算法、**环搜索的候选生成** |
//! | **单个大环** | N 元碳环 | 纯环的**快路径** —— 见下 |
//!
//! 只测一种形状会漏。已经吃过两次亏:
//!
//! - 分量内的立方项在"很多个小环系"上看不出来 —— 那里每个分量只有 6 个原子
//! - 纯环上的退化在"线性并苯"上也看不出来 —— 并苯走的是一般路径
//!
//! # 大环那一档守的**不是**候选生成
//!
//! 这里原先写着"围长大、圈秩只有 1,环搜索的候选生成必须走遍整个分量 ——
//! 这一档专门盯着它"。那是假的:`sssr::component_ring_set` 开头就有快路径
//! (分量内全部顶点度数为 2 ⇒ 分量本身就是一条简单环,直接返回),而**纯大环
//! 正好走这条快路径**。往 `horton_candidates` 里插一句打印实测:
//!
//! ```text
//! --- 大环 500 (500 原子)---        [探针] 调用 0 次
//! --- 多环系 100 (699 原子)---      [探针] 调用 0 次
//! --- 并苯 48 (194 原子)---         [探针] 调用 2 次
//! ```
//!
//! 它实际守的是**快路径还在**,而这有价值:变异实测(把快路径关掉),
//! 这一档从 0.38 µs/原子 涨到 23 / 82 / 322 µs/原子,涨幅 3.48,当场红。
//!
//! **候选生成由并苯那一档守着**,但并苯的围长只有 6。
//! **"围长大 + 走一般路径"这个组合先前没有任何测试覆盖** —— 而它正是原注释
//! 担心的那件事。现在由 [`ring_search_work_is_pinned_per_shape`] 补上,
//! 判据是**数出来的工作量**而不是墙钟。
//!
//! 那一条同时把两件事记了下来:
//!
//! - 这个形状**确实是 O(n²)**(theta 图 122 → 1922 原子,工作量 34 032 →
//!   8 677 880,255 倍 ≈ 16²),而且不是人造的极端形状 —— 大环上稠合一个苯环
//!   就是同一档(56 → 806 原子,222 倍 ≈ 14.4²)。
//! - **超线性本身还没修**,现在只是钉住了它,修好那天会当场红。

use std::time::{Duration, Instant};

use omgkit_chem::{
    assign_radicals, clean_up, kekulize, perceive_rings, ring_set, set_aromaticity,
    sssr::{ring_set_counted, SearchStats},
    update_property_cache,
};

/// 探测规模。必须大到让平方项占主导。
const RING_COUNTS: [usize; 3] = [100, 800, 1600];
const ROUNDS: usize = 20;
/// 每原子耗时的最大允许涨幅
const MAX_GROWTH: f64 = 1.25;

/// N 个苯环用亚甲基串起来:同时压到解析(N 次环闭合)、环感知
/// (N 个独立环系)和 kekulize(N 次匹配)。
fn many_rings(n: usize) -> String {
    vec!["c1ccccc1"; n].join("C")
}

fn run_pipeline(smi: &str) -> usize {
    let mut m = omgkit_io::smiles::parse(smi).expect("语料应能解析");
    clean_up(&mut m);
    update_property_cache(&mut m).expect("价键校验应通过");
    let _ = perceive_rings(&mut m);
    kekulize(&mut m).expect("应能 kekulize");
    assign_radicals(&mut m);
    std::hint::black_box(ring_set(&m));
    set_aromaticity(&mut m);
    m.num_atoms()
}

/// 逐档测量 → (原子数, 每原子微秒数)
fn measure(inputs: &[String], rounds: usize) -> Vec<(usize, f64)> {
    // 充分预热:线程会先落在能效核上再迁到性能核,频率也要爬。预热不足时
    // 最先测的那一档会系统性偏慢,足以把曲线整个翻过来。
    for _ in 0..3 {
        for smi in inputs {
            run_pipeline(smi);
        }
    }

    // **交错**测量:外层轮次、内层规模。各档因而经历同样的 CPU 状态轨迹,
    // 而不是"小的在冷机时测、大的在热机时测"。
    let mut best: Vec<Duration> = vec![Duration::MAX; inputs.len()];
    let mut n_atoms: Vec<usize> = vec![0; inputs.len()];
    for _ in 0..rounds {
        for (i, smi) in inputs.iter().enumerate() {
            let t = Instant::now();
            n_atoms[i] = run_pipeline(smi);
            // 取多轮最小值:噪声只会让测量变慢,不会变快
            best[i] = best[i].min(t.elapsed());
        }
    }

    best.iter()
        .enumerate()
        .map(|(i, &t)| {
            let per_atom_us = t.as_secs_f64() * 1e6 / n_atoms[i] as f64;
            println!(
                "{:>6} 原子  {t:>12?}  {per_atom_us:>6.2} µs/原子",
                n_atoms[i]
            );
            (n_atoms[i], per_atom_us)
        })
        .collect()
}

/// N 元单环 —— 围长大、圈秩为 1,专压环搜索的候选生成。
fn macrocycle(n: usize) -> String {
    let mut s = String::from("C1");
    for _ in 1..n {
        s.push('C');
    }
    s.push('1');
    s
}

/// 逐档检查每原子耗时不上升
fn assert_flat(rows: &[(usize, f64)], what: &str) {
    let (n_min, c_min) = rows[0];
    // 太快时测量会被调度噪声主导。门槛不必很高:计时器分辨率是纳秒级,
    // 而交错测量 + 取多轮最小值已经把调度噪声压掉了 —— 判据真正依赖的是
    // **同一轮内各档之间**的一致性,不是绝对耗时。
    assert!(
        c_min * n_min as f64 > 50.0,
        "{what}:最小规模只跑了 {:.0} µs,比值无意义 —— 请调大规模",
        c_min * n_min as f64
    );
    for &(n, c) in &rows[1..] {
        let growth = c / c_min;
        assert!(
            growth < MAX_GROWTH,
            "{what}:每原子耗时从 {c_min:.2} µs({n_min} 原子)涨到 {c:.2} µs({n} 原子),\
             涨了 {growth:.2} 倍。线性实现应基本持平;涨上去说明某处又混进了平方项 —— \
             典型形状是在按原子/按环系的循环里做了正比于整个分子的事。逐档数据:{rows:?}"
        );
    }
}

/// 线性并苯:n 个六元环稠合成**一个**大双连通分量。
///
/// 与 [`many_rings`] 是两种截然不同的压力形状:后者是很多个小分量,
/// 这里是一个大分量。分量内算法的复杂度问题只在这种形状上暴露。
fn acene(n: usize) -> String {
    assert!(n <= 99, "环闭合标号 %NN 最多到 99");
    let rn = |k: usize| {
        if k < 10 {
            k.to_string()
        } else {
            format!("%{k}")
        }
    };
    let mut s = String::from("c1ccc2");
    for k in 2..n {
        s.push_str(&format!("cc{}", rn(k + 1)));
    }
    s.push_str(&format!("cccc{}", rn(n)));
    for k in (2..n).rev() {
        s.push_str(&format!("cc{}", rn(k)));
    }
    s.push_str("cc1");
    s
}

/// **theta 图**:两个三度顶点之间三条内部不相交的长路径。
///
/// `C12` + k×C + `C(` + k×C + `1)` + k×C + `2` —— 3k+2 个原子、圈秩 2、
/// 三个 2k+2 元环,**只有 2 个原子度数为 3**,所以走不了纯环那条快路径。
///
/// 化学上这不是造出来的形状:穴醚(cryptand)、环芳烃、带一条交联的双环肽
/// 都是它。
fn theta(k: usize) -> String {
    let c = "C".repeat(k);
    format!("C12{c}C({c}1){c}2")
}

/// **大环上稠合一个苯环** —— 比 theta 更接近真实分子的那种大围长形状。
///
/// 环芳烃、带芳环的大环内酰胺都是它。同样有度数 3 的顶点、走不了纯环快路径,
/// 同样是 O(n²):56 → 806 原子(14.4 倍),工作量 6 652 → 1 478 210(**222 倍
/// ≈ 14.4²**)。
fn aryl_fused_macrocycle(n: usize) -> String {
    format!("c1ccc2c(c1){}2", "C".repeat(n))
}

/// 跑一遍环搜索,返回(原子数, 工作量)。
fn ring_work(smi: &str) -> (usize, SearchStats) {
    let mut m = omgkit_io::smiles::parse(smi).expect("语料应能解析");
    clean_up(&mut m);
    update_property_cache(&mut m).expect("价键校验应通过");
    let (_, st) = ring_set_counted(&m);
    (m.num_atoms(), st)
}

/// 环搜索的工作量:**数出来的整数**,不是墙钟。
///
/// # 为什么另起一条,而不是再加一档墙钟
///
/// 上面三条判的是"每原子**耗时**不随规模上升"。两处不够:
///
/// - 墙钟会抖。这个文件里已经记着一次:整个工作区并行跑时大环那一档慢 55%,
///   涨幅越过阈值,**判据在报噪声**。
/// - "涨幅"这个形状看不见**按比例整体变慢**的退化。
///
/// 同一条教训在 `omgkit-match/tests/scaling.rs` 上栽过一次(墙钟把一个改的是
/// 别的 crate 的提交打红了),那里已经换成数工作量并钉死绝对值。这里照办。
/// 墙钟那三条**留着**,它们守的是计数看不见的另一半:每次操作本身变贵了。
///
/// # 三种形状,三种判法
///
/// | 形状 | 判据 | 现值 |
/// |---|---|---|
/// | 单个大环 | 工作量**恒为 0** | 走快路径,一次 BFS 都不做 |
/// | 线性并苯 | 每原子工作量不上升 | 约 15.5 次 BFS/原子 |
/// | **theta 图** | 钉死绝对值 | **O(n²)**,见下 |
///
/// # theta 那一档是**钉住现状**,不是"这样就对了"
///
/// Horton 算法要对**每个**顶点做一次覆盖整个分量的 BFS,而 theta 图的围长
/// 就是 2n/3 量级,迭代加深必须一路加到那么深 —— 于是 O(n²)。实测
/// 122 → 1922 原子(16 倍)工作量 34 032 → 8 677 880(**255 倍 ≈ 16²**)。
///
/// **这个超线性本身没修**,这一条只保证它不再悄悄变坏,而且修好的那天会
/// **当场红**(数变小了),逼着改的人把新数写进来。教科书的修法是
/// **抑制二度顶点**(把每条二度链缩成一条带权边,再在缩图上跑带权 Horton),
/// 那要把 BFS 换成 Dijkstra、处理重边、再把环展开回去 —— 是单独一块的活。
///
/// 顺带量过:**真实语料一次都够不着这一档**。`large.smi` 的 8830 个分子
/// 合计只做 242 637 次 BFS,最坏的一个 1 252 次(90 原子的双环缩肽);
/// `hard.smi` 68 个分子合计 2 896、最坏 864;`bridged.smi` 25 个合计 2 517。
/// 都不到 theta k=40 一个分子(34 032)的零头 —— 也就是说这个形状**在语料里
/// 是缺的**,这条测试补的正是它。
#[test]
fn ring_search_work_is_pinned_per_shape() {
    // 一、单个大环:快路径,一次 BFS 都不做。
    //    这一条比上面那条墙钟版严得多 —— 快路径一旦失效,数就不是 0 了,
    //    而墙钟版要等到慢出 1.25 倍才红。
    for n in [100usize, 400, 1600] {
        let (atoms, st) = ring_work(&macrocycle(n));
        assert_eq!(
            (st.bfs_visits, st.edge_tests, st.path_steps),
            (0, 0, 0),
            "{atoms} 元大环走了一般路径 —— `component_ring_set` 开头那条\
             '全部顶点度数为 2 ⇒ 分量本身就是一条环'的快路径失效了"
        );
    }

    // 二、线性并苯:每原子工作量不上升。
    let mut rows: Vec<(usize, f64)> = Vec::new();
    let mut path_rows: Vec<(usize, f64)> = Vec::new();
    for n in [12usize, 24, 48] {
        let (atoms, st) = ring_work(&acene(n));
        #[allow(clippy::cast_precision_loss)]
        let per = st.bfs_visits as f64 / atoms as f64;
        #[allow(clippy::cast_precision_loss)]
        let per_path = st.path_steps as f64 / atoms as f64;
        println!(
            "并苯 {n:>3}:{atoms:>4} 原子  bfs={:<8} {per:>6.2} 次/原子  \
             回溯 {:<8} {per_path:>6.2} 步/原子",
            st.bfs_visits, st.path_steps
        );
        rows.push((atoms, per));
        path_rows.push((atoms, per_path));
    }
    let base = rows[0].1;
    for &(atoms, per) in &rows[1..] {
        assert!(
            per / base < MAX_GROWTH,
            "稠合体系:每原子 BFS 次数从 {base:.2}({} 原子)涨到 {per:.2}({atoms} 原子)",
            rows[0].0
        );
    }
    // 回溯量也要判 —— 两条剪枝失效时 `bfs_visits` 与 `edge_tests` 纹丝不动,
    // 只有这一项会涨(见 `SearchStats::path_steps`)。
    let base_path = path_rows[0].1;
    for &(atoms, per) in &path_rows[1..] {
        assert!(
            per / base_path < MAX_GROWTH,
            "稠合体系:每原子回溯步数从 {base_path:.2}({} 原子)涨到 {per:.2}({atoms} 原子)\
             —— 平衡剪枝或首步剪枝失效了",
            path_rows[0].0
        );
    }

    // 三、大围长 + 走一般路径:钉死绝对值。**变小也红** —— 见本条文档。
    const THETA: &[(usize, usize, u64, u64, u64)] = &[
        // (k, 原子数, bfs_visits, edge_tests, path_steps)
        (40, 122, 34_032, 35_010, 20_252),
        (80, 242, 136_034, 138_240, 78_892),
        (160, 482, 543_412, 548_310, 311_372),
    ];
    for &(k, want_atoms, want_bfs, want_edge, want_path) in THETA {
        let (atoms, st) = ring_work(&theta(k));
        assert_eq!(atoms, want_atoms, "theta k={k} 的原子数变了 —— 生成器改了?");
        assert_eq!(
            (st.bfs_visits, st.edge_tests, st.path_steps),
            (want_bfs, want_edge, want_path),
            "theta k={k}({atoms} 原子)的环搜索工作量变了。\
             **变大**是回归;**变小**说明超线性被治了 —— 那是好事,\
             把新的数写进来,并在提交信息里说明治的是什么。"
        );
    }

    // 四、苯稠大环 —— 同一个超线性,但形状是真实分子里有的那种。
    //     theta 图容易被读成"人造的极端形状";这一档说明不是。
    const FUSED: &[(usize, usize, u64, u64, u64)] = &[
        // (链长, 原子数, bfs_visits, edge_tests, path_steps)
        (50, 56, 6_652, 6_962, 3_095),
        (200, 206, 95_526, 97_060, 41_879),
    ];
    for &(n, want_atoms, want_bfs, want_edge, want_path) in FUSED {
        let (atoms, st) = ring_work(&aryl_fused_macrocycle(n));
        assert_eq!(atoms, want_atoms, "苯稠大环 {n} 的原子数变了");
        assert_eq!(
            (st.bfs_visits, st.edge_tests, st.path_steps),
            (want_bfs, want_edge, want_path),
            "苯稠大环 {n}({atoms} 原子)的环搜索工作量变了 —— 同上"
        );
    }
}

#[test]
fn per_atom_cost_does_not_grow_with_molecule_size() {
    println!("形状:很多个独立小环系");
    let inputs: Vec<String> = RING_COUNTS.iter().map(|&n| many_rings(n)).collect();
    let rows = measure(&inputs, ROUNDS);
    assert_flat(&rows, "多环系");
}

/// 第三种压力形状:**单个大环** —— 盯的是 `sssr` 里那条纯环快路径。
///
/// (先前这里写的是"盯环搜索的候选生成",而纯大环恰恰**走不到**候选生成 ——
///  实测 `horton_candidates` 在这一档一次都没被调用。详见文件头。)
///
/// 变异实测:把 `component_ring_set` 开头那条快路径关掉,
/// 这一档从 0.38 µs/原子 变成 23.5 / 81.8 / 322 µs/原子,涨幅 3.48,当场红。
#[test]
fn per_atom_cost_does_not_grow_in_a_large_macrocycle() {
    println!("形状:单个大环");
    let inputs: Vec<String> = [500usize, 1000, 2000]
        .iter()
        .map(|&n| macrocycle(n))
        .collect();
    // **轮数必须和另外两条一样是 `ROUNDS`。** 先前这里写的是 3,而
    // [`measure`] 的"取多轮最小值"要靠轮数够多才挡得住噪声 ——
    // 3 轮时,`cargo test --workspace` 并行争用一下就能把某一档的三轮全抬高。
    //
    // 实测:单独跑这条测试八次,每原子耗时是 0.38 / 0.38 / 0.39 µs,
    // 涨幅 **1.02**,稳得很;而在整个工作区一起跑的那一次,1000 原子那一档
    // 测到 588 µs(单独跑时约 380 µs,慢 55%),涨幅 1.55,越过 `MAX_GROWTH = 1.25`
    // 直接落进文档标定的"缺陷版 1.4~3.1"那一带 —— **判据在报噪声,不是在报缺陷**。
    //
    // 代价可以忽略:20 轮 × 三档 ≈ 27 ms。
    let rows = measure(&inputs, ROUNDS);
    assert_flat(&rows, "大环");
}

/// 另一种压力形状:**单个大稠合体系**。
///
/// 分量内算法(相关环搜索、匹配搜索)的复杂度只在这种形状上暴露 ——
/// 上面那条"很多个小环系"里每个分量只有 6 个原子,看不出来。
#[test]
fn per_atom_cost_does_not_grow_in_a_single_large_fused_system() {
    println!("形状:单个大稠合体系(线性并苯)");
    // 环闭合标号 %NN 最多到 99,故并苯环数上限 99
    let inputs: Vec<String> = [48usize, 72, 96].iter().map(|&n| acene(n)).collect();
    let rows = measure(&inputs, ROUNDS);
    assert_flat(&rows, "稠合体系");
}
