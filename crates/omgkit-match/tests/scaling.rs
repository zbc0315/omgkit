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
use omgkit_match::{substructure_matches, substructure_matches_counted, MatchOptions, MolProps};

/// 取多轮最小值的轮数。
///
/// **三个时序判据文件对齐成 20**(`omgkit-chem/tests/scaling.rs` 一直是 20,
/// `omgkit-io/tests/canon_scaling.rs` 从 8 抬上来)。这里先前是 5,是三者里最低的,
/// 也是最爱抽风的:**本机在跑完一轮全量测试之后**单独跑 10 次红 2 次
/// (独立审核在空闲机器上 42 次零红 —— 抽风依赖机器状态,这个数只在
/// 那个条件下成立,别当成稳定的发生率)。签名与另外两处一模一样 ——
/// **中间那档偏高而最大档回落**(279→1.33、559→**2.39**、1119→1.31 µs/单位),
/// 不是平方项。20 轮时 12 次零红,耗时仍是 0.07 秒。
///
/// 这与"失败了就重试"是两回事:重试是**只会让判据变绿**的单向过滤器
/// (那条路走过一次,量出来边缘缺陷的检出率从 47% 掉到 20%,已回退)。
/// 加轮数不动判据规则,只让 min-of-N 更贴近真最小值 —— 对健康值与缺陷值
/// 同样收紧,不牺牲检出力。
const ROUNDS: usize = 20;
/// 每原子**耗时**的最大允许涨幅。
///
/// 只剩 `closing_byproducts_is_linear_in_fragment_size` 在用它 —— 那条路上
/// 还没有确定性的工作量计数。匹配那两条已经换成数 `candidate_tests`,
/// 理由见 [`assert_pinned`]。
const MAX_GROWTH: f64 = 1.6;

/// 每单位**工作量**的最大允许涨幅。
///
/// 比耗时那个紧,因为它不抖:同一份输入永远同一个数。实测健康值 1.000
/// (分子那一维**恰好**是 1.000)与 1.105(模式那一维),取 1.25。
const MAX_COUNT_GROWTH: f64 = 1.25;

/// 每单位耗时的**粗**闸,常数因子崩了才拦。
///
/// 细的那道交给工作量计数。留着它是因为计数看不见"每次候选测试本身变贵"
/// (比如 `atom_feasible` 里多了一遍 O(n) 扫描)—— 那种回归不改变计数,
/// 只改变耗时。
const MAX_TIME_GROWTH_COARSE: f64 = 3.0;

/// N 个苯环用亚甲基串起来 —— 分子规模这一维。
fn many_rings(n: usize) -> String {
    vec!["c1ccccc1"; n].join("C")
}

/// N 元碳链 —— 模式规模这一维用它当靶子,匹配位置多但每处都很浅。
fn chain(n: usize) -> String {
    "C".repeat(n)
}

/// N 个立方烷笼首尾相连 —— **专给早停那条判据用的高分支度语料**。
///
/// 立方烷的八个碳全是三度,而且每个邻居都还能继续往下走(不是死路)。
/// 实测 60 个单元 = 480 原子、单连通、度数只有 3 和 4,平均度 **3.25**。
///
/// 对比现成的那几种形状(**全部是 16 个通配原子**的自避路径全枚举,release 实测):
///
/// | 语料 | 原子 | 平均度 | 全枚举匹配数 | 全枚举耗时 |
/// |---|---|---|---|---|
/// | 主链每个碳挂一个甲基 | 401 | 2.00 | 1 492 | 0.6 ms |
/// | `many_rings(60)`(苯环用亚甲基串起来) | 419 | 2.28 | 13 392 | 6.3 ms |
/// | 线性并苯 | 386 | 2.49 | 103 456 | 13 ms |
/// | **`cubanes(60)`** | 480 | **3.25** | **7 297 154** | **2.2 s** |
///
/// 平均度从 2.00 到 3.25,自避路径数差 **4891 倍(3.7 个数量级)** —— 路径数按
/// (度数−1) 的幂增长,而链状/环状语料的分支因子接近 1。
///
/// (这张表的并苯与挂甲基两行原先填的是 **12** 个通配原子的数,而表头写着 16,
///  独立审核逐条重量才发现。**换了模式长度,整张表都要重量** —— 表是用来论证
///  选哪个语料的,论据错了三倍,结论就是蒙对的。)
fn cubanes(n: usize) -> String {
    "C12C3C4C1C5C4C3C25".repeat(n)
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

/// 一次计数测量:(规模, 工作量, 每单位工作量, 每单位耗时微秒)。
///
/// 工作量取 [`SearchStats::candidate_tests`],**整数、确定、与机器无关** ——
/// 实测 debug 与 release 两档逐位相同。
fn measure_work<F>(sizes: &[usize], mut run: F) -> Vec<(usize, u64, f64, f64)>
where
    F: FnMut(usize) -> (usize, u64),
{
    for _ in 0..2 {
        for &n in sizes {
            run(n);
        }
    }
    let mut best = vec![Duration::MAX; sizes.len()];
    let mut units = vec![0usize; sizes.len()];
    let mut work = vec![0u64; sizes.len()];
    for _ in 0..ROUNDS {
        for (i, &n) in sizes.iter().enumerate() {
            let t = Instant::now();
            let (u, w) = run(n);
            best[i] = best[i].min(t.elapsed());
            units[i] = u;
            work[i] = w;
        }
    }
    (0..sizes.len())
        .map(|i| {
            let per_w = work[i] as f64 / units[i] as f64;
            let per_t = best[i].as_secs_f64() * 1e6 / units[i] as f64;
            println!(
                "{:>6} 单位  候选测试 {:>9}  每单位 {per_w:>9.3}  ({:>10?}, {per_t:>7.3} µs/单位)",
                units[i], work[i], best[i]
            );
            (units[i], work[i], per_w, per_t)
        })
        .collect()
}

/// **钉死工作量,而不是量耗时的涨幅。**
///
/// # 为什么换掉墙钟
///
/// 这两条判据先前判的是"每单位**耗时**的涨幅 < 1.6"。8/16/32 元模式那一档实测
/// 213 µs / 440 µs / 1374 µs —— 微秒级的数放在共享 CI 机器上,2026-08-25
/// 把一个**改的是 `omgkit-conf`** 的提交打红了(`omgkit-match` 根本不依赖它),
/// 涨幅 1.61 差一点点越过 1.6。
///
/// # 为什么"比值"这个形状本身也不够
///
/// 比值只看得见**随规模变化**的退化,看不见按比例整体变差的。实测:
///
/// | 变异 | 分子那一维的比值 | 模式那一维的比值 | 绝对工作量 |
/// |---|---|---|---|
/// | 候选生成退化成全分子扫描 | 1.0 → **4.0** 抓住 | 1.10 → 1.14 **漏** | 5377 → **777 128**(145 倍) |
/// | 起点不按候选数挑(稀有锚点那条) | — | 1.0 → 1.0 **漏** | 316 → **9 181**(29 倍) |
///
/// 两个"漏"都是同一个原因:退化在每一档上按同样倍数发生,比值当然平。
/// 所以主判据是**钉死绝对值**;比值留着,负责另一半(只在大规模上冒头的退化)。
///
/// # 还发现一件事:那条"模式变长"的判据不守它自称守的东西
///
/// 它的文档说压的是**定序**。而两种定序变异(整个倒过来、起点按下标挑)在
/// `chain(k)` 配 `chain(300)` 上把工作量从 5377 只挪到 5351 —— 链上每个碳
/// 都等价,任何顺序都是连通的,**这个语料对定序没有压力**。
/// 所以另加了一条 `Br` 锚点的形状,那一条上同样的变异是 29–99 倍。
fn assert_pinned(rows: &[(usize, u64, f64, f64)], want: &[u64], what: &str) {
    let got: Vec<u64> = rows.iter().map(|r| r.1).collect();
    assert_eq!(
        got, want,
        "{what}:搜索的工作量变了。\n\
         这个数是确定的(debug 与 release 逐位相同),变了就是搜索行为变了。\n\
         **变小**说明剪枝更好了 —— 把新数填进来,并在提交信息里说清为什么;\n\
         **变大**说明剪枝退了,先查再改数。"
    );
    let (n_min, _, w_min, t_min) = rows[0];
    for &(n, _, w, t) in &rows[1..] {
        assert!(
            w / w_min < MAX_COUNT_GROWTH,
            "{what}:每单位工作量从 {w_min:.3}({n_min} 单位)涨到 {w:.3}({n} 单位),\
             涨了 {:.2} 倍",
            w / w_min
        );
        assert!(
            t / t_min < MAX_TIME_GROWTH_COARSE,
            "{what}:每单位耗时涨了 {:.2} 倍(粗闸 {MAX_TIME_GROWTH_COARSE})——\
             工作量没涨而耗时涨了,说明每次候选测试本身变贵了",
            t / t_min
        );
    }
}

/// 分子这一维:模式固定,分子越来越大。
///
/// 压的是候选生成 —— 每层的候选必须来自已映射原子的邻居,而不是全分子。
/// 变异实测:改成全分子扫描,工作量 5301 → 602 919(114 倍),比值 1.0 → 4.0。
#[test]
fn work_per_atom_does_not_grow_with_molecule_size() {
    println!("形状:模式固定,分子变大");
    let q = smarts::parse("c1ccccc1").expect("模式应能解析");
    let opts = MatchOptions::default();
    let rows = measure_work(&[40usize, 80, 160], |n| {
        let mut m = smiles::parse(&many_rings(n)).expect("语料应能解析");
        sanitize(&mut m).expect("应能净化");
        let props = MolProps::compute(&m);
        let (hits, st) = substructure_matches_counted(&q, &m, &props, opts);
        std::hint::black_box(&hits);
        (m.num_atoms(), st.candidate_tests)
    });
    assert_pinned(&rows, &[5301, 10621, 21261], "分子规模");
}

/// 模式这一维:分子固定,模式越来越长。
///
/// 压的是候选生成在**深**搜索下的表现。**它压不到定序** —— 见
/// [`assert_pinned`] 的最后一节,定序那件事由
/// [`work_does_not_blow_up_when_the_anchor_is_rare`] 守。
#[test]
fn work_per_query_atom_does_not_grow_with_pattern_size() {
    println!("形状:分子固定,模式变长");
    let mut m = smiles::parse(&chain(300)).expect("语料应能解析");
    sanitize(&mut m).expect("应能净化");
    let props = MolProps::compute(&m);
    let opts = MatchOptions {
        max_matches: 200,
        uniquify: true,
        use_chirality: true,
    };
    let rows = measure_work(&[8usize, 16, 32], |k| {
        let q = smarts::parse(&chain(k)).expect("模式应能解析");
        let (hits, st) = substructure_matches_counted(&q, &m, &props, opts);
        std::hint::black_box(&hits);
        (k, st.candidate_tests)
    });
    assert_pinned(&rows, &[5377, 11633, 23761], "模式规模");
}

/// 定序这一维:模式一端挂着**稀有原子**,分子里只有一个。
///
/// 压的是"起点挑候选最少的那个"。挑错起点时整条链的每个碳都要各试一遍 ——
/// `search_order` 的文档举的正是 `CCCCCCCCBr` 这个例子,而先前**没有任何判据
/// 压得到它**:`chain(k)` 配 `chain(300)` 上,把起点改成按下标挑,工作量
/// 5377 → 5351,纹丝不动。
///
/// 换成这个形状之后,同一个变异是 316 → 9 181 / 332 → 18 405 / 364 → 36 085
/// (**29 / 55 / 99 倍**)。
#[test]
fn work_does_not_blow_up_when_the_anchor_is_rare() {
    println!("形状:稀有锚点,模式变长");
    let mut m = smiles::parse(&format!("Br{}", chain(300))).expect("语料应能解析");
    sanitize(&mut m).expect("应能净化");
    let props = MolProps::compute(&m);
    let opts = MatchOptions {
        max_matches: 200,
        uniquify: true,
        use_chirality: true,
    };
    let rows = measure_work(&[8usize, 16, 32], |k| {
        let q = smarts::parse(&format!("{}Br", chain(k))).expect("模式应能解析");
        let (hits, st) = substructure_matches_counted(&q, &m, &props, opts);
        assert_eq!(hits.len(), 1, "{k} 元模式该正好命中一次");
        (k, st.candidate_tests)
    });
    // 这一维的每单位工作量是**下降**的(锚点钉死之后,模式越长摊得越薄),
    // 所以只钉绝对值,不判涨幅。
    assert_eq!(
        rows.iter().map(|r| r.1).collect::<Vec<_>>(),
        vec![316u64, 332, 364],
        "稀有锚点:工作量变了 —— 变大说明起点挑错了"
    );
}

/// 早停必须真的截断搜索,而不是先枚举完再截断。
///
/// # 判据要用**枚举不完**的用例,而这一条先前不是
///
/// 头一版用 `cc` 在 640 个苯环上测,不去重也只有约 7680 个匹配,全枚举几微秒
/// 就完了 —— 去掉早停测试照样绿。于是换成 12 个通配原子的链,并在文档里写
/// "匹配数按分支因子的 12 次方增长,全枚举根本跑不完"。
///
/// **那句话是错的,而且这条判据因此一直是空的。** `many_rings(60)` 是苯环用
/// 亚甲基串成的**链**,绝大多数原子度数为 2,分支因子接近 1 —— 12 个通配原子
/// 的自避路径统共只有 7496 条。变异实测(把 `extend` 里两处
/// `out.len() >= opts.max_matches` 关掉、末尾改成 `out.truncate(max_matches)`,
/// 正是本文档点名的"先枚举完再截断"):
///
/// ```text
/// 419 原子  找到首个匹配 4.608042ms      ← 阈值 50 ms
/// test result: ok. 1 passed             ← 照样绿
/// ```
///
/// # 换成高分支度语料 + 更长的模式
///
/// 现在用 [`cubanes`] 这个语料(那里有几种形状的对比表)和 16 个通配原子。
/// release 实测:
///
/// | | 耗时 | 相对 50 ms 阈值 |
/// |---|---|---|
/// | 早停生效(健康) | **50–75 µs** | 阈值是它的 **700~1000 倍** |
/// | **同一个变异**(先枚举完再截断) | **2.2 s**(730 万个匹配) | 阈值的 **1/44** |
///
/// 健康值写成区间是有意的:多次实测 release 47–63 µs、`[profile.test]` 67–72 µs,
/// 机器安静时还能低到 21 µs。**墙钟值别写成单点** —— 写死一个数,下一个人在
/// 别的机器上量到不一样就不知道该信哪个。
///
/// 而**同一个变异**在换语料之前只跑出 4.6 ms —— 判据纹丝不动。
///
/// 两头都留得很宽,所以阈值不用动 —— 这正是墙钟阈值该有的样子:
/// **放在实测的健康值与退化值之间**,而不是紧贴任何一边。
///
/// 断言的是**绝对耗时上界**,不是相对增长 —— 相对比较在"跑不完"面前没有意义。
///
/// # 两条路各测一遍:`uniquify` 开与关
///
/// 早停在这两条路上是**分别**判断的,而判据先前只跑 `uniquify: false`。
/// 独立审核做过这个变异:把早停判断改成只在 `!opts.uniquify` 时生效 ——
/// 判据照样绿,而 `uniquify: true` 那条路从 54 µs 变成 **7.2 秒**。
/// `omgkit-py` 把这两个参数直接交给用户(`Query.match(mol, uniquify=True,
/// max_matches=1)`),那条正是没人守的路。所以下面两种都测。
#[test]
fn max_matches_truncates_the_search() {
    // 通配原子 + 任意键:分支因子最大化
    let q = smarts::parse("*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*").expect("模式应能解析");
    // **语料必须是高分支度的**,见函数文档:链状语料上全枚举只要几毫秒,
    // 这条判据就是空的。
    let mut m = smiles::parse(&cubanes(60)).expect("语料应能解析");
    sanitize(&mut m).expect("应能净化");
    let props = MolProps::compute(&m);

    for uniquify in [false, true] {
        let opts = MatchOptions {
            max_matches: 1,
            uniquify,
            use_chirality: true,
        };
        let t = Instant::now();
        let hits = substructure_matches(&q, &m, &props, opts);
        let el = t.elapsed();

        assert_eq!(hits.len(), 1, "uniquify={uniquify}:早停应当只返回一个匹配");
        assert!(
            el < Duration::from_millis(50),
            "uniquify={uniquify}:早停没有截断搜索,{} 原子上找第一个匹配花了 {el:?}。\n\
             这条语料上全枚举是 730 万个匹配、约 2.2 秒,而健康值是 50–75 µs ——\n\
             耗时上到几十毫秒就说明它在往下枚举。",
            m.num_atoms()
        );
        println!(
            "{:>6} 原子  uniquify={uniquify:<5} 找到首个匹配 {el:?}",
            m.num_atoms()
        );
    }
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
