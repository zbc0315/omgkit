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

use omgkit_conf::smooth::triangle_smooth;
use omgkit_conf::{bounds, chiral, pipeline, smooth, threading};
use omgkit_core::MolBuilder;

/// 拓扑距离,封顶 4。`topo[i*n+j]` ∈ {0,1,2,3,4},4 表示"≥4 或不连通"。
fn topo_dist(mol: &MolBuilder, n: usize) -> Vec<u8> {
    let mut topo = vec![4u8; n * n];
    for start in 0..n {
        let mut d = vec![u8::MAX; n];
        d[start] = 0;
        let mut q = std::collections::VecDeque::from([start]);
        while let Some(x) = q.pop_front() {
            if d[x] >= 3 {
                continue;
            }
            let Ok(xu) = u32::try_from(x) else { continue };
            for (y, _) in mol.neighbors(xu) {
                let y = y as usize;
                if y < n && d[y] == u8::MAX {
                    d[y] = d[x] + 1;
                    q.push_back(y);
                }
            }
        }
        for j in 0..n {
            topo[start * n + j] = d[j].min(4);
        }
    }
    topo
}

/// 逐档统计越界:下标 1..=4 分别是 1-2 键 / 1-3 角 / 1-4 扭转 / 长程,
/// 每档 `(越界 >0.1 Å 的对数, 总对数)`。返回的第 0 档恒为空(拓扑距离 0 是自己)。
///
/// **`d` 是 NaN 时记成越界。** `f64::max` 按 IEEE 忽略 NaN,写成
/// `(lo-d).max(d-hi).max(0.0)` 会让 NaN 算出 `over = 0`,于是一组 NaN 坐标
/// 在这张表上拿满分 —— 那是**最好看的分数**。
///
/// 在**这个**判官里这一支够不着(调用方在它之前就把含非有限数的分子 `continue`
/// 掉了),留着是给别的调用方的:同一个坑在 `conformer_oracle` 里是活的。
fn viol_by_class(coords: &[[f64; 3]], b: &smooth::Bounds, topo: &[u8]) -> [(u64, u64); 5] {
    let n = b.len();
    let mut out = [(0u64, 0u64); 5];
    for i in 0..n {
        for j in (i + 1)..n {
            let d = ((coords[i][0] - coords[j][0]).powi(2)
                + (coords[i][1] - coords[j][1]).powi(2)
                + (coords[i][2] - coords[j][2]).powi(2))
            .sqrt();
            let c = topo[i * n + j] as usize;
            out[c].1 += 1;
            let bad = if d.is_finite() {
                (b.lower(i, j) - d).max(d - b.upper(i, j)) > 0.1
            } else {
                true
            };
            if bad {
                out[c].0 += 1;
            }
        }
    }
    out
}

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
/// 一路降下来:0.65% → 0.54%(修了角的回退顺序)→ 0.34%(修了 1-3 的写入规则)
/// → **0.06%**(1-4 的五个距离改用中点而非上限)。现在是 RDKit 那条线的 **1/9**。
///
/// 闸设在 0.12%(约 10 个分子),贴着现值 5 个分子的棘轮,留了一倍余量。
///
/// 另外这条闸放在**全语料 8831 个分子**上,不放在 400 个分子的判官里:
/// 400 个样本上真实率 0.34% 只对应 1.4 个分子,泊松噪声足以让闸随机红绿。
const MAX_INFEASIBLE_FRAC: f64 = 0.0012;

/// 整条流水线跑完之后,**坐标里有原子完全重合**的分子数上限。必须是 0。
///
/// # 这一条是补上的一个洞
///
/// 对称分子的 Gram 矩阵有重特征值,等价原子会拿到逐位相同的坐标 ——
/// 而完全重合的两个原子**梯度恰好为零**,优化器永远分不开。坐标照样返回,
/// 只是废的:**静默的错**。
///
/// 实测(破对称落地之前):`large.smi` 44 个分子(0.50%)、`hard.smi` 8 个(11.8%)。
/// 0.50% 已经与 RDKit 的整体失败率同量级。
///
/// **先前没有任何判据看得见它。** 端到端判官跑的是带手性中心的分子,
/// 而手性中心本身就破了对称 —— 判据的输入分布系统性地排除了要测的那一档。
/// 所以这一条放在**全语料**这边,它的样本里有的是对称分子。
const MAX_COINCIDENT: u64 = 0;

/// 坐标里有非有限数的分子数上限。必须是 0 —— NaN 一路淌下去,查起来极贵。
const MAX_NONFINITE: u64 = 0;

/// **1-2 键**越界 >0.1 Å 的对占比上限。键是最硬的一档:参数表给的区间宽只有
/// `2×DIST12_TOL = 0.02 Å`,压不下去说明产物在化学上就是错的。
///
/// # 这一条为什么现在才有
///
/// 先前**唯一**守几何的是 `conformer_oracle`,而 CI 只喂它
/// `smoke.chirality.jsonl`(150 个药物样分子);跑全语料的正是这个判官,
/// 而它一条几何都不查。**闸有(那边 2%)、会让它红的数据也有(`hard.smi`),
/// 两者从没见过面。** 接上之后当场实测:
///
/// | | 1-2 键越界 | 至少断一根键的分子 |
/// |---|---|---|
/// | `smoke.chirality`(那边跑的) | 0.000% | 0 |
/// | `large.smi` | **0.439%** | **273(3.09%)** |
/// | `hard.smi` | **4.670%** | **13(19.1%)** |
///
/// 根因是环外 1-4 扭转被按环内的值钉死(见 `bounds::ring_path_torsion`),
/// 修掉之后是 0.026% / 0.074%,断键分子 30(0.34%)/ 1(1.47%)。
const MAX_BOND_VIOL_FRAC: f64 = 0.0015;

/// **1-3 角**越界 >0.1 Å 的对占比上限。修根因前是 2.185% / 9.366%,现在 0.148% / 0.818%。
const MAX_ANGLE_VIOL_FRAC: f64 = 0.015;

/// **至少断一根键的分子**占比上限,外加一个绝对下限(小语料上一个分子就是几个百分点,
/// 按比例设闸会被量化噪声抖红)。
///
/// 现值 `large` 30 个(0.34%)、`hard` 1 个(1.47%)。**这不是终点**:真实构象能满足
/// 全部界(拿 RDKit 构象回量,越界 0 对),所以全局极小是 0,剩下的是优化器停在了
/// **简并驻点**。二苯醚是标本:C–O 两根都压到 1.211、C···C 撑到 2.421,夹角 176.4° ——
/// **醚氧被拉成直线**,而三点共线时把氧侧向挪 `h` 只让键长变 `h²/2a`,弯折方向的梯度
/// 是二阶零。与重合原子(`crate::spread`)是同一类陷阱。从产物原地再起一次 L-BFGS,
/// 能量一步不动,证实是真极小而不是没跑够。
///
/// 这一档归**确定性重试阶梯**清,清完这条闸应当能压到 0。
const MAX_BROKEN_MOL_FRAC: f64 = 0.006;
const MIN_BROKEN_ALLOWANCE: u64 = 3;

/// 有**键交叉**的分子占比上限,外加绝对下限。
///
/// 先前 `conformer_oracle` 把"键交叉 + 环穿刺"列为四件硬事之一、也打印了这个数,
/// 但它的 fatal 段**只判环穿刺** —— 键交叉一处都没人拦。
///
/// 键交叉那一条先前有一例假阳性(四面体烷 `C12C3C1C23` 的三对**对棱**几何上
/// 就相距 `a/√2 = 1.066 Å`,低于 `CROSS_TOL`,而它们被环系锁死无法互穿)。
/// 已经在 `threading::detect` 里按拓扑距离排除掉,所以这里的下限可以取 1。
const MAX_CROSS_MOL_FRAC: f64 = 0.001;
const MIN_CROSS_ALLOWANCE: u64 = 1;

/// **键都没断、却穿了环**的分子数上限。必须是 0。
///
/// # 为什么闸是这个形状,而不是"穿刺分子数 ≤ N"
///
/// 环穿刺先前**故意没有闸**:`threading::detect` 用 `.any()` 数"有没有交点",
/// 非凸环上会假阳性,给一个会撒谎的量装硬闸,红起来查的是判据不是产物。
/// 检测器现在改成数交点奇偶了(见 `omgkit_conf::threading` 的模块文档),
/// 可以上闸。
///
/// 但"穿刺分子数 ≤ N"是个**可以随手往上抬的数**。实测全语料只有 1 个分子
/// 报穿刺,而它**已经在"至少断一根键"那 29 个里面**:
/// `C12([P+](CCC3=CC=CC=C3)(C)C)CC4CC(C1)CC(C2)C4`,金刚烷笼塌了 ——
/// 那根穿过去的 C–C 键距离六元环质心只有 0.15 Å,笼上 13 对里 13 对越界,
/// 连 1-2 键 `19-0` 都只有 1.344 Å(界 `[1.517, 1.537]`),精修跑满 400 步没收敛。
///
/// 所以闸判的是**关系**:穿了环的分子,必须同时是键断了的分子。
/// 这不是定律 —— 自穿判据存在的理由正是"穿过去时每一对距离都可以完全合法"
/// (见 `threading` 的模块文档)。它是一条**回归闸**:钉住"今天每一处穿刺都
/// 由一个已经报出来的缺陷解释得了",哪天冒出一个键完好却穿了环的分子,
/// 那是新的一类,当场红。
const MAX_PIERCE_WITH_INTACT_BONDS: u64 = 0;

/// **漏了把 `/` `\` 折算成 `BondStereo` 的分子数**上限。必须是 0。
///
/// 这是一条**前置条件闸**:折算那一步(`omgkit_io::stereo::perceive_bond_stereo`)
/// 在 `omgkit-io` 里,构型流水线管不着它,只能靠调用方记得调。
/// 而"靠记得"实测就是不行 —— 整条流水线先前压根没调,于是
/// `bounds::stereo_path_torsion` 一次都没发力,双键顺反整档退回
/// "顺式到反式的全程",外部判据上 `large.smi` 里 10 个分子交付的几何站错了边。
///
/// SMILES 层的 `harness/check_ez.py` 一直全绿 —— 它跑的是 `parse → write`,
/// **不经过这条流水线**。判据的输入分布又一次排除了要测的那一档。
///
/// 谓词用 `omgkit-io` 那个(与感知**由构造保证一致**),**不在这边自己判** ——
/// "该不该折算"要用规范秩分辨等价取代基,写第二份实现两边迟早分岔。
const MAX_UNPERCEIVED_STEREO: u64 = 0;

/// **拿到了 `BondStereo` 的双键数下限** —— 语料里有方向键时必须 > 0。
///
/// 上面那条只看得见"调用方漏了调",看不见"折算这件事整体失效":
/// 变异验证过 —— 把 SMILES 解析器改成丢掉所有 `/` `\`,那条闸读 0、全绿退出 0,
/// 而断键分子当场退回修复前的数。**只让判据变绿的东西必须配一道反向闸。**
///
/// 判的是条件式:语料里**一根方向键都没有**时(`hard.smi` 就是这样)这条不生效,
/// 但两个数都会打印出来 —— 恒 0/0 的闸看着在守,实际什么都没守。
const MIN_PERCEIVED_STEREO: u64 = 1;

/// **该出构型却没出**的分子数上限。必须是 0。
///
/// # 这一条堵的是"分母静默归零"
///
/// 上面那几条几何闸读的计数器,全都在 `pipeline::conformer` 成功之后才累加。
/// 先前那一步失败是裸 `continue` —— 不计数、不报告、没有闸。于是**任何让构型
/// 生成失败率上升的回归,都会让几何闸变得更好看**:极限情形下每个分子都失败,
/// 四条几何闸读的是 `0/0`,判官照样打印"都过"并退出 0(变异验证过)。
///
/// 只让判据变绿的东西必须配一道上限闸,否则没人拦得住它。
const MAX_NO_CONFORMER: u64 = 0;

fn main() {
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "harness/corpus/large.smi".into());
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| {
        eprintln!("读不了语料 {path}:{e}");
        std::process::exit(1);
    });

    let (mut n, mut n_parse_fail, mut n_empty, mut n_infeasible) = (0u64, 0u64, 0u64, 0u64);
    let (mut n13_conflict, mut n14_degenerate) = (0u64, 0u64);
    let (mut n_unperceived, mut n_perceived, mut n_directional) = (0u64, 0u64, 0u64);
    let (mut n_conf, mut n_spread, mut n_coincident, mut n_nonfinite) = (0u64, 0u64, 0u64, 0u64);
    let mut coincident_cases: Vec<String> = Vec::new();
    let mut empty_cases: Vec<String> = Vec::new();
    let mut infeasible_cases: Vec<String> = Vec::new();
    // ---- 几何:先前**只有 150 个分子的判官在看**,全语料这边一条都没有 ----
    let mut viol = [(0u64, 0u64); 5];
    let (mut n_broken_bond, mut n_cross_mol, mut n_pierce_mol) = (0u64, 0u64, 0u64);
    let mut n_pierce_intact = 0u64;
    let mut pierce_intact_cases: Vec<String> = Vec::new();
    let (mut n_cross, mut worst_over) = (0u64, 0.0f64);
    let mut broken_cases: Vec<String> = Vec::new();
    let mut thread_cases: Vec<String> = Vec::new();
    let mut no_conf_cases: Vec<String> = Vec::new();
    let mut worst_case = String::new();
    let mut n_no_conf = 0u64;

    for line in text.lines() {
        // 语料格式:`SMILES<TAB>名字`,**`#` 开头是注释、空行忽略**
        // (`large.smi` 里有 24 条被注释掉的分子,`bridged.smi` 与 `hard.smi`
        //  的抬头都是整段注释)。先前这里没认注释,把它们当 SMILES 解析,
        // 于是"解析失败"那个数里混着注释行 —— 分母没受影响(它只数建成功的),
        // 但报出来的数是错的。
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
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
        // **把 `/` `\` 折算成双键自己的 `BondStereo`。** SMILES 里顺反记在相邻单键的
        // `direction` 上,不经这一步的话 `bounds::stereo_path_torsion` 一次都不发力,
        // 而它是 1-4 扭转最硬的那条来源。整条流水线先前少的就是这一环。
        // 折算**之前**先记一笔:有几根方向键。这是纯拓扑事实,
        // 用来判断下面那条反向闸该不该生效。
        n_directional += mol
            .bonds()
            .iter()
            .filter(|b| b.direction != omgkit_core::BondDirection::None)
            .count() as u64;
        omgkit_io::stereo::perceive_bond_stereo(&mut mol);
        n_perceived += mol
            .bonds()
            .iter()
            .filter(|b| b.stereo != omgkit_core::BondStereo::None)
            .count() as u64;
        // 谓词与感知由构造保证一致,所以感知跑过之后它必须闭嘴
        n_unperceived += u64::from(omgkit_io::stereo::directions_not_perceived(&mol));
        // 补氢要给一个与写法无关的秩;这里只关心界可不可行,用恒等秩即可
        let order: Vec<u32> = (0..mol.num_atoms() as u32).collect();
        omgkit_chem::explicit_hs::add_explicit_hs(&mut mol, &order);
        n += 1;
        let (mut b, stats) = bounds::build(&mol);
        // **把 Stats 里那两个"退让计数"报出来。** 它们记的是参数表自相矛盾
        // (`n13_conflict`)与几何退化丢约束(`n14_degenerate`)——
        // 只写不读的计数器等于没有,而这两件事都只会让界更松、判据更容易绿。
        n13_conflict += stats.n13_conflict as u64;
        n14_degenerate += stats.n14_degenerate as u64;

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
            continue;
        }

        // ---- 跑完整条流水线,查两条**硬**不变量 ----
        //
        // 光有"界自洽"是不够的:界自洽的分子照样可能吐出一组废坐标
        // (原子完全重合、或者 NaN)。这两条必须逐分子成立,不是统计。
        let centers: Vec<chiral::Center> = chiral::centers(&mol);
        // **失败必须计数。** 后面几条几何闸读的计数器都在这一行之后才累加,
        // 裸 `continue` 会让"生成失败"变成"几何更好看"。见 `MAX_NO_CONFORMER`。
        let conf = match pipeline::conformer(&mol, &centers) {
            Ok(c) => c,
            Err(e) => {
                n_no_conf += 1;
                if no_conf_cases.len() < 6 {
                    no_conf_cases.push(format!("{smi}({e:?})"));
                }
                continue;
            }
        };
        n_conf += 1;
        n_spread += u64::from(conf.spread > 0);
        if conf.coords.iter().any(|p| p.iter().any(|v| !v.is_finite())) {
            n_nonfinite += 1;
            continue;
        }
        let na = conf.coords.len();
        let mut coincident = false;
        for i in 0..na {
            for j in (i + 1)..na {
                let d2: f64 = (0..3)
                    .map(|t| (conf.coords[i][t] - conf.coords[j][t]).powi(2))
                    .sum();
                if d2 <= 1e-12 {
                    coincident = true;
                }
            }
        }
        if coincident {
            n_coincident += 1;
            if coincident_cases.len() < 6 {
                coincident_cases.push(smi.to_string());
            }
        }

        // ---- 几何:界是"能不能摆",这里量的是"摆得对不对" ----
        //
        // 这两件事**不是一回事**,而先前只有 `conformer_oracle` 在看后者,
        // 它跑的是 150 个药物样分子 —— 闸有(2%)、会让它红的数据也有(hard.smi),
        // 两者从没见过面。
        let topo = topo_dist(&mol, na);
        let v = viol_by_class(&conf.coords, &b, &topo);
        for (k, (bad, tot)) in v.iter().enumerate() {
            viol[k].0 += bad;
            viol[k].1 += tot;
        }
        let bonds_broken = v[1].0 > 0;
        if bonds_broken {
            n_broken_bond += 1;
            if broken_cases.len() < 6 {
                broken_cases.push(smi.to_string());
            }
        }
        for i in 0..na {
            for j in (i + 1)..na {
                let d = ((conf.coords[i][0] - conf.coords[j][0]).powi(2)
                    + (conf.coords[i][1] - conf.coords[j][1]).powi(2)
                    + (conf.coords[i][2] - conf.coords[j][2]).powi(2))
                .sqrt();
                let over = (b.lower(i, j) - d).max(d - b.upper(i, j));
                if over > worst_over {
                    worst_over = over;
                    worst_case = smi.to_string();
                }
            }
        }
        let t = threading::detect(&mol, &conf.coords);
        n_cross += t.crossings as u64;
        n_cross_mol += u64::from(t.crossings > 0);
        n_pierce_mol += u64::from(t.pierces > 0);
        if t.pierces > 0 && !bonds_broken {
            n_pierce_intact += 1;
            if pierce_intact_cases.len() < 6 {
                pierce_intact_cases.push(smi.to_string());
            }
        }
        if (t.crossings > 0 || t.pierces > 0) && thread_cases.len() < 6 {
            thread_cases.push(format!("{smi}(交叉{} 穿刺{})", t.crossings, t.pierces));
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
    println!("  1-3 两条估计交空退并集 {n13_conflict} 次;1-4 几何退化丢约束 {n14_degenerate} 次");
    println!(
        "  顺反:方向键 {n_directional} 根 → 折算出 {n_perceived} 根双键立体\
         (下限 {MIN_PERCEIVED_STEREO},仅当方向键 > 0 时生效);\
         **折算完谓词还在报**的分子 {n_unperceived} 个(上限 {MAX_UNPERCEIVED_STEREO})"
    );
    println!("  ── 整条流水线跑完的硬不变量(逐分子,不是统计)──");
    println!(
        "    出了构型的分子 {n_conf};**该出没出**的 {n_no_conf}(上限 {MAX_NO_CONFORMER});\
         破对称动过的 {n_spread}"
    );
    if !no_conf_cases.is_empty() {
        println!("    没出构型的例子:{}", no_conf_cases.join("  "));
    }
    println!("    **原子完全重合**的 {n_coincident}(上限 {MAX_COINCIDENT});坐标含非有限数的 {n_nonfinite}(上限 {MAX_NONFINITE})");
    if !coincident_cases.is_empty() {
        println!("    重合的例子:{}", coincident_cases.join("  "));
    }
    #[allow(clippy::cast_precision_loss)]
    let pct = |a: u64, b: u64| 100.0 * a as f64 / b.max(1) as f64;
    println!("  ── 几何:出厂坐标越界 >0.1 Å 的对,按拓扑档 ──");
    for (c, name) in [
        (1usize, "1-2 键"),
        (2, "1-3 角"),
        (3, "1-4 扭转"),
        (4, "长程"),
    ] {
        println!(
            "    {name:8} {:9} 对中 {:7} 越界  {:6.3}%",
            viol[c].1,
            viol[c].0,
            pct(viol[c].0, viol[c].1)
        );
    }
    println!(
        "    **至少断一根键**的分子 {n_broken_bond}({:.2}%);最坏越界 {worst_over:.2} Å  {worst_case}",
        pct(n_broken_bond, n_conf)
    );
    if !broken_cases.is_empty() {
        println!("    断键的例子:{}", broken_cases.join("  "));
    }
    println!(
        "    键交叉 {n_cross} 处,分布在 {n_cross_mol} 个分子;有环穿刺的分子 {n_pierce_mol},\
         其中**键都没断**的 {n_pierce_intact}(上限 {MAX_PIERCE_WITH_INTACT_BONDS})"
    );
    if !thread_cases.is_empty() {
        println!("    自穿的例子:{}", thread_cases.join("  "));
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
    if n_coincident > MAX_COINCIDENT {
        eprintln!(
            "\n有 {n_coincident} 个分子的坐标里原子完全重合 —— 那里梯度恰好为零,\
             优化器永远分不开,坐标是废的。破对称那一步(crate::spread)可能没接上"
        );
        fatal = true;
    }
    if n_nonfinite > MAX_NONFINITE {
        eprintln!("\n有 {n_nonfinite} 个分子的坐标含非有限数");
        fatal = true;
    }
    if n_unperceived > MAX_UNPERCEIVED_STEREO {
        eprintln!(
            "\n有 {n_unperceived} 个分子折算完谓词还在报漏 —— \
             要么调用方漏了 `omgkit_io::stereo::perceive_bond_stereo`,\
             要么谓词与感知分岔了。这一档的 1-4 扭转会退回顺式到反式的全程,\
             交付的几何有一半站错边"
        );
        fatal = true;
    }
    if n_directional > 0 && n_perceived < MIN_PERCEIVED_STEREO {
        eprintln!(
            "\n语料里有 {n_directional} 根方向键,却一根双键立体都没折算出来 —— \
             折算这件事整体失效了。上面那条闸看不见这种失败(它只看得见\
             '调用方漏了调'),所以必须有这条反向闸"
        );
        fatal = true;
    }
    if n_no_conf > MAX_NO_CONFORMER {
        eprintln!(
            "\n有 {n_no_conf} 个分子界可行却没出构型(上限 {MAX_NO_CONFORMER})—— \
             下面几条几何闸的分母全在这一步之后累加,放着不管等于给它们开后门"
        );
        fatal = true;
    }
    // ---- 几何:先前全语料这边一条都没有 ----
    let bond_frac = pct(viol[1].0, viol[1].1) / 100.0;
    if bond_frac > MAX_BOND_VIOL_FRAC {
        eprintln!(
            "\n1-2 键越界 {:.3}% > {:.3}% —— 最硬的一档,产物在化学上就是错的",
            100.0 * bond_frac,
            100.0 * MAX_BOND_VIOL_FRAC
        );
        fatal = true;
    }
    let angle_frac = pct(viol[2].0, viol[2].1) / 100.0;
    if angle_frac > MAX_ANGLE_VIOL_FRAC {
        eprintln!(
            "\n1-3 角越界 {:.3}% > {:.3}%",
            100.0 * angle_frac,
            100.0 * MAX_ANGLE_VIOL_FRAC
        );
        fatal = true;
    }
    #[allow(
        clippy::cast_precision_loss,
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss
    )]
    let allow = |frac: f64, floor: u64| ((n_conf as f64 * frac).ceil() as u64).max(floor);
    let broken_allowed = allow(MAX_BROKEN_MOL_FRAC, MIN_BROKEN_ALLOWANCE);
    if n_pierce_intact > MAX_PIERCE_WITH_INTACT_BONDS {
        eprintln!(
            "\n有 {n_pierce_intact} 个分子**键都没断却穿了环** > 上限 \
             {MAX_PIERCE_WITH_INTACT_BONDS} —— 这是新的一类缺陷:\n\
             它不会被越界那几档看见(穿过去时每一对距离都可以完全合法)\n\
             例子:{}",
            pierce_intact_cases.join("  ")
        );
        fatal = true;
    }
    if n_broken_bond > broken_allowed {
        eprintln!(
            "\n有 {n_broken_bond} 个分子至少断一根键 > 允许的 {broken_allowed} 个\
             (共 {n_conf};上限 {:.1}%)",
            100.0 * MAX_BROKEN_MOL_FRAC
        );
        fatal = true;
    }
    let cross_allowed = allow(MAX_CROSS_MOL_FRAC, MIN_CROSS_ALLOWANCE);
    if n_cross_mol > cross_allowed {
        eprintln!(
            "\n有 {n_cross_mol} 个分子出现键交叉 > 允许的 {cross_allowed} 个 —— \
             距离判据看不见自穿,穿过去时每一对距离都可以合法"
        );
        fatal = true;
    }
    if fatal {
        std::process::exit(1);
    }
    // 逐条点名,别只报个数 —— 加了闸忘了改数,或者改了数没加闸,都是这么来的。
    println!(
        "\n十一条都过:空区间 / 界不可行 / 顺反没折算 / 原子重合 / 非有限数 / 该出没出 / \
         1-2 键 / 1-3 角 / 断键分子 / 键交叉分子 / 键没断却穿环的分子。"
    );
}
