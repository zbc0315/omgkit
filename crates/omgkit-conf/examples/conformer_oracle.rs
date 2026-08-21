//! **端到端判官:整条流水线跑完之后,构型到底好不好。**
//!
//! 前面几条判官各守一段(界矩阵、特征分解、嵌入、手性、自穿),这一条守的是
//! **产物**:分子进去、坐标出来,那组坐标满不满足化学。
//!
//! 量四件事,**精修前后各量一遍** —— 只报"之后"看不出精修有没有在干活:
//!
//! | 量 | 为什么它是硬的 |
//! |---|---|
//! | 键长/键角越界(按拓扑档分) | 长程对占 87%,只看总数会被稀释 |
//! | 自穿(键交叉 + 环穿刺) | 距离判据看不见它 —— 穿过去时每一对距离都可以合法 |
//! | 手性号 | 立体化学错了分子就是错的,不是"差一点" |
//! | 耗时 | 要与 RDKit 比的那条线 |
//!
//! ```shell
//! cargo run -p omgkit-conf --release --example conformer_oracle -- harness/baseline/smoke.chirality.jsonl
//! ```

use omgkit_conf::chiral::Center;
use omgkit_conf::{bounds, chiral, embed, pipeline, smooth, spread, threading};
use omgkit_core::{BondOrder, ChiralTag, MolBuilder};

/// 精修之后,**键长**(1-2)越界超过这个数的对占比上限。
///
/// 键是最硬的一档:参数表给的区间宽只有 `2×DIST12_TOL = 0.02 Å`,
/// 而精修的误差函数罚**相对**越界,键上最贵 —— 压不下去说明优化器没干活。
///
/// 先前是 **2%**,而实测早已是 **0.0%**(4202 对里 0 对)—— 那不是闸,是橡皮图章:
/// 84 根键推到越界都拦不住。收到 0.2%(约 9 对),贴着现值的棘轮。
/// 全语料那一档的同名闸在 `feasibility` 里,设在 0.15%。
const MAX_BOND_VIOL_FRAC: f64 = 0.002;

/// 精修之后仍有环穿刺的分子占比上限(小样本给下限,同 `threading_oracle`)。
///
/// **这里保留余量是有理由的**:环穿刺用质心扇形三角化 + `.any()`,
/// 非凸环上偶数次相交会假阳性(见 `threading` 模块文档)。在检测器改成数交点奇偶
/// 之前,不把它收成 0。
const MAX_PIERCE_FRAC: f64 = 0.05;
const MIN_PIERCE_ALLOWANCE: u64 = 2;

/// 精修之后仍有**键交叉**的分子数上限。
///
/// 先前这个数只 `println!`,fatal 段**只判环穿刺** —— 于是"键交叉 1731 → 0"
/// 这行漂亮数字底下一个闸都没有,1731 处一处不减也照样打印"判据全过"。
/// 键交叉与环穿刺不同,已经按拓扑距离排掉了刚性假阳性(见 `threading::RIGID_TOPO`),
/// 所以这一条按硬不变量设成 **0**。
const MAX_CROSS_MOL: u64 = 0;

/// 基准里允许有几个分子**根本没量到**。
///
/// # 这一条堵的是分母
///
/// 上面每一个计数器(键长越界、环穿刺、键交叉、手性真值)都在三处
/// `n_fail += 1; continue;` **之后**才累加,而 `n_fail` 先前只被 `println!`,
/// 一道闸都没有 —— 唯一的下限是 `n_mol == 0`。
///
/// 变异实测(让 `pipeline::conformer` 对超过 25 个原子的分子直接失败):
///
/// ```text
/// 判官:端到端构型,分子 66 个(失败 84)     ← 基准 150 行
/// 端到端判据全过。                          ← 退出码 0
/// ```
///
/// 孤对那一档更狠:15 → 11 个,照样退 0 —— 而 CI 里那一步存在的**全部理由**
/// 就是"上面那份基准里一个三配位中心都没有",它的分母没有任何反向闸。
///
/// 与 `feasibility` 的 `MAX_NO_CONFORMER = 0` 是同一条:几何判据的计数器都在
/// 生成成功之后才累加,不给它配闸,任何让失败率上升的回归都会让判据变好看。
/// 实测两份基准都是 0 个失败,所以钉死 0。
const MAX_FAIL: u64 = 0;

/// 基准里允许有几行**没进比对**(解析失败、字段缺失等)。
///
/// `n_mol + n_fail` 应当等于基准的行数;差额是在读入阶段被 `continue` 掉的。
/// 实测两份基准都是 0,所以同样钉死 0。
const MAX_UNREAD: u64 = 0;

/// 精修之后手性号正确的比例下限。
///
/// 立体化学错了分子就是错的,不是"差一点",所以这是**硬不变量:1.0**。
/// 先前设 0.90,而实测早已是 100.0% —— 247 个中心里 24 个翻号都拦不住。
///
/// # 但要清楚这一条量的是什么
///
/// `conf.chiral_ok` 是拿**驱动 `Field` 的那同一份 `centers`** 数出来的,
/// 所以它答的是"优化器有没有打中自己的靶",**不是**"立体化学对不对" ——
/// `chiral::centers` 若系统性地抽错,这里照样 100%。
/// 真值判据在 `chiral_oracle`(真值取自真实构象),两条要一起看。
const MIN_CHIRAL_OK: f64 = 1.0;

/// 逐档统计越界:`[1-2, 1-3, 1-4, 长程]` 的 `(越界数, 总数)`。
///
/// **`d` 是 NaN 时记成越界。** `f64::max` 按 IEEE 忽略 NaN,先前这里写的是
/// `(lo-d).max(d-hi).max(0.0)` —— NaN 算出 `over = 0`,于是一组 NaN 坐标
/// 四档越界全报 0.0%,拿到的是**最好看的分数**;自穿也报 0,没有手性中心的分子
/// 连那条闸都碰不到。方向是反的:坏输入 → 满分。
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

/// 拓扑距离,封顶 4。
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

fn main() {
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "harness/baseline/smoke.chirality.jsonl".into());
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| {
        eprintln!("读不了判官基准 {path}:{e}");
        std::process::exit(1);
    });

    // `n_lines` 是基准行数(分母),`n_mol` 是真正量到的,`n_fail` 是生成失败的。
    let (mut n_lines, mut n_mol, mut n_fail) = (0u64, 0u64, 0u64);
    let mut before = [(0u64, 0u64); 5];
    let mut after = [(0u64, 0u64); 5];
    let (mut cross_b, mut cross_a) = (0u64, 0u64);
    let (mut pierce_mol_b, mut pierce_mol_a) = (0u64, 0u64);
    let (mut chi_total, mut chi_before, mut chi_after) = (0u64, 0u64, 0u64);
    let (mut truth_total, mut truth_ok, mut truth_declared) = (0u64, 0u64, 0u64);
    let mut truth_bad: Vec<String> = Vec::new();
    let (mut iters, mut e_before, mut e_after) = (0u64, 0.0f64, 0.0f64);
    let start = std::time::Instant::now();

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
        let nat = v["n"].as_u64().unwrap_or(0) as usize;
        let z: Vec<u8> = v["z"]
            .as_array()
            .map(|a| {
                a.iter()
                    .filter_map(|x| x.as_u64().map(|y| y as u8))
                    .collect()
            })
            .unwrap_or_default();
        let chg: Vec<i8> = v["charge"]
            .as_array()
            .map(|a| {
                a.iter()
                    .filter_map(|x| x.as_i64().map(|y| y as i8))
                    .collect()
            })
            .unwrap_or_default();
        if z.len() != nat {
            continue;
        }
        let mut m = MolBuilder::new();
        for (k, &a) in z.iter().enumerate() {
            let mut ad = omgkit_core::AtomData::new(a);
            ad.formal_charge = chg.get(k).copied().unwrap_or(0);
            m.add_atom_data(ad);
        }
        let mut ok = true;
        for e in v["bonds"].as_array().into_iter().flatten() {
            let Some(t) = e.as_array() else { continue };
            let (Some(i), Some(j), Some(o)) = (
                t.first().and_then(serde_json::Value::as_u64),
                t.get(1).and_then(serde_json::Value::as_u64),
                t.get(2).and_then(serde_json::Value::as_u64),
            ) else {
                continue;
            };
            let ord = match o {
                2 => BondOrder::Double,
                3 => BondOrder::Triple,
                4 => BondOrder::Aromatic,
                _ => BondOrder::Single,
            };
            if m.add_bond(i as u32, j as u32, ord).is_err() {
                ok = false;
            }
        }
        if !ok || m.num_atoms() != nat {
            continue;
        }
        for c in v["centers"].as_array().into_iter().flatten() {
            let (Some(atom), Some(tag)) = (c["atom"].as_u64(), c["tag"].as_i64()) else {
                continue;
            };
            if let Some(a) = m.atom_mut(atom as u32) {
                a.chiral_tag = if tag == 2 {
                    ChiralTag::Ccw
                } else {
                    ChiralTag::Cw
                };
            }
        }
        let centers: Vec<Center> = chiral::centers(&m);
        // **真值那一份**:中心原子、配体序、以及真实构象上算出来的号,全部取自基准。
        // 完全绕开 `chiral::centers` —— 这样"抽错中心"与"摆错几何"就不会互相掩盖。
        let truth: Vec<(usize, [usize; 3], f64)> = v["centers"]
            .as_array()
            .into_iter()
            .flatten()
            .filter_map(|c| {
                let atom = c["atom"].as_u64()? as usize;
                let sign = c["sign"].as_f64()?;
                let nb = c["nbrs"].as_array()?;
                let g = |k: usize| nb.get(k)?.as_u64().map(|x| x as usize);
                Some((atom, [g(0)?, g(1)?, g(2)?], sign))
            })
            .collect();
        // **分母不许静默缩小。** `filter_map` 里任一 `?` 落空(基准的 `nbrs`/`sign`
        // 字段缺失或改名)都会悄悄丢掉那个中心,于是"号对的 N/N,100%"照样绿 ——
        // 实测把基准里除第一个外的 `nbrs` 全删掉,判官打印"真值口径 1 个:1(100.0%)"
        // 并退出 0。基准里**声明**了几个,留着待会儿对账。
        let declared = v["centers"].as_array().map_or(0, Vec::len) as u64;

        // ---- 精修前:嵌入 + 全局反射 ----
        let (mut b, _) = bounds::build(&m);
        if smooth::triangle_smooth(&mut b).is_err() {
            n_fail += 1;
            continue;
        }
        let Ok(emb) = embed::embed(&embed::reference_distances(&b), nat) else {
            n_fail += 1;
            continue;
        };
        let mut pre = emb.coords;
        // **破对称也要跑** —— 否则"精修前"量的是流水线里根本不存在的中间态。
        // 而且漏掉它会让 `needs_reflection` 偏:重合配体的有符号体积恒为 0,
        // 那些中心一律计成"号不对"。
        spread::break_coincidence(&mut pre);
        if chiral::needs_reflection(&pre, &centers) {
            chiral::reflect(&mut pre);
        }

        // ---- 精修后 ----
        let Ok(conf) = pipeline::conformer(&m, &centers) else {
            n_fail += 1;
            continue;
        };
        n_mol += 1;
        iters += conf.iterations as u64;
        e_before += conf.energy_before;
        e_after += conf.energy;

        let topo = topo_dist(&m, nat);
        for (k, (bad, tot)) in viol_by_class(&pre, &b, &topo).iter().enumerate() {
            before[k].0 += bad;
            before[k].1 += tot;
        }
        for (k, (bad, tot)) in viol_by_class(&conf.coords, &b, &topo).iter().enumerate() {
            after[k].0 += bad;
            after[k].1 += tot;
        }
        let tb = threading::detect(&m, &pre);
        let ta = threading::detect(&m, &conf.coords);
        cross_b += tb.crossings as u64;
        cross_a += ta.crossings as u64;
        if tb.pierces > 0 {
            pierce_mol_b += 1;
        }
        if ta.pierces > 0 {
            pierce_mol_a += 1;
        }
        chi_total += centers.len() as u64;
        chi_before += chiral::correct_count(&pre, &centers) as u64;
        chi_after += conf.chiral_ok as u64;

        // ---- 真值口径:在**我们交付的坐标**上,按真值的中心与配体序复算 ----
        // 对账放在这里,不放在解析处 —— 中途 `continue` 掉的分子(光滑化/嵌入失败)
        // 两个计数都不该加,否则会假红。
        truth_declared += declared;
        for &(atom, lig, sign) in &truth {
            if atom >= nat || lig.iter().any(|&k| k >= nat) {
                continue;
            }
            let o = conf.coords[atom];
            let d = |k: usize| {
                let p = conf.coords[lig[k]];
                [p[0] - o[0], p[1] - o[1], p[2] - o[2]]
            };
            let (a, b, c) = (d(0), d(1), d(2));
            let vol = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0]);
            truth_total += 1;
            if vol != 0.0 && vol.signum() == sign {
                truth_ok += 1;
            } else if truth_bad.len() < 6 {
                truth_bad.push(format!(
                    "{}#{atom}(V={vol:+.3} 目标号 {sign:+})",
                    v["smiles"]
                ));
            }
        }
    }
    let elapsed = start.elapsed();

    #[allow(clippy::cast_precision_loss)]
    let pct = |a: u64, b: u64| 100.0 * a as f64 / b.max(1) as f64;
    let unread = n_lines.saturating_sub(n_mol + n_fail);
    println!(
        "判官:端到端构型,基准 {n_lines} 行,真正量到 {n_mol} 个(生成失败 {n_fail},没读进来 {unread})"
    );
    println!(
        "  精修:平均 {} 步;误差 {:.3e} → {:.3e}(降 {:.1}%)",
        iters / n_mol.max(1),
        e_before / n_mol.max(1) as f64,
        e_after / n_mol.max(1) as f64,
        100.0 * (1.0 - e_after / e_before.max(1e-30))
    );
    println!("  ── 越界 >0.1 Å 的对,按拓扑档(精修前 → 后)──");
    for (c, name) in [
        (1usize, "1-2 键"),
        (2, "1-3 角"),
        (3, "1-4 扭转"),
        (4, "长程"),
    ] {
        println!(
            "    {name:8} 对数 {:6}   {:5.1}% → {:5.1}%",
            after[c].1,
            pct(before[c].0, before[c].1),
            pct(after[c].0, after[c].1)
        );
    }
    println!("  ── 自穿 ──");
    println!("    键交叉 {cross_b} → {cross_a};有环穿刺的分子 {pierce_mol_b} → {pierce_mol_a}(共 {n_mol})");
    println!("  ── 手性 ──");
    println!(
        "    自洽口径(我们抽的中心)  中心 {chi_total} 个:{:.1}% → {:.1}%",
        pct(chi_before, chi_total),
        pct(chi_after, chi_total)
    );
    println!(
        "    **真值口径**(中心与配体序取自基准,在我们交付的坐标上复算)\
         {truth_total} 个:号对的 {truth_ok}({:.1}%)",
        pct(truth_ok, truth_total)
    );
    if !truth_bad.is_empty() {
        println!("      号不对的:{}", truth_bad.join("  "));
    }
    println!(
        "  ── 耗时 ── 合计 {:.2} 秒,每分子 {:.2} ms",
        elapsed.as_secs_f64(),
        1000.0 * elapsed.as_secs_f64() / n_mol.max(1) as f64
    );

    let mut fatal = false;
    // **分母闸,排在其它闸之前。** 少量到几个分子,下面每一个比例都会变好看。
    if n_fail > MAX_FAIL {
        eprintln!(
            "\n{n_fail} 个分子没能生成构型(上限 {MAX_FAIL})—— 下面每一条几何判据的\n\
             计数器都在生成成功之后才累加,失败率一涨,那些数就会**变好看**"
        );
        std::process::exit(1);
    }
    if unread > MAX_UNREAD {
        eprintln!(
            "\n基准里有 {unread} 行没进比对(上限 {MAX_UNREAD})—— 分母核不上,\n\
             这条判据算出来的比例没有意义"
        );
        std::process::exit(1);
    }
    if n_mol == 0 {
        eprintln!("\n一个分子都没跑成");
        fatal = true;
    }
    let bond_frac = pct(after[1].0, after[1].1) / 100.0;
    if bond_frac > MAX_BOND_VIOL_FRAC {
        eprintln!(
            "\n精修之后键长越界的对占 {:.1}% > {:.0}% —— 最硬的一档没压下去",
            100.0 * bond_frac,
            100.0 * MAX_BOND_VIOL_FRAC
        );
        fatal = true;
    }
    #[allow(
        clippy::cast_precision_loss,
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss
    )]
    let allowed = ((n_mol as f64 * MAX_PIERCE_FRAC).ceil() as u64).max(MIN_PIERCE_ALLOWANCE);
    if pierce_mol_a > allowed {
        eprintln!("\n精修之后仍有 {pierce_mol_a} 个分子环穿刺 > 允许的 {allowed} 个");
        fatal = true;
    }
    if cross_a > MAX_CROSS_MOL {
        eprintln!(
            "\n精修之后还有 {cross_a} 处键交叉(上限 {MAX_CROSS_MOL})—— \
             距离判据看不见自穿,穿过去时每一对距离都可以合法"
        );
        fatal = true;
    }
    if truth_total > 0 && truth_ok < truth_total {
        eprintln!(
            "\n真值口径下有 {} 个中心的号不对 —— 这一条**绕开了我们自己的抽取逻辑**,\
             红了就是交付的坐标里立体化学是错的(整分子对映体)",
            truth_total - truth_ok
        );
        fatal = true;
    }
    if truth_total != truth_declared {
        eprintln!(
            "\n基准里写着 {truth_declared} 个中心,真值口径只读进 {truth_total} 个 —— \
             差的那些被 `filter_map` 静默丢了(字段缺失/改名/下标越界)。\
             分母缩小只会让上面那条闸更好看,必须当场红"
        );
        fatal = true;
    }
    if truth_total == 0 && chi_total > 0 {
        eprintln!("\n基准里一个真值中心都没读到 —— 上面那条闸的分母是 0,等于没在看");
        fatal = true;
    }
    if chi_total > 0 && pct(chi_after, chi_total) / 100.0 < MIN_CHIRAL_OK {
        eprintln!(
            "\n精修之后手性正确率 {:.1}% < {:.0}% —— 个别中心那一档本该由精修救回来",
            pct(chi_after, chi_total),
            100.0 * MIN_CHIRAL_OK
        );
        fatal = true;
    }
    if fatal {
        std::process::exit(1);
    }
    println!("\n端到端判据全过。");
}
