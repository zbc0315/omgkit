//! 手性真值里那两列**有符号体积**的读取方。
//!
//! `smoke.chirality.jsonl` / `smoke.lonepair.jsonl` 的每个中心带四样东西:
//! `nbrs`(RDKit 的邻居序)、`sign`(±1 的真值)、`vol`(以中心为基点的行列式)、
//! `vol_ligand`(四配体行列式,旧口径)。**前两样有五个读取方,后两样一个都没有** ——
//! 从入库到现在,把 `vol` 全改成 0 也没有任何判据会红。
//!
//! # 只比号是看不见公式的
//!
//! 现有判据(`examples/chiral_oracle.rs`)拿真实坐标复算体积,但**只比号**。
//! 号只有两种取值,所以只有"整体反号"那一类错抓得到。实测三种变异:
//!
//! | 变异 `chiral::center_volume` | 号 | `chiral_oracle` | 本文件 |
//! |---|---|---|---|
//! | 基点从中心原子换成配体 0 | 整体反 | **退 1**(248 个号对不上) | 红 |
//! | 结果乘 2 | 不变 | **退 0**,三行数字与未变异时逐字相同 | 红 |
//! | 换一个配体三元组(`d(1),d(2),d(3)`) | 多数不变 | 退 101(三配位那格是孤对,越界 panic) | 红 |
//!
//! 第二行是这条判据存在的理由:**一个差常数因子的体积公式,所有比号的判据
//! 一个数都不会动** —— "符号预测错 0 个"、"号对不上的 0 个"、"全部中心都对
//! 76.7%" 原样照印。
//!
//! 这一类不是纸上假想:`dump_chirality.py` 头一版用的就是四配体口径,注释还
//! 写着"与 `chiral.rs::signed_volume` 同一个式子" —— 真值与待验实现同一条式子,
//! 那个判官在结构上就抓不到"中心被挤出配体四面体"那一档。
//!
//! 所以这里比**数值**:用产品的 [`chiral::center_volume`],按真值给的配体序、
//! 在真值给的坐标上复算,必须与 `vol` 逐个相同。同一份输入、同一个定义,
//! 两个独立实现(一个 Python 在 dump 里,一个 Rust 在产品里)必须给出同一个数。
//!
//! # 另外两条,守的是"真值本身站不站得住"
//!
//! **中心必须在配体四面体里面。** 正四面体上 `V_配体 = −4·V_中心`,两者恒反号;
//! **同号就意味着中心原子被挤到四个配体张成的四面体外面**(伞形翻转),那时
//! `sign` 说的构型与分子实际的构型已经对不上了。`dump_chirality.py` 在导的时候
//! 会数这一档并据此定退出码,但那是**导的那一刻**的事 —— 入库之后没人再看。
//!
//! **号不能压在一个近乎共面的中心上。** 无量纲平度 `|V| / (|a||b||c|)`:
//! 三个单位向量两两 109.47° 时是 `|1−t|·√(1+2t) = 0.7698`(`t = cos109.47° = −1/3`)。
//! 实测这两份基准的中位数是 0.7529 / 0.8904,最小 0.5438 / 0.7331 —— 离共面很远。
//! 真值若哪天出现一个平的中心,它的 `sign` 就是掷硬币,而"我们号对了"也就没有意义。

use std::path::PathBuf;

use omgkit_conf::chiral::{self, Center};

/// 复算与 `vol` 的差多大算不同。实测最大差 2.7e-15(体积量级 2–10),
/// 留到 1e-9 是给不同平台的浮点求和顺序留的余量,不是给公式差异留的。
const VOL_TOL: f64 = 1e-9;

/// 无量纲平度的下限。正四面体是 0.7698,实测最小 0.5438 ——
/// 0.2 这条线两边都远:真值不会误红,而一个真的压平了的中心(0.05 量级)拦得住。
const MIN_FLATNESS: f64 = 0.2;

struct Rec {
    smiles: String,
    coords: Vec<[f64; 3]>,
    centers: Vec<serde_json::Value>,
}

fn read(name: &str) -> Vec<Rec> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../harness/baseline")
        .join(name);
    let text =
        std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("读不到 {}: {e}", path.display()));
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|line| {
            let v: serde_json::Value = serde_json::from_str(line).expect("每行都是合法 JSON");
            let coords = v["coords"]
                .as_array()
                .expect("有 coords")
                .iter()
                .map(|c| {
                    let a = c.as_array().expect("坐标是三元组");
                    [
                        a[0].as_f64().expect("x"),
                        a[1].as_f64().expect("y"),
                        a[2].as_f64().expect("z"),
                    ]
                })
                .collect();
            Rec {
                smiles: v["smiles"].as_str().unwrap_or("?").to_string(),
                coords,
                centers: v["centers"].as_array().cloned().unwrap_or_default(),
            }
        })
        .collect()
}

/// 按真值给的邻居序造一个 `Center`。三配位的最后一格是孤对。
fn center_from_truth(c: &serde_json::Value) -> (Center, Vec<usize>) {
    let nbrs: Vec<usize> = c["nbrs"]
        .as_array()
        .expect("有 nbrs")
        .iter()
        .map(|x| x.as_u64().expect("邻居下标") as usize)
        .collect();
    assert!(
        nbrs.len() == 3 || nbrs.len() == 4,
        "真值里的中心有 {} 个邻居 —— 只应当是 3(带孤对)或 4",
        nbrs.len()
    );
    let mut ligands = [Center::IMPLICIT; 4];
    for (k, &n) in nbrs.iter().enumerate() {
        ligands[k] = u32::try_from(n).expect("下标放得进 u32");
    }
    (
        Center {
            atom: u32::try_from(c["atom"].as_u64().expect("有 atom")).expect("放得进 u32"),
            ligands,
            sign: c["sign"].as_f64().expect("有 sign"),
        },
        nbrs,
    )
}

const BASELINES: [&str; 2] = ["smoke.chirality.jsonl", "smoke.lonepair.jsonl"];

#[test]
fn 真值的体积列与产品的公式在同一坐标上逐个相同() {
    let mut n = 0usize;
    let mut worst = 0f64;
    let mut bad = Vec::new();
    for name in BASELINES {
        for r in read(name) {
            for c in &r.centers {
                let (center, _) = center_from_truth(c);
                let want = c["vol"].as_f64().expect("有 vol");
                let got = chiral::center_volume(&r.coords, &center);
                n += 1;
                worst = worst.max((got - want).abs());
                if (got - want).abs() > VOL_TOL {
                    bad.push(format!(
                        "  {} 原子 {}:真值 {want}, 产品复算 {got}",
                        r.smiles, center.atom
                    ));
                }
            }
        }
    }
    assert!(n > 0, "一个中心都没读到 —— 基准空了或者字段改名了");
    assert!(
        bad.is_empty(),
        "{} 个中心的体积与真值对不上(容差 {VOL_TOL:.0e}):\n{}\n\n\
         同一份坐标、同一个配体序、同一个定义,两个独立实现必须给出同一个数。\n\
         对不上通常意味着基点或配体三元组变了 —— 那会让“只比号”的判据继续全绿。",
        bad.len(),
        bad.join("\n")
    );
    println!("真值体积复算:{n} 个中心,最大差 {worst:.3e}(容差 {VOL_TOL:.0e})");
}

#[test]
fn 真值的中心都落在配体四面体里面() {
    let (mut n4, mut outside) = (0usize, Vec::new());
    for name in BASELINES {
        for r in read(name) {
            for c in &r.centers {
                let (center, nbrs) = center_from_truth(c);
                if nbrs.len() != 4 {
                    continue; // 三配位那一格是孤对,四配体行列式无从谈起
                }
                n4 += 1;
                let vol = c["vol"].as_f64().expect("有 vol");
                let ligand = c["vol_ligand"].as_f64().expect("有 vol_ligand");
                if vol == 0.0 || ligand == 0.0 || vol.signum() == ligand.signum() {
                    outside.push(format!(
                        "  {} 原子 {}:vol={vol}, vol_ligand={ligand}",
                        r.smiles, center.atom
                    ));
                }
            }
        }
    }
    assert!(n4 > 0, "一个四配位中心都没有 —— 这条判据空过了");
    assert!(
        outside.is_empty(),
        "{} 个中心的两种体积同号 —— 中心原子在配体四面体**外面**(伞形翻转):\n{}\n\n\
         正四面体上 V_配体 = −4·V_中心,恒反号。同号时 `sign` 说的构型\n\
         与分子实际的构型对不上,拿它当真值没有意义。",
        outside.len(),
        outside.join("\n")
    );
    println!("配体四面体:{n4} 个四配位中心,全部反号(中心在里面)");
}

#[test]
fn 真值的号不是压在共面的中心上() {
    let mut worst: Option<(f64, String)> = None;
    let mut n = 0usize;
    for name in BASELINES {
        for r in read(name) {
            for c in &r.centers {
                let (center, nbrs) = center_from_truth(c);
                let o = r.coords[center.atom as usize];
                let d = |k: usize| {
                    let p = r.coords[nbrs[k]];
                    [p[0] - o[0], p[1] - o[1], p[2] - o[2]]
                };
                let norm = |v: [f64; 3]| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
                let scale = norm(d(0)) * norm(d(1)) * norm(d(2));
                assert!(
                    scale > 0.0,
                    "{} 原子 {} 有一根零长度的键",
                    r.smiles,
                    center.atom
                );
                let flat = chiral::center_volume(&r.coords, &center).abs() / scale;
                n += 1;
                if worst.as_ref().is_none_or(|(w, _)| flat < *w) {
                    worst = Some((flat, format!("{} 原子 {}", r.smiles, center.atom)));
                }
            }
        }
    }
    let (flat, who) = worst.expect("至少有一个中心");
    assert!(
        flat > MIN_FLATNESS,
        "最平的真值中心是 {who},无量纲平度 {flat:.4} ≤ {MIN_FLATNESS} —— \
         它的 sign 是掷硬币,拿它当真值没有意义。\n\
         (正四面体是 0.7698;这两份基准实测中位数 0.75 / 0.89)"
    );
    println!("真值平度:{n} 个中心,最平的是 {who},{flat:.4}(下限 {MIN_FLATNESS})");
}
