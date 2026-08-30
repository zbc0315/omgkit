//! **把三维图连同它的输入一起导出成 jsonl**,交给 `harness/check_depict3d.py`。
//!
//! ```shell
//! cargo run -p omgkit-depict --release --example dump_depict3d -- harness/corpus/large.smi 400 > /tmp/three.jsonl
//! .venv/bin/python harness/check_depict3d.py /tmp/three.jsonl
//! ```
//!
//! 导出的是**产品真正吐出来的那段 SVG**,不是一份"给判据看的"平行表示。
//! 判官从 SVG 里把圆和线读回来,再拿坐标、旋转矩阵、元素表(它自己那份)
//! 独立算一遍该是什么样。少了这一条,判据比的就是我们自己算出来的中间量。
//!
//! 同一个分子导两遍:一遍原样,一遍**把原子倒序重编号**。两份 SVG 必须逐字节
//! 相同 —— 这是本仓头号契约在三维上的样子。重编号那段与
//! `three.rs` 判据里的 `renumbered` 是同一件事,那边跑十个分子给快速反馈,
//! 这边跑全语料。

use omgkit_core::MolBuilder;
use omgkit_depict::three::{self, Style3D};

/// 原子倒序重排,坐标跟着走。见模块文档。
fn renumbered(mol: &MolBuilder, coords: &[[f64; 3]]) -> (MolBuilder, Vec<[f64; 3]>) {
    let n = mol.num_atoms();
    let mut out = MolBuilder::with_capacity(n, mol.num_bonds());
    for a in mol.atoms().iter().rev() {
        out.add_atom_data(*a);
    }
    let map = |old: u32| n as u32 - 1 - old;
    for b in mol.bonds() {
        let mut nb = *b;
        nb.begin = map(b.begin);
        nb.end = map(b.end);
        out.add_bond_data(nb).expect("重编号不该造出坏键");
    }
    let mut c = vec![[0.0; 3]; n];
    for (i, p) in coords.iter().enumerate() {
        c[n - 1 - i] = *p;
    }
    (out, c)
}

fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("用法:dump_depict3d <语料> [上限]");
    let limit: usize = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(usize::MAX);
    let text = std::fs::read_to_string(&path).expect("读语料");
    let style2d = omgkit_depict::style::Style::ACS_1996;

    let mut n = 0usize;
    for line in text.lines() {
        let smi = line.split_whitespace().next().unwrap_or("");
        if smi.is_empty() || smi.starts_with('#') {
            continue;
        }
        if n >= limit {
            break;
        }
        let Ok(mut mol) = omgkit_io::smiles::parse(smi) else {
            continue;
        };
        let Ok(conf) = omgkit_conf::pipeline::conformer_for(&mut mol) else {
            continue;
        };
        n += 1;

        let (rmol, rcoords) = renumbered(&mol, &conf.coords);
        let mut styles = serde_json::Map::new();
        for style in &Style3D::ALL {
            let d = three::depict(&mol, &conf.coords, style).expect("画不出来");
            let r = three::depict(&rmol, &rcoords, style).expect("画不出来");
            styles.insert(
                style.name.to_string(),
                serde_json::json!({
                    "rot": d.view.rot,
                    "centre": d.view.centre,
                    "degenerate": d.view.degenerate,
                    "width": d.scene.width,
                    "height": d.scene.height,
                    "ball_vdw_frac": style.ball_vdw_frac,
                    "stick_radius_a": style.stick_radius_a,
                    "spacing_a": style.multiple_bond_spacing_a,
                    "scale": style.scale_pt_per_a,
                    "placed": d.placed.iter()
                        .map(|p| [p.at.x, p.at.y, p.radius, p.depth])
                        .collect::<Vec<_>>(),
                    "svg": omgkit_depict::svg::to_svg(&d.scene, &style2d),
                    "svg_renumbered": omgkit_depict::svg::to_svg(&r.scene, &style2d),
                }),
            );
        }

        println!(
            "{}",
            serde_json::json!({
                "smiles": smi,
                "z": mol.atoms().iter().map(|a| a.atomic_num).collect::<Vec<_>>(),
                "bonds": mol.bonds().iter()
                    .map(|b| serde_json::json!([b.begin, b.end, format!("{:?}", b.order)]))
                    .collect::<Vec<_>>(),
                "coords": conf.coords,
                "styles": styles,
            })
        );
    }
    eprintln!("导出 {n} 个分子");
}
