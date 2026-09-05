//! **把二维结构式导出成 jsonl**,交给 `harness/check_depict2d.py`。
//!
//! ```shell
//! cargo run -p omgkit-depict --release --example dump_depict2d -- harness/corpus/large.smi 400 > /tmp/two.jsonl
//! .venv/bin/python harness/check_depict2d.py /tmp/two.jsonl
//! ```
//!
//! 导出的是**产品真正吐出来的那段 SVG**,不是一份"给判据看的"平行表示。
//! 判官从 SVG 里把多边形和渐变读回来,芳香环有几个、各多大则由 RDKit 独立
//! 数一遍 —— 那一侧不经过本库任何一行代码。
//!
//! 每个分子每套规范导**三份**:
//!
//! | 键 | 是什么 | 判官拿它做什么 |
//! |---|---|---|
//! | `plain` | 不铺底色 | 与 `fill` 逐字节比:抠掉多边形与 `<defs>` 之后必须相同,于是"铺底色只多了底色、别的一笔没动"是判出来的,不是声称的 |
//! | `fill` | 默认白 → 浅蓝 | 几何、图层、块数都在这一份上判 |
//! | `custom` | 一组自定义色 | 颜色接没接上。默认色与自定义色写反了这一份会红 |

use omgkit_depict::{
    generate,
    render::scene,
    style::{AromaticFill, Style},
    svg::to_svg,
};

/// 判官与 Python 侧那条判据都要照着这两个颜色再算一遍,所以**写在这里一份**。
const CUSTOM: AromaticFill = AromaticFill {
    centre: [0xff, 0xfb, 0xe6],
    edge: [0xf5, 0xc2, 0x6b],
};

fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("用法:dump_depict2d <语料> [上限]");
    let limit: usize = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(usize::MAX);
    let text = std::fs::read_to_string(&path).expect("读语料");

    let mut n = 0usize;
    for line in text.lines() {
        let smi = line.split_whitespace().next().unwrap_or("");
        if smi.is_empty() || smi.starts_with('#') {
            continue;
        }
        if n >= limit {
            break;
        }
        let Ok(mut m) = omgkit_io::smiles::parse(smi) else {
            continue;
        };
        if omgkit_chem::pipeline::sanitize(&mut m).is_err() {
            continue;
        }
        // 顺反感知不在 sanitize 里,要单独跑 —— 见 stereo 模块文档。漏了的话
        // `generate` 的 debug_assert 会拦,但 release 下不拦,画出来是另一张图。
        omgkit_io::stereo::perceive_bond_stereo(&mut m);
        n += 1;

        let mut styles = serde_json::Map::new();
        for base in &Style::ALL {
            let mut variants = serde_json::Map::new();
            for (tag, fill) in [
                ("plain", None),
                ("fill", Some(AromaticFill::DEFAULT)),
                ("custom", Some(CUSTOM)),
            ] {
                let st = Style {
                    aromatic_fill: fill,
                    ..base.clone()
                };
                let d = generate(&m, &st);
                variants.insert(
                    tag.to_string(),
                    serde_json::Value::String(to_svg(&scene(&m, &d, &st), &st)),
                );
            }
            styles.insert(base.name.to_string(), serde_json::Value::Object(variants));
        }
        println!(
            "{}",
            serde_json::json!({
                "smiles": smi,
                "custom": {
                    "centre": format!("#{:02x}{:02x}{:02x}", CUSTOM.centre[0], CUSTOM.centre[1], CUSTOM.centre[2]),
                    "edge": format!("#{:02x}{:02x}{:02x}", CUSTOM.edge[0], CUSTOM.edge[1], CUSTOM.edge[2]),
                },
                "styles": styles,
            })
        );
    }
    eprintln!("导出 {n} 个分子");
}
