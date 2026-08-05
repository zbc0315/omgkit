//! 把画出来的图导成 V2000 molblock,交给外部实现反读构型。
//!
//! ```shell
//! cargo run -p omgkit-depict --release --example dump_molblock -- harness/corpus/large.smi \
//!     > /tmp/blocks.txt
//! python3 harness/check_wedge_readback.py /tmp/blocks.txt
//! ```
//!
//! # 为什么要导出去让别人读
//!
//! [`assign_wedges`](omgkit_depict::stereo::assign_wedges) 是"试 `Up`/`Down`,
//! 取**反读回来**等于图里那个标记的那一个"构造出来的,而反读用的就是
//! [`read_chirality`](omgkit_depict::stereo::read_chirality) —— 两者共谋,拿它们
//! 的往返去检验是**空过的**。函数自己的文档就写着这一点。
//!
//! 真正要问的是另一个问题:**别人照着这张图读,读出来是不是同一个分子。**
//! 那就必须把图交出去。molblock 是最没有歧义的载体:2D 坐标 + 每根键的立体标记
//! (1 = 实楔,6 = 虚楔,**窄端是键的第一个原子**),外部实现照着它指派手性。
//!
//! # 输出格式
//!
//! 每个分子一段:
//!
//! ```text
//! >>> <语料行号>\t<SMILES>
//! #unwedged <逗号分隔的原子下标>
//! <molblock>
//! $$$$
//! ```
//!
//! `#unwedged` 是**图上没能画出构型的中心**,判官据此把"如实报过的"与"画错了
//! 还说自己对的"分开 —— 两者不是一回事。
//!
//! 没有四面体中心的分子直接跳过。
use omgkit_core::BondOrder;
use omgkit_depict::{generate, render::drawn_orders, stereo::Wedge, style::Style};

fn main() {
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "harness/corpus/large.smi".into());
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("读不了 {path}:{e}"));
    let style = &Style::ACS_1996;

    for (lineno, line) in text.lines().enumerate() {
        let smi = line.split_whitespace().next().unwrap_or("");
        if smi.is_empty() {
            continue;
        }
        let Ok(mut m) = omgkit_io::smiles::parse(smi) else {
            continue;
        };
        if omgkit_chem::pipeline::sanitize(&mut m).is_err() || m.num_atoms() < 2 {
            continue;
        }
        omgkit_io::stereo::perceive_bond_stereo(&mut m);
        let d = generate(&m, style);
        if d.wedges.iter().all(|w| w.narrow().is_none()) && d.unwedged.is_empty() {
            continue;
        }
        let orders = drawn_orders(&m);

        println!(">>> {lineno}\t{smi}");
        println!(
            "#unwedged {}",
            d.unwedged
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(",")
        );
        println!();
        println!("omgkit");
        println!();
        println!(
            "{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000",
            m.num_atoms(),
            m.num_bonds()
        );
        for (i, a) in m.atoms().iter().enumerate() {
            let p = d.coords[i];
            let sym = omgkit_core::element::by_atomic_num(a.atomic_num).map_or("*", |e| e.symbol);
            println!(
                "{:>10.4}{:>10.4}{:>10.4} {sym:<3} 0  0  0  0  0  0  0  0  0  0  0  0",
                p.x, p.y, 0.0
            );
        }
        for (bi, b) in m.bonds().iter().enumerate() {
            // **窄端必须写成键的第一个原子** —— molblock 的立体标记是这么定义的。
            // 写反了,楔形描述的就是另一头那个原子的构型。
            let (first, second, code) = match d.wedges[bi] {
                Wedge::Up { narrow } => {
                    let o = if narrow == b.begin { b.end } else { b.begin };
                    (narrow, o, 1)
                }
                Wedge::Down { narrow } => {
                    let o = if narrow == b.begin { b.end } else { b.begin };
                    (narrow, o, 6)
                }
                Wedge::None => (b.begin, b.end, 0),
            };
            // 芳香键要按凯库勒式写 —— molblock 的 4 号键级各家读法不一
            let ord = match orders[bi] {
                BondOrder::Double => 2,
                BondOrder::Triple => 3,
                _ => 1,
            };
            println!(
                "{:>3}{:>3}{ord:>3}{code:>3}  0  0  0",
                first + 1,
                second + 1
            );
        }
        for (i, a) in m.atoms().iter().enumerate() {
            if a.formal_charge != 0 {
                println!("M  CHG  1{:>4}{:>4}", i + 1, a.formal_charge);
            }
            if a.isotope != 0 {
                println!("M  ISO  1{:>4}{:>4}", i + 1, a.isotope);
            }
        }
        println!("M  END");
        println!("$$$$");
    }
}
