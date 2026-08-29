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
use omgkit_depict::{generate, render::drawn_orders, style::Style};
use omgkit_io::molblock::{write_v2000, Record};

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
        // **导出的必须是被画的那个分子。** 为画出构型补的显式氢也在里面,而
        // 楔形恰恰就打在那根 C–H 上 —— 拿原分子导出的话,那根键根本不存在,
        // 判官看到的是"没有立体信息",149 个中心全部报成"画成 None"。
        let grown = d.drawn(&m);
        let m = &*grown;
        let orders = drawn_orders(m);

        println!(">>> {lineno}\t{smi}");
        println!(
            "#unwedged {}",
            d.unwedged
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(",")
        );
        // **格式化只有一份。** 计数行、原子块的价键字段、`M CHG`/`M ISO`/`M RAD`、
        // "窄端必须是键的第一个原子"这条规则,全在 `omgkit_io::molblock` 里 ——
        // 这里只负责把二维坐标与楔形翻过去。先前这个例子自己写了一份格式化,
        // 而三维那条路要再写一份,两份必然分家。
        let coords: Vec<[f64; 3]> = d.coords.iter().map(|p| [p.x, p.y, 0.0]).collect();
        // 楔形两边是**同一个类型**(`omgkit_io::wedge::Wedge`),不用转换。
        // 作者没写顺反的双键要标成交叉双键 —— 不标的话,图上那个几何会被读成
        // 化学信息。判断的对象是**画出来的那个分子**,键下标要与写出的一致。
        let unknown = omgkit_io::stereo::unspecified_cis_trans(m);
        let rec = Record {
            title: "",
            coords: &coords,
            wedges: &d.wedges,
            orders: &orders,
            unknown_stereo: &unknown,
        };
        match write_v2000(m, &rec) {
            Ok(block) => print!("{block}"),
            Err(e) => {
                eprintln!("第 {lineno} 行 {smi} 写不出来:{e}");
                continue;
            }
        }
        println!("$$$$");
    }
}
