//! 出图看效果。
//!
//! ```shell
//! cargo run -p omgkit-depict --features raster --example draw -- <输出目录>
//! ```
//!
//! 判据能守住"非白像素够多""格式头对",守不住"这张图画得对不对" —— 那要人眼看。
use omgkit_depict::{generate, raster, render::scene, style::Style, svg::to_svg};

/// 文件名里用的短名。**加了新规范就要在这里加一支** —— 否则它会跟别人
/// 撞同一个短名,把前一套的图悄悄覆盖掉。
fn tag_of(st: &Style) -> &'static str {
    match st.name {
        "ACS Document 1996" => "acs",
        "ChemDraw New Document" => "cd",
        other => panic!("规范 {other:?} 还没定短名"),
    }
}

fn main() {
    let mut args = std::env::args().skip(1);
    let dir = args.next().unwrap_or_else(|| ".".into());
    // 命令行上可以直接给 `名字=SMILES`,不给就用下面这组内置的
    let extra: Vec<(String, String)> = args
        .map(|a| match a.split_once('=') {
            // 不合式的参数**不能默默丢掉**:丢掉之后 `extra` 可能空掉,于是
            // 转去画内置的那一批 —— 跑得好好的,画的却不是你要的分子
            None => panic!("参数要写成 `名字=SMILES`,收到的是 {a:?}"),
            Some((n, _)) if n.contains('\t') || n.is_empty() => {
                panic!("名字不能为空、也不能含制表符(要写进 mols.tsv):{n:?}")
            }
            Some((n, s)) => (n.to_string(), s.to_string()),
        })
        .collect();
    // 挑的时候按"能考出什么"分:芳环、稠环、桥环、糖、甾体、β-内酰胺、
    // 对映体、顺反 —— 每一类都有一个代表。
    let mols = [
        ("aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
        ("caffeine", "CN1C=NC2=C1C(=O)N(C)C(=O)N2C"),
        ("naphthalene", "c1ccc2ccccc2c1"),
        ("paracetamol", "CC(=O)Nc1ccc(O)cc1"),
        ("nicotine", "CN1CCC[C@H]1c1cccnc1"),
        ("benzoic-acid", "OC(=O)c1ccccc1"),
        ("L-alanine", "C[C@H](N)C(=O)O"),
        ("D-alanine", "C[C@@H](N)C(=O)O"),
        ("trans-butene", "C/C=C/C"),
        ("ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
        ("glucose", "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O"),
        // 环内 C=C 两端各挂一个 OH —— 按邻居计数会正好抵消,考的是
        // "环上偏哪侧不看取代基"这条
        ("ascorbic-acid", "OC[C@H](O)[C@H]1OC(=O)C(O)=C1O"),
        (
            "penicillin-G",
            "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O",
        ),
        (
            "cholesterol",
            "CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C",
        ),
        // 以下三个是桥环 —— 平面上没有好解,出图会如实报退化
        //
        // 樟脑这个式子先前写成了 `[C@]1`,两个桥头碳的构型组合在三维里**根本
        // 嵌不出来** —— 平面图照画不误,谁也看不出来。桥环体系的桥头构型不是
        // 独立的,写式子时容易各写各的。
        ("camphor", "CC1(C)[C@@H]2CC[C@@]1(C)C(=O)C2"),
        (
            "morphine",
            "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5",
        ),
        (
            "atropine",
            "CN1[C@H]2CC[C@@H]1C[C@@H](C2)OC(=O)C(CO)c1ccccc1",
        ),
    ];
    let all: Vec<(String, String)> = if extra.is_empty() {
        mols.iter()
            .map(|(n, s)| ((*n).to_string(), (*s).to_string()))
            .collect()
    } else {
        extra
    };
    // 把分子表也写出来:对照脚本(compare_rdkit.py)读它,两边不必各维护一份
    // SMILES —— 那种重复迟早会漂移,而漂移之后对照的就是两个不同的分子了。
    let table: String = all.iter().map(|(n, s)| format!("{n}\t{s}\n")).collect();
    std::fs::write(format!("{dir}/mols.tsv"), table).unwrap();
    // 规范表同理:键长在 Python 那边硬编码的话,改了 `Style` 两侧的键长就
    // 不一样了,而并排图看着依旧正常 —— 那种对照会把人引到错的结论上
    let styles: String = Style::ALL
        .iter()
        .map(|st| {
            format!(
                "{}\t{}\t{}\t{}\n",
                tag_of(st),
                st.name,
                st.bond_length_pt,
                st.atom_label_pt
            )
        })
        .collect();
    std::fs::write(format!("{dir}/styles.tsv"), styles).unwrap();

    for (name, smi) in &all {
        let (name, smi) = (name.as_str(), smi.as_str());
        let mut m = omgkit_io::smiles::parse(smi).unwrap();
        omgkit_chem::pipeline::sanitize(&mut m).unwrap();
        // 顺反感知不在 sanitize 里,要单独跑 —— 见 stereo 模块文档
        omgkit_io::stereo::perceive_bond_stereo(&mut m);
        for st in &Style::ALL {
            let tag = tag_of(st);
            let d = generate(&m, st);
            let sc = scene(&m, &d, st);
            std::fs::write(format!("{dir}/{name}.{tag}.svg"), to_svg(&sc, st)).unwrap();
            std::fs::write(
                format!("{dir}/{name}.{tag}.png"),
                raster::to_png(&sc, st, 300.0 / 72.0).unwrap(),
            )
            .unwrap();
            std::fs::write(
                format!("{dir}/{name}.{tag}.jpg"),
                raster::to_jpeg(&sc, st, 300.0 / 72.0, 92).unwrap(),
            )
            .unwrap();
            // `未画手性` 要一起报 —— 立体中心在图上没画出楔形的话,读者读到的是
            // 一个构型未定的分子,而线条本身看着一点毛病没有。
            println!(
                "{name:14} {tag:3} {:>3}×{:<3}pt  退化{} 未解冲突{} 交叉{} 未画手性{}",
                sc.width.round(),
                sc.height.round(),
                d.degraded.len(),
                d.unresolved.len(),
                d.crossings.len(),
                d.unwedged.len()
            );
        }
    }
}
