//! 把生成的三维构象导成 V2000 molblock,交给外部实现读回。
//!
//! ```shell
//! cargo run -q -p omgkit-conf --release --example dump_molblock3d -- \
//!     harness/corpus/large.smi > /tmp/blocks3d.sdf
//! python3 harness/check_molblock.py /tmp/blocks3d.sdf
//! ```
//!
//! # 与 `dump_conformers` 有什么不同
//!
//! 那一条把原子表与坐标导成 JSONL,判官在 Python 里**自己拼**一个 RDKit 分子。
//! 那条路把**文件格式本身**整个绕开了:计数行、原子块的价键字段、`M CHG` /
//! `M ISO` / `M RAD`、键块的字段宽度 —— 全都没有人读过。
//!
//! 这一条走的是文件:我们写,别人读。写错一个字段,读的一方拿到的就是另一个
//! 分子,而这正是 `.mol` / `.sdf` 交出去之后真实会发生的事。
//!
//! 只导**带立体标记**的分子:立体是最容易在格式里丢掉的一档,而没有立体的
//! 分子读回来对不对,`check_write.py` 那条已经在 SMILES 上守着了。
//!
//! # 导的是**真 SDF**,不是自造的包装
//!
//! 先前每条前面加一行 `>>> 行号\t原串`、后面手写 `$$$$`。那样判官读的是我方
//! 自造的格式,于是 `$$$$` 摆在哪、数据字段怎么写,**外部实现一次都没读过**。
//! 现在走 `write_sdf_record`:行号与原串是两个数据字段,记录边界由写出器负责,
//! 判官拿 `ForwardSDMolSupplier` 当普通 SDF 读。
use omgkit_io::molblock::{write_sdf_record, Record};

fn main() {
    let path = std::env::args().nth(1).expect("语料");
    let text = std::fs::read_to_string(&path).expect("读语料");
    let mut n = 0usize;
    for (lineno, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let smi = line.split('\t').next().unwrap_or("").trim();
        if !smi.contains('@') && !smi.contains('/') && !smi.contains('\\') {
            continue;
        }
        let Ok(mut mol) = omgkit_io::smiles::parse(smi) else {
            continue;
        };
        let Ok(conf) = omgkit_conf::pipeline::conformer_for(&mut mol) else {
            continue;
        };
        // **先凯库勒化。** molblock 没有"芳香键"这回事,留着它写出去要么歧义、
        // 要么被写成单键(噻吩读回来成四氢噻吩)。写出器现在会拒绝,这里给它
        // 一份凯库勒化之后的键级。
        let mut kek = mol.clone();
        if omgkit_chem::kekulize(&mut kek).is_err() {
            continue;
        }
        let orders: Vec<_> = kek.bonds().iter().map(|b| b.order).collect();
        // 三维构象**不写楔形**:楔形是二维图上的记号,写在三维坐标旁边等于给出
        // 两个可能互相矛盾的说法。立体从坐标本身读。
        //
        // 但"顺反未知"要标:嵌出来的构象给每根双键一个确定的二面角,作者没写
        // 顺反的那些键不标交叉,读回来就成了"作者说是顺式"。
        let unknown = omgkit_io::stereo::unspecified_cis_trans(&mol);
        let rec = Record {
            title: "",
            coords: &conf.coords,
            wedges: &[],
            orders: &orders,
            unknown_stereo: &unknown,
        };
        let line_no = lineno.to_string();
        let Ok(block) = write_sdf_record(&mol, &rec, &[("行号", &line_no), ("原串", smi)])
        else {
            continue;
        };
        print!("{block}");
        n += 1;
    }
    eprintln!("导出 {n} 个分子");
}
