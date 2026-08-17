//! 把语料逐行解析再写出,输出 `原始<TAB>写出` 两列。
//!
//! 往返测试只用自家解析器验证写出,原理上有一种漏网可能:解析与写出**共享
//! 同一个误解**,互为逆运算却都偏离了 SMILES 的语义。手性尤其危险 ——
//! 标记写反了,原子数、键集合、连通性全都对,只有分子是镜像的。
//!
//! 这个程序把写出结果交给外部实现去判:两边各自规范化,再比字符串。
//! 加 `--canonical` 则写出规范 SMILES(经规范化排序)。
//!
//! 用法见 `harness/README.md`。

use std::io::{BufWriter, Write};

fn main() {
    let mut args = std::env::args().skip(1);
    let Some(path) = args.next() else {
        eprintln!("用法: write_smiles <语料.smi> [--canonical]");
        std::process::exit(2);
    };
    let canonical = args.any(|a| a == "--canonical");
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| {
        eprintln!("读不到 {path}: {e}");
        std::process::exit(2);
    });

    let stdout = std::io::stdout();
    let mut out = BufWriter::new(stdout.lock());
    let (mut ok, mut failed) = (0usize, 0usize);

    for line in text.lines() {
        let smi = line.split_whitespace().next().unwrap_or("");
        if smi.is_empty() || smi.starts_with('#') {
            continue;
        }
        match omgkit_io::smiles::parse(smi) {
            Ok(mol) => {
                let w = if canonical {
                    omgkit_io::canon::canonical_smiles(&mol)
                } else {
                    omgkit_io::smiles::write(&mol)
                };
                let _ = writeln!(out, "{smi}\t{}", w.smiles);
                ok += 1;
            }
            // 解析失败的行不参与比对 —— 解析器的报错另有测试守着
            Err(_) => failed += 1,
        }
    }
    let _ = out.flush();
    eprintln!("写出 {ok} 条,跳过(解析失败){failed} 条");
}
