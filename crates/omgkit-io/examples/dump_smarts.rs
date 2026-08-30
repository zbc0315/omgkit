//! 逐行解析 SMARTS 语料,输出 `SMARTS<TAB>原子数<TAB>键数`,失败输出 `ERR<TAB>原因`。
//!
//! 供与外部实现对拍用。(写出那一档的判官是 `harness/check_smarts_write.py`,
//! 它吃的是 `dump_smarts_written` 的输出,不是这个例子的。)

use std::io::{BufWriter, Write};

fn main() {
    let Some(path) = std::env::args().nth(1) else {
        eprintln!("用法: dump_smarts <语料.txt>");
        std::process::exit(2);
    };
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| {
        eprintln!("读不到 {path}: {e}");
        std::process::exit(2);
    });

    let stdout = std::io::stdout();
    let mut out = BufWriter::new(stdout.lock());
    let (mut ok, mut err) = (0usize, 0usize);

    for line in text.lines() {
        let s = line.trim();
        if s.is_empty() || s.starts_with('#') {
            continue;
        }
        match omgkit_io::smarts::parse(s) {
            Ok(q) => {
                let _ = writeln!(out, "{s}\t{}\t{}", q.num_atoms(), q.num_bonds());
                ok += 1;
            }
            Err(e) => {
                let _ = writeln!(out, "{s}\tERR\t{}", e.kind);
                err += 1;
            }
        }
    }
    let _ = out.flush();
    eprintln!("解析成功 {ok} 条,失败 {err} 条");
}
