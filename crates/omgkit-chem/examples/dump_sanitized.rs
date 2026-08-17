//! 把分子语料解析、**净化**后写出,每行 `原文<TAB>写出的串`。
//!
//! 净化会重排氢的存放位置(第 12 步把隐式氢挪成显式氢),写出器必须照样把它们
//! 写回去。不净化的那一档见 omgkit-io 的 `dump_written` —— 两档要分开跑,
//! 因为它们守的是不同的性质:那一档守"写出不丢立体",这一档守"净化之后写出
//! 仍然不丢氢"。
fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("用法: dump_sanitized <分子.smi>");
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("读不到 {path}: {e}"));
    for tok in text.lines().filter_map(|l| l.split_whitespace().next()) {
        if tok.starts_with('#') {
            continue;
        }
        match omgkit_io::smiles::parse(tok) {
            Ok(mut m) => {
                if omgkit_chem::sanitize(&mut m).is_err() {
                    println!("{tok}\t<sanitize-error>");
                    continue;
                }
                println!("{tok}\t{}", omgkit_io::smiles::write(&m).smiles);
            }
            Err(_) => println!("{tok}\t<parse-error>"),
        }
    }
}
