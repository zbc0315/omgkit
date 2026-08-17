//! 把分子语料解析、净化后写出**规范** SMILES,每行 `原文<TAB>规范串`。
//!
//! 与 `dump_sanitized` 的分工:那一档走 `smiles::write`(按存储顺序),
//! 这一档走 `canon::canonical_smiles`。两条路径不同 —— 规范化那条会打破对称、
//! 枚举取最小,还会**抹掉不携带信息的立体标记**。抹错就把两个分子塌成一个,
//! 所以它需要自己的判据,不能靠另一档代劳。
fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("用法: dump_canonical <分子.smi>");
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
                omgkit_io::stereo::perceive_bond_stereo(&mut m);
                println!("{tok}\t{}", omgkit_io::canon::canonical_smiles(&m).smiles);
            }
            Err(_) => println!("{tok}\t<parse-error>"),
        }
    }
}
