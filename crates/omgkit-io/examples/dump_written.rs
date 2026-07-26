//! 把分子语料解析后原样写出,每行一条,供外部实现比对立体化学是否守恒。
fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("用法: dump_written <分子.smi>");
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("读不到 {path}: {e}"));
    for tok in text.lines().filter_map(|l| l.split_whitespace().next()) {
        if tok.starts_with('#') {
            continue;
        }
        match omgkit_io::smiles::parse(tok) {
            Ok(m) => println!("{tok}\t{}", omgkit_io::smiles::write(&m).smiles),
            Err(_) => println!("{tok}\t<parse-error>"),
        }
    }
}
