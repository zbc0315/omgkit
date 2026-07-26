//! 把 SMARTS 语料解析后写出,每行 `原文<TAB>本实现写出的串`。
//!
//! 交给 harness/check_smarts_write.py 做语义比对:两个 SMARTS 都由外部实现
//! 去匹配同一批分子,匹配到的原子集合应当一致。往返幂等只保证"这一趟没丢
//! 信息",保证不了"语义没变"。
fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("用法: dump_smarts_written <smarts.txt>");
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("读不到 {path}: {e}"));
    for line in text.lines() {
        let src = line.trim();
        if src.is_empty() || src.starts_with('#') {
            continue;
        }
        match omgkit_io::smarts::parse(src) {
            Ok(q) => println!("{src}\t{}", omgkit_io::smarts::write(&q)),
            // 语料里有故意的非法输入,解析器另有测试守着
            Err(_) => println!("{src}\t<parse-error>"),
        }
    }
}
