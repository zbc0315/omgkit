//! 读一批 molblock,每条打印一行我方写回的 SMILES,交给判官与外部实现比。
//!
//! ```shell
//! python3 harness/check_molblock_read.py --write /tmp/in.sdf harness/corpus/large.smi
//! cargo run -q -p omgkit-io --release --example read_molblock -- /tmp/in.sdf > /tmp/ours.txt
//! python3 harness/check_molblock_read.py --compare /tmp/in.sdf /tmp/ours.txt
//! ```
//!
//! 输入的每条前面有一行 `>>> <行号>\t<原始 SMILES>`,后面是 molblock,
//! 以 `$$$$` 收尾 —— 与 `dump_molblock` 那边同一个装法。
//!
//! **不做净化以外的事。** 读出来的分子要净化才谈得上写 SMILES(价键、环、
//! 芳香性都在净化里),净化失败就如实报出来,不跳过 —— 跳过会让判据的分母
//! 悄悄变小。
fn main() {
    let path = std::env::args().nth(1).expect("输入文件");
    let text = std::fs::read_to_string(&path).expect("读输入");
    let mut header = String::new();
    let mut buf = String::new();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix(">>> ") {
            header = rest.to_string();
            buf.clear();
        } else if line == "$$$$" {
            let out = match omgkit_io::molblock::read_v2000(&buf) {
                Err(e) => format!("<读不了:{e}>"),
                Ok(block) => {
                    let mut m = block.mol;
                    match omgkit_chem::pipeline::sanitize(&mut m) {
                        Err(e) => format!("<净化不了:{e}>"),
                        Ok(()) => {
                            // **净化之后**才打手性标记:判一个中心要知道它有几个
                            // 隐式氢,而那一栏是净化算出来的。
                            let _ = omgkit_io::wedge::assign_chirality_2d(
                                &mut m,
                                &block.coords,
                                &block.wedges,
                            );
                            omgkit_io::canon::canonical_smiles(&m).smiles
                        }
                    }
                }
            };
            println!("{header}\t{out}");
            buf.clear();
        } else {
            buf.push_str(line);
            buf.push('\n');
        }
    }
}
