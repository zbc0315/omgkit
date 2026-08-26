//! 逐条读一个 SDF,每条打印一行,交给判官与外部实现比。
//!
//! ```shell
//! python3 harness/check_sdf.py --write /tmp/data.sdf harness/corpus/large.smi
//! cargo run -q -p omgkit-io --release --example read_sdf -- /tmp/data.sdf > /tmp/ours.txt
//! python3 harness/check_sdf.py --compare /tmp/data.sdf /tmp/ours.txt
//! ```
//!
//! 每条一行,制表符分隔:`第几条`、规范 SMILES(或 `<读不了:…>`)、
//! 数据字段(JSON 数组,每项是 `[名字, 值]`)。
//!
//! **数据字段必须用 JSON,不能自己拼转义。** 头一版把换行写成 `\n`、字段之间
//! 用 `\x1f` 隔开,当场就错了:语料里有条 SMILES 是 `[H]/N=c/1\nc[nH]s1`,
//! 里面那两个字符**本来就是**反斜杠加 n。判官解码时把它当成换行,报了两条
//! "数据字段不同" —— 而读取器一点毛病没有,是判据自己的传输把值改了。
//!
//! **一条都不许跳过。** 读不了的那条也要占一行 —— 跳过会让判官那边的条数
//! 对不上,而"对不上"正是它要抓的东西之一。
fn main() {
    let path = std::env::args().nth(1).expect("输入文件");
    let text = std::fs::read_to_string(&path).expect("读输入");
    for (i, rec) in omgkit_io::molblock::read_sdf(&text).enumerate() {
        let (smi, data) = match rec {
            Err(e) => (format!("<读不了:{e}>"), String::new()),
            Ok(rec) => {
                let mut m = rec.block.mol;
                let smi = match omgkit_chem::pipeline::sanitize(&mut m) {
                    Err(e) => format!("<净化不了:{e}>"),
                    Ok(()) => {
                        // 与 `read_molblock` 同一条顺序:净化之后才打立体标记。
                        let _ = omgkit_io::stereo::assign_chirality_2d(
                            &mut m,
                            &rec.block.coords,
                            &rec.block.wedges,
                        );
                        let _ = omgkit_io::stereo::assign_bond_stereo_2d(
                            &mut m,
                            &rec.block.coords,
                            &rec.block.unknown_stereo,
                        );
                        // 三维那两个:四个 `assign_*` 各自认出维数,不合的那一对返回 0。
                        // 由它们自己判、而不是在这里 `if is_3d`,是为了让"什么时候读得出立体"
                        // 只有一个住处。
                        let _ = omgkit_io::stereo::assign_chirality_3d(&mut m, &rec.block.coords);
                        let _ = omgkit_io::stereo::assign_bond_stereo_3d(
                            &mut m,
                            &rec.block.coords,
                            &rec.block.unknown_stereo,
                        );
                        omgkit_io::canon::canonical_smiles(&m).smiles
                    }
                };
                let data = serde_json::to_string(&rec.data).expect("字段转 JSON");
                (smi, data)
            }
        };
        println!("{i}\t{smi}\t{data}");
    }
}
