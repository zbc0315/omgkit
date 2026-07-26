//! 对分子语料跑双键立体感知,输出每根被标注的双键。
//!
//! 每行 `SMILES<TAB>begin,end,顺反,参照a,参照b;...`,顺反写作 CIS/TRANS。
//! 没有标注的分子也出一行(第二列为空),这样"本实现漏标"也能被判官看见。
use omgkit_core::BondStereo;

fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("用法: dump_bond_stereo <分子.smi>");
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("读不到 {path}: {e}"));
    for tok in text.lines().filter_map(|l| l.split_whitespace().next()) {
        if tok.starts_with('#') {
            continue;
        }
        let Ok(mut m) = omgkit_io::smiles::parse(tok) else {
            println!("{tok}\t<parse-error>");
            continue;
        };
        omgkit_io::stereo::perceive_bond_stereo(&mut m);
        let cells: Vec<String> = m
            .bonds()
            .iter()
            .filter(|b| b.stereo != BondStereo::None)
            .map(|b| {
                let s = match b.stereo {
                    BondStereo::Cis => "CIS",
                    BondStereo::Trans => "TRANS",
                    _ => "?",
                };
                format!(
                    "{},{},{s},{},{}",
                    b.begin, b.end, b.stereo_atoms[0], b.stereo_atoms[1]
                )
            })
            .collect();
        println!("{tok}\t{}", cells.join(";"));
    }
}
