//! 出三维图看效果。**四套样式各一份。**
//!
//! ```shell
//! cargo run -p omgkit-depict --features raster --example draw3d -- <输出目录> 名字=SMILES ...
//! ```
//!
//! 每个分子每套样式写两个文件:`<名字>.<样式>.svg` 与 `<名字>.<样式>.png`,
//! 命名与二维那个 `draw` 一致,`docs/figures/make_figures.py` 拿它出文档配图。
//!
//! # `--scale=<磅每埃>`:把四套样式压到同一个比例尺
//!
//! 四套样式各带各的默认比例尺(空间填充 24,其余 36),那是**单独出一张图**时
//! 的合理默认 —— 球取满范德华半径,不缩小画布会大一圈。
//!
//! 但**并排比就成了误导**:同一个分子,空间填充那格看着比球棍小一圈,而分子
//! 并没有变。实测拿吗啡摆四格,第一眼像是换了个视角(其实四格每个原子的深度
//! 逐位相同)。所以要并排的图一律显式压到同一个比例尺。
//!
//! 判据能守住球心、半径、颜色、画序,守不住"这张图看着像不像那个分子" ——
//! 那要人眼看。

use omgkit_depict::three::{self, Style3D};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let dir = std::path::PathBuf::from(args.next().unwrap_or_else(|| ".".into()));
    let rest: Vec<String> = args.collect();
    let mut scale: Option<f64> = None;
    let mut mol_args: Vec<String> = Vec::new();
    for a in rest {
        if let Some(v) = a.strip_prefix("--scale=") {
            scale = Some(v.parse().expect("--scale 后面要跟一个数(磅每埃)"));
        } else {
            mol_args.push(a);
        }
    }
    let mols: Vec<(String, String)> = mol_args
        .into_iter()
        .map(|a| match a.split_once('=') {
            // 不合式的参数**不能默默丢掉**:丢掉之后这一批就少画一个分子,
            // 而排版脚本读不到文件才报错,报的还是别的错。
            None => panic!("参数要写成 `名字=SMILES`,收到的是 {a:?}"),
            Some(("", s)) => panic!("名字不能为空:{s:?}"),
            Some((n, s)) => (n.to_string(), s.to_string()),
        })
        .collect();
    if mols.is_empty() {
        return Err("用法:draw3d <输出目录> 名字=SMILES ...".into());
    }
    std::fs::create_dir_all(&dir)?;
    let style2d = omgkit_depict::style::Style::ACS_1996;

    for (name, smi) in &mols {
        let mut mol = omgkit_io::smiles::parse(smi)?;
        let conf = omgkit_conf::pipeline::conformer_for(&mut mol).map_err(|e| format!("{e}"))?;
        for style in &Style3D::ALL {
            let mut style = style.clone();
            if let Some(s) = scale {
                style.scale_pt_per_a = s;
            }
            let style = &style;
            let d = three::depict(&mol, &conf.coords, style)?;
            // **文件名要自己拼,不能用 `with_extension`。** 样式名里带个连字符
            // (`space-filling`),`Path` 把最后一个点后面的当扩展名 ——
            // `aspirin.space-filling` 加 `.svg` 会得到 `aspirin.svg`,四套样式
            // 互相覆盖,而程序一声不吭。
            let stem = format!("{name}.{}", style.name);
            std::fs::write(
                dir.join(format!("{stem}.svg")),
                omgkit_depict::svg::to_svg(&d.scene, &style2d),
            )?;
            std::fs::write(
                dir.join(format!("{stem}.png")),
                omgkit_depict::raster::to_png(&d.scene, &style2d, 300.0 / 72.0)?,
            )?;
            // 每行末尾报诊断,好让排版脚本把"画得不干净"的挑出来 —— 与二维那个
            // `draw` 是同一个约定。三维目前只有一档:视角定不定得下来。
            println!(
                "{name}\t{}\t{:.0}×{:.0}pt\t图元{}\t视角退化{}",
                style.name,
                d.scene.width,
                d.scene.height,
                d.scene.items.len(),
                u8::from(d.view.degenerate)
            );
        }
    }
    Ok(())
}
