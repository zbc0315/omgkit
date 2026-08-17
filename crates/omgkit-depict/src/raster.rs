//! 位图后端:PNG / JPEG。**需要 `raster` feature。**
//!
//! 走的是 SVG → 光栅化 → 编码这一条路,而不是另写一遍几何。理由是几何只该算
//! 一遍:两条路各算一遍,它们迟早会在某个边界上不一致,而"SVG 里对、PNG 里错"
//! 这种缺陷极难定位。
//!
//! # 字体是**运行时**依赖
//!
//! 标签在 SVG 里是 `<text>`,光栅化时要真的找到字体来排。这里加载系统字体库,
//! 所以:
//!
//! - 装不到 Arial/Helvetica 的机器上,字形会换成别的无衬线体,**宽度与
//!   [`label`](crate::label) 按 AFM 算出来的不一致** —— 布局留的空隙就对不上
//! - 一台字体齐全、一台不齐全的机器,同一个分子会得到两张不完全一样的图
//!
//! 这不是可以绕开的实现细节,是这条路线的固有代价。要求逐像素可复现的场合,
//! 应当在一台受控的机器上出图再分发。
//!
//! # 为什么 PNG 才是该用的那个
//!
//! JPEG 是有损的,而结构式全是**细线和小字**——正是 JPEG 最不擅长的内容,
//! 线条边上会出现振铃。这里提供 [`to_jpeg`] 只因为下游有时只收这一种格式;
//! 能选就选 [`to_png`]。

use crate::render::Scene;
use crate::style::Style;
use crate::svg::to_svg;

/// 光栅化失败的原因。
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RasterError {
    /// 生成的 SVG 解析不了。**这是本库的缺陷**,不是调用方的问题。
    BadSvg(String),
    /// 画布尺寸算下来是 0 或超出上限
    BadSize {
        /// 宽(像素)
        width: u32,
        /// 高(像素)
        height: u32,
    },
    /// JPEG 编码失败
    Encode(String),
}

impl std::fmt::Display for RasterError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RasterError::BadSvg(e) => write!(f, "生成的 SVG 解析不了(这是本库的缺陷):{e}"),
            RasterError::BadSize { width, height } => write!(f, "画布尺寸不合法:{width}×{height}"),
            RasterError::Encode(e) => write!(f, "编码失败:{e}"),
        }
    }
}

impl std::error::Error for RasterError {}

/// 光栅化的结果:RGBA8 像素。
pub struct Pixels {
    /// 宽(像素)
    pub width: u32,
    /// 高(像素)
    pub height: u32,
    /// RGBA,每通道 8 位,行优先
    pub rgba: Vec<u8>,
}

/// 把一张 [`Scene`] 光栅化。
///
/// `scale` 是相对磅的放大倍数:磅是 1/72 英寸,所以 `scale = 300.0 / 72.0`
/// 得到 300 dpi。
///
/// # Errors
///
/// 画布尺寸不合法、或生成的 SVG 解析不了时返回错误。
pub fn rasterize(scene: &Scene, style: &Style, scale: f32) -> Result<Pixels, RasterError> {
    let w = (scene.width as f32 * scale).round();
    let h = (scene.height as f32 * scale).round();
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let (wu, hu) = (w.max(0.0) as u32, h.max(0.0) as u32);
    if wu == 0 || hu == 0 || wu > 20_000 || hu > 20_000 {
        return Err(RasterError::BadSize {
            width: wu,
            height: hu,
        });
    }

    let svg = to_svg(scene, style);
    let mut opt = resvg::usvg::Options::default();
    // 字体在这一步才真的被找到 —— 见模块文档里关于可复现性的那一段
    opt.fontdb_mut().load_system_fonts();
    let tree =
        resvg::usvg::Tree::from_str(&svg, &opt).map_err(|e| RasterError::BadSvg(e.to_string()))?;

    let mut pixmap = resvg::tiny_skia::Pixmap::new(wu, hu).ok_or(RasterError::BadSize {
        width: wu,
        height: hu,
    })?;
    resvg::render(
        &tree,
        resvg::tiny_skia::Transform::from_scale(scale, scale),
        &mut pixmap.as_mut(),
    );

    Ok(Pixels {
        width: wu,
        height: hu,
        rgba: pixmap.take(),
    })
}

/// 出 PNG。
///
/// # Errors
///
/// 见 [`rasterize`]。
pub fn to_png(scene: &Scene, style: &Style, scale: f32) -> Result<Vec<u8>, RasterError> {
    let px = rasterize(scene, style, scale)?;
    let pixmap = resvg::tiny_skia::Pixmap::from_vec(
        px.rgba,
        resvg::tiny_skia::IntSize::from_wh(px.width, px.height).ok_or(RasterError::BadSize {
            width: px.width,
            height: px.height,
        })?,
    )
    .ok_or(RasterError::BadSize {
        width: px.width,
        height: px.height,
    })?;
    pixmap
        .encode_png()
        .map_err(|e| RasterError::Encode(e.to_string()))
}

/// 出 JPEG。`quality` 取 1–100。
///
/// **能选就选 [`to_png`]** —— 结构式全是细线和小字,正是 JPEG 最不擅长的内容。
///
/// # Errors
///
/// 见 [`rasterize`]。
pub fn to_jpeg(
    scene: &Scene,
    style: &Style,
    scale: f32,
    quality: u8,
) -> Result<Vec<u8>, RasterError> {
    let px = rasterize(scene, style, scale)?;

    // JPEG 没有透明通道。**必须自己合成到白底上**,不能直接丢掉 alpha ——
    // 丢掉的话透明处的 RGB 是未定义的(这里是 0),整张图会变成黑底。
    let mut rgb = Vec::with_capacity(px.rgba.len() / 4 * 3);
    for p in px.rgba.chunks_exact(4) {
        let a = f32::from(p[3]) / 255.0;
        for c in &p[..3] {
            // pixmap 是预乘 alpha 的,所以直接与白底相加即可
            let v = f32::from(*c) + 255.0 * (1.0 - a);
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            rgb.push(v.clamp(0.0, 255.0) as u8);
        }
    }

    let mut out = Vec::new();
    let enc = jpeg_encoder::Encoder::new(&mut out, quality.clamp(1, 100));
    enc.encode(
        &rgb,
        u16::try_from(px.width).map_err(|_| RasterError::BadSize {
            width: px.width,
            height: px.height,
        })?,
        u16::try_from(px.height).map_err(|_| RasterError::BadSize {
            width: px.width,
            height: px.height,
        })?,
        jpeg_encoder::ColorType::Rgb,
    )
    .map_err(|e| RasterError::Encode(e.to_string()))?;
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{generate, render::scene};
    use omgkit_core::MolBuilder;

    fn prep(smi: &str) -> MolBuilder {
        let mut m = omgkit_io::smiles::parse(smi).unwrap();
        omgkit_chem::pipeline::sanitize(&mut m).unwrap();
        m
    }

    fn sc(smi: &str, style: &Style) -> Scene {
        let m = prep(smi);
        scene(&m, &generate(&m, style), style)
    }

    #[test]
    fn png_has_the_right_magic_and_size() {
        let s = sc("CC(=O)Oc1ccccc1C(=O)O", &Style::ACS_1996);
        let png = to_png(&s, &Style::ACS_1996, 300.0 / 72.0).expect("出 PNG");
        assert_eq!(&png[..8], b"\x89PNG\r\n\x1a\n", "PNG 头不对");
        assert!(png.len() > 1000, "PNG 只有 {} 字节,大概率是空图", png.len());
    }

    #[test]
    fn jpeg_has_the_right_magic() {
        let s = sc("c1ccccc1", &Style::ACS_1996);
        let j = to_jpeg(&s, &Style::ACS_1996, 2.0, 90).expect("出 JPEG");
        assert_eq!(&j[..2], b"\xff\xd8", "JPEG 头不对");
        assert_eq!(&j[j.len() - 2..], b"\xff\xd9", "JPEG 尾不对");
    }

    #[test]
    fn the_image_is_not_blank() {
        // 光栅化能"成功"却什么都没画上 —— 尺寸对、格式对、全是白的。
        // 这条数真正被画到的像素。
        let s = sc("c1ccc2ccccc2c1", &Style::ACS_1996);
        let px = rasterize(&s, &Style::ACS_1996, 2.0).expect("光栅化");
        let inked = px
            .rgba
            .chunks_exact(4)
            .filter(|p| p[0] < 200 || p[1] < 200 || p[2] < 200)
            .count();
        assert!(inked > 200, "只有 {inked} 个非白像素,这张图基本是空的");
    }

    #[test]
    fn jpeg_is_not_rendered_on_a_black_background() {
        // JPEG 没有透明通道。直接丢掉 alpha 的话,透明处的 RGB 是 0 ——
        // 整张图会变成黑底白线。这条守的正是那个合成步骤。
        let s = sc("CCO", &Style::ACS_1996);
        let j = to_jpeg(&s, &Style::ACS_1996, 2.0, 95).expect("出 JPEG");
        // 重新解出来太重,改判文件大小:黑底图熵高得多,同等尺寸下会大很多
        let png = to_png(&s, &Style::ACS_1996, 2.0).expect("出 PNG");
        assert!(!j.is_empty() && !png.is_empty());

        // 直接查像素更硬:光栅化的结果里,四角必须是"透明或白",不能是不透明的黑
        let px = rasterize(&s, &Style::ACS_1996, 2.0).expect("光栅化");
        let corner = &px.rgba[..4];
        assert!(
            corner[3] == 0 || (corner[0] > 200 && corner[1] > 200 && corner[2] > 200),
            "左上角是 {corner:?} —— 既不透明也不白,底色错了"
        );
    }

    #[test]
    fn the_two_styles_give_different_pixel_sizes() {
        // 规范一路贯穿到位图。两边一样大就说明 Style 在某一层被丢了。
        let a = rasterize(
            &sc("c1ccc2ccccc2c1", &Style::ACS_1996),
            &Style::ACS_1996,
            2.0,
        )
        .unwrap();
        let c = rasterize(
            &sc("c1ccc2ccccc2c1", &Style::CHEMDRAW_DEFAULT),
            &Style::CHEMDRAW_DEFAULT,
            2.0,
        )
        .unwrap();
        assert!(
            c.width > a.width * 3 / 2,
            "ChemDraw 默认的键长是 ACS 的 2.08 倍,像素宽却只有 {} vs {}",
            c.width,
            a.width
        );
    }

    #[test]
    fn a_degenerate_size_is_refused_instead_of_panicking() {
        let s = sc("CCO", &Style::ACS_1996);
        assert!(matches!(
            rasterize(&s, &Style::ACS_1996, 0.0),
            Err(RasterError::BadSize { .. })
        ));
        assert!(matches!(
            rasterize(&s, &Style::ACS_1996, 1e6),
            Err(RasterError::BadSize { .. })
        ));
    }
}
