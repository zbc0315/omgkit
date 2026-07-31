//! 绘图规范 —— **同一份参数同时喂给布局和渲染**。
//!
//! # 为什么规范要参与布局,而不只是渲染
//!
//! 直觉上"排坐标"和"画出来"是两件事(Mayfield 在 RDKit UGM 2016 上正是这么
//! 切的),几何骨架层面也确实如此:苯环是正六边形,与用什么字号无关。
//!
//! 但**"算不算挤在一起"是规范相关的**。原子标签要占地方,而标签尺寸与键长的
//! 比例随规范变:
//!
//! | | ACS 1996 | ChemDraw 默认 |
//! |---|---|---|
//! | 键长 | 14.4 pt | 30 pt |
//! | 原子标签 | 10 pt | 10 pt |
//! | **标签占一个键长的** | **69%** | **33%** |
//!
//! 字号一样、键长差 2.08 倍 —— 按默认规范排得开的图,换成 ACS 就会挤上。所以
//! 碰撞阈值必须来自规范,布局不能对规范一无所知。
//!
//! # 一份配置,两端共用
//!
//! 这里刻意**不**分成"布局配置"和"渲染配置"两个结构。两套平行配置迟早漂移,
//! 而漂移的后果恰恰是最难发现的那种:按 A 的间距排版、按 B 的字号绘制,图上
//! 挤成一团却没有任何一步报错。
//!
//! [`Depiction`](crate::Depiction) 会记下产生它的规范指纹,拿另一套规范去渲染
//! 时可以被查出来。
//!
//! # 数值出处
//!
//! 全部取自 ChemDraw 17.1 用户手册第 4 章 "Preferences and Settings" 的样式表
//! 清单(手册共收录 11 套)。选这两套的理由:ACS 1996 是发表用的事实标准,
//! New Document 是 ChemDraw 打开即用的默认值。其余 9 套(Wiley、Synthesis、
//! Adv. Synth. Catal. 等)键长都在 17–20 pt,夹在这两者之间。

/// 一套绘图规范。
///
/// 长度一律以**磅(pt)** 为单位,与 ChemDraw 的文档设置对齐,便于逐项核对。
/// 布局内部用的是"键长 = 1"的无量纲坐标,换算由 [`Style::bond_length_pt`] 负责。
#[derive(Debug, Clone, PartialEq)]
pub struct Style {
    /// 规范名,出现在 [`Depiction`](crate::Depiction) 的指纹里
    pub name: &'static str,
    /// 键长(ChemDraw 的 Fixed Length)
    pub bond_length_pt: f64,
    /// 普通键的线宽(Line Width)
    pub line_width_pt: f64,
    /// 粗键与楔形键的宽度(Bold Width)。楔形末端宽度是它的 1.5 倍 —— 手册原话
    pub bold_width_pt: f64,
    /// 原子标签四周的留白(Margin Width)。键画到标签边上要留出这一圈
    pub margin_width_pt: f64,
    /// 虚楔形/虚线键的横线间距(Hash Spacing)
    pub hash_spacing_pt: f64,
    /// 双键两条线的间距,**占键长的百分比**(Bond Spacing)
    pub bond_spacing_pct: f64,
    /// 链的键角(Chain Angle)
    pub chain_angle_deg: f64,
    /// 原子标签字号
    pub atom_label_pt: f64,
    /// 图注字号
    pub caption_pt: f64,
    /// 字体族。
    ///
    /// 手册里这一项是以图片呈现的,文本抽不出来,所以这里**不假称是从手册读到
    /// 的**:取无衬线族,与两套规范的实际观感一致。要精确复刻请自行覆盖。
    pub font_family: &'static str,
}

impl Style {
    /// **ACS Document 1996** —— ACS 期刊双栏版式,发表用的事实标准。
    ///
    /// 数值取自 ChemDraw 17.1 手册。
    pub const ACS_1996: Style = Style {
        name: "ACS Document 1996",
        bond_length_pt: 14.4,
        line_width_pt: 0.6,
        bold_width_pt: 2.0,
        margin_width_pt: 1.6,
        hash_spacing_pt: 2.5,
        bond_spacing_pct: 18.0,
        chain_angle_deg: 120.0,
        atom_label_pt: 10.0,
        caption_pt: 10.0,
        font_family: "Arial, Helvetica, sans-serif",
    };

    /// **New Document** —— ChemDraw 新建文档的默认值。
    ///
    /// 键长是 ACS 的 2.08 倍,于是同样 10 pt 的标签只占三分之一个键长 ——
    /// 这套排得开的图换到 ACS 未必排得开。
    pub const CHEMDRAW_DEFAULT: Style = Style {
        name: "ChemDraw New Document",
        bond_length_pt: 30.0,
        line_width_pt: 1.0,
        bold_width_pt: 2.0,
        margin_width_pt: 2.0,
        hash_spacing_pt: 2.7,
        bond_spacing_pct: 12.0,
        chain_angle_deg: 120.0,
        atom_label_pt: 10.0,
        caption_pt: 12.0,
        font_family: "Arial, Helvetica, sans-serif",
    };

    /// 本库内置的全部规范。
    pub const ALL: [Style; 2] = [Style::ACS_1996, Style::CHEMDRAW_DEFAULT];

    /// 原子标签的字号,换算成**键长为单位**。
    ///
    /// 布局全程用无量纲坐标,碰撞判定要的正是这个比值 —— ACS 是 0.69,
    /// ChemDraw 默认是 0.33。
    #[must_use]
    pub fn label_size(&self) -> f64 {
        self.atom_label_pt / self.bond_length_pt
    }

    /// 标签四周留白,换算成键长为单位。
    #[must_use]
    pub fn margin(&self) -> f64 {
        self.margin_width_pt / self.bond_length_pt
    }

    /// 双键第二条线与主线的间距,键长为单位。
    #[must_use]
    pub fn bond_spacing(&self) -> f64 {
        self.bond_spacing_pct / 100.0
    }

    /// 两个**没有成键**的原子,中心至少要离这么远才不算撞上。
    ///
    /// 取标签高度加两侧留白 —— 两个标签各占一半高度,中间还要留出空隙。带标签
    /// 的原子(杂原子、带电原子)才真的占这么大;裸碳可以更近,那由调用方按
    /// 实际标签尺寸细化,这里给的是**上界**。
    #[must_use]
    pub fn clash_distance(&self) -> f64 {
        (self.label_size() + 2.0 * self.margin()).max(0.35)
    }

    /// 影响**布局结果**的那些参数的指纹。
    ///
    /// 字号、线宽这些只影响渲染的项不计入 —— 换个线宽不该让坐标失效。计入的是
    /// 真正会改变坐标的:碰撞阈值与链角。
    ///
    /// [`Depiction`](crate::Depiction) 存下它,于是"用 A 规范排版、拿 B 规范
    /// 渲染"这种静默错配可以被查出来。
    #[must_use]
    pub fn layout_fingerprint(&self) -> u64 {
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        for v in [self.clash_distance(), self.chain_angle_deg] {
            for b in v.to_bits().to_le_bytes() {
                h ^= u64::from(b);
                h = h.wrapping_mul(0x1000_0000_01b3);
            }
        }
        h
    }
}

impl Default for Style {
    fn default() -> Self {
        Style::ACS_1996
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_two_styles_really_do_differ_in_what_counts_as_a_clash() {
        // 这是"规范必须参与布局"的全部依据。两套规范的字号一样(都是 10 pt),
        // 如果碰撞阈值也一样,那把 Style 传进布局就是白传 —— 整个架构决定
        // 就落空了。这里量的正是那个差别。
        let acs = Style::ACS_1996;
        let cd = Style::CHEMDRAW_DEFAULT;
        assert_eq!(
            acs.atom_label_pt, cd.atom_label_pt,
            "两套的字号本来就该一样"
        );
        assert!(
            acs.label_size() > cd.label_size() * 1.9,
            "标签占键长的比例应当差近一倍:ACS {:.2} vs 默认 {:.2}",
            acs.label_size(),
            cd.label_size()
        );
        assert!(
            acs.clash_distance() > cd.clash_distance() * 1.5,
            "碰撞阈值必须显著不同,否则规范参与布局就没有意义"
        );
    }

    #[test]
    fn the_fingerprint_tracks_layout_and_ignores_rendering() {
        let base = Style::ACS_1996;

        // 只改渲染项 —— 指纹必须不变,否则换个线宽就要重排一遍
        let mut render_only = base.clone();
        render_only.line_width_pt = 99.0;
        render_only.hash_spacing_pt = 42.0;
        render_only.caption_pt = 7.0;
        render_only.font_family = "Comic Sans";
        assert_eq!(
            base.layout_fingerprint(),
            render_only.layout_fingerprint(),
            "只改渲染项不该让已有坐标失效"
        );

        // 改布局项 —— 指纹必须变,否则错配查不出来
        let mut layout_changed = base.clone();
        layout_changed.bond_length_pt = 30.0;
        assert_ne!(
            base.layout_fingerprint(),
            layout_changed.layout_fingerprint(),
            "改了键长(进而改了碰撞阈值)却拿到同一个指纹,错配就查不出来了"
        );

        assert_ne!(
            Style::ACS_1996.layout_fingerprint(),
            Style::CHEMDRAW_DEFAULT.layout_fingerprint(),
            "两套内置规范的布局指纹必须不同"
        );
    }

    #[test]
    fn the_numbers_match_the_chemdraw_manual() {
        // 逐项对着 ChemDraw 17.1 手册第 4 章的样式表清单核过。抄错一个数不会让
        // 任何东西报错,只会让图和期刊要求差一点点 —— 而那正是没人会去查的那种错。
        let a = Style::ACS_1996;
        assert_eq!(
            (a.bond_length_pt, a.line_width_pt, a.bold_width_pt),
            (14.4, 0.6, 2.0)
        );
        assert_eq!((a.margin_width_pt, a.hash_spacing_pt), (1.6, 2.5));
        assert_eq!(
            (a.bond_spacing_pct, a.chain_angle_deg, a.atom_label_pt),
            (18.0, 120.0, 10.0)
        );

        let c = Style::CHEMDRAW_DEFAULT;
        assert_eq!(
            (c.bond_length_pt, c.line_width_pt, c.bold_width_pt),
            (30.0, 1.0, 2.0)
        );
        assert_eq!((c.margin_width_pt, c.hash_spacing_pt), (2.0, 2.7));
        assert_eq!(
            (c.bond_spacing_pct, c.chain_angle_deg, c.atom_label_pt),
            (12.0, 120.0, 10.0)
        );
    }
}
