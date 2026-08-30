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
    /// 图注字号。
    ///
    /// **本 crate 目前不画图注,所以生产代码一次都不读它。** 留着是因为两套规范
    /// 都规定了这一项,而 `Style` 是"手册怎么说"的转录;下游自己排图注时用得上。
    ///
    /// 注意 `style.rs` 那道穷举解构闸只逼新字段"在表里露面",分不出"只进渲染"
    /// 与"哪儿都不进" —— 它给不了这个字段任何背书。
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

    /// 影响**布局结果**的那些参数的指纹。
    ///
    /// [`Depiction`](crate::Depiction) 存下它,于是"用 A 规范排版、拿 B 规范
    /// 渲染"这种静默错配可以被查出来。
    ///
    /// # 计入哪两项,是**量出来的**
    ///
    /// 逐字段扰动、比全量语料的坐标(8831 个分子)。**只有两个旋钮**:
    ///
    /// | 旋钮 | 扰动 | 坐标变了的分子 |
    /// |---|---|---:|
    /// | `chain_angle_deg` | 120° → 108° | **6687** |
    /// | `label_size()` | ×1.4(字号 ×1.4,或等价地键长 ÷1.4) | **1409** |
    /// | 同上 | ÷1.4 | 727 |
    /// | 其余**八个**字段 | ×1.4 | **全 0** |
    ///
    /// `atom_label_pt` 与 `bond_length_pt` **不是两个字段,是同一个旋钮的两头**:
    /// 布局只经 [`Style::label_size`] 读它们,而那正是二者之商。实测对称性
    /// 逐值吻合 —— 字号 ×k 与键长 ÷k 在 k = 1.36/1.37/1.38/1.4/2.0 上给出的数
    /// 分别是 1022/1351/1358/1409/4740,**两列一个不差**。
    ///
    /// **这些计数对扰动幅度极敏感,引用时必须连幅度一起写**:k = 1.36 → 1022,
    /// 1.37 → 1351,1.38 → 1358 —— 幅度差 1%,数字挪三成。它们说明的是"这个
    /// 旋钮**动不动**得了布局",不是"它有多敏感"。
    ///
    /// **验死了**:保持 `label_size()` 与 `chain_angle_deg` 不变、把**其余十个
    /// 字段全改掉**(键长与字号双双翻倍,线宽、粗宽、留白、虚线间距、双键间距、
    /// 图注字号、字体、名字全换),8831 个分子的坐标**逐点相同**;
    /// `margin_width_pt` 从 ×0.01 扫到 ×50,一张图都不动。
    ///
    /// 所以计入的正好是 `(label_size(), chain_angle_deg)` 两项 —— 不多不少。
    /// 判据 `the_fingerprint_covers_exactly_what_moves_the_layout` 逐字段钉住
    /// 这个等价关系,**并且用穷举解构逼新加的字段必须在那张表里露面**。
    ///
    /// # 先前打的是 `clash_distance()`,那是个漏
    ///
    /// 那个量是 `(label_size + 2*margin).max(0.35)` —— **布局一处都没用过它**
    /// (全仓只出现在本文件里),它把两个字段有损地揉在一起,还带个截断。
    /// 于是两套 `label_size` 差整整一倍的规范能撞出同一个指纹:实测那样一对
    /// 规范下**8831 个分子里有 787 个坐标不同,而 `matches()` 一个都没查出来**。
    /// 一套专为拦静默错配而生的机制,8.9% 的错配拦不住。
    ///
    /// `clash_distance()` 已随这次一并删掉:它没有任何调用方,而文档却写着
    /// "由调用方按实际标签尺寸细化" —— 那个调用方不存在。布局真正用的碰撞
    /// 半径是 `refine::radii`,按逐原子的标签墨迹盒算。
    ///
    /// # 常数是 FNV-1a 的那两个,别再写错一次
    ///
    /// 乘数**曾经写成 `0x1000_0000_01b3`**(17592186044851)。它既不是 FNV 的
    /// 常数,**也根本不是质数** —— 分解得 `11 × 41 × 199 × 196015399`。
    /// 与真常数 1099511628211 的关系是:两者**置位数都是 7、低 12 位都是
    /// `0x1b3`,只有最高那一位从 2⁴⁰ 挪到了 2⁴⁴** —— 手一滑多写了一格。
    /// (别说成"差整 16 倍":比值是 15.999999994,低位那 `0x1b3` 没跟着挪。)
    ///
    /// 当哈希用**确实看不出差别** —— 两者位型如此接近,雪崩表现自然相当。
    /// 但**别把理由写成"奇数乘子都散得开"**:奇数只保证乘法在 mod 2⁶⁴ 下是
    /// 双射(不丢信息),双射与散得开是两回事 —— `m = 1` 就是奇数,而它翻一个
    /// 输入位只翻一个输出位。
    ///
    /// 改它是安全的,两条理由:
    ///
    /// 1. 这个值只进 [`Depiction::style_fingerprint`](crate::Depiction) 做运行时
    ///    自比,**本仓哪儿都不落盘**(无 serde、模板表不含它、基准文件不含它);
    /// 2. 更硬的一条:**即使下游存了旧值,改动也是 fail-safe 的** —— 存了旧值
    ///    比出来是"不匹配",于是重排一遍,代价是白干,**不会画错**。危险的那个
    ///    方向(两套不同规范撞成同一个指纹)完全不受影响。
    ///
    /// 这个值**不保证跨版本稳定,别存**。
    ///
    /// 改的理由是**两处不一致**:`examples/audit.rs` 里打图元指纹用的是对的,
    /// 以后有人"对齐"这两处,很容易往错的那边改。判据
    /// `the_fingerprint_really_is_fnv_1a` 钉住它 —— 而且钉的方式要能挡住
    /// "两处一起改",见那里。
    #[must_use]
    pub fn layout_fingerprint(&self) -> u64 {
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        for v in [self.label_size(), self.chain_angle_deg] {
            for b in v.to_bits().to_le_bytes() {
                h ^= u64::from(b);
                h = h.wrapping_mul(0x0000_0100_0000_01b3);
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

    /// 指纹用的真是 FNV-1a 的那两个常数。
    ///
    /// # 光"独立重算一遍"挡不住,必须闭到外部常数上
    ///
    /// 先前这条只拿十进制字面量在测试里重算一遍,与实现比。**那挡不住这条判据
    /// 要防的那件事** —— 实现与判据里的常数**一起**改错,它照样绿(变异验过)。
    /// 而"以后有人把两处对齐、却往错的那边对"正是文档里写着的威胁模型:
    /// 判据自己就是那两处之一,一起改就一起倒。
    ///
    /// 所以再钉一环:拿 **FNV 公布的测试向量**。三方闭环 ——
    /// 实现 ←(十进制重算)→ 判据里的字面量 ←(测试向量)→ 外部常数。
    /// 那条向量同时钉住偏移基、乘数、**以及"先异或后乘"这个次序**
    /// (写成 FNV-1 会得 `0x340d_8765_a4dd_a9c2`,当场红)。
    #[test]
    fn the_fingerprint_really_is_fnv_1a() {
        // ① 这对常数确实是 FNV-1a 的 —— 拿外部公布的测试向量验
        let mut h: u64 = 14_695_981_039_346_656_037; // FNV-1a 64 位偏移基
        for b in b"foobar" {
            h ^= u64::from(*b);
            h = h.wrapping_mul(1_099_511_628_211); // FNV-1a 64 位质数
        }
        assert_eq!(
            h, 0x8594_4171_f739_67e8,
            "这对常数不是 FNV-1a 的 —— 官方测试向量 foobar 对不上"
        );

        // ② 实现用的就是这对常数,而且函数体的其余部分也没动
        //    (换 `to_be_bytes`、调换两个值的次序,这一条同样会红)
        let expect = |s: &Style| -> u64 {
            let mut h: u64 = 14_695_981_039_346_656_037;
            for v in [s.label_size(), s.chain_angle_deg] {
                for b in v.to_bits().to_le_bytes() {
                    h ^= u64::from(b);
                    h = h.wrapping_mul(1_099_511_628_211);
                }
            }
            h
        };
        let mut checked = 0usize;
        for s in &Style::ALL {
            assert_eq!(
                s.layout_fingerprint(),
                expect(s),
                "{} 的指纹不是 FNV-1a —— 常数写错了",
                s.name
            );
            checked += 1;
        }
        assert!(checked > 0, "一套规范都没查,这条判据在空过");
    }

    /// 一次逐字段扰动:字段名 + 怎么改它。
    type Tweak = (&'static str, fn(&mut Style));

    /// 指纹覆盖的**正好**是会挪动布局的那些字段 —— 一个不多、一个不少。
    ///
    /// 逐字段扰动一次,比两件事:**指纹变没变**、**坐标变没变**。两者必须一致。
    /// 少了(坐标变而指纹不变)就是静默错配拦不住;多了(指纹变而坐标不变)
    /// 就是换个线宽还要重排一遍。
    ///
    /// 这一条是补一个真漏的:先前指纹打的是 `clash_distance()`,布局一处都没
    /// 用过它 —— 两套 `label_size` 差一倍的规范能撞出同一个指纹,实测那样一对
    /// 规范下 8831 个分子有 787 个坐标不同而 `matches()` 一个都没查出来。
    ///
    /// 分子取自真实语料(第 4134、498 行),而且是**实测挑的** —— 合起来要能
    /// 覆盖标签尺寸与链角两个旋钮,否则"坐标没变"会是假象。
    #[test]
    fn the_fingerprint_covers_exactly_what_moves_the_layout() {
        // 两个都取自 `harness/corpus/large.smi`(第 4134、498 行)。
        // **不许自己造分子** —— 造出来的会把自身缺陷带进结论。
        let sensitive = ["CCC(CC)C(O)(C(CC)CC)C(O)=O", "BrC(Br)CS(=O)(=O)C1=CC=CC=C1"];
        let base = Style::ACS_1996;
        let mols: Vec<_> = sensitive.iter().map(|s| crate::tests_prep(s)).collect();
        let coords: Vec<Vec<crate::geom::Point2>> = mols
            .iter()
            .map(|m| crate::generate(m, &base).coords)
            .collect();

        // `bond_length_pt` 与 `atom_label_pt` 只经它们的**比值**进布局,所以这里
        // 让两者给出同一个比值变化 —— 它们理应同变同不变。
        let cases: [Tweak; 11] = [
            ("bond_length_pt", |s| s.bond_length_pt /= 1.4),
            ("atom_label_pt", |s| s.atom_label_pt *= 1.4),
            ("chain_angle_deg", |s| s.chain_angle_deg = 108.0),
            ("line_width_pt", |s| s.line_width_pt *= 1.4),
            ("bold_width_pt", |s| s.bold_width_pt *= 1.4),
            ("margin_width_pt", |s| s.margin_width_pt *= 1.4),
            ("hash_spacing_pt", |s| s.hash_spacing_pt *= 1.4),
            ("bond_spacing_pct", |s| s.bond_spacing_pct *= 1.4),
            ("caption_pt", |s| s.caption_pt *= 1.4),
            ("font_family", |s| s.font_family = "serif"),
            ("name", |s| s.name = "另一套规范"),
        ];
        // **`Style` 多一个字段,这里就编不过**(`error[E0027]` 会直接点名那个
        // 字段)—— 逼新字段必须在上面那张表里露面。少了这道闸门,一个真影响
        // 布局、指纹却漏掉的新字段会让这条判据一声不吭:实测加过一个
        // `extra_angle_deg` 并让 `chains::ideal_angle` 真的读它,五条判据全绿。
        let Style {
            name: _,
            bond_length_pt: _,
            line_width_pt: _,
            bold_width_pt: _,
            margin_width_pt: _,
            hash_spacing_pt: _,
            bond_spacing_pct: _,
            chain_angle_deg: _,
            atom_label_pt: _,
            caption_pt: _,
            font_family: _,
        } = &base;

        let (mut moved_any, mut still_any) = (false, false);
        for (field, tweak) in cases {
            let mut st = base.clone();
            tweak(&mut st);
            assert_ne!(st, base, "{field}:扰动没生效,这一档在空过");
            let fp_changed = st.layout_fingerprint() != base.layout_fingerprint();
            let moved = mols.iter().zip(&coords).any(|(m, c)| {
                crate::generate(m, &st)
                    .coords
                    .iter()
                    .zip(c)
                    .any(|(u, v)| u.dist(*v) > 1e-9)
            });
            assert_eq!(
                fp_changed, moved,
                "{field}:指纹变={fp_changed} 而坐标变={moved} —— 指纹的覆盖面对不上"
            );
            moved_any |= moved;
            still_any |= !moved;
        }
        // 防空过:两边都要真的出现过,否则上面那句可能一直在比 false == false
        assert!(
            moved_any,
            "没有一个字段挪动过布局 —— 分子挑得不敏感,这条判据在空过"
        );
        assert!(still_any, "没有一个字段是只进渲染的 —— 那也不对");
    }

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
        // 而且要落到实处:同一个分子,两套规范画出来的坐标**确实不同**。
        // 只比参数是不够的 —— 参数差着而布局不看,同样是白传。
        let m = crate::tests_prep("BrC(Br)CS(=O)(=O)C1=CC=CC=C1");
        let (a, b) = (crate::generate(&m, &acs), crate::generate(&m, &cd));
        assert!(
            a.coords
                .iter()
                .zip(&b.coords)
                .any(|(u, v)| u.dist(*v) > 1e-9),
            "两套规范画出来的坐标一模一样 —— 那把 Style 传进布局就是白传"
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
