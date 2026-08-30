//! **CPK 元素配色** —— 三维分子图按元素上色的那张表。
//!
//! # 为什么不放在元素表里
//!
//! [`omgkit_core::element::Element`] 装的是**物性**:共价半径、范德华半径、
//! 默认价、Pauling 电负性。颜色不是物性 —— 它是一条**绘图惯例**,由 Corey、
//! Pauling、Koltun 三个人在实体模型上定下来,后来各家软件各自微调。同一个碳
//! 在 Jmol 是 `#909090`、在 Rasmol 是 `#c8c8c8`,而它的范德华半径两边都是
//! 1.7 Å。把两者混在一张表里,下一个人就会以为颜色也有个"正确值"。
//!
//! # 取 Jmol 那一套
//!
//! 三套通行配色(原始 CPK、Rasmol、Jmol)里取 Jmol:它是覆盖最全的一套
//! (到 Mt,109 个元素;Rasmol 那套一大半元素落进"未知"的粉色),而且
//! PyMOL、VMD、3Dmol.js、ASE 的默认元素色都与它同源。
//!
//! 表的出处、以及"两份互相独立的转录逐字节相同"这件事,写在
//! `harness/params/jmol_colors.tsv` 的头注释里。

use crate::palette_data::CPK;

/// 元素 `atomic_num` 的 CPK 颜色,RGB 各一字节。
///
/// # 表外的元素给的是什么
///
/// Jmol 的表到 Mt(109)为止。110 以上、以及 SMILES 的通配原子 `*`(这里
/// 记作原子序数 0),拿的是 Jmol 的**未知元素色** deeppink `#ff1493` ——
/// 那是"表里没有"的标志色,不是"这个元素是粉的"。它刺眼是故意的。
///
/// 原子序数超出周期表(> 118)时同样给这个颜色。
///
/// ```
/// use omgkit_depict::palette::cpk;
/// assert_eq!(cpk(6), [0x90, 0x90, 0x90]);   // 碳,Jmol 的灰
/// assert_eq!(cpk(8), [0xff, 0x0d, 0x0d]);   // 氧,红
/// assert_eq!(cpk(0), [0xff, 0x14, 0x93]);   // 通配原子:未知色
/// ```
#[must_use]
pub fn cpk(atomic_num: u8) -> [u8; 3] {
    CPK.get(atomic_num as usize).copied().unwrap_or(CPK[0])
}

#[cfg(test)]
mod tests {
    use super::{cpk, CPK};

    /// 表的长度与元素表一致 —— 否则按原子序数索引会在末尾静默落进 `unwrap_or`。
    #[test]
    fn 表长与元素表相同() {
        assert_eq!(CPK.len(), omgkit_core::element::count());
    }

    /// **配色表的元素符号与元素表逐项对上。**
    ///
    /// 两张表来自不同的上游:元素表转自 BODR(经 RDKit),配色表转自 Jmol。
    /// 生成器只按原子序数往下填,谁都没检查过"第 16 行真的是硫"。这条判据
    /// 逐项对号:表里的符号注释若与元素表错位,颜色就整体挪了位 —— 而挪位之后
    /// 图还是画得出来,只是每个元素都穿了别人的衣服。
    ///
    /// 注释在生成的文件里,读不到运行时,所以这里改比**颜色本身**:
    /// 拿几个各家文档都写明的元素当锚点。
    #[test]
    fn 锚点元素的颜色与各家文档一致() {
        // 这几个是 CPK 惯例里最没有争议的,写错一个整张表多半也错位了
        assert_eq!(cpk(1), [0xff, 0xff, 0xff], "H 白");
        assert_eq!(cpk(6), [0x90, 0x90, 0x90], "C 灰");
        assert_eq!(cpk(7), [0x30, 0x50, 0xf8], "N 蓝");
        assert_eq!(cpk(8), [0xff, 0x0d, 0x0d], "O 红");
        assert_eq!(cpk(15), [0xff, 0x80, 0x00], "P 橙");
        assert_eq!(cpk(16), [0xff, 0xff, 0x30], "S 黄");
        assert_eq!(cpk(17), [0x1f, 0xf0, 0x1f], "Cl 绿");
        assert_eq!(cpk(35), [0xa6, 0x29, 0x29], "Br 暗红");
    }

    /// 表里 109 个元素之后必须**全是**未知色,而 109 及以前一个都不许是。
    ///
    /// 反着也要判:只查"110 以后是粉的"的话,生成器把整张表填成粉色也照样绿。
    #[test]
    fn 未知色只出现在表外() {
        let unknown = [0xff, 0x14, 0x93];
        for z in 1..=109u8 {
            assert_ne!(cpk(z), unknown, "Z={z} 在 Jmol 表里,不该是未知色");
        }
        for z in 110..=118u8 {
            assert_eq!(cpk(z), unknown, "Z={z} 不在 Jmol 表里");
        }
        assert_eq!(cpk(0), unknown, "通配原子不在表里");
    }

    /// 越界不 panic,给未知色。
    #[test]
    fn 越界的原子序数给未知色() {
        assert_eq!(cpk(200), [0xff, 0x14, 0x93]);
        assert_eq!(cpk(u8::MAX), [0xff, 0x14, 0x93]);
    }
}
