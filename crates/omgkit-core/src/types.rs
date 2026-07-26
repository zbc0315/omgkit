//! 基础标量类型:键级、手性、立体、以及原子/键的标志位。
//!
//! 这些类型刻意保持为 `Copy` 的小整数包装,以便它们能直接作为
//! [`MolBatch`](crate::MolBatch) 的列元素存储 —— 列在内存中必须是
//! 紧凑的 POD 数组,才能同时满足 SIMD、`memcpy` 到 GPU 和零拷贝
//! 暴露给 numpy 三个要求。

/// 键级。
///
/// `Aromatic` 只在净化(L2)之后出现。解析阶段(L1)由小写原子推断出的
/// 芳香键先记为 `Aromatic`,但此时它只是"作者声称芳香",尚未经过感知验证。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default, PartialOrd, Ord)]
#[repr(u8)]
pub enum BondOrder {
    /// 未指定(查询键、或尚未 kekulize 的占位)
    #[default]
    Unspecified = 0,
    /// 单键
    Single = 1,
    /// 双键
    Double = 2,
    /// 三键
    Triple = 3,
    /// 四重键(金属-金属)
    Quadruple = 4,
    /// 芳香键
    Aromatic = 5,
    /// 配位键。两个电子都由 [`begin`](crate::BondData::begin) 端提供,
    /// 所以端点顺序**有语义**,不能随意对调。
    Dative = 6,
}

impl BondOrder {
    /// 键级的数值形式。
    ///
    /// 配位键在这里算 **1** —— 它是一条 σ 键,判断"是否不饱和"这类问题时
    /// 与单键同级。
    ///
    /// 注意这**不是**价键计算用的量 —— 配位键对两端的价贡献不对称,
    /// 见 [`BondData::valence_contribution_to`](crate::BondData::valence_contribution_to)。
    #[must_use]
    pub fn as_double(self) -> f32 {
        match self {
            Self::Unspecified => 0.0,
            Self::Single | Self::Dative => 1.0,
            Self::Aromatic => 1.5,
            Self::Double => 2.0,
            Self::Triple => 3.0,
            Self::Quadruple => 4.0,
        }
    }

    /// SMILES 中的键符号。单键与芳香键按惯例省略,故返回 `None`。
    #[must_use]
    pub fn smiles_symbol(self) -> Option<char> {
        match self {
            Self::Single | Self::Aromatic | Self::Unspecified => None,
            Self::Double => Some('='),
            Self::Triple => Some('#'),
            Self::Quadruple => Some('$'),
            Self::Dative => Some('>'),
        }
    }
}

/// 原子上的立体标记 —— 记录的是**几何类别**,不是 CIP 的 R/S。
///
/// R/S 需要 CIP 优先级排序,属于 L6;这里只承载书写时给出的排列信息。
///
/// # 两部分信息
///
/// 一个立体标记由"几何类别 + 类内排列序号"两部分构成。类别存在本枚举里,
/// 序号存在 [`AtomData::stereo_perm`](crate::AtomData::stereo_perm)。
///
/// | 类别 | 配位数 | 排列数 | 序号存放 |
/// |---|---|---|---|
/// | [`Cw`](Self::Cw) / [`Ccw`](Self::Ccw) | 4 | 2 | 已由变体本身表达 |
/// | [`SquarePlanar`](Self::SquarePlanar) | 4 | 3 | `stereo_perm` |
/// | [`TrigonalBipyramidal`](Self::TrigonalBipyramidal) | 5 | 20 | `stereo_perm` |
/// | [`Octahedral`](Self::Octahedral) | 6 | 30 | `stereo_perm` |
/// | [`Other`](Self::Other) | — | — | `stereo_perm` |
///
/// 四面体的两种排列直接做成了两个变体,因为"邻居顺序对换一次就翻转"这条
/// 规则只对它成立(见 [`inverted`](Self::inverted)),下游到处都在用。
/// 其余类别的排列在邻居重排下的变换是一张查找表,不能用一个布尔位表达。
///
/// # 四面体两个变体的含义
///
/// 相对于邻居在**存储中**的顺序:从第一个邻居看向该原子,其余邻居依次
/// 是逆时针([`Ccw`](Self::Ccw),写作 `@`)还是顺时针([`Cw`](Self::Cw),
/// 写作 `@@`)。存储顺序与 SMILES 书写顺序不一定一致,L1 解析时已做过补偿。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[repr(u8)]
pub enum ChiralTag {
    /// 无立体标记
    #[default]
    Unspecified = 0,
    /// 四面体,`@@` —— 顺时针
    Cw = 1,
    /// 四面体,`@` —— 逆时针
    Ccw = 2,
    /// 有立体标记但不属于下列任何一种几何。`@AL`(丙二烯轴手性)归入此类:
    /// 它的立体信息属于**一根轴**而非一个中心,与配位几何不是一回事。
    Other = 3,
    /// 平面四方形,`@SP`
    SquarePlanar = 4,
    /// 三角双锥,`@TB`
    TrigonalBipyramidal = 5,
    /// 八面体,`@OH`
    Octahedral = 6,
}

impl ChiralTag {
    /// 反转手性。邻居顺序发生一次对换时使用。
    ///
    /// 只有四面体会变 —— 其余类别的排列序号在邻居重排下按查找表变换,
    /// 不是一个可以就地翻转的布尔量,故原样返回。这些类别的重排要连同
    /// [`AtomData::stereo_perm`](crate::AtomData::stereo_perm) 一起处理。
    #[must_use]
    pub fn inverted(self) -> Self {
        match self {
            Self::Cw => Self::Ccw,
            Self::Ccw => Self::Cw,
            other => other,
        }
    }

    /// 是否是四面体手性(即 [`inverted`](Self::inverted) 对其有效)。
    #[must_use]
    pub fn is_tetrahedral(self) -> bool {
        matches!(self, Self::Cw | Self::Ccw)
    }
}

/// SMILES 中的方向键(`/` 与 `\`),用于表达双键顺反。
///
/// 注意这是**书写方向**而非最终的 Z/E 判定;Z/E 要在净化后由
/// [`BondStereo`] 承载。二者的区别是 SMILES 立体处理里最常见的混淆点。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[repr(u8)]
pub enum BondDirection {
    /// 无方向
    #[default]
    None = 0,
    /// `/`
    UpRight = 1,
    /// `\`
    DownRight = 2,
}

impl BondDirection {
    /// 翻转方向(书写顺序颠倒时使用)。
    #[must_use]
    pub fn flipped(self) -> Self {
        match self {
            Self::UpRight => Self::DownRight,
            Self::DownRight => Self::UpRight,
            Self::None => Self::None,
        }
    }
}

/// 双键立体化学:两个取代基在双键同侧还是异侧。
///
/// 由 [`stereo::perceive_bond_stereo`](../../omgkit_io/stereo/fn.perceive_bond_stereo.html)
/// 从 [`BondDirection`] 感知得出,参照原子记在
/// [`BondData::stereo_atoms`](crate::BondData::stereo_atoms) 里。
///
/// # 与 [`BondDirection`] 的分工
///
/// 方向是**写法**,依附于某根单键 —— 那根键被图编辑删掉,信息就没了。
/// 顺反是双键**自己**的性质,只要两个参照原子还在就一直成立。
/// 反应改图之后仍能写出正确的方向键,靠的正是这一层。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[repr(u8)]
pub enum BondStereo {
    /// 未指定
    #[default]
    None = 0,
    /// Z(CIP 优先级较高的取代基同侧)
    Z = 1,
    /// E
    E = 2,
    /// 顺(相对于记录的参照原子,不涉及 CIP)
    Cis = 3,
    /// 反
    Trans = 4,
}

macro_rules! flag_set {
    ($(#[$meta:meta])* $name:ident { $($(#[$fmeta:meta])* $flag:ident = $bit:expr),+ $(,)? }) => {
        $(#[$meta])*
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
        pub struct $name(u32);

        impl $name {
            /// 空标志集
            pub const NONE: Self = Self(0);
            $($(#[$fmeta])* pub const $flag: Self = Self(1 << $bit);)+

            /// 是否包含 `other` 的全部位
            #[must_use]
            pub const fn contains(self, other: Self) -> bool {
                (self.0 & other.0) == other.0
            }

            /// 置位
            pub fn insert(&mut self, other: Self) {
                self.0 |= other.0;
            }

            /// 清位
            pub fn remove(&mut self, other: Self) {
                self.0 &= !other.0;
            }

            /// 按条件置位或清位
            pub fn set(&mut self, other: Self, on: bool) {
                if on { self.insert(other) } else { self.remove(other) }
            }

            /// 原始位
            #[must_use]
            pub const fn bits(self) -> u32 {
                self.0
            }
        }

        impl core::ops::BitOr for $name {
            type Output = Self;
            fn bitor(self, rhs: Self) -> Self {
                Self(self.0 | rhs.0)
            }
        }

        impl core::ops::BitOrAssign for $name {
            fn bitor_assign(&mut self, rhs: Self) {
                self.0 |= rhs.0;
            }
        }
    };
}

/// 原子的杂化状态。
///
/// 由"成键数 + 孤对数"推出,不是从几何构型测量的。判别值与列式存储的编码
/// 一一对应,改动需同步更新差分基准的编码表。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[repr(u8)]
pub enum Hybridization {
    /// 未确定
    #[default]
    Unspecified = 0,
    /// s —— 孤立原子或只有一个轨道
    S = 1,
    /// sp —— 两个轨道,直线形
    Sp = 2,
    /// sp² —— 三个轨道,平面三角形
    Sp2 = 3,
    /// sp³ —— 四个轨道,四面体
    Sp3 = 4,
    /// sp²d —— 平面四方形
    Sp2d = 5,
    /// sp³d —— 三角双锥
    Sp3d = 6,
    /// sp³d² —— 八面体
    Sp3d2 = 7,
}

flag_set! {
    /// 原子标志位。
    AtomFlags {
        /// 芳香。L1 记录作者声称,L2 净化后为感知结果。
        AROMATIC = 0,
        /// 位于任意环中(L2 环感知后置位)
        IN_RING = 1,
        /// 原子写在方括号中 —— 隐式氢由作者显式指定,不再推断。
        /// 这是 `[CH4]` 与 `C` 语义不同的根源。
        NO_IMPLICIT = 2,
        /// 该原子是查询原子(L4),语义由外挂的查询树承载
        HAS_QUERY = 3,
        /// 共轭体系成员(L2)
        CONJUGATED = 4,
    }
}

flag_set! {
    /// 键标志位。
    BondFlags {
        /// 芳香
        AROMATIC = 0,
        /// 位于任意环中
        IN_RING = 1,
        /// 共轭
        CONJUGATED = 2,
        /// 该键是查询键(L4)
        HAS_QUERY = 3,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bond_order_as_double() {
        assert_eq!(BondOrder::Single.as_double(), 1.0);
        assert_eq!(BondOrder::Aromatic.as_double(), 1.5);
        // 配位键的数值键级是 1;
        // 它对两端价贡献的不对称性由 BondData::valence_contribution_to 处理
        assert_eq!(BondOrder::Dative.as_double(), 1.0);
        assert_eq!(BondOrder::Unspecified.as_double(), 0.0);
    }

    /// 反转必须是对合的,**每一个**变体都要满足 —— 非四面体的变体走的是
    /// "原样返回"那条路,同样是这个契约的一部分。
    #[test]
    fn chiral_tag_inversion_is_involutive() {
        for tag in [
            ChiralTag::Unspecified,
            ChiralTag::Cw,
            ChiralTag::Ccw,
            ChiralTag::Other,
            ChiralTag::SquarePlanar,
            ChiralTag::TrigonalBipyramidal,
            ChiralTag::Octahedral,
        ] {
            assert_eq!(tag.inverted().inverted(), tag, "{tag:?}");
            // 只有四面体会真的变;其余原样返回
            if tag.is_tetrahedral() {
                assert_ne!(tag.inverted(), tag, "{tag:?} 应当翻转");
            } else {
                assert_eq!(tag.inverted(), tag, "{tag:?} 不该翻转");
            }
        }
    }

    #[test]
    fn bond_direction_flip_is_involutive() {
        for d in [
            BondDirection::UpRight,
            BondDirection::DownRight,
            BondDirection::None,
        ] {
            assert_eq!(d.flipped().flipped(), d);
        }
    }

    #[test]
    fn flags_basic_ops() {
        let mut f = AtomFlags::NONE;
        assert!(!f.contains(AtomFlags::AROMATIC));

        f.insert(AtomFlags::AROMATIC);
        f.insert(AtomFlags::IN_RING);
        assert!(f.contains(AtomFlags::AROMATIC));
        assert!(f.contains(AtomFlags::AROMATIC | AtomFlags::IN_RING));

        f.remove(AtomFlags::AROMATIC);
        assert!(!f.contains(AtomFlags::AROMATIC));
        assert!(f.contains(AtomFlags::IN_RING));

        f.set(AtomFlags::NO_IMPLICIT, true);
        assert!(f.contains(AtomFlags::NO_IMPLICIT));
        f.set(AtomFlags::NO_IMPLICIT, false);
        assert!(!f.contains(AtomFlags::NO_IMPLICIT));
    }
}
