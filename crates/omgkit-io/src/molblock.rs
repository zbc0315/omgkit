//! MDL molblock(V2000)的写出。
//!
//! # 只做 V2000
//!
//! V2000 的计数行给原子数和键数各留 3 位,所以**上限 999**。仓库里几份语料
//! 最大的分子是 122 个重原子(补完氢约 250),离上限很远;超过时这里
//! **报错而不是把数字挤出格** —— 挤出格的计数行别人读进去是另一个分子。
//! 真需要更大的分子时补 V3000,那是另一块活。
//!
//! # 立体信息有两条路,别混
//!
//! | 坐标 | 立体从哪儿来 |
//! |---|---|
//! | 2D | 每根键的**楔形码**(1 实楔、6 虚楔) |
//! | 3D | 坐标本身 |
//!
//! 三维构象**不写楔形码**:楔形是二维图上的记号,写在三维坐标旁边会让读的一方
//! 拿到两个可能互相矛盾的说法。
//!
//! # 窄端必须是键的第一个原子
//!
//! molblock 的楔形是有方向的:`1` 表示"从第一个原子指向第二个原子、朝观察者"。
//! 写反了,楔形描述的就是另一头那个原子的构型。这条规则**只写在这里** ——
//! 调用方给的是"窄端在哪个原子",端点谁先谁后由这个模块决定。

use core::fmt::Write as _;

use omgkit_core::{element, AtomFlags, BondOrder, MolBuilder};

/// V2000 计数行给原子数与键数各留 3 位。
const V2000_LIMIT: usize = 999;

/// 写不出来的原因。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WriteError {
    /// 原子数或键数超过 V2000 的 999 上限。
    TooLarge {
        /// 原子数
        atoms: usize,
        /// 键数
        bonds: usize,
    },
    /// 调用方给的逐键数组长度与键数对不上。
    BondArrayLen {
        /// 给了多少
        got: usize,
        /// 该给多少
        want: usize,
    },
    /// 有芳香键,而调用方没给凯库勒化之后的键级。
    ///
    /// molblock 的 4 号键级各家读法不一,写出去等于把歧义交给对方;而按单键写
    /// 是**更坏**的做法 —— 噻吩会被读成四氢噻吩,而且一声不响。所以这里报错,
    /// 由调用方先凯库勒化(`omgkit_chem::kekulize`)再把键级传进来。
    AromaticBond {
        /// 第一根芳香键的下标
        bond: usize,
    },
}

impl core::fmt::Display for WriteError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::TooLarge { atoms, bonds } => write!(
                f,
                "{atoms} 个原子 / {bonds} 根键,超过 V2000 的 {V2000_LIMIT} 上限"
            ),
            Self::BondArrayLen { got, want } => {
                write!(f, "逐键数组给了 {got} 项,而分子有 {want} 根键")
            }
            Self::AromaticBond { bond } => write!(
                f,
                "第 {bond} 根键是芳香键,而调用方没给凯库勒化之后的键级 —— \
                 先跑 `omgkit_chem::kekulize`"
            ),
        }
    }
}

/// 一根键在图上画成什么样。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BondWedge {
    /// 普通实线
    #[default]
    Plain,
    /// 实楔形:窄端那个原子在纸面,另一端朝观察者
    Up {
        /// 窄端所在的原子
        narrow: u32,
    },
    /// 虚楔形:窄端那个原子在纸面,另一端背离观察者
    Down {
        /// 窄端所在的原子
        narrow: u32,
    },
}

/// 写一条 molblock 要的东西。
#[derive(Debug, Clone, Copy)]
pub struct Record<'a> {
    /// 第一行的标题。空串就是空行 —— 合法。
    pub title: &'a str,
    /// 逐原子坐标。二维图把 `z` 全给 0。
    pub coords: &'a [[f64; 3]],
    /// 逐键的楔形。空切片表示全是普通实线(三维构象就该是空的)。
    pub wedges: &'a [BondWedge],
    /// 逐键的键级。空切片表示照分子里存的写。
    ///
    /// **芳香键要先凯库勒化**(`omgkit_chem::kekulize`)。留着芳香键直接写会
    /// 报 [`WriteError::AromaticBond`] —— 不是退回单键:那样噻吩写出去、
    /// 读回来就成了四氢噻吩,而且一声不响。
    pub orders: &'a [BondOrder],
}

impl<'a> Record<'a> {
    /// 三维构象的最简用法:只给坐标,不画楔形,键级照分子里存的。
    #[must_use]
    pub fn from_coords(title: &'a str, coords: &'a [[f64; 3]]) -> Self {
        Self {
            title,
            coords,
            wedges: &[],
            orders: &[],
        }
    }
}

/// 写成一条 V2000 molblock,末尾带 `M  END`(不带 SDF 的 `$$$$`)。
///
/// 多条拼成 SDF 时,每条后面接数据字段(`> <名字>` + 值 + 空行)再接 `$$$$`。
///
/// # Errors
///
/// 原子数或键数超过 999,或逐键数组长度与键数对不上。
pub fn write_v2000(mol: &MolBuilder, rec: &Record) -> Result<String, WriteError> {
    let (na, nb) = (mol.num_atoms(), mol.num_bonds());
    if na > V2000_LIMIT || nb > V2000_LIMIT {
        return Err(WriteError::TooLarge {
            atoms: na,
            bonds: nb,
        });
    }
    for (got, want) in [(rec.wedges.len(), nb), (rec.orders.len(), nb)] {
        if got != 0 && got != want {
            return Err(WriteError::BondArrayLen { got, want });
        }
    }
    // **芳香键必须先凯库勒化。** 不报错的话它会落进下面的"其余按单键",
    // 而那是一个安静的错:噻吩写出去、读回来就是四氢噻吩。
    if rec.orders.is_empty() {
        if let Some(bond) = mol
            .bonds()
            .iter()
            .position(|b| b.order == BondOrder::Aromatic)
        {
            return Err(WriteError::AromaticBond { bond });
        }
    } else if let Some(bond) = rec.orders.iter().position(|o| *o == BondOrder::Aromatic) {
        return Err(WriteError::AromaticBond { bond });
    }

    let mut out = String::with_capacity(64 + na * 70 + nb * 22);
    // 三行头:标题 / 程序行 / 注释行。第二行按规范是程序名与时间戳,
    // **不写时间戳** —— 同一个分子每次写出都该逐字节相同,时间戳会毁掉这一条。
    let _ = writeln!(out, "{}", rec.title);
    let _ = writeln!(out, "  omgkit");
    out.push('\n');
    let _ = writeln!(out, "{na:>3}{nb:>3}  0  0  0  0  0  0  0  0999 V2000");

    for (i, a) in mol.atoms().iter().enumerate() {
        let p = rec.coords.get(i).copied().unwrap_or([0.0; 3]);
        let sym = element::by_atomic_num(a.atomic_num).map_or("*", |e| e.symbol);
        // **价键字段 `vvv`:氢数就靠它钉住。**
        //
        // molblock 的原子块不写氢数,读的一方按元素的默认价自己补 ——
        // 于是 `[CH]` 读回来成了 `[CH3]`,自由基碳读回来成了甲基。写上这个字段
        // 就把总价钉死,补氢无从下手。参照实现也是这么做的(实测 `[CH2]C` 的
        // 自由基碳那一行写的正是 `vvv=3`)。
        //
        // 只在**必要时**写:作者钉过氢数(`NO_IMPLICIT`)或者带自由基电子。
        // 常规原子留 0,让读的一方按默认价补 —— 那正是它该做的。
        let valence = if a.flags.contains(AtomFlags::NO_IMPLICIT) || a.num_radical_electrons != 0 {
            omgkit_core::valence::explicit_valence_nonstrict(mol, u32::try_from(i).unwrap_or(0))
                .clamp(0, 14)
        } else {
            0
        };
        let _ = writeln!(
            out,
            "{:>10.4}{:>10.4}{:>10.4} {sym:<3} 0  0  0  0  0{valence:>3}  0  0  0  0  0  0",
            p[0], p[1], p[2]
        );
    }

    for (bi, b) in mol.bonds().iter().enumerate() {
        let wedge = rec.wedges.get(bi).copied().unwrap_or_default();
        // **窄端必须写成键的第一个原子。** 这条规则只在这里。
        let (first, second, code) = match wedge {
            BondWedge::Plain => (b.begin, b.end, 0),
            BondWedge::Up { narrow } | BondWedge::Down { narrow } => {
                let other = if narrow == b.begin { b.end } else { b.begin };
                let code = if matches!(wedge, BondWedge::Up { .. }) {
                    1
                } else {
                    6
                };
                (narrow, other, code)
            }
        };
        let order = rec.orders.get(bi).copied().unwrap_or(b.order);
        let ord = match order {
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            _ => 1,
        };
        let _ = writeln!(
            out,
            "{:>3}{:>3}{ord:>3}{code:>3}  0  0  0",
            first + 1,
            second + 1
        );
    }

    // 属性块。**电荷、同位素、自由基都要写** —— 少写一样,读的一方按默认补氢,
    // 拿到的是另一个分子。原子块里那两个旧字段(质量差、旧电荷码)一律留 0,
    // 现代读法只认这里的 `M` 行。
    for (i, a) in mol.atoms().iter().enumerate() {
        if a.formal_charge != 0 {
            let _ = writeln!(out, "M  CHG  1{:>4}{:>4}", i + 1, a.formal_charge);
        }
        if a.isotope != 0 {
            let _ = writeln!(out, "M  ISO  1{:>4}{:>4}", i + 1, a.isotope);
        }
        if a.num_radical_electrons != 0 {
            // 规范的编码:1 = 双线态(1 个电子)、2 = 单线态、3 = 三线态(2 个)。
            // 我们只记电子数,按 1 个 → 2(双线态)、2 个 → 3(三线态)写,
            // 与 RDKit 的写法一致。
            let code = match a.num_radical_electrons {
                1 => 2,
                n => i32::from(n) + 1,
            };
            let _ = writeln!(out, "M  RAD  1{:>4}{code:>4}", i + 1);
        }
    }
    out.push_str("M  END\n");
    Ok(out)
}
