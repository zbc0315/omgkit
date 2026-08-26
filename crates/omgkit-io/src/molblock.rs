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

// ---------------------------------------------------------------------------
// 读取
// ---------------------------------------------------------------------------

/// 读不了的原因。
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReadError {
    /// 行数不够(头三行 + 计数行 + 原子块 + 键块)。
    Truncated {
        /// 缺的是哪一段
        what: &'static str,
    },
    /// 计数行读不出原子数/键数。
    BadCounts,
    /// V3000。本模块只做 V2000 —— 认出来并明确报错,好过把它当 V2000 硬读。
    V3000,
    /// 某一行的格式不对。
    BadLine {
        /// 行号(从 0 数)
        line: usize,
        /// 哪儿不对
        what: &'static str,
    },
    /// 不认识的元素符号。
    UnknownElement {
        /// 行号
        line: usize,
        /// 读到的符号
        symbol: String,
    },
    /// 键的类型不在 1..=4。查询用的 5..=8 这里不收 —— 它们不是分子。
    BadBondType {
        /// 行号
        line: usize,
        /// 读到的值
        got: i32,
    },
}

impl core::fmt::Display for ReadError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Truncated { what } => write!(f, "molblock 截断了:缺{what}"),
            Self::BadCounts => write!(f, "计数行读不出原子数与键数"),
            Self::V3000 => write!(f, "这是 V3000,本模块只读 V2000"),
            Self::BadLine { line, what } => write!(f, "第 {line} 行{what}"),
            Self::UnknownElement { line, symbol } => {
                write!(f, "第 {line} 行:不认识的元素符号 `{symbol}`")
            }
            Self::BadBondType { line, got } => {
                write!(f, "第 {line} 行:键类型 {got} 不在 1..=4(5..=8 是查询用的)")
            }
        }
    }
}

/// 读出来的东西。
///
/// **立体化学不在 `mol` 里。** 二维图的立体靠 [`wedges`](Self::wedges),三维的
/// 靠 [`coords`](Self::coords) —— 两者都要在更上一层赋值(赋值要用对称等价类,
/// 那在 L1 之上)。这里把两样都如实交出来,而不是悄悄给一个没有立体的分子。
#[derive(Debug, Clone)]
pub struct Molblock {
    /// 第一行的标题。
    pub title: String,
    /// 分子。**未净化**:价键没算、环没感知、芳香性没感知。
    pub mol: MolBuilder,
    /// 逐原子坐标。二维图的 `z` 是 0。
    pub coords: Vec<[f64; 3]>,
    /// 逐键的楔形。
    pub wedges: Vec<BondWedge>,
    /// 坐标是不是三维的(有任何一个 `z` 不为 0)。
    ///
    /// 二维和三维的立体读法完全不同,而文件里没有哪个字段直说 —— 只能这么判,
    /// 与 RDKit 同法。
    pub is_3d: bool,
}

/// 取固定列的一段并去掉空白。列超出行尾时给空串 —— molblock 允许行尾截断。
fn field(line: &str, from: usize, to: usize) -> &str {
    let b = line.as_bytes();
    if from >= b.len() {
        return "";
    }
    let end = to.min(b.len());
    // 只在 ASCII 边界上切;molblock 是 ASCII 格式,非 ASCII 只可能出现在标题里
    line.get(from..end).unwrap_or("").trim()
}

fn parse_i32(s: &str) -> Option<i32> {
    if s.is_empty() {
        Some(0)
    } else {
        s.parse().ok()
    }
}

/// 读一条 V2000 molblock。到 `M  END` 为止,后面的东西(SDF 的数据字段、
/// `$$$$`)一概不看。
///
/// # Errors
///
/// 截断、计数行读不出、V3000、某行格式不对、元素不认识、键类型不在 1..=4。
#[allow(clippy::too_many_lines)]
pub fn read_v2000(text: &str) -> Result<Molblock, ReadError> {
    // **`\r` 要去掉。** 真实文件多半来自 Windows,留着 `\r` 会让最后一个字段
    // 带上它 —— 元素符号 `C\r` 查不到,而错误信息里看不出多了什么。
    let lines: Vec<&str> = text.lines().map(|l| l.trim_end_matches('\r')).collect();
    if lines.len() < 4 {
        return Err(ReadError::Truncated { what: "头四行" });
    }
    let counts = lines[3];
    if counts.contains("V3000") {
        return Err(ReadError::V3000);
    }
    let na = parse_i32(field(counts, 0, 3)).ok_or(ReadError::BadCounts)?;
    let nb = parse_i32(field(counts, 3, 6)).ok_or(ReadError::BadCounts)?;
    let (na, nb) = (
        usize::try_from(na).map_err(|_| ReadError::BadCounts)?,
        usize::try_from(nb).map_err(|_| ReadError::BadCounts)?,
    );
    if lines.len() < 4 + na + nb {
        return Err(ReadError::Truncated {
            what: "原子块或键块",
        });
    }

    let mut mol = MolBuilder::with_capacity(na, nb);
    let mut coords = Vec::with_capacity(na);
    // 原子块里那个旧电荷码,先记着 —— **只有在没有任何 `M  CHG` 行时才作数**。
    let mut legacy_charge = vec![0i8; na];
    let mut legacy_iso = vec![0i16; na];

    for k in 0..na {
        let ln = 4 + k;
        let line = lines[ln];
        let bad = |what| ReadError::BadLine { line: ln, what };
        let x: f64 = field(line, 0, 10).parse().map_err(|_| bad("x 读不出来"))?;
        let y: f64 = field(line, 10, 20).parse().map_err(|_| bad("y 读不出来"))?;
        let z: f64 = field(line, 20, 30).parse().map_err(|_| bad("z 读不出来"))?;
        let sym = field(line, 31, 34);
        let el = element::by_symbol(sym).ok_or_else(|| ReadError::UnknownElement {
            line: ln,
            symbol: sym.to_string(),
        })?;
        let idx = mol.add_atom(el.atomic_num);

        // 旧电荷码:0 无、1 = +3、2 = +2、3 = +1、4 = 双线态自由基、5 = −1、
        // 6 = −2、7 = −3。这个编码常年被写错,所以只在没有 `M  CHG` 时才用。
        legacy_charge[k] = match parse_i32(field(line, 36, 39)).unwrap_or(0) {
            1 => 3,
            2 => 2,
            3 => 1,
            5 => -1,
            6 => -2,
            7 => -3,
            _ => 0,
        };
        // 质量差:相对该元素**最常见同位素**的差值。0 表示不指定。
        let mass_diff = parse_i32(field(line, 34, 36)).unwrap_or(0);
        if mass_diff != 0 {
            legacy_iso[k] = i16::try_from(i32::from(el.common_isotope) + mass_diff).unwrap_or(0);
        }
        // 价键字段:0 = 按默认价补氢;15 = 零价;1..=14 = 总价钉死。
        let valence = parse_i32(field(line, 48, 51)).unwrap_or(0);
        if valence != 0 {
            if let Some(a) = mol.atom_mut(idx) {
                a.flags.insert(AtomFlags::NO_IMPLICIT);
            }
        }
        coords.push([x, y, z]);
    }

    let mut wedges = Vec::with_capacity(nb);
    for k in 0..nb {
        let ln = 4 + na + k;
        let line = lines[ln];
        let bad = |what| ReadError::BadLine { line: ln, what };
        let a = parse_i32(field(line, 0, 3)).ok_or_else(|| bad("第一个原子号读不出来"))?;
        let b = parse_i32(field(line, 3, 6)).ok_or_else(|| bad("第二个原子号读不出来"))?;
        let t = parse_i32(field(line, 6, 9)).ok_or_else(|| bad("键类型读不出来"))?;
        let stereo = parse_i32(field(line, 9, 12)).unwrap_or(0);
        let (a, b) = (
            usize::try_from(a - 1).map_err(|_| bad("原子号越界"))?,
            usize::try_from(b - 1).map_err(|_| bad("原子号越界"))?,
        );
        if a >= na || b >= na {
            return Err(bad("原子号越界"));
        }
        let order = match t {
            1 => BondOrder::Single,
            2 => BondOrder::Double,
            3 => BondOrder::Triple,
            4 => BondOrder::Aromatic,
            got => return Err(ReadError::BadBondType { line: ln, got }),
        };
        let (a, b) = (
            u32::try_from(a).map_err(|_| bad("原子号越界"))?,
            u32::try_from(b).map_err(|_| bad("原子号越界"))?,
        );
        mol.add_bond(a, b, order)
            .map_err(|_| bad("这根键建不起来"))?;
        // 楔形的**窄端就是第一个原子** —— 与写出侧同一条规则
        wedges.push(match stereo {
            1 => BondWedge::Up { narrow: a },
            6 => BondWedge::Down { narrow: a },
            _ => BondWedge::Plain,
        });
    }

    // 属性块。**`M  CHG` / `M  ISO` 一出现,原子块里那两个旧字段整体作废** ——
    // 规范就是这么定的,而"两处都读、后者覆盖前者"会在只写了一部分原子的文件上
    // 给出错的电荷。
    let mut saw_chg = false;
    let mut saw_iso = false;
    let mut chg: Vec<(usize, i8)> = Vec::new();
    let mut iso: Vec<(usize, i16)> = Vec::new();
    let mut rad: Vec<(usize, u8)> = Vec::new();
    for line in &lines[4 + na + nb..] {
        if line.starts_with("M  END") {
            break;
        }
        let tag = field(line, 0, 6);
        if !matches!(tag, "M  CHG" | "M  ISO" | "M  RAD") {
            continue;
        }
        let n = parse_i32(field(line, 6, 9)).unwrap_or(0).max(0);
        for e in 0..usize::try_from(n).unwrap_or(0) {
            let at = parse_i32(field(line, 9 + e * 8, 13 + e * 8)).unwrap_or(0);
            let v = parse_i32(field(line, 13 + e * 8, 17 + e * 8)).unwrap_or(0);
            let Ok(at) = usize::try_from(at - 1) else {
                continue;
            };
            if at >= na {
                continue;
            }
            match tag {
                "M  CHG" => {
                    saw_chg = true;
                    chg.push((at, i8::try_from(v).unwrap_or(0)));
                }
                "M  ISO" => {
                    saw_iso = true;
                    iso.push((at, i16::try_from(v).unwrap_or(0)));
                }
                _ => rad.push((at, u8::try_from(radical_electrons(v)).unwrap_or(0))),
            }
        }
    }
    for k in 0..na {
        let Some(a) = mol.atom_mut(u32::try_from(k).unwrap_or(0)) else {
            continue;
        };
        if !saw_chg {
            a.formal_charge = legacy_charge[k];
        }
        if !saw_iso {
            a.isotope = u16::try_from(legacy_iso[k]).unwrap_or(0);
        }
    }
    for (at, v) in chg {
        if let Some(a) = mol.atom_mut(u32::try_from(at).unwrap_or(0)) {
            a.formal_charge = v;
        }
    }
    for (at, v) in iso {
        if let Some(a) = mol.atom_mut(u32::try_from(at).unwrap_or(0)) {
            a.isotope = u16::try_from(v).unwrap_or(0);
        }
    }
    for (at, v) in rad {
        if let Some(a) = mol.atom_mut(u32::try_from(at).unwrap_or(0)) {
            a.num_radical_electrons = v;
        }
    }

    let is_3d = coords.iter().any(|p| p[2].abs() > 1e-9);
    Ok(Molblock {
        title: (*lines.first().unwrap_or(&"")).to_string(),
        mol,
        coords,
        wedges,
        is_3d,
    })
}

/// `M  RAD` 的编码 → 未成对电子数。1 = 双线态(1 个)、2 = 单线态(2 个)、
/// 3 = 三线态(2 个)。
fn radical_electrons(code: i32) -> i32 {
    match code {
        1 | 2 => code.min(2),
        3 => 2,
        _ => 0,
    }
}
