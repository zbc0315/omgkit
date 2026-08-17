//! SMILES 解析 —— 手写递归下降。
//!
//! # 为什么是手写而不是生成的解析器
//!
//! SMILES 本质是**带环闭合的 DFS 序列**,LL(1) 足够。手写换来三样东西:
//! 精确到列的错误位置、零拷贝、以及没有构建期依赖。
//!
//! # 本层只做纯解析
//!
//! 不含任何化学语义:
//!
//! - 不做芳香性感知 —— 小写原子只是打上 [`AtomFlags::AROMATIC`],记录"作者
//!   如此声称",是否真芳香由后续的净化判定
//! - 不做价键计算、不推断隐式氢 —— `num_implicit_hs` 恒为 0
//! - 不做环感知 —— [`BondFlags::IN_RING`](omgkit_core::BondFlags::IN_RING) 不置位
//!
//! # 约定一:环闭合键延后追加,按环标号排序
//!
//! 环闭合键**全部推到键表末尾**,不在书写位置就地插入。`C1CC2CCC1CC2` 的
//! 键序因此是 `(0,1)(1,2)(2,3)(3,4)(4,5)(5,6)(6,7)(5,0)(7,2)` —— 两条环键
//! 在最后,即使 `(5,0)` 在书写上早于 `(5,6)`。
//!
//! 环键之间按**环标号升序**排列,而不是闭合先后:
//!
//! ```text
//! C2CC2C1CC1        环键 = [(5,3), (2,0)]   标号 1 在前,尽管标号 2 先闭合
//! C%10CC1CCC%10CC1  环键 = [(7,2), (5,0)]   1 < 10
//! C1CC1C1CC1        环键 = [(2,0), (5,3)]   同标号复用则按出现顺序
//! ```
//!
//! 端点方向与直觉相反:**`begin` 是闭合原子,`end` 是开环原子**
//! (`C1CCCCC1` → `(5,0)`)。唯一的例外是开环端写了键级符号,那时两端对调 ——
//! 见 `Parser::ring_closure`。
//!
//! # 约定二:手性标记相对**存储序**
//!
//! 由于约定一,原子的邻居存储顺序偏离书写顺序,手性标记必须补偿:
//!
//! ```text
//! nSwaps = 置换宇称(存储的键顺序 ↔ 书写的键顺序)   // 只算键,不含隐式氢
//! if needs_tag_inversion(...) { nSwaps += 1 }
//! if nSwaps 为奇 { 翻转标记 }
//! ```
//!
//! 其中 `needs_tag_inversion` 为:
//!
//! ```text
//! degree == 3 && ( (是片段首原子 && 显式氢数 == 1)
//!               || (显式氢数 != 1 && 环闭合数 == 1 && !不饱和) )
//! ```
//!
//! **隐式氢不参与置换**,这一点强烈反直觉。把 H 当成一个槽位塞进置换的两种
//! 自然做法都是错的:
//!
//! | 用例 | 正确结果 | "H 排最后" | "H 在 index 1" |
//! |---|---|---|---|
//! | `[C@H](N)(O)F` | 翻转 | ✓ | ✓ |
//! | `N[C@H](O)F` | 不翻 | ✓ | ✓ |
//! | `[C@@H]1CCCCC1O` | 翻转 | ✓ | ✗ |
//! | `C[P@H]C` | 不翻 | ✗ | ✓ |
//!
//! 两种建模各能解释一半用例。正确的模型是把 H 完全排除在置换之外,改由
//! 上面那个 `degree == 3` 的特判补偿。
//!
//! 书写序按**邻居原子序号**排序重建:非环键照序号排,环闭合键整体插在"自身"
//! 所处的位置。DFS 下前驱原子序号必小于本原子、后续原子序号必大于本原子,
//! 所以序号序等价于书写序。
//!
//! 之所以用"相对存储序"而不是"相对书写序":分子会被反应编辑,那时书写顺序
//! 早已不存在,只有相对存储序的约定能活下来。

mod write;

pub use write::{write, write_with_priority, write_with_priority_styled, WriteStyle, Written};

use std::collections::{HashMap, HashSet};

use omgkit_core::{AtomData, AtomFlags, BondData, BondDirection, BondOrder, ChiralTag, MolBuilder};

use crate::error::{ParseError, ParseErrorKind as K, Result};

/// 解析一条 SMILES。
///
/// 返回**未净化**的分子:芳香性未感知、隐式氢未推断、环未感知。
///
/// # Errors
/// 语法错误时返回带位置的 [`ParseError`];用 [`ParseError::render`] 可得到
/// 带插字号的两行视图。
pub fn parse(input: &str) -> Result<MolBuilder> {
    Parser::new(input.as_bytes()).run()
}

/// 解析 `.smi` 的一行:`SMILES[<空白>名字]`。
///
/// # Errors
/// 同 [`parse`]。
pub fn parse_line(line: &str) -> Result<MolBuilder> {
    let line = line.trim();
    let (smi, name) = match line.find(char::is_whitespace) {
        Some(i) => (&line[..i], line[i..].trim()),
        None => (line, ""),
    };
    let mut mol = parse(smi)?;
    if !name.is_empty() {
        mol.set_name(name);
    }
    Ok(mol)
}

// ---------------------------------------------------------------------------

/// 一个待配对的环闭合。
#[derive(Debug, Clone, Copy)]
struct RingOpen {
    atom: u32,
    order: Option<BondOrder>,
    /// 开环端是否写了**键级**符号(`- = # $ :`)。
    /// `/` `\` 只是方向符号,不算 —— 这个区分决定环键端点的朝向,见模块文档。
    explicit_order: bool,
    /// 开环端写的是 `<-`
    swap_ends: bool,
    dir: BondDirection,
    /// 开环流水号,用于把最终键下标关联回手性记录的 `ring_seqs`
    seq: u32,
    /// 开环处的位置,报错用
    pos: usize,
}

/// 已配对、待追加到键表末尾的环键。
///
/// 追加顺序按 `(number, seq)` 排序 —— 见模块文档约定 1。
#[derive(Debug, Clone, Copy)]
struct RingBond {
    /// 环标号(`1`、`%10`、`%(123)` 的数值)
    number: u32,
    /// 闭合发生的先后,同标号复用时用于稳定排序
    seq: u32,
    /// 开环流水号,用于把最终键下标回填给手性记录
    open_seq: u32,
    /// 闭合原子(放在 begin,见模块文档约定一)
    begin: u32,
    /// 开环原子
    end: u32,
    order: BondOrder,
    dir: BondDirection,
}

/// 尚未消费的键符号。
#[derive(Debug, Clone, Copy)]
struct Pending {
    /// 该符号指定的键级。`/` `\` 是纯方向标记,不指定键级,故为 `None`
    /// (键级仍走默认规则,两端芳香则为芳香键)。
    order: Option<BondOrder>,
    /// 是否是**键级**符号。`/` `\` 为 `false` —— 这决定环键端点朝向。
    explicit_order: bool,
    /// 配位键写成 `<-` 时置位:电子对由**右**侧原子提供,故端点要对调。
    /// `->` 与其余符号为 `false`。
    swap_ends: bool,
    dir: BondDirection,
    pos: usize,
}

/// 需要在解析结束后修正宇称的手性原子。
#[derive(Debug, Clone)]
struct ChiralRec {
    atom: u32,
    tag: ChiralTag,
    /// 是否是某个 `.` 片段的首原子;手性宇称补偿要用
    is_smiles_start: bool,
    /// 方括号中书写的氢数
    num_explicit_hs: u8,
    /// 该原子上的环闭合,按在此原子处书写的先后记录开环流水号。
    /// 开环端与闭合端都要记。
    ring_seqs: Vec<u32>,
}

struct Parser<'a> {
    src: &'a [u8],
    pos: usize,
    mol: MolBuilder,
    prev: Option<u32>,
    /// 分支栈:(返回到的原子, 入栈时的原子数 —— 用于检测空分支)
    branches: Vec<(u32, usize)>,
    rings: HashMap<u32, RingOpen>,
    ring_bonds: Vec<RingBond>,
    ring_seq: u32,
    pending: Option<Pending>,
    chiral: Vec<ChiralRec>,
    chiral_of: HashMap<u32, usize>,
    /// 开环流水号 → 最终键下标,在环键追加后填好
    ring_bond_index: HashMap<u32, u32>,
    /// 待建环键的端点对(已归一为 `(min, max)`),供 [`Parser::bond_exists`] O(1) 查重。
    ///
    /// 环键要到解析末尾才统一追加进分子,所以在那之前它们不在 `mol.bonds()` 里。
    /// 线性扫 `ring_bonds` 的话,每次闭环都是 O(已有环键数),环多的分子直接
    /// 退化成平方。
    pending_ring_pairs: HashSet<(u32, u32)>,
}

/// 把一对端点归一成 `(min, max)`,使查重与书写方向无关。
fn endpoint_pair(a: u32, b: u32) -> (u32, u32) {
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}

impl<'a> Parser<'a> {
    fn new(src: &'a [u8]) -> Self {
        Self {
            src,
            pos: 0,
            mol: MolBuilder::new(),
            prev: None,
            branches: Vec::new(),
            rings: HashMap::new(),
            ring_bonds: Vec::new(),
            ring_seq: 0,
            pending: None,
            chiral: Vec::new(),
            chiral_of: HashMap::new(),
            ring_bond_index: HashMap::new(),
            pending_ring_pairs: HashSet::new(),
        }
    }

    fn err<T>(&self, kind: K, pos: usize) -> Result<T> {
        Err(ParseError::new(kind, pos, self.src))
    }

    fn peek(&self) -> Option<u8> {
        self.src.get(self.pos).copied()
    }

    fn bump(&mut self) -> Option<u8> {
        let b = self.peek();
        if b.is_some() {
            self.pos += 1;
        }
        b
    }

    fn eat(&mut self, b: u8) -> bool {
        if self.peek() == Some(b) {
            self.pos += 1;
            true
        } else {
            false
        }
    }

    // -- 主循环 ------------------------------------------------------------

    fn run(mut self) -> Result<MolBuilder> {
        if self.src.is_empty() {
            return self.err(K::Empty, 0);
        }

        while let Some(b) = self.peek() {
            match b {
                b'(' => self.open_branch()?,
                b')' => self.close_branch()?,
                b'[' => self.bracket_atom()?,
                b'.' => {
                    if let Some(p) = self.pending {
                        return self.err(K::DanglingBond, p.pos);
                    }
                    self.pos += 1;
                    self.prev = None;
                }
                b'-' | b'=' | b'#' | b'$' | b':' | b'/' | b'\\' | b'<' => self.bond_symbol()?,
                b'0'..=b'9' => {
                    let pos = self.pos;
                    let n = u32::from(self.bump().expect("已 peek") - b'0');
                    self.ring_closure(n, pos)?;
                }
                b'%' => {
                    let pos = self.pos;
                    self.pos += 1;
                    let n = self.ring_number()?;
                    self.ring_closure(n, pos)?;
                }
                _ => self.organic_atom()?,
            }
        }

        // -- 收尾校验 --
        if let Some(p) = self.pending {
            return self.err(K::DanglingBond, p.pos);
        }
        if let Some(&(_, _)) = self.branches.last() {
            let pos = self.pos;
            return self.err(K::UnbalancedParen, pos);
        }
        // 注意:取 `HashMap` 的任意一项会让报错随运行而变。多个环都没闭合时,
        // 报**最先出现**的那个 —— 既确定,也最贴近用户视角。
        if let Some((&num, open)) = self.rings.iter().min_by_key(|(_, o)| o.pos) {
            let pos = open.pos;
            return self.err(K::UnclosedRingBond(num), pos);
        }
        if self.mol.num_atoms() == 0 {
            return self.err(K::Empty, 0);
        }

        // 环键在此统一追加,见模块文档约定一。
        // 排序键是 (环标号, 闭合先后)。
        let mut ring_bonds = std::mem::take(&mut self.ring_bonds);
        ring_bonds.sort_by_key(|r| (r.number, r.seq));
        let (pos, src) = (self.pos, self.src);
        for rb in ring_bonds {
            let mut bd = BondData::new(rb.begin, rb.end, rb.order);
            bd.direction = rb.dir;
            mark_aromatic(&mut bd);
            let idx = self
                .mol
                .add_bond_data(bd)
                .map_err(|_| ParseError::new(K::RingBondToSelf(rb.number), pos, src))?;
            self.ring_bond_index.insert(rb.open_seq, idx);
        }

        self.fix_chirality();
        Ok(self.mol)
    }

    // -- 分支 --------------------------------------------------------------

    fn open_branch(&mut self) -> Result<()> {
        let pos = self.pos;
        self.pos += 1;
        match self.prev {
            Some(a) => {
                self.branches.push((a, self.mol.num_atoms()));
                Ok(())
            }
            None => self.err(K::UnbalancedParen, pos),
        }
    }

    fn close_branch(&mut self) -> Result<()> {
        let pos = self.pos;
        self.pos += 1;
        match self.branches.pop() {
            Some((atom, n_at_open)) => {
                if self.mol.num_atoms() == n_at_open {
                    return self.err(K::EmptyBranch, pos);
                }
                if let Some(p) = self.pending {
                    return self.err(K::DanglingBond, p.pos);
                }
                self.prev = Some(atom);
                Ok(())
            }
            None => self.err(K::UnbalancedParen, pos),
        }
    }

    // -- 键符号 ------------------------------------------------------------

    fn bond_symbol(&mut self) -> Result<()> {
        let pos = self.pos;
        if self.pending.is_some() {
            return self.err(K::DanglingBond, pos);
        }
        let b = self.bump().expect("已 peek");
        // `/` `\` 是**纯方向**标记,不指定键级 —— 键级仍走默认规则。
        //
        // OpenSMILES 把方向键描述为"单键",照字面实现会错:
        // `Cc1cs/c(=N\C(C)C)/n1...` 中 `s/c` 这条键应当是 AROMATIC,
        // 因为两端都是芳香原子。误当单键会造成几十条差分分歧。
        //
        // 这个区分还决定环闭合键的端点朝向(见模块文档约定 1):
        // `C=1CCCCC1` 交换端点,`C/1CCCCC1` 不交换。
        let (order, explicit_order, swap_ends, dir) = match b {
            // `->` 是配位键,`-` 单独出现才是单键
            b'-' => {
                if self.eat(b'>') {
                    (Some(BondOrder::Dative), true, false, BondDirection::None)
                } else {
                    (Some(BondOrder::Single), true, false, BondDirection::None)
                }
            }
            // `<` 只可能是 `<-` 的开头
            b'<' => {
                if !self.eat(b'-') {
                    return self.err(K::UnexpectedChar('<'), pos);
                }
                (Some(BondOrder::Dative), true, true, BondDirection::None)
            }
            b'=' => (Some(BondOrder::Double), true, false, BondDirection::None),
            b'#' => (Some(BondOrder::Triple), true, false, BondDirection::None),
            b'$' => (Some(BondOrder::Quadruple), true, false, BondDirection::None),
            b':' => (Some(BondOrder::Aromatic), true, false, BondDirection::None),
            b'/' => (None, false, false, BondDirection::UpRight),
            b'\\' => {
                // `\\`(两个字面反斜杠)当作单个方向键 —— 有些 SMILES 写出
                // 工具会转义反斜杠,故 `F/C=C\\F` 与 `F/C=C\F` 等价。
                // 真实语料里这种写法并不罕见,不兼容会直接
                // 造成上百条差分分歧。
                self.eat(b'\\');
                (None, false, false, BondDirection::DownRight)
            }
            _ => unreachable!("由调用方 match 保证"),
        };
        self.pending = Some(Pending {
            order,
            explicit_order,
            swap_ends,
            dir,
            pos,
        });
        Ok(())
    }

    // -- 原子 --------------------------------------------------------------

    /// 有机子集原子(可不写方括号):B C N O P S F Cl Br I + 芳香 b c n o p s
    fn organic_atom(&mut self) -> Result<()> {
        let pos = self.pos;
        let b = self.peek().expect("已 peek");

        // 双字符元素必须先于单字符匹配,否则 Cl 会被读成 C 再撞上 l
        let two = self.src.get(self.pos + 1).copied();
        let (z, len, aromatic) = match (b, two) {
            (b'C', Some(b'l')) => (17u8, 2usize, false),
            (b'B', Some(b'r')) => (35, 2, false),
            (b'B', _) => (5, 1, false),
            (b'C', _) => (6, 1, false),
            (b'N', _) => (7, 1, false),
            (b'O', _) => (8, 1, false),
            (b'P', _) => (15, 1, false),
            (b'S', _) => (16, 1, false),
            (b'F', _) => (9, 1, false),
            (b'I', _) => (53, 1, false),
            (b'b', _) => (5, 1, true),
            (b'c', _) => (6, 1, true),
            (b'n', _) => (7, 1, true),
            (b'o', _) => (8, 1, true),
            (b'p', _) => (15, 1, true),
            (b's', _) => (16, 1, true),
            (b'*', _) => (0, 1, false),
            _ => {
                let c = char::from(b);
                return self.err(K::UnexpectedChar(c), pos);
            }
        };
        self.pos += len;

        let mut atom = AtomData::new(z);
        if aromatic {
            atom.flags.insert(AtomFlags::AROMATIC);
        }
        // 有机子集原子不写方括号 → 隐式氢由 L2 推断,此处不置 NO_IMPLICIT
        self.push_atom(atom, None, pos)
    }

    /// `[` isotope? symbol chiral? hcount? charge? (`:` class)? `]`
    fn bracket_atom(&mut self) -> Result<()> {
        let open_pos = self.pos;
        self.pos += 1; // '['

        let isotope = self.opt_number(5)?;

        // -- 元素符号 --
        let sym_pos = self.pos;
        let (z, aromatic) = self.bracket_symbol(sym_pos)?;

        // -- 立体 --
        let (tag, stereo_perm) = self.opt_chirality()?;

        // -- 氢数 --
        let mut hcount = 0u8;
        if self.peek() == Some(b'H') {
            self.pos += 1;
            hcount = match self.opt_number(2)? {
                Some(n) => u8::try_from(n)
                    .map_err(|_| ParseError::new(K::NumberOverflow, self.pos, self.src))?,
                None => 1,
            };
        }

        // -- 电荷 --
        let charge = self.opt_charge()?;

        // -- 原子映射号 --
        let mut map = 0u16;
        if self.eat(b':') {
            let n = self.opt_number(5)?.ok_or_else(|| {
                ParseError::new(K::BadBracketAtom("`:` 后缺少映射号"), self.pos, self.src)
            })?;
            map = u16::try_from(n)
                .map_err(|_| ParseError::new(K::NumberOverflow, self.pos, self.src))?;
        }

        if !self.eat(b']') {
            let kind = if self.peek().is_none() {
                K::UnexpectedEnd
            } else {
                K::BadBracketAtom("方括号未正确闭合")
            };
            return self.err(kind, self.pos);
        }

        let mut atom = AtomData::new(z);
        atom.isotope = isotope
            .map(|n| u16::try_from(n).unwrap_or(u16::MAX))
            .unwrap_or(0);
        atom.num_explicit_hs = hcount;
        atom.formal_charge = charge;
        atom.atom_map = map;
        atom.chiral_tag = tag;
        atom.stereo_perm = stereo_perm;
        // 方括号原子的氢数由作者显式给定,后续不再推断。
        // 这对**所有**方括号原子成立,包括 `[CH4]` `[Fe+2]` `[*:1]`
        atom.flags.insert(AtomFlags::NO_IMPLICIT);
        if aromatic {
            atom.flags.insert(AtomFlags::AROMATIC);
        }

        self.push_atom(atom, Some((tag, hcount)), open_pos)
    }

    fn bracket_symbol(&mut self, pos: usize) -> Result<(u8, bool)> {
        let b = match self.peek() {
            Some(b) => b,
            None => return self.err(K::UnexpectedEnd, pos),
        };

        if b == b'*' {
            self.pos += 1;
            return Ok((0, false));
        }

        if b.is_ascii_uppercase() {
            // 贪心取两字符,失败再退回一字符 —— `[Cl]` 与 `[C@]` 都要正确
            let mut sym = String::new();
            sym.push(char::from(b));
            let second = self.src.get(self.pos + 1).copied();
            if let Some(c) = second {
                if c.is_ascii_lowercase() {
                    let two: String = [char::from(b), char::from(c)].iter().collect();
                    if let Some(z) = omgkit_core::element::atomic_num_of(&two) {
                        self.pos += 2;
                        return Ok((z, false));
                    }
                }
            }
            if let Some(z) = omgkit_core::element::atomic_num_of(&sym) {
                self.pos += 1;
                return Ok((z, false));
            }
            if let Some(c) = second {
                if c.is_ascii_lowercase() {
                    sym.push(char::from(c));
                }
            }
            return self.err(K::UnknownElement(sym), pos);
        }

        if b.is_ascii_lowercase() {
            // 芳香形式:单字符 b c n o p s,双字符 se as te
            let second = self.src.get(self.pos + 1).copied();
            if let Some(c) = second {
                if c.is_ascii_lowercase() {
                    let two: String = [char::from(b), char::from(c)].iter().collect();
                    let up = format!("{}{}", two[..1].to_ascii_uppercase(), &two[1..]);
                    if matches!(two.as_str(), "se" | "as" | "te" | "si") {
                        if let Some(z) = omgkit_core::element::atomic_num_of(&up) {
                            self.pos += 2;
                            return Ok((z, true));
                        }
                    }
                }
            }
            let up = char::from(b).to_ascii_uppercase().to_string();
            if let Some(z) = omgkit_core::element::atomic_num_of(&up) {
                if omgkit_core::element::can_be_aromatic_lowercase(z) {
                    self.pos += 1;
                    return Ok((z, true));
                }
            }
            return self.err(K::UnknownElement(char::from(b).to_string()), pos);
        }

        self.err(K::BadBracketAtom("缺少元素符号"), pos)
    }

    /// 读立体标记,返回(几何类别, 类内排列序号)。
    ///
    /// 支持简写 `@` / `@@` 与扩展形式 `@TH1` `@AL1` `@SP1` `@TB15` `@OH25`。
    ///
    /// # 序号的取值范围要当场校验
    ///
    /// 每种几何的排列数是**几何本身决定的**:配体在多面体位置上的排法除以
    /// 该多面体的转动群阶数 —— 平面四方 4!/8 = 3,三角双锥 5!/6 = 20,
    /// 八面体 6!/24 = 30。所以 `@TB21` 不是"某个暂时不认识的排列",而是一个
    /// **不存在**的排列,只能是笔误。放它过去等于把一个无法解释的数字
    /// 塞进分子里,越往后越难查。
    ///
    /// 序号 0 与省略序号(光秃的 `@SP`)都表示"有这个几何但没指定排列",
    /// 是合法的。
    ///
    /// # 四面体的序号不单独存
    ///
    /// `TH1` 就是 `@`、`TH2` 就是 `@@`,排列已由类别本身表达,再存一份序号
    /// 就有了两个可以互相矛盾的真相来源。故四面体一律返回 0。
    fn opt_chirality(&mut self) -> Result<(ChiralTag, u8)> {
        if self.peek() != Some(b'@') {
            return Ok((ChiralTag::Unspecified, 0));
        }
        self.pos += 1;
        if self.eat(b'@') {
            return Ok((ChiralTag::Cw, 0));
        }
        let rest = &self.src[self.pos..];
        let two = rest.get(..2).map(<[u8]>::to_ascii_uppercase);
        // (几何书写形式, 最大序号, 类别)
        let (name, max, class) = match two.as_deref() {
            Some(b"TH") => ("TH", 2, ChiralTag::Ccw),
            // 丙二烯轴手性:立体信息属于一根轴而非一个配位中心,不归入配位几何
            Some(b"AL") => ("AL", 2, ChiralTag::Other),
            Some(b"SP") => ("SP", 3, ChiralTag::SquarePlanar),
            Some(b"TB") => ("TB", 20, ChiralTag::TrigonalBipyramidal),
            Some(b"OH") => ("OH", 30, ChiralTag::Octahedral),
            // 光秃秃的 `@` 就是四面体逆时针
            _ => return Ok((ChiralTag::Ccw, 0)),
        };
        self.pos += 2;
        let pos = self.pos;
        let perm = self.opt_number(2)?.unwrap_or(0);
        if perm > max {
            return self.err(
                K::StereoPermOutOfRange {
                    geometry: name,
                    got: perm,
                    max,
                },
                pos,
            );
        }
        if name == "TH" {
            // TH1 = `@` = 逆时针,TH2 = `@@` = 顺时针
            return Ok((
                if perm == 2 {
                    ChiralTag::Cw
                } else {
                    ChiralTag::Ccw
                },
                0,
            ));
        }
        Ok((class, u8::try_from(perm).unwrap_or(0)))
    }

    fn opt_charge(&mut self) -> Result<i8> {
        let sign = match self.peek() {
            Some(b'+') => 1i8,
            Some(b'-') => -1i8,
            _ => return Ok(0),
        };
        self.pos += 1;

        // `++` / `--` 形式
        let mut n = 1i32;
        while self.peek() == Some(if sign > 0 { b'+' } else { b'-' }) {
            self.pos += 1;
            n += 1;
        }
        if n > 1 {
            return i8::try_from(n * i32::from(sign))
                .map_err(|_| ParseError::new(K::NumberOverflow, self.pos, self.src));
        }

        // `+2` 形式
        if let Some(v) = self.opt_number(2)? {
            n = i32::try_from(v).unwrap_or(i32::MAX);
        }
        i8::try_from(n * i32::from(sign))
            .map_err(|_| ParseError::new(K::NumberOverflow, self.pos, self.src))
    }

    /// 读至多 `max_digits` 位十进制数;无数字时返回 `None`。
    fn opt_number(&mut self, max_digits: usize) -> Result<Option<u32>> {
        let start = self.pos;
        let mut val: u64 = 0;
        let mut n = 0;
        while n < max_digits {
            match self.peek() {
                Some(d @ b'0'..=b'9') => {
                    val = val * 10 + u64::from(d - b'0');
                    self.pos += 1;
                    n += 1;
                }
                _ => break,
            }
        }
        if n == 0 {
            return Ok(None);
        }
        u32::try_from(val)
            .map(Some)
            .map_err(|_| ParseError::new(K::NumberOverflow, start, self.src))
    }

    /// `%NN` 或 `%(N...)`
    fn ring_number(&mut self) -> Result<u32> {
        if self.eat(b'(') {
            let n = self.opt_number(5)?.ok_or_else(|| {
                ParseError::new(K::BadBracketAtom("`%(` 后缺少数字"), self.pos, self.src)
            })?;
            if !self.eat(b')') {
                return self.err(K::BadBracketAtom("`%(` 未闭合"), self.pos);
            }
            return Ok(n);
        }
        let pos = self.pos;
        match self.opt_number(2)? {
            Some(n) if self.pos - pos == 2 => Ok(n),
            _ => self.err(K::BadBracketAtom("`%` 后需要两位数字"), pos),
        }
    }

    // -- 原子入图 ----------------------------------------------------------

    /// 把原子加入分子,并处理与前一个原子之间的键。
    ///
    /// `chiral` 为 `Some((tag, has_h))` 时登记手性记录。
    fn push_atom(
        &mut self,
        atom: AtomData,
        chiral: Option<(ChiralTag, u8)>,
        pos: usize,
    ) -> Result<()> {
        let aromatic = atom.flags.contains(AtomFlags::AROMATIC);
        let idx = self.mol.add_atom_data(atom);

        // 手性记录须在成键前建立 —— `is_smiles_start` 要看的是"此刻还没有前驱原子"
        if let Some((tag, num_explicit_hs)) = chiral {
            if tag != ChiralTag::Unspecified {
                self.chiral_of.insert(idx, self.chiral.len());
                self.chiral.push(ChiralRec {
                    atom: idx,
                    tag,
                    is_smiles_start: self.prev.is_none(),
                    num_explicit_hs,
                    ring_seqs: Vec::new(),
                });
            }
        }

        if let Some(prev) = self.prev {
            let pending = self.pending.take();
            let order = match pending.and_then(|p| p.order) {
                Some(o) => o,
                None => self.default_order(prev, aromatic),
            };
            // 配位键 `<-` 的电子对由右侧原子提供,故 begin 是**后**写的那个
            let (begin, end) = if pending.is_some_and(|p| p.swap_ends) {
                (idx, prev)
            } else {
                (prev, idx)
            };
            let mut bd = BondData::new(begin, end, order);
            bd.direction = pending.map_or(BondDirection::None, |p| p.dir);
            mark_aromatic(&mut bd);
            let src = self.src;
            self.mol
                .add_bond_data(bd)
                .map_err(|_| ParseError::new(K::UnexpectedChar('?'), pos, src))?;
        } else if let Some(p) = self.pending.take() {
            // 键符号后面没有可连接的前驱原子
            return self.err(K::DanglingBond, p.pos);
        }

        self.prev = Some(idx);
        Ok(())
    }

    /// 未显式书写键符号时的默认键级:两端都芳香则芳香键,否则单键。
    ///
    /// `c1ccccc1` 的键全为 AROMATIC,而
    /// `[O-][N+](=O)c1ccccc1` 中 N→c 的键为 SINGLE。
    fn default_order(&self, prev: u32, cur_aromatic: bool) -> BondOrder {
        let prev_aromatic = self.mol.atoms()[prev as usize]
            .flags
            .contains(AtomFlags::AROMATIC);
        if prev_aromatic && cur_aromatic {
            BondOrder::Aromatic
        } else {
            BondOrder::Single
        }
    }

    /// 若 `atom` 是手性原子,记下它这里出现的一次环闭合(按书写先后)。
    ///
    /// 开环端与闭合端都要记 —— 宇称补偿要用到两端各自的环闭合数。
    fn note_ring(&mut self, atom: u32, open_seq: u32) {
        if let Some(&i) = self.chiral_of.get(&atom) {
            self.chiral[i].ring_seqs.push(open_seq);
        }
    }

    // -- 环闭合 ------------------------------------------------------------

    fn ring_closure(&mut self, num: u32, pos: usize) -> Result<()> {
        let cur = match self.prev {
            Some(a) => a,
            None => return self.err(K::UnexpectedChar(char::from(self.src[pos])), pos),
        };
        let pending = self.pending.take();

        match self.rings.remove(&num) {
            // 闭环
            Some(open) => {
                if open.atom == cur {
                    return self.err(K::RingBondToSelf(num), pos);
                }
                if self.bond_exists(open.atom, cur) {
                    return self.err(K::DuplicateRingBond(num), pos);
                }

                let close_order = pending.and_then(|p| p.order);
                let close_explicit = pending.is_some_and(|p| p.explicit_order);

                // 只有两端都写了**键级**符号才可能冲突;`/` 蕴含的单键不参与
                if open.explicit_order && close_explicit && open.order != close_order {
                    return self.err(K::ConflictingRingBondOrder(num), pos);
                }
                let order = match open.order.or(close_order) {
                    Some(o) => o,
                    None => {
                        let cur_arom = self.mol.atoms()[cur as usize]
                            .flags
                            .contains(AtomFlags::AROMATIC);
                        self.default_order(open.atom, cur_arom)
                    }
                };

                // 端点朝向:开环端写了键级符号时,键在开环处即已确定,故开环
                // 原子在 begin;否则键在闭合处建立,闭合原子在 begin:
                //   C1CCCCC1 → (5,0)   C=1CCCCC1 → (0,5)   C/1CCCCC1 → (5,0)
                let open_is_begin = open.explicit_order;
                let (begin, end) = if open_is_begin {
                    (open.atom, cur)
                } else {
                    (cur, open.atom)
                };

                // `<-` 在上面的朝向之上再对调一次。取的是**定了键级的那一端**
                // 写的符号 —— 与 `order` 的取法保持一致(开环端优先),否则
                // 两端都写 `->` 时键级和方向会来自不同的端。
                //   N->1CCCCC->1   → (0,5)   开环端说 N 给电子
                //   [Cu]<-1CCCC1   → (4,0)   开环端说 Cu 收电子
                //   [Cu]1<-NCCN->1 → (4,0)   开环端没定键级,听闭合端的
                let swap_ends = if open.order.is_some() {
                    open.swap_ends
                } else {
                    pending.is_some_and(|p| p.swap_ends)
                };
                // 方向符号处于书写它那一端的参考系;该端若被存为 `end` 就要翻转。
                let close_dir = pending.map_or(BondDirection::None, |p| p.dir);
                let dir = if open_is_begin {
                    if open.dir != BondDirection::None {
                        open.dir
                    } else {
                        close_dir.flipped()
                    }
                } else if close_dir != BondDirection::None {
                    close_dir
                } else {
                    open.dir.flipped()
                };

                // 端点对调后,方向也处在了相反的参考系里,要跟着翻。
                // 一端写箭头、另一端写 `/` 或 `\` 时才看得出来 ——
                // `N/1CCCCC<-1` 的方向应是 `/`,漏了这一翻就成了 `\`。
                let (begin, end, dir) = if swap_ends {
                    (end, begin, dir.flipped())
                } else {
                    (begin, end, dir)
                };

                let seq = u32::try_from(self.ring_bonds.len()).unwrap_or(u32::MAX);
                self.pending_ring_pairs.insert(endpoint_pair(begin, end));
                self.ring_bonds.push(RingBond {
                    number: num,
                    seq,
                    open_seq: open.seq,
                    begin,
                    end,
                    order,
                    dir,
                });

                // 闭合端也要记这次环闭合(开环端在开环时已记)
                self.note_ring(cur, open.seq);
            }
            // 开环
            None => {
                let seq = self.ring_seq;
                self.ring_seq += 1;
                self.rings.insert(
                    num,
                    RingOpen {
                        atom: cur,
                        order: pending.and_then(|p| p.order),
                        explicit_order: pending.is_some_and(|p| p.explicit_order),
                        swap_ends: pending.is_some_and(|p| p.swap_ends),
                        dir: pending.map_or(BondDirection::None, |p| p.dir),
                        seq,
                        pos,
                    },
                );
                self.note_ring(cur, seq);
            }
        }
        Ok(())
    }

    /// `a` 与 `b` 之间是否已有键 —— 既查已建的,也查待建的环键。
    ///
    /// 两边都是 O(1)/O(度数):分子侧走 `MolBuilder` 的邻接索引,
    /// 待建侧走 [`Parser::pending_ring_pairs`]。
    fn bond_exists(&self, a: u32, b: u32) -> bool {
        self.mol.bond_between(a, b).is_some()
            || self.pending_ring_pairs.contains(&endpoint_pair(a, b))
    }

    // -- 手性宇称修正 ------------------------------------------------------

    /// 把每个手性标记从"书写序"转换到"存储序"。
    ///
    /// 算法见模块文档约定二:
    ///
    /// ```text
    /// nSwaps = 置换宇称(存储的键顺序 ↔ 书写的键顺序)   // 只算键,不含隐式氢
    /// if 需要补偿(见模块文档) { nSwaps += 1 }
    /// if nSwaps 为奇 { 翻转标记 }
    /// ```
    ///
    /// 隐式氢**不参与置换**,而是由那条 degree==3 的补偿规则统一处理 ——
    /// 见模块文档里的对照表。
    ///
    /// 只有四面体能这样修:它的两种排列互为镜像,一次对换即翻转。配位几何
    /// (SP/TB/OH)的排列序号在邻居重排下按查找表变换,不是一个可翻转的
    /// 布尔量,处理方式见下方注释。
    fn fix_chirality(&mut self) {
        for rec in &self.chiral {
            let atom = rec.atom;

            // 该原子上的环闭合键下标(按在此原子处书写的先后)
            let ring_bonds: Vec<u32> = rec
                .ring_seqs
                .iter()
                .filter_map(|s| self.ring_bond_index.get(s).copied())
                .collect();

            // 存储序:按键下标递增
            let stored: Vec<u32> = self
                .mol
                .bonds()
                .iter()
                .enumerate()
                .filter(|(_, b)| b.other_end(atom).is_some())
                .map(|(i, _)| i as u32)
                .collect();

            // 书写序:非环键按**邻居原子序号**排序,
            // 环闭合键整体插在"自身"所处的位置。这样能重建 SMILES 的书写顺序,
            // 因为 DFS 下前驱原子序号必小于本原子、后续原子序号必大于本原子。
            let mut entries: Vec<(u32, Option<u32>)> = vec![(atom, None)];
            for (i, b) in self.mol.bonds().iter().enumerate() {
                let i = i as u32;
                if ring_bonds.contains(&i) {
                    continue;
                }
                if let Some(other) = b.other_end(atom) {
                    entries.push((other, Some(i)));
                }
            }
            entries.sort_by_key(|e| e.0);

            let mut written = Vec::with_capacity(stored.len());
            for (_, bond) in entries {
                match bond {
                    None => written.extend(ring_bonds.iter().copied()),
                    Some(i) => written.push(i),
                }
            }
            if written.len() != stored.len() {
                continue; // 结构异常,已在别处报错
            }

            if !rec.tag.is_tetrahedral() {
                // 配位几何的排列序号换参考系要走查找表(见 `AtomData::stereo_perm`
                // 的文档),不是一次奇偶翻转能表达的。这里不动它 —— 它保存的是
                // **书写时的字面值**,而字面值不会因为存储序变了就变错。
                continue;
            }

            let mut odd = permutation_is_odd(&written, &stored).unwrap_or(false);

            if stored.len() == 3 {
                // 此刻价键尚未计算,所以"是否有第四个价位"只能看写出来的氢数
                let has_fourth_valence = rec.num_explicit_hs == 1;
                // 不饱和:任一键的数值键级 > 1(芳香键的 1.5 也算)
                let unsaturated = self
                    .mol
                    .bonds()
                    .iter()
                    .filter(|b| b.other_end(atom).is_some())
                    .any(|b| b.order.as_double() > 1.0);

                if (rec.is_smiles_start && rec.num_explicit_hs == 1)
                    || (!has_fourth_valence && ring_bonds.len() == 1 && !unsaturated)
                {
                    odd = !odd;
                }
            }

            if odd {
                if let Some(a) = self.mol.atom_mut(atom) {
                    a.chiral_tag = rec.tag.inverted();
                }
            }
        }
    }
}

/// 芳香键的标志位与键级 AROMATIC 始终同步。
fn mark_aromatic(bd: &mut BondData) {
    if bd.order == BondOrder::Aromatic {
        bd.flags.insert(omgkit_core::BondFlags::AROMATIC);
    }
}

/// `written` → `storage` 置换的宇称。两者不是同一多重集时返回 `None`。
///
/// n ≤ 6,O(n²) 完全够用。
pub(crate) fn permutation_is_odd(written: &[u32], storage: &[u32]) -> Option<bool> {
    let n = written.len();
    if n != storage.len() {
        return None;
    }
    let mut used = vec![false; n];
    let mut perm = Vec::with_capacity(n);
    for w in written {
        let j = storage
            .iter()
            .enumerate()
            .position(|(j, s)| !used[j] && s == w)?;
        used[j] = true;
        perm.push(j);
    }
    let mut inversions = 0usize;
    for i in 0..n {
        for j in i + 1..n {
            if perm[i] > perm[j] {
                inversions += 1;
            }
        }
    }
    Some(inversions % 2 == 1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::ParseErrorKind as K;

    fn parse_ok(s: &str) -> MolBuilder {
        parse(s).unwrap_or_else(|e| panic!("应能解析 {s:?}:\n{}", e.render()))
    }

    // -- 基本形状 --

    #[test]
    fn shapes() {
        for (smi, na, nb) in [
            ("C", 1, 0),
            ("CC", 2, 1),
            ("CCO", 3, 2),
            ("CC(=O)O", 4, 3),
            ("C1CCCCC1", 6, 6),
            ("CCO.CCN", 6, 4),
            ("c1ccccc1", 6, 6),
            ("C1CC2CCC1CC2", 8, 9),
        ] {
            let m = parse_ok(smi);
            assert_eq!((m.num_atoms(), m.num_bonds()), (na, nb), "{smi}");
        }
    }

    #[test]
    fn dot_keeps_one_molecule_with_two_fragments() {
        // `.` 不切分分子,只是断开成键
        let m = parse_ok("CCO.CCN");
        assert_eq!(m.num_atoms(), 6);
        assert_eq!(m.num_bonds(), 4);
    }

    // -- 方括号原子 --

    #[test]
    fn bracket_atom_fields() {
        let m = parse_ok("[13CH4:7]");
        let a = m.atoms()[0];
        assert_eq!(a.atomic_num, 6);
        assert_eq!(a.isotope, 13);
        assert_eq!(a.num_explicit_hs, 4);
        assert_eq!(a.atom_map, 7);
        assert!(a.flags.contains(AtomFlags::NO_IMPLICIT));
    }

    #[test]
    fn charge_forms_are_equivalent() {
        assert_eq!(parse_ok("[Fe+2]").atoms()[0].formal_charge, 2);
        assert_eq!(parse_ok("[Fe++]").atoms()[0].formal_charge, 2);
        assert_eq!(parse_ok("[O-]").atoms()[0].formal_charge, -1);
        assert_eq!(parse_ok("[O-2]").atoms()[0].formal_charge, -2);
        assert_eq!(parse_ok("[O--]").atoms()[0].formal_charge, -2);
    }

    #[test]
    fn wildcard_is_atomic_number_zero() {
        assert_eq!(parse_ok("*").atoms()[0].atomic_num, 0);
        assert_eq!(parse_ok("[*:1]").atoms()[0].atom_map, 1);
    }

    #[test]
    fn two_char_element_beats_one_char() {
        // Cl 必须先于 C 匹配,否则会读成 C 再撞上 l
        let m = parse_ok("ClCCBr");
        let z: Vec<u8> = m.atoms().iter().map(|a| a.atomic_num).collect();
        assert_eq!(z, vec![17, 6, 6, 35]);
    }

    #[test]
    fn bracket_greedy_symbol_backtracks() {
        // `[C@]` 里 `C@` 不是元素,必须退回单字符 C
        assert_eq!(parse_ok("[C@](N)(O)(F)Cl").atoms()[0].atomic_num, 6);
        assert_eq!(parse_ok("[Cl-]").atoms()[0].atomic_num, 17);
    }

    // -- 芳香 --

    #[test]
    fn lowercase_atoms_are_flagged_aromatic() {
        let m = parse_ok("c1ccccc1");
        assert!(m
            .atoms()
            .iter()
            .all(|a| a.flags.contains(AtomFlags::AROMATIC)));
        assert!(m.bonds().iter().all(|b| b.order == BondOrder::Aromatic));
    }

    #[test]
    fn default_bond_is_single_when_either_end_is_not_aromatic() {
        // N→c 是单键,即使 c 芳香
        let m = parse_ok("[O-][N+](=O)c1ccccc1");
        assert_eq!(m.bonds()[2].order, BondOrder::Single);
    }

    #[test]
    fn direction_symbols_do_not_force_single_bond() {
        // `/` 是纯方向标记;两端芳香时键级仍是芳香
        let m = parse_ok("Cc1cs/c(=N\\C)/n1");
        assert!(
            m.bonds()
                .iter()
                .any(|b| b.order == BondOrder::Aromatic && b.direction != BondDirection::None),
            "应存在既芳香又带方向的键"
        );
    }

    #[test]
    fn double_backslash_is_one_direction_bond() {
        // `\\` 等价于 `\`
        let a = parse_ok(r"F/C=C\F");
        let b = parse_ok(r"F/C=C\\F");
        assert_eq!(a.num_atoms(), b.num_atoms());
        assert_eq!(a.num_bonds(), b.num_bonds());
        assert_eq!(a.bonds()[2].direction, b.bonds()[2].direction);
    }

    // -- 环闭合 --

    #[test]
    fn ring_bond_is_appended_last_with_closer_first() {
        let m = parse_ok("C1CCCCC1");
        let last = m.bonds()[5];
        assert_eq!(
            (last.begin, last.end),
            (5, 0),
            "环键端点应为 (闭合原子, 开环原子)"
        );
    }

    #[test]
    fn ring_bonds_sort_by_ring_number() {
        // ring2 先开先闭,但 ring1 的键要排在前面
        let m = parse_ok("C2CC2C1CC1");
        assert_eq!(m.num_bonds(), 7, "5 条链键 + 2 条环键");
        let ring: Vec<(u32, u32)> = m.bonds()[m.num_bonds() - 2..]
            .iter()
            .map(|b| (b.begin, b.end))
            .collect();
        assert_eq!(ring, vec![(5, 3), (2, 0)]);
    }

    #[test]
    fn explicit_order_at_ring_open_swaps_endpoints() {
        assert_eq!(
            (
                parse_ok("C1CCCCC1").bonds()[5].begin,
                parse_ok("C1CCCCC1").bonds()[5].end
            ),
            (5, 0)
        );
        let m = parse_ok("C=1CCCCC1");
        assert_eq!((m.bonds()[5].begin, m.bonds()[5].end), (0, 5));
        // 方向符号不触发交换
        let m = parse_ok("C/1CCCCC1");
        assert_eq!((m.bonds()[5].begin, m.bonds()[5].end), (5, 0));
    }

    #[test]
    fn multi_digit_ring_numbers() {
        assert_eq!(parse_ok("C%10CCCCC%10").num_bonds(), 6);
        assert_eq!(parse_ok("C%(123)CCCCC%(123)").num_bonds(), 6);
    }

    // -- 错误:位置必须精确 --

    fn err_at(smi: &str) -> ParseError {
        parse(smi).expect_err(&format!("{smi:?} 应当解析失败"))
    }

    #[test]
    fn unclosed_ring_points_at_the_digit() {
        let e = err_at("C1CC");
        assert_eq!(e.kind, K::UnclosedRingBond(1));
        assert_eq!(e.pos, 1, "应指向环标号所在位置");
        assert!(e.render().contains('^'));
    }

    #[test]
    fn ring_to_self() {
        assert_eq!(err_at("C11").kind, K::RingBondToSelf(1));
    }

    #[test]
    fn unbalanced_parens() {
        assert_eq!(err_at("CC(").kind, K::UnbalancedParen);
        assert_eq!(err_at("CC)").kind, K::UnbalancedParen);
    }

    #[test]
    fn empty_branch() {
        assert_eq!(err_at("CC()C").kind, K::EmptyBranch);
    }

    #[test]
    fn dangling_bond() {
        assert_eq!(err_at("C=").kind, K::DanglingBond);
        assert_eq!(err_at("=C").kind, K::DanglingBond);
    }

    #[test]
    fn unknown_element() {
        assert!(matches!(err_at("[Xx]").kind, K::UnknownElement(_)));
    }

    #[test]
    fn unclosed_bracket() {
        assert_eq!(err_at("[C").kind, K::UnexpectedEnd);
    }

    #[test]
    fn empty_input() {
        assert_eq!(err_at("").kind, K::Empty);
    }

    #[test]
    fn error_render_puts_caret_under_the_right_column() {
        let e = err_at("CCC1CC");
        let rendered = e.render();
        let caret_line = rendered.lines().nth(1).unwrap();
        assert_eq!(caret_line.find('^'), Some(e.pos), "插字号列号应等于 pos");
    }

    // -- 配位键 --

    /// `->` 与 `<-` 是同一条键的两种写法,区别只在端点朝向:
    /// **begin 端提供电子对**。
    #[test]
    fn dative_bond_direction_follows_the_arrow() {
        let m = parse_ok("N->[Cu]");
        let b = m.bonds()[0];
        assert_eq!(
            (b.begin, b.end, b.order),
            (0, 1, BondOrder::Dative),
            "N->[Cu]:给电子的 N 应在 begin"
        );

        let m = parse_ok("[Cu]<-N");
        let b = m.bonds()[0];
        assert_eq!(
            (b.begin, b.end, b.order),
            (1, 0, BondOrder::Dative),
            "[Cu]<-N:给电子的 N 是后写的那个,端点要对调"
        );
    }

    /// 单键的 `-` 不能被 `->` 抢走。
    #[test]
    fn single_bond_dash_is_not_swallowed_by_dative() {
        let m = parse_ok("C-C");
        assert_eq!(m.bonds()[0].order, BondOrder::Single);
        // `-` 后面接 `->` 是两个键符号连写,应报悬空键
        assert_eq!(err_at("N-->[Cu]").kind, K::DanglingBond);
    }

    /// 孤立的 `<` 不是合法符号。
    #[test]
    fn lone_angle_bracket_is_an_error() {
        assert_eq!(err_at("C<C").kind, K::UnexpectedChar('<'));
    }

    /// 配位键也能当环闭合键。朝向仍由箭头决定,且**定了键级的那一端**说了算。
    ///
    /// 端点以集合比对:环键统一追加到键表末尾,下标顺序与书写顺序本就不同。
    #[test]
    fn dative_ring_closures() {
        for (smi, expect) in [
            // 开环端写 `->`:开环原子给电子
            ("N->1CCCCC1", &[(0u32, 5u32)][..]),
            // 开环端写 `<-`:开环原子收电子,对方给
            ("[Cu]<-1CCCC1", &[(4, 0)][..]),
            // 两端都定了键级时,以开环端为准
            ("N->1CCCCC<-1", &[(0, 5)][..]),
            // 开环端只写了环标号,没定键级,于是听闭合端的;
            // 另有一条链上的配位键 `[Cu]<-N`
            ("[Cu]1<-NCCN->1", &[(1, 0), (4, 0)][..]),
        ] {
            let m = parse_ok(smi);
            let mut got: Vec<(u32, u32)> = m
                .bonds()
                .iter()
                .filter(|b| b.order == BondOrder::Dative)
                .map(|b| (b.begin, b.end))
                .collect();
            got.sort_unstable();
            let mut want = expect.to_vec();
            want.sort_unstable();
            assert_eq!(got, want, "{smi}");
        }
    }

    // -- 非四面体立体 --

    /// `@SP` / `@TB` / `@OH` 必须保留**几何类别**和**类内排列序号**。
    /// 只留一个"其它"标记的话,序号丢失,写出时就还原不回去了。
    #[test]
    fn coordination_geometry_keeps_class_and_permutation() {
        for (smi, tag, perm) in [
            ("[Pt@SP1](Cl)(Cl)(N)N", ChiralTag::SquarePlanar, 1u8),
            ("[Pt@SP3](Cl)(Cl)(N)N", ChiralTag::SquarePlanar, 3),
            ("F[P@TB15](Cl)(Br)(I)S", ChiralTag::TrigonalBipyramidal, 15),
            ("C[Co@OH25](N)(O)(S)(P)Cl", ChiralTag::Octahedral, 25),
        ] {
            let m = parse_ok(smi);
            let a = m
                .atoms()
                .iter()
                .find(|a| a.chiral_tag != ChiralTag::Unspecified)
                .unwrap_or_else(|| panic!("{smi}:应有立体标记"));
            assert_eq!((a.chiral_tag, a.stereo_perm), (tag, perm), "{smi}");
        }
    }

    /// 四面体的排列由标记自身表达,序号保持 0 —— 两处都记就会有一致性隐患。
    #[test]
    fn tetrahedral_leaves_permutation_zero() {
        for (smi, tag) in [
            ("[C@](N)(O)(F)Cl", ChiralTag::Ccw),
            ("[C@@](N)(O)(F)Cl", ChiralTag::Cw),
            ("[C@TH1](N)(O)(F)Cl", ChiralTag::Ccw),
            ("[C@TH2](N)(O)(F)Cl", ChiralTag::Cw),
        ] {
            let a = parse_ok(smi).atoms()[0];
            assert_eq!((a.chiral_tag, a.stereo_perm), (tag, 0), "{smi}");
        }
    }

    /// 丙二烯的 `@AL` 说的是一根**轴**,不是配位中心,故不归入配位几何。
    #[test]
    fn allene_is_not_a_coordination_geometry() {
        let m = parse_ok("N[C@AL1]=C=C(O)F");
        assert_eq!(m.atoms()[1].chiral_tag, ChiralTag::Other);
        assert_eq!(m.atoms()[1].stereo_perm, 1);
    }

    /// 排列序号原样保管书写值,**不**随邻居重排而改动。
    ///
    /// 换参照系要走查找表(见 `AtomData::stereo_perm`),属于 L6。在那之前
    /// 本字段的语义就是"作者写了几",这个语义不会因为存储序变了而失效。
    #[test]
    fn coordination_permutation_keeps_the_written_literal() {
        // 环闭合键被追加到键表末尾,故存储序 ≠ 书写序
        let m = parse_ok("[Co@OH5]1(N)(O)(S)(P)CCC1");
        let a = m.atoms()[0];
        assert_eq!(a.chiral_tag, ChiralTag::Octahedral, "类别要保住");
        assert_eq!(a.stereo_perm, 5, "序号是字面值,不因重排而变");

        let m = parse_ok("[Co@OH5](N)(O)(S)(P)C");
        assert_eq!(m.atoms()[0].stereo_perm, 5);
    }

    /// 排列序号超出该几何的取值范围时报错。
    ///
    /// 排列数由几何本身决定(SP 3 / TB 20 / OH 30),超出的序号不是"暂时
    /// 不认识",而是不存在 —— 放过去等于把一个无法解释的数字存进分子。
    #[test]
    fn out_of_range_stereo_permutation_is_rejected() {
        for (smi, geometry, got, max) in [
            ("[Pt@SP4](Cl)(Cl)(N)N", "SP", 4u32, 3u32),
            ("F[P@TB21](Cl)(Br)(I)S", "TB", 21, 20),
            ("C[Co@OH31](N)(O)(S)(P)Cl", "OH", 31, 30),
            ("[C@TH3](N)(O)(F)Cl", "TH", 3, 2),
            ("N[C@AL3]=C=C(O)F", "AL", 3, 2),
        ] {
            assert_eq!(
                err_at(smi).kind,
                K::StereoPermOutOfRange { geometry, got, max },
                "{smi} 应被拒绝"
            );
        }
        // 0 与省略序号都表示"有这个几何但没指定排列",合法
        assert_eq!(parse_ok("[Pt@SP0](Cl)(Cl)(N)N").atoms()[0].stereo_perm, 0);
        assert_eq!(parse_ok("[Pt@SP](Cl)(Cl)(N)N").atoms()[0].stereo_perm, 0);
        assert_eq!(
            parse_ok("[Pt@SP](Cl)(Cl)(N)N").atoms()[0].chiral_tag,
            ChiralTag::SquarePlanar
        );
    }

    /// 环闭合两端写了互相矛盾的键级 —— **刻意报错**。
    ///
    /// 这是一处有意为之的分歧:同一条键被写了两个不同的键级,只可能是笔误。
    /// 常见做法是取其中一端(通常是开环端)静默放过,那样错误会一路带到
    /// 分子里去。精确报出位置正是手写解析器的意义所在。
    ///
    /// 配位键把这条规则的触发面扩大了:`->` 与 `-` 现在也算冲突。
    #[test]
    fn conflicting_ring_bond_orders_are_rejected_on_purpose() {
        for smi in [
            "C=1CCCCC#1",
            "N-1CCCCC<-1",
            "N=1CCCCC<-1",
            "N->1CCCCC-1",
            "N->1CCCCC=1",
        ] {
            assert_eq!(
                err_at(smi).kind,
                K::ConflictingRingBondOrder(1),
                "{smi} 的两端键级矛盾,应当报错"
            );
        }
        // 两端一致(含"一端 `->` 一端 `<-`"这种朝向相反但键级相同的写法)不算冲突
        for smi in ["N->1CCCCC->1", "N->1CCCCC<-1", "C=1CCCCC=1"] {
            let _ = parse_ok(smi);
        }
    }

    /// 端点对调时方向符号也要跟着翻。
    ///
    /// 一端写箭头、另一端写 `/` 或 `\` 才看得出来 —— 箭头自身不带方向。
    #[test]
    fn dative_ring_closure_flips_the_direction_too() {
        for (smi, begin, end, dir) in [
            ("N/1CCCCC<-1", 0u32, 5u32, BondDirection::UpRight),
            ("N<-1CCCCC/1", 5, 0, BondDirection::UpRight),
            ("N\\1CCCCC<-1", 0, 5, BondDirection::DownRight),
        ] {
            let m = parse_ok(smi);
            let b = m
                .bonds()
                .iter()
                .find(|b| b.order == BondOrder::Dative)
                .unwrap_or_else(|| panic!("{smi}:应有一条配位键"));
            assert_eq!((b.begin, b.end, b.direction), (begin, end, dir), "{smi}");
        }
    }

    // -- .smi 行 --

    #[test]
    fn parse_line_takes_name() {
        let m = parse_line("CCO\tethanol").unwrap();
        assert_eq!(m.num_atoms(), 3);
        assert_eq!(m.name(), Some("ethanol"));

        let m = parse_line("CCO").unwrap();
        assert_eq!(m.name(), None);
    }
}
