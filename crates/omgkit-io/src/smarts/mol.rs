//! 完整 SMARTS 模式的解析:拓扑 + 逐原子逐键的查询树。
//!
//! 拓扑部分(分支、环闭合、片段分隔)与 SMILES 完全同构 —— SMARTS 在语法上
//! 就是 SMILES 的超集,差别集中在原子和键的位置上可以写表达式。
//!
//! # 与 SMILES 解析共享哪些约定
//!
//! 环闭合键同样**延后追加到键表末尾**,端点朝向也照 SMILES 的规矩。这不是
//! 图省事:上层拿到 [`QueryMol`] 之后做匹配、做写出,都指望"查询分子"和
//! "被查分子"的拓扑约定一致。两边各搞一套的话,每个用到端点顺序的地方都要
//! 分情况。
//!
//! # 有机子集的裸原子含义不同
//!
//! SMILES 里 `C` 就是碳。SMARTS 里 `C` 是**查询**"脂肪碳",`c` 是"芳香碳",
//! `*` 是"任意原子"。方括号外只允许有机子集与 `*`,和 SMILES 一样。

use omgkit_core::{AtomData, AtomFlags, BondData, BondFlags, BondOrder, MolBuilder};

use super::bond::{parse_bond_expr, starts_bond_expr};
use super::expr::{AtomExpr, AtomPrim, BondExpr, BondPrim};
use super::parse::parse_atom_expr;
use super::QueryMol;
use crate::error::{ParseError, ParseErrorKind as K, Result};

/// 解析一条 SMARTS。
///
/// # Errors
/// 语法错误时返回带位置的 [`ParseError`];用 [`ParseError::render`] 可得到
/// 带插字号的两行视图。
pub fn parse(input: &str) -> Result<QueryMol> {
    Parser::new(input.as_bytes()).run()
}

/// 一个待配对的环闭合。
#[derive(Debug, Clone)]
struct RingOpen {
    atom: u32,
    expr: Option<BondExpr>,
    pos: usize,
    seq: u32,
}

/// 已配对、待追加到键表末尾的环键。
#[derive(Debug, Clone)]
struct RingBond {
    number: u32,
    seq: u32,
    begin: u32,
    end: u32,
    expr: BondExpr,
    open_seq: u32,
}

struct Parser<'a> {
    src: &'a [u8],
    pos: usize,
    topology: MolBuilder,
    atoms: Vec<AtomExpr>,
    bonds: Vec<BondExpr>,
    prev: Option<u32>,
    branches: Vec<(u32, usize)>,
    rings: std::collections::HashMap<u32, RingOpen>,
    ring_bonds: Vec<RingBond>,
    pending: Option<(BondExpr, usize)>,
    /// 每处环闭合的唯一 id,两端共用
    next_ring_seq: u32,
    /// 原子 → 在**该原子处**书写的环闭合 id,按书写先后
    ring_at_atom: std::collections::HashMap<u32, Vec<u32>>,
    /// 是不是某个 `.` 片段的首原子。括号氢的位置随它而变。
    fragment_start: Vec<bool>,
}

impl<'a> Parser<'a> {
    fn new(src: &'a [u8]) -> Self {
        Self {
            src,
            pos: 0,
            topology: MolBuilder::new(),
            atoms: Vec::new(),
            bonds: Vec::new(),
            prev: None,
            branches: Vec::new(),
            rings: std::collections::HashMap::new(),
            ring_bonds: Vec::new(),
            pending: None,
            next_ring_seq: 0,
            ring_at_atom: std::collections::HashMap::new(),
            fragment_start: Vec::new(),
        }
    }

    fn err<T>(&self, kind: K, pos: usize) -> Result<T> {
        Err(ParseError::new(kind, pos, self.src))
    }

    fn peek(&self) -> Option<u8> {
        self.src.get(self.pos).copied()
    }

    fn run(mut self) -> Result<QueryMol> {
        if self.src.is_empty() {
            return self.err(K::Empty, 0);
        }
        while let Some(b) = self.peek() {
            match b {
                b'(' => self.open_branch()?,
                b')' => self.close_branch()?,
                b'[' => self.bracket_atom()?,
                b'.' => {
                    if let Some((_, p)) = &self.pending {
                        return self.err(K::DanglingBond, *p);
                    }
                    self.pos += 1;
                    self.prev = None;
                }
                b'0'..=b'9' => {
                    let pos = self.pos;
                    let n = u32::from(b - b'0');
                    self.pos += 1;
                    self.ring_closure(n, pos)?;
                }
                b'%' => {
                    let pos = self.pos;
                    self.pos += 1;
                    let n = self.ring_number()?;
                    self.ring_closure(n, pos)?;
                }
                _ if starts_bond_expr(b) => self.bond_symbol()?,
                _ => self.organic_atom()?,
            }
        }

        if let Some((_, p)) = &self.pending {
            return self.err(K::DanglingBond, *p);
        }
        if !self.branches.is_empty() {
            let pos = self.pos;
            return self.err(K::UnbalancedParen, pos);
        }
        if let Some((&num, open)) = self.rings.iter().min_by_key(|(_, o)| o.pos) {
            let pos = open.pos;
            return self.err(K::UnclosedRingBond(num), pos);
        }
        if self.topology.num_atoms() == 0 {
            return self.err(K::Empty, 0);
        }

        // 环键统一追加到末尾,与 SMILES 解析同一套约定
        let mut ring_bonds = std::mem::take(&mut self.ring_bonds);
        ring_bonds.sort_by_key(|r| (r.number, r.seq));
        let (pos, src) = (self.pos, self.src);
        let mut ring_bond_index: std::collections::HashMap<u32, u32> =
            std::collections::HashMap::new();
        for rb in ring_bonds {
            let mut bd = BondData::new(rb.begin, rb.end, BondOrder::Unspecified);
            bd.flags.insert(BondFlags::HAS_QUERY);
            let idx = self
                .topology
                .add_bond_data(bd)
                .map_err(|_| ParseError::new(K::RingBondToSelf(rb.number), pos, src))?;
            ring_bond_index.insert(rb.open_seq, idx);
            self.bonds.push(rb.expr);
        }

        self.fix_chirality(&ring_bond_index);

        let q = QueryMol {
            topology: self.topology,
            atoms: self.atoms,
            bonds: self.bonds,
        };
        debug_assert!(q.is_consistent(), "查询树与拓扑长度不一致");
        Ok(q)
    }

    /// 把带定值四面体手性的标记从"书写序"换算到"存储序"。
    ///
    /// 两项相加,与 SMILES 解析里那套同构:
    ///
    /// **一、环闭合造成的置换。** 环键在串里写在紧跟原子之后,却统一追加到
    /// 键表末尾。按键下标排出来的存储序与书写序差一个置换,算它的奇偶。
    ///
    /// **二、括号氢的位置。** 氢不是图上的节点,它在四元组里占哪一位由书写
    /// 规则定:紧跟前一个原子。原子是片段首原子时没有前一个,氢落在**第一位**;
    /// 否则落在第二位。两者差一次对换,再翻一次。
    ///
    /// 少任何一项都不对:只补第一项会翻掉本来正确的写法,只补第二项则连环闭合
    /// 带来的置换都没消掉。
    ///
    /// # 第二项只在写了三根键时补
    ///
    /// 氢要占哪一位,前提是它**在四元组里**。查询只写了一两根键时,四个位置
    /// 凑不齐,标记在查询自身的范围内还定不下构型 —— 那时氢的位置无从谈起,
    /// 补一次反而凭空翻掉一次。
    ///
    /// 这个界限在匹配那一路看不出来:匹配器对写了不到三根键的查询原子只要求
    /// "底物这里有手性",不判是哪一个(见 `chirality_ok`)。看得出来的是产物
    /// 构建 —— 反应模板里 `[C@H:1]-[Cl:2]` 这种只写一根键的中心极常见,标记
    /// 会原样落到产物上,翻错就是造出对映体。
    fn fix_chirality(&mut self, ring_bond_index: &std::collections::HashMap<u32, u32>) {
        for atom in 0..self.topology.num_atoms() as u32 {
            let expr = &self.atoms[atom as usize];
            let Some(tag) = super::required_chirality(expr) else {
                continue;
            };
            if !tag.is_tetrahedral() {
                continue;
            }

            let ring_here: Vec<u32> = self
                .ring_at_atom
                .get(&atom)
                .map(|seqs| {
                    seqs.iter()
                        .filter_map(|s| ring_bond_index.get(s).copied())
                        .collect()
                })
                .unwrap_or_default();

            // 存储序:按键下标递增
            let stored: Vec<u32> = self
                .topology
                .bonds()
                .iter()
                .enumerate()
                .filter(|(_, b)| b.other_end(atom).is_some())
                .map(|(i, _)| u32::try_from(i).unwrap_or(u32::MAX))
                .collect();

            // 书写序:非环键按邻居原子序号排,环闭合整体插在"自身"的位置。
            // DFS 下前驱序号必小于本原子、后继必大于本原子,所以这样排得出串里的先后。
            let mut entries: Vec<(u32, Option<u32>)> = vec![(atom, None)];
            for (i, b) in self.topology.bonds().iter().enumerate() {
                let i = u32::try_from(i).unwrap_or(u32::MAX);
                if ring_here.contains(&i) {
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
                    None => written.extend(ring_here.iter().copied()),
                    Some(i) => written.push(i),
                }
            }
            if written.len() != stored.len() {
                continue; // 结构异常,已在别处报错
            }

            let mut odd = permutation_is_odd(&written, &stored);
            if stored.len() == 3
                && needs_h_compensation(
                    self.fragment_start[atom as usize],
                    expr,
                    ring_here.len(),
                    has_unsaturated_bond(&self.topology, &self.bonds, atom),
                )
            {
                odd = !odd;
            }
            if odd {
                invert_chirality(&mut self.atoms[atom as usize]);
            }
        }
    }

    // -- 分支 --------------------------------------------------------------

    fn open_branch(&mut self) -> Result<()> {
        let pos = self.pos;
        self.pos += 1;
        match self.prev {
            Some(a) => {
                self.branches.push((a, self.topology.num_atoms()));
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
                if self.topology.num_atoms() == n_at_open {
                    return self.err(K::EmptyBranch, pos);
                }
                if let Some((_, p)) = &self.pending {
                    return self.err(K::DanglingBond, *p);
                }
                self.prev = Some(atom);
                Ok(())
            }
            None => self.err(K::UnbalancedParen, pos),
        }
    }

    // -- 键 ----------------------------------------------------------------

    /// 读一段键表达式。表达式的边界是"下一个不属于键表达式的字节"。
    fn bond_symbol(&mut self) -> Result<()> {
        let pos = self.pos;
        if self.pending.is_some() {
            return self.err(K::DanglingBond, pos);
        }
        let start = self.pos;
        let mut prev = 0u8;
        while let Some(b) = self.peek() {
            // `>` 只在紧跟 `-` 时才属于键表达式(配位键 `->`)。
            //
            // 不能把 `>` 直接并入"键表达式字符"的集合:反应 SMARTS 里
            // `>` 是组分分隔符(`反应物>试剂>产物`),吞掉它会把整条反应
            // 读成一个分子。
            let part_of_expr = starts_bond_expr(b)
                || matches!(b, b'&' | b',' | b';')
                || (b == b'>' && prev == b'-');
            if !part_of_expr {
                break;
            }
            prev = b;
            self.pos += 1;
        }
        let expr = parse_bond_expr(&self.src[start..self.pos])
            .map_err(|e| ParseError::new(e.kind, start + e.pos, self.src))?;
        self.pending = Some((expr, pos));
        Ok(())
    }

    // -- 原子 --------------------------------------------------------------

    /// 方括号外的原子:有机子集与 `*`,含义是查询而非具体元素。
    fn organic_atom(&mut self) -> Result<()> {
        let pos = self.pos;
        let b = self.peek().expect("已 peek");
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
            (b'a', _) => {
                self.pos += 1;
                return self.push_atom(AtomExpr::Prim(AtomPrim::Aromatic), pos);
            }
            (b'A', _) => {
                self.pos += 1;
                return self.push_atom(AtomExpr::Prim(AtomPrim::Aliphatic), pos);
            }
            (b'*', _) => {
                self.pos += 1;
                return self.push_atom(AtomExpr::Prim(AtomPrim::Any), pos);
            }
            _ => {
                return self.err(K::UnknownElement(char::from(b).to_string()), pos);
            }
        };
        self.pos += len;
        self.push_atom(
            AtomExpr::Prim(AtomPrim::Element {
                z,
                aromatic: Some(aromatic),
            }),
            pos,
        )
    }

    /// `[` ... `]`。内容整体交给表达式解析器,但要先过一遍 `[H]` 特例表。
    fn bracket_atom(&mut self) -> Result<()> {
        let open_pos = self.pos;
        self.pos += 1;
        let start = self.pos;

        // 找到配对的 `]`。递归 SMARTS `$(...)` 里可能嵌套方括号,要计数。
        let mut depth = 0usize;
        loop {
            match self.peek() {
                None => return self.err(K::UnexpectedEnd, self.pos),
                Some(b'[') => {
                    depth += 1;
                    self.pos += 1;
                }
                Some(b']') if depth == 0 => break,
                Some(b']') => {
                    depth -= 1;
                    self.pos += 1;
                }
                Some(_) => self.pos += 1,
            }
        }
        debug_assert_eq!(depth, 0, "方括号计数没归零");
        let inner = &self.src[start..self.pos];
        self.pos += 1; // ']'

        let expr = match hydrogen_special_case(inner) {
            Some(e) => e,
            None => parse_atom_expr(inner)
                .map_err(|e| ParseError::new(e.kind, start + e.pos, self.src))?,
        };
        self.push_atom(expr, open_pos)
    }

    fn push_atom(&mut self, expr: AtomExpr, pos: usize) -> Result<()> {
        let mut atom = AtomData::new(0);
        atom.flags.insert(AtomFlags::HAS_QUERY);
        let idx = self.topology.add_atom_data(atom);
        self.fragment_start.push(self.prev.is_none());
        self.atoms.push(expr);

        if let Some(prev) = self.prev {
            let expr = self
                .pending
                .take()
                .map_or_else(BondExpr::default_bond, |(e, _)| e);
            let mut bd = BondData::new(prev, idx, BondOrder::Unspecified);
            bd.flags.insert(BondFlags::HAS_QUERY);
            let src = self.src;
            self.topology
                .add_bond_data(bd)
                .map_err(|_| ParseError::new(K::UnexpectedChar('?'), pos, src))?;
            self.bonds.push(expr);
        } else if let Some((_, p)) = self.pending.take() {
            return self.err(K::DanglingBond, p);
        }
        self.prev = Some(idx);
        Ok(())
    }

    // -- 环闭合 ------------------------------------------------------------

    fn ring_number(&mut self) -> Result<u32> {
        let pos = self.pos;
        let mut v = 0u32;
        let mut n = 0;
        while n < 2 {
            match self.peek() {
                Some(d @ b'0'..=b'9') => {
                    v = v * 10 + u32::from(d - b'0');
                    self.pos += 1;
                    n += 1;
                }
                _ => break,
            }
        }
        if n == 2 {
            Ok(v)
        } else {
            self.err(K::BadBracketAtom("`%` 后需要两位数字"), pos)
        }
    }

    fn ring_closure(&mut self, num: u32, pos: usize) -> Result<()> {
        let Some(cur) = self.prev else {
            return self.err(K::UnexpectedChar(char::from(self.src[pos])), pos);
        };
        let pending = self.pending.take().map(|(e, _)| e);

        match self.rings.remove(&num) {
            Some(open) => {
                if open.atom == cur {
                    return self.err(K::RingBondToSelf(num), pos);
                }
                if open.expr.is_some() && pending.is_some() && open.expr != pending {
                    return self.err(K::ConflictingRingBondOrder(num), pos);
                }
                // 键表达式可以写在**任一端**,而有朝向的基元(`->` `<-` `/` `\`)
                // 说的是"从我这一端看过去"。所以端点顺序要跟着表达式的出处走:
                // 写在开环端就存成 开环→闭环,写在闭环端就存成 闭环→开环。
                //
                // 固定成一种顺序的话,`[O]->1...[Fe]1` 与 `[O]1...[Fe]<-1` 这两种
                // 写同一件事的形式会存出相反的朝向,而下游(产物构建、匹配)
                // 只看得到端点,分不出是哪种写法留下的。
                let (expr, from_open) = match (open.expr, pending) {
                    (Some(e), _) => (e, true),
                    (None, Some(e)) => (e, false),
                    (None, None) => (BondExpr::default_bond(), false),
                };
                let (begin, end) = if from_open {
                    (open.atom, cur)
                } else {
                    (cur, open.atom)
                };
                let seq = u32::try_from(self.ring_bonds.len()).unwrap_or(u32::MAX);
                self.ring_at_atom.entry(cur).or_default().push(open.seq);
                self.ring_bonds.push(RingBond {
                    number: num,
                    seq,
                    begin,
                    end,
                    expr,
                    open_seq: open.seq,
                });
            }
            None => {
                let seq = self.next_ring_seq;
                self.next_ring_seq += 1;
                self.ring_at_atom.entry(cur).or_default().push(seq);
                self.rings.insert(
                    num,
                    RingOpen {
                        atom: cur,
                        expr: pending,
                        pos,
                        seq,
                    },
                );
            }
        }
        Ok(())
    }
}

/// `[H]` 特例表:可选同位素 + `H` + 可选电荷 + 可选映射号,整个括号到此为止。
///
/// 匹配上的才是**氢元素**;其余一律走一般规则,那里 `H` 是氢计数。
/// 详见 [`mod@super::parse`] 的模块文档。
///
/// 链接要写 `mod@`:`smarts` 底下 `parse` 既是一个私有模块、又是 `mol::parse`
/// 这个函数的再导出名,不消歧的话 rustdoc 报"既是函数又是模块"。
fn hydrogen_special_case(inner: &[u8]) -> Option<AtomExpr> {
    let mut i = 0;
    let mut parts: Vec<AtomExpr> = Vec::new();

    // 可选同位素
    let start = i;
    while i < inner.len() && inner[i].is_ascii_digit() {
        i += 1;
    }
    if i > start {
        let n: u32 = std::str::from_utf8(&inner[start..i]).ok()?.parse().ok()?;
        parts.push(AtomExpr::Prim(AtomPrim::Isotope(u16::try_from(n).ok()?)));
    }

    // 必须是 H
    if inner.get(i) != Some(&b'H') {
        return None;
    }
    i += 1;
    parts.push(AtomExpr::Prim(AtomPrim::Element {
        z: 1,
        aromatic: None,
    }));

    // 可选电荷
    if let Some(&c @ (b'+' | b'-')) = inner.get(i) {
        i += 1;
        let sign: i32 = if c == b'+' { 1 } else { -1 };
        let mut n = 1i32;
        while inner.get(i) == Some(&c) {
            i += 1;
            n += 1;
        }
        if n == 1 {
            let s = i;
            while i < inner.len() && inner[i].is_ascii_digit() {
                i += 1;
            }
            if i > s {
                n = std::str::from_utf8(&inner[s..i]).ok()?.parse().ok()?;
            }
        }
        parts.push(AtomExpr::Prim(AtomPrim::Charge(n * sign)));
    }

    // 可选映射号
    if inner.get(i) == Some(&b':') {
        i += 1;
        let s = i;
        while i < inner.len() && inner[i].is_ascii_digit() {
            i += 1;
        }
        if i == s {
            return None;
        }
        let n: u32 = std::str::from_utf8(&inner[s..i]).ok()?.parse().ok()?;
        parts.push(AtomExpr::Prim(AtomPrim::AtomMap(u16::try_from(n).ok()?)));
    }

    // 必须正好读完 —— 有残余就说明不在表里
    if i != inner.len() {
        return None;
    }
    Some(if parts.len() == 1 {
        parts.pop().expect("非空")
    } else {
        AtomExpr::And(parts)
    })
}

/// `from` → `to` 的置换是否为奇。两者不是同一多重集时按偶处理(不翻)。
pub(super) fn permutation_is_odd(from: &[u32], to: &[u32]) -> bool {
    let mut cur = from.to_vec();
    let mut swaps = 0usize;
    for i in 0..to.len() {
        if cur[i] == to[i] {
            continue;
        }
        let Some(j) = (i + 1..cur.len()).find(|&j| cur[j] == to[i]) else {
            return false;
        };
        cur.swap(i, j);
        swaps += 1;
    }
    swaps % 2 == 1
}

/// 括号氢那一项该不该补。解析与写出**共用这一份**。
///
/// 两处必须互逆:解析把标记从书写序换到存储序,写出换回来。各写一份的话,
/// 改动只落到一边就会让往返连翻两次 —— 而分子看上去毫无异常,只是成了镜像。
///
/// 两种情形都是"氢占了四元组里的一位,而书写序把它放在了别处":
///
/// - **片段首原子 + 一个括号氢** —— 氢本该紧跟前一个原子,首原子没有前一个,
///   于是氢落到第一位而不是第二位,差一次对换。
/// - **没有括号氢 + 恰好一个环闭合 + 不含不饱和键** —— 第四位由环闭合占,
///   而环闭合键被追加到了键表末尾,这一支在只算键的置换里补不回来。
///
/// 调用方还要自己判"这个原子写了三根键"。写得更少时四个位置凑不齐,标记在
/// 查询自身的范围内定不下构型,氢的位置无从谈起。
pub(super) fn needs_h_compensation(
    is_root: bool,
    expr: &AtomExpr,
    n_ring: usize,
    unsaturated: bool,
) -> bool {
    let h = bracket_h_count(expr);
    if is_root && h == Some(1) {
        return true;
    }
    h != Some(1) && n_ring == 1 && !unsaturated
}

/// 这个原子上有没有键**要求**键级大于一。
///
/// 析取与否定下面的键级说不定,一律当作单键 —— 与产物构建里取键级的走法一致。
pub(super) fn has_unsaturated_bond(topology: &MolBuilder, bonds: &[BondExpr], atom: u32) -> bool {
    fn multiple(e: &BondExpr) -> bool {
        match e {
            BondExpr::Prim(p) => matches!(
                p,
                BondPrim::Double | BondPrim::Triple | BondPrim::Quadruple | BondPrim::Aromatic
            ),
            BondExpr::And(parts) => parts.iter().any(multiple),
            _ => false,
        }
    }
    topology
        .bonds()
        .iter()
        .enumerate()
        .filter(|(_, b)| b.other_end(atom).is_some())
        .any(|(i, _)| bonds.get(i).is_some_and(multiple))
}

/// 方括号里写的氢数。没写返回 `None`。
pub(super) fn bracket_h_count(expr: &AtomExpr) -> Option<u32> {
    match expr {
        AtomExpr::Prim(AtomPrim::TotalHs(n)) => Some(*n),
        AtomExpr::And(parts) => parts.iter().find_map(bracket_h_count),
        _ => None,
    }
}

/// 把表达式里的四面体手性基元翻个面。
///
/// 遍历范围与 [`required_chirality`](super::required_chirality) 一致 ——
/// 只走合取,不进析取与否定。那两处的手性没有定值。
pub(super) fn invert_chirality(expr: &mut AtomExpr) {
    match expr {
        AtomExpr::Prim(AtomPrim::Chirality(t)) => *t = t.inverted(),
        AtomExpr::And(parts) => {
            for p in parts {
                invert_chirality(p);
            }
        }
        _ => {}
    }
}
