//! SMARTS 写出。
//!
//! # 判据是往返恒等,不是逐字节相同
//!
//! 同一个查询有很多写法(`[CH1]`、`[C&H1]`、`[H1&C]` 都一样),逐字节比对
//! 比的是写法而不是语义。这里守的是 **parse → write → parse 得到同一棵
//! 表达式树**。
//!
//! # 原子表达式里**没有括号**,所以运算符不能随便挑
//!
//! `&` 与 `;` 都对应 [`AtomExpr::And`],但优先级不同(`&` > `,` > `;`)。
//! 没有括号可用,选错就改语义:
//!
//! | 表达式树 | 只能写成 | 写成 `&` 会变成 |
//! |---|---|---|
//! | `And[Or[C,N], H1]` | `[C,N;H1]` | `[C,N&H1]` = `C , (N&H1)` |
//!
//! 规则:`And` 的孩子里出现 `Or` 就必须用 `;`,否则 `&` 够用。
//!
//! 解析器产出的树保证写得出来 —— `Or` 的孩子来自 `&` 那一档,不可能再含
//! `,` 或 `;`,所以不会出现"要在 `,` 里嵌 `;`"这种无解的形状。

use std::fmt::Write as _;

use omgkit_core::ChiralTag;

use super::{AtomExpr, AtomPrim, BondExpr, BondPrim, QueryMol, Reaction};

/// 把查询分子写成 SMARTS。
#[must_use]
pub fn write(q: &QueryMol) -> String {
    let mut out = String::new();
    write_into(&mut out, q);
    out
}

/// 把反应写成 SMARTS 的三段式。
///
/// 试剂段为空时仍写出两个 `>` —— `A>>B` 与 `A>>B` 的区别不在这里,
/// 而在于**段数固定是三段**,少写一个 `>` 就成了另一条反应。
#[must_use]
pub fn write_reaction(r: &Reaction) -> String {
    let join = |v: &[QueryMol]| v.iter().map(write).collect::<Vec<_>>().join(".");
    format!(
        "{}>{}>{}",
        join(&r.reactants),
        join(&r.agents),
        join(&r.products)
    )
}

fn write_into(out: &mut String, q: &QueryMol) {
    let mol = &q.topology;
    let n = mol.num_atoms();
    if n == 0 {
        return;
    }

    // 两趟。第一趟只定生成树:哪些边是树边、哪些是环闭合边。
    //
    // 不能一趟写完 —— 环闭合是从**后**一个原子发现的回边,而**前**一个原子
    // 那时早已写进串里,没法回头补标号。SMILES 写出器同样是两趟,原因一样。
    let mut visited = vec![false; n];
    let mut edge_used = vec![false; mol.num_bonds()];
    let mut children: Vec<Vec<(u32, u32)>> = vec![Vec::new(); n];
    let mut closures: Vec<Vec<u32>> = vec![Vec::new(); n];
    let mut roots: Vec<u32> = Vec::new();

    for root in 0..n as u32 {
        if visited[root as usize] {
            continue;
        }
        visited[root as usize] = true;
        roots.push(root);
        let mut stack = vec![root];
        while let Some(a) = stack.pop() {
            for (other, bond) in mol.neighbors(a) {
                if edge_used[bond as usize] {
                    continue;
                }
                edge_used[bond as usize] = true;
                if visited[other as usize] {
                    // 环闭合:两端都要记,先写到的那一端开环
                    closures[a as usize].push(bond);
                    closures[other as usize].push(bond);
                } else {
                    visited[other as usize] = true;
                    children[a as usize].push((other, bond));
                    stack.push(other);
                }
            }
        }
    }

    // 第二趟:按树写串
    let mut open: Vec<Option<u32>> = vec![None; mol.num_bonds()];
    let mut next_label = 1u32;
    for (i, &root) in roots.iter().enumerate() {
        if i > 0 {
            out.push('.');
        }
        emit(
            out,
            q,
            root,
            None,
            &children,
            &closures,
            &mut open,
            &mut next_label,
        );
    }
}

#[allow(clippy::too_many_arguments)]
fn emit(
    out: &mut String,
    q: &QueryMol,
    atom: u32,
    via: Option<u32>,
    children: &[Vec<(u32, u32)>],
    closures: &[Vec<u32>],
    open: &mut [Option<u32>],
    next_label: &mut u32,
) {
    if let Some(bond) = via {
        out.push_str(&bond_expr_string(&q.bonds[bond as usize]));
    }
    out.push('[');
    out.push_str(&atom_expr_string(&for_writing(
        q, atom, via, children, closures,
    )));
    out.push(']');

    for &bond in &closures[atom as usize] {
        match open[bond as usize] {
            Some(label) => {
                // 闭环端只写标号,键表达式已经写在开环端了
                out.push_str(&label_string(label));
                open[bond as usize] = None;
            }
            None => {
                let label = *next_label;
                *next_label += 1;
                open[bond as usize] = Some(label);
                out.push_str(&bond_expr_string(&q.bonds[bond as usize]));
                out.push_str(&label_string(label));
            }
        }
    }

    let kids = &children[atom as usize];
    for (i, &(other, bond)) in kids.iter().enumerate() {
        let last = i + 1 == kids.len();
        if !last {
            out.push('(');
        }
        emit(
            out,
            q,
            other,
            Some(bond),
            children,
            closures,
            open,
            next_label,
        );
        if !last {
            out.push(')');
        }
    }
}

/// 把手性标记从**存储序**换回这次写出的**书写序**。
///
/// 解析时做过一次相反的换算(见 `mol::Parser::fix_chirality`),写出必须与它
/// 互逆,否则 `解析 → 写出 → 解析` 会连翻两次,往返不再幂等 —— 而分子看上去
/// 完全正常,只是成了镜像。
///
/// 书写序由这次遍历决定:前驱键在最前,然后是本原子上的环闭合,最后是分支。
/// 括号氢那一项同理:原子是片段根(没有前驱)时氢排在最前,否则排在第二位。
fn for_writing(
    q: &QueryMol,
    atom: u32,
    via: Option<u32>,
    children: &[Vec<(u32, u32)>],
    closures: &[Vec<u32>],
) -> AtomExpr {
    let expr = &q.atoms[atom as usize];
    let Some(tag) = super::required_chirality(expr) else {
        return expr.clone();
    };
    if !tag.is_tetrahedral() {
        return expr.clone();
    }

    let mut written: Vec<u32> = Vec::new();
    if let Some(b) = via {
        written.push(b);
    }
    written.extend(closures[atom as usize].iter().copied());
    written.extend(children[atom as usize].iter().map(|&(_, b)| b));

    let stored: Vec<u32> = q
        .topology
        .bonds()
        .iter()
        .enumerate()
        .filter(|(_, b)| b.other_end(atom).is_some())
        .map(|(i, _)| u32::try_from(i).unwrap_or(u32::MAX))
        .collect();
    if written.len() != stored.len() {
        return expr.clone();
    }

    // 置换与其逆的奇偶相同,所以这里与解析侧用同一段计算
    let mut odd = super::mol::permutation_is_odd(&stored, &written);
    // 括号氢那一项也共用解析侧那份判定,两边必须互逆
    if stored.len() == 3
        && super::mol::needs_h_compensation(
            via.is_none(),
            expr,
            closures[atom as usize].len(),
            super::mol::has_unsaturated_bond(&q.topology, &q.bonds, atom),
        )
    {
        odd = !odd;
    }
    let mut out = expr.clone();
    if odd {
        super::mol::invert_chirality(&mut out);
    }
    out
}

fn label_string(label: u32) -> String {
    if label < 10 {
        label.to_string()
    } else {
        format!("%{label:02}")
    }
}

/// 原子表达式渲染。返回**方括号内部**的内容,不含括号本身。
#[must_use]
pub fn atom_expr_string(e: &AtomExpr) -> String {
    match e {
        AtomExpr::Prim(p) => prim_string(p),
        AtomExpr::Not(inner) => format!("!{}", atom_expr_string(inner)),
        AtomExpr::Or(parts) => parts
            .iter()
            .map(atom_expr_string)
            .collect::<Vec<_>>()
            .join(","),
        AtomExpr::And(parts) => {
            // 孩子里有 `Or` 就必须用 `;` —— `&` 比 `,` 紧,会把括号不存在的
            // 这门语言里唯一的分组手段弄丢
            let sep = if parts.iter().any(|p| matches!(p, AtomExpr::Or(_))) {
                ";"
            } else {
                "&"
            };
            parts
                .iter()
                .map(atom_expr_string)
                .collect::<Vec<_>>()
                .join(sep)
        }
    }
}

fn prim_string(p: &AtomPrim) -> String {
    match p {
        AtomPrim::Any => "*".into(),
        AtomPrim::Aromatic => "a".into(),
        AtomPrim::Aliphatic => "A".into(),
        AtomPrim::Element { z, aromatic } => match aromatic {
            // 不限芳香性只能写 `#n` —— 元素符号本身就带着芳香性的含义
            None => format!("#{z}"),
            Some(arom) => {
                let sym = omgkit_core::element::by_atomic_num(*z).map_or("*", |e| e.symbol);
                if *arom {
                    sym.to_ascii_lowercase()
                } else {
                    sym.to_string()
                }
            }
        },
        AtomPrim::Degree(n) => format!("D{n}"),
        AtomPrim::TotalDegree(n) => format!("X{n}"),
        AtomPrim::TotalHs(n) => format!("H{n}"),
        AtomPrim::ImplicitHs(n) => format!("h{n}"),
        AtomPrim::RingCount(n) => opt_num("R", *n),
        AtomPrim::RingSize(n) => opt_num("r", *n),
        AtomPrim::RingBondCount(n) => opt_num("x", *n),
        AtomPrim::Valence(n) => format!("v{n}"),
        AtomPrim::Charge(c) => charge_string(*c),
        AtomPrim::Isotope(m) => m.to_string(),
        AtomPrim::AtomMap(m) => format!(":{m}"),
        AtomPrim::Chirality(t) => chirality_string(*t),
        AtomPrim::Recursive(q) => format!("$({})", write(q)),
    }
}

fn opt_num(sym: &str, n: Option<u32>) -> String {
    match n {
        Some(v) => format!("{sym}{v}"),
        None => sym.to_string(),
    }
}

fn charge_string(c: i32) -> String {
    let mut s = String::new();
    let sign = if c < 0 { '-' } else { '+' };
    let mag = c.unsigned_abs();
    // `+` 与 `++` 都合法,但 `+2` 更不容易读错,统一写数字形式
    let _ = write!(s, "{sign}{mag}");
    s
}

fn chirality_string(t: ChiralTag) -> String {
    match t {
        ChiralTag::Ccw => "@".into(),
        ChiralTag::Cw => "@@".into(),
        // 非四面体的几何在 SMARTS 里没有对应写法,退回"有手性"
        _ => "@".into(),
    }
}

/// 键表达式渲染。
///
/// **默认键写空串**。SMARTS 里"没写键符号"的默认是析取式`单键或芳香键`
/// (两端可能是通配,当场定不下来),照着渲染就成了 `-,:` —— 语义没错,
/// 但每根普通键都多两个字符,而且再解析回来还是同一棵树。省掉更合适。
#[must_use]
pub fn bond_expr_string(e: &BondExpr) -> String {
    if *e == BondExpr::default_bond() {
        return String::new();
    }
    match e {
        BondExpr::Prim(p) => bond_prim_string(*p),
        BondExpr::Not(inner) => format!("!{}", bond_expr_string(inner)),
        BondExpr::Or(parts) => parts
            .iter()
            .map(bond_expr_string)
            .collect::<Vec<_>>()
            .join(","),
        BondExpr::And(parts) => {
            let sep = if parts.iter().any(|p| matches!(p, BondExpr::Or(_))) {
                ";"
            } else {
                "&"
            };
            parts
                .iter()
                .map(bond_expr_string)
                .collect::<Vec<_>>()
                .join(sep)
        }
    }
}

fn bond_prim_string(p: BondPrim) -> String {
    match p {
        BondPrim::Any => "~",
        BondPrim::Single => "-",
        BondPrim::Double => "=",
        BondPrim::Triple => "#",
        BondPrim::Quadruple => "$",
        BondPrim::Aromatic => ":",
        BondPrim::InRing => "@",
        BondPrim::UpRight => "/",
        BondPrim::DownRight => "\\",
        BondPrim::Dative => "->",
        BondPrim::DativeReversed => "<-",
    }
    .to_string()
}
