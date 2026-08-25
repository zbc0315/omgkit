//! SMILES 写出。
//!
//! # 输出顺序由调用方决定
//!
//! 写出本身不挑顺序 —— 给定一个优先级数组,它就照着走 DFS。这样"怎么排"
//! (规范化排序)和"怎么写"(本模块)是两件可以分别验证的事:
//!
//! - 写出的判据是**往返恒等**:解析 → 写出 → 再解析,得到同一个分子。
//!   这条判据不需要任何外部参照。
//! - 排序的判据是**重排不变**:原子编号任意重排,规范秩不变。
//!
//! 把两者揉在一起的话,一个失败就分不清是谁的锅。
//!
//! # 写法的取舍由调用方选
//!
//! 两种取舍服务于两条互相冲突的性质,见 [`WriteStyle`]:
//!
//! - **往返恒等**要照原样再现方括号 —— `[CH3][CH2][OH]` 与 `CCO` 是同一个分子,
//!   却不是同一份原子表示,再解析要拿回原来那一份
//! - **规范**要抹掉这个差别 —— 同一个分子只能有一串
//!
//! "能省则省"用不上完整的价键模型:判据只要键级和、总氢数与该元素的首位默认价,
//! 三样都在 L0。判据保守,算不准的一律留框 —— 留框只是啰嗦,去错了会改掉分子。
//!
//! # 立体化学
//!
//! 四面体手性(`@` / `@@`)会写出,并且做过与解析器互逆的宇称换算 ——
//! 输出会重排邻居,标记必须跟着换参照系,否则写出的是镜像分子,而拓扑
//! 完全正确、看不出来。换算式见 [`output_chiral_tag`]。
//!
//! 双键方向键(`/` `\`)也会写出。存储的方向一律相对键的 `begin → end`,
//! 而 DFS 从哪一端进入这条键不受存储顺序约束,所以写出时要按遍历方向换算
//! (见 [`bond_symbol`])。少了这次换算,顺式会写成反式。
//!
//! 配位几何(`@SP`/`@TB`/`@OH`)会写出。序号是相对"配体按什么顺序列出"的,
//! 所以要按该多面体的转动群从存储序换算到本次输出的顺序 ——
//! 见 [`omgkit_core::polyhedron::renumber`]。方括号里的氢、或者一个空的配位
//! 位置,也占一个顶点而不在键序列里;它落在"自身位置",两侧的配体序列由
//! `smiles::coordination_ligands` 统一补出来。缺两个及以上顶点时那几个顶点
//! 彼此分不开,序号换算不唯一,整个丢掉那个标记(丢掉是老实的,瞎写一个
//! 序号是撒谎)。
//!
//! 丙二烯的轴手性(`@AL`)**尚未写出**:它的立体信息属于一根轴而不是一个中心,
//! 没有多面体排列表可言。
//!
//! # 环闭合标号
//!
//! 分配的是**当前空闲的最小标号**,闭合后立即回收。标号超过 9 时写成 `%NN`,
//! 超过 99 时写成 `%(NNN)`。后者在同时打开的环超过 99 个时才会出现 ——
//! 稠合体系确实做得到。

use std::collections::BTreeMap;
use std::fmt::Write as _;

use omgkit_core::{element, AtomFlags, BondData, BondDirection, BondOrder, ChiralTag, MolBuilder};

/// 写出的结果。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Written {
    /// SMILES 字符串
    pub smiles: String,
    /// 输出中第 `i` 个原子对应原分子的哪个原子下标。
    ///
    /// 往返比对要用它:再解析得到的分子里,原子 `i` 对应原分子的
    /// `atom_order[i]`。没有它就只能比"分子是否同构",那要跑图同构,
    /// 既慢又会把写出的错误和匹配的错误混在一起。
    pub atom_order: Vec<u32>,
}

/// 写出的取舍:忠实回写,还是规范。
///
/// 两者服务于两条互相冲突的性质,所以必须由调用方选:
///
/// | | 要的性质 | 判据 |
/// |---|---|---|
/// | [`Faithful`](Self::Faithful) | **往返恒等** —— 再解析得到的原子逐字段相同 | `roundtrip_smiles.rs` |
/// | [`Canonical`](Self::Canonical) | **规范** —— 同一个分子只有一串 | `canonical_invariance.rs` |
///
/// 冲突是实打实的:`[CH3][CH2][OH]` 与 `CCO` 是同一个**分子**,却不是同一份
/// **原子表示**(前者 `NO_IMPLICIT` 置位、氢记在显式一侧)。往返要保住这个差别,
/// 规范化要抹掉它。
///
/// # 规范式要抹掉的"输入写法痕迹"有两处
///
/// 1. **方括号**:常见原子多写了框,去掉之后读回来氢数不变的就去掉
/// 2. **方向键的整体翻转**:同一约束片段内 `/` 与 `\` 可以全体互换而不改变
///    任何一对取代基的相对位置。这个自由度原先由**键的存储下标**定,而存储
///    下标是输入写法留下的痕迹;规范写法改成由**输出顺序**定 —— 每个片段第一个
///    写出来的方向符号一律取 `/`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WriteStyle {
    /// 照原样再现:方括号、方向键的写法都沿用分子里存着的那一份。
    Faithful,
    /// 抹掉输入写法的痕迹,只留分子本身决定的东西。
    ///
    /// 方括号**只会去、不会加** —— 本来没框的原子不受影响。加框需要知道简写形式
    /// 会推出几个氢,而未净化的分子上那个数还没算出来,猜了就会把分子写坏。
    Canonical,
}

/// 按原子存储顺序写出 SMILES。
///
/// 四面体手性与双键方向键会写出;配位几何尚未 —— 见 [`write_with_priority`]。
#[must_use]
pub fn write(mol: &MolBuilder) -> Written {
    let priority: Vec<u32> = (0..mol.num_atoms() as u32).collect();
    write_with_priority(mol, &priority)
}

/// 按给定优先级写出 SMILES。`priority[a]` 越小,原子 `a` 越早被访问。
///
/// 优先级同时决定两件事:每个连通片段的起点(片段内优先级最小者),以及
/// 每个原子处分支的先后。
///
/// # 立体化学
///
/// 四面体手性与双键方向键(`/` `\`)会写出,两者都做过与解析器互逆的换算。
///
/// 配位几何(`@SP`/`@TB`/`@OH`)目前**不输出**:那要一张排列换算表,属于 L6。
/// 在那之前**宁可不写**,也好过写出一个可能是错的立体信息 —— 立体写错了
/// 拓扑还是对的,只有分子是镜像的,极难发现。
///
/// # Panics
/// `priority` 长度与原子数不符时 panic —— 这是调用方的编程错误。
#[must_use]
pub fn write_with_priority(mol: &MolBuilder, priority: &[u32]) -> Written {
    write_with_priority_styled(mol, priority, WriteStyle::Faithful)
}

/// 同 [`write_with_priority`],但由调用方选写法,见 [`WriteStyle`]。
///
/// # Panics
/// `priority` 长度与原子数不符时 panic —— 这是调用方的编程错误。
#[must_use]
pub fn write_with_priority_styled(
    mol: &MolBuilder,
    priority: &[u32],
    style: WriteStyle,
) -> Written {
    let n = mol.num_atoms();
    assert_eq!(
        priority.len(),
        n,
        "优先级数组长度 {} 与原子数 {n} 不符",
        priority.len()
    );
    if n == 0 {
        return Written {
            smiles: String::new(),
            atom_order: Vec::new(),
        };
    }

    // 规范写法下,双键的参照原子改按规范秩挑,而不是沿用感知留下的那一个 ——
    // 感知挑的是"存储顺序里第一个带方向的邻居",那带着输入写法的痕迹,会让
    // 规范串不成不动点。理由与做法见 `stereo::normalized_stereo_refs`。
    //
    // 忠实写法不动:那边要的是原样再现输入,包括方向符号落在哪根键上。
    let normalized = if style == WriteStyle::Canonical {
        crate::stereo::normalized_stereo_refs(mol, priority)
    } else {
        None
    };
    let mol = normalized.as_ref().unwrap_or(mol);

    let tree = build_tree(mol, priority);
    emit(mol, &tree, style)
}

// ---------------------------------------------------------------------------
// 第一趟:DFS 生成树
// ---------------------------------------------------------------------------

/// DFS 的产物:哪些键是树边、哪些是环闭合边,以及每个原子的孩子顺序。
struct Dfs {
    /// 每个片段的根原子,按片段被访问的先后
    roots: Vec<u32>,
    /// 每个原子的孩子(存**键**下标),按访问先后
    children: Vec<Vec<u32>>,
    /// 每个原子处的环闭合(存**键**下标),按发现先后
    ring_closures: Vec<Vec<u32>>,
}

fn build_tree(mol: &MolBuilder, priority: &[u32]) -> Dfs {
    let n = mol.num_atoms();

    // 每个原子的邻居按优先级排好。分支先后与起点选择都只看优先级,
    // 与邻居的存储顺序无关 —— 存储顺序是建图留下的痕迹,不是分子的性质。
    let mut nbrs: Vec<Vec<(u32, u32)>> = Vec::with_capacity(n);
    for a in 0..n as u32 {
        let mut v: Vec<(u32, u32)> = mol
            .neighbors(a)
            .map(|(other, bond)| (bond, other))
            .collect();
        v.sort_unstable_by_key(|&(bond, other)| (priority[other as usize], bond));
        nbrs.push(v);
    }

    let mut order: Vec<u32> = (0..n as u32).collect();
    order.sort_unstable_by_key(|&a| priority[a as usize]);

    let mut visited = vec![false; n];
    let mut edge_used = vec![false; mol.num_bonds()];
    let mut dfs = Dfs {
        roots: Vec::new(),
        children: vec![Vec::new(); n],
        ring_closures: vec![Vec::new(); n],
    };

    // 显式栈而不是递归:大环语料里有几千个原子的分子,递归深度等于原子数。
    let mut stack: Vec<(u32, usize)> = Vec::new();
    for &root in &order {
        if visited[root as usize] {
            continue;
        }
        visited[root as usize] = true;
        dfs.roots.push(root);
        stack.push((root, 0));

        while let Some(&mut (a, ref mut cursor)) = stack.last_mut() {
            let Some(&(bond, other)) = nbrs[a as usize].get(*cursor) else {
                stack.pop();
                continue;
            };
            *cursor += 1;

            if edge_used[bond as usize] {
                continue; // 父边,或已记过的环闭合边
            }
            edge_used[bond as usize] = true;

            if visited[other as usize] {
                // 环闭合:两端都要记,先被写出的那一端开环
                dfs.ring_closures[a as usize].push(bond);
                dfs.ring_closures[other as usize].push(bond);
            } else {
                visited[other as usize] = true;
                dfs.children[a as usize].push(bond);
                stack.push((other, 0));
            }
        }
    }

    dfs
}

// ---------------------------------------------------------------------------
// 第二趟:按生成树写字符串
// ---------------------------------------------------------------------------

/// 写出期间要发生的事。用显式栈跑,理由同 [`build_tree`]。
enum Step {
    /// 写一个原子;`via` 是从父原子过来的那条键
    Atom { atom: u32, via: Option<u32> },
    /// 写一个字面量(分支括号、片段分隔符)
    Literal(&'static str),
}

fn emit(mol: &MolBuilder, dfs: &Dfs, style: WriteStyle) -> Written {
    // 每根键该写什么方向。感知过顺反的双键由它重新生成方向,没感知过的
    // 沿用存储的写法 —— 见 stereo::directions_for_writing
    let written = crate::stereo::directions_for_writing(mol);
    let (dirs, comps) = (written.dirs, written.component);
    // 片段 → 是否整体翻转。规范写法下由**第一个写出来的**方向符号定死,见 WriteStyle。
    let mut gauge: BTreeMap<u32, bool> = BTreeMap::new();
    let mut out = String::new();
    let mut atom_order = Vec::with_capacity(mol.num_atoms());
    // 键下标 → 已分配的环闭合标号。有值即表示该环已开、等着闭合。
    let mut open_label: Vec<Option<u32>> = vec![None; mol.num_bonds()];
    let mut label_in_use: Vec<bool> = Vec::new();

    let mut stack: Vec<Step> = Vec::new();
    for (i, &root) in dfs.roots.iter().enumerate() {
        if i > 0 {
            stack.push(Step::Literal("."));
        }
        stack.push(Step::Atom {
            atom: root,
            via: None,
        });
    }
    stack.reverse();

    while let Some(step) = stack.pop() {
        let (atom, via) = match step {
            Step::Literal(s) => {
                out.push_str(s);
                continue;
            }
            Step::Atom { atom, via } => (atom, via),
        };

        if let Some(bond) = via {
            // 箭头方向要从**父**原子看过去
            let parent = other_end(mol, bond, atom);
            out.push_str(gauged_symbol(
                bond, parent, mol, &dirs, &comps, style, &mut gauge,
            ));
        }

        // 该原子的邻居在输出串里出现的顺序:父键、环闭合键、子键。
        // 立体标记要相对这个顺序写,所以必须在写原子**之前**就定下来。
        let mut written_bonds: Vec<u32> = Vec::with_capacity(mol.degree(atom));
        written_bonds.extend(via);
        written_bonds.extend(dfs.ring_closures[atom as usize].iter().copied());
        written_bonds.extend(dfs.children[atom as usize].iter().copied());

        let tag = output_chiral_tag(
            mol,
            atom,
            &written_bonds,
            via.is_none(),
            dfs.ring_closures[atom as usize].len(),
        );
        write_atom(&mut out, mol, atom, tag, style);
        atom_order.push(atom);

        for &bond in &dfs.ring_closures[atom as usize] {
            match open_label[bond as usize] {
                // 第二次遇到:闭合并回收标号
                Some(label) => {
                    out.push_str(&ring_label(label));
                    open_label[bond as usize] = None;
                    // 标号从 1 起,槽位从 0 起
                    label_in_use[label as usize - 1] = false;
                }
                // 第一次遇到:开环。键级符号写在开环端 —— 与解析器的端点
                // 约定配套,配位键的箭头方向才能原样还原。
                None => {
                    let label = alloc_label(&mut label_in_use);
                    open_label[bond as usize] = Some(label);
                    out.push_str(gauged_symbol(
                        bond, atom, mol, &dirs, &comps, style, &mut gauge,
                    ));
                    out.push_str(&ring_label(label));
                }
            }
        }

        // 最后一个孩子不套括号
        let kids = &dfs.children[atom as usize];
        if let Some((&last, rest)) = kids.split_last() {
            stack.push(Step::Atom {
                atom: other_end(mol, last, atom),
                via: Some(last),
            });
            for &bond in rest.iter().rev() {
                stack.push(Step::Literal(")"));
                stack.push(Step::Atom {
                    atom: other_end(mol, bond, atom),
                    via: Some(bond),
                });
                stack.push(Step::Literal("("));
            }
        }
    }

    Written {
        smiles: out,
        atom_order,
    }
}

/// 把存储序上的四面体标记换算成**要写进串里**的标记。
///
/// # 换算式
///
/// 解析时做的是(见 `smiles` 模块文档约定二):
///
/// ```text
/// 存储标记 = 翻转^(p ⊕ c)(串里的标记)
/// p = 置换宇称(串里的邻居顺序 → 存储顺序)
/// c = 那条 degree==3 的补偿规则
/// ```
///
/// 写出要反过来解出"串里的标记"。翻转是对合的,宇称是可加的,于是
/// **同一个式子倒着用**即可:
///
/// ```text
/// 串里的标记 = 翻转^(p' ⊕ c')(存储标记)
/// p' = 置换宇称(本次输出的邻居顺序 → 本分子的存储顺序)
/// c' = 补偿规则,按**本次输出**的形态算
/// ```
///
/// 关键在于 `p'` 用的是**当前分子**的存储序,而不是重新解析之后的存储序 ——
/// 后者未知,但两处的差值恰好与 `p'` 相消。
///
/// # 只处理四面体
///
/// 配位几何(`@SP`/`@TB`/`@OH`)的排列序号换参照系要走查找表,属于 L6;
/// 在那之前它们不写出。丙二烯的轴手性同理。
fn output_chiral_tag(
    mol: &MolBuilder,
    atom: u32,
    written_bonds: &[u32],
    is_fragment_start: bool,
    ring_closures: usize,
) -> (ChiralTag, u8) {
    let a = mol.atoms()[atom as usize];

    // 配位几何(`@SP`/`@TB`/`@OH`):把序号从存储序换算到本次输出的顺序。
    //
    // 方括号里的氢、或者一个空的配位位置,也占一个顶点而不在键序列里:输出串里
    // 它落在"自身位置"(紧跟前驱原子之后、环闭合之前 —— 方括号本来就写在那儿),
    // 存储序里排最前,见 [`super::coordination_ligands`]。
    //
    // 换算不了(缺两个及以上顶点、序号越界)时写出**不带序号的类别**是不行的
    // —— `[Pt@SP]` 读回来是错的。所以给 0,由 `write_atom` 整个略过。
    // 解析侧用的是同一个条件(见 `fix_chirality`),两侧因此对齐。
    if omgkit_core::polyhedron::ligand_count(a.chiral_tag).is_some() {
        let perm = super::coordination_ligands(mol, atom, a.chiral_tag)
            .and_then(|stored| {
                let mut written = written_bonds.to_vec();
                if stored.len() == written.len() + 1 {
                    written.insert(usize::from(!is_fragment_start), super::VACANT_LIGAND);
                }
                omgkit_core::polyhedron::renumber(a.chiral_tag, a.stereo_perm, &stored, &written)
            })
            .unwrap_or(0);
        return (a.chiral_tag, perm);
    }

    if !a.chiral_tag.is_tetrahedral() {
        return (ChiralTag::Unspecified, 0);
    }

    let stored: Vec<u32> = mol.neighbors(atom).map(|(_, bond)| bond).collect();
    let Some(mut odd) = super::permutation_is_odd(written_bonds, &stored) else {
        // 两个序列不是同一个多重集 —— 只可能是本模块自己算错了邻居
        debug_assert!(false, "输出的邻居顺序与存储顺序不是同一组键");
        return (ChiralTag::Unspecified, 0);
    };

    // 补偿规则:隐式/显式氢不参与置换,由这条特判统一处理。
    // 判据要按**输出串**的形态算 —— 是不是片段首原子、有几个环闭合,
    // 都随写出方式而变。
    if stored.len() == 3 {
        let hs = total_hs(&a);
        let unsaturated = mol
            .neighbors(atom)
            .any(|(_, bond)| mol.bonds()[bond as usize].order.as_double() > 1.0);
        if (is_fragment_start && hs == 1) || (hs != 1 && ring_closures == 1 && !unsaturated) {
            odd = !odd;
        }
    }

    let tag = if odd {
        a.chiral_tag.inverted()
    } else {
        a.chiral_tag
    };
    (tag, 0)
}

/// 方括号里要写的氢数。
///
/// 未置 [`AtomFlags::NO_IMPLICIT`] 的原子把氢记在 `num_implicit_hs` 里,
/// 一旦要给它加方括号,氢数就必须显式写出来 —— 否则 `[C]` 会被读成零个氢。
/// 两个字段互斥(置位的那类隐式氢恒为 0),相加即总数。
fn total_hs(a: &omgkit_core::AtomData) -> u8 {
    a.num_explicit_hs.saturating_add(a.num_implicit_hs)
}

/// 取键 `bond` 上 `from` 的对端。`from` 必是端点之一,否则是本模块的逻辑错误。
fn other_end(mol: &MolBuilder, bond: u32, from: u32) -> u32 {
    mol.bonds()[bond as usize]
        .other_end(from)
        .expect("遍历产生的键必以当前原子为端点")
}

/// 取当前空闲的最小标号(从 1 起)。
fn alloc_label(in_use: &mut Vec<bool>) -> u32 {
    match in_use.iter().position(|&used| !used) {
        Some(i) => {
            in_use[i] = true;
            i as u32 + 1
        }
        None => {
            in_use.push(true);
            in_use.len() as u32
        }
    }
}

/// 环闭合标号的字面形式。
fn ring_label(label: u32) -> String {
    if label < 10 {
        label.to_string()
    } else if label < 100 {
        format!("%{label}")
    } else {
        // 同时打开的环超过 99 个 —— 稠合体系做得到
        format!("%({label})")
    }
}

/// 从 `from` 端看这条键该写什么符号。
///
/// # `from` 为什么是必需的
///
/// 有两类键的符号取决于**从哪一端看**:
///
/// - 配位键的箭头要指向受体(`->` / `<-`)
/// - 方向键 `/` 与 `\` 表达的是"从这一端走向另一端时是上行还是下行"
///
/// 存储里的 `direction` 一律相对 `begin → end`。写出时的遍历方向可能相反
/// (DFS 从哪个原子进入这条键不受存储顺序约束),那时必须翻转 —— 与解析器
/// 处理环闭合端点交换时是同一个变换。
///
/// 少了这次翻转,顺式会写成反式:分子变了,而且变得静悄悄。
/// 写一根键的符号,顺带把方向键的**整体翻转自由度**按输出顺序定死。
///
/// 同一约束片段内各键的方向互相锁死,整体翻转不改变任何一对取代基的相对位置 ——
/// 对分子而言那是个真自由度。可 [`directions_for_writing`] 是按**键的存储下标**
/// 取的种子,而存储下标是输入写法留下的痕迹:同一个分子换一种写法读进来,规范串
/// 里的 `/` 与 `\` 就整体互换了。实测语料 8831 条里有 118 条这样。
///
/// 所以 [`WriteStyle::Canonical`] 下:每个片段**第一个写出来的**方向符号一律取
/// `/`,同片段其余的跟着它走。忠实写法不动,那边要的是原样再现。
///
/// [`directions_for_writing`]: crate::stereo::directions_for_writing
fn gauged_symbol(
    bond: u32,
    from: u32,
    mol: &MolBuilder,
    dirs: &[BondDirection],
    comps: &[Option<u32>],
    style: WriteStyle,
    gauge: &mut BTreeMap<u32, bool>,
) -> &'static str {
    let sym = bond_symbol(bond, from, mol, dirs);
    if style != WriteStyle::Canonical || (sym != "/" && sym != "\\") {
        return sym;
    }
    // 只有由约束定下方向的键才有可翻的自由度
    let Some(comp) = comps.get(bond as usize).copied().flatten() else {
        return sym;
    };
    // 该片段第一次露面时定调:让它写成 `/`
    let flip = *gauge.entry(comp).or_insert(sym == "\\");
    match (flip, sym) {
        (true, "/") => "\\",
        (true, _) => "/",
        (false, s) => s,
    }
}

fn bond_symbol(bond: u32, from: u32, mol: &MolBuilder, dirs: &[BondDirection]) -> &'static str {
    let b = mol.bonds()[bond as usize];
    match b.order {
        BondOrder::Double => "=",
        BondOrder::Triple => "#",
        BondOrder::Quadruple => "$",
        BondOrder::Dative => {
            if b.begin == from {
                "->"
            } else {
                "<-"
            }
        }
        // 芳香键在两端都芳香时是默认,不必写 —— 除非它还带着方向。
        // 双键挂在芳香环外时,指方向的正是环上的芳香键。
        BondOrder::Aromatic => match direction_from(b, from, dirs[bond as usize]) {
            BondDirection::UpRight => "/",
            BondDirection::DownRight => "\\",
            BondDirection::None => {
                if both_aromatic(mol, b.begin, b.end) {
                    ""
                } else {
                    ":"
                }
            }
        },
        // 单键通常省略,但有两种情形非写不可:带方向的,以及两端都芳香的
        // (那时默认是芳香键 —— 联苯的两个环之间就是这种情形)
        BondOrder::Single | BondOrder::Unspecified => {
            match direction_from(b, from, dirs[bond as usize]) {
                BondDirection::UpRight => "/",
                BondDirection::DownRight => "\\",
                BondDirection::None => {
                    if both_aromatic(mol, b.begin, b.end) {
                        "-"
                    } else {
                        ""
                    }
                }
            }
        }
    }
}

/// 把存储的方向换算到"从 `from` 走向另一端"的参照系。
fn direction_from(b: BondData, from: u32, stored: BondDirection) -> BondDirection {
    if b.begin == from {
        stored
    } else {
        stored.flipped()
    }
}

fn both_aromatic(mol: &MolBuilder, a: u32, b: u32) -> bool {
    let at = mol.atoms();
    at[a as usize].flags.contains(AtomFlags::AROMATIC)
        && at[b as usize].flags.contains(AtomFlags::AROMATIC)
}

/// 写一个原子。`tag` 是已经换算到**输出顺序**的立体标记。
fn write_atom(
    out: &mut String,
    mol: &MolBuilder,
    idx: u32,
    stereo: (ChiralTag, u8),
    style: WriteStyle,
) {
    let (tag, perm) = stereo;
    let a = mol.atoms()[idx as usize];
    let aromatic = a.flags.contains(AtomFlags::AROMATIC);

    // 通配原子:只要没别的要说,`*` 就够
    if a.atomic_num == 0 && !needs_brackets(mol, idx, style) {
        out.push('*');
        return;
    }

    if !needs_brackets(mol, idx, style) {
        let sym = element::by_atomic_num(a.atomic_num).map_or("*", |e| e.symbol);
        if aromatic {
            out.push_str(&sym.to_ascii_lowercase());
        } else {
            out.push_str(sym);
        }
        return;
    }

    out.push('[');
    if a.isotope != 0 {
        let _ = write!(out, "{}", a.isotope);
    }
    if a.atomic_num == 0 {
        out.push('*');
    } else {
        let sym = element::by_atomic_num(a.atomic_num).map_or("*", |e| e.symbol);
        if aromatic {
            out.push_str(&sym.to_ascii_lowercase());
        } else {
            out.push_str(sym);
        }
    }
    // 立体标记写在元素符号之后、氢数之前
    match tag {
        ChiralTag::Ccw => out.push('@'),
        ChiralTag::Cw => out.push_str("@@"),
        // 序号为 0 表示"这个标记表达不出来",见 `output_chiral_tag`
        ChiralTag::SquarePlanar if perm != 0 => {
            let _ = write!(out, "@SP{perm}");
        }
        ChiralTag::TrigonalBipyramidal if perm != 0 => {
            let _ = write!(out, "@TB{perm}");
        }
        ChiralTag::Octahedral if perm != 0 => {
            let _ = write!(out, "@OH{perm}");
        }
        // 丙二烯轴手性(`@AL`)尚不写出,见 output_chiral_tag
        _ => {}
    }
    match total_hs(&a) {
        0 => {}
        1 => out.push('H'),
        k => {
            let _ = write!(out, "H{k}");
        }
    }
    match a.formal_charge.cmp(&0) {
        std::cmp::Ordering::Greater => {
            out.push('+');
            if a.formal_charge > 1 {
                let _ = write!(out, "{}", a.formal_charge);
            }
        }
        std::cmp::Ordering::Less => {
            out.push('-');
            if a.formal_charge < -1 {
                let _ = write!(out, "{}", -i32::from(a.formal_charge));
            }
        }
        std::cmp::Ordering::Equal => {}
    }
    if a.atom_map != 0 {
        let _ = write!(out, ":{}", a.atom_map);
    }
    out.push(']');
}

/// 该原子是否必须写成方括号形式。
///
/// 方括号一旦出现,氢数就由字面决定、不再推断,所以这个判断同时决定了
/// 氢数怎么表达。判据是"简写形式表达不了这个原子":
fn needs_brackets(mol: &MolBuilder, idx: u32, style: WriteStyle) -> bool {
    if hard_bracket(mol, idx) {
        return true;
    }
    let a = mol.atoms()[idx as usize];
    // 作者钉死过氢数的原子才谈得上要不要留框。
    //
    // 光看 `NO_IMPLICIT` 是不够的 —— 净化会把这个标志清掉,同时把氢挪进
    // `num_explicit_hs`(第 12 步)。于是"净化之后写出"会把吡咯型氮的 `[nH]`
    // 写成裸 `n`,氢凭空消失,写出的串连凯库勒化都做不到。
    // 实测:8839 条语料净化后写出,633 条因此坏掉。
    let author_fixed_hs = a.flags.contains(AtomFlags::NO_IMPLICIT) || a.num_explicit_hs != 0;
    match style {
        WriteStyle::Faithful => author_fixed_hs,
        WriteStyle::Canonical => author_fixed_hs && !hs_survive_without_brackets(mol, idx),
    }
}

/// 简写形式表达不了这个原子 —— 与氢数无关的那几条。
fn hard_bracket(mol: &MolBuilder, idx: u32) -> bool {
    let a = mol.atoms()[idx as usize];
    a.isotope != 0
        || a.formal_charge != 0
        || a.atom_map != 0
        || a.num_radical_electrons != 0
        || a.chiral_tag != ChiralTag::Unspecified
        // 有机子集之外的元素没有简写形式
        || (a.atomic_num != 0 && !element::is_organic_subset(a.atomic_num))
        // 小写形式只对少数几个元素有定义
        || (a.flags.contains(AtomFlags::AROMATIC)
            && !element::can_be_aromatic_lowercase(a.atomic_num))
}

/// 去掉方括号之后,再读回来氢数还是不是原来那个。
///
/// # 判据必须与**解析侧**同口径 —— 现在是**同一份代码**
///
/// 简写形式的氢数由价反推,而反推的规则住在
/// [`omgkit_core::valence::implicit_hs_for_bare_form`]。所以这里问的是一个很具体
/// 的问题:**按那条规则,裸写形式会补出几个氢?** 与本原子实际的氢数相等才能去框。
///
/// 先前这里自己写了一份近似(芳香分支只看首位默认价、非芳香分支用
/// `default_valence_for`),注释里明写着"一处已知的不同步:
/// `explicit_valence_of` 还有一步芳香价回落,这边没有 …… 两处规则分处两个
/// crate 各写一遍、靠人同步"。现在两处走同一份代码,不同步这件事从结构上没有了。
///
/// 实测两条规则在全语料上**确有分歧**:`large.smi` 的 133 537 个(中性、无同位素、
/// 无自由基的)原子里有 5 746 处不同,全是"旧的说算不准、新的说补 0 个氢",
/// 典型是并环芳香碳(三根芳香键,键级和 4.5;旧规则拿首位默认价 4 一比就放弃,
/// 新规则先做芳香价回落到 4、再算出 0 个氢 —— 而那正是对的)。
///
/// 但**写出的结果一行没变**(大语料两个方向 + 冒烟语料,逐行相同):
/// 那些原子本来就走不到这个函数。也就是说这次是把一个**没在发生**的分岔
/// 从结构上堵掉,不是修一个正在漏的洞。
///
/// 先前这里只看**首位**默认价(`bonds >= valences[0]` / `bonds + hs == valences[0]`),
/// 而多价元素的首位默认价根本不是读者会选的那个:
///
/// | | 键级和 | 价表 | 读者补的氢 | 先前判定 |
/// |---|---|---|---|---|
/// | `Cl[I]Cl` 的 I | 2 | `[1, 3, 5]` | 1(补到 3) | `2 >= 1` → 去框 |
/// | `NC[S](=O)=O` 的 S | 5 | `[2, 4, 6]` | 1(补到 6) | `5 >= 2` → 去框 |
///
/// 去框之后写出的是 `ClICl` / `NCS(=O)=O`,任何读者(**包括我们自己的净化**)
/// 读回来都多一个氢 —— 写出来的是另一个分子。全语料 9 条如此,
/// 外部判据 `harness/check_write.py --canonical` 抓的就是这个。
///
/// 往返测试抓不住:本模块的解析器**不推隐式氢**(`num_implicit_hs` 恒为 0),
/// 氢要到净化才出现,于是解析与写出在这一点上是共谋的。
///
/// # 超价一律留框
///
/// 已用价超过该元素**全部**允许价时(`default_valence_for` 给 `None`),
/// 读者补几个氢**取决于它那份价表有多长**,而价表各版本各实现都不一样。
/// 实测(全语料只有一个这样的原子:六根键的中性 P,`large.smi` 第 4554 行的
/// PF₆):
///
/// | 价表来源 | P 的价表 | 裸写 `FP(F)(F)(F)(F)F` 读回来 |
/// |---|---|---|
/// | 本仓 `element_data.rs`(转录自 RDKit **2025.09.2**) | `[3, 5]` | 补 0 个氢 |
/// | RDKit **2025.09.2** | `[3, 5]` | **净化直接失败**(6 > 5) |
/// | RDKit **2022.09.5** | `[3, 5, 7]` | 补 1 个氢 → 是另一个分子 |
///
/// 三种结果两两不同,所以**算不准就留框**:留框在三种读者下都是对的。
///
/// (这段先前写的是"我们的 P 价表是 `[3,5]`,RDKit 的是 `[3,5,7]`"。那个
///  `[3,5,7]` 是**开发机 `.venv` 里那个 2022.09.5** 的,而本仓的元素表恰恰就是
///  从 2025.09.2 转录的 —— 两边其实一样。把没核过的版本差写成"我们 vs 别人"
///  的分歧,正是 `harness/requirements.lock` 里记过的那个坑。)
///
/// # 其余的保守取舍
///
/// 一个中性碳(`[C]`)氢数是 0,可去掉方括号写成 `C` 再读回来就补上了
/// 四个氢 —— 所以"没氢"绝不等于"能去框"。
///
/// # 配位键先前是无条件留框的,那让同一个分子有了两串
///
/// 先前这里有一条"邻居里有配位键就直接返回 false",理由写的是"配位键的给体端
/// 不计价,这一带的价本就不好谈"。它的后果是:`N->[Cu]` 与 `[NH3]->[Cu]` 是同一个
/// 分子(氮都是 3 个氢),而净化把前者的氢放在 `num_implicit_hs`、后者放在
/// `num_explicit_hs`,于是**只有后者走到这个函数**,被无条件留框 ——
/// 一个写成 `N->[Cu]`,一个写成 `[NH3]->[Cu]`。方括号是书写习惯,不是分子的性质,
/// 规范串不能跟着它变(`tests/canonical_invariance.rs` 有一条判据专门守这个)。
///
/// 那条兜底现在没有了:氢数由 [`omgkit_core::valence::implicit_hs_for_bare_form`]
/// 算,而那条规则本来就把配位键的价贡献算对了(给体 0、受体 1,见
/// `BondData::valence_contribution_to`),与下面那个 `bonds` 求和同一个约定。
/// 读者按同一条规则反推,裸写形式补出来的氢数一样。
///
/// 抓住它的是 `tests/differential_l3.rs`(拿 RDKit 的规范串当"另一种写法"):
/// 冒烟语料 149 条里 5 条因此不收敛,全是配位键那一带。
///
/// # 前置条件:调用方已经挡掉了带电 / 自由基 / 同位素的原子
///
/// 本函数**不做**形式电荷与自由基的调整,而 `implicit_hs_of` 是做的
/// (`effective_atomic_num`、`can_be_hypervalent`、`+ n_radicals`)。
/// 全语料实测有 1330 多个原子两边算出的氢数不同,**无一例外全是带电或自由基**
/// (O⁻ 1121、N⁺ 135、S⁺ 25、N⁻ 25…)。它们全部被
/// [`hard_bracket`] 提前判成"必须留框",走不到这里。
///
/// 这是**调用方保证的前置条件,不是本函数的性质** —— 所以下面有一句
/// `debug_assert!`。哪天有人把调用顺序换了,那 1330 个原子会静默地按错规则算。
fn hs_survive_without_brackets(mol: &MolBuilder, idx: u32) -> bool {
    debug_assert!(
        !hard_bracket(mol, idx),
        "hs_survive_without_brackets 的前置条件被破坏了:带电/自由基/同位素的原子\
         必须由 hard_bracket 提前挡掉,本函数不做那几项调整 —— 见函数文档"
    );
    let a = mol.atoms()[idx as usize];
    let Some(e) = element::by_atomic_num(a.atomic_num) else {
        return false;
    };
    if !e.has_valence_constraint() {
        return false;
    }
    let bonds: f32 = mol
        .neighbors(idx)
        .map(|(_, bi)| mol.bonds()[bi as usize].valence_contribution_to(idx))
        .sum();
    // x.5 向上取整(芳香键各计 1.5),与价键计算同一约定
    #[allow(clippy::cast_possible_truncation)]
    let bonds = (bonds + 0.1).round() as i32;
    let Ok(used) = i8::try_from(bonds) else {
        return false; // 键级和大得离谱,算不准 → 留框
    };
    let total_hs = i32::from(a.num_explicit_hs) + i32::from(a.num_implicit_hs);

    // 超价一律留框 —— 读者补几个氢取决于它那份价表有多长,见函数文档。
    // 这问的是"算不算得准",不是"补几个氢",所以留在这一侧。
    if e.default_valence_for(used).is_none() {
        return false;
    }

    // **氢数由那条唯一的规则算**,不在这里另写一份。
    //
    // 先前这里是自己写的一份近似:芳香分支只看首位默认价、非芳香分支用
    // `default_valence_for`。注释里还明写着"一处已知的不同步:
    // `explicit_valence_of` 还有一步芳香价回落(差在 1.5 以内就取该价态),
    // 这边没有 …… 两处规则分处两个 crate 各写一遍、靠人同步"。
    //
    // 现在两处走同一份代码([`omgkit_core::valence`]),不同步这件事从结构上
    // 没有了。那条芳香保守分支也一并没了 —— 它守的
    // `c1cc[sH]c1`(S 两根芳香键、1 个氢)现在由规则本身给出正确答案:
    // 裸写形式补 0 个氢 ≠ 实际的 1 个,照样留框。
    let Some(bare) = omgkit_core::valence::implicit_hs_for_bare_form(mol, idx) else {
        return false; // 电荷/同位素/自由基裸写表达不出来 —— 本来就该留框
    };
    i32::from(bare) == total_hs
}
