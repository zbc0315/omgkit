//! 规范化排序的判据:**原子任意重排,规范 SMILES 不变**。
//!
//! 这是整条管线里唯一不需要外部参照就能验证的性质,而且抓得住绝大多数排序
//! 错误 —— 一个依赖输入编号的实现照样能产出一个"看起来像样"的全序,只有
//! 换个编号才露馅。
//!
//! # 重排必须连键序一起换
//!
//! 只打乱原子编号是不够的。邻居的存储顺序等于**建键顺序**,规范化若不小心
//! 读了它,单纯换原子编号的测试可能仍然通过。所以本测试同时打乱两者。
//!
//! # 判据是字符串,不是秩
//!
//! 比"重排前后的秩数组"要先把秩映射回去,绕一圈还容易把比对本身写错。
//! 直接比最终产物 —— 规范 SMILES 字符串 —— 既简单又是真正关心的东西。

use std::path::{Path, PathBuf};

use omgkit_core::{BondData, MolBuilder};
use omgkit_io::{canon, smiles};

fn corpus(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../harness/corpus")
        .join(name)
}

/// 定死种子的 xorshift。测试要可复现 —— 随机失败没法查。
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    /// Fisher-Yates
    fn shuffle<T>(&mut self, v: &mut [T]) {
        for i in (1..v.len()).rev() {
            let j = (self.next() % (i as u64 + 1)) as usize;
            v.swap(i, j);
        }
    }
}

/// 置换的宇称。两者不是同一多重集时返回 `None`。
fn permutation_is_odd(from: &[usize], to: &[usize]) -> Option<bool> {
    if from.len() != to.len() {
        return None;
    }
    let mut cur: Vec<usize> = from.to_vec();
    let mut swaps = 0usize;
    for i in 0..to.len() {
        if cur[i] == to[i] {
            continue;
        }
        let j = (i + 1..cur.len()).find(|&j| cur[j] == to[i])?;
        cur.swap(i, j);
        swaps += 1;
    }
    Some(swaps % 2 == 1)
}

/// 按给定的原子置换与键序重排,造一个**同构但编号不同**的分子。
///
/// # 手性标记必须跟着换参照系
///
/// 标记的含义是相对**邻居存储顺序**的,而存储顺序等于建键顺序 —— 打乱键序
/// 就改变了参照系。照抄 `chiral_tag` 造出来的不是同一个分子,而是它在部分
/// 手性中心上的镜像。
///
/// 这个坑很实在:重排测试本来就是用来抓"实现读了输入编号"的,若重排本身
/// 就把分子改了,测试会报出一堆并非缺陷的失败,进而诱使人去"修"正确的实现。
fn renumber(mol: &MolBuilder, atom_perm: &[u32], bond_perm: &[usize]) -> MolBuilder {
    // atom_perm[new] = old
    let mut old_to_new = vec![0u32; mol.num_atoms()];
    for (new, &old) in atom_perm.iter().enumerate() {
        old_to_new[old as usize] = new as u32;
    }
    // 新键下标 → 旧键下标
    let new_to_old_bond: Vec<usize> = bond_perm.to_vec();

    let mut out = MolBuilder::new();
    for &old in atom_perm {
        out.add_atom_data(mol.atoms()[old as usize]);
    }
    for &bi in bond_perm {
        let b = mol.bonds()[bi];
        let mut nb = BondData::new(
            old_to_new[b.begin as usize],
            old_to_new[b.end as usize],
            b.order,
        );
        nb.direction = b.direction;
        nb.stereo = b.stereo;
        nb.flags = b.flags;
        out.add_bond_data(nb).expect("端点合法");
    }

    // 逐个立体中心把标记搬到新的存储序上
    for new_a in 0..out.num_atoms() as u32 {
        let old_a = atom_perm[new_a as usize];
        let a = mol.atoms()[old_a as usize];
        let old_seq: Vec<usize> = mol.neighbors(old_a).map(|(_, b)| b as usize).collect();
        let new_seq: Vec<usize> = out
            .neighbors(new_a)
            .map(|(_, b)| new_to_old_bond[b as usize])
            .collect();
        if a.chiral_tag.is_tetrahedral() {
            if permutation_is_odd(&old_seq, &new_seq).expect("同一组键") {
                out.atom_mut(new_a).expect("原子存在").chiral_tag = a.chiral_tag.inverted();
            }
            continue;
        }
        // 配位几何的序号同样是**相对存储序**的,照抄就是换了个分子 ——
        // 与上面那条四面体的规则同源,只是变换不是"翻一下"而是按该多面体的
        // 转动群换一个排法。
        if omgkit_core::polyhedron::ligand_count(a.chiral_tag).is_some() {
            let old_ids: Vec<u32> = old_seq.iter().map(|&b| b as u32).collect();
            let new_ids: Vec<u32> = new_seq.iter().map(|&b| b as u32).collect();
            let p =
                omgkit_core::polyhedron::renumber(a.chiral_tag, a.stereo_perm, &old_ids, &new_ids)
                    .unwrap_or(a.stereo_perm);
            out.atom_mut(new_a).expect("原子存在").stereo_perm = p;
        }
    }
    out
}

/// 解析 + 净化 + 感知顺反 —— 调用方真正走的那条路。
///
/// 感知这一步不能省:`stereo_atoms` 是它填的,而"参照原子挑在哪一侧"正是
/// 规范串不动点的要害。少了它,相关判据全部空过。
fn perceived(smi: &str) -> MolBuilder {
    let mut m = smiles::parse(smi).unwrap_or_else(|e| panic!("{smi}: {}", e.render()));
    omgkit_chem::sanitize(&mut m).unwrap_or_else(|e| panic!("{smi}: {e}"));
    omgkit_io::stereo::perceive_bond_stereo(&mut m);
    m
}

fn canonical_smiles(mol: &MolBuilder) -> String {
    canon::canonical_smiles(mol).smiles
}

struct Failure {
    smi: String,
    base: String,
    other: String,
    round: usize,
}

#[derive(Default)]
struct Stats {
    molecules: usize,
    /// 细化后仍有对称等价原子的分子数 —— 打破对称那一步真正生效的用例
    needed_tie_breaking: usize,
}

/// 对语料里每个分子做 `rounds` 次随机重排,断言规范 SMILES 恒等。
fn check_corpus(path: &Path, rounds: usize, seed: u64) -> (Stats, Vec<Failure>) {
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("读不到语料 {}: {e}", path.display()));
    let mut rng = Rng(seed);
    let mut stats = Stats::default();
    let mut bad: Vec<Failure> = Vec::new();

    for line in text.lines() {
        let smi = line.split_whitespace().next().unwrap_or("");
        if smi.is_empty() || smi.starts_with('#') {
            continue;
        }
        let Ok(mol) = smiles::parse(smi) else {
            continue;
        };
        if mol.num_atoms() < 2 {
            continue; // 单原子无从重排
        }
        stats.molecules += 1;

        let base = canonical_smiles(&mol);
        // 秩数少于原子数说明细化没能分开所有原子,是打破对称在起作用
        if has_symmetry(&mol) {
            stats.needed_tie_breaking += 1;
        }

        for round in 0..rounds {
            let mut ap: Vec<u32> = (0..mol.num_atoms() as u32).collect();
            let mut bp: Vec<usize> = (0..mol.num_bonds()).collect();
            rng.shuffle(&mut ap);
            rng.shuffle(&mut bp);
            let other = canonical_smiles(&renumber(&mol, &ap, &bp));
            if other != base {
                bad.push(Failure {
                    smi: smi.to_string(),
                    base: base.clone(),
                    other,
                    round,
                });
                break; // 一个分子报一次就够
            }
        }
    }
    (stats, bad)
}

/// 分子里是否存在互相等价的原子(据"重排后规范串仍相同但原子可互换"判断
/// 代价太高,这里用一个便宜的充分条件:存在两个原子的初始不变量与邻居
/// 元素多重集完全相同)。
///
/// 只用来统计覆盖,不参与判定,所以宁可漏报也不要拖慢主循环。
fn has_symmetry(mol: &MolBuilder) -> bool {
    let mut keys: Vec<(u8, usize, Vec<u8>)> = (0..mol.num_atoms() as u32)
        .map(|a| {
            let mut nbrs: Vec<u8> = mol
                .neighbors(a)
                .map(|(o, _)| mol.atoms()[o as usize].atomic_num)
                .collect();
            nbrs.sort_unstable();
            (mol.atoms()[a as usize].atomic_num, mol.degree(a), nbrs)
        })
        .collect();
    keys.sort();
    keys.windows(2).any(|w| w[0] == w[1])
}

fn report(bad: &[Failure], limit: usize) -> String {
    let mut out = format!("\n规范 SMILES 不随重排恒定,{} 条分子失败:\n\n", bad.len());
    for f in bad.iter().take(limit) {
        out.push_str(&format!(
            "  原文:{}\n  重排前:{}\n  第 {} 次重排后:{}\n\n",
            f.smi, f.base, f.round, f.other
        ));
    }
    if bad.len() > limit {
        out.push_str(&format!("  ...(另有 {} 条)\n", bad.len() - limit));
    }
    out
}

#[test]
fn canonical_smiles_is_invariant_under_renumbering_smoke() {
    let (stats, bad) = check_corpus(&corpus("smoke.smi"), 20, 0x9E37_79B9_7F4A_7C15);
    assert!(bad.is_empty(), "{}", report(&bad, 10));
    assert!(stats.molecules > 100, "只测到 {} 条分子", stats.molecules);
    assert!(
        stats.needed_tie_breaking > 0,
        "语料里没有含等价原子的分子,打破对称那一步是空过的"
    );
    println!(
        "规范化重排不变(冒烟):{} 条分子 × 20 次重排全部恒等;其中含等价原子 {} 条",
        stats.molecules, stats.needed_tie_breaking
    );
}

#[test]
// **不再 `#[ignore]`。** 实测**几秒量级**(不同机器/负载下 2~8 s 都量到过 ——
// 墙钟值别写成单点,见 `omgkit-match/tests/scaling.rs` 的同一条)。见 `adjacency_index.rs`
// 里同一处说明:标着 `#[ignore]` 等于这条判据从来没在 CI 里跑过。
fn canonical_smiles_is_invariant_under_renumbering_corpus() {
    let (stats, bad) = check_corpus(&corpus("large.smi"), 5, 0x243F_6A88_85A3_08D3);
    assert!(bad.is_empty(), "{}", report(&bad, 15));
    assert!(stats.molecules > 8000, "只测到 {} 条分子", stats.molecules);
    // **覆盖断言,与冒烟档同一套。** 冒烟档断了 `needed_tie_breaking > 0`,
    // 这一档先前一条都没有 —— 而它现在进了 CI:分母敞着的话,"全部恒等"
    // 可能只说明打破对称那条路根本没走到。现值 8717 条。
    assert!(
        stats.needed_tie_breaking > 8000,
        "只有 {} 条分子含等价原子 —— 打破对称那条路几乎没走到",
        stats.needed_tie_breaking
    );
    println!(
        "规范化重排不变(大语料):{} 条分子 × 5 次重排全部恒等;其中含等价原子 {} 条",
        stats.molecules, stats.needed_tie_breaking
    );
}

/// 同一个分子的不同写法必须收敛到同一个规范串。
///
/// 这条比重排测试更强:重排只换编号,这里连**书写结构**都变了
/// (起点不同、分支顺序不同、环闭合标号不同)。
///
/// 注意各组内必须是**同一个 L1 图**。`c1ccccc1` 与 `C1=CC=CC=C1` 化学上是
/// 同一个分子,但在 L1 一个是芳香标志、一个是交替单双键,要经过 L2 的芳香性
/// 感知才等价 —— 那不是规范化排序该负责的事。
#[test]
fn different_spellings_converge_to_one_canonical_form() {
    for group in [
        &["OCC", "CCO", "C(O)C"][..],
        &["OC(=O)c1ccccc1N", "Nc1ccccc1C(O)=O", "c1cc(N)c(cc1)C(=O)O"][..],
        &["C1CCCCC1", "C1(CCCCC1)"][..],
        &["CC(C)C", "C(C)(C)C"][..],
        &["N[C@@H](C)C(=O)O", "C([C@@H](N)C)(=O)O"][..],
    ] {
        let mut canonical: Option<String> = None;
        for smi in group {
            let m = smiles::parse(smi).unwrap_or_else(|e| panic!("{smi}: {}", e.render()));
            let c = canonical_smiles(&m);
            match &canonical {
                None => canonical = Some(c),
                Some(first) => {
                    assert_eq!(first, &c, "{smi} 与 {} 是同一个分子,规范串却不同", group[0])
                }
            }
        }
    }
}

/// 顺反的**整体翻转**也是书写习惯:`/C=C/` 与 `\C=C\` 是同一个几何。
///
/// 一根双键两侧的方向键同时取反,任何一对取代基的相对位置都不变。对分子而言
/// 那是个真自由度,对规范串而言必须定死,否则同一个分子有两串。
///
/// 关键在于**没跑过 `perceive_bond_stereo` 的分子也要定死**。跑过感知的分子,
/// 方向是从 `stereo` 重新生成的,写出器早就按约束片段把翻转定死了;没跑过的
/// 分子只剩沿用来的方向,一个约束都没有,那一段自由度会悬空。
///
/// `run_reactants` 交出来的产物正是后者 —— 这一条最初就是在反应产物上露的馅:
/// 同一个肉桂酸,一条从 SMILES 直接读进来、一条由反应生成,两串规范式互为
/// 整体翻转。
#[test]
fn canonical_is_independent_of_the_global_direction_flip() {
    for group in [
        &["F/C=C/F", "F\\C=C\\F", "C(\\F)=C/F"][..],
        &[
            "OC(=O)/C=C/c1ccccc1",
            "OC(=O)\\C=C\\c1ccccc1",
            "c1ccccc1/C=C/C(O)=O",
            "C(=C\\C(O)=O)/c1ccccc1",
        ][..],
        // 共轭链:整条链一起翻
        &["C/C=C/C=C/C", "C\\C=C\\C=C\\C"][..],
        // 顺式那一支同样要定死
        &["C/C=C\\C", "C\\C=C/C"][..],
    ] {
        let mut first: Option<String> = None;
        for smi in group {
            let m = smiles::parse(smi).unwrap_or_else(|e| panic!("{smi}: {}", e.render()));
            let c = canonical_smiles(&m);
            match &first {
                None => first = Some(c),
                Some(f) => assert_eq!(f, &c, "{smi} 与 {} 是同一个几何,规范串却不同", group[0]),
            }
        }
    }
}

/// 方括号是**书写习惯**,不是分子的性质 —— 规范串不能跟着它变。
///
/// `[CH3][CH2][OH]` 与 `CCO` 是同一个分子,只是前者把氢写在方括号里。写出器要
/// 往返恒等,所以照原样再现方括号;规范化要的却是"同一个分子只有一串",两者
/// 冲突,由 `BracketStyle` 分开。这一条守规范化那一侧。
///
/// 上面那条 `different_spellings_converge_to_one_canonical_form` 盖不到:它换的
/// 是**遍历起点与分支顺序**,原子的表示没变过。判据的洞不是判据松,是语料没这个
/// 形态 —— 修之前那一条是绿的,而 `CCO` 与 `[CH3][CH2][OH]` 给出两串。
#[test]
fn canonical_is_independent_of_bracket_notation() {
    for group in [
        // 常见原子多写了方括号 —— 去掉框读回来氢数不变
        &["CCO", "[CH3][CH2][OH]", "[CH3]C[OH]"][..],
        &["c1ccccc1", "[cH]1[cH][cH][cH][cH][cH]1"][..],
        &["CC(=O)O", "[CH3][C](=[O])[OH]"][..],
        &["CS(=O)(=O)O", "[CH3][S](=[O])(=[O])[OH]"][..],
        &[
            "CN(C)c1ccc(N(C)C)cc1",
            "CN(C)[c]1[cH][cH][c]([cH][cH]1)[N](C)C",
        ][..],
        // 框**必须**留着的:去掉之后氢数会变
        &["c1cc[nH]c1", "[cH]1[cH][cH][nH][cH]1"][..],
        // 吡啶氮无氢,框可去
        &["c1ccncc1", "[cH]1[cH][cH][n][cH][cH]1"][..],
        // 平面四方:序号说的是“哪两对互为反位”,而两个氯彼此可换 ——
        // `@SP1` 与 `@SP3` 在这个分子上表达同一件事,必须收敛到一串。
        // 判据是参照自己给的:RDKit 把它重排原子再规范化,会在这两串之间跳。
        &["[Pt@SP1](Cl)(Cl)(N)N", "[Pt@SP3](Cl)(Cl)(N)N"][..],
        // 四个配体互不相同时三种排法各是一个分子,这里只验其中一种的两种写法:
        // 交换列出顺序里的两个配体,序号跟着换,分子不变。
        &["[Pt@SP1](Cl)(Br)(N)P", "[Pt@SP2](Cl)(N)(Br)P"][..],
        // 配位键的给体端:写不写框都是同一个分子(氮都是 3 个氢)。
        // 写出器先前碰到配位键无条件留框,于是这一组给出两串 ——
        // 是 `differential_l3.rs` 拿 RDKit 的规范串比出来的。
        &["N->[Cu]", "[NH3]->[Cu]", "[Cu]<-N"][..],
        &[
            "O->[Fe](<-O)(<-O)<-O",
            "[OH2]->[Fe](<-[OH2])(<-[OH2])<-[OH2]",
        ][..],
    ] {
        let mut canonical: Option<String> = None;
        for smi in group {
            // 不净化:方括号里的氢本来就记在显式一侧,判据不依赖净化。
            // 少一层依赖,这条判据坏了就一定是写出器的事。
            let m = smiles::parse(smi).unwrap_or_else(|e| panic!("{smi}: {}", e.render()));
            let c = canonical_smiles(&m);
            match &canonical {
                None => canonical = Some(c),
                Some(first) => {
                    assert_eq!(first, &c, "{smi} 与 {} 是同一个分子,规范串却不同", group[0])
                }
            }
        }
    }
}

/// 该留的框不能被"能省则省"去掉 —— 去掉之后读回来就是另一个分子。
///
/// 与上一条是一对:那条守"同一个分子一串",这条守"别把分子写坏"。只有前者的话,
/// 一个把所有框都去掉的实现也能通过。
#[test]
fn brackets_that_carry_hydrogen_counts_are_kept() {
    for smi in [
        "c1cc[nH]c1", // 吡咯氮:去框读回来氢没了,连凯库勒化都做不到
        "[CH2]C",     // 卡宾/双自由基碳:去框会补满氢
        "[nH]1cccc1", // 同吡咯,起点不同
    ] {
        let m = smiles::parse(smi).unwrap_or_else(|e| panic!("{smi}: {}", e.render()));
        let c = canonical_smiles(&m);
        assert!(
            c.contains('['),
            "{smi} 的规范串 {c} 丢了方括号 —— 氢数会跟着变"
        );
    }
}

/// 环闭合落在平面四方中心上 —— **书写序 ≠ 存储序**,归一那一步只有在这里才做事。
///
/// 语料里的 `@SP` 分子(`[Pt@SP1](Cl)(Cl)(N)N` 那几条)书写序恰好等于存储序,
/// 于是解析侧那一步归一是空操作:变异实测**把它整个删掉,全套判据一条都不红**。
/// 环闭合键是在闭合的那一刻才建的,排在存储序末尾,而它在串里写在最前面 ——
/// 这个分子把两者岔开。
///
/// # 这个分子上"三种排法"只有两个分子
///
/// 两个环碳等价(都只连着 Pt 与彼此),交换它们是一个自同构。按列出顺序
/// `(C_闭环, Cl, N, C_链)`:
///
/// | 写法 | 反位配对 | 自同构映到 |
/// |---|---|---|
/// | `@SP1` | C_闭环–N、Cl–C_链 | `@SP2` 的那个 |
/// | `@SP2` | C_闭环–Cl、N–C_链 | `@SP1` 的那个 |
/// | `@SP3` | C_闭环–C_链、Cl–N | 自己 |
///
/// 所以 `@SP1` 与 `@SP2` 是同一个分子,`@SP3` 是另一个。我方给两串。
/// **参照(RDKit 2025.09.2)给三串** —— `[NH2][Pt@SP2]1([Cl])[CH2][CH2]1`、
/// `[NH2][Pt@SP3]1(...)`、`[NH2][Pt@SP1]1(...)` 各一 —— 它的 `@SP` 规范化不考虑
/// 这个自同构,与它在 `[Pt@SP1](Cl)(Cl)(N)N` 上重排原子会跳串是同一件事。
#[test]
fn square_planar_with_a_ring_closure_on_the_centre() {
    let one = canonical_smiles(&perceived("[Pt@SP1]1(Cl)(N)CC1"));
    let two = canonical_smiles(&perceived("[Pt@SP2]1(Cl)(N)CC1"));
    let three = canonical_smiles(&perceived("[Pt@SP3]1(Cl)(N)CC1"));
    assert_eq!(
        one, two,
        "@SP1 与 @SP2 在这个分子上是同一个分子(两个环碳可换),规范串却不同"
    );
    assert_ne!(
        one, three,
        "@SP3 是另一个分子(它把两个环碳排成反位),却与 @SP1 塌成了一串"
    );
    // 换个起笔位置再写一遍同样的三个分子,必须落到同一组串上。
    // 起笔一换,书写序与存储序的错位方式也跟着变。
    assert_eq!(one, canonical_smiles(&perceived("C1C[Pt@SP1]1(Cl)N")));
    assert_eq!(three, canonical_smiles(&perceived("C1C[Pt@SP2]1(Cl)N")));
    assert_eq!(one, canonical_smiles(&perceived("C1C[Pt@SP3]1(Cl)N")));

    // **另一种写法由外部实现给出。** 上面几条全是我方自己的写法,而
    // `SQUARE_PLANAR_TRANS` 那张表若整体错位(比如两行对调),解析与写出两侧
    // 会一起错、正好抵消 —— 变异实测:把表的第 2、3 行对调,全套外部判据
    // 照样退 0(语料里那几个 `@SP` 分子书写序恰好等于存储序,错表抵消掉了)。
    //
    // 这三串是 RDKit 2025.09.2 对上面三个分子给出的规范串。它与我方的配体顺序
    // 不同,错表在这里就抵消不掉了。语料里没有"环闭合落在 `@SP` 中心上"的分子,
    // 所以这一对是手写的 —— 但**第二种写法不是我方自己写的**。
    for (ours, reference) in [
        ("[Pt@SP1]1(Cl)(N)CC1", "[NH2][Pt@SP2]1([Cl])[CH2][CH2]1"),
        ("[Pt@SP2]1(Cl)(N)CC1", "[NH2][Pt@SP3]1([Cl])[CH2][CH2]1"),
        ("[Pt@SP3]1(Cl)(N)CC1", "[NH2][Pt@SP1]1([Cl])[CH2][CH2]1"),
    ] {
        assert_eq!(
            canonical_smiles(&perceived(ours)),
            canonical_smiles(&perceived(reference)),
            "{ours} 与外部实现写的 {reference} 是同一个分子,规范串却不同"
        );
    }
}

/// 丙二烯型轴手性,**分组由外部实现钉住**。
///
/// 仓库钉的 RDKit 2025.09.2 在这一档上完全没有能力:SMILES 读写、
/// `FindPotentialStereo`、`rdCIPLabeler`、从 3D 坐标反推、molblock 往返、
/// 带手性的子结构匹配 —— 六条路都把 `@AL1` 与 `@AL2` 当成同一个东西。
/// 所以这里的参照是 **Indigo 1.46.0**,它认这个区别而且自洽
/// (`@AL1` ≡ `@`,读→写链上每一步都还是同一个分子)。
///
/// 下表里"第几组"全部由 Indigo 的 `exactMatch(…, "ALL")` 判定,每组还带一串
/// **它自己写出来的**写法(把 `@` / `@@` 写在中心上的那种形式)。
///
/// **一族里必须放几种把配体角色拆开的写法。** 只放"同一个骨架换个序号"是不够的
/// —— 那时参照那一侧配体的排布与我方相同,插错位置两边一起错、正好抵消。
/// 变异实测:把端上那个氢排到该端配体的末尾(而不是"紧跟前驱原子"),
/// 只有 `NC=[C@AL{i}]=C(O)F` 一种形状时**全绿**;加上 `[C@AL{i}](=CN)=C(O)F`
/// (氢端从中心那边写)才红。环那一族同理,要把三元环**反着写**一遍。
#[test]
fn allene_stereo_matches_the_reference() {
    for (family, cases) in [
        (
            "两端各两个取代基",
            &[
                ("NC(Br)=[C@AL1]=C(O)F", 0),
                ("[C@AL1](=C(N)Br)=C(O)F", 0),
                ("OC(F)=[C@AL1]=C(N)Br", 0),
                ("BrC(N)=[C@AL2]=C(O)F", 0),
                ("NC(=[C@@]=C(F)O)Br", 0),
                ("BrC(N)=[C@AL1]=C(O)F", 1),
                ("NC(Br)=[C@AL2]=C(O)F", 1),
                ("[C@AL2](=C(N)Br)=C(O)F", 1),
                ("OC(F)=[C@AL2]=C(N)Br", 1),
                ("BrC(=[C@@]=C(F)O)N", 1),
            ][..],
        ),
        (
            "一端带一个氢",
            &[
                ("NC=[C@AL1]=C(O)F", 0),
                ("N[CH]=[C@AL1]=C(O)F", 0),
                ("[C@AL2](=CN)=C(O)F", 0),
                ("NC=[C@@]=C(F)O", 0),
                ("[C@AL1](=CN)=C(O)F", 1),
                ("NC=[C@AL2]=C(O)F", 1),
                ("N[CH]=[C@AL2]=C(O)F", 1),
                ("[C@@](=C(F)O)=CN", 1),
            ][..],
        ),
        (
            "一端的配体在环里",
            &[
                ("C1(=[C@AL1]=C(N)Br)OC1", 0),
                ("O1CC1=[C@AL1]=C(N)Br", 0),
                ("C1(=[C@AL2]=C(N)Br)CO1", 0),
                ("C1OC1=[C@AL2]=C(N)Br", 0),
                ("C1(CO1)=[C@]=C(Br)N", 0),
                ("C1(=[C@AL1]=C(N)Br)CO1", 1),
                ("C1OC1=[C@AL1]=C(N)Br", 1),
                ("C1(=[C@AL2]=C(N)Br)OC1", 1),
                ("O1CC1=[C@AL2]=C(N)Br", 1),
                ("C1(OC1)=[C@]=C(Br)N", 1),
            ][..],
        ),
    ] {
        let mut canonical: Vec<(String, usize)> = Vec::new();
        for &(smi, group) in cases {
            canonical.push((canonical_smiles(&perceived(smi)), group));
        }
        for (i, (a, ga)) in canonical.iter().enumerate() {
            for (b, gb) in &canonical[i + 1..] {
                if ga == gb {
                    assert_eq!(a, b, "{family}:外部实现说这两条是同一个分子,规范串却不同");
                } else {
                    assert_ne!(a, b, "{family}:外部实现说这两条是两个分子,规范串却塌成一串");
                }
            }
        }
    }
}

/// 缺一个顶点时的序号,**标号由参照实现钉住**。
///
/// 方括号里的氢、或者一个空的配位位置,也占多面体的一个顶点,而它不在键序列
/// 里。它落在哪一位是量出来的(见 `smiles.rs::coordination_ligands`):
/// "自身位置" —— 紧跟前驱原子之后、环闭合之前。
///
/// 这条判据必须由**外部**给第二种写法。自家的写法钉不住标号:插错位置时
/// 解析与写出两侧会一起错、正好抵消 —— 往返照样恒等,三个序号照样给三串。
/// 下面每一对里,右边那串是 RDKit 2025.09.2 对左边那个分子给出的规范串,
/// 它的配体顺序与我方不同,错位在这里就抵消不掉了。
///
/// 六族分别压住:空位 / 方括号里的氢、手性原子居首 / 前面有原子、
/// 中心上有没有环闭合。
#[test]
fn coordination_stereo_with_a_missing_vertex_matches_the_reference() {
    for (family, cases) in [
        (
            "空位·平面四方",
            [
                ("[Pt@SP1](Cl)(N)O", "[NH2][Pt@SP2]([OH])[Cl]"),
                ("[Pt@SP2](Cl)(N)O", "[NH2][Pt@SP1]([OH])[Cl]"),
                ("[Pt@SP3](Cl)(N)O", "[NH2][Pt@SP3]([OH])[Cl]"),
            ],
        ),
        (
            "空位·三角双锥",
            [
                ("[P@TB1](F)(Cl)(Br)I", "F[P@TB9](Cl)(Br)I"),
                ("[P@TB5](F)(Cl)(Br)I", "F[P@TB13](Cl)(Br)I"),
                ("[P@TB12](F)(Cl)(Br)I", "F[P@TB4](Cl)(Br)I"),
            ],
        ),
        (
            "空位·八面体",
            [
                (
                    "[Co@OH1](N)(O)(S)(P)F",
                    "[NH2][Co@OH2]([OH])([F])([PH2])[SH]",
                ),
                (
                    "[Co@OH11](N)(O)(S)(P)F",
                    "[NH2][Co@OH24]([OH])([F])([PH2])[SH]",
                ),
                (
                    "[Co@OH20](N)(O)(S)(P)F",
                    "[NH2][Co@OH12]([OH])([F])([PH2])[SH]",
                ),
            ],
        ),
        (
            "方括号里的氢·平面四方",
            [
                ("[Pt@SP1H](Cl)(N)O", "[NH2][Pt@SP2H]([OH])[Cl]"),
                ("[Pt@SP2H](Cl)(N)O", "[NH2][Pt@SP1H]([OH])[Cl]"),
                ("[Pt@SP3H](Cl)(N)O", "[NH2][Pt@SP3H]([OH])[Cl]"),
            ],
        ),
        (
            "方括号里的氢·八面体",
            [
                (
                    "[Co@OH1H](N)(O)(S)(P)F",
                    "[NH2][Co@OH2H]([OH])([F])([PH2])[SH]",
                ),
                (
                    "[Co@OH11H](N)(O)(S)(P)F",
                    "[NH2][Co@OH24H]([OH])([F])([PH2])[SH]",
                ),
                (
                    "[Co@OH20H](N)(O)(S)(P)F",
                    "[NH2][Co@OH12H]([OH])([F])([PH2])[SH]",
                ),
            ],
        ),
        (
            "前面有原子·八面体",
            [
                ("N[Co@OH1](O)(S)(P)F", "[NH2][Co@OH7]([OH])([F])([PH2])[SH]"),
                (
                    "N[Co@OH11](O)(S)(P)F",
                    "[NH2][Co@OH9]([OH])([F])([PH2])[SH]",
                ),
                (
                    "N[Co@OH20](O)(S)(P)F",
                    "[NH2][Co@OH22]([OH])([F])([PH2])[SH]",
                ),
            ],
        ),
        (
            // 环里的两个原子必须不一样。都是碳的话交换它们是这个分子的自同构,
            // 三种排法只剩两个异构体 —— 而参照实现的规范串在这种情形下**不合并**
            // (见 `square_planar_with_a_ring_closure_on_the_centre`),
            // 拿它当"互不相同"的判据会判错。
            "环闭合·平面四方",
            [
                ("[Pt@SP1]1(N)CO1", "[NH2][Pt@SP2]1[CH2][O]1"),
                ("[Pt@SP2]1(N)CO1", "[NH2][Pt@SP3]1[CH2][O]1"),
                ("[Pt@SP3]1(N)CO1", "[NH2][Pt@SP1]1[CH2][O]1"),
            ],
        ),
        (
            "环闭合·八面体",
            [
                (
                    "[Co@OH1]1(N)(P)(S)CO1",
                    "[NH2][Co@OH5]1([PH2])([SH])[CH2][O]1",
                ),
                (
                    "[Co@OH11]1(N)(P)(S)CO1",
                    "[NH2][Co@OH24]1([PH2])([SH])[CH2][O]1",
                ),
                (
                    "[Co@OH20]1(N)(P)(S)CO1",
                    "[NH2][Co@OH29]1([PH2])([SH])[CH2][O]1",
                ),
            ],
        ),
    ] {
        let mut seen = std::collections::BTreeSet::new();
        for (ours, reference) in cases {
            let a = canonical_smiles(&perceived(ours));
            assert_eq!(
                a,
                canonical_smiles(&perceived(reference)),
                "{family}:{ours} 与外部实现写的 {reference} 是同一个分子,规范串却不同"
            );
            seen.insert(a);
        }
        // 少了这一条,一个把序号整个丢光的实现也能过上面那些 —— 三串会塌成一串,
        // 而参照给的三串同样塌成一串。
        assert_eq!(seen.len(), 3, "{family}:三种排法塌成了 {} 串", seen.len());
    }
}

/// 三角双锥与八面体也一样:环闭合落在中心上时,**书写序 ≠ 存储序**。
///
/// 与上一条同源。第二种写法同样由外部实现给出 —— 语料里只有一个这种形态的
/// 分子(`[Co@OH5]1(N)(O)(S)(P)CCC1`,由 `differential_l3.rs` 覆盖),
/// 三角双锥一个都没有。
#[test]
fn coordination_stereo_with_a_ring_closure_matches_the_reference() {
    for (ours, reference) in [
        // 三角双锥
        ("[P@TB1]1(F)(Cl)(Br)CC1", "F[P@TB9]1(Cl)(Br)CC1"),
        ("[P@TB5]1(F)(Cl)(Br)CC1", "F[P@TB13]1(Cl)(Br)CC1"),
        ("[P@TB12]1(F)(Cl)(Br)CC1", "F[P@TB4]1(Cl)(Br)CC1"),
        // 八面体
        (
            "[Co@OH11]1(N)(O)(S)(P)CCC1",
            "[NH2][Co@OH21]1([OH])([PH2])([SH])[CH2]C[CH2]1",
        ),
        (
            "[Co@OH20]1(N)(O)(S)(P)CCC1",
            "[NH2][Co@OH8]1([OH])([PH2])([SH])[CH2]C[CH2]1",
        ),
    ] {
        assert_eq!(
            canonical_smiles(&perceived(ours)),
            canonical_smiles(&perceived(reference)),
            "{ours} 与外部实现写的 {reference} 是同一个分子,规范串却不同"
        );
    }
    // 三个三角双锥的排法互不相同 —— 少了这一条,一个把序号丢光的实现也能过上面那些
    let a = canonical_smiles(&perceived("[P@TB1]1(F)(Cl)(Br)CC1"));
    let b = canonical_smiles(&perceived("[P@TB5]1(F)(Cl)(Br)CC1"));
    let c = canonical_smiles(&perceived("[P@TB12]1(F)(Cl)(Br)CC1"));
    assert_ne!(a, b);
    assert_ne!(b, c);
    assert_ne!(a, c);
}

/// 不同的分子不能塌到同一个规范串。
///
/// 只测"相同的要相同"是不够的 —— 一个恒返回空串的实现也能通过。
#[test]
fn different_molecules_get_different_canonical_forms() {
    let smis = [
        "CCO",
        "CCN",
        "CCC",
        "C1CC1",
        "c1ccccc1",
        "C1CCCCC1",
        "OC(=O)c1ccccc1N",
        "OC(=O)c1ccccc1O",
        "CC(C)C",
        "CCCC",
        "C1CC2CCC1CC2",
        "C1CC2CCC(C1)CC2",
        // 平面四方:两个几何异构体只差一个排列序号,图完全相同。
        // 写不出 `@SP` 的时候这两条塌成一串 —— 那正是补上写出之前的状态。
        "[Pt@SP1](Cl)(Cl)(N)N",
        "[Pt@SP2](Cl)(Cl)(N)N",
    ];
    let mut seen: std::collections::HashMap<String, &str> = std::collections::HashMap::new();
    for smi in smis {
        let m = smiles::parse(smi).unwrap();
        let c = canonical_smiles(&m);
        if let Some(prev) = seen.insert(c.clone(), smi) {
            panic!("{smi} 与 {prev} 是不同的分子,却都规范成了 {c}");
        }
    }
}

/// 规范串必须是**不动点**:读回来再规范化,还是同一串。
///
/// 这条与重排恒定不是一回事,重排判据换的是分子对象的编号、不经过"写出再解析"
/// 那一趟,而那一趟会同时改掉原子次序、键序、氢的表示,以及**方向符号落在哪根
/// 键上**。最后这一项正是这条判据要守的。
///
/// # 为什么会不成立
///
/// 顺反记的是"相对某两个参照原子"。同一根双键换一个参照、把顺反翻一次,说的
/// 是同一件事;可写出器是**按参照键**放方向符号的,换个参照符号就落到另一根
/// 键上。而感知挑参照用的是"存储顺序里第一个带方向的邻居" —— 输入写法留下的
/// 痕迹。于是同一个分子换种写法读进来,规范串就差一个不携带信息的方向符号。
///
/// 下面两条取自语料(`harness/corpus/large.smi`),是全语料 8831 条里仅有的
/// 两条:双键一端挂着两个取代基,而其中一个又是另一根双键的端点,那个端点
/// 于是带了两根有方向的键,读回去时感知只认下其中一根。
#[test]
fn canonical_smiles_is_a_fixed_point() {
    for smi in [
        "CC1=C/C(=N\\O)/C(=N\\O)/N=C1",
        "CN1CCC\\2=C1/C(=N\\O)/S/C2=N\\c3ccc(cc3)F",
        // 对照:普通顺反、共轭链、环上双键 —— 本来就该成立,防止判据只盖到特例
        "F/C=C/F",
        "OC(=O)/C=C/c1ccccc1",
        "F/C=C/C=C/F",
    ] {
        // **必须走调用方真正走的那条路**:净化之后再感知顺反。少了感知,
        // `stereo_atoms` 根本没被填,参照原子的挑法压根不参与,判据当场空过 ——
        // 实测只 parse 的话撤掉修复也是绿的。
        let once = canonical_smiles(&perceived(smi));
        let twice = canonical_smiles(&perceived(&once));
        assert_eq!(once, twice, "{smi} 的规范串不是不动点");
    }
}
