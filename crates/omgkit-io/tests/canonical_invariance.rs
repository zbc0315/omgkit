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

    // 逐个手性原子把标记搬到新的存储序上
    for new_a in 0..out.num_atoms() as u32 {
        let old_a = atom_perm[new_a as usize];
        let tag = mol.atoms()[old_a as usize].chiral_tag;
        if !tag.is_tetrahedral() {
            continue;
        }
        let old_seq: Vec<usize> = mol.neighbors(old_a).map(|(_, b)| b as usize).collect();
        let new_seq: Vec<usize> = out
            .neighbors(new_a)
            .map(|(_, b)| new_to_old_bond[b as usize])
            .collect();
        if permutation_is_odd(&old_seq, &new_seq).expect("同一组键") {
            out.atom_mut(new_a).expect("原子存在").chiral_tag = tag.inverted();
        }
    }
    out
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
#[ignore = "语料大,用 cargo test -- --ignored 运行"]
fn canonical_smiles_is_invariant_under_renumbering_corpus() {
    let (stats, bad) = check_corpus(&corpus("large.smi"), 5, 0x243F_6A88_85A3_08D3);
    assert!(bad.is_empty(), "{}", report(&bad, 15));
    assert!(stats.molecules > 8000, "只测到 {} 条分子", stats.molecules);
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
