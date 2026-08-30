//! 按反应模板生成产物。
//!
//! # 映射号定义了一切
//!
//! 反应物模板与产物模板之间**只有映射号**这一条纽带:
//!
//! | 原子在哪 | 有映射号 | 无映射号 |
//! |---|---|---|
//! | 只在反应物模板 | 不可能(号在两侧都要找得到才叫"有") | **删掉**匹配到的那个原子 |
//! | 两侧都有 | **保留**匹配到的原子,按产物模板改写属性 |  |
//! | 只在产物模板 | 同上,视作新建 | **新建**一个原子 |
//!
//! 键同理:产物模板里有的键就建,反应物模板里有而产物模板里没有的键就断。
//!
//! # 产物分子数由连通性决定
//!
//! 产物模板描述的是反应中心的片段,**不是**一个片段一个分子。所有产物模板建进
//! 同一张图,模板之外的原子只搬一次,最后按连通分量切开。
//!
//! 逐产物各建一张图会出大问题:两个片段的锚点若仍通过未匹配的原子相连
//! (分子内反应、断环键),那批原子会被复制进每一个产物,**质量当场不守恒**
//! 而没有任何东西报错。
//!
//! # 模板之外的部分原样带过来
//!
//! 反应物分子里没被模板匹配到的原子和键要**照原样搬进产物**。这一条决定了
//! 反应模板可以写得很小 —— `[C:1][OH:2]>>[C:1][Cl:2]` 只描述羟基变氯,
//! 分子其余部分自动跟着走。
//!
//! # 产物不做净化
//!
//! 产物是"按模板改写出来的图",可能价键不合法(模板本身就能写出不合法的东西)。
//! 净化与否交给调用方 —— 在这里净化会把"模板写错了"和"这条反应本来就不该
//! 用在这个底物上"混成同一种失败。
//!
//! # 原子映射号是可选产出
//!
//! 模板里的映射号只在**模板内部**成立,它连的是"反应物模板的这个原子"与
//! "产物模板的那个原子",与底物无关。底物上真正的原子对应关系是运行时才
//! 产生的:模板匹配定下一部分,模板之外原样搬运的部分定下其余的。
//!
//! [`run_reactants`] 的 `atom_mapping` 参数控制要不要把这份运行时对应关系
//! 固化成映射号。开启时返回的 [`Outcome`] 里反应物与产物**两侧都带号**,
//! 同一个号出现在两侧就表示是同一个原子。

use std::collections::{BTreeMap, BTreeSet};

use omgkit_core::{
    AtomData, BondData, BondDirection, BondOrder, BondStereo, ChiralTag, MolBuilder,
};
use omgkit_io::smarts::{
    map_number, required_chirality, AtomExpr, AtomPrim, BondExpr, BondPrim, QueryMol, Reaction,
};

use crate::matcher::{substructure_matches, MatchOptions};
use crate::props::MolProps;

/// 一次反应的产物。
///
/// **长度不等于产物模板数。** 产物模板描述的是反应中心的片段;这些片段在底物里
/// 断没断开,要看模板之外的原子还连不连着。分子内环化的逆向就是这样:模板写成
/// 两个片段,可打开一个环并不产生两个分子。
pub type ProductSet = Vec<MolBuilder>;

/// 一个匹配组合跑出来的结果。
#[derive(Debug, Clone)]
pub struct Outcome {
    /// 产物。数目由**连通性**决定,不由产物模板数决定,见 [`ProductSet`]。
    pub products: ProductSet,
    /// 带映射号的反应物副本。
    ///
    /// **只在开启 `atom_mapping` 时非空。** 关掉时反应物一个字节都不会被改动,
    /// 调用方手上的输入就是答案,复制一份纯属浪费 —— 而一次反应可以产出上百组
    /// 结果,每组都复制一遍反应物不是小开销。
    pub reactants: Vec<MolBuilder>,
    /// `discarded[i]` = 第 i 个输入分子里**没有进入任何产物**的原子下标,升序。
    ///
    /// 这是一条**事实记录,不含任何推断**:模板明说要删的原子、以及只挂在它们
    /// 身上因而失去落脚点的原子,都在这里。产物侧看不到它们,所以不记的话这批
    /// 原子就是凭空消失 —— 而"消失"与"被判定为不该存在"是两回事。
    ///
    /// 把它们收口成分子是**另一件事**,由 [`crate::byproduct::reconstruct`] 做,
    /// 而且是推断:模板里没有"离去基团变成了什么"这条信息。两者分开是有意的 ——
    /// 本字段永远可信,那边的结论要看它自己给出的档次。
    pub discarded: Vec<Vec<u32>>,
}

/// 对一组反应物跑反应,返回所有结果组。
///
/// 每个反应物模板配一个**互不相同**的输入分子;分子数与模板数不等时返回空
/// (那一档是 [`run_on_substrate`] 的形状:多个片段落在同一个分子上)。
///
/// 每个匹配组合产出一组结果 —— 底物上有几处能反应就有几组,内容可能重复
/// (对称位点)。去重是调用方的事:要按什么去重取决于用途,规范 SMILES
/// 多重集只是其中一种。
///
/// # 递入顺序不影响出不出产物
///
/// **位置不是化学。** 谁先谁后是调用方敲键盘的顺序,不是分子的性质,所以
/// 它不该决定这条反应跑不跑得起来。
///
/// 实现上仍然**先试恒等分配**(第 i 个模板配第 i 个分子)—— 顺序本来就对得上
/// 时开销与只试这一种完全相同;恒等分配一个产物都给不出,才去找别的一一对应,
/// 按字典序取第一个能出产物的。所以:
///
/// - 顺序对得上:行为与耗时都不变
/// - 顺序不对:照样出产物,而不是交白卷
/// - **返回空只剩一个意思:这批分子上没有反应位点**
///
/// 这一条是量出来的,不是想出来的。USPTO-50k 正向语料里,按记录自带的分子
/// 顺序直接调用,约 **689 条**交白卷;抽样 4000 条逐条核过,其中 **59 条全部**
/// 只是顺序对不上 —— 换个顺序就出产物,**没有一条**是真的匹配不上。
/// 而调用方拿到的是同一个空列表,分不出这两件事。
///
/// 回退那条路要多算最多 n²−n 次子结构搜索(n 是反应物模板数,现实中 1–3),
/// 而且**只在本来就要返回空的时候才走** —— 拿它换的是一个静默的错答案。
///
/// 开了 `atom_mapping` 时,哪个分子担了哪个角色可以从映射号读回来:
/// [`Outcome::reactants`] 里的副本按**输入顺序**排,号是按底物原子发的。
///
/// # 哪些原子会进产物:只有从保留下来的原子**走得到**的
///
/// 产物 = 模板产物侧建出来的原子,加上从它们出发在底物里能走到的部分。走不到
/// 的原子不进产物,**不报错**。这条规则有两个看得见的后果,都是刻意的:
///
/// **一、模板删掉一个原子时,只挂在它身上的东西跟着走。** 叔丁酯水解的模板写
/// `C-C-[O:1]-[C:2]=[O:3]`,删掉的 `C-C` 是叔丁基的一个甲基加季碳;季碳上另外
/// 两个甲基没有别的路连回保留部分,于是一并消失 —— 一次少掉 4 个重原子而不是 2 个。
/// 真实反应语料里这一档数以千计,不是边角情形。
///
/// **二、完全不连通的旁观组分原样交回来,不丢。** 底物写成 `内酰胺.HCl` 或
/// `[Na+].[O-]CC(=O)OCC[O-].[K+]` 时,反离子与任何匹配到的原子都不连通,遍历
/// 走不到它们。**但走不到不等于该丢** —— 丢了产物的重原子数就少于底物,是引擎
/// 自己在破坏质量守恒,而且不报错。逆合成正是把模板作用到任意分子上,盐是常态。
///
/// 所以这些组分会被原样搬进产物图,按连通分量切开之后各自成为一个产物分子。
/// **引擎不替调用方决定归属** —— "这个反离子该跟哪一半走"没有普遍答案,模板里
/// 也没有这条信息;交回去,由调用方定。
///
/// 与上一条的分界是"这个组分里**有没有**原子被模板匹配到":有,留下还是删掉
/// 是模板的表态;没有,模板压根没提到它。
///
/// # `atom_mapping`
///
/// 开启后,[`Outcome::reactants`] 填上带映射号的反应物副本,产物侧的对应原子
/// 打上同一个号 —— 两侧合起来就是一条完整的原子映射反应。关闭时两侧都不带号,
/// `reactants` 留空。
///
/// 发号的规则:
///
/// - **只给两侧都在的原子发。** 被反应删掉的、产物侧新建的都不发 —— 一个在
///   另一侧找不到的号,读的人只能理解成"这个原子凭空消失/出现",而那正是号
///   要表达的反面。
/// - **号是新发的,不沿用模板里的。** 模板的 `[C:1]` 连的是两个模板而非底物,
///   而且模板只覆盖分子的一小块,搬运过来的部分本来就没有号可沿用。反应物
///   副本上原有的号会先清掉,免得留下在产物侧找不到的悬空号。
/// - **顺序**:按 `(第几个反应物, 原子下标)` 升序,从 1 连续发。同一份输入
///   永远得到同一套号;换一种原子编号写同一个分子,号会跟着变 —— 映射号本就
///   是贴在某一种画法上的标签,不是分子的不变量。
/// - 一个号在同一侧只出现一次。产物是按连通分量切出来的,每个底物原子只进
///   一个产物,所以这一条自然成立。
///
/// ## 写出带号的产物之前要先净化
///
/// 与"产物不做净化"那一条配套:**隐式氢数是派生量,图改完就过期了**。模板删掉
/// 一个邻居之后,那个原子该补几个氢要重算,而重算在净化里。
///
/// 平时看不出来 —— 裸写形式(`N`、`C`)把氢数交给读的人按价规则去推,推出来
/// 的正是重算后的值。可**带映射号的原子必须写进方括号**,而方括号里的氢数是
/// 显式的,写出去的就是那个过期的缓存值。于是同一个产物,开不开映射号会写出
/// **氢数不同**的两串。
///
/// 所以带号写出之前先 [`sanitize`](omgkit_chem::sanitize)。实测语料上:不净化
/// 直接写,31 万个 outcome 里有 176 个两串对不上;先净化再写,**0 个**。
///
/// # 调用前要先感知双键顺反
///
/// 反应物应当先跑
/// [`perceive_bond_stereo`](omgkit_io::stereo::perceive_bond_stereo) ——
/// 净化那 12 步里**没有**它(感知要用对称等价类,那在净化的上一层,调不到)。
///
/// 漏了这一步不会报错,会**静默丢几何**:方向键(`/` `\`)依附在某根单键上,
/// 反应把那根键删掉,几何就跟着没了 —— 哪怕双键本身根本没被碰过。产物照样
/// 合法、原子数照样对,只有顺反悄悄少了。感知一次之后信息记在双键自己身上,
/// 只要参照原子还在就活得下来。
///
/// ```no_run
/// # use omgkit_core::MolBuilder;
/// # use omgkit_match::{run_reactants, MolProps};
/// # fn demo(mut mol: MolBuilder, rxn: &omgkit_io::smarts::Reaction) {
/// omgkit_chem::sanitize(&mut mol).unwrap();
/// omgkit_io::stereo::perceive_bond_stereo(&mut mol); // ← 别漏
/// let props = MolProps::compute(&mol);
/// let out = run_reactants(rxn, &[(mol, props)], 0, false);
/// # let _ = out;
/// # }
/// ```
///
/// debug 构建下漏了会被 [`debug_assert`] 当场拦住;release 下不做这个检查。
#[must_use]
pub fn run_reactants(
    reaction: &Reaction,
    reactants: &[(MolBuilder, MolProps)],
    max_products: usize,
    atom_mapping: bool,
) -> Vec<Outcome> {
    debug_assert!(
        !reactants
            .iter()
            .any(|(m, _)| omgkit_io::stereo::directions_not_perceived(m)),
        "反应物里有双键的几何**方向键已经写明**、却没有感知过顺反 —— \
         漏了 omgkit_io::stereo::perceive_bond_stereo。这样跑不会报错,\
         但反应一旦删掉承载方向的那根单键,几何会静默丢失"
    );
    if reactants.len() != reaction.reactants.len() || reaction.products.is_empty() {
        return Vec::new();
    }

    // 逐个反应物模板找匹配,再取笛卡尔积
    let opts = MatchOptions {
        max_matches: 0,
        uniquify: false,
        // 反应侧**不判**立体。
        //
        // 反应模板是跨工具流通的东西,读得更严会让现成的模板不再出产物,
        // 而"少了产物"比"多了产物"难发现得多。子结构匹配那边默认判,
        // 因为那里作者写什么就该算什么。
        use_chirality: false,
    };
    let n = reaction.reactants.len();
    let per_template: Vec<Vec<Vec<u32>>> = reaction
        .reactants
        .iter()
        .zip(reactants)
        .map(|(t, (mol, props))| substructure_matches(t, mol, props, opts))
        .collect();
    // 连通分量按分子算一次就够 —— 它与匹配到哪个位点无关,更与分配无关
    let comps: Vec<Vec<u32>> = reactants.iter().map(|(m, _)| components(m)).collect();

    // **恒等分配先试**:第 i 个模板配第 i 个分子。它只要 n 次子结构搜索,
    // 所以顺序本来就对得上时,这条路的开销与先前一模一样。
    if per_template.iter().all(|m| !m.is_empty()) {
        let identity: Vec<usize> = (0..n).collect();
        let out = outcomes_under(
            reaction,
            reactants,
            &per_template,
            &identity,
            &comps,
            max_products,
            atom_mapping,
        );
        if !out.is_empty() {
            return out;
        }
    }

    // 恒等分配一个产物都给不出 —— 再看**别的一一对应**行不行。
    // 到这里才把匹配表补满(最多再 n²−n 次搜索),恒等那条路一分钱不多花。
    let mut table: Vec<Vec<Vec<Vec<u32>>>> = Vec::with_capacity(n);
    for (t, tpl) in reaction.reactants.iter().enumerate() {
        let mut row = Vec::with_capacity(n);
        for (m, (mol, props)) in reactants.iter().enumerate() {
            row.push(if m == t {
                per_template[t].clone()
            } else {
                substructure_matches(tpl, mol, props, opts)
            });
        }
        table.push(row);
    }
    let mut assign = vec![0usize; n];
    let mut used = vec![false; n];
    search_assignment(
        reaction,
        reactants,
        &table,
        &comps,
        max_products,
        atom_mapping,
        0,
        &mut assign,
        &mut used,
    )
    .unwrap_or_default()
}

/// 在一个**确定的分配**下枚举匹配的笛卡尔积、造产物。
///
/// `assign[i]` 是第 i 个反应物模板落在第几个输入分子上;`per_template[i]` 是
/// 那个模板在**那个分子**里的全部匹配。
fn outcomes_under(
    reaction: &Reaction,
    reactants: &[(MolBuilder, MolProps)],
    per_template: &[Vec<Vec<u32>>],
    assign: &[usize],
    comps: &[Vec<u32>],
    max_products: usize,
    atom_mapping: bool,
) -> Vec<Outcome> {
    let mut out = Vec::new();
    let mut combo: Vec<usize> = vec![0; per_template.len()];
    loop {
        let mapping: Vec<&Vec<u32>> = combo
            .iter()
            .enumerate()
            .map(|(i, &j)| &per_template[i][j])
            .collect();
        let built = build_products(reaction, reactants, &mapping, assign, comps);
        out.push(stamp_atom_maps(reactants, built, atom_mapping));
        if max_products != 0 && out.len() >= max_products {
            return out;
        }
        // 进位
        let mut i = 0;
        loop {
            if i == combo.len() {
                return out;
            }
            combo[i] += 1;
            if combo[i] < per_template[i].len() {
                break;
            }
            combo[i] = 0;
            i += 1;
        }
    }
}

/// 按字典序找**第一个能出产物的一一对应**。
///
/// 搜索的是"模板 ↔ 分子"的完美匹配,匹配表已经算好,所以每一步只是查表 ——
/// 某个模板在剩下的分子里一个都匹配不上时当场剪掉,不往下展。
#[allow(clippy::too_many_arguments)]
fn search_assignment(
    reaction: &Reaction,
    reactants: &[(MolBuilder, MolProps)],
    table: &[Vec<Vec<Vec<u32>>>],
    comps: &[Vec<u32>],
    max_products: usize,
    atom_mapping: bool,
    depth: usize,
    assign: &mut Vec<usize>,
    used: &mut Vec<bool>,
) -> Option<Vec<Outcome>> {
    if depth == assign.len() {
        let per: Vec<Vec<Vec<u32>>> = assign
            .iter()
            .enumerate()
            .map(|(t, &m)| table[t][m].clone())
            .collect();
        let out = outcomes_under(
            reaction,
            reactants,
            &per,
            assign,
            comps,
            max_products,
            atom_mapping,
        );
        return if out.is_empty() { None } else { Some(out) };
    }
    for m in 0..used.len() {
        if used[m] || table[depth][m].is_empty() {
            continue;
        }
        used[m] = true;
        assign[depth] = m;
        if let Some(out) = search_assignment(
            reaction,
            reactants,
            table,
            comps,
            max_products,
            atom_mapping,
            depth + 1,
            assign,
            used,
        ) {
            return Some(out);
        }
        used[m] = false;
    }
    None
}

/// 把若干个分子拼成一张图(不加任何键),返回拼好的图。
///
/// 只有 [`BondData`] 的 `begin`、`end`、`stereo_atoms` 带原子下标,要加偏移;
/// [`AtomData`] 一个下标都不带 —— 手性是**相对邻居顺序**说的,不是相对下标。
/// 正因如此,原子与键都必须**按原顺序**逐个搬:顺序一乱,每个原子的邻居序
/// 跟着乱,手性标记的含义就变了,而拓扑、原子数、电荷全对,只有构型悄悄反了。
fn concat(mols: &[(MolBuilder, MolProps)]) -> MolBuilder {
    let n_atoms = mols.iter().map(|(m, _)| m.num_atoms()).sum();
    let n_bonds = mols.iter().map(|(m, _)| m.num_bonds()).sum();
    let mut out = MolBuilder::with_capacity(n_atoms, n_bonds);
    for (m, _) in mols {
        let base = u32::try_from(out.num_atoms()).unwrap_or(u32::MAX);
        for a in m.atoms() {
            out.add_atom_data(*a);
        }
        for b in m.bonds() {
            let mut nb = *b;
            nb.begin += base;
            nb.end += base;
            for s in &mut nb.stereo_atoms {
                if *s != BondData::NO_STEREO_ATOM {
                    *s += base;
                }
            }
            let _ = out.add_bond_data(nb);
        }
    }
    out
}

/// 把整个反应物侧当作**一张图**上的查询来跑,而不是按位置配对。
///
/// # 与 [`run_reactants`] 的分工
///
/// [`run_reactants`] 的契约是"N 个反应物模板 ↔ N 个输入分子,**一一对应**" ——
/// 先恒等分配,给不出产物再搜别的对应关系,所以递入顺序不决定它跑不跑得起来。
/// 但它仍然是**一对一**的:模板的片段数比分子数多时直接交白卷,而那正是
/// **分子内反应**的形状 —— 两个片段落在同一个分子上。
///
/// 本函数把输入拼成一张图,让每个反应物模板在整张图上自由找位置,只要求
/// 各模板匹配到的原子**两两不重叠**。于是
///
/// - 分子间:片段落在不同的连通分量上,与位置式的结果一致(不必再枚举排列)
/// - 分子内:片段落在同一个分量的不同部位 —— 位置式表达不了的那一档
/// - 盐:阳离子与阴离子是同一个分子的两个组分,模板可以同时碰到它们
///
/// 这与产物侧是同一条原则:**片数是(模板, 底物)共同的性质,不是模板的性质**。
/// 产物侧早就这么做了(建进同一张图、按连通分量切开),这里只是把同一条原则
/// 补到反应物侧。
///
/// # 代价
///
/// 匹配的搜索空间变大:位置式下第 i 个模板只在第 i 个分子里找,这里在整张图里
/// 找,再靠不相交筛掉大部分组合。片段多、分子大时组合数会涨得很快,`max_products`
/// 只截输出、不截枚举。要可预测的耗时就用 [`run_reactants`]。
///
/// # 调用前同样要先感知双键顺反
///
/// 理由与 [`run_reactants`] 完全相同,见那里。
#[must_use]
pub fn run_on_substrate(
    reaction: &Reaction,
    substrate: &[(MolBuilder, MolProps)],
    max_products: usize,
    atom_mapping: bool,
) -> Vec<Outcome> {
    debug_assert!(
        !substrate
            .iter()
            .any(|(m, _)| omgkit_io::stereo::directions_not_perceived(m)),
        "底物里有双键的几何**方向键已经写明**、却没有感知过顺反 —— \
         漏了 omgkit_io::stereo::perceive_bond_stereo。理由见 run_reactants"
    );
    if substrate.is_empty() || reaction.reactants.is_empty() || reaction.products.is_empty() {
        return Vec::new();
    }

    // 拼图之前先记下各分子的原子数 —— `discarded` 要按这个切回去,见
    // `regroup_discarded`
    let sizes: Vec<usize> = substrate.iter().map(|(m, _)| m.num_atoms()).collect();
    let mol = concat(substrate);
    let props = MolProps::compute(&mol);
    let inputs = [(mol, props)];

    // 反应侧不判立体,理由见 `run_reactants`
    let opts = MatchOptions {
        max_matches: 0,
        uniquify: false,
        use_chirality: false,
    };
    let per_template: Vec<Vec<Vec<u32>>> = reaction
        .reactants
        .iter()
        .map(|t| substructure_matches(t, &inputs[0].0, &inputs[0].1, opts))
        .collect();
    if per_template.iter().any(Vec::is_empty) {
        return Vec::new();
    }

    // 所有模板都落在这唯一一张图上
    let home = vec![0usize; reaction.reactants.len()];
    let n_atoms = inputs[0].0.num_atoms();
    let comps: Vec<Vec<u32>> = inputs.iter().map(|(m, _)| components(m)).collect();

    let mut out = Vec::new();
    let mut combo: Vec<usize> = vec![0; per_template.len()];
    let mut used = vec![false; n_atoms];
    loop {
        let mapping: Vec<&Vec<u32>> = combo
            .iter()
            .enumerate()
            .map(|(i, &j)| &per_template[i][j])
            .collect();
        // 两个模板抢同一个原子是不合法的:位置式契约靠"分子各不相同"天然
        // 保证了这一点,拼成一张图之后必须自己判。
        used.iter_mut().for_each(|u| *u = false);
        let disjoint = mapping.iter().all(|m| {
            m.iter().all(|&a| {
                let fresh = !used[a as usize];
                used[a as usize] = true;
                fresh
            })
        });
        if disjoint {
            let built = build_products(reaction, &inputs, &mapping, &home, &comps);
            let mut outcome = stamp_atom_maps(&inputs, built, atom_mapping);
            outcome.discarded = regroup_discarded(&outcome.discarded, &sizes);
            out.push(outcome);
            if max_products != 0 && out.len() >= max_products {
                break;
            }
        }
        // 进位
        let mut i = 0;
        loop {
            if i == combo.len() {
                return out;
            }
            combo[i] += 1;
            if combo[i] < per_template[i].len() {
                break;
            }
            combo[i] = 0;
            i += 1;
        }
    }
    out
}

/// 反应物侧一个映射号对应的**具体原子**:(第几个反应物, 原子下标)。
type Anchor = (usize, u32);

/// 反应物侧提前算好的三张表,产物构建全程只读。
struct ReactantFacts {
    /// 映射号 → 反应物里的那个原子
    anchors: BTreeMap<u16, Anchor>,
    /// 映射号 → 该原子在**反应物模板**里的度数。产物侧要拿它比,判断
    /// "这个原子的连接有没有变" —— 变了的话氢数要重算,不能照抄。
    degree: BTreeMap<u16, usize>,
    /// 映射号 → **反应物模板**在这个原子上写的手性。产物侧要拿它比,
    /// 见 [`ChiralityPlan`]。
    chirality: BTreeMap<u16, Option<ChiralTag>>,
    /// 映射号 → 该原子在**反应物模板**里的邻居次序,按映射号记(没号的记
    /// `None`)。手性标记是相对邻居顺序的,两侧比标记之前要先比这个 ——
    /// 见 [`template_order_is_odd`]。
    neighbors: BTreeMap<u16, Vec<Option<u16>>>,
}

/// 一个模板原子的邻居次序,按映射号记。没有映射号的邻居记 `None`。
fn neighbor_maps(template: &QueryMol, qi: u32) -> Vec<Option<u16>> {
    template
        .topology
        .neighbors(qi)
        .map(|(other, _)| map_number(&template.atoms[other as usize]))
        .collect()
}

/// 反应物模板与产物模板在同一个原子上的邻居次序,置换是不是奇的。
///
/// 手性标记说的是"按**这张模板自己**的邻居顺序看过去"的构型。两侧的顺序不同时,
/// 同一个 `@` 说的是两种构型 —— 所以"两侧写得一样不一样"要在同一个顺序下问,
/// 不能直接比标记。典型形状:
///
/// ```text
/// [F:2][C@H:1]([Cl:3])[Br:4]>>[Br:4][C@H:1]([Cl:3])[F:2]
/// ```
///
/// 两侧都是 `@`,取代基一个没换,可产物侧把 F 与 Br 对调着写 —— 置换是奇的,
/// 所以这条模板说的是**翻转**,不是保留。
///
/// # 只容忍每侧一个对不上的邻居
///
/// 一侧有、另一侧没有的邻居各至多一个时,它们互相顶替,对应关系唯一。多于一个
/// 就不唯一了(两种配对宇称相反),这时返回 `None` 表示**不调整** —— 猜一个
/// 只会把未定义变成另一个未定义。
fn template_order_is_odd(react: &[Option<u16>], prod: &[Option<u16>]) -> Option<bool> {
    // 度数不足 3 的中心谈不上四面体手性;两侧差出一个以上也没法对应
    if react.len() < 3 || prod.len() < 3 || react.len().abs_diff(prod.len()) > 1 {
        return None;
    }
    // 短的那侧补一个空位,代表"对面多出来的那个邻居"
    let mut r: Vec<Option<u16>> = react.to_vec();
    let mut p: Vec<Option<u16>> = prod.to_vec();
    if r.len() < p.len() {
        r.push(None);
    } else if p.len() < r.len() {
        p.push(None);
    }
    if r.iter().filter(|x| x.is_none()).count() > 1 || p.iter().filter(|x| x.is_none()).count() > 1
    {
        return None;
    }
    fill_missing(&mut r, &p)?;
    fill_missing(&mut p, &r)?;
    let enc = |v: &[Option<u16>]| -> Vec<u32> {
        v.iter().map(|x| x.map_or(u32::MAX, u32::from)).collect()
    };
    omgkit_core::permutation_is_odd(&enc(&r), &enc(&p))
}

/// `want` 里有而 `have` 里没有的映射号,填进 `have` 唯一的那个空位。
///
/// 空位不够用就说明两侧的邻居对应不上,返回 `None`。
fn fill_missing(have: &mut [Option<u16>], want: &[Option<u16>]) -> Option<()> {
    for &elem in want.iter().flatten() {
        if have.contains(&Some(elem)) {
            continue;
        }
        let slot = have.iter().position(Option::is_none)?;
        have[slot] = Some(elem);
    }
    Some(())
}

/// 产物这个原子的手性该怎么定。由**模板两侧写没写**决定,四种组合各有各的含义。
///
/// 模板作者用"写不写手性"表达反应对构型做了什么,这是一套约定,不是可以自由
/// 发挥的地方:
///
/// | 反应物侧 | 产物侧 | 含义 |
/// |---|---|---|
/// | 没写 | 没写 | 模板没管这件事 —— 底物的构型原样带过来 |
/// | 写了 | 没写 | 构型被破坏 —— 清掉 |
/// | 没写 | 写了 | 构型是**新建**的 —— 照模板写死 |
/// | 写了 | 写了 | 相对底物**保留**(两个标记相同)或**翻转**(不同) |
///
/// 最后一行最容易做错:两侧都写时,产物侧那个标记不是"要建成这个构型",而是
/// "与反应物侧那个标记比一比"。照字面写死的话,产物的构型就与底物无关了 ——
/// 同一个模板作用在一对对映体上会给出同一个产物,而正确答案是一对对映体。
#[derive(Clone, Copy, PartialEq, Eq)]
enum ChiralityPlan {
    /// 两侧都没写 —— 继承底物,之后由 [`rebase_chirality`] 换参照系
    Inherit,
    /// 只有反应物侧写了 —— 清掉
    Drop,
    /// 只有产物侧写了 —— 照模板写死
    Set,
    /// 两侧都写了且相同 —— 保留底物的构型
    Retain,
    /// 两侧都写了且不同 —— 底物的构型翻一次
    Invert,
}

impl ChiralityPlan {
    /// `order_is_odd` 是两侧模板邻居次序的置换宇称,见
    /// [`template_order_is_odd`]。次序对不上时给 `None`,那就只比标记。
    fn decide(
        reactant: Option<ChiralTag>,
        product: Option<ChiralTag>,
        order_is_odd: Option<bool>,
    ) -> Self {
        match (reactant, product) {
            (None, None) => Self::Inherit,
            (Some(_), None) => Self::Drop,
            (None, Some(_)) => Self::Set,
            // 标记一样、次序也一样 → 保留;两者恰好有一个反了 → 翻转。
            // 只比标记的话,产物侧把两个取代基对调着写就会被当成"保留"。
            (Some(r), Some(p)) => {
                if (r == p) != order_is_odd.unwrap_or(false) {
                    Self::Retain
                } else {
                    Self::Invert
                }
            }
        }
    }
}

/// 一个产物,连同"它的每个原子从哪个反应物原子来"。
///
/// 第二项按反应物下标分组,每组是 `反应物原子 → 产物原子`。产物侧新建的原子
/// 不在表里 —— 它们没有反应物出处。
type BuiltProduct = (MolBuilder, Vec<BTreeMap<u32, u32>>);

/// `home[ti]` = 第 ti 个反应物模板匹配进了第几个输入分子。
///
/// 拆出这个间接层是为了让"一个模板 ↔ 一个分子"不再是写死的假设:分子内反应
/// 里好几个模板落在同一个分子上,`home` 全指向同一个下标。位置式的
/// [`run_reactants`] 传的是恒等映射,行为一字不变。
/// 连通分量编号,每个原子一个。**每个分子只算一次**。
///
/// 分量结构是分子自己的性质,与匹配到哪个位点无关;而 `build_products` 每出一个
/// outcome 就要判一次"哪些分量没有被匹配到",一个底物上有几十处匹配就要重算几十遍。
/// 算一次传下去,省下的是常数乘以 outcome 数。
pub(crate) fn components(mol: &MolBuilder) -> Vec<u32> {
    let n = mol.num_atoms();
    let mut comp = vec![u32::MAX; n];
    let mut stack: Vec<u32> = Vec::new();
    let mut next = 0u32;
    for s in 0..n as u32 {
        if comp[s as usize] != u32::MAX {
            continue;
        }
        comp[s as usize] = next;
        stack.push(s);
        while let Some(a) = stack.pop() {
            for (other, _) in mol.neighbors(a) {
                if comp[other as usize] == u32::MAX {
                    comp[other as usize] = next;
                    stack.push(other);
                }
            }
        }
        next += 1;
    }
    comp
}

fn build_products(
    reaction: &Reaction,
    reactants: &[(MolBuilder, MolProps)],
    matches: &[&Vec<u32>],
    home: &[usize],
    comps: &[Vec<u32>],
) -> Vec<BuiltProduct> {
    // 被模板匹配到的原子(不论有没有映射号),这些不再作为"模板外的部分"搬运
    let mut matched: Vec<Vec<bool>> = reactants
        .iter()
        .map(|(m, _)| vec![false; m.num_atoms()])
        .collect();
    // 反应物模板**亲自匹配到**的那些底物键,按底物的键下标索引。
    //
    // 只有这些键归产物模板管;模板没看见的键它无权删,理由见 `carry_over`。
    let mut template_bonds: Vec<Vec<bool>> = reactants
        .iter()
        .map(|(m, _)| vec![false; m.num_bonds()])
        .collect();

    let mut facts = ReactantFacts {
        anchors: BTreeMap::new(),
        degree: BTreeMap::new(),
        chirality: BTreeMap::new(),
        neighbors: BTreeMap::new(),
    };

    for (ti, template) in reaction.reactants.iter().enumerate() {
        // 第 ti 个模板落在第几个输入分子上。位置式契约下就是 ti 自己;把两个
        // 模板放到同一个分子上(分子内反应)时,好几个 ti 会指向同一个下标。
        let ri = home[ti];
        for qb in template.topology.bonds() {
            let (a, b) = (matches[ti][qb.begin as usize], matches[ti][qb.end as usize]);
            if let Some(bi) = reactants[ri].0.bond_between(a, b) {
                template_bonds[ri][bi as usize] = true;
            }
        }
        for (qi, &target) in matches[ti].iter().enumerate() {
            matched[ri][target as usize] = true;
            if let Some(n) = map_number(&template.atoms[qi]) {
                facts.anchors.entry(n).or_insert((ri, target));
                facts
                    .degree
                    .entry(n)
                    .or_insert_with(|| template.topology.degree(qi as u32));
                facts
                    .chirality
                    .entry(n)
                    .or_insert_with(|| required_chirality(&template.atoms[qi]));
                facts
                    .neighbors
                    .entry(n)
                    .or_insert_with(|| neighbor_maps(template, qi as u32));
            }
        }
    }

    // 所有产物模板建进**同一张图**,未匹配的部分只搬一次,最后按连通分量切开。
    //
    // 逐产物各建一张图是错的:两个产物模板的锚点若仍通过未匹配的原子相连
    // (分子内反应、断环键都是这样),那批原子会被**复制**进每一个产物 ——
    // 质量当场不守恒。分子内环化的逆向尤其明显:打开一个环并不产生两个分子。
    let mut out = MolBuilder::new();
    let mut from_reactant: Vec<BTreeMap<u32, u32>> =
        reactants.iter().map(|_| BTreeMap::new()).collect();
    // 手性已经由模板定死、不再定基的那些原子,见 `rebase_chirality`。
    let mut settled_chirality: BTreeSet<u32> = BTreeSet::new();

    for pt in &reaction.products {
        emit_template(
            pt,
            reactants,
            &facts,
            &mut out,
            &mut from_reactant,
            &mut settled_chirality,
        );
    }

    // 模板之外的部分:每个原子只搬一次。
    //
    // 先埋旁观组分的种子,再走遍历 —— 两步共用同一趟 `carry_over`,键的排序
    // 纪律与重定基因此完全一致。
    for (ti, (mol, _)) in reactants.iter().enumerate() {
        seed_spectators(
            mol,
            &comps[ti],
            &matched[ti],
            &mut from_reactant[ti],
            &mut out,
        );
        carry_over(
            mol,
            &matched[ti],
            &template_bonds[ti],
            &mut from_reactant[ti],
            &mut out,
        );
    }
    // 手性与顺反都要在切分**之前**定基 —— 切分保持邻居的相对顺序,定基却要看全图
    for (ti, (mol, _)) in reactants.iter().enumerate() {
        rebase_chirality(mol, &from_reactant[ti], &settled_chirality, &mut out);
        rebase_bond_stereo(mol, &from_reactant[ti], &mut out);
    }

    split_components(&out, &from_reactant)
}

/// 模板里哪些方向键**真的表了态**,按模板键下标索引。
///
/// # 一根方向键单独出现时什么也没说
///
/// 顺反要靠双键**两端各一根**方向键才定得下来。`F/C=C/F` 是反式;而 `F/C=CF`
/// 里那根 `/` 定不了任何东西 —— 它只说了"这根单键画在右上",另一端没画,
/// 两个取代基的相对位置仍然未知。SMILES 与 SMARTS 在这一点上是同一套规则。
///
/// # 照抄一根孤立方向键会把参照系换掉
///
/// 孤立的那根抄进产物之后,并不会孤立地待着:双键另一侧的键多半是从底物
/// **继承**来的,它带着底物的方向。于是产物里凑出了一对方向 —— 一根来自
/// 模板的书写顺序,一根来自底物的书写顺序,两个互不相干的参照系。凑出来的
/// 几何是任意的:底物明明是反式,产物可以变成顺式,而拓扑、原子数、电荷
/// 全对,只有几何被悄悄换掉。这正是本项目最难发现的那一类错。
///
/// rdchiral 从 USPTO-50k 抽出的模板大量落在这一档:两侧各写一根孤立方向键,
/// 且书写朝向一正一反(反应物侧写 `[C:2]/[C:4]=`,产物侧写 `=[C:4]/[C:2]`),
/// 于是每一条都恰好翻一次。实测正向 100 条以上因此翻错。
///
/// # 判据
///
/// 双键"被模板定死"= 它两端**各自**还挂着至少一根写了方向的键。方向键"算数"
/// = 它挨着至少一根被定死的双键。两条都不满足时按没写方向处理,几何交回给
/// 继承那一支 —— 那一支两侧都取自同一个底物,参照系是一致的。
fn honoured_directions(template: &QueryMol) -> Vec<bool> {
    let bonds = template.topology.bonds();
    let has_dir: Vec<bool> = template
        .bonds
        .iter()
        .map(|e| bond_direction_from(e) != BondDirection::None)
        .collect();
    // 这个原子上,除了 `skip` 那根键之外还有没有写了方向的键
    let flanked = |atom: u32, skip: usize| {
        template
            .topology
            .neighbors(atom)
            .any(|(_, bi)| bi as usize != skip && has_dir[bi as usize])
    };
    let determined: Vec<bool> = (0..bonds.len())
        .map(|bi| {
            product_bond_from(&template.bonds[bi]) == ProductBond::Fixed(BondOrder::Double)
                && flanked(bonds[bi].begin, bi)
                && flanked(bonds[bi].end, bi)
        })
        .collect();
    (0..bonds.len())
        .map(|bi| {
            has_dir[bi]
                && [bonds[bi].begin, bonds[bi].end].iter().any(|&a| {
                    template
                        .topology
                        .neighbors(a)
                        .any(|(_, ob)| ob as usize != bi && determined[ob as usize])
                })
        })
        .collect()
}

/// 把一个产物模板的原子与键建进共享图。
fn emit_template(
    template: &QueryMol,
    reactants: &[(MolBuilder, MolProps)],
    facts: &ReactantFacts,
    out: &mut MolBuilder,
    from_reactant: &mut [BTreeMap<u32, u32>],
    settled_chirality: &mut BTreeSet<u32>,
) {
    let mut from_template: Vec<u32> = Vec::with_capacity(template.num_atoms());
    let mut anchor_of: Vec<Option<Anchor>> = Vec::with_capacity(template.num_atoms());

    // 一、产物模板里的原子
    for (qi, expr) in template.atoms.iter().enumerate() {
        let anchor = map_number(expr)
            .and_then(|n| facts.anchors.get(&n))
            .copied();
        anchor_of.push(anchor);
        let base = match anchor {
            // 有映射号且反应物侧找得到 —— 继承原子,再按模板改写
            Some((ti, ai)) => reactants[ti].0.atoms()[ai as usize],
            // 没有 —— 新建
            None => AtomData::new(0),
        };
        // 该原子在产物模板里的度数与它在反应物模板里的度数是否一致
        let degree_kept = map_number(expr)
            .and_then(|n| facts.degree.get(&n).copied())
            .is_some_and(|d| d == template.topology.degree(qi as u32));
        let plan = ChiralityPlan::decide(
            map_number(expr)
                .and_then(|n| facts.chirality.get(&n).copied())
                .flatten(),
            required_chirality(expr),
            map_number(expr)
                .and_then(|n| facts.neighbors.get(&n))
                .and_then(|r| template_order_is_odd(r, &neighbor_maps(template, qi as u32))),
        );
        let idx = out.add_atom_data(apply_template(base, expr, degree_kept, plan));
        // 模板在这个原子上表过态的,标记就此定死,不再定基。
        //
        // 定基换的是"反应物邻居序 → 产物邻居序",只对**继承来的**标记成立。
        // 模板一旦发话,产物侧那个标记说的就是产物自己参照系里的构型,再套一次
        // 反应物侧的置换等于凭空多翻一道。`Drop` 那档没有标记可言,也不必进来。
        if plan == ChiralityPlan::Set {
            settled_chirality.insert(idx);
        }
        from_template.push(idx);
        if let Some((ti, ai)) = anchor {
            from_reactant[ti].insert(ai, idx);
        }
    }

    // 二、产物模板里的键
    let honoured = honoured_directions(template);
    for (bi, expr) in template.bonds.iter().enumerate() {
        let b = template.topology.bonds()[bi];
        // 键级有三种来源,见 [`ProductBond`]。原子已经在上一段建完了,所以
        // "两端都芳香吗"此刻问得出来。
        let order = match product_bond_from(expr) {
            ProductBond::Fixed(o) => o,
            ProductBond::FollowAromaticity => {
                let aromatic = |ti: u32| {
                    out.atoms()[from_template[ti as usize] as usize]
                        .flags
                        .contains(omgkit_core::AtomFlags::AROMATIC)
                };
                if aromatic(b.begin) && aromatic(b.end) {
                    BondOrder::Aromatic
                } else {
                    BondOrder::Single
                }
            }
            ProductBond::Inherit => {
                match (anchor_of[b.begin as usize], anchor_of[b.end as usize]) {
                    (Some((t1, a1)), Some((t2, a2))) if t1 == t2 => {
                        inherited_order(&reactants[t1].0, a1, a2)
                    }
                    // 底物里没有这根键,模板又没说建什么 —— 谁都没表过态。
                    _ => BondOrder::Unspecified,
                }
            }
        };
        // 配位键的**朝向靠端点顺序表达**:`begin` 必须是给电子的一端。
        //
        // 而查询侧不是这么存的:`A<-B` 的端点按书写顺序记成 (A, B),由基元
        // `DativeReversed` 去区分朝向 —— 匹配时那样最直接。照搬端点建产物,
        // 箭头就反了:模板写 `>>[Fe:1]<-O(C)C`(氧给铁)会建成铁给氧,
        // 而"接受"的配位键计入受体的价,那个氧当场超价。
        let (tb, te) = if is_dative_reversed(expr) {
            (b.end, b.begin)
        } else {
            (b.begin, b.end)
        };
        let mut bd = BondData::new(
            from_template[tb as usize],
            from_template[te as usize],
            order,
        );
        // 芳香键的标志位与键级始终同步
        bd.flags.set(
            omgkit_core::BondFlags::AROMATIC,
            order == BondOrder::Aromatic,
        );
        // 方向有两个可能的来源,**模板真的表了态时**它优先。
        //
        // 一、模板自己写了 `/` 或 `\`,而且写成了**能定下几何的那种**:
        //    `>>C/[C:1]=[C:2]/C` 说的就是"新生成的双键是反式"。它相对模板键的
        //    begin → end,而新键的两端正是模板两端的像,朝向一致,直接照抄。
        //
        //    孤零零一根方向键不算表态,判据见 `honoured_directions` —— 它
        //    什么几何也没定,却会把另一侧继承来的方向拽进另一个参照系。
        //
        // 二、模板没写方向,但两端都源自同一个反应物、那里本来就有这根键。
        //    方向表达的是**旁边那根双键**两侧取代基的相对位置 —— 反应没碰那根
        //    双键的话,这个关系就该原样留着。丢了它,产物从确定的顺式(或反式)
        //    退化成未指定,是实打实的丢信息。
        //
        //    这一支的朝向要对齐:存的方向相对源键的 begin → end,而新键的 begin
        //    对应的是模板端点 b.begin 的出处,两者不一定同向。
        let from_template_dir = if honoured[bi] {
            bond_direction_from(expr)
        } else {
            BondDirection::None
        };
        bd.direction = if from_template_dir != BondDirection::None {
            from_template_dir
        } else if let (Some((t1, a1)), Some((t2, a2))) =
            (anchor_of[b.begin as usize], anchor_of[b.end as usize])
        {
            if t1 == t2 {
                inherited_direction(&reactants[t1].0, a1, a2)
            } else {
                BondDirection::None
            }
        } else {
            BondDirection::None
        };
        let _ = out.add_bond_data(bd);
    }
}

/// 把共享图按连通分量切成一个个产物分子。
///
/// # 产物分子数由**连通性**决定,不等于产物模板数
///
/// 产物模板描述的是反应中心的片段。片段之间断没断开,要看底物 —— 模板之外的
/// 原子可能仍把它们连着。分子内环化的逆向就是这样:模板写成两个片段,可打开
/// 一个环并不产生两个分子,原子一个不多不少。
///
/// # 邻居的相对顺序必须保住
///
/// 手性标记相对邻居存储顺序,而切分是在定基**之后**做的。按原下标升序放原子、
/// 按原键序建键,每个原子的邻居相对顺序就与切分前一致,标记照样成立。
fn split_components(
    shared: &MolBuilder,
    from_reactant: &[BTreeMap<u32, u32>],
) -> Vec<BuiltProduct> {
    let n = shared.num_atoms();
    let mut comp = vec![usize::MAX; n];
    let mut n_comp = 0usize;
    let mut stack: Vec<u32> = Vec::new();
    for s in 0..n as u32 {
        if comp[s as usize] != usize::MAX {
            continue;
        }
        comp[s as usize] = n_comp;
        stack.push(s);
        while let Some(a) = stack.pop() {
            for (other, _) in shared.neighbors(a) {
                if comp[other as usize] == usize::MAX {
                    comp[other as usize] = n_comp;
                    stack.push(other);
                }
            }
        }
        n_comp += 1;
    }

    let mut mols: Vec<MolBuilder> = (0..n_comp).map(|_| MolBuilder::new()).collect();
    // 共享图下标 → 该分量里的下标
    let mut local = vec![u32::MAX; n];
    for a in 0..n as u32 {
        let c = comp[a as usize];
        local[a as usize] = mols[c].add_atom_data(shared.atoms()[a as usize]);
    }
    for b in shared.bonds() {
        let c = comp[b.begin as usize];
        let mut nb = *b;
        nb.begin = local[b.begin as usize];
        nb.end = local[b.end as usize];
        nb.stereo_atoms = [
            translate_stereo_atom(b.stereo_atoms[0], &local),
            translate_stereo_atom(b.stereo_atoms[1], &local),
        ];
        let _ = mols[c].add_bond_data(nb);
    }

    // 出处表也要按分量拆开
    let mut tables: Vec<Vec<BTreeMap<u32, u32>>> = (0..n_comp)
        .map(|_| from_reactant.iter().map(|_| BTreeMap::new()).collect())
        .collect();
    for (ti, table) in from_reactant.iter().enumerate() {
        for (&src, &dst) in table {
            let c = comp[dst as usize];
            tables[c][ti].insert(src, local[dst as usize]);
        }
    }

    mols.into_iter().zip(tables).collect()
}

/// 参照原子的下标换算。哨兵值不换算 —— 它不是下标。
fn translate_stereo_atom(idx: u32, local: &[u32]) -> u32 {
    if idx == BondData::NO_STEREO_ATOM {
        return BondData::NO_STEREO_ATOM;
    }
    local
        .get(idx as usize)
        .copied()
        .filter(|&v| v != u32::MAX)
        .unwrap_or(BondData::NO_STEREO_ATOM)
}

/// 按运行时的原子对应关系,给反应物副本与产物打上映射号。
///
/// 规则连同理由写在 [`run_reactants`] 的文档里(那是调用方看得到的地方)。
/// 这里只补两处实现上的要点:
///
/// - 对应关系来自 [`BuiltProduct`] 的第二项,它由模板匹配与搬运共同填出;
///   产物侧新建的原子不在表里,自然拿不到号
/// - `first_home` 用 `BTreeMap` 而不是 `HashMap`:它的迭代顺序正好是发号要的
///   `(反应物, 原子下标)` 升序,顺手把"号必须确定"这条要求解决掉
///
/// `atom_mapping` 为假时直接把产物取出来,反应物侧留空,一次复制都不做。
fn stamp_atom_maps(
    reactants: &[(MolBuilder, MolProps)],
    built: Vec<BuiltProduct>,
    atom_mapping: bool,
) -> Outcome {
    let discarded = discarded_atoms(reactants, &built);
    if !atom_mapping {
        return Outcome {
            products: built.into_iter().map(|(m, _)| m).collect(),
            reactants: Vec::new(),
            discarded,
        };
    }

    let mut products: ProductSet = Vec::with_capacity(built.len());
    // (反应物下标, 反应物原子) → (第几个产物, 产物原子)。BTreeMap 的迭代顺序
    // 正好是发号要的顺序;`or_insert` 让先到的产物赢。
    let mut first_home: BTreeMap<(usize, u32), (usize, u32)> = BTreeMap::new();
    for (pi, (mol, per_reactant)) in built.into_iter().enumerate() {
        for (ti, table) in per_reactant.iter().enumerate() {
            for (&src, &dst) in table {
                first_home.entry((ti, src)).or_insert((pi, dst));
            }
        }
        products.push(mol);
    }

    let mut mapped: Vec<MolBuilder> = reactants.iter().map(|(m, _)| m.clone()).collect();
    for m in &mut mapped {
        for i in 0..m.num_atoms() as u32 {
            if let Some(a) = m.atom_mut(i) {
                a.atom_map = 0;
            }
        }
    }

    let mut next: u32 = 1;
    for (&(ti, src), &(pi, dst)) in &first_home {
        // u16 装不下更多号了。继续发会绕回去,把两个不同的原子说成同一个,
        // 那比留几个原子无号糟得多。
        let Ok(n) = u16::try_from(next) else { break };
        // 两侧都写得进去才发 —— 单边的号正是本函数要避免的东西
        if mapped[ti].atoms().get(src as usize).is_none()
            || products[pi].atoms().get(dst as usize).is_none()
        {
            continue;
        }
        if let Some(a) = mapped[ti].atom_mut(src) {
            a.atom_map = n;
        }
        if let Some(a) = products[pi].atom_mut(dst) {
            a.atom_map = n;
        }
        next += 1;
    }

    Outcome {
        products,
        reactants: mapped,
        discarded,
    }
}

/// 把拼接图上的丢弃原子下标切回**各个输入分子**的下标。
///
/// [`run_on_substrate`] 把输入拼成一张图跑,于是 `discarded` 只有一条、下标是
/// 拼接图的。可 [`Outcome::discarded`] 的契约是"第 i 个输入分子的原子下标" ——
/// **契约不该随入口而变**:调用方拿到 `discarded` 时不该还要先知道上游走的是
/// 哪个函数。
///
/// 不切回去的后果不是报错,是**静默算错**:下游拿拼接图的下标去索引原始分子,
/// 越界的被悄悄跳过,片段少了原子,账跟着错,而每一步都跑得通。
fn regroup_discarded(flat: &[Vec<u32>], sizes: &[usize]) -> Vec<Vec<u32>> {
    let mut out: Vec<Vec<u32>> = sizes.iter().map(|_| Vec::new()).collect();
    for a in flat.iter().flatten() {
        let mut rest = *a as usize;
        for (i, &n) in sizes.iter().enumerate() {
            if rest < n {
                out[i].push(u32::try_from(rest).unwrap_or(u32::MAX));
                break;
            }
            rest -= n;
        }
    }
    out
}

/// 每个输入分子里没有进入任何产物的原子。
///
/// 出处表是**按产物分量**拆开的,所以"进了产物"要对所有产物取并集再补 ——
/// 只看某一个产物会把搬进别的片段的原子误记成丢弃。
fn discarded_atoms(reactants: &[(MolBuilder, MolProps)], built: &[BuiltProduct]) -> Vec<Vec<u32>> {
    let mut kept: Vec<Vec<bool>> = reactants
        .iter()
        .map(|(m, _)| vec![false; m.num_atoms()])
        .collect();
    for (_, per_reactant) in built {
        for (ti, table) in per_reactant.iter().enumerate() {
            for &src in table.keys() {
                if let Some(slot) = kept[ti].get_mut(src as usize) {
                    *slot = true;
                }
            }
        }
    }
    kept.iter()
        .map(|flags| {
            flags
                .iter()
                .enumerate()
                .filter(|&(_, &k)| !k)
                .map(|(i, _)| u32::try_from(i).unwrap_or(u32::MAX))
                .collect()
        })
        .collect()
}

/// 手性标记是相对**邻居存储顺序**的,而产物的存储顺序与反应物不同 ——
/// 模板的键先建、搬过来的键后建,顺序被打乱了。
///
/// 不重新定基的话产物会是镜像分子:原子数、键集合、连通性全对,只有手性反了,
/// 纯拓扑比对永远发现不了。这与写出器里那次宇称换算是同一类问题。
///
/// 邻居**数目变了**的中心不处理:取代基增减之后手性本就没有定义,
/// 硬翻一次只会把一个未定义的值变成另一个未定义的值。
///
/// # 配位几何(`@SP`/`@TB`/`@OH`)走另一套换算
///
/// 四面体只有两种排列,换参照系就是"置换是奇是偶"。配位几何有 3 / 20 / 30 种,
/// 奇偶说不清它 —— 换算表在 [`omgkit_core::polyhedron::renumber`],写出器与
/// 规范化早就在用它,这里先前没接。
///
/// 后果:产物的邻居存储顺序**一定**变了(模板的键先建、搬运来的键后建),
/// 而 `@SP1` 原样照抄进一个不同的参照系,指的是另一个几何异构体。
///
/// # 要定基的是**产物**当前的标记,不是反应物的
///
/// 产物原子的标记未必等于它继承来的那个:模板可以写死一个构型,也可以刻意
/// 不写(那时 [`apply_template`] 会把它清掉)。拿反应物的标记来写,等于把
/// 模板刚做的决定又覆盖回去 —— 清掉的会被恢复,写死的会被换成继承值。
///
/// 换参照系用的是**邻居的对应关系**,与标记取值无关,所以这两件事可以分开:
/// 置换从两侧的邻居顺序算,取值从产物当前的标记取。
///
/// # 模板表过态的原子不定基
///
/// 定基换的是"反应物的邻居序 → 产物的邻居序",这只对模板**没管**的原子成立 ——
/// 它们的标记原样继承自反应物,记在反应物的参照系里。
///
/// 模板一旦在这个原子上写了手性([`ChiralityPlan`] 的后三档),标记说的就是
/// 产物自己参照系里的构型:`Set` 直接来自产物模板,而产物模板的键正是按模板
/// 顺序先建进产物图的;`Retain`/`Invert` 表达的是"与底物同构型/反构型",
/// 也是就产物这张图而言。再套一次反应物侧的置换,等于凭空多翻一道。
///
/// 这一档的触发面在小语料上是 0(22 条模板没有一条在产物侧写手性),
/// 真实模板里却很常见:立体专一的反应就是靠产物侧的 `@` 表达构型的。
/// 判据见 `harness/check_product_chirality.py`。
fn rebase_chirality(
    mol: &MolBuilder,
    kept: &BTreeMap<u32, u32>,
    settled_chirality: &BTreeSet<u32>,
    out: &mut MolBuilder,
) {
    for (&src, &dst) in kept {
        if settled_chirality.contains(&dst) {
            continue;
        }
        let tag = out.atoms()[dst as usize].chiral_tag;
        if tag == ChiralTag::Unspecified {
            continue;
        }
        let after: Vec<u32> = out.neighbors(dst).map(|(other, _)| other).collect();
        if !tag.is_tetrahedral() {
            rebase_coordination(mol, src, dst, kept, &after, out);
            continue;
        }
        // 反应物侧的邻居,按原顺序换算成产物下标。空出来的槽位下面补。
        //
        // 槽位空不空要看"**还连不连在这个中心上**",不是看那个原子有没有进产物。
        // 反应可以把一个邻居挪到别的产物片段去、同时接上一个新的:
        // `[C@@H:1]-[O:2]>>C-C(=O)-O-[C@@H:1].[OH:2]` 里 O:2 活着,只是不再连着
        // :1 了。只判"进没进产物"的话这个槽位算被占着,于是 before 里有个 after
        // 里没有的原子,置换算不出来(多重集都不同),重定基被**静默跳过** ——
        // 标记原样照抄,而产物的邻居顺序早变了,得到的是镜像。
        let slots: Vec<Option<u32>> = mol
            .neighbors(src)
            .map(|(other, _)| kept.get(&other).copied().filter(|p| after.contains(p)))
            .collect();
        let Some((before, after)) = align_for_rebase(&slots, &after) else {
            continue;
        };
        if omgkit_core::permutation_is_odd(&before, &after) == Some(true) {
            if let Some(a) = out.atom_mut(dst) {
                a.chiral_tag = tag.inverted();
            }
        }
    }
}

/// 配位几何的重定基:把排列序号从反应物的邻居序换算到产物的邻居序。
///
/// 两侧的配体必须是同一组(不增不减)—— 增减了的话这个标记本就没有意义。
/// 换算不出来就把标记整个丢掉:**表达不出来要说出来,照抄一个错的序号是撒谎。**
fn rebase_coordination(
    mol: &MolBuilder,
    src: u32,
    dst: u32,
    kept: &BTreeMap<u32, u32>,
    after: &[u32],
    out: &mut MolBuilder,
) {
    let tag = out.atoms()[dst as usize].chiral_tag;
    let perm = out.atoms()[dst as usize].stereo_perm;
    let before: Vec<u32> = mol
        .neighbors(src)
        .filter_map(|(other, _)| kept.get(&other).copied())
        .filter(|p| after.contains(p))
        .collect();
    let renumbered = if perm == 0 || before.len() != after.len() {
        None
    } else {
        omgkit_core::polyhedron::renumber(tag, perm, &before, after)
    };
    if let Some(a) = out.atom_mut(dst) {
        match renumbered {
            Some(p) => a.stereo_perm = p,
            None => {
                a.stereo_perm = 0;
                a.chiral_tag = ChiralTag::Unspecified;
            }
        }
    }
}

/// 代表**隐式氢**的哨兵。氢不是图里的节点,可它占着四面体的一个位置,
/// 换参照系时必须算进去。取 `u32::MAX` 是因为真实原子下标不可能是它。
pub(crate) const IMPLICIT_H: u32 = u32::MAX;

/// 把反应物侧与产物侧的邻居顺序对齐成两条等长、同多重集的序列,供求宇称用。
///
/// # 度数变了不等于放弃
///
/// 取代基被**隐式氢**接管是常事 —— 脱保护、脱羧、脱卤都是:
///
/// ```text
/// C-C(-C)(-C)-O-C(=O)-[C@@;H0;D4:1](-[C:2])(-[N:3])-[C:4]=[O:6]
///                  >> [C:4](=[O:6])-[C@H;D3:1](-[C:2])-[N:3]
/// ```
///
/// 中心从 D4H0 变成 D3H1。因为"长度对不上"就跳过重定基的话,标记会原样留在
/// **反应物的**参照系里,而产物的邻居顺序已经变了 —— 于是拿到镜像。拓扑、
/// 原子数、电荷全对,只有构型反了。USPTO-50k 上实测有这一档。
///
/// # 氢占着腾出来的那个位置
///
/// 所以做法是:少邻居的那一侧把哨兵**插在下标 1**,也就是紧跟第一个邻居 ——
/// 这与解析器留下的存储约定一致(`N[C@H](O)F` 的存储序是 `(N, O, F)`,标记
/// 说的是"从 N 看过去,H、O、F 依次逆时针")。
///
/// 插在哪一头其实**不影响结果**:四个邻居时把一个元素从首挪到尾是个三轮换,
/// 宇称不变。选下标 1 只是为了与约定对得上,读代码的人不必再推一遍。
///
/// 两个方向对称:取代基被氢接管时哨兵进**产物侧**,反应物侧原本的隐式氢被
/// 新邻居顶替时哨兵进**反应物侧**。
///
/// 只处理恰好差一个、且对齐后凑满四个邻居的情形 —— "对换一次就翻转"这条规则
/// 只对四面体成立。
pub(crate) fn align_for_rebase(
    slots: &[Option<u32>],
    after: &[u32],
) -> Option<(Vec<u32>, Vec<u32>)> {
    if let Some(before) = fill_replaced_slots(slots, after) {
        if before.len() == after.len() {
            return Some((before, after.to_vec()));
        }
    }
    let vacated = slots.iter().filter(|s| s.is_none()).count();
    let occupied = slots.len() - vacated;
    if vacated == 1 && occupied == after.len() && slots.len() == 4 {
        // 反应物侧多一个邻居:它在产物里被隐式氢接管
        let before: Vec<u32> = slots.iter().map(|s| s.unwrap_or(IMPLICIT_H)).collect();
        let mut aligned = after.to_vec();
        aligned.insert(1, IMPLICIT_H);
        return Some((before, aligned));
    }
    if vacated == 0 && after.len() == slots.len() + 1 && after.len() == 4 {
        // 产物侧多一个邻居:它顶了反应物侧原本的隐式氢
        let taken: BTreeSet<u32> = slots.iter().flatten().copied().collect();
        let mut fresh = after.iter().filter(|a| !taken.contains(a));
        let new = *fresh.next()?;
        if fresh.next().is_some() {
            return None;
        }
        let mut before: Vec<u32> = slots.iter().flatten().copied().collect();
        before.insert(1, new);
        return Some((before, after.to_vec()));
    }
    None
}

/// 把"被替换掉的取代基"那些空槽用产物侧新建的原子填上,保持原顺序。
///
/// # 为什么不能直接把空槽丢掉
///
/// `[C:1][OH]>>[C:1]Cl` 会删掉氧、新建一个氯。中心的度数没变,变的是邻居的
/// **身份**。把氧丢掉的话,反应物侧只剩 2 个邻居而产物侧有 3 个,长度对不上,
/// 重定基整个被跳过 —— 标记原样照抄,而产物侧的邻居顺序已经变了,于是手性反了。
///
/// 取代反应的几何含义是新取代基**占据被替换者原来的空间位置**,中心的构型
/// 不因此改变。所以这里按原顺序把空槽填上,新原子依次顶替。
///
/// # 一次换掉两个取代基:这里给的是一个**约定**
///
/// "新的顶替旧的"在只换一个时是唯一的。同一个中心上一次换掉两个就不唯一了:
///
/// ```text
/// [C@:1]([F:2])([Cl:3])([Br:4])[I:5]>>[C@:1]([N])([O])([I:5])[Br:4]
/// ```
///
/// N 顶 F 还是顶 Cl,两种配对宇称相反,得到的是一对对映体,而拓扑本身说不出是
/// 哪一个。本实现按**产物侧的邻居顺序**依次顶替。
///
/// 要强调的是:另一种做法("对应关系不唯一就干脆不定基")并不是"不猜",而是
/// 猜了恒等置换 —— 标记原样留在底物的参照系里,同样是一个约定。两者都推不出来,
/// 只能拿记录裁。真实语料上这两种约定分歧 3 条,记录判 **2:1** 支持这里的做法
/// (反应 6855、11651 本实现对,13618 另一种对);判据见
/// `harness/check_product_chirality.py`。
///
/// 空槽数与新原子数对不上时返回 `None` —— 那时连"依次顶替"都排不出来。
fn fill_replaced_slots(slots: &[Option<u32>], after: &[u32]) -> Option<Vec<u32>> {
    if slots.iter().all(Option::is_some) {
        return Some(slots.iter().flatten().copied().collect());
    }
    // 顶上来的是产物侧**还没占住槽位**的邻居。
    //
    // 不能只挑"没有反应物出处"的:接上来的那个可以是别处搬来的、有出处的原子,
    // 那样会一个都挑不到而白白放弃重定基。
    let taken: BTreeSet<u32> = slots.iter().flatten().copied().collect();
    let mut fresh = after.iter().filter(|a| !taken.contains(a));
    let filled: Option<Vec<u32>> = slots
        .iter()
        .map(|s| match s {
            Some(x) => Some(*x),
            None => fresh.next().copied(),
        })
        .collect();
    let filled = filled?;
    // 还剩下新原子没用上 —— 连接变了不止"替换"这么简单,不猜
    if fresh.next().is_some() {
        return None;
    }
    Some(filled)
}

/// 把反应物里未被模板占用的部分接到产物上。
/// 反应物里 `a` 与 `b` 之间那根键的方向,换算到"从 `a` 走向 `b`"的参照系。
///
/// 没有这根键、或它没有方向时返回 `None` 对应的无方向值。
/// `~` 那一档要沿用的键级:底物里 `a`–`b` 那根键的键级。没有这根键就是"未指定"。
fn inherited_order(mol: &MolBuilder, a: u32, b: u32) -> BondOrder {
    mol.neighbors(a)
        .find(|&(other, _)| other == b)
        .map_or(BondOrder::Unspecified, |(_, bi)| {
            mol.bonds()[bi as usize].order
        })
}

fn inherited_direction(mol: &MolBuilder, a: u32, b: u32) -> BondDirection {
    let Some((_, bi)) = mol.neighbors(a).find(|&(other, _)| other == b) else {
        return BondDirection::None;
    };
    let src = mol.bonds()[bi as usize];
    if src.begin == a {
        src.direction
    } else {
        src.direction.flipped()
    }
}

/// 旁观组分:模板**一个原子都没匹配到**的那些连通分量,原样搬进产物。
///
/// # 不搬就是丢原子
///
/// 底物写成 `内酰胺.HCl` 或 `[Na+].[O-]CC(=O)OCC[O-].[K+]` 时,反离子与任何
/// 匹配到的原子都不连通,[`carry_over`] 的遍历永远走不到它们。丢掉它们,产物的
/// 重原子数就少于底物 —— **引擎自己违反了质量守恒**,而且不报错。
///
/// 逆合成正是把模板作用到**任意**分子上,盐是常态而非边角,所以这一条不做成开关。
///
/// # 做法:埋一个种子,剩下的交给同一趟遍历
///
/// 给每个这样的组分往产物图里放一个种子原子、登记进 `kept`,`carry_over` 随后
/// 会从种子出发把整个组分连同它的键搬过来。这样键的排序纪律(按源键下标,不按
/// 遍历发现顺序)、手性与顺反的重定基、按连通分量切分、原子映射发号,全都与别处
/// 走同一条路径,不必再开一条。
///
/// # 判据是"有没有被匹配",不是"在不在 `kept` 里"
///
/// 组分里只要有**一个**原子被模板匹配到,它就归模板管:留下还是删掉是模板的表态,
/// 与旁观无关。模板删掉某个原子、挂在它身上的东西跟着失去落脚点,那是另一条约定
/// (见本模块的契约),这里不碰。
fn seed_spectators(
    mol: &MolBuilder,
    comp: &[u32],
    matched: &[bool],
    kept: &mut BTreeMap<u32, u32>,
    out: &mut MolBuilder,
) {
    // 单分量的分子不可能有旁观组分 —— 模板既然匹配上了,那唯一的分量就有匹配原子。
    // 绝大多数底物是这一档,先短路掉。
    let Some(&n_comp) = comp.iter().max() else {
        return;
    };
    if n_comp == 0 {
        return;
    }
    let n_comp = n_comp as usize + 1;

    // 哪些分量里有被匹配到的原子 —— 它们归模板管,不是旁观
    let mut has_match = vec![false; n_comp];
    for (a, &hit) in matched.iter().enumerate() {
        if hit {
            has_match[comp[a] as usize] = true;
        }
    }
    // 每个旁观分量埋一个种子:取分量里**下标最小**的原子。遍历的发现顺序不能用,
    // 它会让同一个分子因为写法不同而搬出不同的邻居次序。
    let mut seeded = vec![false; n_comp];
    for (a, &c) in comp.iter().enumerate() {
        let c = c as usize;
        if has_match[c] || seeded[c] {
            continue;
        }
        seeded[c] = true;
        let mut carried = mol.atoms()[a];
        carried.atom_map = 0;
        let idx = out.add_atom_data(carried);
        kept.insert(a as u32, idx);
    }
}

fn carry_over(
    mol: &MolBuilder,
    matched: &[bool],
    template_bonds: &[bool],
    kept: &mut BTreeMap<u32, u32>,
    out: &mut MolBuilder,
) {
    // 从已保留的原子出发做一次遍历,沿途把未匹配的原子拉进来
    let mut stack: Vec<u32> = kept.keys().copied().collect();
    let mut seen: Vec<bool> = vec![false; mol.num_atoms()];
    for &a in kept.keys() {
        seen[a as usize] = true;
    }
    // 每条边连它在**反应物里的键下标**一起记,建进产物时按这个下标排 ——
    // 搬过来的键因此保持反应物里的相对顺序。
    //
    // 遍历的发现顺序不能直接用:它取决于从哪个原子起步、栈怎么弹,与反应物
    // 自身的顺序无关。手性标记正是相对邻居顺序的,顺序一乱,产物就可能是镜像。
    // 模板只钉住一两根键的中心尤其明显 —— 余下的位置全由搬过来的键填,
    // 它们的先后直接定构型。
    let mut edges: Vec<(u32, u32, u32, BondData)> = Vec::new();

    while let Some(a) = stack.pop() {
        for (other, bi) in mol.neighbors(a) {
            let b = mol.bonds()[bi as usize];
            // `other` 被反应物模板匹配到、却不在**本**产物模板里 —— 它归别的
            // 产物片段(或压根没有映射号、该被删掉)。遍历到此为止。
            //
            // 少了这一条,断环键就出错:`[C:1][N:2]>>[C:1].[N:2]` 作用在
            // 皮啶上时,从 [C:1] 往外走会绕过整个环、从另一侧撞上那个氮,
            // 于是氮被拉进碳那一片,环又接了回去。断键反应在**环上**与在
            // 链上的区别正在这里 —— 链上走不回去,环上走得回去。
            if matched[other as usize] && !kept.contains_key(&other) {
                continue;
            }
            // 反应物模板**亲自匹配到**的键才归产物模板负责,不搬。
            //
            // 判据不能是"两端都被匹配就不搬"。子结构匹配只要求模板的每根键
            // 在底物里找得到,并不要求底物在这些原子之间没有别的键 —— 模板
            // 把一个环写成**开链路径**时,环闭合的那根键两端确实都被匹配了,
            // 可模板从没看见它。按"两端都匹配"判,这根键就被当成模板的地盘
            // 删掉,环被撕开。撕开之后报出来的是**芳香**错误(原子不在环中
            // 却带着芳香标志),病因与症状隔着一层,极难回溯。
            //
            // rdchiral 抽模板时只沿反应中心走一趟,稠环因此普遍写成路径,
            // 这不是边角情形。
            //
            // 反过来也不能一律搬:模板匹配到、产物侧又不写的键,正是**断键**
            // 反应要表达的东西。两条判据缺一不可。
            if template_bonds[bi as usize] {
                continue;
            }
            if !seen[other as usize] {
                seen[other as usize] = true;
                // 搬过来的原子同样清掉映射号,理由同 apply_template
                let mut carried = mol.atoms()[other as usize];
                carried.atom_map = 0;
                let idx = out.add_atom_data(carried);
                kept.insert(other, idx);
                stack.push(other);
            }
            edges.push((bi, a, other, b));
        }
    }
    edges.sort_by_key(|&(bi, ..)| bi);

    let mut done: std::collections::HashSet<(u32, u32)> = std::collections::HashSet::new();
    for (_, a, b, src) in edges {
        let key = if a <= b { (a, b) } else { (b, a) };
        if !done.insert(key) {
            continue;
        }
        // 新键要沿用**源键的朝向**,而不是遍历到它的方向。
        //
        // `a` 是遍历的出发端,与 `src.begin` 不一定是同一个原子。两者相反时
        // 若按遍历方向建键,朝向就翻了 —— 而朝向是有语义的:
        //
        // - `direction`(`/` `\`)相对 `begin → end`,翻转即把顺式写成反式
        // - 配位键的箭头必须从给电子的一端指向受体
        //
        // 照源键的 begin/end 映射过去,这些量就能原样抄,不必逐个换参照系;
        // 少一次换算就少一处会悄悄写反的地方。
        let (Some(&na), Some(&nb)) = (kept.get(&src.begin), kept.get(&src.end)) else {
            continue;
        };
        // 产物模板已经建过这根键 —— 成环模板作用在**已经成环**的底物上就是
        // 这样:闭合键不在反应物模板里(反应物侧写的是一条链),却被产物模板
        // 明写了出来。模板说了算,再搬一遍就是重复键:不报错,价数却翻倍。
        if out.bond_between(na, nb).is_some() {
            continue;
        }
        // 整条键的属性都要跟过来:只抄键级的话,芳香键的标志位会丢,
        // 而"键级为 Aromatic 时标志位必须同步"是全局不变量 —— 破了它,
        // 写出来的芳香环会变成"大写原子 + 冒号键"这种半吊子形式
        let mut nb_data = BondData::new(na, nb, src.order);
        nb_data.direction = src.direction;
        // 顺反照抄,参照原子先留空 —— 挑参照要看**产物**的连接,而这会儿
        // 图还没建完。留空是给 `rebase_bond_stereo` 认的记号,见那里。
        nb_data.stereo = src.stereo;
        nb_data.stereo_atoms = [BondData::NO_STEREO_ATOM; 2];
        nb_data.flags = src.flags;
        let _ = out.add_bond_data(nb_data);
    }
}

/// 给搬运过来的双键在**产物**里重挑顺反的参照原子。
///
/// # 参照必须在产物里还连着这根双键
///
/// 顺反记在双键上,可它是**相对两个参照原子**说的。反应能动到参照原子的方式
/// 有两种,判据都不能只看"那个原子进没进产物":
///
/// - 参照被**删掉**(模板里没有映射号)
/// - 参照活着,却被挪到了**别的产物片段**去
///
/// 后一种最阴:`kept` 里查得到,于是参照被换成一个属于另一个分子的下标。切分
/// 之后那个下标在本分子里要么越界(几何静默丢失),要么正好落在某个真邻居上
/// (几何**静默变错**)。两种都不报错。
///
/// # 顶替者占的是被顶替者的**位置**
///
/// 挑替代不能随便找同端的另一个取代基 —— 那个在双键的另一侧,顺反的含义跟着
/// 反。取代的几何含义是新取代基占据被替换者原来的空间位置,所以按**槽位**挑:
/// 把该端的取代基按反应物里的顺序排开,走掉的留空,再拿产物侧多出来的邻居依次
/// 顶上。顶替者与被顶替者同侧,顺反因此原样成立,一次都不用翻。
///
/// # 为什么放在这里而不是搬运的时候
///
/// 挑参照要看产物侧该端还连着谁,而搬运时图还没建完。`carry_over` 因此只抄
/// `stereo`、把参照留成哨兵值,由本函数认这个记号再填。模板自己建的键不带
/// `stereo`,所以不会被误认。
fn rebase_bond_stereo(mol: &MolBuilder, kept: &BTreeMap<u32, u32>, out: &mut MolBuilder) {
    for src in mol.bonds() {
        if src.stereo == BondStereo::None
            || src.stereo_atoms[0] == BondData::NO_STEREO_ATOM
            || src.stereo_atoms[1] == BondData::NO_STEREO_ATOM
        {
            continue;
        }
        let (Some(&pb), Some(&pe)) = (kept.get(&src.begin), kept.get(&src.end)) else {
            continue;
        };
        let Some(bi) = out.bond_between(pb, pe) else {
            continue;
        };
        // 只认 `carry_over` 留下的记号 —— 模板建的键不走这一路
        let cur = out.bonds()[bi as usize];
        if cur.stereo == BondStereo::None || cur.stereo_atoms[0] != BondData::NO_STEREO_ATOM {
            continue;
        }

        let mut refs = [BondData::NO_STEREO_ATOM; 2];
        // 换到"另一侧那个取代基"上的次数。奇数次就要把顺反翻过来。
        let mut flips = 0usize;
        for (i, (end, other, p_end, p_other)) in
            [(src.begin, src.end, pb, pe), (src.end, src.begin, pe, pb)]
                .into_iter()
                .enumerate()
        {
            let want = src.stereo_atoms[i];
            // 产物侧该端的取代基(不含双键另一端)
            let p_subs: Vec<u32> = out
                .neighbors(p_end)
                .map(|(o, _)| o)
                .filter(|&o| o != p_other)
                .collect();
            // 参照原封不动地还连着 —— 直接用
            if let Some(&p) = kept.get(&want) {
                if p_subs.contains(&p) {
                    refs[i] = p;
                    continue;
                }
            }
            // 走掉了(被删,或挪去了别的片段)—— 找占了它那个槽位的
            let subs: Vec<u32> = mol
                .neighbors(end)
                .map(|(o, _)| o)
                .filter(|&o| o != other)
                .collect();
            let Some(pos) = subs.iter().position(|&o| o == want) else {
                break;
            };
            let slots: Vec<Option<u32>> = subs
                .iter()
                .map(|o| kept.get(o).copied().filter(|p| p_subs.contains(p)))
                .collect();
            let Some(filled) = fill_replaced_slots(&slots, &p_subs) else {
                // **没人顶这个槽位** —— 走掉的位置由隐式氢补上,而隐式氢没有
                // 下标,做不了参照原子。这时候要**改参照到该端另一个取代基,
                // 并把顺反翻一次**:它在双键的另一侧。
                //
                // 先前这里直接放弃,整根双键的顺反跟着作废。实测大语料上
                // 硝基烯烃 `[N+](/C(=C/Ar)C)([O-])=O` 被 `[C:1][N:2]>>[C:1].[N:2]`
                // 打掉硝基之后,交付的是 `CC=Cc1ccccc1` —— **E/Z 整个没了**,
                // 而双键还在、两端各自还有取代基,构型依旧成立。
                //
                // 谁对不是推的:把底物嵌成真实三维构象、把离去的氮**原地**换成
                // 氢再读回构型,三个分子五个 seed 都给出 `C/C=C\Ar` 这一类
                // (判据先自校准过:同一条路读底物本身,五个 seed 都还原输入)。
                let Some(&alt) = subs.iter().find(|&&o| o != want) else {
                    // 该端只有走掉的那一个取代基 —— 换成两个隐式氢,
                    // 这根双键**真的**没有构型可言了,作废是对的
                    break;
                };
                let Some(&p_alt) = kept.get(&alt) else {
                    break;
                };
                if !p_subs.contains(&p_alt) {
                    break;
                }
                refs[i] = p_alt;
                flips += 1;
                continue;
            };
            refs[i] = filled[pos];
        }

        if let Some(mut b) = out.bond_mut(bi) {
            if refs[0] == BondData::NO_STEREO_ATOM || refs[1] == BondData::NO_STEREO_ATOM {
                // 挑不出参照就作废 —— 留一个指向别人的下标比没有更糟
                b.set_stereo(BondStereo::None);
            } else if flips % 2 == 0 {
                b.set_stereo_atoms(refs);
            } else {
                match src.stereo {
                    BondStereo::Cis => {
                        b.set_stereo(BondStereo::Trans);
                        b.set_stereo_atoms(refs);
                    }
                    BondStereo::Trans => {
                        b.set_stereo(BondStereo::Cis);
                        b.set_stereo_atoms(refs);
                    }
                    // Z/E 是按 **CIP 优先级**定的,与记录的参照原子无关 ——
                    // 换参照不该翻它;而取代基换掉之后 CIP 排序本身也可能变,
                    // 那要重新定优先级,不是翻个号能解决的。不猜,作废。
                    _ => b.set_stereo(BondStereo::None),
                }
            }
        }
    }
}

/// 把产物模板里写死的属性盖到原子上。
///
/// 模板里没写的属性**保持继承来的值** —— 这正是"分子其余部分自动跟着走"
/// 的原子级体现:`[C:1]` 只说"这里是个碳",电荷、同位素都不动。
///
/// 模板里的映射号**不**写进产物:它连的是两个模板,不是分子的属性。留在产物
/// 里的话,写出的 SMILES 会带上 `[CH2:1]` 这种本不该有的标注,而且下一次拿这个
/// 产物当底物时,那些号会跟新模板的号撞上。
///
/// 要的是底物层面的原子对应关系时,用 [`run_reactants`] 的 `atom_mapping`
/// 参数 —— 那套号是运行时另发的,见 [`stamp_atom_maps`]。
fn apply_template(
    mut base: AtomData,
    expr: &AtomExpr,
    degree_kept: bool,
    plan: ChiralityPlan,
) -> AtomData {
    // 自由基电子数是**派生量**,不是原子的固有属性:它由具体的 Kekulé 结构、
    // 电荷与氢数一起定下来,而模板恰恰会把这三样都改掉。继承过来就是一个陈旧值。
    //
    // 陈旧在哪不显眼:净化里 kekulize 排在自由基重算**之前**(自由基数要等键级
    // 定下来才算得出),于是那个陈旧值会被 kekulize 当真。一个三价碳
    // (`[C]`,带一个自由基)被模板改写成芳香碳之后,kekulize 认为它不能再要
    // 双键,整个芳香环就配不出 Kekulé 结构 —— 报的错落在环上某个无辜的原子身上,
    // 离根因很远。
    //
    // 清成 0 之后,产物走的路与"把这个分子写出来再读回去"完全一致:
    // 新解析的分子这个字段本来就是 0,由 `assign_radicals` 在 kekulize 之后重算。

    base.num_radical_electrons = 0;

    // 元素一旦被模板改掉,继承来的一切都失去意义 —— `[OH:2]` 变成 `[Cl:2]`
    // 时若把氧的那个氢留下,得到的是 ClH,直接超价。电荷、同位素同理。
    let element_changed = template_element(expr).is_some_and(|z| z != base.atomic_num);

    // 氢数只在"元素没变**且**连接没变"时才继承。断了一条键的原子要补氢:
    // `[C:1][N:2]>>[C:1].[N:2]` 把丙氨酸的 C—N 断开,那个碳应当从 CH 变成
    // CH2。照抄氢数的话会得到一个凭空少了个氢的自由基。
    if element_changed || !degree_kept {
        base.num_explicit_hs = 0;
        base.num_implicit_hs = 0;
        base.flags.remove(omgkit_core::AtomFlags::NO_IMPLICIT);
    }
    if element_changed {
        base.formal_charge = 0;
        base.isotope = 0;
        base.chiral_tag = ChiralTag::Unspecified;
    }
    // 继承来的构型,`apply_expr` 有可能把它盖掉,所以先存下来
    let inherited = base.chiral_tag;
    apply_expr(&mut base, expr);
    base.chiral_tag = match plan {
        ChiralityPlan::Inherit | ChiralityPlan::Set => base.chiral_tag,
        ChiralityPlan::Drop => ChiralTag::Unspecified,
        ChiralityPlan::Retain => inherited,
        ChiralityPlan::Invert => inherited.inverted(),
    };
    base.atom_map = 0;
    base
}

/// 产物模板给这个原子指定的元素。写成析取式(`[C,N:1]`)时说不出是哪个,
/// 返回 `None`。
fn template_element(expr: &AtomExpr) -> Option<u8> {
    match expr {
        AtomExpr::Prim(AtomPrim::Element { z, .. }) => Some(*z),
        AtomExpr::And(parts) => parts.iter().find_map(template_element),
        _ => None,
    }
}

fn apply_expr(a: &mut AtomData, expr: &AtomExpr) {
    match expr {
        AtomExpr::Prim(p) => apply_prim(a, p),
        AtomExpr::And(parts) => {
            for p in parts {
                apply_expr(a, p);
            }
        }
        // 析取与否定在产物侧没有确定含义 —— 写 `[C,N:1]` 说不出该建哪个,
        // 所以一律忽略,保留继承来的值
        AtomExpr::Or(_) | AtomExpr::Not(_) => {}
    }
}

fn apply_prim(a: &mut AtomData, p: &AtomPrim) {
    match p {
        AtomPrim::Element { z, aromatic } => {
            a.atomic_num = *z;
            if let Some(arom) = aromatic {
                a.flags.set(omgkit_core::AtomFlags::AROMATIC, *arom);
            }
        }
        AtomPrim::Charge(c) => a.formal_charge = i8::try_from(*c).unwrap_or(0),
        AtomPrim::Isotope(i) => a.isotope = *i,
        AtomPrim::TotalHs(n) => {
            a.num_explicit_hs = u8::try_from(*n).unwrap_or(0);
            a.num_implicit_hs = 0;
            a.flags.insert(omgkit_core::AtomFlags::NO_IMPLICIT);
        }
        AtomPrim::Chirality(t) => a.chiral_tag = *t,
        // 其余基元是**筛选**条件,不是构建指令:`[C;R1:1]` 里的 R1 说的是
        // "只匹配环上的碳",不是"把产物做成环"
        _ => {}
    }
}

/// 这条键表达式写的是 `<-` 吗。
///
/// 查询侧用两个基元区分配位键的朝向,端点则按书写顺序存;产物侧靠端点顺序
/// 表达朝向。两种表示之间要换算,靠的就是这个判断。
fn is_dative_reversed(expr: &BondExpr) -> bool {
    match expr {
        BondExpr::Prim(BondPrim::DativeReversed) => true,
        BondExpr::And(parts) => parts.iter().any(is_dative_reversed),
        // 析取与否定说不出确定的朝向,按写法原样建
        _ => false,
    }
}

/// 产物模板里的键表达式指定的方向(`/` `\`)。没写方向时返回 `None` 对应的值。
///
/// 方向与键级是**两件事**:`/` 既说"这是单键",也说"取代基在双键的哪一侧"。
/// [`product_bond_from`] 只取前者,后者要靠这里取,否则模板里写的几何会被
/// 悄悄丢掉 —— 产物从确定的顺反退化成未指定。
fn bond_direction_from(expr: &BondExpr) -> BondDirection {
    match expr {
        BondExpr::Prim(BondPrim::UpRight) => BondDirection::UpRight,
        BondExpr::Prim(BondPrim::DownRight) => BondDirection::DownRight,
        // 合取式里任一支写了方向就算数(`/&!@` 这类)
        BondExpr::And(parts) => parts
            .iter()
            .map(bond_direction_from)
            .find(|d| *d != BondDirection::None)
            .unwrap_or(BondDirection::None),
        // 析取与否定说不出确定的方向:`/,\` 是"两侧都行",不是某一侧
        _ => BondDirection::None,
    }
}

/// 产物模板里的一根键该建成什么 —— 三种情形,不是一种。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ProductBond {
    /// 模板写死了键级
    Fixed(BondOrder),
    /// **省略了键符号。** 两端都是芳香原子就建芳香键,否则单键。
    FollowAromaticity,
    /// `~` —— 沿用底物那根键;底物没有这根键时是"未指定"。
    Inherit,
}

/// 产物模板里的键表达式该建成什么键。
///
/// # 省略键符号**不是**单键
///
/// SMARTS 里省略键符号的默认值是析取 `单键 或 芳香键`([`BondExpr::default_bond`]),
/// 而这正是产物模板最常见的写法。先前这里把析取整个退回单键,于是
/// `>>[c:1]1[c:2][c:3][c:4][n:5]1` 建出来的是吡咯**烷** —— 五个芳香原子之间连
/// 五根单键,而这样的产物净化得过、不报任何错,调用方拿到一个结构良好的错分子。
///
/// 参照实现的规矩在 `ReactionRunner.cpp:391-406`:析取默认值时看两端原子的
/// 芳香性 —— 都芳香就建芳香键,否则单键;`~` 另算,标记成"沿用底物"。这里照办。
///
/// 其余析取/否定说不出确定的键级,取第一个说得出的子式;一个都没有就按省略处理。
fn product_bond_from(expr: &BondExpr) -> ProductBond {
    match expr {
        BondExpr::Prim(BondPrim::Any) => ProductBond::Inherit,
        BondExpr::Prim(p) => ProductBond::Fixed(match p {
            BondPrim::Double => BondOrder::Double,
            BondPrim::Triple => BondOrder::Triple,
            BondPrim::Quadruple => BondOrder::Quadruple,
            BondPrim::Aromatic => BondOrder::Aromatic,
            BondPrim::Dative | BondPrim::DativeReversed => BondOrder::Dative,
            _ => BondOrder::Single,
        }),
        BondExpr::And(parts) => parts
            .iter()
            .map(product_bond_from)
            .find(|o| !matches!(o, ProductBond::Fixed(BondOrder::Single)))
            .unwrap_or(ProductBond::Fixed(BondOrder::Single)),
        BondExpr::Or(_) | BondExpr::Not(_) => {
            if *expr == BondExpr::default_bond() {
                return ProductBond::FollowAromaticity;
            }
            let parts = match expr {
                BondExpr::Or(parts) => parts.as_slice(),
                _ => &[],
            };
            parts
                .iter()
                .map(product_bond_from)
                .find(|o| matches!(o, ProductBond::Fixed(_)))
                .unwrap_or(ProductBond::FollowAromaticity)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 两侧模板的邻居次序怎么算宇称,以及什么时候该放弃。
    ///
    /// 端到端那条判据在 `tests/reaction.rs`
    /// (`both_sides_written_are_compared_in_a_common_neighbour_order`);
    /// 这里守的是边界:对应关系不唯一时必须给 `None`,而不是随便算一个宇称
    /// 出来 —— 算出来的那个会被当真,把构型悄悄写反。
    #[test]
    fn template_order_parity_gives_up_when_the_correspondence_is_not_unique() {
        let n = |v: &[u16]| -> Vec<Option<u16>> { v.iter().map(|&x| Some(x)).collect() };

        // 次序相同 —— 偶
        assert_eq!(
            template_order_is_odd(&n(&[2, 3, 4]), &n(&[2, 3, 4])),
            Some(false)
        );
        // 对调一对 —— 奇
        assert_eq!(
            template_order_is_odd(&n(&[2, 3, 4]), &n(&[4, 3, 2])),
            Some(true)
        );
        // 轮换一圈是两次对调 —— 偶
        assert_eq!(
            template_order_is_odd(&n(&[2, 3, 4]), &n(&[3, 4, 2])),
            Some(false)
        );

        // 每侧各有一个对方没有的邻居:互相顶替,对应唯一
        let mut react = n(&[2, 3, 4]);
        react[0] = None;
        let mut prod = n(&[2, 3, 4]);
        prod[2] = None;
        assert!(
            template_order_is_odd(&react, &prod).is_some(),
            "各有一个对不上时该顶替得起来"
        );

        // 一侧两个对不上 —— 两种配对宇称相反,不能挑一个
        assert_eq!(
            template_order_is_odd(&n(&[2, 3, 4]), &n(&[5, 6, 4])),
            None,
            "产物侧有两个邻居在反应物侧找不到,对应关系不唯一"
        );
        // 度数不足 3 谈不上四面体手性
        assert_eq!(template_order_is_odd(&n(&[2, 3]), &n(&[2, 3])), None);
        // 两侧差出一个以上
        assert_eq!(
            template_order_is_odd(&n(&[2, 3, 4]), &n(&[2, 3, 4, 5, 6])),
            None
        );
    }
}
