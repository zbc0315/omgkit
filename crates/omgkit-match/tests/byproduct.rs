//! 副产物收口。
//!
//! 判据分两层,**下面那一层比上面那一层硬**:
//!
//! 1. 收口出来的分子的规范 SMILES —— 看给得对不对
//! 2. **质量守恒**:底物的重原子数 = 产物 + 副产物的重原子数,逐条精确相等。
//!    这一条不依赖任何记录,也不依赖判据作者认得几个副产物 —— 它是本模块唯一
//!    的正确性来源,所以每个能闭合的用例都要过一遍。
//!
//! 收口不了的用例同样要测:那时**不能给分子**。编一个出来是最坏的结果 ——
//! 它拓扑合法、能净化、看不出破绽,只是错的。

use omgkit_chem::sanitize;
use omgkit_core::MolBuilder;
use omgkit_io::{canon, smarts, smiles};
use omgkit_match::byproduct::{self, Unresolved, Verdict};
use omgkit_match::{run_on_substrate, run_reactants, MolProps};

fn sanitized(smi: &str) -> MolBuilder {
    let mut m = smiles::parse(smi).unwrap_or_else(|e| panic!("{smi}: {}", e.render()));
    sanitize(&mut m).unwrap_or_else(|e| panic!("{smi}: {e}"));
    omgkit_io::stereo::perceive_bond_stereo(&mut m);
    m
}

fn heavy(mol: &MolBuilder) -> usize {
    mol.atoms().iter().filter(|a| a.atomic_num != 1).count()
}

/// 跑一条反应,返回第一组结果的(产物规范式, 副产物规范式, 结论)。
fn run(rxn_smarts: &str, reactant_smis: &[&str]) -> (Vec<String>, Vec<String>, Verdict) {
    let rxn = smarts::parse_reaction(rxn_smarts)
        .unwrap_or_else(|e| panic!("{rxn_smarts}:\n{}", e.render()));
    let mols: Vec<MolBuilder> = reactant_smis.iter().map(|s| sanitized(s)).collect();
    let inputs: Vec<(MolBuilder, MolProps)> = mols
        .iter()
        .map(|m| (m.clone(), MolProps::compute(m)))
        .collect();

    let outs = run_reactants(&rxn, &inputs, 0, false);
    let outcome = outs.first().expect("这条反应应当出产物");
    let by = byproduct::reconstruct(&mols, outcome);

    let canonical = |m: &MolBuilder| {
        let mut c = m.clone();
        match sanitize(&mut c) {
            Ok(()) => canon::canonical_smiles(&c).smiles,
            Err(e) => format!("<净化失败: {e}>"),
        }
    };
    let products: Vec<String> = outcome.products.iter().map(canonical).collect();
    let byproducts: Vec<String> = by.molecules.iter().map(canonical).collect();

    // 能闭合就必须守恒 —— 这是本文件最硬的一条,每个用例都过
    if by.verdict.is_closed() {
        let got: usize = outcome.products.iter().map(heavy).sum::<usize>()
            + by.molecules.iter().map(heavy).sum::<usize>();
        let want: usize = mols.iter().map(heavy).sum();
        assert_eq!(
            got, want,
            "质量不守恒:底物 {want} 个重原子,产物+副产物 {got} 个\n\
             模板 {rxn_smarts}\n产物 {products:?}\n副产物 {byproducts:?}"
        );
    }
    (products, byproducts, by.verdict)
}

/// 酯化丢掉的是水,不是"什么都没有"。
///
/// 这是整件事最小的用例:酸的羟基氧被删掉(它带着 1 个氢、欠 1 处价),
/// 醇的氧则少了 1 个氢。两笔加起来正好是 H2O。
#[test]
fn esterification_gives_back_water() {
    let (products, byproducts, verdict) = run(
        "[C:1](=[O:2])[OH:3].[OH:4][C:5]>>[C:1](=[O:2])[O:4][C:5]",
        &["CC(=O)O", "CCO"],
    );
    assert_eq!(products, ["CC(=O)OCC"]);
    assert_eq!(byproducts, ["O"]);
    assert_eq!(verdict, Verdict::Capped);
}

/// 酰胺缩合同样掉一个水。
#[test]
fn amide_coupling_gives_back_water() {
    let (_, byproducts, verdict) = run(
        "[C:1](=[O:2])[OH:3].[N:4]>>[C:1](=[O:2])[N:4]",
        &["CC(=O)O", "NCc1ccccc1"],
    );
    assert_eq!(byproducts, ["O"]);
    assert_eq!(verdict, Verdict::Capped);
}

/// 离去的卤素拿到氢,变成卤化氢 —— 而不是留下一个裸原子。
#[test]
fn a_leaving_halide_comes_back_as_the_hydrogen_halide() {
    let (_, byproducts, verdict) = run("[C:1][Cl:2].[O:3][C:4]>>[C:1][O:3][C:4]", &["CCCl", "OCC"]);
    assert_eq!(byproducts, ["Cl"]);
    assert_eq!(verdict, Verdict::Capped);
}

/// 叔丁酯水解:丢掉的叔丁基带的氢**比预算多一个**,于是摘一个氢、成一根 π 键,
/// 给出异丁烯而不是异丁烷。
///
/// 这一条守的是 `borrow_hydrogens`:预算为负时不能一味补氢,得反过来摘。摘完
/// 空价必须落在**相邻**的原子上,两处空价才成得了双键 —— 落在远处就只能成环
/// 或者干脆收不了口。
#[test]
fn a_tert_butyl_ester_leaves_as_the_alkene_not_the_alkane() {
    let (products, byproducts, verdict) = run(
        "[C:1](=[O:2])[O:3][C](C)(C)C>>[C:1](=[O:2])[OH:3]",
        &["CC(=O)OC(C)(C)C"],
    );
    assert_eq!(products, ["CC(O)=O"]);
    assert_eq!(byproducts, ["CC(C)=C"]);
    assert_eq!(verdict, Verdict::Bonded { bonds: 1 });
}

/// 乙酯掉下来的乙基收成乙烯。
///
/// **这一条守的是消除的形状本身,不守 `has_spare_hydrogen`** —— 乙基上与断点
/// 相邻的原子恰好都带氢,新旧两种写法给同一个答案,判据碰不到分歧点。那一处
/// 由 [`a_fragment_whose_neighbours_carry_no_hydrogen`] 守,两条不可互相替代。
#[test]
fn an_ethyl_ester_leaves_as_ethylene() {
    let (_, byproducts, verdict) = run(
        "[C:1](=[O:2])[O:3][CH2]C>>[C:1](=[O:2])[OH:3]",
        &["CC(=O)OCC"],
    );
    assert_eq!(byproducts, ["C=C"]);
    assert_eq!(verdict, Verdict::Bonded { bonds: 1 });
}

/// 与断点相邻的原子**一个氢都没有**时,摘氢要摘到更远处去。
///
/// 这一条守 `has_spare_hydrogen`:先前按"非氢非卤就当它出得起氢"来估,于是
/// 醚氧、酯羰基碳这类不带氢的原子也会被选中,空价落上去、成键之后当场超价。
///
/// Cbz 正是这个形状 —— 断点在氨基甲酸酯的羰基碳上,它的两个邻居(羰基氧、
/// 酯氧)都不带氢。**乙酯那条判据抓不到它**:那里相邻原子恰好都带氢,新旧写法
/// 结论相同。实测把老写法塞回去,乙酯那条照样绿,只有这一条会红。
///
/// 顺带记一件事:这里给出的 `c1ccccc1C1C(=O)O1` 是个三元 α-内酯 —— 账是平的,
/// 化学上却不是实际分离到的那个(实际是 CO₂ + 甲苯)。判据守的是**账**与**收口
/// 走对了哪条路**,不是"这个分子合不合成化学直觉";后者要靠分解规则表,
/// 见模块文档"形式副产物"一节。
#[test]
fn a_fragment_whose_neighbours_carry_no_hydrogen() {
    let (_, byproducts, verdict) = run("[N:1]C(=O)OCc1ccccc1>>[NH2:1]", &["CNC(=O)OCc1ccccc1"]);
    assert_eq!(byproducts, ["c1ccccc1C1C(=O)O1"]);
    assert_eq!(verdict, Verdict::Bonded { bonds: 1 });
}

/// 方括号写法的离去基团也要补上氢。
///
/// 这一条与 [`esterification_gives_back_water`] 看着重复,其实走的是**另一条**
/// 代码路径。裸写的 `O` 没有 `NO_IMPLICIT`,净化会按价规则自己把氢补齐,收口
/// 什么都不用做;而 `[OH]` 的氢数是写死的,净化不会动它,必须由
/// `settle_hydrogens` 手工加。
///
/// 少了这一条,把手工加的那一步删掉,上面那些用例**照样全绿** —— 实测过。
/// 同一个化学结果由两条路径产出,判据必须两条都走到。
#[test]
fn a_leaving_group_written_in_brackets_is_also_capped() {
    let (_, byproducts, verdict) = run(
        "[C:1](=[O:2])[OH:3].[OH:4][C:5]>>[C:1](=[O:2])[O:4][C:5]",
        &["CC(=O)[OH]", "CCO"],
    );
    assert_eq!(byproducts, ["O"]);
    assert_eq!(verdict, Verdict::Capped);
}

/// 季铵化:溴以 **Br⁻** 离去 —— 一处空价由电荷填掉,一个氢都不需要。
///
/// 只按氢记账的话这一档永远配不成对(空价 1、要补的氢 0,剩 1 处配不成对),
/// 会被误判成"记录不平"而拒答。而记录本身是完整自洽的:产物带 +1,副产物就该
/// 带 −1,电荷守恒把这一处空价交代得清清楚楚。
///
/// 这是**归因分析在真实语料上翻出来的**一处判据缺陷,不是构造出来的边角情形:
/// 亲核取代、季铵化这些反应里,离去基团本来就以阴离子离去。
#[test]
fn a_leaving_group_that_departs_as_an_anion_closes_on_charge() {
    let (products, byproducts, verdict) =
        run("[C:1][Br:2].[N:3]>>[C:1][N+:3]", &["CCBr", "CN(C)C"]);
    assert_eq!(products, ["C[N+](C)(C)CC"]);
    assert_eq!(byproducts, ["[Br-]"]);
    assert_eq!(verdict, Verdict::Capped);
}

/// 上一条的防空过:电荷不是随便加的,总量由**电荷守恒**定死。
///
/// 产物不带电时,副产物也不该带电 —— 同一个卤素这时要拿氢变成 HBr。两条用例
/// 走同一段代码,给出的却是带电与不带电两种结果,说明电荷那一项确实是算出来的,
/// 不是写死的。
#[test]
fn the_same_halide_stays_neutral_when_the_product_is_neutral() {
    let (_, byproducts, verdict) = run("[C:1][Br:2].[O:3]>>[C:1][O:3]", &["CCBr", "OC"]);
    assert_eq!(byproducts, ["Br"]);
    assert_eq!(verdict, Verdict::Capped);
}

/// 两个入口必须给出同一个副产物。
///
/// `run_on_substrate` 把输入**拼成一张图**跑,于是 `discarded` 只有一条、下标是
/// 拼接图的;而 [`Outcome::discarded`](omgkit_match::Outcome::discarded) 的契约是
/// "第 i 个输入分子的原子下标"。不切回去的后果**不是报错,是静默算错**:收口时
/// 拿拼接图的下标去索引原始分子,越界的被悄悄跳过,片段少了原子,账跟着错,
/// 而每一步都跑得通。
///
/// 判据挑"分子间底物"是有意的:那时两个入口本该等价,任何差别都是实现的问题。
/// 这也是唯一能把这类错照出来的形状 —— 单分子底物上两条路径的下标恰好重合。
#[test]
fn both_entry_points_agree_on_the_byproduct() {
    let rxn = smarts::parse_reaction("[OH:3][C:4].[C:1][Cl]>>[C:1][O:3][C:4]").unwrap();
    let mols = [sanitized("OCC"), sanitized("CCCl")];
    let inputs: Vec<(MolBuilder, MolProps)> = mols
        .iter()
        .map(|m| (m.clone(), MolProps::compute(m)))
        .collect();

    let by_positional = byproduct::reconstruct(&mols, &run_reactants(&rxn, &inputs, 0, false)[0]);
    let by_graph = byproduct::reconstruct(&mols, &run_on_substrate(&rxn, &inputs, 0, false)[0]);

    let names = |b: &byproduct::Byproducts| -> Vec<String> {
        b.molecules
            .iter()
            .map(|m| {
                let mut c = m.clone();
                sanitize(&mut c).unwrap();
                canon::canonical_smiles(&c).smiles
            })
            .collect()
    };
    assert_eq!(
        by_positional.verdict, by_graph.verdict,
        "两个入口的结论应当相同"
    );
    assert_eq!(
        names(&by_positional),
        names(&by_graph),
        "两个入口的副产物应当相同"
    );
    assert_eq!(names(&by_graph), ["Cl"]);
    // 契约:两个入口都按**输入分子**给下标,长度等于输入分子数
    assert_eq!(
        run_on_substrate(&rxn, &inputs, 0, false)[0].discarded.len(),
        2
    );
}

/// CDI 的两个咪唑要收得回来 —— 片段的键级得按**凯库勒式**数。
///
/// 底物与模板取自语料第 4635 行(`US20100105723A1`),不是造的。CDI
/// (羰基二咪唑)把羰基交出去,两个咪唑离去;它们各自是完整的芳香环。
///
/// 先前直接从芳香式的分子上切,芳香键与芳香标志被原样搬进片段,净化当场报
/// "原子不在环中却带着芳香标志",于是报成 `FragmentUnsanitizable` —— 那个理由
/// 是**错的**:不是收口路线不成立,是切的时候把标志搬错了地方。
///
/// 全量对拍过:这个修在 52261 个 outcome 里只改动 3 个,这是其中 2 个。
#[test]
fn an_aromatic_leaving_group_survives_the_cut() {
    let (_, byproducts, verdict) = run(
        "[#7;a:1]:[c:2](-[NH2;D1;+0:3]):[c:4]-[NH2;D1;+0:5].\
         [O;D1;H0:6]=[C;H0;D3;+0:7](-n1:c:c:n:c:1)-n1:c:c:n:c:1>>\
         [#7;a:1]:[c:2]1:[c:4]:[nH;D2;+0:5]:[c;H0;D3;+0:7](=[O;D1;H0:6]):[nH;D2;+0:3]:1",
        &[
            "Nc1cc(-c2c(F)cncc2F)c(-c2ccccc2F)nc1N",
            "O=C(n1ccnc1)n1ccnc1",
        ],
    );
    assert_eq!(byproducts, ["c1ncc[nH]1", "c1ncc[nH]1"], "CDI 的两个咪唑");
    assert_eq!(verdict, Verdict::Capped);
}

/// 收口不许把三键塞进小环 —— 那是几何上不可能的东西。
///
/// 底物与模板取自语料第 1446 行(`US20100279990A1`)。模板把一个碘原子交出去,
/// 而账要求片段再成一根键;唯一配得上的两处空价落在同一个苯环上,成出来的是
/// 苯炔式的 `c#c`。
///
/// **这一档躲得过前面每一道判据**:原子账平、电荷账平、净化也过得去。价规则管
/// 的是"几根键",管不到"这几根键摆得下摆不下"。不单独拦的话,输出的是一个
/// 配平、合法、下游任何检查都看不出问题的**错分子**。
///
/// 修 kekulize 那一处**顺带把这条路打开了**(修之前它报 `BudgetMismatch`),
/// 所以两处必须一起做 —— 只修前者等于把一条诚实的失败换成一个错答案。
#[test]
fn a_triple_bond_is_never_closed_into_a_small_ring() {
    let (_, byproducts, verdict) = run(
        "C-N-[C@@H]1-C-C-[C@@H](-c2:c:c:c(-Cl):c(-Cl):c:2)-c2:c:c:c(-[I;H0;D1;+0:1]):c:c:2-1\
         >>[IH;D0;+0:1]",
        &["CN[C@H]1CC[C@@H](c2ccc(Cl)c(Cl)c2)c2ccc(I)cc21"],
    );
    assert!(
        byproducts.is_empty(),
        "几何上不可能的东西不该交出去,给了 {byproducts:?}"
    );
    assert_eq!(verdict, Verdict::Unresolved(Unresolved::StrainedClosure));
}

/// 切点落在手性中心上时,标记要跟着换参照系。
///
/// 标记是相对邻居**存储顺序**说的。原子被切下来时那个邻居从列表里消失,顶上来的
/// 隐式氢按本库的约定占下标 1 —— 被切邻居原在下标 0 或 2 时这个置换是**奇**的,
/// 标记必须翻;在 1 或 3 时是偶的,不该翻。
///
/// # 判据必须同时抓住两个方向
///
/// 下面四条**故意各写各的邻居次序**:前两条的被切邻居落在需要翻的槽位,后两条
/// 落在不该翻的槽位。只测前两条的话,一个"一律翻转"的实现照样全绿;只测后两条,
/// 一个"一律不翻"的实现全绿 —— 而后者正是修之前的样子。
///
/// # 这一组是**拼装**的,不是从语料里取的
///
/// 语料给不出:两万条里副产物带手性标记的有 7 条,但它们的切点**全部**落在氧上
/// (酯的酰氧断裂),手性中心的邻居一个没少,这条路一次都走不到。所以这里按
/// 最小形状拼一个脱羧 —— 拼的是模板与底物,**判据本身不是拼的**:期望值由
/// RDKit 独立算出(`FragmentOnBonds` 切开、哑原子换氢),八种写法逐条对过。
#[test]
fn a_stereocentre_at_the_cut_is_rebased() {
    let tpl = "[C](-C)(-C-C)(-F)-[C:1](=[O:2])-[OH:3]>>[O:3]=[C:1]=[O:2]";
    // (底物, 期望的副产物) —— 期望值取自 RDKit 的独立计算
    for (substrate, want) in [
        // 被切邻居在下标 0:置换为奇,必须翻
        ("OC(=O)[C@](C)(CC)F", "CC[C@H](C)F"),
        ("OC(=O)[C@@](C)(CC)F", "CC[C@@H](C)F"),
        // 被切邻居在下标 3:置换为偶,不该翻
        ("CC[C@](C)(F)C(=O)O", "CC[C@H](C)F"),
        ("CC[C@@](C)(F)C(=O)O", "CC[C@@H](C)F"),
    ] {
        let (_, byproducts, verdict) = run(tpl, &[substrate]);
        assert_eq!(verdict, Verdict::Capped, "{substrate}");
        // 两侧都过一遍本库的规范化,免得比的是写法差异
        let want_canon = canon::canonical_smiles(&sanitized(want)).smiles;
        assert_eq!(byproducts, [want_canon], "{substrate} 的副产物构型不对");
    }
}

/// 配位键的**给体端不占价** —— 断开之后它不欠任何东西。
///
/// 电子对是给体自己出的,所以配位键对给体端的价贡献是 0(全库的唯一真相来源是
/// `BondData::valence_contribution_to`)。按对称的键级算会给它凭空记一处空价,
/// 于是一条本来完整的记录被判成收不平。
///
/// 判据盯的是 `open_valence` 而不是最终副产物:这一处修的正是**空价怎么数**,
/// 拿最终分子去比会把它与别的因素混在一起。
#[test]
fn the_donor_end_of_a_dative_bond_owes_nothing() {
    let rxn = smarts::parse_reaction("N(C)(C)->[B:1]>>[BH3:1]").unwrap();
    let mols = [sanitized("CN(C)->B")];
    let inputs: Vec<(MolBuilder, MolProps)> = mols
        .iter()
        .map(|m| (m.clone(), MolProps::compute(m)))
        .collect();
    let outs = run_reactants(&rxn, &inputs, 0, false);
    let by = byproduct::reconstruct(&mols, &outs[0]);
    assert_eq!(
        by.budget.open_valence, 0,
        "配位键断开,给体端不该欠价(按对称键级算会记成 1)"
    );
}

/// "配不成对"要与"键太多"分开报。
///
/// 甲苯脱甲基:丢下来的甲基只有**一个**位点,剩余空价再多也没有第二个位点跟它
/// 配对。这与"要成的键超过上限"是两回事 —— 后者是本实现主动不找了,前者是物理上
/// 配不起来。混成一个出口的话,归因脚本会把后者全算成"搜索爆掉"。
#[test]
fn valences_that_cannot_pair_up_say_so() {
    let (_, byproducts, verdict) = run("[c:1][CH3]>>[cH:1]", &["Cc1ccccc1"]);
    assert!(byproducts.is_empty());
    assert_eq!(verdict, Verdict::Unresolved(Unresolved::NoPairing));
}

/// 电荷对空价的作用**由元素定**,不是"负减正加"。
///
/// 两条走同一段代码,却要求相反的方向:
///
/// - 叔丁基氯离解:片段**得到正电荷**。按"正加一处价"算的话账多出两处,收不平;
///   而 C⁺ 的价是 3(不是 4),正电荷在这里**填掉**一处空价。
/// - 溴以阴离子离去:片段得到负电荷,Br⁻ 的价是 0,负电荷同样填掉一处。
///
/// 只测后者的话,"负减正加"那种写死的实现照样全绿 —— 而它在前者上是错的。
/// 真相来源是 `omgkit_chem::valence_shift`,与净化的隐式氢推断共用同一张价表。
#[test]
fn charge_changes_valence_by_element_not_by_sign() {
    // 片段拿到**正**电荷:碳正离子的价是 3,所以这一处电荷填掉一处空价
    let (_, byproducts, verdict) = run("[C:1][Cl:2]>>[Cl-:2]", &["CC(C)(C)Cl"]);
    assert_eq!(byproducts, ["C[C+](C)C"], "叔丁基正离子");
    assert_eq!(verdict, Verdict::Capped);

    // 片段拿到**负**电荷:方向相反,但同样是"填掉一处"
    let (_, byproducts, verdict) = run("[C:1][Br:2].[N:3]>>[C:1][N+:3]", &["CCBr", "CN(C)C"]);
    assert_eq!(byproducts, ["[Br-]"]);
    assert_eq!(verdict, Verdict::Capped);
}

/// 副产物不带底物的映射号。
///
/// 映射号是**模板内部**的东西 —— 它连的是模板两侧,与副产物无关。底物身上带号
/// 是常事(反应数据库导出的 SMILES 普遍带号),原样搬过来的话副产物会写成
/// `[OH2:10]` 而不是 `O`:分子是对的,而串不是。这类串进了下游会当成"带映射
/// 的反应",而那个号在这里没有任何一侧与之对应。
///
/// 别的用例的底物都不带号,所以这条路一次也走不到 —— 撤掉清零它们照样全绿。
#[test]
fn byproducts_do_not_inherit_substrate_atom_maps() {
    let (_, byproducts, verdict) = run(
        "[C:1](=[O:2])[OH:3].[OH:4][C:5]>>[C:1](=[O:2])[O:4][C:5]",
        &["[CH3:7][C:8](=[O:9])[OH:10]", "CCO"],
    );
    assert_eq!(byproducts, ["O"], "带号会写成 [OH2:10],分子对而串不对");
    assert_eq!(verdict, Verdict::Capped);
}

/// 一个原子都没丢的反应,结论是"没有副产物",不是"收口失败"。
///
/// 两者要分得开:前者是正常情形,后者要引起注意。混在一起的话,一个总是失败的
/// 实现看起来会像"大多数反应本来就没有副产物"。
#[test]
fn a_reaction_that_discards_nothing_reports_nothing() {
    let (_, byproducts, verdict) = run("[C:1][OH:2]>>[C:1][O:2]C", &["CCO"]);
    assert!(byproducts.is_empty());
    assert_eq!(verdict, Verdict::Nothing);
}

/// 不连通的旁观组分(盐的反离子)不算被丢弃 —— 它们**原样进产物**,
/// 由 `seed_spectators` 管,不该在这里再冒出来一次。
#[test]
fn a_spectator_counter_ion_is_not_a_byproduct() {
    let (products, byproducts, verdict) = run(
        "[C:1](=[O:2])[OH:3].[OH:4][C:5]>>[C:1](=[O:2])[O:4][C:5]",
        &["CC(=O)O.[Na+].[Cl-]", "CCO"],
    );
    assert_eq!(verdict, Verdict::Capped);
    assert_eq!(byproducts, ["O"]);
    // 反离子在产物侧,不在副产物侧
    assert!(products.iter().any(|p| p.contains("Na")));
    assert!(products.iter().any(|p| p.contains("Cl")));
}

/// 收不了口的时候**不给分子**。
///
/// 这里构造的是"氧被拿走了,而记录里没有任何东西给它配对"(还原剂没写)——
/// 丢掉的是一个欠着两处价的氧,预算给不出两个氢,剩余空价为奇。引擎必须说
/// 答不了,而不是编一个水或者双氧水出来。
#[test]
fn an_unbalanced_record_is_reported_not_guessed() {
    let (_, byproducts, verdict) = run(
        "[N+:1](=[O:2])[O-:3]>>[NH2:1]",
        &["Cc1ccc(cc1)[N+](=O)[O-]"],
    );
    assert!(
        byproducts.is_empty(),
        "收不了口就不该给分子,给了 {byproducts:?}"
    );
    assert!(
        matches!(verdict, Verdict::Unresolved(_)),
        "应当报未决,实际 {verdict:?}"
    );
}

/// 记录漏了供氢的试剂时,要报得**具体**,不能只说"配不成对"。
///
/// 硝基还原成胺:产物比底物多两个氢,而那两个氢来自还原剂,记录里没有它。
/// `delta_h` 就是副产物应有的氢数,这时它是负的 —— 物理上讲不通,所以是一条
/// 硬拒。单独成档是因为它指向的东西很具体(缺供氢试剂),混进"配不成对"里
/// 就只剩一句没信息的话。
#[test]
fn a_missing_hydrogen_donor_is_named_as_such() {
    let (_, byproducts, verdict) = run(
        "[N+:1](=[O:2])[O-:3]>>[NH2:1]",
        &["Cc1ccc(cc1)[N+](=O)[O-]"],
    );
    assert!(byproducts.is_empty());
    assert_eq!(
        verdict,
        Verdict::Unresolved(Unresolved::HydrogenBudgetNegative)
    );
}

/// 判据不能只看"给出来的分子对不对" —— 还要看**账**。
///
/// 这一条直接查 `Budget`:三个量互相牵制,任何一个算错都会让等式不成立,
/// 而分子本身可能照样"看着像水"。
#[test]
fn the_budget_is_internally_consistent() {
    let rxn =
        smarts::parse_reaction("[C:1](=[O:2])[OH:3].[OH:4][C:5]>>[C:1](=[O:2])[O:4][C:5]").unwrap();
    let mols = [sanitized("CC(=O)O"), sanitized("CCO")];
    let inputs: Vec<(MolBuilder, MolProps)> = mols
        .iter()
        .map(|m| (m.clone(), MolProps::compute(m)))
        .collect();
    let outs = run_reactants(&rxn, &inputs, 0, false);
    let by = byproduct::reconstruct(&mols, &outs[0]);

    let b = by.budget;
    assert_eq!(b.open_valence, 1, "只断了一根单键");
    assert_eq!(b.fragment_hydrogens, 1, "被删的羟基氧自带一个氢");
    assert_eq!(b.delta_h, 2, "水要两个氢");
    assert_eq!(b.need, 1, "还差一个");
    assert_eq!(b.remaining, 0, "补上就闭合,不用再成键");
    assert_eq!(
        b.remaining,
        i32::try_from(b.open_valence).unwrap() - b.need,
        "remaining 的定义"
    );
}

/// `discarded` 是**事实**,与收不收得了口无关。
///
/// 收口是推断,可能答不了;而"哪些原子没进产物"永远答得出来。两者混在一起的话,
/// 一条收不了口的反应会连"丢了哪些原子"都报不出来 —— 那是最需要这条信息的时候。
#[test]
fn discarded_is_recorded_even_when_closure_fails() {
    let rxn = smarts::parse_reaction("[N+:1](=[O:2])[O-:3]>>[NH2:1]").unwrap();
    let mols = [sanitized("Cc1ccc(cc1)[N+](=O)[O-]")];
    let inputs: Vec<(MolBuilder, MolProps)> = mols
        .iter()
        .map(|m| (m.clone(), MolProps::compute(m)))
        .collect();
    let outs = run_reactants(&rxn, &inputs, 0, false);
    let by = byproduct::reconstruct(&mols, &outs[0]);

    assert!(matches!(by.verdict, Verdict::Unresolved(_)));
    assert_eq!(
        outs[0].discarded[0].len(),
        2,
        "两个氧没进产物,这一点无论收不收得了口都成立"
    );
}

/// 产物净化不过时,氢预算无从算起 —— 要说清是这个原因,不能报成别的。
///
/// 产物的隐式氢是净化才填的派生量。拿一个没净化过的产物去算总氢,得到的数
/// 没有意义,据此收口出来的副产物会稳定地错,而且**看不出**错在哪。
#[test]
fn unsanitizable_products_are_reported_as_such() {
    // 五价碳:模板明写了一个填不满的东西,产物必然净化不过
    let rxn = smarts::parse_reaction("[C:1][OH:2]>>[C:1](C)(C)(C)(C)C").unwrap();
    let mols = [sanitized("CCO")];
    let inputs: Vec<(MolBuilder, MolProps)> = mols
        .iter()
        .map(|m| (m.clone(), MolProps::compute(m)))
        .collect();
    let outs = run_reactants(&rxn, &inputs, 0, false);
    let by = byproduct::reconstruct(&mols, &outs[0]);
    assert_eq!(
        by.verdict,
        Verdict::Unresolved(Unresolved::ProductsUnsanitizable)
    );
    assert!(by.molecules.is_empty());
}
