//! 产物生成的正确性。
//!
//! 判据是**产物的规范 SMILES 多重集** —— 产物原子的编号是构建顺序留下的
//! 痕迹,不是语义量。同一条反应在同一个底物上可能有多处反应位点,产物之间
//! 可能重复,所以比的是多重集而不是集合。

use omgkit_chem::sanitize;
use omgkit_core::MolBuilder;
use omgkit_io::{canon, smarts, smiles};
use omgkit_match::{run_reactants, MolProps};

/// 净化,再把方向键换算成双键自己的顺反属性。
///
/// 第二步不能省。方向依附在某根**单键**上,反应把那根键删掉,几何就跟着没了 ——
/// 哪怕双键本身根本没被碰过。感知成双键自己的属性之后,只要参照原子还在,信息
/// 就还在。这也是 `run_reactants` 的调用方该走的路子(净化那 12 步里没有它,
/// 理由见 `omgkit_io::stereo` 的模块文档)。
fn sanitized(smi: &str) -> MolBuilder {
    let mut m = smiles::parse(smi).unwrap_or_else(|e| panic!("{smi}: {}", e.render()));
    sanitize(&mut m).unwrap_or_else(|e| panic!("{smi}: {e}"));
    omgkit_io::stereo::perceive_bond_stereo(&mut m);
    m
}

/// 跑反应,把每组产物规范化成 SMILES;净化失败的原样标出来。
fn products(rxn_smarts: &str, reactant_smis: &[&str]) -> Vec<Vec<String>> {
    let rxn = smarts::parse_reaction(rxn_smarts)
        .unwrap_or_else(|e| panic!("{rxn_smarts}:\n{}", e.render()));
    let inputs: Vec<(MolBuilder, MolProps)> = reactant_smis
        .iter()
        .map(|s| {
            let m = sanitized(s);
            let p = MolProps::compute(&m);
            (m, p)
        })
        .collect();

    run_reactants(&rxn, &inputs, 0, false)
        .into_iter()
        .map(|outcome| {
            outcome
                .products
                .into_iter()
                .map(|mut p| match sanitize(&mut p) {
                    Ok(()) => canon::canonical_smiles(&p).smiles,
                    Err(e) => format!("<净化失败: {e}>"),
                })
                .collect()
        })
        .collect()
}

/// 同上,但产物走**调用方真正走的那条路**:净化之后再感知一次顺反。
///
/// `sanitize` 那 12 步里没有顺反感知(理由见 `omgkit_io::stereo` 的模块文档),
/// 所以 [`products`] 得到的产物身上只有沿用来的方向键,`stereo` 还没算。真实
/// 调用方不会停在那里 —— Python 绑定的 `Mol.sanitize` 就是净化完接着感知。
///
/// 差别不只是"多算一步":感知会**读**方向键。产物身上的方向键可能来自模板,
/// 而搬运过来的 `stereo` 来自底物,两者参照系不同。只有走完这一步,那类冲突
/// 才会显形 —— 停在 [`products`] 的判据对它一律空过。
fn perceived_products(rxn: &str, mols: &[&str]) -> Vec<String> {
    let rxn_parsed =
        smarts::parse_reaction(rxn).unwrap_or_else(|e| panic!("{rxn}:\n{}", e.render()));
    let inputs: Vec<(MolBuilder, MolProps)> = mols
        .iter()
        .map(|s| {
            let m = sanitized(s);
            let p = MolProps::compute(&m);
            (m, p)
        })
        .collect();
    let mut v: Vec<String> = run_reactants(&rxn_parsed, &inputs, 0, false)
        .into_iter()
        .flat_map(|outcome| outcome.products)
        .map(|mut p| match sanitize(&mut p) {
            Ok(()) => {
                omgkit_io::stereo::perceive_bond_stereo(&mut p);
                canon::canonical_smiles(&p).smiles
            }
            Err(e) => format!("<净化失败: {e}>"),
        })
        .collect();
    v.sort();
    v
}

/// 展平成"产物 SMILES 的多重集",排序后便于比对。
fn flat(rxn: &str, mols: &[&str]) -> Vec<String> {
    let mut v: Vec<String> = products(rxn, mols).into_iter().flatten().collect();
    v.sort();
    v
}

/// 最基本的一条:羟基换成氯。
#[test]
fn hydroxyl_to_chloride() {
    assert_eq!(flat("[C:1][OH:2]>>[C:1][Cl:2]", &["CCO"]), vec!["CCCl"]);
}

/// 底物上有两处反应位点就出两组产物 —— 即使产物相同。
///
/// 去重是调用方的事:关心"能得到哪些产物"时要去重,关心"有几条路径"时不能去。
#[test]
fn two_sites_give_two_product_sets() {
    let v = flat("[C:1][OH:2]>>[C:1][Cl:2]", &["OCCO"]);
    assert_eq!(v.len(), 2, "乙二醇的两个羟基各反应一次");
    assert_eq!(v[0], v[1], "两条路径给出同一个产物");
}

/// 模板之外的部分要原样带过来。
///
/// 模板只描述"羟基变氯",苯环、支链都得自动跟着走。这一条错了的话产物会
/// 只剩模板里画出来的那几个原子。
#[test]
fn untouched_parts_are_carried_over() {
    assert_eq!(
        flat("[C:1][OH:2]>>[C:1][Cl:2]", &["c1ccccc1CO"]),
        vec![canonical("ClCc1ccccc1")]
    );
    assert_eq!(
        flat("[C:1][OH:2]>>[C:1][Cl:2]", &["CC(C)(C)CO"]),
        vec![canonical("CC(C)(C)CCl")]
    );
}

/// 反应物模板里**没有映射号**的原子会被删掉。
#[test]
fn unmapped_reactant_atoms_are_deleted() {
    // OH 没有映射号 → 氧被删;产物侧的 Cl 没有映射号 → 新建
    assert_eq!(flat("[C:1][OH]>>[C:1]Cl", &["CCO"]), vec!["CCCl"]);
}

/// 键级的改写。
///
/// 对称模板会给出**两组**产物:`[C:1]=[C:2]` 配到 `C=C` 上有两种映射方向。
/// 内容相同不代表只有一条路径 —— 去重是调用方的事,见
/// [`two_sites_give_two_product_sets`]。
#[test]
fn bond_orders_are_rewritten() {
    assert_eq!(flat("[C:1]=[C:2]>>[C:1]-[C:2]", &["C=C"]), vec!["CC", "CC"]);
    assert_eq!(
        flat("[C:1]-[C:2]>>[C:1]=[C:2]", &["CC"]),
        vec!["C=C", "C=C"]
    );
}

/// 两个反应物模板:各自匹配,再取组合。
#[test]
fn two_reactant_templates() {
    let v = flat("[C:1][OH:2].[N:3]>>[C:1][N:3]", &["CO", "N"]);
    assert_eq!(v.len(), 1);
    assert_eq!(v[0], canonical("CN"));
}

/// 产物模板可以有多个,一次给出多个分子。
#[test]
fn multiple_product_templates() {
    let sets = products("[C:1][O:2][C:3]>>[C:1][O:2].[C:3]", &["COC"]);
    assert!(!sets.is_empty());
    assert_eq!(sets[0].len(), 2, "两个产物模板 → 每组两个分子");
}

/// 电荷与同位素这类**写死在产物模板里**的属性要盖上去;
/// 没写的属性保持继承。
#[test]
fn template_attributes_are_applied() {
    // 产物模板写了 +1,要盖上去
    let v = flat("[N:1]>>[N+:1]", &["CN"]);
    assert_eq!(v.len(), 1);
    assert!(v[0].contains('+'), "产物应带正电荷,实际 {}", v[0]);
}

/// 匹配不上就没有产物,而不是崩或者给出空分子。
#[test]
fn no_match_gives_no_products() {
    assert!(flat("[C:1][OH:2]>>[C:1][Cl:2]", &["CC"]).is_empty());
    // 反应物数目不符也一样
    assert!(flat("[C:1].[N:2]>>[C:1][N:2]", &["CC"]).is_empty());
}

fn canonical(smi: &str) -> String {
    let m = sanitized(smi);
    canon::canonical_smiles(&m).smiles
}

/// 重原子数,用来把"少了几个原子"写成数字而不是感觉。
fn heavy(smi: &str) -> usize {
    sanitized(smi).num_atoms()
}

/// 去重后的产物多重集,排序;判据要的是"输出了哪些",不是"输出了几次"。
fn distinct(rxn: &str, mols: &[&str]) -> Vec<String> {
    let mut v = flat(rxn, mols);
    v.dedup();
    v
}

/// 语料 US05026856(USPTO-50k 第 49427 行)的 N-脱苄模板。
///
/// `-C-c1:c:c:c:c:c:1` 那一支整支没有映射号,是要删掉的。
const DEBENZYL: &str =
    "[C:1]-[N;H0;D3;+0:2](-C-c1:c:c:c:c:c:1)-[C:3]>>[C:1]-[NH;D2;+0:2]-[C:3]";
/// 同一条记录的底物:N-苄基异吲哚啉,苯环上带一个甲氧基。
const DEBENZYL_SUB: &str = "COc1ccc2c(c1)CN(Cc1ccccc1)C2";

/// 产物只包含从**保留下来的原子走得到**的部分 —— 之一:模板删掉某个原子时,
/// 只挂在它身上的东西跟着一起走。
///
/// 上面这条真模板在这个真底物上有**两处**匹配,第二处正好把约定暴露出来:
///
/// - 一处:`-C-c1ccccc1` 那一支匹配到真正的苄基 → 脱苄,产物 11 个重原子
/// - 另一处:匹配到**稠合的苯环**那一支 → 苯环被删,挂在苯环上的甲氧基
///   再没有别的路连回保留部分,于是一并消失 —— 产物只剩 9 个重原子
///
/// 后者不是缺陷:模板说了要删那 7 个原子,而甲氧基是通过它们才连着的。换任何
/// 别的做法都得回答"这两个原子该接到哪里去",而那没有普遍答案。RDKit 在同一
/// 条模板同一个底物上给出**一模一样**的两组产物;全语料 1493 个 outcome 落在
/// 这一档,两个引擎的数字逐条相同。
#[test]
fn atoms_reachable_only_through_a_deleted_atom_go_with_it() {
    let got = distinct(DEBENZYL, &[DEBENZYL_SUB]);
    let mut want = vec![
        canonical("COc1ccc2c(c1)CNC2"), // 脱苄
        canonical("CNCc1ccccc1"),       // 苯环那一支被删,甲氧基跟着走
    ];
    want.sort();
    assert_eq!(got, want);

    // 把"少掉的是甲氧基那 2 个原子"写成数字:底物 18,模板删 7,若什么都不带走
    // 该剩 11;实际 9。
    assert_eq!(heavy(DEBENZYL_SUB), 18);
    assert_eq!(heavy("CNCc1ccccc1"), 9);
}

/// 与上一条相反的一档:**完全不连通的旁观组分要原样交回来**,不能丢。
///
/// 搬运是从匹配到的原子出发做遍历的,与任何匹配原子都不连通的组分本来没有路
/// 可走。曾经因此把它们丢掉 —— 而那是**引擎自己在丢原子**:产物的重原子数少于
/// 底物,不报错。逆合成正是把模板作用到任意分子上,盐是常态,所以这一条不做成
/// 开关,直接改成默认搬过来(见 `seed_spectators`)。
///
/// 与上一条的分界是"这个组分里**有没有**原子被模板匹配到":有,留下还是删掉就是
/// 模板的表态;没有,模板压根没提到它,不能替它做主。
///
/// # 底物为什么是拼的
///
/// 模板与主体底物都取自语料(同上一条),但**后面那个 HCl 是拼上去的**。
/// USPTO-50k 抽模板时把多组分物种拆成了独立分子,50016 条记录的输入分子全是
/// 单组分,旁观反离子一条都没有 —— 语料给不出这一档的底物。拼上去的是**书写
/// 形式**,不是化学:盐在真实数据里到处都是。
#[test]
fn a_component_with_no_matched_atom_comes_back_untouched() {
    let with_salt = format!("{DEBENZYL_SUB}.Cl");
    // 先证明 HCl 确实进了引擎 —— 否则这条判据可能是在"解析时就丢了"上空过,
    // 而那与搬运时走不到完全是两回事
    assert_eq!(heavy(&with_salt), heavy(DEBENZYL_SUB) + 1, "HCl 没被读进来");

    let plain = distinct(DEBENZYL, &[DEBENZYL_SUB]);
    let salted = distinct(DEBENZYL, &[&with_salt]);
    assert!(!plain.is_empty(), "对照组没有产物,判据空过");

    // 旁观组分作为自己的一个分子回来,别的产物一个字不变
    let mut want = plain.clone();
    want.push(canonical("Cl"));
    want.sort();
    assert_eq!(salted, want, "旁观的 HCl 没被原样交回来");

    // 逐 outcome 数重原子:这条模板会删掉 7 个原子,不新建;带不带 HCl,
    // 产物总数就该差**恰好一个** HCl 的重原子数。数字才是判据,SMILES 只是佐证。
    for (a, b) in products(DEBENZYL, &[DEBENZYL_SUB])
        .iter()
        .zip(products(DEBENZYL, &[&with_salt]).iter())
    {
        let ha: usize = a.iter().map(|s| heavy(s)).sum();
        let hb: usize = b.iter().map(|s| heavy(s)).sum();
        assert_eq!(hb, ha + 1, "带上 HCl 之后产物少了原子:{a:?} → {b:?}");
    }
}

/// 盐:两个反离子各自不连通,断键之后两个都要回来。
///
/// 底物是使用者给的一条实例。它把上一条判据推到更难的形状:**两个**旁观组分,
/// 而且反应把主体切成了两半 —— 引擎既不能丢掉反离子,也不该替使用者决定
/// 哪个反离子跟哪一半走(那没有普遍答案,模板里也没有这条信息)。
///
/// 判据是**质量守恒**:模板新建一个氧,所以产物重原子数 = 底物 + 1。
#[test]
fn both_counter_ions_of_a_salt_come_back() {
    const SALT: &str = "[Na+].[O-]CC(=O)OCC[O-].[K+]";
    const ESTER: &str = "[C:1](=[O:2])[O:3][C:4]>>[C:1](=[O:2])[O:3].[O:5][C:4]";

    let got = flat(ESTER, &[SALT]);
    assert!(!got.is_empty(), "一个产物都没有,判据空过");

    let total: usize = got.iter().map(|s| heavy(s)).sum();
    assert_eq!(
        total,
        heavy(SALT) + 1,
        "质量不守恒 —— 模板只新建一个氧,产物却是 {got:?}"
    );
    for ion in ["[Na+]", "[K+]"] {
        assert!(
            got.contains(&canonical(ion)),
            "{ion} 没回来:{got:?}"
        );
    }
}

/// 搬运未匹配部分时,键的**朝向**必须沿用源键,不能沿用遍历方向。
///
/// 遍历是从已保留的原子往外走的,走到某条键时的出发端与该键存储的 `begin`
/// 不一定是同一个原子。两者相反时若按遍历方向建键,朝向就翻了 —— 而朝向
/// 有语义:`direction`(`/` `\`)相对 `begin → end`,翻转即把顺式写成反式。
///
/// 这类错误不改变拓扑,只把分子换成几何异构体,肉眼极难发现,而带方向键的
/// 底物一多就成批出现。
#[test]
fn carried_over_bonds_keep_their_orientation() {
    // 反应只把羰基变成醇,C=N 那根双键与它两侧的方向键都只是被搬运
    let got = flat(
        "[C:1][C:2]=[O:3]>>[C:1][C:2][OH:3]",
        &["CC(=O)/C(Cl)=N/Nc1ccccc1"],
    );
    // 羰基碳两侧都能当 [C:1],所以有两处位点
    assert!(!got.is_empty(), "一个产物都没有");
    for p in &got {
        assert!(
            p.contains('/') || p.contains('\\'),
            "产物 {p} 丢了双键立体 —— 搬运时朝向翻了,方向键就成了噪声被丢弃"
        );
    }
}

/// 配位键的箭头同样靠朝向表达,搬运时不能翻。
#[test]
fn carried_over_dative_bonds_keep_their_arrow() {
    let rxn = "[C:1][OH:2]>>[C:1][Cl:2]";
    let got = flat(rxn, &["OCC[N+](C)(C)C"]);
    assert_eq!(got.len(), 1);
    // 拓扑对齐即可 —— 这里要守的是"搬运不改朝向",不是具体写法
    assert!(got[0].contains("Cl"), "产物 {}", got[0]);
}

/// 断**环**键给出的是一个开环分子,不是两片。
///
/// `[C:1][N:2]>>[C:1].[N:2]` 作用在哌啶上:断掉 C—N 之后,环上其余的原子仍把
/// 两端连着,所以结果是**一个**开环分子。产物模板有两个,产物分子只有一个 ——
/// 分子数由连通性决定,不由模板数决定。
///
/// 这条最要紧的判据是**质量守恒**:逐产物各搬一次"模板之外的部分",环上那几个
/// 原子会被复制进两片,原子凭空变多而没有任何东西报错。
#[test]
fn breaking_a_ring_bond_gives_one_ring_opened_molecule() {
    // 按**组**看:哌啶上有多处 C—N 可匹配,每处给出一组产物
    let sets = products("[C:1][N:2]>>[C:1].[N:2]", &["C1CCNCC1"]);
    assert!(!sets.is_empty(), "一组产物都没有");
    for set in &sets {
        assert_eq!(set.len(), 1, "环上还连着,每组应当只有一个分子,实际 {set:?}");
        let p = &set[0];
        let heavy = p.chars().filter(|c| *c == 'C' || *c == 'N').count();
        assert_eq!(
            heavy, 6,
            "哌啶 6 个重原子,产物 {p} 有 {heavy} 个 —— 质量不守恒"
        );
        assert!(
            !p.contains('1') && !p.contains('2'),
            "产物 {p} 里还有环闭合标号 —— 环没断开"
        );
    }
}

/// 断链键的老行为不能被上面那条改坏。
///
/// 两个产物模板给出两个分子,`flat` 把它们展平成两条 —— 不是一条带 `.` 的串。
#[test]
fn breaking_a_chain_bond_still_splits_cleanly() {
    assert_eq!(flat("[C:1][N:2]>>[C:1].[N:2]", &["CCN"]), vec!["CC", "N"]);
}

/// 反应物模板**没匹配到**的键,产物模板无权删。
///
/// 子结构匹配只要求模板的每根键在底物里找得到,**不**要求这些原子之间没有
/// 别的键。模板把环写成一条**开链路径**时,环闭合的那根键两端确实都被匹配了,
/// 可模板从没看见它。判据若只看"两端都匹配就不搬",这根键就被当成模板的
/// 地盘删掉,环被撕开。
///
/// 撕开之后报出来的是**芳香**错误(原子不在环中却带着芳香标志),病因与症状
/// 隔着一层,极容易把它误记成"芳香标志没清干净"。
///
/// # 用语料里的真模板,不自己造
///
/// 下面这条是 rdchiral 从 USPTO-50k 的 US05965555(语料第 5351 行)抽出的
/// 逆向模板,底物与真值也照抄该记录 —— 嘌呤的 N-脱烷基。模板两侧都把稠环
/// 写成六个原子的开链路径,连接两端的那根键不在模板里。
///
/// 自己造的模板可能根本不是它看上去的意思(产物侧不写键级时该继承什么、
/// 芳香原子之间的默认键级是什么,都不是想当然的),拿它当判据会把判据本身
/// 变成错的。全语料扫下来,这一档正向 11 条、逆向 67 条。
#[test]
fn a_ring_closure_the_template_never_matched_survives() {
    let rxn = "[#7;a:3]:[c:4]:[c:5]:[n;H0;D3;+0:6](-[CH2;D2;+0:1]-[C:2]):[c:7]:[#7;a:8]\
               >>Cl-[CH2;D2;+0:1]-[C:2].[#7;a:3]:[c:4]:[c:5]:[nH;D2;+0:6]:[c:7]:[#7;a:8]";
    let substrate = "CCCn1cnc2c1c(=O)n(CCCCC(C)(O)C#CCN(CC)CC)c(=O)n2C";
    let sets = products(rxn, &[substrate]);
    assert!(!sets.is_empty(), "一组产物都没有");
    for set in &sets {
        for p in set {
            assert!(!p.starts_with("<净化失败"), "产物净化不过:{p}");
        }
    }
    // 记录里的真值:氯代物 + 完整的嘌呤
    let mut want = vec![
        canonical("CCN(CC)CC#CC(C)(O)CCCCCl"),
        canonical("CCCn1cnc2c1c(=O)[nH]c(=O)n2C"),
    ];
    want.sort();
    assert!(
        sets.iter().any(|s| {
            let mut got = s.clone();
            got.sort();
            got == want
        }),
        "没有一组给出记录里的反应物 {want:?},实际 {sets:?}"
    );
}

/// 取代基被**替换**时,新原子占据被替换者的空间位置,手性几何不变。
///
/// `[C:1][OH]>>[C:1]Cl` 删掉氧、新建一个氯。中心的度数没变,变的是邻居身份。
/// 重定基时若把没进产物的邻居直接丢掉,反应物侧就比产物侧少一个,长度对不上
/// 整个重定基被跳过 —— 标记原样照抄,而产物侧的邻居顺序已经变了,手性就反了。
///
/// 判据:结果应当与"把 O 直接写成 Cl"的分子相同。
#[test]
fn replaced_substituents_keep_the_geometry() {
    for (substrate, expected) in [
        ("[C@H](N)(O)F", "[C@H](N)(Cl)F"),
        ("N[C@@H](O)F", "N[C@@H](Cl)F"),
        ("O[C@H](N)F", "Cl[C@H](N)F"),
        ("C[C@H](O)F", "C[C@H](Cl)F"),
    ] {
        let got = flat("[C:1][OH]>>[C:1]Cl", &[substrate]);
        assert_eq!(
            got,
            vec![canonical(expected)],
            "{substrate}:羟基换成氯之后手性该保持,不该翻"
        );
    }
}

/// 取代基被**隐式氢**接管时,构型同样要换参照系。
///
/// 脱保护、脱羧、脱卤都是这个形状:中心从 D4H0 变成 D3H1,少了一个图上的邻居。
/// 氢不是图里的节点,可它占着四面体的一个位置 —— 因为"两侧邻居个数对不上"
/// 就跳过重定基的话,标记会原样留在**反应物的**参照系里,而产物的邻居顺序已经
/// 变了,于是拿到镜像。
///
/// 这一档极难发现:拓扑、原子数、电荷、分子式全对,只有构型反了。
/// 反面(反应物侧的隐式氢被新邻居顶替)一并守在下面一条。
#[test]
fn a_substituent_replaced_by_a_hydrogen_still_rebases() {
    // Boc 保护的季碳脱去 Boc:两侧都不写手性,走的是继承 + 重定基那条路
    let rxn = "C-C(-C)(-C)-O-C(=O)-[C;H0;D4;+0:1](-[C:2])(-[N;D1;H2:3])-[C:4](-[#7:5])=[O;D1;H0:6]\
               >>[#7:5]-[C:4](=[O;D1;H0:6])-[CH;D3;+0:1](-[C:2])-[N;D1;H2:3]";
    // 叔丁基的三个甲基都配得上 `C-C(-C)(-C)`,同一个产物会出现好几组
    let got = flat(rxn, &["CNC(=O)[C@@](N)(CC)C(=O)OC(C)(C)C"]);
    assert!(!got.is_empty(), "一个产物都没有");
    for p in &got {
        assert_eq!(
            *p,
            canonical("CC[C@H](N)C(=O)NC"),
            "取代基换成隐式氢之后构型反了 —— 重定基被长度对不上挡掉了"
        );
    }
}

/// 反面:反应物侧的**隐式氢**被新邻居顶替时,产物不能随模板的书写顺序而变。
///
/// 新邻居占的正是氢原来的位置。模板把同样这几个邻居换个次序写,描述的仍是同一个
/// 产物 —— 结果必须一样。这条判据比"结果等于某个值"更强,它同时守住两件事:
/// 换参照系这一步**真的做了**,而且**没做过头**。
///
/// 下面六种写法里前三种撤掉那一支就会给出镜像(枚举全部 24 种写法时是 12 种),
/// 后三种撤掉也照样对 —— 只挑后三种当判据就会空过。
#[test]
fn the_product_does_not_depend_on_how_the_template_orders_neighbours() {
    for tail in [
        "[N:3]-[C;H0;D4;+0:1](-[C:2])(-[O:4])-Cl",
        "[C:2]-[C;H0;D4;+0:1](-[N:3])(-Cl)-[O:4]",
        "Cl-[C;H0;D4;+0:1](-[O:4])(-[C:2])-[N:3]",
        "[C:2]-[C;H0;D4;+0:1](-[N:3])(-[O:4])-Cl",
        "[N:3]-[C;H0;D4;+0:1](-[O:4])(-[C:2])-Cl",
        "[O:4]-[C;H0;D4;+0:1](-Cl)(-[C:2])-[N:3]",
    ] {
        let rxn = format!("[C:2]-[CH;D3;+0:1](-[N:3])-[O:4]>>{tail}");
        assert_eq!(
            flat(&rxn, &["C[C@H](N)O"]),
            vec![canonical("C[C@](N)(O)Cl")],
            "{tail}:同一个产物,换个书写次序就变了构型"
        );
    }
}

/// 邻居**挪到别的产物片段**去了,也算腾出了槽位。
///
/// `[C@@H:1]-[O:2]>>C-C(=O)-O-[C@@H:1].[OH:2]`:氧带着映射号活到了另一个产物里,
/// 可它不再连着 `:1` 了,同时有个新的氧接了上来。判"这个槽位空没空"只看
/// "那个原子进没进产物"的话,槽位会被算成占着 —— 于是反应物侧的邻居集合里
/// 有个产物侧没有的原子,置换的多重集对不上,重定基被**静默跳过**。
///
/// 跳过的后果不是报错,是标记原样照抄,而产物的邻居顺序早变了 —— 得到镜像。
/// 反应物原子数、产物原子数、连通性全对,只有构型反了。
///
/// # 判据必须是**绝对**的
///
/// 不能拿"模板不写手性"那一档当参照:它走的是同一段重定基,同一个缺陷把两边
/// 一起带偏,比出来永远相等。第一版就这么空过了 —— 缺陷还在,测试却是绿的。
///
/// 所以直接写出应得的产物。醚氧换成乙酸酯氧、构型保留:
/// `CO[C@@H](C)CC` → `CC(=O)O[C@@H](C)CC`。两串里那个中心的书写序都是
/// (氧, 氢, 甲基, 乙基),标记的含义没变,肉眼就能核对。
#[test]
fn a_neighbour_moving_to_another_fragment_frees_its_slot() {
    // `[O;H0;D2]` 钉住醚氧 —— 不限定的话会配到别的氧上,多出来的产物组会把
    // 判据搅浑
    let retain = "[C@@H;D3:1]-[O;H0;D2:2]>>C-C(=O)-O-[C@@H;D3:1].[OH:2]";
    let invert = "[C@@H;D3:1]-[O;H0;D2:2]>>C-C(=O)-O-[C@H;D3:1].[OH:2]";

    // 中心的四个取代基:氧、氢、甲基、乙基 —— 都不同,是真中心。
    //
    // 醚氧要写在**中间**那一位:走掉的氧与接上来的氧若正好占同一个位置,置换是
    // 偶的,漏掉重定基也看不出来 —— `CO[C@@H](C)CC` 就是这样,拿它当用例是空过的。
    for (substrate, kept, flipped) in [
        ("C[C@@H](OC)CC", "C[C@@H](OC(C)=O)CC", "C[C@H](OC(C)=O)CC"),
        ("C[C@H](OC)CC", "C[C@H](OC(C)=O)CC", "C[C@@H](OC(C)=O)CC"),
    ] {
        // 只比带立体中心的那一片。另一片是甲醇,它的氢怎么写由模板的 `[OH:2]`
        // 定(写死氢数即不再补隐式氢),与这一档要守的东西无关。
        let got = flat(retain, &[substrate]);
        assert_eq!(got.len(), 2, "{substrate}:该出两个产物,实得 {got:?}");
        assert!(
            got.contains(&canonical(kept)),
            "{substrate}:两侧标记相同该保留构型,应当得到 {kept},实得 {got:?} —— \
             挪到另一个产物片段去的那个氧若还把槽位占着,重定基会被静默跳过,\
             得到的是镜像"
        );

        let got = flat(invert, &[substrate]);
        assert!(
            got.contains(&canonical(flipped)),
            "{substrate}:两侧标记不同该翻转构型,应当得到 {flipped},实得 {got:?}"
        );
    }
}

/// 顺反的参照原子被反应动掉时,几何要靠**顶替者**继续描述。
///
/// 参照原子有两种走法,判据都不能只看"那个原子进没进产物":
///
/// - 被**删掉**(模板里没有映射号)
/// - 活着,却被挪到了**别的产物片段**去 —— `kept` 里查得到,于是参照被换成一个
///   属于另一个分子的下标。切分之后那个下标在本分子里要么越界(几何静默丢失),
///   要么正好落在某个真邻居上(几何**静默变错**)
///
/// 两种都不报错。判据是一对顺反底物:产物必须跟着变,而且要与外部实现一致。
#[test]
fn bond_stereo_survives_losing_its_reference_atom() {
    // 醚氧承载着 `/`,反应把它换成碘。第一条把氧挪进另一个产物,第二条直接删掉。
    for rxn in [
        "[O;H0;D2:1]-[C;D3:2]>>[OH:1].I-[C;D3:2]",
        "[O;H0;D2]-[C;D3:2]>>I-[C;D3:2]",
    ] {
        let cis = flat(rxn, &["CO/C(F)=C(/Cl)C"]);
        let trans = flat(rxn, &["CO/C(F)=C(\\Cl)C"]);
        for (got, src) in [(&cis, "顺"), (&trans, "反")] {
            assert!(
                got.iter().any(|p| p.contains('/') || p.contains('\\')),
                "{rxn} 在{src}式底物上给出 {got:?} —— 双键没被碰过,几何却丢了"
            );
        }
        // 顺与反必须给出不同的产物 —— 否则上面那条是空过的:
        // 两边都丢了几何时,它们同样"都不含斜杠",而这里会双双通过
        assert_ne!(
            cis, trans,
            "{rxn}:一对顺反底物给出了同一批产物,几何被抹平了"
        );
    }
}

// ---------------------------------------------------------------------------
// 原子映射号
// ---------------------------------------------------------------------------

/// 开启映射号跑反应。
///
/// 这里**不净化也不规范化**:两者都会重排原子,而这一档要看的正是号贴在哪个
/// 原子上,重排之后就对不上了。
fn run_mapped(rxn_smarts: &str, reactant_smis: &[&str]) -> Vec<omgkit_match::Outcome> {
    let rxn = smarts::parse_reaction(rxn_smarts)
        .unwrap_or_else(|e| panic!("{rxn_smarts}:\n{}", e.render()));
    let inputs: Vec<(MolBuilder, MolProps)> = reactant_smis
        .iter()
        .map(|s| {
            let m = sanitized(s);
            let p = MolProps::compute(&m);
            (m, p)
        })
        .collect();
    run_reactants(&rxn, &inputs, 0, true)
}

/// `映射号 → 原子序数`,把一侧的若干分子合起来看。
///
/// 顺带断言**同一个号在一侧只出现一次** —— 出现两次就不再是映射,而这正是
/// "一个反应物原子被复制进多个产物"时会犯的错。
fn by_map(mols: &[MolBuilder]) -> std::collections::BTreeMap<u16, u8> {
    let mut out = std::collections::BTreeMap::new();
    for m in mols {
        for a in m.atoms() {
            if a.atom_map != 0 {
                assert!(
                    out.insert(a.atom_map, a.atomic_num).is_none(),
                    "映射号 {} 在同一侧出现了两次",
                    a.atom_map
                );
            }
        }
    }
    out
}

/// 两侧的号必须一一配对 —— 这是"映射"这个词的全部含义。
#[test]
fn atom_maps_pair_up_across_the_arrow() {
    for (rxn, subs) in [
        ("[C:1][OH:2]>>[C:1][Cl:2]", &["CCO"][..]),
        ("[C:1][OH:2].[N:3]>>[C:1][N:3]", &["CO", "N"][..]),
        (
            "[C:1][C:2]=[O:3]>>[C:1][C:2][OH:3]",
            &["CC(=O)Cc1ccccc1"][..],
        ),
        ("[C:1][N:2]>>[C:1].[N:2]", &["CCN"][..]),
    ] {
        let outs = run_mapped(rxn, subs);
        assert!(!outs.is_empty(), "{rxn}:一组产物都没有");
        for o in &outs {
            let left = by_map(&o.reactants);
            let right = by_map(&o.products);
            assert!(!left.is_empty(), "{rxn}:一个号都没发出来");
            assert_eq!(
                left.keys().copied().collect::<Vec<_>>(),
                right.keys().copied().collect::<Vec<_>>(),
                "{rxn}:两侧的号对不上"
            );
        }
    }
}

/// 号不只要配对,还得贴在**真正对应**的那个原子上。
///
/// `[C:1][OH:2]>>[C:1][Cl:2]` 把氧就地改成氯:带同一个号的两个原子,左边是氧
/// 右边就该是氯。光比号的集合发现不了"号发给了别的原子"。
#[test]
fn atom_maps_point_at_the_atoms_that_really_correspond() {
    let outs = run_mapped("[C:1][OH:2]>>[C:1][Cl:2]", &["CCO"]);
    assert_eq!(outs.len(), 1);
    let left = by_map(&outs[0].reactants);
    let right = by_map(&outs[0].products);
    assert_eq!(left.len(), 3, "乙醇三个重原子都该有号,实际 {left:?}");

    // 元素没被模板改的原子,两侧元素必须一致
    for (n, &lz) in &left {
        let rz = right[n];
        if lz != 8 {
            assert_eq!(lz, rz, "号 {n} 两侧元素不同:{lz} vs {rz}");
        }
    }
    // 唯一被改写的那个:左边是氧,右边是氯
    let oxygen = left
        .iter()
        .find(|(_, &z)| z == 8)
        .map(|(n, _)| *n)
        .expect("反应物侧应有一个带号的氧");
    assert_eq!(right[&oxygen], 17, "号 {oxygen} 那个氧应当变成氯");
}

/// 只在一侧出现的原子不发号。
///
/// 发了的话写出来就是一条两侧对不上的反应:左边有 `:5` 右边没有,读的人只能
/// 理解成"这个原子凭空消失",而那正是号要表达的反面。
#[test]
fn atoms_without_a_counterpart_get_no_number() {
    // OH 在反应物模板里没有映射号 → 氧被删;产物侧的 Cl 是新建的
    let outs = run_mapped("[C:1][OH]>>[C:1]Cl", &["CCO"]);
    assert_eq!(outs.len(), 1);
    let o = &outs[0];
    let left = by_map(&o.reactants);
    assert_eq!(left.len(), 2, "只有两个碳有对应者,实际 {left:?}");
    assert!(left.values().all(|&z| z == 6), "带号的都该是碳");

    let oxygen = o.reactants[0]
        .atoms()
        .iter()
        .find(|a| a.atomic_num == 8)
        .expect("反应物里有氧");
    assert_eq!(oxygen.atom_map, 0, "被删掉的氧不该有号");

    let chlorine = o.products[0]
        .atoms()
        .iter()
        .find(|a| a.atomic_num == 17)
        .expect("产物里有氯");
    assert_eq!(chlorine.atom_map, 0, "新建的氯不该有号");
}

/// 跨反应物的号也不能撞 —— 一条反应 SMILES 是一个整体。
#[test]
fn numbers_are_unique_across_reactants() {
    let outs = run_mapped("[C:1][OH:2].[N:3]>>[C:1][N:3]", &["CO", "N"]);
    assert_eq!(outs.len(), 1);
    let o = &outs[0];
    assert_eq!(o.reactants.len(), 2, "两个反应物都要还回来");
    // by_map 会在号重复时直接失败
    let left = by_map(&o.reactants);
    assert_eq!(left.len(), 2, "碳和氮各一个号,实际 {left:?}");
    for (i, m) in o.reactants.iter().enumerate() {
        assert!(
            m.atoms().iter().any(|a| a.atom_map != 0),
            "第 {i} 个反应物一个号都没拿到"
        );
    }
}

/// 底物自带的映射号会被清掉重发。
///
/// 那些号与本次反应无关,留着就会出现"左边有号、右边找不到"的悬空号。
#[test]
fn map_numbers_in_the_input_are_replaced() {
    let outs = run_mapped("[C:1][OH:2]>>[C:1][Cl:2]", &["[CH3:7]CO"]);
    assert_eq!(outs.len(), 1);
    let left = by_map(&outs[0].reactants);
    assert_eq!(
        left.keys().copied().collect::<Vec<_>>(),
        vec![1, 2, 3],
        "应当重发成连续的 1..3,实际 {left:?}"
    );
}

/// 关掉时:产物不带号,反应物一份都不复制。
#[test]
fn mapping_off_leaves_no_numbers_and_no_reactant_copies() {
    let rxn = smarts::parse_reaction("[C:1][OH:2]>>[C:1][Cl:2]").unwrap();
    let m = sanitized("CCO");
    let inputs = vec![(m.clone(), MolProps::compute(&m))];
    let outs = run_reactants(&rxn, &inputs, 0, false);
    assert_eq!(outs.len(), 1);
    assert!(outs[0].reactants.is_empty(), "关掉时不该复制反应物");
    for p in &outs[0].products {
        assert!(
            p.atoms().iter().all(|a| a.atom_map == 0),
            "关掉时产物不该带号"
        );
    }
}

/// 端到端:写出来就是一条带原子映射号的反应 SMILES,再读回去号还在。
#[test]
fn mapped_reaction_survives_a_write_read_roundtrip() {
    let outs = run_mapped("[C:1][OH:2]>>[C:1][Cl:2]", &["CCO"]);
    let o = &outs[0];
    let side = |ms: &[MolBuilder]| {
        ms.iter()
            .map(|m| smiles::write(m).smiles)
            .collect::<Vec<_>>()
            .join(".")
    };
    let (left, right) = (side(&o.reactants), side(&o.products));

    // 号只到 3,`:1` 不会误配上 `:12` 这种更长的号
    for n in [":1", ":2", ":3"] {
        assert_eq!(left.matches(n).count(), 1, "{left}:{n} 应恰好出现一次");
        assert_eq!(right.matches(n).count(), 1, "{right}:{n} 应恰好出现一次");
    }

    let back = smiles::parse(&left).unwrap_or_else(|e| panic!("{left} 读不回:{}", e.render()));
    let mut got: Vec<u16> = back
        .atoms()
        .iter()
        .map(|a| a.atom_map)
        .filter(|&n| n != 0)
        .collect();
    got.sort_unstable();
    assert_eq!(got, vec![1, 2, 3], "往返之后号丢了:{left}");
}

/// 断环键时质量要守恒,而且每个原子恰好拿一个映射号。
///
/// 这里守的是"模板之外的部分只搬一次"。逐产物各搬一次的话,环上那几个原子会被
/// 复制进两片:原子数凭空变多,同一个原子在产物侧出现两次,映射号也就不再是映射。
#[test]
fn breaking_a_ring_bond_conserves_mass_and_maps_every_atom_once() {
    let outs = run_mapped("[C:1][N:2]>>[C:1].[N:2]", &["C1CCNCC1"]);
    assert!(!outs.is_empty(), "一组产物都没有");
    let o = &outs[0];
    assert_eq!(o.products.len(), 1, "环上还连着,应当只有一个产物分子");

    let total: usize = o.products.iter().map(MolBuilder::num_atoms).sum();
    assert_eq!(
        total,
        o.reactants[0].num_atoms(),
        "产物共 {total} 原子,底物 {} 个 —— 质量不守恒",
        o.reactants[0].num_atoms()
    );

    // by_map 内部会在号重复时失败
    let left = by_map(&o.reactants);
    let right = by_map(&o.products);
    assert_eq!(
        left.keys().copied().collect::<Vec<_>>(),
        right.keys().copied().collect::<Vec<_>>(),
        "两侧的号对不上"
    );
    assert_eq!(left.len(), 6, "哌啶六个重原子都该有号,实际 {left:?}");
}

/// 产物模板里写的方向键(`/` `\`)要落到产物上。
///
/// `>>C/[C:1]=[C:2]/C` 说的是"新生成的双键是反式" —— 那是作者对产物几何的
/// 明确指定,不是可有可无的书写痕迹。丢掉它,产物就从确定的顺反退化成未指定,
/// 而下游拿到"未指定"不会报错,只会当成两种几何都行。
///
/// 键级与方向是**两件事**:`/` 既说"这是单键"也说"在双键的哪一侧"。
/// 只取键级的话拓扑完全正确,只有立体没了。
#[test]
fn directions_written_in_the_product_template_are_applied() {
    for (rxn, want) in [
        ("[C:1]=[C:2]>>C/[C:1]=[C:2]/C", "C/C=C/C"),
        ("[C:1]=[C:2]>>C/[C:1]=[C:2]\\C", "C/C=C\\C"),
    ] {
        let got = flat(rxn, &["C=C"]);
        assert!(!got.is_empty(), "{rxn}:一个产物都没有");
        assert_eq!(got[0], canonical(want), "{rxn}:模板指定的双键几何没被应用");
    }
}

/// 模板没写方向时,仍然沿用反应物那边继承来的方向。
///
/// 两个来源有优先级:模板写了就听模板的,没写才继承。这条守的是"加了模板
/// 优先之后,继承那一支没被改坏"。
#[test]
fn inherited_directions_still_work_when_the_template_is_silent() {
    let got = flat(
        "[C:1][C:2]=[O:3]>>[C:1][C:2][OH:3]",
        &["CC(=O)/C(Cl)=N/Nc1ccccc1"],
    );
    assert!(!got.is_empty(), "一个产物都没有");
    for p in &got {
        assert!(
            p.contains('/') || p.contains('\\'),
            "产物 {p} 丢了继承来的双键立体"
        );
    }
}

/// 模板里**孤零零**一根方向键定不了任何几何,不能拿它去顶掉继承来的方向。
///
/// 顺反要双键两端各一根方向键才定得下来。`F/C=CF` 里那根 `/` 只说了"这根单键
/// 画在右上",另一端没画,取代基的相对位置仍然未知。
///
/// 危险在于:孤立的那根抄进产物之后并不孤立 —— 双键另一侧的键是从底物继承的,
/// 带着底物的方向。两根凑成一对,一根来自模板的书写顺序、一根来自底物的书写
/// 顺序,几何就成了两个参照系拼出来的任意结果。
///
/// 下面这条模板是 rdchiral 从 USPTO-50k 抽出来的酯水解,原样照抄:两侧各写一根
/// 孤立方向键,而且书写朝向恰好相反(反应物侧 `[C:2]/[C:4]=`,产物侧
/// `=[C:4]/[C:2]`)。照抄的话每一条都会翻一次,实测正向有 100 条以上栽在这里。
#[test]
fn a_lone_direction_in_the_template_does_not_override_the_substrate() {
    let rxn = "C-C-[O;H0;D2;+0:1]-[C:2](=[O;D1;H0:3])/[C:4]=[C:5]\
               >>[C:5]=[C:4]/[C:2](=[O;D1;H0:3])-[OH;D1;+0:1]";
    for (sub, want) in [
        ("CCOC(=O)/C=C/c1ccccc1", "OC(=O)/C=C/c1ccccc1"),
        ("CCOC(=O)/C=C\\c1ccccc1", "OC(=O)/C=C\\c1ccccc1"),
    ] {
        let got = flat(rxn, &[sub]);
        assert_eq!(
            got,
            vec![canonical(want)],
            "{sub}:模板那根孤立方向键把底物的几何顶掉了"
        );
    }
}

/// 共轭链上,模板写的方向键**同时贴着两根双键**,不能让它去重算另一根。
///
/// 模板给 `[C:3]=[C:4]` 写了一对方向,其中 `[C:4]-[C:5]` 那根同时是**下一根**
/// 双键的侧翼,而下一根的另一侧是从底物搬过来的。搬运用的是双键自己的参照原子
/// (与写法无关、天然正确),模板那根用的是模板的书写顺序 —— 两个参照系。
///
/// 拿它们凑成一对去感知下一根双键,就把一个正确答案覆盖成任意值。判据要走
/// **调用方真正走的那条路**:净化之后再感知一次顺反(Python 绑定的 `sanitize`
/// 就是这么做的),否则这一步压根不会执行,判据空过。
///
/// 实测 USPTO-50k 上 12 条,全是共轭多烯的酰胺化/酯化。
#[test]
fn a_template_direction_does_not_reach_into_the_neighbouring_double_bond() {
    let rxn = "O-[C;H0;D3;+0:1](=[O;D1;H0:2])/[C:3]=[C:4]/[C:5].[C:6]-[NH2;D1;+0:7]\
               >>[C:6]-[NH;D2;+0:7]-[C;H0;D3;+0:1](=[O;D1;H0:2])/[C:3]=[C:4]/[C:5]";
    // 两条底物在**模板管着的**那根双键(挨着羧基的)上都是反式 —— 模板要求如此;
    // 区别在**隔壁**那根:一条反式一条顺式,而它归继承管。
    for (diene, want) in [
        ("C/C=C/C=C/C(=O)O", "C/C=C/C=C/C(=O)NCC"),
        ("C/C=C\\C=C\\C(=O)O", "C/C=C\\C=C\\C(=O)NCC"),
    ] {
        let got = perceived_products(rxn, &[diene, "CCN"]);
        assert_eq!(
            got,
            vec![canonical(want)],
            "{diene}:模板那根方向键伸到隔壁双键上去了"
        );
    }
}

/// 反面:模板把几何**写全**了的时候,它说了算 —— 哪怕与底物相反。
///
/// 守的是上一条别修过头。"孤立的不算数"不能滑成"模板写的都不算数":产物模板
/// 双键两端各写一根方向键时,那是作者对产物几何的明确指定,底物继承来的要让位。
#[test]
fn a_determined_pair_in_the_template_still_overrides_the_substrate() {
    // 丁烯两端等价,模板在同一个底物上配得上两次,两组产物相同
    let got = flat(
        "[C:1][C:2]=[C:3][C:4]>>[C:1]/[C:2]=[C:3]\\[C:4]",
        &["C/C=C/C"],
    );
    assert!(!got.is_empty(), "一个产物都没有");
    for p in &got {
        assert_eq!(
            *p,
            canonical("C/C=C\\C"),
            "模板两端都写了方向,产物却还是底物的几何"
        );
    }
}

/// 反应物模板写了手性、产物模板没写 —— 那是模板在说"这个中心的构型不再确定"。
///
/// 反应物侧写 `[C@;...]`、产物侧写 `[C;...]` 是有意的不对称:作者在同一个映射号
/// 上先要求了构型、又不再给出构型。继承着不放,产物就凭空多出一个原始底物才有的
/// 构型 —— 而下游不会怀疑它。
///
/// 反面见 [`replaced_substituents_keep_the_geometry`]:两侧都没写手性时,
/// 模板压根没管这件事,继承来的构型要保住。
#[test]
fn chirality_is_dropped_when_the_product_template_stops_specifying_it() {
    let got = flat(
        "[C@;H0;D4;+0:1]-[NH;D2;+0:2]>>Cl-[C;H0;D4;+0:1].[NH2;D1;+0:2]",
        &["CC[C@](Nc1ccc(C)cc1)(C#N)c1cccc(F)c1"],
    );
    assert!(!got.is_empty(), "一个产物都没有");
    for p in &got {
        assert!(
            !p.contains('@'),
            "产物 {p} 还带着手性 —— 产物模板已经不写 `@` 了,继承来的构型该清掉"
        );
    }
}

/// 模板两侧写没写手性,是**四种不同的指令**,不是一个字面值。
///
/// | 反应物侧 | 产物侧 | 含义 |
/// |---|---|---|
/// | 没写 | 没写 | 模板没管 —— 底物的构型带过来 |
/// | 写了 | 没写 | 构型被破坏(另见 [`chirality_is_dropped_when_the_product_template_stops_specifying_it`]) |
/// | 没写 | 写了 | 构型是**新建**的,与底物无关 |
/// | 写了 | 写了 | 相对底物**保留**(两标记相同)或**翻转**(不同) |
///
/// 最后一行最容易做错。把产物侧那个标记照字面写死的话,同一个模板作用在一对
/// 对映体上会给出**同一个**产物 —— 而正确答案是一对对映体。这个错误不改拓扑、
/// 不改原子数、不改价键,只有拿一对对映体分别跑一遍才看得见。
///
/// 所以判据是**产物随不随底物变**,不是某个具体的 SMILES:
///
/// - "新建"那一档:两个对映体底物必须给出同一个产物
/// - "保留"那一档:必须与"两侧都没写"给出同一个产物
/// - "翻转"那一档:必须与"保留"给出互为对映体的产物
#[test]
fn template_chirality_is_an_instruction_not_a_literal() {
    // 一对对映体,差别只在那一个中心上
    const R: &str = "C[C@H](O)CC";
    const S: &str = "C[C@@H](O)CC";

    let one = |rxn: &str, smi: &str| {
        let got = flat(rxn, &[smi]);
        assert_eq!(
            got.len(),
            1,
            "{rxn} 在 {smi} 上该恰好出一个产物,实得 {got:?}"
        );
        got.into_iter().next().unwrap()
    };

    // 两侧都没写:构型跟着底物走。这一条同时守着后面几档的判据不空过 ——
    // 底物的构型若压根传不过来,下面的比较全都成立而毫无意义。
    let inherit = "[C:1]-[OH:2]>>[C:1]-[Cl:2]";
    let (inherit_r, inherit_s) = (one(inherit, R), one(inherit, S));
    assert_ne!(
        inherit_r, inherit_s,
        "模板没管手性,一对对映体底物却给出同一个产物 —— 构型没传过来"
    );

    // 只有产物侧写:构型是新建的,与底物无关
    for rxn in [
        "[C:1]-[OH:2]>>[C@:1]-[Cl:2]",
        "[CH:1]-[OH:2]>>[C@H:1]-[Cl:2]",
    ] {
        assert_eq!(
            one(rxn, R),
            one(rxn, S),
            "{rxn}:产物侧写死了构型,两个对映体底物却给出不同产物"
        );
    }

    // 两侧都写:相对底物保留或翻转。括号氢写不写都一样 —— 它只改参照系,
    // 不改这条约定。
    for (retain, invert) in [
        (
            "[C@:1]-[OH:2]>>[C@:1]-[Cl:2]",
            "[C@:1]-[OH:2]>>[C@@:1]-[Cl:2]",
        ),
        (
            "[C@H:1]-[OH:2]>>[C@H:1]-[Cl:2]",
            "[C@H:1]-[OH:2]>>[C@@H:1]-[Cl:2]",
        ),
    ] {
        assert_eq!(
            one(retain, R),
            inherit_r,
            "{retain}:两侧标记相同该是保留构型,却与两侧都不写给出的不一样"
        );
        assert_eq!(one(retain, S), inherit_s, "{retain}:S 底物那一侧同理");
        assert_eq!(
            one(invert, R),
            inherit_s,
            "{invert}:两侧标记不同该是翻转构型,产物应当是保留那一档的对映体"
        );
        assert_eq!(one(invert, S), inherit_r, "{invert}:S 底物那一侧同理");
    }
}

/// 自由基电子数是**派生量**,不能从底物照抄进产物。
///
/// 底物那个 `[C]` 是三价中性碳,净化会给它记一个自由基电子。模板把它改写成芳香碳,
/// 自由基数若照抄,kekulize 会认为它不能再要双键 —— 整个苯环配不出 Kekulé 结构,
/// 而报错落在环上另一个原子身上,离根因很远。
///
/// 净化里 kekulize 排在自由基重算**之前**(自由基数要等键级定下来才算得出),所以
/// 这个陈旧值一定会被 kekulize 撞上。把同一个产物写出来再读回去却是好的 —— 新解析
/// 的分子该字段本来就是 0。这种"只在反应路径上犯"的毛病最难发现。
#[test]
fn radical_count_is_not_inherited_from_the_substrate() {
    let rxn = "[N+;H0;D3:1]=[C;H0;D3;+0:2]1-[CH;D2;+0:3]=[CH;D2;+0:4]-[C;H0;D3;+0:5]\
               -[CH;D2;+0:6]=[CH;D2;+0:7]-1\
               >>[N;H0;D3;+0:1]-[c;H0;D3;+0:2]1:[cH;D2;+0:3]:[cH;D2;+0:4]\
               :[c;H0;D3;+0:5]:[cH;D2;+0:6]:[cH;D2;+0:7]:1";
    let got = flat(rxn, &["CN(C)[C]1C=CC(C=C1)=[N+](C)C"]);
    assert_eq!(got.len(), 2, "环的对称性给出两组产物,实得 {got:?}");
    for p in &got {
        assert!(
            !p.starts_with("<净化失败"),
            "醌式环芳构化之后应当净化得过,实得 {p}"
        );
    }
    let want = "CN(C)c1ccc(cc1)N(C)C";
    assert_eq!(got, vec![want, want], "产物应当是对苯二胺");
}

/// 两侧都写手性时,比标记之前要先把两张模板的**邻居次序**摆到一起。
///
/// 标记说的是"按这张模板自己的邻居顺序看过去"的构型。产物模板把取代基对调着
/// 写时,同一个 `@` 说的已经是另一种构型 —— 只比标记就会把"翻转"读成"保留",
/// 而产物照样合法、原子数照样对,纯拓扑比对永远发现不了。
///
/// 这里取代基一个没换,变的只有次序,所以守的纯是次序这一维,与"新取代基
/// 顶替谁"那一档无关。
#[test]
fn both_sides_written_are_compared_in_a_common_neighbour_order() {
    // 一对对映体。四个取代基各不相同,中心是货真价实的手性中心。
    const R: &str = "F[C@@H](Cl)Br";
    const S: &str = "F[C@H](Cl)Br";

    let one = |rxn: &str, smi: &str| {
        let got = flat(rxn, &[smi]);
        assert_eq!(
            got.len(),
            1,
            "{rxn} 在 {smi} 上该恰好出一个产物,实得 {got:?}"
        );
        got.into_iter().next().unwrap()
    };

    // 次序不变、标记相同 —— 保留。这一条同时守着下面两条不空过:构型若压根
    // 传不过来,两个对映体会给出同一个产物,后面比什么都成立。
    let keep = "[F:2][C@H:1]([Cl:3])[Br:4]>>[F:2][C@H:1]([Cl:3])[Br:4]";
    assert_eq!(
        one(keep, R),
        "[C@H](F)(Cl)Br",
        "次序标记都没变,构型该原样留着"
    );
    assert_eq!(one(keep, S), "[C@@H](F)(Cl)Br", "S 底物那一侧同理");

    // 产物模板把 F 与 Br 对调着写 —— 置换是奇的。标记仍然相同,可这条模板
    // 说的已经是**翻转**。
    let swapped = "[F:2][C@H:1]([Cl:3])[Br:4]>>[Br:4][C@H:1]([Cl:3])[F:2]";
    assert_eq!(
        one(swapped, R),
        "[C@@H](F)(Cl)Br",
        "{swapped}:两侧标记相同但次序对调了,该翻转 —— 只比标记会读成保留"
    );
    assert_eq!(one(swapped, S), "[C@H](F)(Cl)Br", "S 底物那一侧同理");

    // 次序对调 + 标记不同 —— 两处各反一次,合起来是保留
    let swapped_flipped = "[F:2][C@H:1]([Cl:3])[Br:4]>>[Br:4][C@@H:1]([Cl:3])[F:2]";
    assert_eq!(
        one(swapped_flipped, R),
        "[C@H](F)(Cl)Br",
        "{swapped_flipped}:次序与标记各反一次,合起来该是保留"
    );
    assert_eq!(
        one(swapped_flipped, S),
        "[C@@H](F)(Cl)Br",
        "S 底物那一侧同理"
    );
}

/// 手性重定基要基于**产物当前**的标记,不是反应物的。
///
/// 产物的标记可能被模板定死、也可能被刻意清掉。定基时若回头去读反应物的标记,
/// 就把模板刚做的决定覆盖回去了 —— 清掉的会被恢复,定死的会被换成继承值。
///
/// 这里只守"构型没被丢掉"。构型具体该是哪一个由
/// [`template_chirality_is_an_instruction_not_a_literal`] 守。
#[test]
fn rebasing_uses_the_products_own_tag() {
    // 两侧都写了手性且标记不同 —— 相对底物翻转,不管怎样都得留下构型
    let got = flat("[C@H:1]([OH])>>[C@@H:1]([Cl])", &["C[C@H](O)CC"]);
    assert_eq!(got.len(), 1, "应当恰好一组产物");
    assert!(got[0].contains('@'), "产物 {} 丢了模板定下的构型", got[0]);
}

/// 配位键的**朝向**在两种表示之间要换算。
///
/// 查询侧用 `Dative` / `DativeReversed` 两个基元区分朝向,端点按书写顺序存;
/// 产物侧靠端点顺序表达(`begin` 是给电子的一端)。照搬端点建产物,`<-` 那一支
/// 的箭头就反了。
///
/// 反了不只是写法难看:"接受"的配位键计入受体的价,给电子的那个氧当场变成
/// 三价而被净化拒绝 —— 症状出现在离根因很远的地方。
#[test]
fn dative_arrows_keep_their_direction_through_the_template() {
    for (rxn, why) in [
        ("[Fe:1]>>CO(C)->[Fe:1]", "箭头向右"),
        ("[Fe:1]>>[Fe:1]<-O(C)C", "箭头向左"),
    ] {
        let got = flat(rxn, &["[Fe]"]);
        assert_eq!(got.len(), 1, "{rxn}({why}):应当恰好一个产物");
        assert_eq!(
            got[0],
            canonical("CO(C)->[Fe]"),
            "{rxn}({why}):配位键的朝向变了 —— 两种写法说的是同一件事"
        );
    }
}

/// 环闭合上的配位键同样不能翻。
#[test]
fn dative_arrows_on_ring_closures_keep_their_direction() {
    // 两条模板写的是同一个东西:两个氧都向铁给电子
    let a = flat("[Fe:1]>>[O]1CC[O](C)->[Fe:1]<-1", &["[Fe]"]);
    let b = flat("[Fe:1]>>[O]->1CC[O](C)->[Fe:1]1", &["[Fe]"]);
    assert_eq!(a.len(), 1, "闭环端写箭头:应当恰好一个产物");
    assert_eq!(a, b, "同一个分子的两种写法给出了不同的产物");
    assert!(
        !a[0].contains("<-"),
        "产物 {} 里有指向氧的箭头 —— 模板说的是氧给出电子",
        a[0]
    );
}

/// 质量守恒:产物的重原子总数不该多于底物。
///
/// 这是产物构建最基本的一条,而它**很容易被悄悄破掉** —— 模板之外的部分若按
/// 产物各搬一次,共享的那一段就被复制进每个产物,原子凭空变多。拓扑、价键、
/// 立体全都自洽,只有原子数不对,而没有任何一档判据在看这个数。
#[test]
fn products_never_invent_atoms() {
    // 数字母是不行的:`Cl` 会被数成两个原子。老老实实解析再数。
    let heavy = |s: &str| {
        smiles::parse(s)
            .unwrap_or_else(|e| panic!("{s} 读不回:{}", e.render()))
            .num_atoms()
    };
    for (rxn, sub) in [
        ("[C:1][N:2]>>[C:1].[N:2]", "C1CCNCC1"),
        ("[C:1][N:2]>>[C:1].[N:2]", "CCN"),
        ("[C:1][O:2][C:3]>>[C:1][O:2].[C:3]", "COC"),
        ("[C:1][O:2][C:3]>>[C:1][O:2].[C:3]", "C1CCOCC1"),
        // 取代:氧被删、氯新建,原子数不变
        ("[C:1][OH]>>[C:1]Cl", "OCCO"),
    ] {
        let want = heavy(sub);
        for set in products(rxn, &[sub]) {
            if set.iter().any(|p| p.starts_with('<')) {
                continue; // 净化失败的那组数不了,另有判据管它
            }
            let got: usize = set.iter().map(|p| heavy(p)).sum();
            assert!(
                got <= want,
                "{rxn} 作用在 {sub}({want} 个重原子)上给出 {set:?},\
                 共 {got} 个 —— 产物凭空多出原子"
            );
        }
    }
}

// ---------------------------------------------------------------- 一张图上的匹配

/// 把整个反应物侧当作一张图上的查询来跑,产物规范化成 SMILES 的多重集。
fn on_substrate(rxn_smarts: &str, mols: &[&str]) -> Vec<Vec<String>> {
    let rxn = smarts::parse_reaction(rxn_smarts)
        .unwrap_or_else(|e| panic!("{rxn_smarts}:\n{}", e.render()));
    let inputs: Vec<(MolBuilder, MolProps)> = mols
        .iter()
        .map(|s| {
            let m = sanitized(s);
            let p = MolProps::compute(&m);
            (m, p)
        })
        .collect();
    omgkit_match::run_on_substrate(&rxn, &inputs, 0, false)
        .into_iter()
        .map(|outcome| {
            outcome
                .products
                .into_iter()
                .map(|mut p| match sanitize(&mut p) {
                    Ok(()) => canon::canonical_smiles(&p).smiles,
                    Err(e) => format!("<净化失败: {e}>"),
                })
                .collect()
        })
        .collect()
}

/// 分子内反应:模板的两个片段落在**同一个分子**上。
///
/// 位置式的 `run_reactants` 表达不了这一档 —— 它要求模板片段数等于输入分子数,
/// 这里是 2 对 1,直接交白卷。RDKit 的 `RunReactants` 与 rdchiral 同样跑不了
/// (两者都会抛 "Number of reactants provided does not match")。USPTO-50k 上
/// 这类有 301 条。
///
/// 下面是语料 US08058045B2(第 304 行)的正向模板与底物:分子内 Williamson
/// 成醚,氯与羟基在同一条链上,关成吗啉酮。真值也照抄记录。
#[test]
fn two_template_fragments_can_land_on_one_molecule() {
    let rxn = "Cl-[CH2;D2;+0:1]-[C:2](-[#7:3])=[O;D1;H0:4].[C:5]-[OH;D1;+0:6]\
               >>[#7:3]-[C:2](=[O;D1;H0:4])-[CH2;D2;+0:1]-[O;H0;D2;+0:6]-[C:5]";
    let sub = "O=C(CCl)NC[C@H](O)CO";

    // 位置式:2 个模板片段 vs 1 个分子 —— 契约表达不了
    assert!(
        products(rxn, &[sub]).is_empty(),
        "位置式契约本不该给出产物"
    );

    let sets = on_substrate(rxn, &[sub]);
    assert!(!sets.is_empty(), "一张图上也没找到 —— 分子内没跑通");
    let want = canonical("O=C1CO[C@H](CO)CN1");
    assert!(
        sets.iter().any(|s| s.len() == 1 && s[0] == want),
        "没有给出记录里的产物 {want},实际 {sets:?}"
    );
}

/// 同一条规则在分子间底物上不能改变结果。
///
/// 拼成一张图之后,两个片段落在不同的连通分量上 —— 这与位置式配对是同一件事,
/// 结果必须一致。差别只在调用方不必再把输入的所有排列试一遍。
#[test]
fn one_graph_matching_agrees_with_positional_on_separate_molecules() {
    let rxn = "Cl-[CH2;D2;+0:1]-[C:2](-[#7:3])=[O;D1;H0:4].[C:5]-[OH;D1;+0:6]\
               >>[#7:3]-[C:2](=[O;D1;H0:4])-[CH2;D2;+0:1]-[O;H0;D2;+0:6]-[C:5]";
    let mols = ["O=C(CCl)NC", "CCO"];

    let mut positional: Vec<String> = products(rxn, &mols).into_iter().flatten().collect();
    positional.sort();
    positional.dedup();
    assert!(!positional.is_empty(), "对照组没有产物,判据空过");

    let mut combined: Vec<String> = on_substrate(rxn, &mols).into_iter().flatten().collect();
    combined.sort();
    combined.dedup();
    assert_eq!(combined, positional);
}

/// 两个模板片段不能抢同一个原子。
///
/// 位置式契约靠"分子各不相同"天然保证了这一点;拼成一张图之后必须自己判,
/// 否则同一个羟基会既当亲核试剂又当离去基团那一侧,凭空多出一堆产物。
///
/// 判据:底物只有**一个**可当 `[C:5]-[OH:6]` 的羟基,而它同时也是模板另一个
/// 片段里 `[C:2]` 能碰到的碳;不判不相交的话产物数会多出来。这里比的是
/// "每一组产物里的原子数都不超过底物" —— 抢原子必然导致重复建原子。
#[test]
fn two_template_fragments_never_share_an_atom() {
    let rxn = "Cl-[CH2;D2;+0:1]-[C:2](-[#7:3])=[O;D1;H0:4].[C:5]-[OH;D1;+0:6]\
               >>[#7:3]-[C:2](=[O;D1;H0:4])-[CH2;D2;+0:1]-[O;H0;D2;+0:6]-[C:5]";
    let sub = "O=C(CCl)NC[C@H](O)CO";
    let n_in = heavy(sub);
    for set in on_substrate(rxn, &[sub]) {
        let n_out: usize = set.iter().map(|s| heavy(s)).sum();
        assert!(
            n_out <= n_in,
            "产物 {set:?} 有 {n_out} 个重原子,底物只有 {n_in} —— 原子被抢重了"
        );
    }
}
