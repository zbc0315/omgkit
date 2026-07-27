//! 查询求值:拿一棵查询树去判断一个具体的原子/键是否匹配。
//!
//! # 性质由调用方喂进来
//!
//! 求值需要环成员数、价、隐式氢这些**化学**信息,而它们属于 L2。本模块不去
//! 依赖 L2,而是让调用方填一份 [`AtomProps`] 递进来。
//!
//! 这样分有两个好处:求值逻辑可以脱离化学单独测(手工构造性质,想造什么形状
//! 造什么形状);上层做子结构匹配时本来就要把这些性质**预计算并缓存**起来 ——
//! 匹配过程中同一个原子会被反复求值,每次现算环信息会直接变成性能灾难。
//!
//! # 几个基元的定义,实测得来,不能想当然
//!
//! | 基元 | 含义 | 容易写错成 |
//! |---|---|---|
//! | `D<n>` | 重原子邻居数 | 含氢的总连接数 |
//! | `X<n>` | 总连接数 = 重原子邻居 + 总氢 | 只数重原子 |
//! | `H<n>` | **总**氢数(显式 + 隐式) | 只数隐式 |
//! | `h<n>` | **隐式**氢数 | 总氢 |
//! | `R<n>` | 该原子属于几个环 | 环的大小 |
//! | `r<n>` | **最小**环的大小恰为 n | "在某个 n 元环里" |
//! | `x<n>` | 该原子上有几条环键 | 环的个数 |
//!
//! `r` 那条最隐蔽。以氢化茚(5 元环与 6 元环稠合,共 9 个环原子)为例:
//! `[r5]` 命中 5 个、`[r6]` 命中 4 个,加起来正好 9 —— 每个原子只算一次,
//! 归到它**最小**的那个环。若按"在某个 n 元环里"理解,两个稠合原子会被
//! 重复计入,`[r5]` 就该命中 7 个。

use omgkit_core::{BondDirection, BondOrder, ChiralTag};

use super::expr::{AtomExpr, AtomPrim, BondExpr, BondPrim};
use super::QueryMol;

/// 求值一个原子查询所需的全部性质。
///
/// 字段名对应 SMARTS 基元,含义见模块文档的表。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct AtomProps {
    /// 原子序数。0 表示通配原子。
    pub atomic_num: u8,
    /// 芳香标志(净化后的感知结果)
    pub aromatic: bool,
    /// 形式电荷
    pub charge: i32,
    /// 同位素。0 表示未指定。
    pub isotope: u16,
    /// `D` —— 重原子邻居数
    pub degree: u32,
    /// `H` —— 总氢数(显式 + 隐式)
    pub total_hs: u32,
    /// `h` —— 隐式氢数
    pub implicit_hs: u32,
    /// `v` —— 总价
    pub valence: u32,
    /// `R` —— 该原子属于几个环
    pub ring_count: u32,
    /// `r` —— 最小环的大小;0 表示不在任何环中
    pub min_ring_size: u32,
    /// `x` —— 该原子上的环键数
    pub ring_bonds: u32,
    /// 立体标记
    pub chiral_tag: ChiralTag,
    /// 反应原子映射号。0 表示无。
    pub atom_map: u16,
}

impl AtomProps {
    /// `X` —— 总连接数 = 重原子邻居 + 总氢
    #[must_use]
    pub fn total_degree(&self) -> u32 {
        self.degree + self.total_hs
    }
}

/// 求值一个键查询所需的性质。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BondProps {
    /// 键级
    pub order: BondOrder,
    /// 是否在环中
    pub in_ring: bool,
    /// 方向键 `/` `\`
    pub direction: BondDirection,
    /// 配位键的给体是否是**查询模式里写在前面**的那一端。
    ///
    /// 配位键的方向有语义,而 `->` 与 `<-` 是两个不同的基元,所以求值时
    /// 必须知道"当前这条键在模式里的朝向"。非配位键时该字段无意义。
    pub dative_forward: bool,
}

impl Default for BondProps {
    fn default() -> Self {
        Self {
            order: BondOrder::Unspecified,
            in_ring: false,
            direction: BondDirection::None,
            dative_forward: true,
        }
    }
}

/// 判断原子是否满足查询。
///
/// `recursive` 用来求值 `$(...)`:给它一个子模式,它回答"当前原子能否作为
/// 该模式的首原子匹配上"。子结构匹配器要到 L5 才有,在那之前可以传一个
/// 恒假的闭包 —— 但那会让含递归的模式**静默不匹配**,所以调用方必须自己
/// 清楚这一点。
pub fn atom_matches(
    expr: &AtomExpr,
    props: &AtomProps,
    recursive: &mut dyn FnMut(&QueryMol) -> bool,
) -> bool {
    match expr {
        AtomExpr::Prim(p) => prim_matches(p, props, recursive),
        AtomExpr::Not(e) => !atom_matches(e, props, recursive),
        AtomExpr::And(parts) => parts.iter().all(|e| atom_matches(e, props, recursive)),
        AtomExpr::Or(parts) => parts.iter().any(|e| atom_matches(e, props, recursive)),
    }
}

fn prim_matches(p: &AtomPrim, a: &AtomProps, recursive: &mut dyn FnMut(&QueryMol) -> bool) -> bool {
    match p {
        AtomPrim::Any => true,
        AtomPrim::Aromatic => a.aromatic,
        AtomPrim::Aliphatic => !a.aromatic,
        AtomPrim::Element { z, aromatic } => {
            a.atomic_num == *z && aromatic.map_or(true, |want| want == a.aromatic)
        }
        AtomPrim::Degree(n) => a.degree == *n,
        AtomPrim::TotalDegree(n) => a.total_degree() == *n,
        AtomPrim::TotalHs(n) => a.total_hs == *n,
        AtomPrim::ImplicitHs(n) => a.implicit_hs == *n,
        AtomPrim::Valence(n) => a.valence == *n,
        // 裸写的 `R` / `x` 是"任意非零",带数字才比具体值
        AtomPrim::RingCount(n) => n.map_or(a.ring_count > 0, |k| a.ring_count == k),
        AtomPrim::RingBondCount(n) => n.map_or(a.ring_bonds > 0, |k| a.ring_bonds == k),
        AtomPrim::RingSize(n) => n.map_or(a.min_ring_size > 0, |k| a.min_ring_size == k),
        AtomPrim::Charge(c) => a.charge == *c,
        AtomPrim::Isotope(i) => a.isotope == *i,
        // 映射号是**标签不是条件**。`[C:99]` 匹配任何碳,与目标原子自己的
        // 映射号无关 —— 它的用途是在反应里把反应物原子与产物原子对应起来。
        //
        // 当成条件的话,所有反应模板都会一个也匹配不上(底物的映射号通常是 0),
        // 而错误表现为"反应没有产物",完全指不到根因。
        AtomPrim::AtomMap(_) => true,
        // 手性**不是逐原子能判的性质**,它是**映射**的性质:标记相对各自分子的
        // 邻居存储顺序,要知道查询的邻居分别映到底物的哪些原子才能算宇称。
        // 这里拿不到映射,所以一律放行,由匹配器在映射完成后校验。
        //
        // 直接比原始标记(`a.chiral_tag == *t`)是拿两个参照系里的值去比,
        // 得到的构型可以正好相反。
        AtomPrim::Chirality(_) => true,
        AtomPrim::Recursive(sub) => recursive(sub),
    }
}

/// 从原子表达式里推出**允许的元素集合**。返回 `None` 表示推不出约束。
///
/// 这是给匹配器做候选估算用的:一条 `CCCCCCCCBr` 的模式,从溴那一端起头
/// 候选只有一个,从碳起头候选是几百个。要挑对起点就得先知道"这个查询原子
/// 大概能匹配多少个目标原子",而元素是最强也最便宜的那一维。
///
/// # 只能放宽,不能收紧
///
/// 推不出来时必须返回 `None`(视作"什么都可能"),宁可估得保守。估紧了会让
/// 匹配器**漏掉**候选 —— 那不是慢,是错。所以:
///
/// - `And` 取交集:两个条件都要满足
/// - `Or` 取并集;任一支推不出约束,整体就推不出
/// - `Not` 一律返回 `None`:补集要枚举整张周期表,而且对
///   `[!C]` 这种"非脂肪碳"来说,芳香碳仍然是允许的,补集算错就会漏候选
#[must_use]
pub fn allowed_elements(expr: &AtomExpr) -> Option<std::collections::BTreeSet<u8>> {
    match expr {
        AtomExpr::Prim(AtomPrim::Element { z, .. }) => Some([*z].into_iter().collect()),
        AtomExpr::Prim(_) => None,
        AtomExpr::Not(_) => None,
        AtomExpr::And(parts) => {
            let mut acc: Option<std::collections::BTreeSet<u8>> = None;
            for p in parts {
                if let Some(s) = allowed_elements(p) {
                    acc = Some(match acc {
                        None => s,
                        Some(a) => a.intersection(&s).copied().collect(),
                    });
                }
            }
            acc
        }
        AtomExpr::Or(parts) => {
            let mut acc = std::collections::BTreeSet::new();
            for p in parts {
                // 任一支没有元素约束,整个析取式就没有
                acc.extend(allowed_elements(p)?);
            }
            Some(acc)
        }
    }
}

/// 判断键是否满足查询。
#[must_use]
pub fn bond_matches(expr: &BondExpr, props: &BondProps) -> bool {
    match expr {
        BondExpr::Prim(p) => bond_prim_matches(*p, props),
        BondExpr::Not(e) => !bond_matches(e, props),
        BondExpr::And(parts) => parts.iter().all(|e| bond_matches(e, props)),
        BondExpr::Or(parts) => parts.iter().any(|e| bond_matches(e, props)),
    }
}

fn bond_prim_matches(p: BondPrim, b: &BondProps) -> bool {
    match p {
        BondPrim::Any => true,
        BondPrim::Single => b.order == BondOrder::Single,
        BondPrim::Double => b.order == BondOrder::Double,
        BondPrim::Triple => b.order == BondOrder::Triple,
        BondPrim::Quadruple => b.order == BondOrder::Quadruple,
        BondPrim::Aromatic => b.order == BondOrder::Aromatic,
        BondPrim::InRing => b.in_ring,
        // 方向**不是逐键能判的性质**,与手性同源:`/` 相对键自己的
        // begin → end 朝向,查询与底物的朝向不同,直接比就成了"比写法"。
        //
        // `F/C=C/F` 与 `C(\F)=C/F` 是同一个分子,直接比会一个匹配、一个不
        // 匹配。而且单独一条 `/` 本身什么都不表示 —— 顺反是双键加两条参照键
        // 的性质,同样要等映射齐了才判。
        //
        // 但**键级约束要留着**:整个返回 true 会让 `/` 匹配上双键,于是
        // `[C]/[C]=[C][C]` 比 `[C][C]=[C][C]` 匹配得**更多** —— 加了方向
        // 只该收紧不该放宽。
        //
        // 而这里的约束就是[默认键]那一个:**单键或芳香键**。写成"仅单键"会比
        // 默认键**更严**,于是在芳香键上写方向的模板一处也匹配不上 —— 稠环模板
        // 里这是常态,方向落在环内某根芳香键上并不罕见。方向该改的是朝向,不是
        // 键级;能不能带方向由写模板的人负责,不由匹配代为否决。
        //
        // SMILES 解析那边早就是这么做的(见 `smiles.rs` 的 `bond_symbol`:
        // "`/` `\` 是**纯方向**标记,不指定键级")—— OpenSMILES 把方向键描述为
        // "单键",照字面实现会错。同一条结论要在 SMARTS 这边一并成立,否则
        // 同一根键读进来是芳香的,查询却按单键去配。
        //
        // [默认键]: BondExpr::default_bond
        BondPrim::UpRight | BondPrim::DownRight => {
            matches!(b.order, BondOrder::Single | BondOrder::Aromatic)
        }
        // 配位键的两个基元靠朝向区分
        BondPrim::Dative => b.order == BondOrder::Dative && b.dative_forward,
        BondPrim::DativeReversed => b.order == BondOrder::Dative && !b.dative_forward,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smarts;

    /// 不含递归的模式用它求值 —— 遇到 `$(...)` 直接 panic,免得静默不匹配。
    fn no_recursion(_: &QueryMol) -> bool {
        panic!("这条模式含递归 SMARTS,但测试没有提供求值器");
    }

    fn matches(pat: &str, props: &AtomProps) -> bool {
        let q = smarts::parse(pat).unwrap_or_else(|e| panic!("{pat}: {}", e.render()));
        atom_matches(&q.atoms[0], props, &mut no_recursion)
    }

    /// 萘的 CH 原子(非稠合)
    fn naphthalene_ch() -> AtomProps {
        AtomProps {
            atomic_num: 6,
            aromatic: true,
            degree: 2,
            total_hs: 1,
            implicit_hs: 1,
            valence: 4,
            ring_count: 1,
            min_ring_size: 6,
            ring_bonds: 2,
            ..AtomProps::default()
        }
    }

    /// 萘的稠合碳
    fn naphthalene_fusion() -> AtomProps {
        AtomProps {
            atomic_num: 6,
            aromatic: true,
            degree: 3,
            total_hs: 0,
            implicit_hs: 0,
            valence: 4,
            ring_count: 2,
            min_ring_size: 6,
            ring_bonds: 3,
            ..AtomProps::default()
        }
    }

    /// 各基元的定义。数字取自对萘的实测。
    #[test]
    fn primitive_semantics() {
        let ch = naphthalene_ch();
        let fu = naphthalene_fusion();

        assert!(matches("[c]", &ch) && matches("[c]", &fu));
        assert!(!matches("[C]", &ch), "大写 C 要求脂肪碳");
        assert!(matches("[a]", &ch) && !matches("[A]", &ch));
        assert!(matches("[#6]", &ch), "#6 不限芳香性");

        // R:属于几个环
        assert!(matches("[R1]", &ch) && !matches("[R1]", &fu));
        assert!(matches("[R2]", &fu) && !matches("[R2]", &ch));
        assert!(
            matches("[R]", &ch) && matches("[R]", &fu),
            "裸 R = 在任意环中"
        );

        // r:最小环的大小
        assert!(matches("[r6]", &ch) && matches("[r6]", &fu));
        assert!(!matches("[r5]", &ch));

        // x:环键数
        assert!(matches("[x2]", &ch) && matches("[x3]", &fu));

        // D / X / H / h / v
        assert!(matches("[D2]", &ch) && matches("[D3]", &fu));
        assert!(matches("[X3]", &ch) && matches("[X3]", &fu), "度 + 总氢");
        assert!(matches("[H1]", &ch) && matches("[H0]", &fu));
        assert!(matches("[h1]", &ch) && matches("[h0]", &fu));
        assert!(matches("[v4]", &ch));
    }

    /// 逻辑运算的求值,以及优先级在求值上确实产生了不同结果。
    #[test]
    fn logic_evaluation() {
        let ch = naphthalene_ch();
        assert!(matches("[c,n]", &ch));
        assert!(!matches("[n,o]", &ch));
        assert!(matches("[c;R1]", &ch));
        assert!(!matches("[c;R2]", &ch));
        assert!(matches("[!n]", &ch));
        assert!(matches("[!C;!N]", &ch), "既非脂肪碳也非脂肪氮");

        // `[c,n;H0]` = (c 或 n) 且 H0 —— 萘的 CH 有 1 个氢,不匹配
        assert!(!matches("[c,n;H0]", &ch));
        // `[c,n&H0]` = c 或 (n 且 H0) —— 它是 c,匹配
        assert!(matches("[c,n&H0]", &ch), "优先级不同,结果就不同");
    }

    /// 裸写与带数字的取值区别。
    #[test]
    fn bare_versus_numbered() {
        let acyclic = AtomProps {
            atomic_num: 6,
            degree: 1,
            total_hs: 3,
            implicit_hs: 3,
            valence: 4,
            ..AtomProps::default()
        };
        assert!(!matches("[R]", &acyclic), "不在环中");
        assert!(!matches("[r]", &acyclic));
        assert!(!matches("[x]", &acyclic));
        assert!(matches("[R0]", &acyclic), "R0 = 不属于任何环");
        assert!(matches("[D1]", &acyclic) && matches("[X4]", &acyclic));
    }

    /// 键的求值,包括配位键靠朝向区分的两个基元。
    #[test]
    fn bond_evaluation() {
        let single = BondProps {
            order: BondOrder::Single,
            ..BondProps::default()
        };
        let ring_double = BondProps {
            order: BondOrder::Double,
            in_ring: true,
            ..BondProps::default()
        };
        let q = |pat: &str| {
            let m = smarts::parse(pat).unwrap_or_else(|e| panic!("{pat}: {}", e.render()));
            m.bonds[0].clone()
        };

        assert!(bond_matches(&q("C-C"), &single));
        assert!(!bond_matches(&q("C=C"), &single));
        assert!(bond_matches(&q("C~C"), &single), "~ 匹配一切");
        assert!(bond_matches(&q("C=@C"), &ring_double), "环内双键");
        assert!(!bond_matches(&q("C=!@C"), &ring_double));
        assert!(
            bond_matches(
                &q("C=!@C"),
                &BondProps {
                    order: BondOrder::Double,
                    in_ring: false,
                    ..BondProps::default()
                }
            ),
            "非环双键"
        );

        // 省略键符号 = 单键或芳香键
        assert!(bond_matches(&q("CC"), &single));
        assert!(bond_matches(
            &q("CC"),
            &BondProps {
                order: BondOrder::Aromatic,
                ..BondProps::default()
            }
        ));
        assert!(!bond_matches(&q("CC"), &ring_double));

        // 配位键:同一条键,朝向不同,匹配的基元也不同
        let fwd = BondProps {
            order: BondOrder::Dative,
            dative_forward: true,
            ..BondProps::default()
        };
        let rev = BondProps {
            dative_forward: false,
            ..fwd
        };
        assert!(bond_matches(&q("N->[Cu]"), &fwd));
        assert!(!bond_matches(&q("N->[Cu]"), &rev));
        assert!(bond_matches(&q("[Cu]<-N"), &rev));
        assert!(!bond_matches(&q("[Cu]<-N"), &fwd));
    }

    /// 映射号不参与匹配。
    #[test]
    fn atom_map_is_a_label_not_a_condition() {
        let plain = AtomProps {
            atomic_num: 6,
            ..AtomProps::default()
        };
        let mapped = AtomProps {
            atom_map: 7,
            ..plain
        };
        // 目标没有映射号,模板写了任意映射号,都该命中
        assert!(matches("[C:1]", &plain));
        assert!(matches("[C:99]", &plain));
        // 目标带映射号,模板写不同的号,也该命中
        assert!(matches("[C:1]", &mapped));
        assert!(matches("[C:0]", &mapped));
    }

    /// 元素约束的推导。**只能放宽不能收紧** —— 估紧了会漏候选。
    #[test]
    fn element_constraint_inference() {
        let el = |s: &str| {
            let q = smarts::parse(s).unwrap_or_else(|e| panic!("{s}: {}", e.render()));
            allowed_elements(&q.atoms[0])
        };
        assert_eq!(el("[C]"), Some([6].into_iter().collect()));
        assert_eq!(el("[c]"), Some([6].into_iter().collect()), "芳香碳也是碳");
        assert_eq!(el("[#7]"), Some([7].into_iter().collect()));
        assert_eq!(el("[C,N]"), Some([6, 7].into_iter().collect()), "或取并");
        assert_eq!(el("[C;H3]"), Some([6].into_iter().collect()), "与取交");
        assert_eq!(
            el("[C,N;H3]"),
            Some([6, 7].into_iter().collect()),
            "与的一支没有元素约束时,取另一支的"
        );
        assert_eq!(el("[#6;#7]"), Some([].into_iter().collect()), "矛盾即空集");

        // 推不出来的一律 None
        assert_eq!(el("[*]"), None);
        assert_eq!(el("[R1]"), None);
        assert_eq!(el("[!C]"), None, "补集不推 —— 芳香碳仍可能被 [!C] 接受");
        assert_eq!(el("[C,R1]"), None, "或的一支没约束,整体就没有");
        assert_eq!(el("[a]"), None);
    }

    /// 递归 SMARTS 走的是调用方给的闭包。
    #[test]
    fn recursive_goes_through_the_callback() {
        let q = smarts::parse("[$(CC)]").unwrap();
        let props = AtomProps::default();

        let mut called = 0;
        {
            // 闭包借着 `called`,读它之前要先让借用结束
            let mut yes = |sub: &QueryMol| {
                called += 1;
                assert_eq!(sub.num_atoms(), 2, "子模式应当是 CC");
                true
            };
            assert!(atom_matches(&q.atoms[0], &props, &mut yes));
        }
        assert_eq!(called, 1, "闭包应当被调用一次");

        let mut no = |_: &QueryMol| false;
        assert!(!atom_matches(&q.atoms[0], &props, &mut no));
    }
}
