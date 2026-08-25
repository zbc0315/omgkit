//! 把反应丢弃的片段收口成分子。
//!
//! 反应记录普遍只写主产物:酸与醇成酯,记录里只有酯,水没了。模板是从记录抽的,
//! 于是模板也只描述主产物,[`run_reactants`](crate::run_reactants) 照模板办事,
//! 那部分原子就**从输出里消失**。消失不报错,而它恰恰是引擎自己在破坏质量守恒。
//!
//! 本模块把那批原子接回来:[`Outcome::discarded`](crate::Outcome::discarded) 是
//! 事实记录(哪些原子没进产物),本模块负责**推断**它们收口成了什么分子。
//!
//! # 收口靠的是原子账,不是化学直觉
//!
//! 模板里没有"离去基团变成了什么"这条信息 —— 它只说哪些原子不要了。但**账是
//! 可以算的**:
//!
//! | 量 | 含义 |
//! |---|---|
//! | `open_valence` | 片段与保留部分之间被切断的键的**键级和**,每一处都是一个未闭合的价 |
//! | `fragment_hydrogens` | 片段自身已经带着的氢 |
//! | `delta_h` | 氢预算 = 底物总氢 − 产物总氢 |
//! | `need` | `delta_h − fragment_hydrogens`,要往片段上补几个氢 |
//! | `delta_charge` | 电荷预算 = 底物总电荷 − 产物总电荷 |
//! | `charge_shift` | `delta_charge − 片段自带电荷`,片段要**净得到**多少电荷 |
//! | `remaining` | `open_valence + charge_shift − need`,收完之后还剩几处空价 |
//!
//! 补一个氢填掉一处空价;`need` 为负则要**摘**氢,每摘一个反而多出一处空价 ——
//! 所以两种情形下是同一个式子。
//!
//! # 填空价的不只有氢,还有**形式电荷**
//!
//! 溴从 C–Br 上断下来时带走那对电子,成的是 Br⁻ —— 一处空价被一个负电荷填掉,
//! **一个氢都不需要**。只按氢记账的话这一档永远配不成对,会被误判成"记录不平"。
//! 实测这是最大的一处误判:季铵化、亲核取代这些反应里离去基团本来就以阴离子离去。
//!
//! 反过来,正电荷**多出**一处价(铵氮有四根键),所以电荷是带符号进账的。
//!
//! 剩余空价必须两两成键消化掉:
//!
//! - `remaining == 0` → 只补氢就闭合([`Verdict::Capped`])
//! - `remaining` 为正偶数 → 还要成 `remaining/2` 根键([`Verdict::Bonded`])
//! - `remaining` 为奇数或负数 → **闭不上**([`Verdict::Unresolved`]),不猜
//!
//! 最后一档不是实现能力不足,是**这条记录本身给不出答案**:多半是记录漏写了
//! 贡献原子的试剂(水、酸、碱、氧化剂、还原剂),缺的原子不在手里,收口就无从
//! 谈起。引擎在这里应当明说"答不了",而不是编一个分子出来 —— 编出来的分子
//! 拓扑合法、能净化、看不出破绽,只是错的。
//!
//! # 形式副产物,不是分离得到的副产物
//!
//! 本模块给的是**账平且价键填满**的那个分子。多数时候它就是实际副产物(水、
//! 氯化氢、醇),但对会自发分解的就不是:Boc 脱保护的形式副产物是叔丁基碳酸,
//! 实际拿到的是二氧化碳加异丁烯。分解要靠一张规则表,那是另一件事,不在这里 ——
//! 形式副产物有硬判据(原子账与电荷账精确闭合,不依赖任何记录),分解规则没有,
//! 混在一个输出里就再也分不清哪个是证出来的、哪个是猜的。
//!
//! ```no_run
//! # use omgkit_core::MolBuilder;
//! # use omgkit_match::{byproduct, run_reactants, MolProps};
//! # fn demo(rxn: &omgkit_io::smarts::Reaction, acid: MolBuilder, alcohol: MolBuilder) {
//! let inputs = [acid, alcohol];
//! let props: Vec<_> = inputs.iter().map(MolProps::compute).collect();
//! let pairs: Vec<_> = inputs.iter().cloned().zip(props).collect();
//! for outcome in run_reactants(rxn, &pairs, 0, false) {
//!     let by = byproduct::reconstruct(&inputs, &outcome);
//!     // by.molecules 是收口出来的副产物;by.verdict 说明它是哪一档给的
//! }
//! # }
//! ```

use omgkit_core::{AtomFlags, BondData, BondOrder, MolBuilder};

use crate::react::{align_for_rebase, components, permutation_is_odd, Outcome};

/// 收口最多肯成几根键。
///
/// 语料实测:要成键的那一档里 90.6% 只需一根、9.4% 两根、0.1% 三根,再多的一条
/// 没有。上限摆在这里是为了让"搜索爆掉"变成一条**说得出口的结论**
/// ([`Unresolved::TooManyBonds`]),而不是一个跑很久之后给出的可疑答案。
const MAX_BONDS: u32 = 4;

/// 收口用的原子账。每一项都可以由调用方自己重算,用来核对结论。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Budget {
    /// 被切断的键的键级和
    pub open_valence: u32,
    /// 片段自身带着的氢
    pub fragment_hydrogens: u32,
    /// 氢预算:底物总氢 − 产物总氢
    pub delta_h: i32,
    /// 要补的氢数(负数表示要摘)
    pub need: i32,
    /// 收完氢与电荷之后剩下几处空价
    pub remaining: i32,
    /// 电荷预算:底物总电荷 − 产物总电荷
    pub delta_charge: i32,
    /// 片段自身带着的形式电荷
    pub fragment_charge: i32,
    /// 片段要**净得到**的电荷 = `delta_charge − fragment_charge`。
    ///
    /// 负数表示片段要拿到负电荷(离去基团以阴离子离去),每一个填掉一处空价;
    /// 正数表示要拿到正电荷,每一个反而**多出**一处空价。
    pub charge_shift: i32,
}

/// 闭不上账的原因。**每一档都意味着"不猜"**,不是"暂时没做"。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Unresolved {
    /// 剩余空价是奇数 —— 配不成对。多半是记录漏了贡献原子的试剂。
    OddValence,
    /// 氢预算比空价还多 —— 产物侧的氢比底物能给的还少,记录本身不平。
    BudgetExceedsValence,
    /// 氢预算是**负数** —— 产物的氢比全部贡献反应物加起来还多。
    ///
    /// `delta_h` 就是副产物应有的氢数,负数在物理上讲不通,只可能是记录里少写了
    /// 供氢的试剂。硝基还原成胺是最典型的一档:还原剂不在记录里,那两个氢凭空
    /// 出现在产物上。
    ///
    /// 单独成档是因为它**指向的东西很具体**(缺供氢试剂),混进
    /// [`OddValence`](Self::OddValence) 里就只剩"配不成对"这一句没信息的话。
    HydrogenBudgetNegative,
    /// 要成的键超过本实现肯找的上限(4 根)。
    ///
    /// 上限是实测定的:语料里要成键的那一档,90.6% 只需一根、9.4% 两根、
    /// 0.1% 三根,再多的一条没有。摆一个上限是为了让"搜索爆掉"变成一条
    /// **说得出口的结论**,而不是一个跑很久之后给出的可疑答案。
    TooManyBonds,
    /// 产物净化不过,氢预算无从算起。
    ///
    /// 产物的隐式氢是**净化才填**的派生量,不净化就没有可比的总氢数。
    ProductsUnsanitizable,
    /// 剩余空价配不成对:只剩一个位点,或只剩两个卤素(它们该拿氢,不该互相成键)。
    ///
    /// 与 [`TooManyBonds`](Self::TooManyBonds) 分开是因为**指向的东西不同**:
    /// 那一档是"要成的键太多、本实现不找了",这一档是"这些空价物理上配不起来"。
    /// 混成一个出口的话,归因脚本会把后者全算成"搜索爆掉",而它其实和记录漏试剂
    /// 是同一类问题。
    NoPairing,
    /// 收口出来的东西净化不过 —— 说明这条收口路线不成立。
    FragmentUnsanitizable,
    /// 收口成了一个**几何上不可能**的结构:小环里的三键。
    ///
    /// 三键要求两端与各自的取代基共线,塞进小环里排不下 —— 苯炔那类只在
    /// 瞬态中存在,不会作为副产物被分离到。
    ///
    /// **这一档必须单独拦**,因为它躲得过前面每一道:原子账平、电荷账平、
    /// 净化也过得去。价规则管的是"几根键",管不到"这几根键摆得下摆不下"。
    /// 拦不住的话,输出的是一个配平、合法、下游任何检查都看不出问题的**错分子**。
    StrainedClosure,
    /// 底物 kekulize 不了,断口的键级因而定不下来。
    ///
    /// 本模块要求底物是**净化过**的,而净化里就跑过一次 kekulize,所以正常
    /// 调用路径上到不了这一档。它守的是"调用方递进来一个没净化的分子" ——
    /// 那时若默默按芳香键算,给出的账会稳定地错而毫无迹象。
    SubstrateUnkekulizable,
    /// 收口之后账没对上。**这是本模块自己的缺陷**,不是数据的问题,
    /// 出现即应当当作 bug 追。
    BudgetMismatch,
}

/// 收口的结论。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Verdict {
    /// 一个原子都没被丢弃 —— 本来就守恒,没有副产物
    Nothing,
    /// 只补氢就闭合了
    Capped,
    /// 补氢之后还成了 `bonds` 根键
    Bonded {
        /// 成了几根键
        bonds: u32,
    },
    /// 闭不上,不给分子
    Unresolved(Unresolved),
}

impl Verdict {
    /// 账有没有精确闭合。闭合才有 [`Byproducts::molecules`]。
    #[must_use]
    pub fn is_closed(self) -> bool {
        matches!(self, Verdict::Capped | Verdict::Bonded { .. })
    }
}

/// 收口的产出。
#[derive(Debug, Clone)]
pub struct Byproducts {
    /// 收口出来的分子,按连通分量切开。**只在 [`Verdict::is_closed`] 为真时非空。**
    pub molecules: Vec<MolBuilder>,
    /// 这个结论是哪一档给的
    pub verdict: Verdict,
    /// 算账的中间量,供调用方复核
    pub budget: Budget,
}

/// 片段里的一个原子还欠几处价。
///
/// # `sites[i]` 说的就是片段里的第 i 个原子
///
/// 片段的每个原子都建一条,而且建的顺序与原子入图的顺序完全一致,所以下标本身
/// 就是对应关系,**不再单独存一个 `idx`**。
///
/// 这不只是省一个字段:存了 `idx` 就得靠线性查找把它找回来,而那三处查找都在
/// "按原子"或"按键"的循环里 —— 正是本仓库反复警告的那个形状(在按原子的循环里
/// 做一件正比于整个分子的事)。去掉字段之后查找变成下标索引,顺带也没有了
/// "两者对不上"这种可能。
struct Site {
    /// 还没闭合的价
    opens: u32,
    /// 因为"摘氢"而额外多出来的空价,写回时要从显式氢里扣掉
    borrowed_h: u32,
}

/// 收口**做了多少事**。整数、确定,debug 与 release 逐位相同。
///
/// # 为什么要数,而不是量耗时
///
/// 收口的一批操作天然是"按位点"的,而位点表是**每个片段原子一条** ——
/// 于是"遍历位点"这件事的代价正比于整个片段。在按键或按位点的循环里再做
/// 一次这样的遍历,就是本仓库反复警告的那个形状。
///
/// 判它的判据先前是**每原子耗时**,而墙钟会抖:同一个文件里已经记着两次
/// (一次抽风打红了改别的 crate 的提交,一次把"离散度"错判成"增长")。
/// 匹配那两条早已换成数 `candidate_tests`,这里补上最后一条。
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct CloseStats {
    /// 位点表被访问的次数合计 —— 定电荷、摘氢、配对、写回四处循环。
    ///
    /// 记在各条早退判断**之前**:记在之后的话,被跳过的那些数不到,
    /// 而"跳过得够不够早"正是要守的东西。
    pub site_visits: u64,
    /// 整个片段被走一遍的次数(连通分量)。
    ///
    /// 单独数是因为它藏在 `form_bonds` 的 `while formed < to_bond` **里面** ——
    /// 每成一根键就重算一次连通分量。眼下 `to_bond` 只有两三,看不出来;
    /// 它一旦随片段规模走,这就是一个平方项,而位点计数看不见。
    pub fragment_scans: u64,
}

/// 与 [`reconstruct`] 同一条路,外加工作量计数。
#[must_use]
pub fn reconstruct_counted(
    reactants: &[MolBuilder],
    outcome: &Outcome,
) -> (Byproducts, CloseStats) {
    let mut stats = CloseStats::default();
    let by = reconstruct_inner(reactants, outcome, &mut stats);
    (by, stats)
}

/// 把 [`Outcome`] 丢弃的原子收口成分子。
///
/// `reactants` 必须是**传给引擎的那一组原始分子**(已净化),顺序一致 ——
/// [`Outcome::discarded`](crate::Outcome::discarded) 里的下标是相对它们说的。
///
/// # 为什么要在这里净化产物
///
/// [`run_reactants`](crate::run_reactants) 有意不净化产物,而氢预算要拿产物的
/// 总氢来算,隐式氢又是净化才填的。所以本函数在**副本**上净化一次:被测对象
/// 的契约一个字没改,而算账拿到的是有意义的数。产物净化不过时如实报
/// [`Unresolved::ProductsUnsanitizable`],不去猜。
#[must_use]
pub fn reconstruct(reactants: &[MolBuilder], outcome: &Outcome) -> Byproducts {
    reconstruct_counted(reactants, outcome).0
}

fn reconstruct_inner(
    reactants: &[MolBuilder],
    outcome: &Outcome,
    stats: &mut CloseStats,
) -> Byproducts {
    let empty = Budget {
        open_valence: 0,
        fragment_hydrogens: 0,
        delta_h: 0,
        need: 0,
        remaining: 0,
        delta_charge: 0,
        fragment_charge: 0,
        charge_shift: 0,
    };
    if outcome.discarded.iter().all(Vec::is_empty) {
        return Byproducts {
            molecules: Vec::new(),
            verdict: Verdict::Nothing,
            budget: empty,
        };
    }

    // 产物侧的总氢与总电荷 —— 在副本上净化,不动调用方手上的东西
    let mut h_products: i32 = 0;
    let mut q_products: i32 = 0;
    for p in &outcome.products {
        let mut copy = p.clone();
        if omgkit_chem::sanitize(&mut copy).is_err() {
            return Byproducts {
                molecules: Vec::new(),
                verdict: Verdict::Unresolved(Unresolved::ProductsUnsanitizable),
                budget: empty,
            };
        }
        h_products += total_hydrogens(&copy);
        q_products += total_charge(&copy);
    }

    // 片段从**凯库勒式**的副本上切,而不是从芳香式的原分子上切。
    //
    // 两件事都取决于这一步:
    //
    // - **断口的空价数**。断口的键级是底物的性质,苯环上一个碳的两根环键按凯库
    //   勒式是一单一双、合计 3;按"芳香键算 1"数成 2,奇偶就翻了,整条被误判成
    //   收不平。
    // - **片段内部的键**。芳香键与芳香标志一旦搬进不再成环的片段,净化必然报
    //   "原子不在环中却带着芳香标志"——那时报出来的理由(收口路线不成立)是
    //   错的,真正的原因是这里把标志搬错了地方。
    //
    // `kekulize` 会连标志位一起清掉,所以一次调用两件事都办了。
    let mut kekulized: Vec<MolBuilder> = Vec::with_capacity(reactants.len());
    for m in reactants {
        let mut k = m.clone();
        if omgkit_chem::kekulize(&mut k).is_err() {
            return Byproducts {
                molecules: Vec::new(),
                verdict: Verdict::Unresolved(Unresolved::SubstrateUnkekulizable),
                budget: empty,
            };
        }
        kekulized.push(k);
    }
    let (mut frag, mut sites) = build_fragment(&kekulized, &outcome.discarded);
    let open_valence: u32 = sites.iter().map(|s| s.opens).sum();
    let fragment_hydrogens = u32::try_from(total_hydrogens(&frag)).unwrap_or(0);
    let fragment_charge = total_charge(&frag);

    let h_inputs: i32 = reactants.iter().map(total_hydrogens).sum();
    let q_inputs: i32 = reactants.iter().map(total_charge).sum();
    let delta_h = h_inputs - h_products;
    let delta_charge = q_inputs - q_products;
    let need = delta_h - i32::try_from(fragment_hydrogens).unwrap_or(0);
    let charge_shift = delta_charge - fragment_charge;

    // 电荷要**先落定**再算剩余空价。落在哪个元素上决定了它是填掉一处价还是多出
    // 一处(见 `apply_charges`),所以"落定前先估一个 remaining"是估不准的。
    let charges_ok = charge_shift == 0 || apply_charges(&mut frag, &mut sites, charge_shift, stats);
    let opens_after: i32 = sites
        .iter()
        .map(|s| i32::from(u16::try_from(s.opens).unwrap_or(0)))
        .sum();
    let remaining = opens_after - need;

    let budget = Budget {
        open_valence,
        fragment_hydrogens,
        delta_h,
        need,
        remaining,
        delta_charge,
        fragment_charge,
        charge_shift,
    };
    let bail = |why: Unresolved| Byproducts {
        molecules: Vec::new(),
        verdict: Verdict::Unresolved(why),
        budget,
    };

    // `delta_h` 就是副产物应有的氢数,负数在物理上讲不通。排在最前面判是因为
    // 它**指向的东西最具体**:记录里少写了供氢的试剂。让它掉进后面那些档,
    // 报出来的就只剩"配不成对",线索没了。
    if delta_h < 0 {
        return bail(Unresolved::HydrogenBudgetNegative);
    }
    if !charges_ok {
        return bail(Unresolved::NoPairing);
    }
    if remaining < 0 {
        return bail(Unresolved::BudgetExceedsValence);
    }
    if remaining % 2 != 0 {
        return bail(Unresolved::OddValence);
    }
    let to_bond = u32::try_from(remaining).unwrap_or(0) / 2;
    if to_bond > MAX_BONDS {
        return bail(Unresolved::TooManyBonds);
    }

    // `need` 为负 = 片段带的氢比预算多,得摘掉。摘一个氢就腾出一处空价,
    // 它与原有的空价成 π 键 —— 这正是消除的形状(叔丁酯的叔丁基 → 异丁烯)。
    if need < 0 {
        let extra = u32::try_from(-need).unwrap_or(0);
        if !borrow_hydrogens(&frag, &mut sites, extra, stats) {
            return bail(Unresolved::OddValence);
        }
    }

    let mut closed = frag;
    // 这里失败**不是**"键太多"——上面 `to_bond > MAX_BONDS` 已经拦过了。
    // 走到这里说明剩下的空价物理上配不起来:只剩一个位点,或只剩两个卤素。
    let heavy_before = heavy_atoms(&closed);
    if !form_bonds(&mut closed, &mut sites, to_bond, stats) {
        return bail(Unresolved::NoPairing);
    }
    settle_hydrogens(&mut closed, &sites, stats);

    if omgkit_chem::sanitize(&mut closed).is_err() {
        return bail(Unresolved::FragmentUnsanitizable);
    }
    if let Some(size) = strained_triple_bond(&mut closed) {
        let _ = size;
        return bail(Unresolved::StrainedClosure);
    }
    // 判据:账必须**精确**闭合。这一条不依赖任何记录,是本模块唯一的正确性来源。
    //
    // 重原子那一项要单列。收口只该补氢、落电荷、成键 —— **一个重原子都不该增减**。
    // 前两项(氢、电荷)盯不住它:凭空多一个重原子的同时,氢与电荷完全可以照样
    // 配平。文档一直把"重原子守恒"写成本模块的判据,可它此前只在基准脚本里查,
    // 而 `reconstruct` 是公开 API,调用方拿到的东西没人替他查。
    if heavy_atoms(&closed) != heavy_before
        || total_hydrogens(&closed) != delta_h
        || total_charge(&closed) != budget.delta_charge
    {
        return bail(Unresolved::BudgetMismatch);
    }

    Byproducts {
        molecules: split(&closed),
        verdict: if to_bond == 0 {
            Verdict::Capped
        } else {
            Verdict::Bonded { bonds: to_bond }
        },
        budget,
    }
}

/// 收口有没有把三键塞进小环。有就返回那个环的大小。
///
/// 判据是**几何**,不是化学品味:三键的两端与各自的取代基必须共线,小环里排不下。
/// 门限取 8 —— 环辛炔是已知最小的可分离环炔,再小的(苯炔、环己炔)都是瞬态。
///
/// 为什么非得单独拦:价规则只管"这个原子接了几根键",管不到"这几根键摆得下
/// 摆不下"。所以这类结构原子账平、电荷账平、净化也过得去,前面每一道都放行。
fn strained_triple_bond(mol: &mut MolBuilder) -> Option<u8> {
    const MIN_RING_FOR_TRIPLE: u8 = 8;
    if !mol.bonds().iter().any(|b| b.order == BondOrder::Triple) {
        return None; // 绝大多数分子在这里就走了
    }
    let rings = omgkit_chem::perceive_rings(mol);
    for b in mol.bonds() {
        if b.order != BondOrder::Triple {
            continue;
        }
        // 环大小取两端的最小环里更小的那个;0 表示不在环中,不算
        let sizes = [
            rings.atom_min_ring_size[b.begin as usize],
            rings.atom_min_ring_size[b.end as usize],
        ];
        let in_ring = sizes.iter().copied().filter(|&s| s > 0).min();
        if let Some(size) = in_ring {
            if size < MIN_RING_FOR_TRIPLE {
                return Some(size);
            }
        }
    }
    None
}

/// 分子的总氢:每个原子记着的氢,加上图里作为独立节点存在的 `[H]`。
///
/// 两项都要 —— 本库把 `removeHs` 划在净化之外,显式氢原子是留在图里的,
/// 只数 `num_*_hs` 会漏掉它们。
fn total_hydrogens(mol: &MolBuilder) -> i32 {
    mol.atoms()
        .iter()
        .map(|a| {
            i32::from(a.num_explicit_hs)
                + i32::from(a.num_implicit_hs)
                + i32::from(a.atomic_num == 1)
        })
        .sum()
}

/// 重原子数。氢不算 —— 收口本来就要动氢,拿它当守恒量没有意义。
fn heavy_atoms(mol: &MolBuilder) -> usize {
    mol.atoms().iter().filter(|a| a.atomic_num != 1).count()
}

fn total_charge(mol: &MolBuilder) -> i32 {
    mol.atoms().iter().map(|a| i32::from(a.formal_charge)).sum()
}

/// 把丢弃的原子拼成一张图,并记下每个原子欠了几处价。
///
/// 片段内部的键原样搬来;跨界的键(一端进了产物、一端没进)不搬,它在丢弃这
/// 一侧留下的正是空价。空价按**键级**记 —— 断一根双键欠的是两处,不是一处。
fn build_fragment(reactants: &[MolBuilder], discarded: &[Vec<u32>]) -> (MolBuilder, Vec<Site>) {
    let mut out = MolBuilder::new();
    let mut sites: Vec<Site> = Vec::new();
    // (第几个反应物, 原子下标) → 片段图下标
    let mut index: Vec<Vec<u32>> = reactants
        .iter()
        .map(|m| vec![u32::MAX; m.num_atoms()])
        .collect();

    for (ti, drop_list) in discarded.iter().enumerate() {
        let Some(mol) = reactants.get(ti) else {
            continue;
        };
        for &a in drop_list {
            let Some(&data) = mol.atoms().get(a as usize) else {
                continue;
            };
            let mut carried = data;
            // 映射号是模板内部的东西,搬到副产物上没有意义
            carried.atom_map = 0;
            let idx = out.add_atom_data(carried);
            index[ti][a as usize] = idx;
            debug_assert_eq!(
                idx as usize,
                sites.len(),
                "sites[i] 必须对应片段的第 i 个原子"
            );
            sites.push(Site {
                opens: 0,
                borrowed_h: 0,
            });
        }
    }

    for (ti, drop_list) in discarded.iter().enumerate() {
        let Some(mol) = reactants.get(ti) else {
            continue;
        };
        let gone: Vec<bool> = {
            let mut v = vec![false; mol.num_atoms()];
            for &a in drop_list {
                if let Some(slot) = v.get_mut(a as usize) {
                    *slot = true;
                }
            }
            v
        };
        for b in mol.bonds() {
            let (i, j) = (b.begin as usize, b.end as usize);
            match (gone.get(i), gone.get(j)) {
                (Some(true), Some(true)) => {
                    let mut nb = *b;
                    nb.begin = index[ti][i];
                    nb.end = index[ti][j];
                    // 参照原子的下标在片段里无效,清成哨兵 —— 留着会指到别的原子
                    nb.stereo_atoms = [BondData::NO_STEREO_ATOM; 2];
                    let _ = out.add_bond_data(nb);
                }
                // 跨界:丢弃的那一端欠下这根键的键级
                (Some(true), Some(false)) | (Some(false), Some(true)) => {
                    let inside = if gone[i] { i } else { j };
                    // 配位键的**给体端不占价** —— 电子对是它自己出的,断开之后
                    // 它不欠任何东西。按对称的键级算会给它凭空记一处空价,
                    // 于是一条本来完整的记录被判成收不平。
                    // `valence_contribution_to` 是全库对这件事的唯一真相来源。
                    let owed = if b.order == BondOrder::Dative {
                        let idx = u32::try_from(inside).unwrap_or(u32::MAX);
                        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
                        let v = b.valence_contribution_to(idx).round() as u32;
                        v
                    } else {
                        order_valence(b.order)
                    };
                    let target = index[ti][inside];
                    if let Some(s) = sites.get_mut(target as usize) {
                        s.opens += owed;
                    }
                }
                _ => {}
            }
        }
    }
    rebase_fragment_chirality(reactants, discarded, &index, &mut out);
    (out, sites)
}

/// 切断一根键会改变中心的邻居顺序,手性标记必须跟着换参照系。
///
/// # 不换的后果是**镜像**,而且查不出来
///
/// 标记是相对邻居**存储顺序**说的。原子被切下来时那个邻居从列表里消失,顶上来的
/// 隐式氢按本库的约定占**下标 1**(与解析器的存储约定一致,见 `align_for_rebase`)。
/// 被切邻居原在下标 0 或 2 时这个置换是**奇**的,标记必须翻转;在 1 或 3 时是偶的,
/// 不该翻。一律不翻的话,一半的写法给出镜像分子 —— 而原子数、键、连通性、
/// 原子账全都对,纯拓扑比对与质量守恒判据**都发现不了**。
///
/// 实测:同一个分子的八种等价写法,不重定基时错 4 条。
///
/// # 与产物侧共用同一套机制
///
/// `rebase_chirality` 为产物侧做的是同一件事,所以这里直接复用
/// [`align_for_rebase`] 与 [`permutation_is_odd`],不另写一套 —— 换参照系这件事
/// 只该有一个真相来源,两处各写一份迟早会分叉。
fn rebase_fragment_chirality(
    reactants: &[MolBuilder],
    discarded: &[Vec<u32>],
    index: &[Vec<u32>],
    out: &mut MolBuilder,
) {
    for (ti, drop_list) in discarded.iter().enumerate() {
        let Some(mol) = reactants.get(ti) else {
            continue;
        };
        for &a in drop_list {
            let Some(&dst) = index[ti].get(a as usize) else {
                continue;
            };
            if dst == u32::MAX {
                continue;
            }
            let tag = out.atoms()[dst as usize].chiral_tag;
            if !tag.is_tetrahedral() {
                continue;
            }
            let after: Vec<u32> = out.neighbors(dst).map(|(other, _)| other).collect();
            // 槽位空不空看"**在片段里还连不连着这个中心**",不是看那个原子存不存在
            let slots: Vec<Option<u32>> = mol
                .neighbors(a)
                .map(|(other, _)| {
                    index[ti]
                        .get(other as usize)
                        .copied()
                        .filter(|&p| p != u32::MAX && after.contains(&p))
                })
                .collect();
            let Some((before, aligned)) = align_for_rebase(&slots, &after) else {
                continue;
            };
            if permutation_is_odd(&before, &aligned) == Some(true) {
                if let Some(at) = out.atom_mut(dst) {
                    at.chiral_tag = tag.inverted();
                }
            }
        }
    }
}

/// 键级折算成它占掉的价数。
///
/// 走到这里时键级已经是**凯库勒式**的(见 [`reconstruct`] 里的 kekulize 那一步),
/// 所以不会遇到 `Aromatic`。留着那一支只是兜底 —— 真撞上说明上游漏了 kekulize,
/// 按 1 折至少不会 panic,而账目复核会把它抓出来。
fn order_valence(order: BondOrder) -> u32 {
    match order {
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        _ => 1,
    }
}

/// 把片段该拿到的电荷落到具体原子上。
///
/// # 电荷对空价的作用**由元素定**,不是"负减正加"
///
/// 价由有效原子序数定,所以同样加一个负电荷,氧少一处价(O⁻ 只连一根键)而硼
/// **多**一处(BH₄⁻ 连四根);同样加一个正电荷,氮多一处而碳少一处。按"负减正加"
/// 写死的话,碳正离子与硼酸根这两类会算反 —— 算反的后果不是报错,是给出一个
/// 价数不对的分子,或者把一条完整的记录判成收不平。
///
/// 真相来源是 [`omgkit_chem::valence_shift`],与净化的隐式氢推断同一张价表。
///
/// # 挑哪个原子
///
/// 先按"这个元素拿得住这个电荷吗"排(卤素/氧/硫 > 氮 > 其它),再要求它**确实
/// 需要**这次改动:填空价的那一档只落在还欠着价的原子上。挑错原子只影响写法,
/// 总量由预算定死,所以账不会因此不平。
fn apply_charges(
    frag: &mut MolBuilder,
    sites: &mut [Site],
    shift: i32,
    stats: &mut CloseStats,
) -> bool {
    let want_negative = shift < 0;
    let step: i8 = if want_negative { -1 } else { 1 };
    let mut left = shift.unsigned_abs();
    while left > 0 {
        let mut best: Option<(u8, usize)> = None;
        for (k, site) in sites.iter().enumerate() {
            stats.site_visits += 1;
            let Some(a) = frag.atoms().get(k) else {
                continue;
            };
            let dv = omgkit_chem::valence_shift(
                a.atomic_num,
                a.formal_charge,
                a.formal_charge.saturating_add(step),
            );
            // 价降一格 = 填掉一处空价,那就得真有一处空价可填
            if dv < 0 && site.opens == 0 {
                continue;
            }
            if dv == 0 {
                continue; // 这个元素没有价约束,落上去说明不了什么
            }
            let rank = if want_negative {
                match a.atomic_num {
                    9 | 17 | 35 | 53 => 0,
                    8 | 16 => 1,
                    7 => 2,
                    _ => 3,
                }
            } else {
                match a.atomic_num {
                    7 => 0,
                    8 | 16 => 1,
                    _ => 2,
                }
            };
            let better = match best {
                None => true,
                Some((r, _)) => rank < r,
            };
            if better {
                best = Some((rank, k));
            }
        }
        let Some((_, k)) = best else {
            return false;
        };
        let Some(a) = frag.atom_mut(u32::try_from(k).unwrap_or(u32::MAX)) else {
            return false;
        };
        let dv = omgkit_chem::valence_shift(
            a.atomic_num,
            a.formal_charge,
            a.formal_charge.saturating_add(step),
        );
        a.formal_charge = a.formal_charge.saturating_add(step);
        if dv < 0 {
            sites[k].opens = sites[k].opens.saturating_sub(1);
        } else {
            sites[k].opens += 1;
        }
        left -= 1;
    }
    true
}

/// 摘掉 `extra` 个氢,腾出同样多的空价。
///
/// 优先摘在**与已有空价相邻**的原子上:摘出来的空价与原有空价正好成一根 π 键,
/// 得到的是烯烃而不是一对远隔的自由价。叔丁酯水解就是这一档 —— 叔丁基欠一处价、
/// 带的氢比预算多一个,摘掉相邻甲基上的一个氢,两处空价成双键,给出异丁烯。
fn borrow_hydrogens(
    frag: &MolBuilder,
    sites: &mut [Site],
    extra: u32,
    stats: &mut CloseStats,
) -> bool {
    let mut left = extra;
    for adjacent_first in [true, false] {
        for k in 0..sites.len() {
            stats.site_visits += 1;
            while left > 0 && has_spare_hydrogen(frag, sites, k) {
                let near = frag
                    .neighbors(u32::try_from(k).unwrap_or(u32::MAX))
                    .any(|(other, _)| sites.get(other as usize).is_some_and(|s| s.opens > 0));
                if adjacent_first != near {
                    break;
                }
                sites[k].opens += 1;
                sites[k].borrowed_h += 1;
                left -= 1;
            }
        }
        if left == 0 {
            return true;
        }
    }
    left == 0
}

/// 这个原子还有没有氢可摘。
///
/// 判据是它**从底物带过来的氢**:底物是净化过的,所以显式与隐式两栏加起来就是
/// 真实氢数,两类原子共用同一条判据。摘氢的实现是让新成的键去占掉那处价 ——
/// 隐式氢原子净化时自然少补一个,方括号原子由 [`settle_hydrogens`] 手工扣。
///
/// # 不能按"是不是氢或卤素"来估
///
/// 先前的写法是"非氢非卤就当它填得下"。这对**一个氢都没有**的原子是错的:
/// 醚氧、酯羰基碳、季碳都会被当成能出氢,于是空价落在它们身上,成键之后当场
/// 超价,收口整个失败。实测这一条占未决档的大头 —— Cbz/Boc 这类保护基的离去
/// 片段里,与断点相邻的正好都是不带氢的原子。
fn has_spare_hydrogen(frag: &MolBuilder, sites: &[Site], k: usize) -> bool {
    let Some(a) = frag.atoms().get(k) else {
        return false;
    };
    if a.atomic_num == 1 {
        return false;
    }
    let carried = u32::from(a.num_explicit_hs) + u32::from(a.num_implicit_hs);
    carried > sites[k].borrowed_h
}

/// 把剩余空价两两成键。
///
/// 挑对的优先级,由高到低:
///
/// 1. **两个原子在片段里已经成键** → 提一级键级,得到 π 键(消除)
/// 2. **两个原子分属不同连通分量** → 新建单键(离去基团互相结合,如 Wittig 的 P=O)
/// 3. 同一分量、原先不相邻 → 成环,最不优先
///
/// 卤素之间不成键:它们是离去基团,该拿氢变成卤化氢,而不是两两配成卤素单质。
fn form_bonds(
    frag: &mut MolBuilder,
    sites: &mut [Site],
    to_bond: u32,
    stats: &mut CloseStats,
) -> bool {
    let mut formed = 0;
    // 只有**还欠着价**的位点参与配对。片段可以有上百个原子,而欠价的通常只有
    // 两三个 —— 每成一根键就把所有原子两两扫一遍,正是"在按键的循环里做一件
    // 正比于整个片段的事"那个形状。
    let mut open_sites: Vec<usize> = (0..sites.len()).filter(|&k| sites[k].opens > 0).collect();
    while formed < to_bond {
        stats.fragment_scans += 1;
        let comp = components(frag);
        let mut best: Option<(u32, usize, usize)> = None;
        for (x, &i) in open_sites.iter().enumerate() {
            for &j in &open_sites[x + 1..] {
                stats.site_visits += 1;
                if sites[i].opens == 0 || sites[j].opens == 0 {
                    continue;
                }
                let (a, b) = (
                    u32::try_from(i).unwrap_or(u32::MAX),
                    u32::try_from(j).unwrap_or(u32::MAX),
                );
                if is_halogen(frag, a) && is_halogen(frag, b) {
                    continue;
                }
                let score = if frag.bond_between(a, b).is_some() {
                    0
                } else if comp.get(a as usize) != comp.get(b as usize) {
                    1
                } else {
                    2
                };
                // 写成 match 而不是 `is_none_or` —— 后者要 Rust 1.82,本工作区
                // 的 MSRV 是 1.75,clippy 有闸门盯着
                let better = match best {
                    None => true,
                    Some((s, ..)) => score < s,
                };
                if better {
                    best = Some((score, i, j));
                }
            }
        }
        let Some((_, i, j)) = best else {
            return false;
        };
        let (a, b) = (
            u32::try_from(i).unwrap_or(u32::MAX),
            u32::try_from(j).unwrap_or(u32::MAX),
        );
        if let Some(bi) = frag.bond_between(a, b) {
            let Some(mut edge) = frag.bond_mut(bi) else {
                return false;
            };
            let Some(up) = raise(edge.get().order) else {
                return false;
            };
            edge.set_order(up);
        } else if frag.add_bond(a, b, BondOrder::Single).is_err() {
            return false;
        }
        sites[i].opens -= 1;
        sites[j].opens -= 1;
        open_sites.retain(|&k| sites[k].opens > 0);
        formed += 1;
    }
    true
}

fn is_halogen(frag: &MolBuilder, idx: u32) -> bool {
    frag.atoms()
        .get(idx as usize)
        .is_some_and(|a| matches!(a.atomic_num, 9 | 17 | 35 | 53 | 85))
}

fn raise(order: BondOrder) -> Option<BondOrder> {
    match order {
        BondOrder::Single | BondOrder::Aromatic => Some(BondOrder::Double),
        BondOrder::Double => Some(BondOrder::Triple),
        _ => None,
    }
}

/// 把剩下的空价交给氢。
///
/// 两类原子走两条路,不能混:
///
/// - **隐式氢原子**(没有 `NO_IMPLICIT`):什么都不用做。净化会按价规则把空缺
///   补成隐式氢 —— 断掉一根键,补一个氢,正是要的结果。
/// - **方括号原子**(有 `NO_IMPLICIT`):氢数是写死的,净化不会动它,必须在这里
///   自己加。`[OH:3]` 被切下来时带着 1 个显式氢、欠 1 处价,加成 2 个就是水。
///
/// 摘氢那一侧同理:隐式氢原子靠新成的键占掉价、净化自然少补;方括号原子要
/// 手工扣掉。
fn settle_hydrogens(frag: &mut MolBuilder, sites: &[Site], stats: &mut CloseStats) {
    for (i, s) in sites.iter().enumerate() {
        stats.site_visits += 1;
        let Some(a) = frag.atom_mut(u32::try_from(i).unwrap_or(u32::MAX)) else {
            continue;
        };
        if !a.flags.contains(AtomFlags::NO_IMPLICIT) {
            continue;
        }
        let add = u8::try_from(s.opens).unwrap_or(0);
        a.num_explicit_hs = a.num_explicit_hs.saturating_add(add);
        let take = u8::try_from(s.borrowed_h).unwrap_or(0);
        a.num_explicit_hs = a.num_explicit_hs.saturating_sub(take);
    }
}

/// 按连通分量切开 —— 与产物侧同一条原则:分子数由连通性决定。
fn split(mol: &MolBuilder) -> Vec<MolBuilder> {
    let comp = components(mol);
    let n_comp = comp.iter().copied().max().map_or(0, |m| m as usize + 1);
    let mut out: Vec<MolBuilder> = (0..n_comp).map(|_| MolBuilder::new()).collect();
    let mut local = vec![u32::MAX; mol.num_atoms()];
    for (a, &c) in comp.iter().enumerate() {
        local[a] = out[c as usize].add_atom_data(mol.atoms()[a]);
    }
    for b in mol.bonds() {
        let c = comp[b.begin as usize] as usize;
        let mut nb = *b;
        nb.begin = local[b.begin as usize];
        nb.end = local[b.end as usize];
        nb.stereo_atoms = [BondData::NO_STEREO_ATOM; 2];
        let _ = out[c].add_bond_data(nb);
    }
    out
}
