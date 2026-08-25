//! **手性:有符号四点体积,以及嵌入之后那一次全局定向。**
//!
//! # 为什么全局定向必须在这里定,不能交给精修
//!
//! 嵌入给出的坐标系定向是**任意的** —— 一半的概率拿到镜像(见
//! [`crate::linalg::symmetric_eigen`] 里符号约定那一段)。而翻转整体手性是一次
//! **反射**(`det = −1`),不在 `SO(3)` 的连通分支里:要用连续下降从一个手性走到
//! 它的镜像,必须经过"所有手性体积同时为 0"的构型,也就是把整个分子压平。
//! **下降法不会付这个势垒**,手性罚项的权重调多大都没用。
//!
//! 所以定向只能**离散地**定一次:两个镜像各算一遍总罚,取小的那个。
//! 这一步确定、便宜(只与手性中心数有关),而且必须做。
//!
//! (RDKit 一有手性中心就把分子嵌到四维,`Embedder.cpp:1632`,就是为了绕开这件事:
//! 四维里在 `(x₃, x₄)` 平面转 π 就把 `x₃` 送到 `−x₃`,而**四维两两距离精确不变** ——
//! 三维里的反射在四维里是一次免费的连续旋转。)
//!
//! # 符号约定 —— **有两个体积,号相反,别混**
//!
//! | 名字 | 式子 | 谁在用 |
//! |---|---|---|
//! | 四配体 | `det[l₁−l₀, l₂−l₀, l₃−l₀]` | [`signed_volume`];`omgkit-depict` 的 `read_chirality` |
//! | 中心基点 | `det[l₀−c, l₁−c, l₂−c]` | [`center_volume`]、[`Center::sign`]、[`correct_count`] |
//!
//! **[`Center::sign`] 说的是后者的号:`@` → 正、`@@` → 负。**
//! 前者正好相反(`@` → 负、`@@` → 正)—— 正四面体上 `V_配体 = −4·V_中心`。
//!
//! 两个体积**不是同一个量**:四配体那个完全不看中心原子在哪,所以中心被挤到
//! 配体四面体外面(伞形翻转)时它一点变化都没有,而真实立体化学已经翻了。
//! 判"号对不对"必须用中心基点那个,理由见 [`center_volume`]。
//!
//! `omgkit-depict` 的 `read_chirality` 仍用四配体口径,那边由
//! `the_reference_tetrahedron_pins_the_sign` 钉住,并有外部判官
//! (`harness/check_wedge_readback.py`,拿 RDKit 从导出的 molblock 反读)验过。
//! 它与这里**不冲突**,因为它读的是 2D + 楔形、产出的是 `ChiralTag` 而不是号 ——
//! 但两边都叫"有符号体积"而号相反,这一段就是为了不让下一个人踩进去。
//!
//! # 抽中心那一半:**槽位约定是实测出来的,不是推的**

//!
//! "四个配体按什么槽位排"这个约定,排错的话**错法整批一致** —— 于是"符号正确率"
//! 要么 0% 要么 100%,两个数看起来都像"约定定死了"。这种错推不出来,只能拿真值比。
//!
//! 真值取自 `harness/dump_chirality.py`:每个中心的有符号体积**在真实三维构象上
//! 的实际符号**(不是标记推出来的号 —— 那正是待验的东西)。实测:
//!
//! | 标记 | 中心基点体积 | 样本 |
//! |---|---|---|
//! | `@`(Ccw) | **正** | 127 / 127 |
//! | `@@`(Cw) | **负** | 120 / 120 |
//!
//! (头一版这张表是四配体口径的 `@`→负 22/22、`@@`→正 17/17。换基点之后
//! 号整体反过来,上面这一组是**同一批基准换个式子重算**出来的,不是重新猜的。)
//!
//! 与 `omgkit-depict` 的约定一致,而且那一致性写成了机器可验的断言。
//! 全量判官见 `examples/chiral_oracle.rs`:247 个中心,符号错 0、漏抽 0。
//! 变异验证过它抓得住整体翻号(那一下是 475/475 全错,不是"差不多对")。
//!
//! # 调用方要负责的那一条
//!
//! [`centers`] 要求**标记与当前键序一致**:按邻居迭代顺序取四个配体,正好是
//! `chiral_tag` 所指的槽位顺序。
//!
//! ## "解析 SMILES → 补氢 → 直接调 `centers`" 是**对**的
//!
//! 这里先前写着那条路"缺一步槽位重排",**那是错的,而且是有害的错** ——
//! 谁照着加一步重排,就会把每个带隐式氢的中心翻成对映体。
//!
//! 理由:`omgkit_io::smiles::parse` **已经**把 `chiral_tag` 归一化成"相对存储序,
//! 且**隐式氢不参与置换**"(见 `omgkit-io/src/smiles.rs` 里那张表与
//! `needs_tag_inversion` 的两个特判:片段首原子带一个显式氢、以及带一个环闭合数)。
//! `add_explicit_hs` 把氢追加到邻居表末尾,正好就是归一化后的标记所期望的位置。
//!
//! 实测,拿 RDKit 从我们交付的坐标读回立体化学
//! (`examples/dump_conformers.rs` + `harness/verify_stereo.py`):
//!
//! | 语料 | 一致 |
//! |---|---|
//! | `harness/corpus/large.smi` 里 **642** 个带立体标记的分子 | **640 / 640** |
//! | `harness/corpus/stereo_edge.smi`(专挑槽位边界) | **21 / 21** |
//!
//! `stereo_edge.smi` 那一份是为这条专门造的:手性原子写在**片段开头**
//! (隐式氢排书写序第一)、写在中间、四个显式配体、以及带环闭合数的,
//! 每档都给出同一分子的两种写法。
//!
//! 另有 2 个分子**判官够不着**(三配位磷,RDKit 的 `AssignStereochemistryFrom3D`
//! 不给它赋手性,连 RDKit 自己嵌出来的构象都读不回)—— 那一档单独计数、
//! 单独设上限闸,不混进失配。
//!
//! 按一张四配位齐全的连接表建分子则天然满足这两条(判官走的正是那条路)。

/// 一个四面体手性中心:四个配体的槽位顺序,以及有符号体积**该是什么号**。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Center {
    /// 中心原子。
    pub atom: u32,
    /// 四个配体原子,**按槽位顺序**(见模块文档)。
    pub ligands: [u32; 4],
    /// 目标符号,说的是 [`center_volume`]
    /// (**以中心原子为基点**)该有的号:`+1.0` 对应 `@`、`-1.0` 对应 `@@`。
    ///
    /// **别按四配体行列式理解 —— 那个的号正好相反。** 这一行先前写的就是
    /// 四配体口径,与代码相反;照它手写 `Center` 会拿到对映体,而所有判据
    /// 都会报"号对"。
    pub sign: f64,
}

impl Center {
    /// 槽位上是**一对孤对电子**(或隐式氢),没有对应的原子。
    ///
    /// 只会出现在 `ligands[3]`:三配位中心按约定把孤对存在最后一格,
    /// 见 [`centers`] 里那段推导。
    pub const IMPLICIT: u32 = u32::MAX;

    /// 这是个**三配位**中心(第四个"取代基"是孤对电子)。
    #[must_use]
    pub fn is_three_coordinate(&self) -> bool {
        self.ligands[3] == Self::IMPLICIT
    }

    /// 真正落在原子上的配体 —— 三配位时只有三个。
    #[must_use]
    pub fn real_ligands(&self) -> &[u32] {
        if self.is_three_coordinate() {
            &self.ligands[..3]
        } else {
            &self.ligands
        }
    }
}

/// 从分子里抽出四面体手性中心。
///
/// # 前置条件 —— 调用方负责,判官验的就是这两条
///
/// 1. **氢必须已经是显式原子。** 中心得凑够四个邻居才谈得上四面体。
/// 2. **`chiral_tag` 必须与当前的键序一致**,也就是"按邻居迭代顺序取四个配体"
///    正好是标记所指的槽位顺序。
///
/// 这两条**从 SMILES 走过来是满足的**:解析器已经把 `chiral_tag` 归一化成
/// "相对存储序、隐式氢不参与置换",而 `add_explicit_hs` 把氢追加到邻居表末尾,
/// 正好对上。详细理由与实测见模块文档那一节 —— 那里先前写着这条路"缺一步
/// 槽位重排",**是错的**,照着加一步重排会把每个带隐式氢的中心翻成对映体。
///
/// 按一张**四配位齐全的连接表**建出来的分子同样满足:没有隐式氢,键序就是槽位序。
/// 判官 `examples/chiral_oracle.rs` 走的是后一条路(真值取自真实构象的有符号体积),
/// `examples/dump_conformers.rs` 走的是前一条。
///
/// 抽不出四个邻居、或者标记不是四面体的原子,直接跳过 —— 不猜。
#[must_use]
pub fn centers(mol: &omgkit_core::MolBuilder) -> Vec<Center> {
    use omgkit_core::ChiralTag;
    let mut out = Vec::new();
    for (idx, a) in mol.atoms().iter().enumerate() {
        // **约定随基点一起翻了。** 头一版按四配体行列式标定:`@` → 负、`@@` → 正
        // (语料实测 22/22、17/17)。[`center_volume`] 改成以中心原子为基点之后,
        // `V_配体 = −4·V_中心`,于是符号整体反过来:`@` → **正**、`@@` → **负**。
        // 这不是重新猜的,是同一批标定样本换个式子重算出来的。
        let sign = match a.chiral_tag {
            ChiralTag::Ccw => 1.0,
            ChiralTag::Cw => -1.0,
            _ => continue,
        };
        let Ok(id) = u32::try_from(idx) else { continue };
        let nb: Vec<u32> = mol.neighbors(id).map(|(y, _)| y).collect();
        let ligands = match nb.len() {
            4 => [nb[0], nb[1], nb[2], nb[3]],
            // **三配位 + 一对孤对**:亚砜/亚磺酰胺的 S、膦的 P …… 构型照样确定。
            // 先前这里 `<[u32; 4]>::try_from` 凑不够四个就整个 `continue`,
            // 于是这些中心的构型是掷硬币 —— 语料 13 个分子、16 个中心。
            //
            // 孤对按约定占**槽位 1**(与 `omgkit-depict::read_chirality` 同一条),
            // 所以四元组是 `[n₀, 孤对, n₁, n₂]`。这里存成 `[n₀, n₁, n₂, IMPLICIT]`:
            // 把孤对从槽位 1 挪到槽位 3 是个 **3-轮换,偶置换,号不变**,
            // 于是 `center_volume` 用前三个槽位算出来的号与四配位那一档**同一个约定**。
            //
            // (另一条等价的算法:四配体里省掉第 `k` 个,行列式的号正比于
            //  `(−1)^(k+1)` —— 正四面体上省 0/1/2/3 实测号是 **− / + / − / +**,交替。
            //  省掉槽位 1 与省掉槽位 3 同号,与上面的结论一致。
            //  头一版这里写的是"+4 / +4 / −4 / −4",与同一句的 `(−1)^(k+1)`
            //  自相矛盾 —— 结论对,四个数是错的,独立审核复现时拆穿的。)
            3 if omgkit_core::element::has_stereogenic_lone_pair(a.atomic_num, a.formal_charge) => {
                [nb[0], nb[1], nb[2], Center::IMPLICIT]
            }
            _ => continue, // 氢没补、或者根本不是四面体中心
        };
        out.push(Center {
            atom: id,
            ligands,
            sign,
        });
    }
    out
}

/// 四个点围出的**有符号**体积(实为体积的 6 倍,只用它的号与相对大小)。
///
/// `det[p₁−p₀, p₂−p₀, p₃−p₀]`。只吃前三维 —— 即便将来在四维里精修,
/// 手性也只由前三维定(RDKit 的 `calcChiralVolume` 同此)。
#[must_use]
pub fn signed_volume(p0: [f64; 3], p1: [f64; 3], p2: [f64; 3], p3: [f64; 3]) -> f64 {
    let d = |p: [f64; 3]| [p[0] - p0[0], p[1] - p0[1], p[2] - p0[2]];
    let (a, b, c) = (d(p1), d(p2), d(p3));
    a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
}

/// 一个中心在给定坐标下的有符号体积:**以中心原子为基点**,
/// `det[l₀−c, l₁−c, l₂−c]`。
///
/// # 为什么基点必须是中心原子,不是第一个配体
///
/// 先前这里用的是 `signed_volume(l₀, l₁, l₂, l₃)` —— 四个配体的行列式,
/// **完全不看中心原子在哪**(`Center::atom` 那时在整个 crate 里从未被读过)。
///
/// 中心取四配体质心时两者恰好差 `−4` 倍(`V_配体 = −4·V_中心`,正四面体上可手算),
/// 所以平时看不出区别。但中心原子被挤到配体四面体**外面**去的时候(伞形翻转),
/// `V_中心` 变号而 `V_配体` **一点变化都没有** —— 于是判据说"号对",
/// 而那组坐标是对映体。
///
/// RDKit 的 `assignChiralTypesFrom3D` 用的正是中心基点。
///
/// # 这一档现在**没有在发生**,换过来是把洞堵上,不是修一个正在漏的洞
///
/// 实测(`large.smi`,484 个中心,在**我们交付的坐标**上):
/// 中心原子在配体四面体外的 **0 个**、号与目标不符的 **0 个**。
/// 距离项其实已经间接挡住了大半 —— 中心跑到外面就必然拉长某条中心–配体键
/// 或压缩某个 1-3 距离,那两档罚得很贵。
///
/// **所以别把这段写成"修好了 N 个对映体"。** 先前的注释写过"484 个里 21 个在
/// 外面、2 个号已翻",那是从一份没有复核的报告里转述的,复现不出来。
/// 换基点的理由是:四配体行列式对这一档**在数学上就是瞎的**,而代价实测为零。
/// 变异验证说明这一档真出事就是灾难:给 300 个分子各翻一个中心的伞,
/// 外部判据 `verify_stereo.py` 从 290/301 掉到 **2/301**。
///
/// 判据先前也看不见它,因为 `harness/dump_chirality.py` 的真值用的**是同一个
/// 四配体行列式**(注释里自己写着"与 `chiral.rs::signed_volume` 同一个式子")。
/// 那已经一并换成中心基点。
///
/// # Panics
///
/// 配体下标越界时 panic。
#[must_use]
pub fn center_volume(coords: &[[f64; 3]], c: &Center) -> f64 {
    let o = coords[c.atom as usize];
    let d = |k: usize| {
        let p = coords[c.ligands[k] as usize];
        [p[0] - o[0], p[1] - o[1], p[2] - o[2]]
    };
    let (a, b, e) = (d(0), d(1), d(2));
    a[0] * (b[1] * e[2] - b[2] * e[1]) - a[1] * (b[0] * e[2] - b[2] * e[0])
        + a[2] * (b[0] * e[1] - b[1] * e[0])
}

/// 一组中心里有几个的号是对的。
///
/// # `v != 0.0` 是精确等零,而且**必须**是精确等零
///
/// 这一行看着像典型的浮点味道:体积算成 1e-17 就不等于 0,于是拿到一个
/// 噪声号还照样计数。距离几何**认不出镜像**,嵌出来的坐标里大批中心恰好压平,
/// 实测 `large.smi`(311 个分子、500 个中心)在 [`needs_reflection`] 看到的
/// 那一刻:**98 个中心 `|V| < 1e-6`**,最小 1.7e-30 —— 看上去近两成的票是
/// 掷硬币。(交付坐标上一个都没有,最小 `|V| = 0.795`。)
///
/// **但那个"修法"是退化。** 实测扫过无量纲的平度阈值
/// `|V| / (|a||b||c|)`(标准四面体中心约 0.27),在 `smoke.chirality`
/// 的 247 个中心上:
///
/// | 阈值 | 嵌完直接对 | 全局反射之后 | 翻了几个分子 |
/// |---|---|---|---|
/// | **0(现状)** | 129(52.2%) | **207(83.8%)** | 59 |
/// | 1e-30 / 1e-18 | 129(52.2%) | **207(83.8%)** | 59 |
/// | 1e-12 ~ 1e-3 | 108(43.7%) | 170(68.8%) | 43 |
/// | 1e-2 | 105(42.5%) | 163(66.0%) | 39 |
///
/// 也就是说**那些微小的号不是噪声**:它是嵌出来的结构已经在往哪边偏,
/// 精修会顺着放大它,而不是把它翻过去。把它们排除出投票,只是让票更容易
/// 打平 —— 而平局的规则是"不翻"。
///
/// 这不是纸上推的:`chiral_oracle` 的 80% 闸正好卡在两者之间
/// (68.8% < 80% < 83.8%),加容差当场退 1。要动这一行,先把上表重跑一遍。
#[must_use]
pub fn correct_count(coords: &[[f64; 3]], centers: &[Center]) -> usize {
    centers
        .iter()
        .filter(|c| {
            let v = center_volume(coords, c);
            v != 0.0 && v.signum() == c.sign
        })
        .count()
}

/// **嵌入之后定一次全局定向:要不要把整个结构镜像过来。**
///
/// 判据是"镜像之后号对的中心更多"。镜像把每个中心的有符号体积整体变号,
/// 所以这一步只是数一遍,不需要真的把坐标翻过去再算。
///
/// 平局时返回 `false`(不翻)——**必须有个确定的平局规则**,否则同一个分子
/// 两次跑可能给出互为镜像的答案,而这个 crate 承诺确定性。
#[must_use]
pub fn needs_reflection(coords: &[[f64; 3]], centers: &[Center]) -> bool {
    let ok = correct_count(coords, centers);
    // 镜像后每个体积变号,所以"号对的"恰好换成"号错且非零的"。
    // 这里的"非零"同样是精确等零,理由见 [`correct_count`] 的那张扫描表。
    let nonzero = centers
        .iter()
        .filter(|c| center_volume(coords, c) != 0.0)
        .count();
    nonzero - ok > ok
}

/// 把整个结构镜像过来(翻第一根坐标轴)。
///
/// 翻哪一根都一样 —— 差别只是一次旋转,而旋转不改变分子。
pub fn reflect(coords: &mut [[f64; 3]]) {
    for p in coords.iter_mut() {
        p[0] = -p[0];
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 一个标准四面体:`p0` 在 +z,另外三个在下方平面上,按 0°/120°/240° 排。
    ///
    /// 从 `p0` 往中心看(即沿 −z 方向看),`p1 → p2 → p3` 是**逆时针**,
    /// 所以这就是 `@`([`ChiralTag::Ccw`](omgkit_core::ChiralTag::Ccw))。
    fn reference_tetrahedron() -> [[f64; 3]; 4] {
        let r = (8.0_f64).sqrt() / 3.0;
        let z = -1.0 / 3.0;
        [
            [0.0, 0.0, 1.0],
            [r, 0.0, z],
            [
                r * (2.0 * std::f64::consts::PI / 3.0).cos(),
                r * (2.0 * std::f64::consts::PI / 3.0).sin(),
                z,
            ],
            [
                r * (4.0 * std::f64::consts::PI / 3.0).cos(),
                r * (4.0 * std::f64::consts::PI / 3.0).sin(),
                z,
            ],
        ]
    }

    #[test]
    fn 参照四面体把符号钉死() {
        // **这条测试是整个模块的符号约定本身。** 它必须与 omgkit-depict 的
        // `the_reference_tetrahedron_pins_the_sign` 给出同一个号 —— 两套约定
        // 并存的话,迟早在某个交界处翻出对映体。
        let t = reference_tetrahedron();
        let v = signed_volume(t[0], t[1], t[2], t[3]);
        assert!(v < 0.0, "逆时针(@ / Ccw)的有符号体积应当是负的,实际 {v}");
        // 换两个配体 = 一次对换 = 手性翻转
        let v2 = signed_volume(t[0], t[2], t[1], t[3]);
        assert!(v2 > 0.0, "对换两个配体之后应当变号,实际 {v2}");
        assert!((v + v2).abs() < 1e-12, "一次对换只该变号,大小不变");

        // **把与 omgkit-depict 的一致性变成机器可验的,而不是一句注释。**
        // 下面这四个点是 depict 的 `the_reference_tetrahedron_pins_the_sign`
        // 里逐字写着的那一组(它断言 `det < 0` ⇒ Ccw)。两个 crate 各有一套
        // 约定却没人比过,迟早在交界处翻出对映体 —— 所以这里直接比。
        let depict_pts = [
            [0.0, 0.0, 1.0],
            [1.0, 0.0, -0.33],
            [-0.5, 0.866, -0.33],
            [-0.5, -0.866, -0.33],
        ];
        let dv = signed_volume(depict_pts[0], depict_pts[1], depict_pts[2], depict_pts[3]);
        assert!(
            dv < 0.0,
            "与 omgkit-depict 的符号约定对不上了:同一组点它判 Ccw(det<0),这里得 {dv}"
        );
    }

    #[test]
    fn 镜像把每个体积都变号() {
        let t = reference_tetrahedron();
        let mut m = t;
        reflect(&mut m);
        let (a, b) = (
            signed_volume(t[0], t[1], t[2], t[3]),
            signed_volume(m[0], m[1], m[2], m[3]),
        );
        assert!((a + b).abs() < 1e-12, "镜像后应当恰好变号:{a} vs {b}");
    }

    #[test]
    fn 旋转不改变符号() {
        // 有符号体积在**旋转**下不变、在**反射**下变号 —— 这正是它能当手性判据的理由。
        let t = reference_tetrahedron();
        let (c, s) = (0.6_f64, 0.8_f64);
        let rot = |p: [f64; 3]| [c * p[0] - s * p[1], s * p[0] + c * p[1], p[2]];
        let r = [rot(t[0]), rot(t[1]), rot(t[2]), rot(t[3])];
        let (a, b) = (
            signed_volume(t[0], t[1], t[2], t[3]),
            signed_volume(r[0], r[1], r[2], r[3]),
        );
        assert!((a - b).abs() < 1e-12, "旋转不该改变有符号体积:{a} vs {b}");
    }

    /// 一个**完整**的中心:下标 0 是中心原子(在质心),1..=4 是四个配体。
    ///
    /// 先前这里的夹具写的是 `atom: 0, ligands: [0,1,2,3]` —— 中心原子就是第一个
    /// 配体,是个不存在的构型。旧公式不看中心原子,所以这个荒唐取值一直没暴露;
    /// 它也正是"全模块没有一条测试碰过 `Center::atom`"的由来。
    fn reference_center() -> Vec<[f64; 3]> {
        let mut v = vec![[0.0, 0.0, 0.0]];
        v.extend_from_slice(&reference_tetrahedron());
        v
    }

    fn ctr(sign: f64) -> Center {
        Center {
            atom: 0,
            ligands: [1, 2, 3, 4],
            sign,
        }
    }

    #[test]
    fn 参照四面体的中心基点体积号为正() {
        // 这一条钉住约定本身:`reference_tetrahedron` 的**四配体**行列式为负
        // (旧口径),而**中心基点**的行列式为正 —— 两者反号,不是同一个量。
        let c = reference_center();
        let vl = signed_volume(c[1], c[2], c[3], c[4]);
        let vc = center_volume(&c, &ctr(1.0));
        assert!(vl < 0.0, "四配体行列式该是负的,实得 {vl}");
        assert!(vc > 0.0, "中心基点行列式该是正的,实得 {vc}");
        // 正四面体上恰好差 −4 倍(一般构型**不**成立,见 `center_volume` 的文档)
        assert!(
            (vl / vc + 4.0).abs() < 1e-9,
            "正四面体上 V_配体 该是 −4·V_中心:{vl} / {vc}"
        );
    }

    #[test]
    fn 全局定向按多数决() {
        let coords = reference_center();
        // 参照四面体的**中心基点**体积是正的
        assert!(!needs_reflection(&coords, &[ctr(1.0)]), "号已经对了,不该翻");
        assert!(needs_reflection(&coords, &[ctr(-1.0)]), "号反了,应当翻");
        // 两个中心一对一错 —— 平局,规则是**不翻**(必须确定,不能随实现摆动)
        assert!(
            !needs_reflection(&coords, &[ctr(-1.0), ctr(1.0)]),
            "平局时的规则是不翻"
        );
        // 二比一
        assert!(
            needs_reflection(&coords, &[ctr(-1.0), ctr(-1.0), ctr(1.0)]),
            "二比一应当翻"
        );
    }

    /// 把中心原子**沿三个配体所在的平面镜像**过去 —— 伞形翻转的干净构造。
    ///
    /// `V = det[l₁−c, l₂−c, l₃−c]` 是 `c` 的仿射函数,且在那张平面上恒为零,
    /// 所以它正比于 `c` 到平面的**有号距离**。镜像把有号距离取反,于是
    /// `V` **精确变号** —— 这是算出来的,不是试出来的。
    fn flip_center(coords: &mut [[f64; 3]], lig: [usize; 3]) {
        let (p, q, r) = (coords[lig[0]], coords[lig[1]], coords[lig[2]]);
        let sub = |u: [f64; 3], v: [f64; 3]| [u[0] - v[0], u[1] - v[1], u[2] - v[2]];
        let (a, b) = (sub(q, p), sub(r, p));
        let n = [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ];
        let nn = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
        let d = sub(coords[0], p);
        let t = (d[0] * n[0] + d[1] * n[1] + d[2] * n[2]) / nn;
        for k in 0..3 {
            coords[0][k] -= 2.0 * t * n[k];
        }
    }

    #[test]
    fn 中心原子翻到配体四面体外面时号必须跟着翻() {
        // **这就是先前那个洞。** 四个配体一动不动,只把中心原子从四面体内部
        // 镜像到外面 —— 真实立体化学翻了,而四配体行列式**一点变化都没有**。
        let mut coords = reference_center();
        let before = center_volume(&coords, &ctr(1.0));
        let vl_before = signed_volume(coords[1], coords[2], coords[3], coords[4]);
        flip_center(&mut coords, [1, 2, 3]);
        let after = center_volume(&coords, &ctr(1.0));
        let vl_after = signed_volume(coords[1], coords[2], coords[3], coords[4]);
        assert!(
            (vl_before - vl_after).abs() < 1e-12,
            "四配体行列式**不该**变(它正是看不见翻伞的原因):{vl_before} → {vl_after}"
        );
        assert!(
            before * after < 0.0,
            "中心基点体积必须变号:{before} → {after}"
        );
        // 于是判据也跟着说"错了"
        assert_eq!(correct_count(&coords, &[ctr(1.0)]), 0, "翻伞之后号该判错");
    }

    #[test]
    fn 压平的中心不算数() {
        // 四点共面 → 体积为 0 → 这个中心的号无从谈起,两边都不该把它计进去。
        // **这一条正是"判据只查符号会被压平刷绿"的根** —— 号为 0 时
        // `signum()` 给 0,与 ±1 都不等,所以它既不算对也不算错。
        let flat = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        ];
        assert_eq!(correct_count(&flat, &[ctr(1.0)]), 0);
        assert_eq!(correct_count(&flat, &[ctr(-1.0)]), 0);
        assert!(
            !needs_reflection(&flat, &[ctr(1.0)]),
            "全是零体积,翻了也没用"
        );
    }

    #[test]
    fn 没有手性中心时不翻() {
        let coords = vec![[0.0; 3]; 4];
        assert!(!needs_reflection(&coords, &[]));
    }
}
