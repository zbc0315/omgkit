//! 拿整份语料过一遍硬性质,把违例逐条报出来。
//!
//! ```shell
//! cargo run -p omgkit-depict --release --example audit -- harness/corpus/large.smi
//! ```
//!
//! # 改动半径怎么量
//!
//! 第三个参数给一个路径,就会另外落一份**逐图指纹**:
//!
//! ```shell
//! cargo run -p omgkit-depict --release --example audit -- harness/corpus/large.smi 1 /tmp/before.tsv
//! # …改代码…
//! cargo run -p omgkit-depict --release --example audit -- harness/corpus/large.smi 1 /tmp/after.tsv
//! join /tmp/before.tsv /tmp/after.tsv |
//!   awk '{n++} $2 != $3 {d++} END {print "合上", n+0, "张,其中", d+0, "张变了"}'
//! ```
//!
//! **一定要用 `join`,不许用 `diff | grep -c '^<'`。** 后者会多报。实测(ACS 的
//! `bond_spacing_pct` 18.0 → 18.5,两份 dump 的键集完全相同):
//!
//! | | `join` | `diff \| grep -c '^<'` |
//! |---|---:|---:|
//! | 报出来的 | **8393** | 8414 |
//!
//! 多出来的 21 行,逐条查过**两边逐字节相同**(例如 `0000637:ChemDraw_…`),
//! 而真变了的 8393 行一条不漏。**机制是 `diff` 把邻近的改动并进同一个 hunk、
//! 连中间没变的行一起打印** —— 不是"一处差异把后面的行整体错开、于是被数很多遍"
//! (那是先前写在这里的说法,**错的**:LCS 会重新同步,每处真改动只数一遍;
//! 真按那个说法误差该是全文件级的,而不是 21/8393 = 0.25%)。
//!
//! **"合上几张"必须一起打出来,不能只打"几张变了"。** `join` 在两边行数不等时
//! **只比交集,而且一声不吭** —— 一份被截断或写错的文件会让改动半径静静地少报。
//! 实测踩过:一份 240 行的文件顶掉了 17662 行的基线,比出来"240 张变了",
//! 而真相是那两份**没有一张可比**。上面那行 `awk` 把 `n` 一起打出来正是为此:
//! `n` 不等于语料张数就说明这次比较根本不作数。`n+0` 那个 `+0` 不能省 ——
//! 一行都没合上时不写它打出来是个空白,而那正是最该看清的一档。
//!
//! # 它量不到什么
//!
//! 指纹用的是「写法无关」那份图元多重集,而那份东西**故意丢掉了只进渲染参数
//! 的量** —— 线宽、字号、虚楔形的间距,这些不随写法变,判据不该看它们。
//! 代价是:改这类参数**改动半径恒为 0**。实测 ACS 的 `line_width_pt` 0.6 → 0.9
//! (每一张 ACS 图的每一根线都变粗了),两份 dump 逐字节相同。
//!
//! 这不是缺陷,是口径的边界:要量的是**图元摆在哪儿**变没变。改渲染参数请直接
//! 看图,别问这个工具。
//!
//! 指纹是把「写法无关」那条判据用的**图元多重集**([`fingerprint`])压成
//! FNV-1a,不另立口径 —— 见 [`scene_digest`]。
//!
//! # 为什么单元判据不够
//!
//! 单元判据里的分子是**手挑的**,而手挑的时候心里已经有一个模型 —— 挑出来的
//! 正好是那个模型覆盖得到的。先前"环内双键不许画到环外"那条列了八个分子,
//! 全是单环或邻稠的,于是桥环共用键那一整类漏了出去;"楔形窄端"那条列的三个
//! 分子都不含 P(V) 中心,于是"楔形落到双键上被静默丢弃"漏了出去。
//!
//! 全量过一遍不需要预先知道会错在哪。
//!
//! # 查哪几条
//!
//! | | 性质 | 前提 |
//! |---|---|---|
//! | 写法无关 | 换一种 SMILES 写法,画出来的图元一模一样 | 无 |
//! | 环内双键 | 环上双键的两条线都落在某个含它的环里 | 布局没退化 |
//! | 楔形落地 | `Depiction` 里记了几个楔形,画布上就有几个 | 无 |
//! | 楔形可读 | 画出楔形的中心,反读回来就是它该有的构型 | 无 |
//! | 键长全等 | 所有键画出来一样长 | 布局没退化 |
//! | 不出画布 | 每个图元都在画布里 | 无 |
//! | 线端不压字 | 没有线的端点落在标签的字形盒里 | 两端标签塞得下 |
//!
//! 带"布局没退化"前提的两条,是因为退化的坐标本身就不成形状 —— 那时要求图
//! 画得对没有意义,而**退化这件事已经报在 [`Depiction::degraded`] 里了**。

use std::collections::BTreeMap;

use omgkit_core::MolBuilder;
use omgkit_depict::{
    generate,
    geom::{point_in_polygon, Point2},
    label::Label,
    render::{is_squeezed, label_at, scene, touches_glyphs, Primitive, Scene},
    style::Style,
};

/// 换几种写法比对。
///
/// **每一种都必须真的换了存储序,否则这条判据是空过的。** 先前用乘法哈希凑
/// 优先序,实测有 **10.85%** 的改写原样返回(全量语料 2874/26493 次),苯、
/// 乙醇这类分子三次"改写"产出的 SMILES 字符串**一模一样** —— 本 crate 的头号
/// 契约就是被这样一个九分之一空过的搅拌器在验的。
///
/// 换成货真价实的置换之后每一次都算数,于是可以多试几种。
///
/// # 这个数是量出来的,不是拍的
///
/// 先前是 5。**5 种只是抽样,漏掉的写法依赖不会自己冒出来** —— 而写法无关是
/// 本 crate 的头号契约,判据自己的**灵敏度**必须先量。全量语料上把它拉开:
///
/// | 写法数 | 5 | 10 | 30 | 60 | 100 |
/// |---|---:|---:|---:|---:|---:|
/// | 违例 | **201** | 221 | **223** | 223 | 223 |
///
/// **5 种漏报了 22 处(10%)**,30 种就饱和,60 与 100 一个不多。所以取 30。
/// 代价是全量一遍从 ~45s 变成 ~55s。
///
/// 命令行第二个参数可以再调:
///
/// ```shell
/// cargo run -p omgkit-depict --release --example audit -- harness/corpus/large.smi 100
/// ```
const WRITINGS: usize = 30;

/// 每种写法最多试几个种子,去找一个**确实换了存储序**的改写。
///
/// 找不到不算失败:苯那样的分子,所有写法产出的字符串本来就相同 —— 那是
/// 分子太对称,不是判据不给力。这种情形要如实计数,不能悄悄当成通过。
const SEED_TRIES: usize = 8;

/// 开跑前先验搅拌器本身。
///
/// 这条判据守的是**判据的判据**:头号契约靠改写写法来验,改写器要是退化成恒等,
/// 整条判据就静悄悄地空过了 —— 而审计照样打印一个漂亮的数字。先前那个乘法哈希
/// 正是这样,10.85% 的改写原样返回。
///
/// 放在 `main` 开头而不是 `#[cfg(test)]` 里:`cargo test` 默认不跑 example 的
/// 测试,写在那儿等于没写。
fn check_the_shuffler() {
    for n in [2usize, 3, 5, 8, 13, 30, 100] {
        let mut identical = 0usize;
        for seed in 0..64u64 {
            let p = shuffled(n, seed);
            let mut sorted = p.clone();
            sorted.sort_unstable();
            assert!(
                sorted.into_iter().eq(0..u32::try_from(n).expect("n 不大")),
                "n={n} seed={seed}:搅出来的根本不是一个置换"
            );
            if p.iter().enumerate().all(|(i, x)| *x as usize == i) {
                identical += 1;
            }
        }
        // n=2 只有两种置换,恒等占一半是应该的;再大就必须极少见
        let allowed = if n <= 3 { 40 } else { 2 };
        assert!(
            identical <= allowed,
            "n={n}:64 个种子里有 {identical} 个搅出恒等置换 —— 这样的搅拌器验不了写法无关"
        );
    }
}

fn main() {
    check_the_shuffler();
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "harness/corpus/large.smi".into());
    let writings: usize = std::env::args().nth(2).map_or(WRITINGS, |x| {
        x.parse().expect("第二个参数是写法数,要是个整数")
    });
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("读不了 {path}:{e}"));
    // 逐图指纹落到这里 —— 给出路径才算,不给就一个字节都不写。用法见模块文档。
    let dump = std::env::args().nth(3);
    let mut prints: Vec<String> = Vec::new();

    let mut fails: BTreeMap<&'static str, Vec<String>> = BTreeMap::new();
    let mut n_ok = 0usize;
    let mut n_skip = 0usize;
    let mut checked: BTreeMap<&'static str, usize> = BTreeMap::new();
    let mut quality: BTreeMap<&'static str, usize> = BTreeMap::new();

    // **行号从 1 数起。** 先前直接报 `enumerate()` 的 0 基下标,而报出来的东西
    // 是给人拿编辑器去翻的 —— 差一行翻出来就是**另一个分子**。实测代码注释里
    // 那七处"语料第 N 行"全部指错:第 7879 行是个内酯,真正的 Ni 四齿配合物在
    // 7880;第 573 行是个肟,内消旋环己醇在 574。已逐条改正。
    for (lineno, line) in text.lines().enumerate().map(|(i, l)| (i + 1, l)) {
        let smi = line.split_whitespace().next().unwrap_or("");
        // `#` 开头是注释。先前没跳,靠"解析失败"混过去 —— 于是 `跳过 N` 这个数
        // 把注释行和真正解析不了的分子混在了一起,后者就藏进去了。
        if smi.is_empty() || smi.starts_with('#') {
            continue;
        }
        let Some(m) = prep(smi) else {
            n_skip += 1;
            continue;
        };
        if m.num_atoms() < 2 {
            n_skip += 1;
            continue;
        }
        n_ok += 1;
        for style in &Style::ALL {
            let d = generate(&m, style);
            let s = scene(&m, &d, style);
            // **几何判据要拿被画的那个分子。** 为画出构型补的显式氢也在里面,
            // 而 `coords`/`wedges` 的下标是相对它的 —— 拿原分子索引会错位甚至
            // 越界。写法无关那一条例外,见 `checks` 的注释。
            let grown = d.drawn(&m);
            let tag = format!("{}:{lineno}:{smi}", style.name);
            if dump.is_some() {
                // 键排在前、指纹在后,`join` 默认按第一列合。三处细节:
                //
                // - **行号补零**:`join` 按字典序合,不补零 `10` 会排在 `9` 前面;
                // - **规范名里的空格换掉**:有空格就不是一个字段了;
                // - **SMILES 也进键**:只用"行号+规范"的话,拿**另一份语料**的
                //   指纹去比也能配上一堆,比出来的数看着像模像样。带上 SMILES
                //   就在结构上配不上,`n` 直接掉到 0,不必靠人去核对。
                prints.push(format!(
                    "{lineno:07}:{}:{smi} {:016x}",
                    style.name.replace(' ', "_"),
                    scene_digest(&s)
                ));
            }
            let clean = d.degraded.is_empty() && d.unresolved.is_empty();
            // 画得干不干净也要报 —— 判据只说"没画错",不说"画得好"。
            // 交叉单独记一笔:下面那串 if-else 里它排在退化和未解冲突后面,
            // 只统计"只有交叉"的话,真实的交叉数会被前两档吃掉。
            if !d.crossings.is_empty() {
                *quality.entry("—— 其中有键交叉").or_default() += 1usize;
                // 交叉的键里有没有端基键 —— 端基键消冲突翻不动(翻一个端基
                // 等于没翻),要修得另想办法。先量清楚够不够本。
                let terminal = d.crossings.iter().any(|(b1, b2)| {
                    [b1, b2].iter().any(|b| {
                        let bd = &grown.bonds()[**b as usize];
                        grown.degree(bd.begin) == 1 || grown.degree(bd.end) == 1
                    })
                });
                *quality
                    .entry(if d.degraded.is_empty() {
                        "——   布局没退化的"
                    } else {
                        "——   布局已退化的"
                    })
                    .or_default() += 1usize;
                *quality
                    .entry(if terminal {
                        "——   涉及端基键的"
                    } else {
                        "——   两端都不是端基的"
                    })
                    .or_default() += 1usize;
            }
            // **「未解冲突」是拿外接圆量的,而外接圆是上界。** 单独记一档。
            //
            // 消冲突那一步用圆是**对的**,而且是论证过的:它跑在 `orient` 之前,
            // 那时最终朝向还没定,轴对齐的盒六成情形下方向就是错的,圆是"对朝向
            // 一无所知"时唯一诚实的形状(见 `refine::radii`)。
            //
            // 但**报出来**的时候朝向已经定了,这时可以问一句更准的:那两个标签盒
            // 到底叠没叠。全量语料 2839 对未解冲突,两端都有标签的 1430 对里
            // **910 对(64%)盒根本不重叠** —— 那是圆模型的过估,不是画得不好。
            //
            // 只加这一档,不改消冲突的判据 —— 决策端的模型该保守,报告端该准确,
            // 两者本来就不必是同一个。
            if labels_really_overlap(&grown, &d, style) > 0 {
                *quality.entry("——   其中确有两个标签盒叠上").or_default() += 1usize;
            }
            // **「退化」这一档里有多少只是个诚实的标签。**
            //
            // 桥环用了模板表就记一笔退化 —— 那是"这不是从头算出来的坐标",不是
            // "画坏了"。它排在分档的第一位,于是别的毛病会被它吃掉(未解冲突
            // 那一档就少算了 118 张),反过来它自己也会被读成 340 张烂图。
            //
            // 所以单报一笔:退化的图里,除了这个标签之外**一点别的毛病都没有**
            // 的有几张。实测 340 张里 **302 张(88.8%)** 是这样。
            if !d.degraded.is_empty()
                && d.crossings.is_empty()
                && labels_really_overlap(&grown, &d, style) == 0
                && no_atom_sits_on_another(&grown, &d).2.is_none()
            {
                *quality
                    .entry("——   其中只是标着退化,没有别的毛病")
                    .or_default() += 1usize;
            }
            // **两端标签加起来比一个键还长的键**,单独记一档。
            //
            // 它不是"画错了",是 ACS 规范下标签本来就占 0.69 个键长 —— `O⁻—N⁺`
            // 两端要 1.375 个键长的净空,一个键长塞不下,`trim` 只能压缩兜底,
            // 于是线端点落进字里。翻转动不了它(键长是定死的),所以不进违例列;
            // 但也不能不报 —— 那正是"画不好要说出来"该覆盖的东西。
            if squeezed_bonds(&grown, &d, style) > 0 {
                *quality.entry("—— 有标签在键上塞不下").or_default() += 1usize;
            }
            // **骨架原子被摆成 180°**,单独记一档。
            //
            // 渲染这边会给它补一个符号,否则图上根本看不见它(两根键连成一条
            // 直线,顶点没有拐角)。补符号是对的 —— 但那 154 处两根**单键**碰巧
            // 共线的骨架碳,坐标本身就不对:sp3 碳该是 120°,是取代基避让一档
            // 一档挪出来的。补了符号图能读了,布局的毛病却被盖住了,所以这里
            // 单独报一笔。真正的累积双键(sp,本来就该 180°)只占 42/300。
            if accidental_collinear(&grown, &d) > 0 {
                *quality.entry("—— 有骨架原子被摆成 180°").or_default() += 1usize;
            }
            if cramped_substituents(&grown, &d) > 0 {
                *quality.entry("—— 有取代基挤到另一根键上").or_default() += 1usize;
            }
            // **立体中心的取代基挤在一侧**,单独记一档。
            //
            // 这一档**曾经**是正确性问题:隐式氢先前被摆在中心上,而挤在一侧时
            // 它其实落在对面的空扇区 —— 外部判官量到 21 个中心因此画成了对映体。
            // 那个模型已经修了,现在本实现与 RDKit 读出来一致。
            //
            // 留着这一档是因为它仍是**布局质量**的信号:读者得自己想明白"第四个
            // 配体在空扇区里",比取代基摊开的图费劲。
            // **楔形落在环键上**,单独记一档,并分清是不是没得选。
            //
            // IUPAC 的图示建议说立体键该画向取代基;环键的两个原子在读者眼里
            // 都躺在环平面里。补显式氢把这一档从 159 压到 18,剩下的**不是
            // 缺口**:15 根的中心四根键全在环上、又没有氢可补 —— 补不了也没有
            // 合法楔形。**RDKit 对同一批中心也画在环键上、也不补氢**(逐个核过
            // 它的 molblock),所以这是理论下限。
            //
            // 另外 3 根是"两个相邻立体中心抢同一根共用的非环键"那笔取舍:让出去
            // 就有一个中心画不出构型,信息比样式重要,见 `stereo::assign_wedges`。
            let (forced, avoidable) = ring_wedges(&grown, &d);
            if forced > 0 {
                *quality.entry("—— 有楔形只能画在环键上").or_default() += 1usize;
            }
            if avoidable > 0 {
                *quality
                    .entry("——   其中本可避开(抢共用键输了)")
                    .or_default() += 1usize;
            }
            if crowded_centres(&grown, &d) > 0 {
                *quality.entry("—— 有立体中心的取代基挤在一侧").or_default() += 1usize;
            }
            // 桥环退化里,查表命中与没命中要分开数 —— 命中的是一次昂贵搜索
            // 的结果,通常不自交;没命中只有运行时那 5 个初值,实测常常自交。
            //
            // **"没命中"还要再分两种,它们指向完全不同的动作。** 先前混成一档,
            // 而实测全量语料里那 4 例**全是**指纹算不出来 —— 文案让人去补语料,
            // 是条走不通的路。
            for g in &d.degraded {
                let omgkit_depict::rings::Degradation::BridgedRingSystem { template, .. } = g
                else {
                    continue;
                };
                match template {
                    omgkit_depict::templates::Status::Hit => {}
                    omgkit_depict::templates::Status::NotInTable => {
                        *quality
                            .entry("——   其中表里没有(骨架该补进 bridged.smi)")
                            .or_default() += 1usize;
                    }
                    omgkit_depict::templates::Status::NoFingerprint => {
                        *quality
                            .entry("——   其中骨架指纹算不出来(要查那个分子本身)")
                            .or_default() += 1usize;
                    }
                }
            }
            *quality
                .entry(if !d.degraded.is_empty() {
                    "退化(桥环等)"
                } else if !d.unresolved.is_empty() {
                    "有未解冲突"
                } else if !d.crossings.is_empty() {
                    "有键交叉"
                } else {
                    "干净"
                })
                .or_default() += 1usize;

            for (name, hit, bad) in checks(&m, &grown, &d, &s, style, clean, writings) {
                *checked.entry(name).or_default() += usize::from(hit);
                if let Some(why) = bad {
                    fails
                        .entry(name)
                        .or_default()
                        .push(format!("{tag} —— {why}"));
                }
            }
        }
    }

    if let Some(p) = &dump {
        // **`join` 要两边同序,这一行把"有序"从巧合变成不变量。** 现在的键恰好
        // 是按语料行号、再按 `Style::ALL` 的声明序落下来的,而 `ACS_…` 正好
        // 小于 `ChemDraw_…` —— 哪天加一套名字排在 ACS 前面的规范,文件就不再
        // 有序,而 BSD 的 `join` 遇到乱序**不报警**(man 页写着 "join may not
        // report all field matches"),半径会静静地少报。
        prints.sort();
        prints.push(String::new()); // 末尾换行
        std::fs::write(p, prints.join("\n")).unwrap_or_else(|e| panic!("写不了 {p}:{e}"));
        println!("逐图指纹落在 {p}({} 张)\n", prints.len() - 1);
    }
    println!("语料 {path}:解析成功 {n_ok},跳过 {n_skip};每分子比 {writings} 种写法\n");
    let tot: usize = quality
        .iter()
        .filter(|(k, _)| !k.starts_with('—'))
        .map(|(_, v)| *v)
        .sum();
    println!("出图质量({tot} 个分子×规范):");
    for (k, v) in &quality {
        println!("  {k:<16} {v:>6}  {:>5.1}%", 100.0 * *v as f64 / tot as f64);
    }
    println!();
    println!("{:<14} {:>10} {:>8}", "性质", "查到", "违例");
    let mut total = 0usize;
    for (name, n) in &checked {
        let bad = fails.get(name).map_or(0, Vec::len);
        total += bad;
        println!("{name:<14} {n:>10} {bad:>8}");
    }
    for (name, list) in &fails {
        println!("\n=== {name} 的前几例 ===");
        for x in list.iter().take(300) {
            println!("  {x}");
        }
        if list.len() > 300 {
            println!("  …… 另有 {} 例", list.len() - 300);
        }
    }
    println!(
        "\n{}",
        if total == 0 {
            "全部通过".to_string()
        } else {
            format!("共 {total} 处违例")
        }
    );
    if total > 0 {
        std::process::exit(1);
    }
}

/// 一张图的短指纹:把 [`fingerprint`] 那份**图元多重集**打成 FNV-1a。
///
/// # 为什么不另立一套口径
///
/// "两张图算不算同一张"这件事,本仓已经有一个定义了 —— 头号契约「写法无关」
/// 比的就是 [`fingerprint`] 给的图元多重集。改动半径要是另按坐标去数,就会
/// 与那条判据说的不是一回事:坐标变了而图元没变(例如整体平移被画布归一吃掉)
/// 会虚报,楔形换了而坐标没变会漏报。这里只是把同一份东西压成 8 个字节,好写
/// 进文件用 `join` 比。
fn scene_digest(s: &Scene) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for item in fingerprint(s) {
        for b in item.as_bytes() {
            h ^= u64::from(*b);
            h = h.wrapping_mul(0x0000_0100_0000_01b3);
        }
        // **条目之间要有真分隔,否则 `["AB","C"]` 与 `["A","BC"]` 必然同码。**
        // 偏偏取 `0xff` 是有道理的:`fingerprint` 返回的是 `String`,保证合法
        // UTF-8,而 **0xFF 在 UTF-8 里永远不可能出现** —— 于是"条目序列 → 字节流"
        // 是单射,那一类必然碰撞被彻底排除,只剩 2⁻⁶⁴ 量级的一般碰撞。
        // **改成 `0x00` 或空格就悄悄失效了**(NUL 是合法 UTF-8 字符)。
        h ^= 0xff;
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

fn prep(smi: &str) -> Option<MolBuilder> {
    let mut m = omgkit_io::smiles::parse(smi).ok()?;
    omgkit_chem::pipeline::sanitize(&mut m).ok()?;
    omgkit_io::stereo::perceive_bond_stereo(&mut m);
    Some(m)
}

/// 有几处取代基被挤到离另一根键 **< `CRAMPED`** 的地方。
///
/// # 为什么要单独数
///
/// 两根键只差几度,画出来就是一根 —— 读者会**整个漏掉一个取代基**。那是读错
/// 结构,和键交叉同一档,可现有的指标一条都报不出来:`未解冲突` 按**原子间
/// 距离**判,而挤成 6° 的两个原子相距 0.10 个键长,够不着阈值;`键角不过窄`
/// 只查度数 ≤ 3 的原子。
///
/// **这一档是看图看出来的**:樟脑的偕二甲基塌成一根,而所有指标都说没事。
fn cramped_substituents(m: &MolBuilder, d: &omgkit_depict::Depiction) -> usize {
    /// 小于这个角就算挤在一起了。取 15° —— 30° 栅格的一半,正常的两根键
    /// 至少差一档。
    const CRAMPED: f64 = 15.0;

    let mut n = 0usize;
    for a in 0..u32::try_from(m.num_atoms()).expect("原子数超出 u32") {
        let here = d.coords[a as usize];
        let mut angs: Vec<f64> = m
            .neighbors(a)
            .map(|(nb, _)| {
                (d.coords[nb as usize] - here)
                    .angle()
                    .to_degrees()
                    .rem_euclid(360.0)
            })
            .collect();
        if angs.len() < 3 {
            continue; // 度 ≤ 2 的两根键之间只有一个夹角,窄了归 `键角不过窄` 管
        }
        angs.sort_by(|x, y| x.partial_cmp(y).expect("角度不会是 NaN"));
        for k in 0..angs.len() {
            let gap = (angs[(k + 1) % angs.len()] - angs[k]).rem_euclid(360.0);
            if gap < CRAMPED {
                n += 1;
            }
        }
    }
    n
}

/// 报出来的未解冲突里,**两个标签盒在最终朝向下真的叠上**的有几对。
///
/// 见调用处的注释:决策端(`refine`)只能用外接圆,报告端可以用盒。
///
/// 两端都得有标签才问这一句 —— 裸骨架碳根本不画符号,"它俩离得近"是另一档
/// 事(顶点分不分得开),不是"字压到字"。
fn labels_really_overlap(m: &MolBuilder, d: &omgkit_depict::Depiction, style: &Style) -> usize {
    d.unresolved
        .iter()
        .filter(|(a, b)| {
            let (Some(la), Some(lb)) = (
                label_at(m, *a, style, &d.coords),
                label_at(m, *b, style, &d.coords),
            ) else {
                return false;
            };
            let dv = d.coords[*a as usize] - d.coords[*b as usize];
            dv.x.abs() < la.half_w + lb.half_w && dv.y.abs() < la.half_h + lb.half_h
        })
        .count()
}

/// 被摆成 180° 而**本来不该是 180°** 的骨架原子有几个。
///
/// sp 原子(有三键、或两根双键)本来就该 180°,不算。剩下的是布局把 sp3 摆直了
/// —— 渲染那边会补个符号让它看得见,但坐标本身的毛病要单独报,不许被补符号盖住。
fn accidental_collinear(m: &MolBuilder, d: &omgkit_depict::Depiction) -> usize {
    use omgkit_core::BondOrder;
    use omgkit_depict::render::is_collinear;

    (0..u32::try_from(m.num_atoms()).expect("原子数超出 u32"))
        .filter(|a| {
            if !is_collinear(m, *a, &d.coords) {
                return false;
            }
            // sp 判据与 `chains::ideal_angle` 同源:有三键、或两根双键
            let mut doubles = 0usize;
            let mut triple = false;
            for (_, bi) in m.neighbors(*a) {
                match m.bonds()[bi as usize].order {
                    BondOrder::Triple => triple = true,
                    BondOrder::Double => doubles += 1,
                    _ => {}
                }
            }
            !(triple || doubles >= 2)
        })
        .count()
}

/// 落在环键上的楔形:`(没有别的合法选择的, 有别的选择却没用上的)`。
///
/// 见调用处的注释:前者是理论下限(RDKit 同样如此),后者是"两个相邻中心抢
/// 同一根共用非环键"那笔取舍的结果。
fn ring_wedges(m: &MolBuilder, d: &omgkit_depict::Depiction) -> (usize, usize) {
    use omgkit_core::BondOrder;
    let rings = omgkit_chem::sssr::ring_set(m);
    let (mut forced, mut avoidable) = (0usize, 0usize);
    for (bi, w) in d.wedges.iter().enumerate() {
        let Some(narrow) = w.narrow() else { continue };
        if !rings.iter().any(|r| r.bonds.contains(&(bi as u32))) {
            continue;
        }
        let has_alt = m.neighbors(narrow).any(|(_, b)| {
            m.bonds()[b as usize].order == BondOrder::Single
                && !rings.iter().any(|r| r.bonds.contains(&b))
        });
        if has_alt {
            avoidable += 1;
        } else {
            forced += 1;
        }
    }
    (forced, avoidable)
}

/// 取代基挤在一侧的立体中心有几个。
///
/// 三个画出来的邻居全落在中心的**同一侧**(最大空隙 > 180°),中心就在它们围出
/// 的三角形**外面**。
///
/// # 这一档曾经是正确性问题
///
/// 隐式氢先前被摆在**中心上**(只在楔形的反面)。那只在"三个邻居把中心围住"时
/// 成立:四面体的四个键方向之和为零,所以 2D 投影之和也为零 —— 三个都在半平面
/// 里时第四个必然在对面的空扇区,不可能在中心。外部判官
/// (`harness/check_wedge_readback.py`)量到 **21 个中心因此画成了对映体,而且
/// 没进 `unwedged`**。模型已经改成"面内取三个邻居单位方向之和的负值",本实现与
/// RDKit 现在读出来一致(480 个中心一致、0 处不一致)。
///
/// # 现在留着它做什么
///
/// 它仍是**布局质量**的信号:这样的图读者得自己想明白"第四个配体在空扇区里",
/// 比取代基摊开的图费劲。不进违例,只计数。
///
/// 已经报进 `unwedged` 的中心不算 —— 那是如实说过"没画出来"的。
fn crowded_centres(m: &MolBuilder, d: &omgkit_depict::Depiction) -> usize {
    (0..u32::try_from(m.num_atoms()).expect("原子数超出 u32"))
        .filter(|a| {
            // 画出了楔形才谈得上读法
            if !m.neighbors(*a).any(|(_, bi)| {
                d.wedges
                    .get(bi as usize)
                    .and_then(|w| w.narrow())
                    .is_some_and(|n| n == *a)
            }) {
                return false;
            }
            let c = d.coords[*a as usize];
            let mut angs: Vec<f64> = m
                .neighbors(*a)
                .map(|(n, _)| {
                    let v = d.coords[n as usize] - c;
                    v.y.atan2(v.x).to_degrees().rem_euclid(360.0)
                })
                .collect();
            // 四个邻居全画出来的中心不需要摆隐式氢,读法不含糊
            if angs.len() != 3 {
                return false;
            }
            angs.sort_by(|x, y| x.partial_cmp(y).expect("坐标非 NaN"));
            let mut gap = 360.0 - (angs[2] - angs[0]);
            for w in angs.windows(2) {
                gap = f64::max(gap, w[1] - w[0]);
            }
            gap > 180.0
        })
        .count()
}

/// 画出来的标签,一个原子一个 —— 与 `scene` 同源。
fn labels_of(m: &MolBuilder, d: &omgkit_depict::Depiction, style: &Style) -> Vec<Option<Label>> {
    (0..u32::try_from(m.num_atoms()).expect("原子数超出 u32"))
        .map(|a| label_at(m, a, style, &d.coords))
        .collect()
}

/// 两端标签加起来比这根键还长的键 —— 每根键一个标志。
///
/// **口径向实现要,不在这里重写。** 先前这里手抄了一份居中盒的算法,而实现
/// 早已改成偏心盒(整串朝一侧挪了 `dx`,好让元素符号落在原子上),两边算出来
/// 的净空对不上 —— 多报了 547 例(1818 → 1271)。`render::is_squeezed` 因此
/// 是公开的。
fn squeezed(
    m: &MolBuilder,
    d: &omgkit_depict::Depiction,
    style: &Style,
    labels: &[Option<Label>],
) -> Vec<bool> {
    m.bonds()
        .iter()
        .map(|b| {
            is_squeezed(
                d.coords[b.begin as usize],
                d.coords[b.end as usize],
                labels[b.begin as usize].as_ref(),
                labels[b.end as usize].as_ref(),
                style,
            )
        })
        .collect()
}

fn squeezed_bonds(m: &MolBuilder, d: &omgkit_depict::Depiction, style: &Style) -> usize {
    squeezed(m, d, style, &labels_of(m, d, style))
        .iter()
        .filter(|x| **x)
        .count()
}

/// 形状指纹:两两距离排序后的多重集。与原子编号、平移、旋转、镜像都无关。
fn shape(c: &[Point2]) -> Vec<i64> {
    let mut v: Vec<i64> = (0..c.len())
        .flat_map(|i| ((i + 1)..c.len()).map(move |j| (i, j)))
        .map(|(i, j)| (c[i].dist(c[j]) * 1e4).round() as i64)
        .collect();
    v.sort_unstable();
    v
}

/// 量化后的坐标多重集 —— 与原子编号无关。
fn quantised(c: &[Point2]) -> Vec<(i64, i64)> {
    let mut v: Vec<(i64, i64)> = c
        .iter()
        .map(|p| ((p.x * 1e4).round() as i64, (p.y * 1e4).round() as i64))
        .collect();
    v.sort_unstable();
    v
}

/// 图元的多重集。线段不分方向 —— `L A B` 与 `L B A` 是同一条线,而谁是起点
/// 取决于键的 begin/end,本来就随写法变。楔形分方向:窄端宽端不是一回事。
fn fingerprint(s: &Scene) -> Vec<String> {
    let q = |p: Point2| format!("{:.3},{:.3}", p.x, p.y);
    let mut v: Vec<String> = s
        .items
        .iter()
        .map(|it| match it {
            Primitive::Line { from, to, .. } => {
                let (x, y) = (q(*from), q(*to));
                if x <= y {
                    format!("L {x} {y}")
                } else {
                    format!("L {y} {x}")
                }
            }
            Primitive::Wedge { from, to, .. } => format!("W {} {}", q(*from), q(*to)),
            Primitive::Hash { from, to, .. } => format!("H {} {}", q(*from), q(*to)),
            Primitive::Text { at, runs, .. } => format!("T {} {runs:?}", q(*at)),
        })
        .collect();
    v.sort();
    v
}

type Check = (&'static str, bool, Option<String>);

/// `orig` 是**调用方传进来的**分子,`drawn` 是**真正被画的**那个(可能多几个
/// 为画出构型补的显式氢,见 `Depiction::drawn`)。
///
/// **两者不能混用。** 几何判据的下标相对 `drawn`;而写法无关那一条要把分子
/// **重写成 SMILES** 再画一遍,那必须用 `orig` —— 拿补完的去写,显式氢会进到
/// SMILES 里,改写出来的就不是同一个分子了。
fn checks(
    orig: &MolBuilder,
    drawn: &MolBuilder,
    d: &omgkit_depict::Depiction,
    s: &Scene,
    style: &Style,
    clean: bool,
    writings: usize,
) -> Vec<Check> {
    let mut v = vec![
        ring_double_bonds(drawn, d, s, style, clean),
        wedges_reach_canvas(d, s),
        wedges_read_back(drawn, d),
        bond_lengths_equal(drawn, d, clean),
        inside_canvas(drawn, d, s, style),
        lines_clear_of_labels(drawn, d, s, style),
        no_atom_sits_on_another(drawn, d),
        no_angle_is_pinched(drawn, d, clean),
    ];
    // 写法无关出三行:判据本身、比满没有、以及有没有查成 —— 见其文档注释
    v.extend(writing_independent(orig, d, s, style, clean, writings));
    v
}

/// 键角不许被压到 90° 以下。
///
/// **60° 不只是难看** —— 链上出现一个 60° 的拐角,看着像旁边有个三元环,那是
/// 让人读错结构。取代基避让(`chains::free_direction`)按 30° 一档挪,挪两档
/// 就成了 60°,所以这条必须守着。
///
/// 只对没退化的布局下判断:桥环松弛出来的坐标本来就不成形状。
///
/// # 三元环的内角不算,而且不是"网开一面"
///
/// 这一条**曾经**把三元环也算进去,于是全量语料 181 处违例里有 **156 处正好是
/// 60.0°**,全是环丙烷/环氧/氮丙啶的内角。那不是画错了,是**这条判据与
/// [`bond_lengths_equal`] 自相矛盾**:
///
/// 三元环只有三根键,键长全等就是等边三角形,内角**必然恰好 60°**。也就是说
/// 含三元环的分子无论怎么摆都过不了这一条 —— 除非先违反键长全等。两条判据
/// 不能同时满足,而键长全等是画图规范定死的,所以让路的是这一条。
///
/// 判别方法:`a` 的两个邻居**彼此成键**,等价于 `a-i-j` 是个三元环。四元环
/// 不受影响 —— 正方形内角 90°,在 `FLOOR = 89°` 之上。
fn no_angle_is_pinched(m: &MolBuilder, d: &omgkit_depict::Depiction, clean: bool) -> Check {
    if !clean {
        return ("键角不过窄", false, None);
    }
    const FLOOR: f64 = 89.0;
    for a in 0..u32::try_from(m.num_atoms()).expect("原子数超出 u32") {
        let nbrs: Vec<u32> = m.neighbors(a).map(|(n, _)| n).collect();
        // 四配位的理想角就是 90°,五配位更小 —— 这条只管度数 ≤ 3 的
        if nbrs.len() < 2 || nbrs.len() > 3 {
            continue;
        }
        let c = d.coords[a as usize];
        for i in 0..nbrs.len() {
            for j in (i + 1)..nbrs.len() {
                // 三元环的内角,见上面的论证
                if m.neighbors(nbrs[i]).any(|(n, _)| n == nbrs[j]) {
                    continue;
                }
                let u = (d.coords[nbrs[i] as usize] - c).normalized();
                let v = (d.coords[nbrs[j] as usize] - c).normalized();
                let deg = u.dot(v).clamp(-1.0, 1.0).acos().to_degrees();
                if deg < FLOOR {
                    return (
                        "键角不过窄",
                        true,
                        Some(format!(
                            "原子 {a} 处 {}–{a}–{} 的夹角只有 {deg:.1}°",
                            nbrs[i], nbrs[j]
                        )),
                    );
                }
            }
        }
    }
    ("键角不过窄", true, None)
}

/// 两个原子不许画在同一点上。
///
/// 重合本身会被消冲突报成"未解冲突",但那个说法太轻:两个原子叠在一起时,
/// 它们各自的键看起来首尾相接,**图上就多出一个分子里没有的环**。读者没有
/// 任何办法看出那个环是假的 —— 这比"挤了一点"严重得多。
fn no_atom_sits_on_another(m: &MolBuilder, d: &omgkit_depict::Depiction) -> Check {
    const TOL: f64 = 0.05; // 单位是键长
    for i in 0..d.coords.len() {
        for j in (i + 1)..d.coords.len() {
            let dist = d.coords[i].dist(d.coords[j]);
            if dist < TOL {
                let reported = d
                    .unresolved
                    .iter()
                    .any(|(a, b)| (*a as usize, *b as usize) == (i, j));
                return (
                    "原子不重合",
                    true,
                    Some(format!(
                        "原子 {i} 与 {j} 相距 {dist:.4} 个键长{}",
                        if reported {
                            "(已报未解冲突)"
                        } else {
                            "**而且没报出来**"
                        }
                    )),
                );
            }
        }
    }
    let _ = m;
    ("原子不重合", true, None)
}

/// 确定性的伪随机置换(splitmix64 + Fisher–Yates)。
///
/// **不许拿乘法哈希凑。** `(i * K * k) % M` 那种写法在很多 `n` 上根本不是置换:
/// 一堆原子挤到同一个优先级,排序退化成恒等,改写出来的 SMILES 与原文一字不差。
/// 这里给的是货真价实的均匀置换,而且只由 `seed` 决定 —— 与进程、与迭代顺序
/// 都无关,同一份语料每次运行搅出同一批写法。
fn shuffled(n: usize, seed: u64) -> Vec<u32> {
    let mut state = seed;
    let mut next = move || -> u64 {
        state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    };
    let mut p: Vec<u32> = (0..u32::try_from(n).unwrap_or(u32::MAX)).collect();
    for i in (1..n).rev() {
        let j = usize::try_from(next() % (i as u64 + 1)).unwrap_or(0);
        p.swap(i, j);
    }
    p
}

/// 字符串的确定性哈希(FNV-1a),让不同分子拿到不同的置换。
///
/// 不用 `DefaultHasher`:它的种子由标准库决定,跨版本不保证稳定,而这里要的是
/// **换一台机器、换一个版本都搅出同一批写法**。
fn seed_of(s: &str) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for b in s.bytes() {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

/// 换写法画出来必须一模一样。
///
/// 返回三行,因为**"查出来的缺陷"和"根本没查成"是两回事,不许混进一个数**:
///
/// | 行 | 「查到」列的含义 |
/// |---|---|
/// | `写法无关` | 至少比成过一次的 case;「违例」列才是真缺陷 |
/// | `写法无关·比满` | 凑够了 [`WRITINGS`] 次真比较的 case |
/// | `写法无关·没查成` | **一次都没比成**的 case |
///
/// 先前把"没查成"记进违例列,于是换一个更差的搅拌器反而让违例数从 259 涨到
/// 1293 —— 涨的那 1156 全是查不动的 case。数字越大看着越像在认真查,实际上
/// 恰恰相反。
fn writing_independent(
    m: &MolBuilder,
    d0: &omgkit_depict::Depiction,
    s: &Scene,
    style: &Style,
    clean: bool,
    writings: usize,
) -> [Check; 3] {
    let want = fingerprint(s);
    let n = m.num_atoms();
    let canon = omgkit_io::canon::canonical_smiles(m).smiles;
    // 原分子的规范标号。改写之后若逐位相同,说明**存储序压根没变**,这一次
    // 比较是白比的 —— 必须查出来,不能记进 compared。
    let ranks0 = omgkit_io::canon::canonical_ranks(m);
    let base = seed_of(&canon);

    // 真正比过几次、因为"换出来不是同一个分子"跳过了几次、以及试遍种子都换不
    // 出新写法几次。一次都没比过的话这一条是空过的,必须看得见。
    let (mut compared, mut skipped, mut unshuffled) = (0usize, 0usize, 0usize);
    for k in 1..=writings {
        // 试几个种子,直到搅出一个确实换了存储序的写法
        let mut found: Option<(String, MolBuilder)> = None;
        for t in 0..SEED_TRIES {
            let seed = base
                .wrapping_add((k as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15))
                .wrapping_add((t as u64).wrapping_mul(0xD1B5_4A32_D192_ED03));
            let w = omgkit_io::smiles::write_with_priority(m, &shuffled(n, seed));
            let Some(m2) = prep(&w.smiles) else { continue };
            if omgkit_io::canon::canonical_ranks(&m2) == ranks0 {
                continue; // 存储序没变,换个种子再来
            }
            found = Some((w.smiles, m2));
            break;
        }
        let Some((smiles, m2)) = found else {
            unshuffled += 1;
            continue;
        };
        // 换出来的必须还是同一个分子,否则比的是两个东西。
        //
        // **要用规范式比。** `write` 按存储序写,而 m2 的存储序正是被打乱过的
        // 那一套 —— 拿它去比,几乎每个分子都会判成"不是同一个"而跳过,这一条
        // 就静悄悄地什么都没查。
        if canon != omgkit_io::canon::canonical_smiles(&m2).smiles {
            skipped += 1;
            continue;
        }
        compared += 1;
        let d2 = generate(&m2, style);
        let got = fingerprint(&scene(&m2, &d2, style));
        if got != want {
            let diff = want.iter().zip(&got).filter(|(a, b)| a != b).count();
            // **差在哪一层要分开说。** 坐标就变了是布局的问题;坐标一样只是
            // 楔形或键级不同,那是指派的问题 —— 两者的修法完全不同。
            let coords_same = quantised(&d0.coords) == quantised(&d2.coords);
            let wedge_kinds = |d: &omgkit_depict::Depiction| {
                let mut v: Vec<String> = d
                    .wedges
                    .iter()
                    .filter_map(|w| {
                        w.narrow()
                            .map(|_| format!("{w:?}").split('{').next().unwrap_or("").to_string())
                    })
                    .collect();
                v.sort();
                v
            };
            let layer = if !coords_same {
                // 形状指纹:两两距离的多重集,与平移旋转镜像都无关。形状一样
                // 只是摆位不同,病在规范朝向;形状本身变了,病在布局。
                if shape(&d0.coords) == shape(&d2.coords) {
                    "形状相同,摆位不同"
                } else {
                    "形状就变了"
                }
            } else if wedge_kinds(d0) != wedge_kinds(&d2) {
                "坐标相同,楔形不同"
            } else {
                "坐标相同,键级或落点不同"
            };
            return [
                (
                    "写法无关",
                    true,
                    Some(format!(
                        "[{layer}|{}] 写成 {smiles} 之后有 {diff}/{} 处图元不同",
                        if clean {
                            "布局干净"
                        } else {
                            "布局已退化"
                        },
                        want.len()
                    )),
                ),
                ("写法无关·比满", compared == writings, None),
                ("写法无关·没查成", false, None),
            ];
        }
    }
    // **"没查成"单独记一行,不进违例列。** 混进违例列会让"换个更差的搅拌器"
    // 看着像查得更狠 —— 实测旧哈希的 1293 里有 1156 就是这么来的。
    //
    // 但两种"没查成"的分量完全不同,只有后一种是良性的:
    //
    // - `skipped`:改写出来**不是同一个分子** —— 那是改写器坏了,要吵。
    // - `unshuffled`:分子太对称,所有写法产出同一个串(苯就是) —— 判据
    //   无能为力,如实计数即可。
    let broke_the_molecule = compared == 0 && skipped > 0;
    [
        (
            "写法无关",
            compared > 0,
            broke_the_molecule.then(|| {
                format!(
                    "{skipped} 次改写出来不是同一个分子(另有 {unshuffled} 次换不出新存储序)\
                     —— 改写器坏了,不是画错了"
                )
            }),
        ),
        ("写法无关·比满", compared == writings, None),
        ("写法无关·没查成", compared == 0, None),
    ]
}

/// 环上双键的两条线都要落在某个含它的环里。
fn ring_double_bonds(
    m: &MolBuilder,
    d: &omgkit_depict::Depiction,
    s: &Scene,
    style: &Style,
    clean: bool,
) -> Check {
    if !clean {
        return ("环内双键", false, None);
    }
    let pts = canvas_pts(m, d, style);
    let rings = omgkit_chem::sssr::ring_set(m);
    // **要按图上画成什么来筛,不是按分子里记的是什么** —— 芳香键在图上是交替
    // 单双,而单键根本没有第二条线可查
    let orders = omgkit_depict::render::drawn_orders(m);
    let mut hit = false;
    for (bi, b) in m.bonds().iter().enumerate() {
        let bond_no = u32::try_from(bi).expect("键数超出 u32");
        let mine: Vec<_> = rings
            .iter()
            .filter(|r| r.bonds.contains(&bond_no))
            .collect();
        if mine.is_empty() || orders[bi] != omgkit_core::BondOrder::Double {
            continue;
        }
        let (pa, pb) = (pts[b.begin as usize], pts[b.end as usize]);
        let len = pa.dist(pb);
        if len < 1e-9 {
            continue;
        }
        let mid = (pa + pb) * 0.5;
        let axis = (pb - pa) * (1.0 / len);
        let normal = Point2::new(-axis.y, axis.x);
        // 这根键画出来的线:两端都贴着键轴、落在跨度内
        let lines: Vec<(Point2, Point2)> = s
            .items
            .iter()
            .filter_map(|it| match it {
                Primitive::Line { from, to, .. } => Some((*from, *to)),
                _ => None,
            })
            .filter(|(f, t)| {
                [f, t].iter().all(|p| {
                    let v = **p - pa;
                    v.dot(normal).abs() < 0.30 * len
                        && v.dot(axis) > -0.10 * len
                        && v.dot(axis) < 1.10 * len
                })
            })
            .collect();
        if lines.len() != 2 {
            continue; // 拥挤处会多选到别的线,那时不下判断
        }
        hit = true;
        for (f, t) in lines {
            let lm = (f + t) * 0.5;
            if (lm - mid).dot(normal).abs() < 0.02 * len {
                continue; // 骑在轴上那条
            }
            let ok = mine.iter().any(|r| {
                let poly: Vec<Point2> = r.atoms.iter().map(|a| pts[*a as usize]).collect();
                point_in_polygon(lm, &poly)
            });
            if !ok {
                return (
                    "环内双键",
                    true,
                    Some(format!("键 {bi}({}–{})有一条线画在环外", b.begin, b.end)),
                );
            }
        }
    }
    ("环内双键", hit, None)
}

/// 记了几个楔形,画布上就该有几个。
fn wedges_reach_canvas(d: &omgkit_depict::Depiction, s: &Scene) -> Check {
    let recorded = d.wedges.iter().filter(|w| w.narrow().is_some()).count();
    let drawn = s
        .items
        .iter()
        .filter(|it| matches!(it, Primitive::Wedge { .. } | Primitive::Hash { .. }))
        .count();
    if recorded != drawn {
        return (
            "楔形落地",
            true,
            Some(format!("记了 {recorded} 个,画出来 {drawn} 个")),
        );
    }
    ("楔形落地", recorded > 0, None)
}

/// 画出楔形的中心,反读回来必须还是它该有的构型。
fn wedges_read_back(m: &MolBuilder, d: &omgkit_depict::Depiction) -> Check {
    let genuine = omgkit_io::stereo::genuine_tetrahedral(m);
    let mut hit = false;
    for (i, a) in m.atoms().iter().enumerate() {
        let at = u32::try_from(i).expect("原子数超出 u32");
        if !genuine[i]
            || !matches!(
                a.chiral_tag,
                omgkit_core::ChiralTag::Cw | omgkit_core::ChiralTag::Ccw
            )
            || d.unwedged.contains(&at)
        {
            continue;
        }
        hit = true;
        let got = omgkit_depict::stereo::read_chirality(m, &d.coords, &d.wedges, at);
        if got != Some(a.chiral_tag) {
            return (
                "楔形可读",
                true,
                Some(format!(
                    "中心 {at} 画出来了,反读是 {got:?},该是 {:?}",
                    a.chiral_tag
                )),
            );
        }
    }
    ("楔形可读", hit, None)
}

/// 所有键画出来一样长。
fn bond_lengths_equal(m: &MolBuilder, d: &omgkit_depict::Depiction, clean: bool) -> Check {
    if !clean || m.num_bonds() == 0 {
        return ("键长全等", false, None);
    }
    let first = d.coords[m.bonds()[0].begin as usize].dist(d.coords[m.bonds()[0].end as usize]);
    for (bi, b) in m.bonds().iter().enumerate() {
        let l = d.coords[b.begin as usize].dist(d.coords[b.end as usize]);
        if (l - first).abs() > 1e-6 {
            return (
                "键长全等",
                true,
                Some(format!("键 {bi} 长 {l:.4},第一根长 {first:.4}")),
            );
        }
    }
    ("键长全等", true, None)
}

/// 图元不许伸到画布外。
fn inside_canvas(m: &MolBuilder, d: &omgkit_depict::Depiction, s: &Scene, style: &Style) -> Check {
    // 线、楔形按线宽/宽端算半径;**文字要按真正的字形盒算**。
    //
    // 先前一律拿 `size / 2` 当半径。字号 10pt 给出 ±5.00,而单字符标签的真实
    // 半宽只有 3.34pt(`S`)—— 判据比实现留的白还宽,于是只要有一个单字符标签
    // 落在画布边缘上就误报。实测踩到过:`S1C2=C(...)` 的硫,真盒左边缘在
    // 1.34pt(稳稳在内),判据却按 −0.33 报了违例。
    //
    // **这是判据估粗了,不是实现出界。** 改成用 `label_for` 给的盒,与
    // `render::bounds` 留白用的是同一个来源。
    for it in &s.items {
        let pts: Vec<(Point2, f64, f64)> = match it {
            Primitive::Line { from, to, width } => {
                vec![
                    (*from, *width / 2.0, *width / 2.0),
                    (*to, *width / 2.0, *width / 2.0),
                ]
            }
            Primitive::Wedge { from, to, wide } | Primitive::Hash { from, to, wide, .. } => {
                vec![
                    (*from, *wide / 2.0, *wide / 2.0),
                    (*to, *wide / 2.0, *wide / 2.0),
                ]
            }
            Primitive::Text { .. } => continue, // 文字单独走下面那段
        };
        for (p, rx, ry) in pts {
            if p.x - rx < -0.01
                || p.x + rx > s.width + 0.01
                || p.y - ry < -0.01
                || p.y + ry > s.height
            {
                return (
                    "不出画布",
                    true,
                    Some(format!(
                        "图元在 ({:.2},{:.2})±({rx:.2},{ry:.2}),画布 {:.2}×{:.2}",
                        p.x, p.y, s.width, s.height
                    )),
                );
            }
        }
    }

    let scale = style.bond_length_pt;
    let pts = canvas_pts(m, d, style);
    for a in 0..u32::try_from(m.num_atoms()).expect("原子数超出 u32") {
        let Some(l) = label_at(m, a, style, &d.coords) else {
            continue;
        };
        // **盒心不在原子上**,横排偏 `dx`、竖排还偏 `dy`。口径向实现要:
        // `Label::offset_canvas` 就是 `scene` 摆盒用的那一份。
        let c = pts[a as usize] + l.offset_canvas() * scale;
        let (rx, ry) = (l.half_w * scale, l.half_h * scale);
        if c.x - rx < -0.01 || c.x + rx > s.width + 0.01 || c.y - ry < -0.01 || c.y + ry > s.height
        {
            return (
                "不出画布",
                true,
                Some(format!(
                    "标签 {} 在 ({:.2},{:.2})±({rx:.2},{ry:.2}),画布 {:.2}×{:.2}",
                    l.plain(),
                    c.x,
                    c.y,
                    s.width,
                    s.height
                )),
            );
        }
    }
    ("不出画布", true, None)
}

/// 一根键画出来的每条线,端点都要停在它自己两端标签的字形盒外。
///
/// # 这条判据是怎么来的
///
/// `render::trim` 只管**主线**:它从原子中心出发,按盒边加 margin 切一刀。
/// 双键的第二条线当时不走这条路 —— 它的端点是斜切算出来的角平分线交点,而
/// 角平分线是从**原子中心**量的,端原子有标签时那个交点就落在字里面。
///
/// 实测咖啡因咪唑环上的 C=N(ACS):`N` 的半宽 3.61pt、半高 3.59pt,主线停在离
/// 中心 5.56pt 处,内侧线却停在 3.20pt —— 压在字上,而且比主线还伸得靠前,
/// 看上去像一个指着 N 的楔子。全量 3044 处(17.2%)。
///
/// **实现已经改了**(端原子有标签就不斜切,再加一道 `escape_boxes` 把横向
/// 平移出来的线端点推出盒),这条判据留下来守着它。
///
/// # 只查这根键自己的两端
///
/// 别的键的线压到这个标签上,那是**布局**没摆开,已经报在 `未解冲突`、
/// `有键交叉`、`原子不重合` 里了,不该在这条上再报一遍 —— 实测头几例正是
/// 这样:两个羧基挤到一起,C=O 的线压到了另一个羧基的 O 上。
///
/// # 塞不下的键跳过
///
/// 两端标签加起来比键还长时 `trim` 走压缩兜底,端点落进盒里是明知故犯 ——
/// ACS 规范下 `O⁻—N⁺` 两端要 1.375 个键长的净空,一个键长塞不下。口径由
/// `render::is_squeezed` 给,不在这边重写。
///
/// # 两处判据比实现严,记在这里
///
/// 实现只把 `from` 推出**近端**那个盒、`to` 推出远端那个;判据拿每条线的两个
/// 端点比两端的盒。于是两类构型实现满足不了:(a) 短键上一条线的**远**端点落进
/// **近**端原子的盒;(b) 刚好卡在 `trim` 的 `SQUEEZE` 阈值下方、却触发
/// `escape_boxes` 自己那档压缩的键。扫过真实标签集合 × 两套规范 × 全角度:
/// (b) 几何上可达(5464 种构型),但压缩后端点仍落在这条判据留了线宽余量的盒
/// 之外,报不出来;(a) 语料里不出现。**这条哪天红了,断层多半就在这两处。**
fn lines_clear_of_labels(
    m: &MolBuilder,
    d: &omgkit_depict::Depiction,
    s: &Scene,
    style: &Style,
) -> Check {
    use omgkit_core::BondOrder;
    use omgkit_depict::render::drawn_orders;

    let scale = style.bond_length_pt;
    let pts = canvas_pts(m, d, style);
    let labels = labels_of(m, d, style);
    let tight = squeezed(m, d, style, &labels);
    let orders = drawn_orders(m);

    // `scene` 按键的次序发图元,一根键发几个由画出来的键级定 —— 与那边同源。
    let mut it = s.items.iter();
    // 真比过几对(线端点 × 标签)。**不能拿"过了一个分子"当查到** —— 全量语料
    // 里有 86 个分子×规范一个端点都没得比(没有落在非 tight 键上的带标签端),
    // 那时这条判据是空过的,得如实反映在"查到"那一列里。
    let mut compared = 0usize;
    for (bi, b) in m.bonds().iter().enumerate() {
        let n = match orders[bi] {
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            _ => 1,
        };
        let mine: Vec<&Primitive> = (0..n).filter_map(|_| it.next()).collect();
        if tight[bi] {
            continue;
        }
        for a in [b.begin, b.end] {
            let Some(l) = &labels[a as usize] else {
                continue;
            };
            // **口径向实现要**,见 `render::touches_glyphs`。
            //
            // 这里踩过两次坑:先是按"盒心在原子上"算(整串其实朝一侧挪了 `dx`,
            // 竖排还上下挪了 `dy`),当场多报 4 处假阳;再是拿**整串外接盒**当
            // 字在的地方 —— 上标把盒顶撑高、盒又上下对称,于是符号正下方那一片
            // 纯空白也算进盒里,一条恰好停在符号底下一个 margin 处的线会被报成
            // 压字。字在哪由 `Label::ink` 说了算。
            for p in mine.iter().filter_map(|x| match x {
                Primitive::Line { from, to, .. } => Some([*from, *to]),
                _ => None,
            }) {
                for p in p {
                    compared += 1;
                    if touches_glyphs(l, pts[a as usize], p, scale, style.line_width_pt) {
                        return (
                            "线端不压字",
                            true,
                            Some(format!(
                                "键 {bi} 的线端点 ({:.2},{:.2}) 压在原子 {a} 的标签 {} 的字上",
                                p.x,
                                p.y,
                                l.plain(),
                            )),
                        );
                    }
                }
            }
        }
    }
    // **图元与键的对应必须严丝合缝。** 上面按 `drawn_orders` 数着消耗图元;
    // 哪根键将来多发或少发一个(芳香圈、配位键、跨过交叉点断成两段),整条
    // 判据就静默错位,开始拿这根键的线去比另一根键的标签 —— 而它照样报绿。
    assert!(
        it.all(|x| matches!(x, Primitive::Text { .. })),
        "图元与键对不上号:按键消耗完之后剩下的不全是标签"
    );
    ("线端不压字", compared > 0, None)
}

fn canvas_pts(m: &MolBuilder, d: &omgkit_depict::Depiction, style: &Style) -> Vec<Point2> {
    // **用 render 自己的那两个函数,不许再抄一遍。**
    //
    // 先前这里按"同样的规则"手抄了包围盒的算法。等实现改成"用真正的 `h_side`、
    // 把标签的横向偏移 `dx` 算进去"之后,这份副本没跟着改 —— 判据算出来的原子
    // 位置与 `scene` 画出来的线对不上号,`环内双键` 从 0 违例变成 **1403**,
    // 全是定位错造成的假阳,而画出来的图一点毛病没有。
    let bnd = omgkit_depict::render::bounds(&d.coords, m, style);
    d.coords
        .iter()
        .map(|p| omgkit_depict::render::to_canvas(*p, bnd, style.bond_length_pt))
        .collect()
}
