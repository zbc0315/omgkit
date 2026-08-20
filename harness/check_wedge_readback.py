#!/usr/bin/env python3
"""用外部实现裁判:画出来的楔形,别人读出来是不是同一个构型。

    cargo run -p omgkit-depict --release --example dump_molblock -- \\
        harness/corpus/large.smi > /tmp/blocks.txt
    python3 harness/check_wedge_readback.py /tmp/blocks.txt harness/corpus/large.smi

第二个参数是**同一份语料**:判据要拿它核分母,见下面"分母也要有闸"。

# 为什么判官必须是外部实现

`stereo::assign_wedges` 是"试 Up/Down,取**反读回来**对的那一个"构造出来的,而
反读用的就是 `read_chirality` —— 两者共谋,拿它们的往返去检验是**空过的**,
函数自己的文档就写着这一点。

真正要问的是另一个问题:**别人照着这张图读,读出来是不是同一个分子。** 楔形
读法在一处是**有约定分歧**的:三个画出来的邻居若全挤在中心的同一侧(最大空隙
> 180°,中心落在它们围出的三角形**外面**),"隐式氢在楔形反面"这条读法与四面体
的读法不再等价,不同实现读出**对映体**。这类错误拓扑完全正确、线条毫无毛病,
只有分子是镜像的 —— 自己验自己永远发现不了。

# 按**中心**比,不按分子比

按分子比太粗:一个分子里可能既有如实报了 `unwedged` 的中心、又有没报却画错的
中心,按分子会把后者藏在前者后面。所以逐个原子比 CIP 码。

(下面那张表里,每一档标着它算的是**中心**还是**分子** —— 四档"整条分子进不了
比对"的按分子计,其余按中心计。先前表头统一写"中心数",两种单位混在一列里。)

# 逐原子比就必须先证明两边的下标是同一套

这条判据先前**默认**"molblock 按本实现的原子序写,外部实现按序读",而那句话
是错的:两边的解析器对**显式氢**的取舍不同。

- `MolFromMolBlock` 默认 `removeHs=True`,会把 SMILES 开头写出来的 `[H]` 删掉;
- `MolFromSmiles` 默认也删氢,但**句中**的 `[H]` 与句首的待遇又不一样。

错位是**成对**产生假象的:同一个中心在参考侧落在 k、在读回侧落在 k−1,于是它
同时进"读成了相反的构型"(k 处参考没有码、读回有码)和"外部实现读不出"
(k−1 处参考有码、读回没有码)两档。实测 `large.smi`:

| | 构型一致 | 判官读不出 | 读成了相反的构型 |
|---|---|---|---|
| 错位时(两边都用默认删氢) | 486 | 12 | **10** |
| 对齐后 | **498** | 2 | **0** |

那 10 例一个都不是画错了 —— 它们的 CIP 码两边完全相同,只是下标差一。而"判官
读不出"那一档里 12 个有 10 个也是同一批错位的另一半:那一档本来是拿 RDKit 的
覆盖面立论的,却被判官自己的 bug 灌了五倍的水。

所以现在**两边都关掉删氢**(`removeHs=False` / `SmilesParserParams`),再**逐个
原子核对元素号 + 同位素 + 形式电荷**(见 `aligned`),对不齐直接判违例。

只关 molblock 那一侧是不够的:参考侧留着默认删氢的话,`[H][C@](C)(N)O` 这种
**句中带 `[H]`** 的写法会当场错位,而图一点没画错。`large.smi` 躲过去纯属侥幸
(它带 `[H]` 的行几乎全是 `[H]/N=` 形式,RDKit 因双键顺反留着那个氢),
换一份用 `[H][C@]` 写法的语料就会整片假红。`harness/README.md` 里
`check_write.py` 那一节早就记过同一条("比对时**关掉 `removeHs`**")。

# CIP 码要在**去氢之后**算

RDKit 的传统标注器(`AssignStereochemistry`)在**带显式氢**时会翻标号。实测
`large.smi` 上有 3 个分子如此,例如 `CNC(=O)O[C@@H]1CC[C@H]2[C@@H]1O2`:

- 参考(无显式氢) → 原子 9 = `S`
- 同一个 molblock,带上为画构型补的 2 个氢 → 原子 9 = `R`
- 同一个 molblock,去氢后 → 原子 9 = `S`
- 换成合规的 `rdCIPLabeler`,带不带氢都是 `S`

而**我们导出的 molblock 必然带补出来的氢**(楔形就打在那根 C–H 上),参考侧却
没有 —— 两边喂给标注器的图形状不同,标出来的号自然可能不同。所以两边一律
先去氢再标号,并且改用合规的 `rdCIPLabeler`:一条修图形,一条换掉那个已知会
翻号的部件,两条都不省。下标靠原子属性 `_oi` 穿过去氢带回来。

# 六种结局要分开

- **一致**:图画对了。
- **不一致,但已报 `unwedged`**:图上本来就没画出那个中心,如实说过了。丢信息,
  不好,但下游拿到的是"未指定",不会当真。
  **这一档要求读回侧确实没有码** —— 见下面"`unwedged` 不能当免罪符"。
- **外部实现读不出**:图上画了,而判官给不出 CIP 码。**这不等于画错了** ——
  RDKit 的 2D→手性只认 S(16) 与 Se(34) 这两种三配位中心
  (`Chirality.cpp`:`anum != 16 && anum != 34 && tnzDegree != 4` 就跳过),
  三配位的**磷**不在它的支持范围。这一档单独计数,不当违例,也不藏起来。
- **输入没指定,图上却定死了**:参考侧那个原子没有 CIP 码,我们却画出了构型。
  图在**替分子做主**。**这一档必须是 0。**
- **读成了相反的构型**:图说自己画对了,别人读出来是**对映体**。
  **这一档必须是 0。**
- **该出现却没出现在 dump 里的分子**:见下面"分母也要有闸"。

**把"读不出"和"读反了"混成一档是不行的**:前者是读者的覆盖面,后者是我们画错了,
危害差着量级。

## `unwedged` 不能当免罪符

分支次序是有讲究的:先前只要中心出现在 `#unwedged` 里就一律归"如实报过",
**哪怕图上确实画了楔形、读回来是对映体**。实测(把一个分子整图镜像 = 交付
对映体,同时把它的原子号全塞进 `#unwedged`),那几个中心被静默豁免。

`unwedged` 的语义是"图上没画出这个中心",它该有的后果就是**读回侧没有码**。
读回侧给得出码,就说明我们**画了** —— 那就按画的算,该红照红。

# 分母也要有闸

判据算的是比例式的东西,而它的**分母是别人喂的**:`dump_molblock` 只导出
"有楔形或有 `unwedged`"的分子(源码里那句 `continue`),解析/净化失败的也直接
跳过。任何让 depict **整批看不见立体中心**的回归,都会让这些分子从 dump 里
消失,判据分母缩小、照样全绿。极端情形实测过:**空文件进去,打印"全部通过"
退 0**。

所以现在要第二个参数(语料),拿 RDKit 独立算一遍"哪些分子带四面体标记",
逐个核对 dump 有没有覆盖。实测 `large.smi`:两边都是 311 个分子,**不多不少**。

这与 `feasibility` 那条"该出构型没出"守的是同一件事:几何判据的计数器都在
生成成功之后才累加,不给分母配闸,任何让失败率上升的回归都会让判据变好看。

# 只会让判据变绿的那几档都配了上限,而且档名不许手抄

`unwedged`、"判官读不出"、两档"跳过"都是**单向过滤器** —— 它们只把失配变成
不计数,没有任何力量让判据变红。这种东西必须配上限闸,否则没人拦得住它:
同一个坑在 `harness/verify_stereo.py` 的自校准上栽过。

档名一律用模块级常量,**不在 `main` 里再抄一遍字面量**:`LIMITS` 是按名字查的,
抄错一个字那一档的上限就永久变成"永不触发",而表里还会打印一个 `—`,
看上去像"这一档本来就没上限"。实测过:把档名改一个字再喂整图镜像的 dump,
496 个中心画成对映体,判据打印"全部通过"并退 0。所以末尾还有一道
**闭合校验**:凡是记过数的档,要么在 `LIMITS` 里,要么在 `NO_LIMIT` 里。

双键顺反两边都先清掉:RDKit 会从 2D 坐标读出顺反,而输入 SMILES 里未指定顺反的
双键会因此"凭空多出"立体信息 —— 那是系统性差异,不是画错。

# 实测(`large.smi`,8863 个分子里 311 个带四面体中心,2026-08-20)

    构型一致                        498(中心)
    画了但判官读不出                   2(中心)   ← 同一个分子的两个三配位磷

那 2 个是 `C[P@H]C[P@@H]CS(=O)(=O)[O-]` 的两个磷,与上面 `Chirality.cpp` 那条
限制吻合(逐个原子查过元素号与配位数,不是照着文档推的)。

换 RDKit 2025.09.2 跑同一份 dump:构型一致 495、判官读不出 2、四档硬违例全 0,
照样退 0。少的 3 个是**参考侧**丢了 CIP 码(2023.09 起 `useLegacyStereoPerception`
默认关掉),不是我们画的变了 —— 这条判据两边喂的是同一个 RDKit,版本变化会对消,
所以它断的是上限而不是精确计数。

**`unwedged` 那一档现在是空的**:311 个分子一个 `unwedged` 中心都没有,也就是说
"不一致,但已报 unwedged" 这一档**从未被数据碰到过**。它不是判据在放水,是
depict 眼下把每个中心都画出来了;但空判据要写明白,别当成"验过了"。
"""

import argparse
import collections
import pathlib
import sys

import rdkit
from rdkit import Chem, RDLogger
from rdkit.Chem import rdCIPLabeler

RDLogger.DisableLog("rdApp.*")

#: 去氢之前把原下标记在原子上,去氢之后再取回来。
ORIG_IDX = "_oi"

# **档名当常量用,别在别处再抄字面量。** 见模块文档最后一节:抄错一个字,
# `LIMITS` 就查不到它,那一档的上限永久变成"永不触发"。
OK = "构型一致"
UNWEDGED = "不一致,但已报 unwedged"
UNREADABLE = "画了但判官读不出(见模块文档)"
INVENTED = "**输入没指定,图上却定死了**"
OPPOSITE = "**读成了相反的构型**"
BAD_BLOCK = "**导出的 molblock 读不了**"
MISALIGNED = "**下标对不齐(判据失效)**"
MISSING = "**该出现在 dump 里却没有的分子**"
SKIP_SMILES = "输入外部实现读不了(跳过)"
SKIP_CIP = "指派手性失败(跳过)"

#: 唯一没有上限的一档 —— 它多了只会更好。闭合校验拿它当白名单。
NO_LIMIT = {OK}

#: 各档的上限。**只会让判据变绿的档,一律要有一个**;硬违例的上限是 0。
#:
#: 现值全部是实测值加余量,不是拍的:`large.smi` 上"判官读不出"是 2
#: (一个分子的两个三配位磷),其余全是 0。涨上去要当场查是哪一档退化了,
#: 不是调大这个数。
LIMITS = {
    BAD_BLOCK: 0,
    MISALIGNED: 0,
    MISSING: 0,
    INVENTED: 0,
    OPPOSITE: 0,
    UNREADABLE: 5,
    UNWEDGED: 5,
    SKIP_SMILES: 5,
    SKIP_CIP: 5,
}

#: 每一档算的是中心还是分子。混在一列里会让"跳过 ≤ 5"被读成"最多放过 5 个中心",
#: 而它实际是"最多放过 5 个**分子**的全部中心"(本语料平均 1.6 中心/分子)。
UNIT = {
    OK: "中心",
    UNWEDGED: "中心",
    UNREADABLE: "中心",
    INVENTED: "中心",
    OPPOSITE: "中心",
    BAD_BLOCK: "分子",
    MISALIGNED: "分子",
    MISSING: "分子",
    SKIP_SMILES: "分子",
    SKIP_CIP: "分子",
}

#: 参考侧也要关掉删氢 —— 见模块文档。只关 molblock 那一侧,句中带 `[H]` 的
#: 写法(`[H][C@](C)(N)O`)会当场错位,而图一点没画错。
SMILES_PARAMS = Chem.SmilesParserParams()
SMILES_PARAMS.removeHs = False


def cip_codes(mol):
    """`{原子下标: 'R'/'S'}`,没定的不进表。读不了返回 None。

    下标是**传进来那个分子的**下标 —— 去氢在内部做,靠原子属性 `_oi` 带回来。
    见模块文档"CIP 码要在去氢之后算":两边喂给标注器的图必须是同一个形状,
    否则 RDKit 的传统标注器会因为一边有显式氢而翻号。
    """
    if mol is None:
        return None
    m = Chem.Mol(mol)
    for a in m.GetAtoms():
        a.SetIntProp(ORIG_IDX, a.GetIdx())
    # 只留四面体 —— 双键顺反见模块文档
    for b in m.GetBonds():
        b.SetStereo(Chem.BondStereo.STEREONONE)
        b.SetBondDir(Chem.BondDir.NONE)
    try:
        m = Chem.RemoveHs(m)
        Chem.AssignStereochemistry(m, cleanIt=True, force=True)
        # 合规实现,覆盖上一行留下的传统标号。传统那个在带显式氢时会翻号。
        rdCIPLabeler.AssignCIPLabels(m)
    except Exception:  # noqa: BLE001
        return None
    return {
        a.GetIntProp(ORIG_IDX): a.GetProp("_CIPCode")
        for a in m.GetAtoms()
        if a.HasProp("_CIPCode")
    }


def _fingerprint(mol):
    """逐原子的 (元素号, 同位素, 形式电荷)。"""
    return [(a.GetAtomicNum(), a.GetIsotope(), a.GetFormalCharge()) for a in mol.GetAtoms()]


def aligned(ref, got):
    """两边前 `len(ref)` 个原子是不是同一批(元素号 + 同位素 + 形式电荷)。

    `got` 是从我们导出的 molblock 读回来的,前面是原分子的原子(次序即本实现的
    原子序),末尾可能多出几个**为画构型补的氢**(`with_stereo_hs` 是追加,
    原下标一概不变)。所以判据是:`got` 至少和 `ref` 一样长,而且前 `len(ref)`
    个位置逐个相同。

    这是判据自己的前提,**必须核而不是假设** —— 见模块文档,先前假设它成立,
    结果 10 个分子的下标整体错了一位,判官凭空报出 10 个"画成了对映体"。

    只比元素号是不够的:`C[C@](F)([H])[2H]` 里错开一位之后两边都是氢,元素号
    照样对得上,而 3 号位一边是氘一边是普通氢。所以连同位素和形式电荷一起比。
    """
    fr, fg = _fingerprint(ref), _fingerprint(got)
    return len(fg) >= len(fr) and fg[: len(fr)] == fr


def blocks(text):
    """切成 (行号, SMILES, unwedged 集合, molblock)。"""
    for chunk in text.split("$$$$\n"):
        if not chunk.strip():
            continue
        lines = chunk.splitlines()
        head = next((i for i, l in enumerate(lines) if l.startswith(">>> ")), None)
        if head is None:
            continue
        lineno, smi = lines[head][4:].split("\t", 1)
        unw = set()
        body = head + 1
        if lines[body].startswith("#unwedged"):
            rest = lines[body][len("#unwedged") :].strip()
            unw = {int(x) for x in rest.split(",") if x}
            body += 1
        yield lineno, smi, unw, "\n".join(lines[body:]) + "\n"


def expected_lines(corpus):
    """语料里**该**出现在 dump 中的行号(0 基)。

    口径:RDKit 认为分子里有至少一个原子带四面体标记。实测 `large.smi` 上这个
    集合与 `dump_molblock` 实际导出的 311 个分子**逐个相同**,不多不少 ——
    所以它可以当硬判据用,而不只是个下限。
    """
    want = set()
    for i, line in enumerate(corpus.read_text(encoding="utf-8").splitlines()):
        fields = line.split()
        if not fields:
            continue
        m = Chem.MolFromSmiles(fields[0])
        if m is None:
            continue
        if any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in m.GetAtoms()):
            want.add(i)
    return want


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("dump", type=pathlib.Path, help="dump_molblock 的输出")
    ap.add_argument("corpus", type=pathlib.Path, help="喂给 dump_molblock 的那份语料(核分母用)")
    ap.add_argument("--show", type=int, default=20, help="最多列几例")
    args = ap.parse_args()

    print(f"外部实现:RDKit {rdkit.__version__}")

    tally = collections.Counter()
    bad = []
    seen = set()
    for lineno, smi, unw, block in blocks(args.dump.read_text()):
        seen.add(int(lineno))
        ref = Chem.MolFromSmiles(smi, SMILES_PARAMS)
        # **两边都关掉删氢**,否则下标错位。见模块文档那张表。
        got = Chem.MolFromMolBlock(block, removeHs=False)
        if ref is None:
            tally[SKIP_SMILES] += 1
            continue
        if got is None:
            tally[BAD_BLOCK] += 1
            bad.append((lineno, smi, "molblock 解析失败"))
            continue
        if not aligned(ref, got):
            tally[MISALIGNED] += 1
            bad.append(
                (
                    lineno,
                    smi,
                    f"参考 {ref.GetNumAtoms()} 个原子、读回 {got.GetNumAtoms()} 个,"
                    "元素/同位素/电荷对不上 —— 逐原子比 CIP 的前提不成立",
                )
            )
            continue
        a, b = cip_codes(ref), cip_codes(got)
        if a is None or b is None:
            tally[SKIP_CIP] += 1
            continue
        for k in set(a) | set(b):
            if a.get(k) == b.get(k):
                tally[OK] += 1
            elif b.get(k) is None and k in unw:
                # 图上确实没画出这个中心,而且如实报过了
                tally[UNWEDGED] += 1
            elif b.get(k) is None:
                # 画了,但判官给不出 CIP 码 —— 是它的覆盖面,不是我们画错了
                tally[UNREADABLE] += 1
            elif a.get(k) is None:
                # 反过来:输入没指定那个中心,图却把它定死了 —— 图在替分子做主
                tally[INVENTED] += 1
                bad.append(
                    (lineno, smi, f"原子 {k}:输入没指定构型,本实现画成了 {b.get(k)}")
                )
            else:
                # 读回侧给得出码 —— 说明我们**画了**。哪怕它进了 `unwedged`,
                # 也按画的算:`unwedged` 的语义是"没画",不是免罪符。
                tally[OPPOSITE] += 1
                why = f"原子 {k}:本实现画成 {b.get(k)},应当是 {a.get(k)}"
                if k in unw:
                    why += "(而且它还自报了 unwedged —— 图上明明画了)"
                bad.append((lineno, smi, why))

    # **分母也要有闸。** 见模块文档:dump 少喂几个分子,上面每一档都会变好看。
    missing = sorted(expected_lines(args.corpus) - seen)
    tally[MISSING] += len(missing)
    for i in missing[:20]:
        bad.append((i, "(不在 dump 里)", "语料里这个分子带四面体标记,dump 却没有它"))

    print(f"{'档':<32}{'数量':>8}{'单位':>6}{'上限':>8}")
    for k, v in tally.most_common():
        cap = LIMITS.get(k)
        print(f"{k:<32}{v:>8}{UNIT.get(k, '?'):>6}{'—' if cap is None else cap:>8}")
    if bad:
        print(f"\n=== 前 {args.show} 例 ===")
        for lineno, smi, why in bad[: args.show]:
            print(f"  行 {lineno}: {smi}\n      {why}")
        if len(bad) > args.show:
            print(f"  …… 另有 {len(bad) - args.show} 例")

    print()
    # **闭合校验。** 记过数的档必须有归属,否则改一个档名就能悄悄关掉一道闸。
    orphan = sorted(set(tally) - set(LIMITS) - NO_LIMIT)
    if orphan:
        print(f"有档既不在 LIMITS 也不在 NO_LIMIT 里:{orphan}")
        print("  —— 档名是按字符串查上限的,漏一个等于那一档永不触发。")
        return 1

    over = [(k, tally[k], cap) for k, cap in LIMITS.items() if tally[k] > cap]
    if not over:
        print("全部通过。")
        return 0
    for k, v, cap in over:
        print(f"{k}:{v} > 上限 {cap}")
    print(
        "\n只会让判据变绿的那几档超了上限的话,别调大它 —— 那几档是单向过滤器,"
        "涨上去说明有一类中心正在被悄悄豁免,要当场查是哪一类。"
    )
    return 1


if __name__ == "__main__":
    sys.exit(main())
