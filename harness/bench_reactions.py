#!/usr/bin/env python3
"""在真实反应语料上比对 omgkit 与外部实现的产物生成:正确性 + 耗时。

输入是三列 CSV `反应id,反应SMILES,逆向模板`,其中反应 SMILES 形如
`反应物.反应物>>产物`,模板是 `产物模式>>反应物模式` 的逆合成写法。

# 两个方向都有基准答案

逆向模板作用在**产物**上应当给出反应物;把模板两侧对调就成了正向模板,
作用在**反应物**上应当给回产物。真实语料自带这两侧的正确答案,所以不必
靠人工构造用例。

# 产物比对不能各自规范化

两边的规范化是两套不同的算法,同一个分子写出来的串本来就不一样。所以两边
的产物都交给**同一个读者**读进来再规范化,比的是多重集。这条与
`check_reactions.py` 同一个道理。

# 耗时要分段量

"跑一条反应"包含三件事:编译模板、准备底物、真正做图改写。三者的比例差得很远,
混在一起量出来的数字说明不了任何问题。所以分开计时,并且两边都**预编译**模板、
**预准备**底物,只让被比较的那一段进入计时。

用法:

    python3 harness/bench_reactions.py <语料.csv> [--limit N] [--direction both]
"""

import argparse
import collections
import csv
import pathlib
import sys
import time

try:
    import omgkit
except ImportError:
    sys.exit("导入不到 omgkit。先构建并安装 wheel,见 harness/test_python.py 的说明。")

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")


# ---------------------------------------------------------------------------
# 归一:两侧必须同口径
# ---------------------------------------------------------------------------
#
# 产物生成出来的图**都是未净化的**(两边都一样,模板本就能写出不合法的东西)。
# 直接把未净化的图写成 SMILES 再比很不公平:两边未净化输出的可读性不同,
# 一侧写得出能读回去的串、另一侧写不出,比出来的"分歧"其实是净化时机的差别。
#
# 所以两侧都**先净化再写出**,净化不了的一律记成 `<unsanitizable>` —— 这样
# "净化不了"本身成为可比的一类结果,而不是伪装成结构差异。


def norm_omgkit(mol):
    """omgkit 的一个产物 → 规范 SMILES。"""
    try:
        mol.sanitize()
    except ValueError:
        return "<unsanitizable>"
    m = Chem.MolFromSmiles(mol.to_smiles())
    return Chem.MolToSmiles(m) if m is not None else "<unreadable>"


def norm_rdkit(mol):
    """外部实现的一个产物 → 规范 SMILES,口径同上。"""
    try:
        Chem.SanitizeMol(mol)
    except Exception:  # noqa: BLE001
        return "<unsanitizable>"
    try:
        return Chem.MolToSmiles(mol)
    except Exception:  # noqa: BLE001
        return "<unreadable>"


# 生成与归一必须**分开**:归一里含净化与规范化,把它算进计时,量到的就不再是
# "产物生成"这一段了。所以下面两个函数只做生成,归一在计时之外单独跑。


def omgkit_produce(rxn, mols, cap):
    return [o.products for o in rxn.run(mols, max_products=cap)]


def rdkit_produce(rxn, mols, cap):
    out = []
    for i, prods in enumerate(rxn.RunReactants(tuple(mols))):
        if cap and i >= cap:
            break
        out.append(prods)
    return out


def omgkit_sets(groups):
    return collections.Counter(
        tuple(sorted(norm_omgkit(p) for p in g)) for g in groups
    )


def rdkit_sets(groups):
    return collections.Counter(tuple(sorted(norm_rdkit(p) for p in g)) for g in groups)


def norm_set(smiles_list):
    """一组 SMILES 串 → 规范形式元组。用来把语料里的真实答案归到同一口径。"""
    out = []
    for s in smiles_list:
        m = Chem.MolFromSmiles(s) if s else None
        out.append(Chem.MolToSmiles(m) if m is not None else "<unreadable>")
    return tuple(sorted(out))


# ---------------------------------------------------------------------------
# 一趟测量
# ---------------------------------------------------------------------------


def reverse_template(t):
    """把 `A>>B` 对调成 `B>>A`。"""
    lhs, _, rhs = t.partition(">>")
    return f"{rhs}>>{lhs}"


def run_direction(rows, direction, cap, verbose, merge_hs):
    """跑一个方向,返回 (统计, 计时)。

    `direction` 取 "retro"(模板作用于产物)或 "forward"(对调模板,作用于反应物)。
    """
    stat = collections.Counter()
    timing = collections.defaultdict(float)
    examples = []

    for rid, rxn_smiles, template in rows:
        lhs, _, product = rxn_smiles.partition(">>")
        reactants = lhs.split(".")

        if direction == "retro":
            tmpl, substrates = template, [product]
        else:
            tmpl, substrates = reverse_template(template), reactants

        # --- 编译模板(两边各自计时)---
        t0 = time.perf_counter()
        try:
            om_rxn = omgkit.parse_reaction(tmpl)
        except ValueError:
            om_rxn = None
        timing["om_compile"] += time.perf_counter() - t0

        t0 = time.perf_counter()
        try:
            rd_rxn = AllChem.ReactionFromSmarts(tmpl)
        except Exception:
            rd_rxn = None
        timing["rd_compile"] += time.perf_counter() - t0

        if om_rxn is None or rd_rxn is None:
            stat["模板两边至少一侧编译不了" if (om_rxn is None) != (rd_rxn is None)
                 else "模板两边都编译不了"] += 1
            continue

        # 反应物模板个数必须与底物个数对得上,否则这条用例本身没法比
        if om_rxn.num_reactant_templates != len(substrates):
            stat["底物数与模板对不上(跳过)"] += 1
            continue

        # --- 准备底物 ---
        t0 = time.perf_counter()
        om_mols = []
        try:
            for s in substrates:
                m = omgkit.parse_smiles(s)
                # 显式氢会把邻接原子的度数撑大,写着 D3 的模板就配不上它。
                # 外部实现在解析时就并掉了,这里显式做同一件事,两边的图才可比。
                if merge_hs:
                    m.remove_hs()
                m.sanitize()
                om_mols.append(m)
        except ValueError:
            om_mols = None
        timing["om_prep"] += time.perf_counter() - t0

        t0 = time.perf_counter()
        rd_mols = [Chem.MolFromSmiles(s) for s in substrates]
        if any(m is None for m in rd_mols):
            rd_mols = None
        timing["rd_prep"] += time.perf_counter() - t0

        if om_mols is None or rd_mols is None:
            stat["底物至少一侧读不了(跳过)"] += 1
            continue

        # --- 真正的产物生成(只有这一段进计时)---
        t0 = time.perf_counter()
        try:
            om_raw = omgkit_produce(om_rxn, om_mols, cap)
        except Exception as e:  # noqa: BLE001
            stat[f"omgkit 抛异常: {type(e).__name__}"] += 1
            continue
        timing["om_run"] += time.perf_counter() - t0

        t0 = time.perf_counter()
        try:
            rd_raw = rdkit_produce(rd_rxn, rd_mols, cap)
        except Exception as e:  # noqa: BLE001
            stat[f"外部实现抛异常: {type(e).__name__}"] += 1
            continue
        timing["rd_run"] += time.perf_counter() - t0

        # --- 归一在计时之外 ---
        om_out = omgkit_sets(om_raw)
        rd_out = rdkit_sets(rd_raw)

        # 底物带显式 [H] 的那一档单独统计。
        #
        # omgkit 把 `removeHs` 划在净化之外,显式氢原子留在图里;别的实现通常在
        # 解析时就把它们并掉。于是同一条 SMILES 在两边得到的**图本来就不同** ——
        # 显式氢会把邻接的原子从 D3 撑成 D4,而模板里写着 `D3`,匹配自然不同。
        # 这是已知的架构差异,不该混进总体一致率里,否则那个数说的是两件事。
        has_h = any("[H]" in s for s in substrates)
        bucket = "含显式H" if has_h else "无显式H"
        stat[f"比对过的反应({bucket})"] += 1
        stat["比对过的反应"] += 1
        if om_out == rd_out:
            stat["产物完全一致"] += 1
            stat[f"产物完全一致({bucket})"] += 1
            if om_out:
                stat["其中确实有产物"] += 1
        else:
            stat["产物不同"] += 1
            stat[f"产物不同({bucket})"] += 1
            if not has_h and len(examples) < 6:
                examples.append((rid, tmpl, substrates, om_out, rd_out))

        # 方向性的正确性:结果里应当出现真实的另一侧
        truth = norm_set(reactants if direction == "retro" else [product])
        if truth in om_out:
            stat["omgkit 命中真实答案"] += 1
            stat[f"omgkit 命中真实({bucket})"] += 1
        if truth in rd_out:
            stat["外部实现命中真实答案"] += 1
            stat[f"外部实现命中真实({bucket})"] += 1

    if verbose:
        for rid, tmpl, subs, a, b in examples:
            print(f"\n  反应 {rid}")
            print(f"    模板 {tmpl[:110]}")
            print(f"    底物 {'.'.join(subs)[:110]}")
            print(f"    omgkit  {list(a)[:2]}")
            print(f"    外部实现 {list(b)[:2]}")

    return stat, timing


def report(name, stat, timing, n):
    print(f"\n===== {name} =====")
    for k, v in stat.most_common():
        print(f"  {k:<34} {v}")
    compared = stat["比对过的反应"]
    if compared:
        print(f"  {'一致率(全部)':<34} {stat['产物完全一致'] / compared * 100:.2f}%")
    for bucket in ("无显式H", "含显式H"):
        # 不要用 `n` 当循环变量 —— 它是本函数的参数,遮蔽之后下面的耗时表会印错条数
        cnt = stat[f"比对过的反应({bucket})"]
        if cnt:
            pct = stat[f"产物完全一致({bucket})"] / cnt * 100
            hit_om = stat[f"omgkit 命中真实({bucket})"]
            hit_rd = stat[f"外部实现命中真实({bucket})"]
            print(
                f"  {'一致率(' + bucket + ')':<34} {pct:.2f}%   "
                f"({cnt} 条;命中真实答案 omgkit {hit_om} / 外部 {hit_rd})"
            )
    print(f"\n  耗时(共 {n} 条,单位毫秒,越小越好)")
    print(f"    {'':<12}{'omgkit':>12}{'外部实现':>14}{'比值':>10}")
    for stage, label in [("compile", "编译模板"), ("prep", "准备底物"), ("run", "产物生成")]:
        om = timing[f"om_{stage}"] * 1000
        rd = timing[f"rd_{stage}"] * 1000
        ratio = f"{rd / om:.2f}×" if om > 0 else "—"
        print(f"    {label:<12}{om:>12.1f}{rd:>14.1f}{ratio:>10}")
    om_tot = sum(timing[f"om_{s}"] for s in ("compile", "prep", "run")) * 1000
    rd_tot = sum(timing[f"rd_{s}"] for s in ("compile", "prep", "run")) * 1000
    ratio = f"{rd_tot / om_tot:.2f}×" if om_tot > 0 else "—"
    print(f"    {'合计':<12}{om_tot:>12.1f}{rd_tot:>14.1f}{ratio:>10}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("csv", type=pathlib.Path, help="三列 CSV:id,反应SMILES,逆向模板")
    ap.add_argument("--limit", type=int, default=2000, help="只跑前几条")
    ap.add_argument("--cap", type=int, default=50, help="单条反应最多取几组产物")
    ap.add_argument(
        "--direction",
        choices=["retro", "forward", "both"],
        default="both",
        help="retro=模板作用于产物;forward=对调模板作用于反应物",
    )
    ap.add_argument("--verbose", action="store_true", help="打印不一致的例子")
    ap.add_argument("--merge-hs", action="store_true",
                    help="omgkit 侧先合并显式氢,让两边的图可比")
    args = ap.parse_args()

    rows = []
    with args.csv.open(newline="") as f:
        for rec in csv.reader(f):
            if len(rec) == 3 and ">>" in rec[1] and ">>" in rec[2]:
                rows.append(rec)
            if len(rows) >= args.limit:
                break
    if not rows:
        sys.exit("语料一条都没读到")
    print(f"语料 {len(rows)} 条 | omgkit {omgkit.__version__} | 外部实现 {Chem.rdBase.rdkitVersion}")

    todo = ["retro", "forward"] if args.direction == "both" else [args.direction]
    for d in todo:
        stat, timing = run_direction(rows, d, args.cap, args.verbose, args.merge_hs)
        report("逆向(模板作用于产物)" if d == "retro" else "正向(对调模板,作用于反应物)",
               stat, timing, len(rows))
    return 0


if __name__ == "__main__":
    sys.exit(main())
