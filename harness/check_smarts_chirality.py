#!/usr/bin/env python3
"""SMARTS 手性的判据,自带**区分力检查**。

# 空过是这一档的头号风险

查询 `[C@H]...` 与它的镜像 `[C@@H]...` 若在探针集上匹配到同样的东西,这组
探针就分不出两者 —— 拿它验证任何手性改动都会得到"通过",而标记可能早就被
翻了。手性错误**不改拓扑、不改原子数、不改价键**,除了这一档没有别的判据
看得见它。

所以每条查询先过一道区分力检查:

    matches(Q) != matches(镜像 Q)      两侧各自都要满足

不满足的不计入结论,单独报成**语料缺口**。缺口超过上限直接失败 —— 一个几乎
全是空过的判据比没有判据更危险,它会给出令人安心的绿色。

# 查询是从探针**生成**的,不是手写的

手写查询很容易写出根本匹配不上、甚至不合法的东西(例如给手性碳接上五个
邻居),那种查询在两侧都是 0 命中,看起来"一致"。这里改成:对每个探针的每个
立体中心,用 `rootedAtAtom` 生成以该中心为首、或以别的原子为首的 SMILES,
拿它当 SMARTS。这样查询一定匹配得上,而"首原子 / 括号氢 / 环闭合数"这三维
由生成过程自然铺开。

# 三处建这个判据时踩到的坑

它们都会让判据静默变空:

| 坑 | 症状 |
|---|---|
| 用 `"@@" in s` 判断该往哪翻 | 串里后面的中心写成 `@@` 时命中它,翻的是**另一个**原子 |
| 命中只记原子下标 | 探针集里成对放着对映体,两边的下标元组一样,区分力全线归零 |
| 只比命中**数** | 数目相同但命中的是不同分子,那仍是分歧 |

# 跑之前必须重建 wheel

这一档量的是 Rust 侧的解析行为,而它是**通过 wheel** 看到的。改完 Rust 代码
不重新 `maturin build` + `pip install`,量到的就是上一次的产物 —— 结论会稳稳
地指向错误的方向,而且没有任何迹象。

实测踩过:同一份源码前后量出 48/0 与 52/40 两组数,查下来是撤回源码后没重建
wheel,那次"改动前"的基线其实跑的是带补偿的旧构建。

    maturin build --release -m crates/omgkit-py/Cargo.toml --out <目录>
    pip install --force-reinstall <目录>/omgkit-*.whl
    python3 harness/check_smarts_chirality.py

用法:

    python3 harness/check_smarts_chirality.py [--max-vacuous N]
"""

import argparse
import collections
import pathlib
import re
import sys

import omgkit
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

REPO = pathlib.Path(__file__).resolve().parent.parent


def assert_wheel_is_fresh():
    """装好的扩展模块不能比 Rust 源码旧。

    这一档量的是 Rust 侧的行为,而它只能通过 wheel 看到。源码改了却没重建,
    量到的是上一次的产物 —— 结论会稳稳指向错误的方向,而且**毫无迹象**:
    数字照常出来,只是描述的是另一份代码。

    靠文档提醒守不住,得让判据自己发现。比 mtime 是最直接的办法。
    """
    mod = pathlib.Path(omgkit.__file__).resolve()
    # wheel 里模块是 omgkit/__init__.py + 同目录的 .so,取整个目录的最新时间
    built = max(f.stat().st_mtime for f in mod.parent.rglob("*") if f.is_file())
    sources = [f for f in (REPO / "crates").rglob("*.rs")] + [
        f for f in (REPO / "crates").rglob("Cargo.toml")
    ]
    if not sources:
        return
    newest = max(sources, key=lambda f: f.stat().st_mtime)
    if newest.stat().st_mtime > built:
        sys.exit(
            f"装好的 omgkit 比源码旧 —— 量到的会是上一次构建的行为。\n"
            f"  最新源码 {newest.relative_to(REPO)}\n"
            f"  先重建:maturin build --release -m crates/omgkit-py/Cargo.toml --out <目录>\n"
            f"          pip install --force-reinstall <目录>/omgkit-*.whl"
        )


# 探针分子。每一条在关注的位置都是**真**立体中心,并成对给出不同构型 ——
# 只放一种构型的话,查询与它的镜像里总有一个在整个探针集上颗粒无收,
# 那种"零命中对零命中"分不出任何东西。
PROBES = [
    # 无环,三取代 + 氢
    "N[C@H](O)C", "N[C@@H](O)C",
    "N[C@H](O)F", "N[C@@H](O)F",
    "C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O",
    # 无环,四取代(该中心没有氢)
    "N[C@](O)(F)C", "N[C@@](O)(F)C",
    # 五元环:两条环臂不同,中心才算数
    "O[C@H]1CCC[C@@H]1C", "O[C@@H]1CCC[C@@H]1C",
    "C[C@H]1CC[C@@H](O)C1", "C[C@@H]1CC[C@@H](O)C1",
    "C[C@H]1CC[C@H](O)C1", "C[C@@H]1CC[C@H](O)C1",
    # 六元环
    "N[C@H]1CCCC[C@@H]1C", "N[C@@H]1CCCC[C@@H]1C",
    "O[C@H]1CCCC[C@@H]1N", "O[C@@H]1CCCC[C@@H]1N",
    # 环上季碳
    "C[C@]1(N)CC[C@@H](O)C1", "C[C@@]1(N)CC[C@@H](O)C1",
    # 稠双环:融合处的原子带**两个**环闭合
    "C[C@H]1CC[C@H]2CCCC[C@H]12", "C[C@H]1CC[C@@H]2CCCC[C@H]12",
    "O[C@H]1CC[C@H]2CCCC[C@H]12", "O[C@H]1CC[C@@H]2CCCC[C@H]12",
    # 糖环:多中心密集
    "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H]1O",
    "OC[C@@H]1O[C@@H](O)[C@H](O)[C@@H]1O",
]


def mirror(smarts):
    """把 SMARTS 里**第一个**手性标记翻个面。

    不能用 `"@@" in s` 分支:串里后面的中心若写成 `@@`,那个判断会命中它,
    于是翻的是另一个原子 —— 得到的"镜像"根本不是镜像,区分力检查会全线空过。
    """
    i = smarts.index("@")
    if smarts[i : i + 2] == "@@":
        return smarts[:i] + "@" + smarts[i + 2 :]
    return smarts[:i] + "@@" + smarts[i + 1 :]


def first_chiral_shape(smarts):
    """按"第一个手性原子长什么样"归类:(是不是首原子, 有没有括号氢, 几个环闭合)。"""
    i = smarts.index("@")
    lb = smarts.rindex("[", 0, i)
    rb = smarts.index("]", i)
    is_first = lb == 0
    has_h = re.search(r"@+H", smarts[lb : rb + 1]) is not None
    # 紧跟方括号之后的环闭合标号:一位数字,或 `%` 加两位
    tail, n_ring, k = smarts[rb + 1 :], 0, 0
    while k < len(tail):
        if tail[k].isdigit():
            n_ring += 1
            k += 1
        elif tail[k] == "%" and tail[k + 1 : k + 3].isdigit():
            n_ring += 1
            k += 3
        else:
            break
    return is_first, has_h, n_ring


def make_cases():
    """对每个探针的每个立体中心,生成以它为首/不为首的查询。"""
    seen, out = set(), []
    for smi in PROBES:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        centers = [
            i
            for i, _ in Chem.FindMolChiralCenters(
                m, includeUnassigned=False, useLegacyImplementation=False
            )
        ]
        for i in centers:
            for root in {i, (i + 1) % m.GetNumAtoms()}:
                q = Chem.MolToSmiles(m, rootedAtAtom=int(root), isomericSmiles=True)
                if "@" in q and q not in seen:
                    seen.add(q)
                    out.append(q)
    return out


def load():
    om, rd = [], []
    for s in PROBES:
        mo = omgkit.parse_smiles(s)
        mo.sanitize()
        om.append(mo)
        rd.append(Chem.MolFromSmiles(s))
    return om, rd


# 命中必须带上**分子编号**。探针集里成对放着对映体,只记原子下标的话,
# "在分子 A 命中"与"在它的对映体 B 命中"会给出一模一样的元组。


def hits_om(q, probes):
    return sorted(
        (i, tuple(sorted(h))) for i, m in enumerate(probes) for h in q.match(m)
    )


def hits_rd(q, probes):
    return sorted(
        (i, tuple(sorted(t)))
        for i, m in enumerate(probes)
        for t in m.GetSubstructMatches(q, useChirality=True)
    )


def check(smarts, om_probes, rd_probes):
    """返回 (有区分力?, 两侧是否一致, omgkit 命中数, 外部命中数)。"""
    try:
        qo = omgkit.parse_smarts(smarts)
        qo_m = omgkit.parse_smarts(mirror(smarts))
    except ValueError:
        return None
    qr = Chem.MolFromSmarts(smarts)
    qr_m = Chem.MolFromSmarts(mirror(smarts))
    if qr is None or qr_m is None:
        return None
    a_om, b_om = hits_om(qo, om_probes), hits_om(qo_m, om_probes)
    a_rd, b_rd = hits_rd(qr, rd_probes), hits_rd(qr_m, rd_probes)
    return (a_om != b_om) and (a_rd != b_rd), a_om == a_rd, len(a_om), len(a_rd)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--max-vacuous",
        type=int,
        default=0,
        help="最多容忍几条空过。超了说明探针集盖不住生成出来的查询,判据在变空",
    )
    ap.add_argument("--limit", type=int, default=8, help="最多打印几条")
    args = ap.parse_args()

    assert_wheel_is_fresh()
    om_probes, rd_probes = load()
    cases = make_cases()
    print(f"探针 {len(PROBES)} 条,生成查询 {len(cases)} 条\n")

    stat = collections.defaultdict(lambda: {"一致": 0, "反了": 0, "空过": 0})
    vacuous, wrong = [], []
    for smarts in cases:
        r = check(smarts, om_probes, rd_probes)
        if r is None:
            continue
        disc, same, n_om, n_rd = r
        is_first, has_h, n_ring = first_chiral_shape(smarts)
        name = f"{'首' if is_first else '非首'} {'有H' if has_h else '无H'} {n_ring}环"
        if not disc:
            stat[name]["空过"] += 1
            vacuous.append(smarts)
            continue
        stat[name]["一致" if same else "反了"] += 1
        if not same:
            wrong.append((smarts, n_om, n_rd))

    print(f"{'情形':<16}{'一致':>6}{'反了':>6}{'空过':>6}")
    for k in sorted(stat):
        s = stat[k]
        print(f"  {k:<14}{s['一致']:>6}{s['反了']:>6}{s['空过']:>6}")

    print(f"\n有区分力且**反了**的 {len(wrong)} 条:")
    for s, a, b in wrong[: args.limit]:
        print(f"  {s:<46} omgkit {a} / 外部 {b}")
    if len(wrong) > args.limit:
        print(f"  ...(另有 {len(wrong) - args.limit} 条)")

    if vacuous:
        print(f"\n空过的 {len(vacuous)} 条(探针分不出正反,不计入结论):")
        for s in vacuous[: args.limit]:
            print(f"  {s}")
        if len(vacuous) > args.limit:
            print(f"  ...(另有 {len(vacuous) - args.limit} 条)")
        if len(vacuous) > args.max_vacuous:
            print(
                f"\n空过 {len(vacuous)} 条,超过上限 {args.max_vacuous} —— "
                "探针集盖不住生成出来的查询,这条判据正在变空"
            )
            return 1
    return 1 if wrong else 0


if __name__ == "__main__":
    sys.exit(main())
