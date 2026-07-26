#!/usr/bin/env python3
"""模板在产物侧写手性时,产物构建有没有把参照系换错。带**区分力检查**。

# 这一档在仓库语料上的触发面是 0

`harness/corpus/reactions.txt` 的 22 条模板里,没有一条在产物侧写手性。也就是
说,反应差分跑得再多也照不到这段代码 —— 它一直是个盲区。真实模板里这却是
常态:立体专一的反应正是靠产物侧的 `@` 表达构型的。

所以这一档的用例写死在文件里,不从语料取。

# 模板两侧写没写手性,是四种不同的指令

| 反应物侧 | 产物侧 | 含义 |
|---|---|---|
| 没写 | 没写 | 模板没管 —— 底物的构型带过来 |
| 写了 | 没写 | 构型被破坏 |
| 没写 | 写了 | 构型是新建的,与底物无关 |
| 写了 | 写了 | 相对底物**保留**(两标记相同)或**翻转**(不同) |

最后一行最容易做错。照字面把产物侧那个标记写死的话,同一个模板作用在一对
对映体上会给出**同一个**产物 —— 而正确答案是一对对映体。这个错误不改拓扑、
不改原子数,除了逐对映体比对没有别的判据看得见。

用例因此成对给出对映体底物,并检查:

- `Set` 一档:两个底物必须给出**同一个**产物
- `Retain`/`Invert` 一档:两个底物必须给出**不同**的产物

不满足就是空过 —— 那组用例分不出正反,拿它验证任何改动都会得到"通过"。

# 跑之前必须重建 wheel

这一档量的是 Rust 侧的行为,而它是**通过 wheel** 看到的。改完 Rust 代码不重建,
量到的是上一次的产物,而且毫无迹象。判据自己会检查,见 `assert_wheel_is_fresh`。

    maturin build --release -m crates/omgkit-py/Cargo.toml --out <目录>
    pip install --force-reinstall <目录>/omgkit-*.whl
    python3 harness/check_product_chirality.py
"""

import argparse
import collections
import pathlib
import sys

import omgkit
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

REPO = pathlib.Path(__file__).resolve().parent.parent

# 一对对映体。四种指令都拿它们跑,靠"产物随不随底物变"分辨。
PAIR = ("C[C@H](O)CC", "C[C@@H](O)CC")

# (名字, 模板, 产物该不该随底物变)
INSTRUCTIONS = [
    ("两侧都没写 —— 带过来", "[C:1]-[OH:2]>>[C:1]-[Cl:2]", True),
    ("只有产物侧写 —— 新建", "[C:1]-[OH:2]>>[C@:1]-[Cl:2]", False),
    ("只有产物侧写 + 括号氢", "[CH:1]-[OH:2]>>[C@H:1]-[Cl:2]", False),
    ("两侧同标记 —— 保留", "[C@:1]-[OH:2]>>[C@:1]-[Cl:2]", True),
    ("两侧同标记 + 括号氢", "[C@H:1]-[OH:2]>>[C@H:1]-[Cl:2]", True),
    ("两侧异标记 —— 翻转", "[C@:1]-[OH:2]>>[C@@:1]-[Cl:2]", True),
    ("两侧异标记 + 括号氢", "[C@H:1]-[OH:2]>>[C@@H:1]-[Cl:2]", True),
    # 标记是相对**这张模板自己**的邻居顺序的。产物模板把取代基对调着写时,
    # 同一个 `@` 说的已经是另一种构型 —— 只比标记就会把翻转读成保留。
    # 三条要一起看:次序对调那条给出的产物,应当与次序不变那条正好相反。
    (
        "次序不变 同标记 —— 保留",
        "[CH3:3]-[C@H:1](-[OH:2])-[CH2:4]>>[CH3:3]-[C@H:1](-[Cl:2])-[CH2:4]",
        True,
    ),
    (
        "次序对调 同标记 —— 其实是翻转",
        "[CH3:3]-[C@H:1](-[OH:2])-[CH2:4]>>[CH2:4]-[C@H:1](-[Cl:2])-[CH3:3]",
        True,
    ),
    (
        "次序对调 异标记 —— 各反一次即保留",
        "[CH3:3]-[C@H:1](-[OH:2])-[CH2:4]>>[CH2:4]-[C@@H:1](-[Cl:2])-[CH3:3]",
        True,
    ),
]

# 真实模板的形态。手性原子在产物模板里长什么样,决定要补哪几项参照系换算:
# 是不是片段首原子、有没有括号氢、写了几根键、带不带环闭合。
SHAPES = [
    ("首原子 有H 一根模板键",
     "[CH;D3;+0:1]-[C;H0;D3;+0:2]=[O;H0;D1;+0:3]"
     ">>[C@@H;D3;+0:1]-[C@H;D3;+0:2]-[OH;D1;+0:3]",
     "O=C1CCCC1CCC1=CC=CC2=C1N=CC=C2"),
    ("非首 有H 一根模板键 + 断成两片",
     "[C@@H;D3;+0:1]-[n;H0;D3;+0:2]"
     ">>Br-[C@@H;D3;+0:1].C-[Si](-C)(-C)-[n;H0;D3;+0:2]",
     "C[C@H](n1cccc1)CC"),
    ("非首 有H 两根模板键 + 新建两个中心",
     "[C;H0;D3;+0:2]=[CH;D2;+0:1]"
     ">>C-S(=O)(=O)-O-[C@H;D3;+0:1]-[C@@H;D3;+0:2]",
     "COC(=O)C1=C[C@H](O)[C@@H](O)[C@@H](O)O1"),
    ("非首 无H 三根模板键 + 环闭合",
     "[CH;D2;+0:3]=[C;H0;D3;+0:2]-[C@H;D3;+0:7]-[CH2;D2;+0:4]"
     "-[C;H0;D3;+0:5](=[O;H0;D1;+0:6])-[O;H0;D2;+0:1]"
     ">>[O;H0;D2;+0:1]-[C@;H0;D4;+0:2]1-[C@@H;D3;+0:3]-[CH2;D2;+0:4]"
     "-[C@@H;D3;+0:5](-[OH;D1;+0:6])-[C@H;D3;+0:7]-1",
     "[H][C@@]12C[C@H](CC=C1COC(=O)C2)C(C)(C)C"),
    ("四取代中心全在模板里",
     "[C:1](-[N:2])(-[O:3])-[F:4]>>[C@@:1](-[N:2])(-[O:3])-[F:4]",
     "N[C@H](O)F"),
]


def assert_wheel_is_fresh():
    """装好的扩展模块不能比 Rust 源码旧,理由同 `check_smarts_chirality.py`。"""
    mod = pathlib.Path(omgkit.__file__).resolve()
    built = max(f.stat().st_mtime for f in mod.parent.rglob("*") if f.is_file())
    sources = list((REPO / "crates").rglob("*.rs")) + list(
        (REPO / "crates").rglob("Cargo.toml")
    )
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


def canonical(smi):
    """两侧的产物都交给外部实现规范化 —— 各自规范化再比,比出来的全是噪声。"""
    m = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(m) if m is not None else "<invalid>"


def run_om(tmpl, smi, limit=200):
    rxn = omgkit.parse_reaction(tmpl)
    mol = omgkit.parse_smiles(smi)
    mol.remove_hs()
    mol.sanitize()
    out = []
    for outcome in rxn.run([mol], max_products=limit):
        group = []
        for p in outcome.products:
            try:
                p.sanitize()
            except ValueError:
                group.append("<invalid>")
                continue
            group.append(canonical(p.to_smiles()))
        out.append(tuple(sorted(group)))
    return collections.Counter(out)


def run_rd(tmpl, smi, limit=200):
    rxn = AllChem.ReactionFromSmarts(tmpl)
    mol = Chem.MolFromSmiles(smi)
    out = []
    for i, group in enumerate(rxn.RunReactants((mol,))):
        if i >= limit:
            break
        row = []
        for p in group:
            try:
                Chem.SanitizeMol(p)
            except Exception:  # noqa: BLE001
                row.append("<invalid>")
                continue
            row.append(Chem.MolToSmiles(p))
        out.append(tuple(sorted(row)))
    return collections.Counter(out)


def mirror_products(tmpl):
    """把模板**产物侧**第一个手性标记翻面,用来做形态那一组的区分力检查。"""
    lhs, sep, rhs = tmpl.partition(">>")
    i = rhs.index("@")
    if rhs[i : i + 2] == "@@":
        rhs = rhs[:i] + "@" + rhs[i + 2 :]
    else:
        rhs = rhs[:i] + "@@" + rhs[i + 1 :]
    return lhs + sep + rhs


def check_instructions(limit):
    """四种指令各自的语义:产物该不该随底物变。"""
    bad = vacuous = 0
    print("四种指令")
    for name, tmpl, varies in INSTRUCTIONS:
        a, b = (run_om(tmpl, s) for s in PAIR)
        ra, rb = (run_rd(tmpl, s) for s in PAIR)
        # 区分力:外部实现自己要先分得出这两个底物,否则这条用例说明不了什么
        if varies and ra == rb:
            print(f"  空过  {name} —— 外部实现在这对对映体上也不变")
            vacuous += 1
            continue
        got_varies = a != b
        if got_varies != varies:
            want = "随底物变" if varies else "与底物无关"
            print(f"  错    {name}:应当{want},实际不是")
            bad += 1
            continue
        if (a, b) != (ra, rb):
            print(f"  不同  {name}")
            for label, x, y in (("R", a, ra), ("S", b, rb)):
                if x != y:
                    print(f"          {label} 底物 omgkit {sorted(x)} / 外部 {sorted(y)}")
            bad += 1
        else:
            print(f"  一致  {name}")
    return bad, vacuous


def check_shapes(limit):
    """真实模板的形态:手性原子在产物模板里怎么写,参照系就差几项。"""
    bad = vacuous = 0
    print("\n产物模板里手性原子的形态")
    for name, tmpl, smi in SHAPES:
        a, r = run_om(tmpl, smi), run_rd(tmpl, smi)
        mt = mirror_products(tmpl)
        am, rm = run_om(mt, smi), run_rd(mt, smi)
        if a == am or r == rm:
            print(f"  空过  {name} —— 翻掉模板里的 `@` 后产物不变")
            vacuous += 1
            continue
        if a == r:
            print(f"  一致  {name}")
            continue
        print(f"  不同  {name}")
        for g, n in list((a - r).items())[:limit]:
            print(f"          omgkit 独有 x{n}: {' . '.join(g)}")
        for g, n in list((r - a).items())[:limit]:
            print(f"          外部   独有 x{n}: {' . '.join(g)}")
        bad += 1
    return bad, vacuous


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--limit", type=int, default=2, help="每条最多打印几组产物")
    ap.add_argument(
        "--max-bad",
        type=int,
        default=0,
        help="容忍几条不一致。形态那一组还有已知缺口时用得上",
    )
    args = ap.parse_args()

    assert_wheel_is_fresh()
    b1, v1 = check_instructions(args.limit)
    b2, v2 = check_shapes(args.limit)
    bad, vacuous = b1 + b2, v1 + v2
    print(f"\n不一致 {bad} 条,空过 {vacuous} 条(上限 {args.max_bad})")
    if vacuous:
        print("空过说明用例分不出正反 —— 那样的判据比没有判据更危险")
        return 1
    return 1 if bad > args.max_bad else 0


if __name__ == "__main__":
    sys.exit(main())
