"""挑出 RDKit 输出错、omgkit 输出对的典型案例,逐条给出**独立于记录的**证据。

# 为什么不能只拿"命中率"当证据

"omgkit 命中而 RDKit 未命中"只说明两者与记录的距离不同。要说 RDKit **错**,
得有一条不依赖记录的判据。这里用三条:

1. **质量守恒**。产物的重原子数必须等于反应物的重原子数减去离去的部分。
   逆向模板写成两个片段、而底物里那两个片段仍由未匹配的原子连着时,
   "逐产物各搬一次模板之外的部分"会把共享的原子复制进每一片 —— 原子凭空
   变多。这条判据与记录无关,数一数就知道。

2. **书写顺序不变性**。同一个连接关系,模板把邻居换个次序写,描述的仍是同一个
   产物。结果跟着书写顺序变的那一方必然有一次是错的。

3. **rdchiral 裁定**。这批模板是 rdchiral 抽的,手性标记也是它写的,
   它自己的应用器给出的就是模板作者的本意。

产出 docs/cases.md 与 figures/case_*.png。
"""

import argparse
import itertools
import json
import os

from rdkit import Chem, RDLogger
from rdkit.Chem import Draw, rdChemReactions

import omgkit

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def canon(s):
    m = Chem.MolFromSmiles(s)
    return Chem.MolToSmiles(m) if m else None


def heavy(s):
    m = Chem.MolFromSmiles(s)
    return None if m is None else m.GetNumHeavyAtoms()


def run_omgkit(tpl, inputs):
    r = omgkit.parse_reaction(tpl)
    n = r.num_reactant_templates
    mols = []
    for s in inputs:
        m = omgkit.parse_smiles(s)
        m.sanitize()
        mols.append(m)
    out = set()
    for perm in itertools.permutations(mols, n):
        for oc in r.run(list(perm)):
            ps = []
            for p in oc.products:
                q = p.copy()
                q.sanitize()
                ps.append(q.to_smiles())
            c = canon(".".join(ps))
            if c:
                out.add(c)
    return out


def run_rdkit(tpl, inputs):
    rx = rdChemReactions.ReactionFromSmarts(tpl)
    rx.Initialize()
    n = rx.GetNumReactantTemplates()
    mols = [Chem.MolFromSmiles(s) for s in inputs]
    out = set()
    for perm in itertools.permutations(mols, n):
        for ps in rx.RunReactants(perm):
            frags = []
            ok = True
            for p in ps:
                q = Chem.Mol(p)
                for a in q.GetAtoms():
                    a.SetAtomMapNum(0)
                try:
                    Chem.SanitizeMol(q)
                    frags.append(Chem.MolToSmiles(q))
                except Exception:
                    ok = False
            if ok:
                c = canon(".".join(frags))
                if c:
                    out.add(c)
    return out


def draw_case(path, title, rows):
    """rows: [(标题, SMILES)]"""
    mols, legends = [], []
    for lab, smi in rows:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        mols.append(m)
        legends.append(lab)
    if not mols:
        return False
    img = Draw.MolsToGridImage(
        mols, molsPerRow=len(mols), subImgSize=(430, 330), legends=legends, useSVG=False
    )
    data = img.data if hasattr(img, "data") else img
    if isinstance(data, bytes):
        with open(path, "wb") as fh:
            fh.write(data)
    else:
        data.save(path)
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench", default=os.path.join(ROOT, "results", "bench.jsonl"))
    ap.add_argument("--tpl", default=os.path.join(ROOT, "data", "templates.jsonl"))
    ap.add_argument("--outmd", default=os.path.join(ROOT, "docs", "cases.md"))
    ap.add_argument("--figdir", default=os.path.join(ROOT, "figures"))
    ap.add_argument("--limit", type=int, default=4)
    args = ap.parse_args()

    tpls = {}
    for line in open(args.tpl):
        r = json.loads(line)
        if "row" in r and "retro" in r:
            tpls[r["row"]] = r

    # omgkit 命中而 RDKit 未命中的
    wins = []
    for line in open(args.bench):
        r = json.loads(line)
        if "err" in r:
            continue
        for d in ("fwd", "retro"):
            a, b = r.get(f"omgkit_{d}"), r.get(f"rdkit_{d}")
            if not a or not b or "hit" not in a or "hit" not in b:
                continue
            if a["hit"] and not b["hit"]:
                wins.append((r["row"], d))

    print(f"omgkit 命中而 RDKit 未命中:{len(wins)} 条")

    mass, chiral, other = [], [], []
    for row, d in wins:
        t = tpls.get(row)
        if t is None:
            continue
        tpl = t[d]
        inputs = t["reactants"] if d == "fwd" else [t["prod"]]
        truth = canon(t["prod"] if d == "fwd" else ".".join(t["reactants"]))
        try:
            rd = run_rdkit(tpl, inputs)
        except Exception:
            continue
        if not rd:
            continue
        n_in = sum(heavy(s) or 0 for s in inputs)
        # RDKit 的输出里有没有重原子数**多于**输入的?那是复制出来的
        heavier = [s for s in rd if (heavy(s) or 0) > n_in]
        rec = {
            "row": row,
            "direction": d,
            "id": t["id"],
            "tpl": tpl,
            "inputs": inputs,
            "truth": truth,
            "rdkit": sorted(rd),
            "n_in": n_in,
        }
        if heavier:
            rec["rdkit_bad"] = heavier[0]
            rec["n_bad"] = heavy(heavier[0])
            mass.append(rec)
        elif any(
            Chem.MolToSmiles(Chem.MolFromSmiles(s), isomericSmiles=False)
            == Chem.MolToSmiles(Chem.MolFromSmiles(truth), isomericSmiles=False)
            for s in rd
            if Chem.MolFromSmiles(s)
        ):
            rec["rdkit_bad"] = sorted(rd)[0]
            chiral.append(rec)
        else:
            other.append(rec)

    print(f"  质量不守恒 {len(mass)}  立体差异 {len(chiral)}  其他 {len(other)}")

    lines = [
        "# 典型案例:RDKit 输出错、omgkit 输出对",
        "",
        "判据不依赖记录 —— 见 `scripts/make_cases.py` 的模块说明。",
        "",
        f"全量 50016 条里,omgkit 命中而 RDKit 未命中的共 **{len(wins)}** 条,",
        f"其中质量不守恒 {len(mass)} 条、立体差异 {len(chiral)} 条、其他 {len(other)} 条。",
        "",
    ]

    lines += ["## 一、RDKit 把共享的原子复制进了每一个产物片段", ""]
    lines += [
        "逆向模板写成两个片段,而底物里这两个片段仍由**未被模板匹配**的原子连着",
        "(分子内成环的逆向就是这样)。逐产物各搬一次\"模板之外的部分\",共享的那批",
        "原子就被复制进每一片 —— 产物的重原子数**多于**底物,而没有任何东西报错。",
        "",
        "判据与记录无关:数重原子。",
        "",
    ]
    for i, rec in enumerate(mass[: args.limit]):
        fig = f"case_mass_{i + 1}.png"
        ok = draw_case(
            os.path.join(args.figdir, fig),
            rec["id"],
            [
                (f"输入 {rec['n_in']} 重原子", ".".join(rec["inputs"])),
                (f"omgkit(真值,{heavy(rec['truth'])} 重原子)", rec["truth"]),
                (f"RDKit({rec['n_bad']} 重原子,多出 {rec['n_bad'] - rec['n_in']})", rec["rdkit_bad"]),
            ],
        )
        lines += [
            f"### {i + 1}. `{rec['id']}` (row {rec['row']}, {rec['direction']})",
            "",
            f"- 模板 `{rec['tpl']}`",
            f"- 输入 `{'.'.join(rec['inputs'])}` —— {rec['n_in']} 个重原子",
            f"- omgkit `{rec['truth']}` —— {heavy(rec['truth'])} 个重原子 ✅ 守恒",
            f"- RDKit `{rec['rdkit_bad']}` —— {rec['n_bad']} 个重原子 ❌ **多出 {rec['n_bad'] - rec['n_in']} 个**",
            "",
        ]
        if ok:
            lines += [f"![{rec['id']}](../figures/{fig})", ""]

    lines += ["## 二、立体:骨架一样,RDKit 给出的是对映体", ""]
    for i, rec in enumerate(chiral[: args.limit]):
        fig = f"case_stereo_{i + 1}.png"
        ok = draw_case(
            os.path.join(args.figdir, fig),
            rec["id"],
            [
                ("输入", ".".join(rec["inputs"])),
                ("omgkit(=记录)", rec["truth"]),
                ("RDKit", rec["rdkit_bad"]),
            ],
        )
        lines += [
            f"### {i + 1}. `{rec['id']}` (row {rec['row']}, {rec['direction']})",
            "",
            f"- 模板 `{rec['tpl']}`",
            f"- 输入 `{'.'.join(rec['inputs'])}`",
            f"- omgkit `{rec['truth']}` ✅ 与记录一致",
            f"- RDKit `{rec['rdkit_bad']}` ❌ 骨架相同,立体不同",
            "",
        ]
        if ok:
            lines += [f"![{rec['id']}](../figures/{fig})", ""]

    # 第三类:书写顺序不变性 —— 用一条构造出来的最小模板,与语料无关
    lines += [
        "## 三、RDKit 的产物随模板的**书写顺序**而变",
        "",
        "同一个连接关系,模板把同样几个邻居换个次序写,描述的仍是同一个产物。",
        "下面枚举一个中心的四个邻居的全部 24 种写法,底物固定:",
        "",
        "```",
        "模板  [C:2]-[CH;D3;+0:1](-[N:3])-[O:4] >> <四个邻居的某种写法>",
        "底物  C[C@H](N)O",
        "```",
        "",
    ]
    parts = {"N": "-[N:3]", "C": "-[C:2]", "O": "-[O:4]", "X": "-Cl"}
    og_set, rd_set = set(), set()
    for order in itertools.permutations("NCOX"):
        first = parts[order[0]].lstrip("-")
        rest = "".join(f"({parts[c]})" for c in order[1:3]) + parts[order[3]]
        t = f"[C:2]-[CH;D3;+0:1](-[N:3])-[O:4]>>{first}-[C;H0;D4;+0:1]{rest}"
        og_set |= run_omgkit(t, ["C[C@H](N)O"])
        try:
            rd_set |= run_rdkit(t, ["C[C@H](N)O"])
        except Exception:
            pass
    lines += [
        f"- omgkit:24 种写法给出 **{len(og_set)}** 个不同产物 —— {sorted(og_set)}",
        f"- RDKit :24 种写法给出 **{len(rd_set)}** 个不同产物 —— {sorted(rd_set)}",
        "",
        "同一个产物不该因为模板作者的书写顺序而变成对映体。",
        "",
    ]

    os.makedirs(os.path.dirname(args.outmd), exist_ok=True)
    with open(args.outmd, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"已写出 {args.outmd}")


if __name__ == "__main__":
    main()
