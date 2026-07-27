"""为"分子内不可预知"与"盐"两个论点测出定量支撑,写成 results/intra_salt.json。

# 测什么

**环大小扫描**:同一条断酰胺的逆向模板,作用在 n 元内酰胺上(n = 5…25)。
正确答案永远是**一个**开环分子。记两个引擎在两种模板写法下的产物分子数与
重原子数,看误差怎么随 n 变。

**盐的普遍程度**:USPTO-50k 里带电片段、单原子/双原子离子、产物侧多片段的
记录各有多少。盐若是罕见情形,论点 B 就不值一提;若是常态,契约缺口就得认真对待。
"""

import argparse
import csv
import itertools
import json
import os
from collections import Counter

from rdkit import Chem, RDLogger
from rdkit.Chem import rdChemReactions

import omgkit

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# 断酰胺的逆向模板。两种写法:片段式与组分括号式。
PLAIN = (
    "[C:2]-[C;H0;D3;+0:1](=[O;D1;H0:3])-[NH;D2;+0:4]-[C:5]"
    ">>O-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3].[NH2;D1;+0:4]-[C:5]"
)
GROUPED = (
    "[C:2]-[C;H0;D3;+0:1](=[O;D1;H0:3])-[NH;D2;+0:4]-[C:5]"
    ">>(O-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3].[NH2;D1;+0:4]-[C:5])"
)


def heavy(s):
    m = Chem.MolFromSmiles(s)
    return 0 if m is None else m.GetNumHeavyAtoms()


def lactam(n):
    """n 元内酰胺:环上一个 C=O、一个 N、其余是 CH2。"""
    return "O=C1" + "C" * (n - 2) + "N1"


def run_og(tpl, smi):
    r = omgkit.parse_reaction(tpl)
    m = omgkit.parse_smiles(smi)
    m.sanitize()
    best = None
    for oc in r.run([m]):
        frags = []
        for p in oc.products:
            q = p.copy()
            q.sanitize()
            frags.append(q.to_smiles())
        cand = (len(frags), sum(heavy(f) for f in frags))
        if best is None or cand < best:
            best = cand
    return best


def run_rd(tpl, smi):
    rx = rdChemReactions.ReactionFromSmarts(tpl)
    rx.Initialize()
    best = None
    for ps in rx.RunReactants((Chem.MolFromSmiles(smi),)):
        frags, comps = [], 0
        for p in ps:
            q = Chem.Mol(p)
            for a in q.GetAtoms():
                a.SetAtomMapNum(0)
            try:
                Chem.SanitizeMol(q)
                s = Chem.MolToSmiles(q)
            except Exception:
                s = Chem.MolToSmiles(q)
            frags.append(s)
            mm = Chem.MolFromSmiles(s)
            comps += len(Chem.GetMolFrags(mm)) if mm else 1
        cand = (len(frags), sum(heavy(f) for f in frags), comps)
        if best is None or cand < best:
            best = cand
    return best


def ring_sweep():
    out = []
    for n in [5, 6, 7, 8, 10, 12, 15, 20, 25]:
        smi = lactam(n)
        n_in = heavy(smi)
        want = n_in + 1  # 水解加进来一个氧;正确答案永远是一个分子
        row = {"n": n, "smiles": smi, "n_in": n_in, "want_heavy": want, "want_mols": 1}
        for tag, tpl in (("plain", PLAIN), ("grouped", GROUPED)):
            og = run_og(tpl, smi)
            rd = run_rd(tpl, smi)
            row[f"og_{tag}"] = {"mols": og[0], "heavy": og[1]} if og else None
            row[f"rd_{tag}"] = (
                {"mols": rd[0], "heavy": rd[1], "components": rd[2]} if rd else None
            )
        out.append(row)
    return out


def intermolecular_check():
    """开链底物上,两种写法各给出几个产物分子。正确答案是 2。"""
    sub = "CCC(=O)NCc1ccccc1"
    res = {"smiles": sub, "n_in": heavy(sub), "want_heavy": heavy(sub) + 1, "want_mols": 2}
    for tag, tpl in (("plain", PLAIN), ("grouped", GROUPED)):
        og, rd = run_og(tpl, sub), run_rd(tpl, sub)
        res[f"og_{tag}"] = {"mols": og[0], "heavy": og[1]}
        res[f"rd_{tag}"] = {"mols": rd[0], "heavy": rd[1], "components": rd[2]}
    return res


def salt_prevalence(path):
    """语料里盐有多常见。

    USPTO 的记录用 `.` 同时表示"另一个反应物"和"同一个盐的另一半",文本上分不开。
    所以这里量的是**必然伴随反离子**的证据:
      带电片段  —— 净电荷非零的片段,一定有个反离子在同一条记录里
      小离子    —— 重原子数 ≤ 2 的片段(Cl⁻、Na⁺、HCl、TFA 的一半…)
    """
    c = Counter()
    with open(path) as fh:
        for row in csv.DictReader(fh):
            rxn = row["rxn_smiles"]
            if rxn.count(">") != 2:
                continue
            c["records"] += 1
            lhs, _, rhs = rxn.split(">")
            for side, frags in (("r", lhs.split(".")), ("p", rhs.split("."))):
                charged = small = 0
                for f in frags:
                    m = Chem.MolFromSmiles(f)
                    if m is None:
                        continue
                    q = sum(a.GetFormalCharge() for a in m.GetAtoms())
                    if q != 0:
                        charged += 1
                    if m.GetNumHeavyAtoms() <= 2:
                        small += 1
                if charged:
                    c[f"{side}_has_charged_fragment"] += 1
                if small:
                    c[f"{side}_has_small_ion"] += 1
                if charged or small:
                    c[f"{side}_saltlike"] += 1
            if len(rhs.split(".")) > 1:
                c["product_multi_fragment"] += 1
    return dict(c)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default=os.path.join(ROOT, "data", "uspto50k_raw.csv"))
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "intra_salt.json"))
    args = ap.parse_args()

    data = {
        "ring_sweep": ring_sweep(),
        "intermolecular": intermolecular_check(),
        "salt_prevalence": salt_prevalence(args.raw),
        "templates": {"plain": PLAIN, "grouped": GROUPED},
    }
    with open(args.out, "w") as fh:
        json.dump(data, fh, ensure_ascii=False, indent=2)

    print("== 环大小扫描(正确答案恒为 1 个分子)==")
    print("  n | 输入 | 应得 |  omgkit 片段式 |  RDKit 片段式  | RDKit 括号式")
    for r in data["ring_sweep"]:
        og, rdp, rdg = r["og_plain"], r["rd_plain"], r["rd_grouped"]
        print(
            f" {r['n']:2d} | {r['n_in']:4d} | {r['want_heavy']:4d} |"
            f" {og['mols']}片/{og['heavy']:3d} |"
            f" {rdp['mols']}片/{rdp['heavy']:3d} ({rdp['heavy'] - r['want_heavy']:+d}) |"
            f" {rdg['mols']}片/{rdg['heavy']:3d}"
        )
    im = data["intermolecular"]
    print("\n== 开链底物(正确答案 2 个分子)==")
    for k in ("og_plain", "rd_plain", "og_grouped", "rd_grouped"):
        v = im[k]
        print(f"  {k:12s} {v['mols']} 个产物对象,{v['heavy']} 重原子")

    print("\n== 盐的普遍程度 ==")
    n = data["salt_prevalence"]["records"]
    for k, v in sorted(data["salt_prevalence"].items()):
        if k == "records":
            continue
        print(f"  {k:28s} {v:6d}  ({100 * v / n:.1f}%)")
    print(f"\n已写出 {args.out}")


if __name__ == "__main__":
    main()
