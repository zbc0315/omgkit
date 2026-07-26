"""小样探路:先把三件事量出来,再决定整套基准怎么搭。

1. 语料到底带不带原子映射,带得全不全。
2. rdchiral 抽出来的模板,omgkit 和 rdkit **是不是都能解析** ——
   解析不了就没得比,得先知道缺口有多大。
3. 正向那一侧,模板的反应物模板个数 N 与真实反应物分子个数 M 的关系。
   RunReactants 要求两者相等,N != M 时怎么配对决定了整套基准的形状。

这个脚本只量不判,输出直接进 docs/。
"""

import csv
import re
import sys
from collections import Counter

from rdkit import Chem, RDLogger
from rdchiral.template_extractor import extract_from_reaction

import omgkit

RDLogger.DisableLog("rdApp.*")

MAP_RE = re.compile(r":\d+\]")
N = int(sys.argv[1]) if len(sys.argv) > 1 else 200
PATH = sys.argv[2] if len(sys.argv) > 2 else "../data/uspto50k_raw.csv"


def strip_maps(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    for a in m.GetAtoms():
        a.SetAtomMapNum(0)
    return Chem.MolToSmiles(m)


def main():
    stats = Counter()
    n_by_m = Counter()
    parse_fail_examples = []
    templates = []

    with open(PATH) as fh:
        rd = csv.DictReader(fh)
        for i, row in enumerate(rd):
            if i >= N:
                break
            stats["rows"] += 1
            rxn = row["rxn_smiles"]
            lhs, _, rhs = rxn.split(">")

            # ---- 1. 原子映射 ----
            pm = Chem.MolFromSmiles(rhs)
            rm = Chem.MolFromSmiles(lhs)
            if pm is None or rm is None:
                stats["rdkit_parse_fail"] += 1
                continue
            pmaps = {a.GetAtomMapNum() for a in pm.GetAtoms() if a.GetAtomMapNum()}
            rmaps = {a.GetAtomMapNum() for a in rm.GetAtoms() if a.GetAtomMapNum()}
            if not pmaps:
                stats["product_unmapped"] += 1
            if pmaps and pmaps <= rmaps:
                stats["product_maps_covered"] += 1
            heavy_unmapped_p = sum(1 for a in pm.GetAtoms() if not a.GetAtomMapNum())
            if heavy_unmapped_p:
                stats["product_has_unmapped_atom"] += 1

            # ---- 2. 模板抽取 + 双方解析 ----
            try:
                res = extract_from_reaction(
                    {"reactants": lhs, "products": rhs, "_id": row["id"]}
                )
            except Exception:
                res = None
            if not res or "reaction_smarts" not in res:
                stats["template_extract_fail"] += 1
                continue
            retro = res["reaction_smarts"]
            stats["template_ok"] += 1
            templates.append(retro)

            fwd = ">>".join(reversed(retro.split(">>")))

            for tag, tpl in (("retro", retro), ("fwd", fwd)):
                try:
                    r = Chem.rdChemReactions.ReactionFromSmarts(tpl)
                    ok_rd = r is not None and r.GetNumReactantTemplates() > 0
                except Exception:
                    ok_rd = False
                try:
                    o = omgkit.parse_reaction(tpl)
                    ok_og = o.num_reactant_templates > 0
                except Exception as e:
                    ok_og = False
                    if len(parse_fail_examples) < 8:
                        parse_fail_examples.append((tag, tpl, str(e).split("\n")[0]))
                stats[f"rdkit_parse_{tag}_ok"] += int(ok_rd)
                stats[f"omgkit_parse_{tag}_ok"] += int(ok_og)
                if ok_rd and not ok_og:
                    stats[f"only_rdkit_{tag}"] += 1
                if ok_og and not ok_rd:
                    stats[f"only_omgkit_{tag}"] += 1

            # ---- 3. N vs M ----
            try:
                rr = Chem.rdChemReactions.ReactionFromSmarts(fwd)
                n_tpl = rr.GetNumReactantTemplates()
            except Exception:
                n_tpl = -1
            # 参与反应的反应物:含有出现在产物里的映射号
            contributing = 0
            for frag in lhs.split("."):
                fm = Chem.MolFromSmiles(frag)
                if fm is None:
                    continue
                if any(a.GetAtomMapNum() in pmaps for a in fm.GetAtoms()):
                    contributing += 1
            n_by_m[(n_tpl, len(lhs.split(".")), contributing)] += 1

    print("== 计数 ==")
    for k, v in sorted(stats.items()):
        print(f"{k:34s} {v}")
    print("\n== (模板反应物数 N, 反应物分子数 M, 参与反应的 M') 分布 ==")
    for k, v in n_by_m.most_common(15):
        print(f"  N={k[0]} M={k[1]} M'={k[2]}  ->  {v}")
    print("\n== omgkit 解析失败样例 ==")
    for tag, tpl, err in parse_fail_examples:
        print(f"  [{tag}] {tpl}\n        {err}")
    print("\n== 模板样例 ==")
    for t in templates[:3]:
        print("  ", t)


if __name__ == "__main__":
    main()
