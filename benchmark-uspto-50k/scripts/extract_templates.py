"""从 USPTO-50k 抽反应模板,顺带把基准要用的真值和输入都固化下来。

语料自带原子映射(探路脚本量过:200/200 完整,产物映射号全部能在反应物侧找到),
所以不需要再跑一遍 atom-atom mapping。**但不能假定**:这里逐条查,缺映射的
单独记一类,真出现了再补标注。

产出 data/templates.jsonl,每行:
  id            记录号
  cls           反应类别(语料自带 1..10)
  retro         逆向模板 产物>>反应物 (rdchiral 抽取的原始方向)
  fwd           正向模板 反应物>>产物 (retro 两侧对调)
  prod          产物 SMILES,已去映射号,RDKit 规范式
  reactants     **参与反应的**反应物 SMILES 列表,已去映射号,RDKit 规范式
                判据:该分子里有原子的映射号出现在产物中。旁观的试剂不算 ——
                模板本来就造不出它们,把它们计进真值等于给两边都判死刑。
  spectators    被上面这条判据排除掉的分子,留着备查
  n_heavy_r     参与反应的反应物重原子数合计
  n_heavy_p     产物重原子数

时间开销主要在 rdchiral 的抽取(每条几十毫秒),所以抽取与跑基准分成两步,
中间结果落盘,基准可以反复重跑而不必重抽。
"""

import argparse
import csv
import json
import os
import sys
import time
from collections import Counter

from rdkit import Chem, RDLogger
from rdchiral.template_extractor import extract_from_reaction

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def canon(smi):
    """去映射号 + RDKit 规范式。解析不了返回 None。"""
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    for a in m.GetAtoms():
        a.SetAtomMapNum(0)
    return Chem.MolToSmiles(m)


def heavy(smi):
    m = Chem.MolFromSmiles(smi)
    return 0 if m is None else m.GetNumHeavyAtoms()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default=os.path.join(ROOT, "data", "uspto50k_raw.csv"))
    ap.add_argument("--out", default=os.path.join(ROOT, "data", "templates.jsonl"))
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    stats = Counter()
    t0 = time.time()
    with open(args.src) as fh, open(args.out, "w") as out:
        for i, row in enumerate(csv.DictReader(fh)):
            if args.limit and i >= args.limit:
                break
            stats["rows"] += 1
            rxn = row["rxn_smiles"]
            if rxn.count(">") != 2:
                stats["bad_rxn_field"] += 1
                continue
            lhs, agents, rhs = rxn.split(">")

            pm = Chem.MolFromSmiles(rhs)
            if pm is None:
                stats["product_unparsable"] += 1
                continue
            pmaps = {a.GetAtomMapNum() for a in pm.GetAtoms() if a.GetAtomMapNum()}
            if not pmaps:
                # 真没有原子映射 —— 记下来,后面单独补标注
                stats["needs_atom_mapping"] += 1
                out.write(
                    json.dumps({"id": row["id"], "row": i, "err": "unmapped"}) + "\n"
                )
                continue
            stats["mapped"] += 1

            frags = lhs.split(".")
            reactants, spectators = [], []
            for frag in frags:
                fm = Chem.MolFromSmiles(frag)
                if fm is None:
                    stats["reactant_frag_unparsable"] += 1
                    continue
                if any(a.GetAtomMapNum() in pmaps for a in fm.GetAtoms()):
                    reactants.append(frag)
                else:
                    spectators.append(frag)
            if not reactants:
                stats["no_contributing_reactant"] += 1
                continue

            try:
                res = extract_from_reaction(
                    {"reactants": lhs, "products": rhs, "_id": row["id"]}
                )
            except Exception as e:  # rdchiral 偶尔会在畸形记录上抛
                stats["extract_exception"] += 1
                out.write(
                    json.dumps(
                        {"id": row["id"], "row": i, "err": f"extract:{type(e).__name__}"}
                    )
                    + "\n"
                )
                continue
            if not res or not res.get("reaction_smarts"):
                stats["extract_empty"] += 1
                out.write(
                    json.dumps({"id": row["id"], "row": i, "err": "extract:empty"})
                    + "\n"
                )
                continue

            retro = res["reaction_smarts"]
            if retro.count(">>") != 1:
                stats["template_shape"] += 1
                continue
            a, b = retro.split(">>")
            fwd = f"{b}>>{a}"

            cp = canon(rhs)
            cr = [canon(x) for x in reactants]
            if cp is None or any(x is None for x in cr):
                stats["canon_fail"] += 1
                continue

            rec = {
                "id": row["id"],
                "row": i,
                "cls": row.get("class", ""),
                "retro": retro,
                "fwd": fwd,
                "prod": cp,
                "reactants": cr,
                "spectators": [canon(x) for x in spectators],
                "agents": agents,
                "n_heavy_r": sum(heavy(x) for x in cr),
                "n_heavy_p": heavy(cp),
            }
            out.write(json.dumps(rec) + "\n")
            stats["ok"] += 1

            if stats["rows"] % 2000 == 0:
                el = time.time() - t0
                print(
                    f"  {stats['rows']} 条 / {el:.0f}s  ok={stats['ok']}",
                    file=sys.stderr,
                    flush=True,
                )

    print("== 抽取统计 ==")
    for k, v in sorted(stats.items()):
        print(f"{k:30s} {v}")
    print(f"耗时 {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
