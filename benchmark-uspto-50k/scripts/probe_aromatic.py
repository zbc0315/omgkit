"""产物净化不过时,错误信息指向哪里 —— 兼一条被推翻的猜测的记录。

# 被推翻的猜测

最初的猜测是:产物是**未净化**的,芳香标志逐原子逐键从底物照抄过来;模板一旦
断掉环上的一根键,剩下的原子还带着 `aromatic` 标志却已不构成芳香体系,于是
kekulize 找不到交替单双键的指派,当场失败。按这个猜测,病根在"芳香标志没清干净"。

**猜错了。** 芳香标志确实是照抄的,但那不是病根 —— 环压根不该断。真正的错在
**键的归属规则**:子结构匹配只要求模板的每根键在底物里找得到,不要求这些原子
之间没有别的键;模板把环写成开链路径时,环闭合的那根键两端都被匹配、模板却
没匹配到它,当时的判据("两端都被匹配就不搬")把它当成模板的地盘删掉了。

芳香错误是**症状**,病因在两层之外。改掉键归属规则之后,omgkit 在全语料上
产物净化不过的反应数从 8 条降到 0 条,一行芳香相关的代码都没动。

# 怎么会错到那里去

猜测是拿**自己造的**模板验证的,其中一条是 `[c:1][c:2]>>[C:1].[C:2]` ——
它描述的是"断苯环上一根 C–C 又不补氢",化学上根本不存在的变换,得到的必然是
双自由基。拿它当探针,两个实现都会"失败",于是得出"两个实现共有的缺陷"这个
错误结论。

**判据要用语料里的真模板。** 自己造的模板可能根本不是它看上去的意思,拿它当
判据会把判据本身变成错的 —— 而错判据会把人引到本来正确的代码上去改。

# 还有用的部分

下面这段语料扫描仍然有用:它按**错误信息**给净化失败分类,而不同的错误信息
指向完全不同的成因,不能混在一起数。

    Can't kekulize    芳香标志与环的实际情况对不上
    Explicit valence  价键算错
"""

import argparse
import itertools
import json
import os
from collections import Counter

from rdkit import Chem, RDLogger
from rdkit.Chem import rdChemReactions

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def err_kind(msg):
    m = msg.lower()
    if "kekulize" in m:
        return "无法 kekulize(芳香标志与环对不上)"
    if "valence" in m:
        return "价键超限"
    if "aromatic" in m:
        return "芳香性其他"
    return "其他:" + msg.split("\n")[0][:60]


def rdkit_failures(tpl, inputs):
    """返回净化失败的 (错误类型, 未净化产物串) 列表。"""
    rx = rdChemReactions.ReactionFromSmarts(tpl)
    rx.Initialize()
    n = rx.GetNumReactantTemplates()
    mols = [Chem.MolFromSmiles(s) for s in inputs]
    if any(m is None for m in mols) or len(mols) < n:
        return []
    out = []
    for perm in itertools.permutations(mols, n):
        for ps in rx.RunReactants(perm):
            for p in ps:
                q = Chem.Mol(p)
                for a in q.GetAtoms():
                    a.SetAtomMapNum(0)
                try:
                    Chem.SanitizeMol(q)
                except Exception as e:
                    out.append((err_kind(str(e)), Chem.MolToSmiles(q)))
    return out


def corpus_scan(bench, tpls, limit):
    """语料上净化不过的那些反应,失败原因分类;顺带看反应中心是不是芳香的。"""
    rows = []
    for line in open(bench):
        r = json.loads(line)
        if "err" in r:
            continue
        for d in ("fwd", "retro"):
            v = r.get(f"rdkit_{d}")
            if v and v.get("n_bad"):
                rows.append((r["row"], d))
    print(f"RDKit 有产物净化不过的反应:{len(rows)} 条(逐条重跑取错误信息)")

    kinds, arom = Counter(), Counter()
    samples = {}
    for row, d in rows[:limit]:
        t = tpls.get(row)
        if t is None:
            continue
        inputs = t["reactants"] if d == "fwd" else [t["prod"]]
        smarts = t[d]
        lhs = smarts.split(">>")[0]
        touches_aromatic = any(c.islower() for c in lhs.replace("H", "")) or ";a" in lhs
        for kind, smi in rdkit_failures(smarts, inputs):
            kinds[kind] += 1
            arom[(kind, "反应中心含芳香" if touches_aromatic else "反应中心不含芳香")] += 1
            samples.setdefault(kind, []).append((t["id"], d, smi, smarts))
    print("\n== RDKit 净化失败的原因 ==")
    for k, v in kinds.most_common():
        print(f"  {k:34s} {v}")
    print("\n== 与反应中心是否芳香的交叉 ==")
    for k, v in arom.most_common():
        print(f"  {k[0]:34s} {k[1]:16s} {v}")
    for k, rows_ in samples.items():
        print(f"\n---- {k} 举例 ----")
        for i, dd, smi, sm in rows_[:2]:
            print(f"  {i} [{dd}]")
            print(f"    未净化产物 {smi[:120]}")
            print(f"    模板       {sm[:130]}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench", default=os.path.join(ROOT, "results", "bench.jsonl"))
    ap.add_argument("--tpl", default=os.path.join(ROOT, "data", "templates.jsonl"))
    ap.add_argument("--limit", type=int, default=200)
    args = ap.parse_args()

    tpls = {}
    for line in open(args.tpl):
        r = json.loads(line)
        if "row" in r and "retro" in r:
            tpls[r["row"]] = r

    corpus_scan(args.bench, tpls, args.limit)


if __name__ == "__main__":
    main()
