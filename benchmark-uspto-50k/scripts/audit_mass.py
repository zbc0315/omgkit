"""对**全部**输出查质量守恒,而不是只查未命中的那些。

# 为什么必须重查

"骨架 50016 条无一出错"这个结论是在 1288 条**未命中**上归的因 —— 它回答的是
"没中的时候为什么没中",而不是"输出有没有错"。一条反应完全可以既输出了正确
产物(于是算命中)、又同时输出了撕坏的产物;按未命中归因,后者永远看不见。

撕环缺陷正是这个形状:模板把环写成开链路径时,某些匹配位点撕环、另一些不撕,
命中率一点没掉。所以要按**每一个 outcome** 查,不按反应查。

# 判据:重原子数

    应得 = 输入重原子数 + 模板净增(产物模板里没有映射号的重原子
                                   − 反应物模板里没有映射号的重原子)

多出来 = **复制**,这是最严重的一档:拓扑错了而什么都不报错。
少掉了 = 丢原子,再分两种:

    整个不连通组分被丢掉  —— 已知约定(搬运是从匹配到的原子出发遍历,
                              完全不连通的组分永远走不到)
    别的                  —— 见下

# "别的"那一档不能直接算成缺陷

模板删掉一个原子时,**挂在它身上、又不通过别的路连回保留部分**的那些原子会一并
失去落脚点。叔丁酯水解是标准例子:模板写 `C-C-[O:1]-[C:2]=[O:3]`,删掉的 `C-C`
是叔丁基的一个甲基加季碳,而季碳上另外两个甲基随之孤立 —— 一次少掉 4 个而不是 2 个。

判据能不能认出这一档,不影响结论,因为**两个引擎在这一档上数字完全相同**。
凡是两边数字一模一样的档,成因一定在模板或约定,不可能是某一方的实现差异 ——
实现差异不会撞出相同的数。真正要看的是**不对称**的那些档。

两个引擎都跑,口径完全一致 —— 否则"多出 1823 个重原子"这种话就不可比。
"""

import argparse
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
MAX_PRODUCTS = 1000


def heavy_smiles(s):
    m = Chem.MolFromSmiles(s, sanitize=False)
    return 0 if m is None else m.GetNumHeavyAtoms()


def template_delta(smarts):
    """产物模板新建的重原子数 − 反应物模板删掉的重原子数。

    "被删掉"有**两种**,漏掉任何一种都会把好反应误报成丢原子:

      没有映射号            产物侧根本没有它的位置
      有映射号、但产物侧没有  映射号在反应物侧出现、产物侧不出现,一样是删掉

    第二种在正向模板里很常见(脱保护基:保护基上的原子带着映射号却不进产物)。
    只按第一种算,实测 800 条里就有 16 条被误报 —— 而且**两个引擎数字完全相同**,
    这正是"错在判据不在实现"的signature:引擎差异不可能一模一样。
    """
    rx = rdChemReactions.ReactionFromSmarts(smarts)
    rx.Initialize()

    def maps(getter, n):
        out = set()
        for i in range(n):
            for a in getter(i).GetAtoms():
                if a.GetAtomMapNum():
                    out.add(a.GetAtomMapNum())
        return out

    r_maps = maps(rx.GetReactantTemplate, rx.GetNumReactantTemplates())
    p_maps = maps(rx.GetProductTemplate, rx.GetNumProductTemplates())

    def count(getter, n, keep):
        """`keep` 是对侧出现的映射号集合;不在其中的原子就是新建/删掉的。"""
        c = 0
        for i in range(n):
            for a in getter(i).GetAtoms():
                if a.GetAtomMapNum() not in keep:
                    c += 1
        return c

    made = count(rx.GetProductTemplate, rx.GetNumProductTemplates(), r_maps)
    gone = count(rx.GetReactantTemplate, rx.GetNumReactantTemplates(), p_maps)
    return made - gone, rx


def component_sizes(smis):
    """输入里每个不连通组分的重原子数 —— 用来认"整个组分被丢掉"。"""
    out = []
    for s in smis:
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        for frag in Chem.GetMolFrags(m, asMols=True, sanitizeFrags=False):
            out.append(frag.GetNumHeavyAtoms())
    return out


def subset_sum_reachable(sizes, target):
    """target 能不能由若干个组分大小加出来 —— 判"丢掉的正好是整组分"。"""
    if target == 0:
        return True
    reach = {0}
    for s in sizes:
        reach |= {r + s for r in reach if r + s <= target}
    return target in reach


def run_engine(engine, smarts, inputs, rx, delta):
    """返回每个 outcome 的 (应得重原子数, 实得重原子数, 本次排列的组分大小)。

    应得数必须**按每个排列**算:正向记录里参与反应的分子可能多于模板片段数,
    引擎每次只吃 n 个,没被吃进去的那些当然不该计进应得数。按全部输入算会
    把每一条这样的反应都误报成"丢原子"。实得数净化不过时给 None。
    """
    out = []
    if engine == "rdkit":
        n = rx.GetNumReactantTemplates()
        mols = [Chem.MolFromSmiles(s) for s in inputs]
        if any(m is None for m in mols) or len(mols) < n:
            return None
        for perm in itertools.permutations(mols, n):
            want = sum(m.GetNumHeavyAtoms() for m in perm) + delta
            sizes = [
                f.GetNumHeavyAtoms()
                for m in perm
                for f in Chem.GetMolFrags(m, asMols=True, sanitizeFrags=False)
            ]
            for ps in rx.RunReactants(perm, MAX_PRODUCTS):
                tot, ok = 0, True
                for p in ps:
                    q = Chem.Mol(p)
                    for a in q.GetAtoms():
                        a.SetAtomMapNum(0)
                    try:
                        Chem.SanitizeMol(q)
                    except Exception:
                        ok = False
                    tot += q.GetNumHeavyAtoms()
                out.append((want, tot if ok else None, sizes))
    else:
        r = omgkit.parse_reaction(smarts)
        n = r.num_reactant_templates
        mols, smis = [], []
        for s in inputs:
            m = omgkit.parse_smiles(s)
            m.sanitize()
            mols.append(m)
            smis.append(s)
        if len(mols) < n:
            return None
        for perm in itertools.permutations(range(len(mols)), n):
            want = sum(heavy_smiles(smis[i]) for i in perm) + delta
            sizes = component_sizes([smis[i] for i in perm])
            for oc in r.run([mols[i] for i in perm], max_products=MAX_PRODUCTS):
                tot, ok = 0, True
                for p in oc.products:
                    q = p.copy()
                    try:
                        q.sanitize()
                        tot += heavy_smiles(q.to_smiles())
                    except Exception:
                        ok = False
                out.append((want, tot if ok else None, sizes))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tpl", default=os.path.join(ROOT, "data", "templates.jsonl"))
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--engines", default="omgkit,rdkit")
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "mass_audit.jsonl"))
    ap.add_argument("--summary", default=os.path.join(ROOT, "results", "mass_audit.json"))
    args = ap.parse_args()

    engines = args.engines.split(",")
    tally = Counter()
    sink = open(args.out, "w")
    n = 0
    with open(args.tpl) as fh:
        for line in fh:
            rec = json.loads(line)
            if args.limit and n >= args.limit:
                break
            n += 1
            if n % 5000 == 0:
                print(f"  {n} 条…", flush=True)
            for d in ("fwd", "retro"):
                inputs = rec["reactants"] if d == "fwd" else [rec["prod"]]
                try:
                    delta, rx = template_delta(rec[d])
                except Exception:
                    tally[f"{d}_tpl_unparsable"] += 1
                    continue
                for e in engines:
                    try:
                        got = run_engine(e, rec[d], inputs, rx, delta)
                    except Exception:
                        tally[f"{e}_{d}_run_error"] += 1
                        continue
                    if got is None:
                        continue
                    tally[f"{e}_{d}_outcomes"] += len(got)
                    for want, tot, sizes in got:
                        if tot is None:
                            tally[f"{e}_{d}_unsanitizable"] += 1
                            continue
                        if tot == want:
                            tally[f"{e}_{d}_exact"] += 1
                        elif tot > want:
                            tally[f"{e}_{d}_duplicated"] += 1
                            tally[f"{e}_{d}_extra_atoms"] += tot - want
                            sink.write(json.dumps({"row": rec["row"], "id": rec["id"],
                                                   "dir": d, "engine": e, "kind": "多",
                                                   "want": want, "got": tot}) + "\n")
                        else:
                            lost = want - tot
                            if subset_sum_reachable(sizes, lost):
                                tally[f"{e}_{d}_lost_whole_component"] += 1
                            else:
                                tally[f"{e}_{d}_lost_other"] += 1
                                tally[f"{e}_{d}_lost_atoms"] += lost
                                sink.write(json.dumps({"row": rec["row"], "id": rec["id"],
                                                       "dir": d, "engine": e, "kind": "少",
                                                       "want": want, "got": tot}) + "\n")
    sink.close()
    print(f"\n扫了 {n} 条记录\n")
    for k, v in sorted(tally.items()):
        print(f"  {k:36s} {v}")

    # 两个引擎对照:相同的档是模板/约定,不对称的档才是实现差异
    print("\n== 两个引擎逐档对照 ==")
    print(f"  {'档':28s} {'omgkit':>10s} {'RDKit':>10s}  {'':4s}")
    kinds = ["exact", "duplicated", "extra_atoms", "unsanitizable",
             "lost_whole_component", "lost_other", "lost_atoms"]
    for d in ("fwd", "retro"):
        for kind in kinds:
            a, b = tally.get(f"omgkit_{d}_{kind}", 0), tally.get(f"rdkit_{d}_{kind}", 0)
            if not a and not b:
                continue
            mark = "相同" if a == b else "**不对称**"
            print(f"  {d + '·' + kind:28s} {a:10d} {b:10d}  {mark}")
    with open(args.summary, "w") as fh:
        json.dump({"n_records": n, "tally": dict(tally)}, fh, ensure_ascii=False, indent=2)
    print(f"\n已写出 {args.summary} 与 {args.out}")


if __name__ == "__main__":
    main()
