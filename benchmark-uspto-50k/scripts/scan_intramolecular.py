"""两个引擎都跑不了的那批反应,到底是什么,换一种契约能拿回多少。

# 缺口

`run_reactants` / `RunReactants` 的契约是"N 个反应物模板 ↔ N 个输入分子,按位置
对应"。模板有 2-3 个片段、而记录里参与反应的只有 1-2 个分子时,这个契约表达
不了 —— 那些片段落在**同一个分子**上(分子内成环),或者落在同一个多组分物种
上(盐)。两个引擎都只能交白卷。

# 要量三件事

  一、这批反应到底有多少、长什么样(片段数 vs 分子数的分布)
  二、换成"把输入合成一张图、整条反应物侧作为一个不连通的查询去匹配"之后,
      能不能给出记录里的真值 —— 用 rdchiral 量,它正是这么做的
  三、拿回来的那些,是不是真的分子内(产物比反应物少一个分子)

第二件是决定要不要动 API 的依据:如果换了契约也拿不回来,那就只是把白卷换成
错卷,不值得为它改接口。
"""

import argparse
import json
import os
from collections import Counter

from rdkit import Chem, RDLogger
from rdkit.Chem import rdChemReactions

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def canon(smi):
    m = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(m) if m else None


def truth_of(rec, direction):
    """与基准同一套真值口径:正向取产物,逆向取参与反应的反应物集合。"""
    smis = [rec["prod"]] if direction == "fwd" else rec["reactants"]
    cs = [canon(s) for s in smis]
    if any(c is None for c in cs):
        return None
    return ".".join(sorted(cs))


def n_fragments(smarts, direction):
    """反应物侧有几个模板片段。"""
    try:
        rx = rdChemReactions.ReactionFromSmarts(smarts)
        rx.Initialize()
    except Exception:
        return None
    return rx.GetNumReactantTemplates()


def rdchiral_recovers(smarts, inputs, truth):
    """把输入合成一张图交给 rdchiral,看能不能给出真值。

    结论:**跑不了**。`rdchiralRun` 内部是 `rxn.RunReactants((合并后的分子,))`
    —— 传的是一个 1 元组,所以 RDKit 仍然要求模板只有一个反应物片段。片段数
    多于 1 时同样抛 "Number of reactants provided does not match"。
    """
    try:
        from rdchiral.main import rdchiralReaction, rdchiralReactants, rdchiralRun
    except Exception:
        return None, "rdchiral 不可用"
    try:
        rxn = rdchiralReaction(smarts)
        reactants = rdchiralReactants(".".join(inputs))
        outs = rdchiralRun(rxn, reactants)
    except Exception as e:
        return None, f"{type(e).__name__}:{str(e)[:60]}"
    got = set()
    for o in outs:
        c = canon(o)
        if c:
            got.add(".".join(sorted(c.split("."))))
    return (truth in got), f"{len(got)} 组"


def omgkit_recovers(smarts, inputs, truth):
    """omgkit 的 run_on_substrate:输入拼成一张图,各片段自由找位置且不相交。"""
    import omgkit

    try:
        rxn = omgkit.parse_reaction(smarts)
        mols = []
        for s in inputs:
            m = omgkit.parse_smiles(s)
            m.sanitize()
            mols.append(m)
        outs = rxn.run_on_substrate(mols, max_products=2000)
    except Exception as e:
        return None, f"{type(e).__name__}:{str(e)[:60]}"
    got = set()
    for oc in outs:
        parts = []
        ok = True
        for p in oc.products:
            q = p.copy()
            try:
                q.sanitize()
                c = canon(q.to_smiles())
                if c is None:
                    ok = False
                else:
                    parts.append(c)
            except Exception:
                ok = False
        if ok:
            got.add(".".join(sorted(parts)))
    return (truth in got), f"{len(got)} 组"


def attribute_misses(rows, raw_csv):
    """没命中的那些逐个立体中心定责,判据与主基准同一套。

    **不能按 SMILES 标记判。** `@`/`@@` 是相对**邻居书写顺序**说的,而反应会
    改变邻居顺序 —— 同一个空间构型,反应前后的标记完全可以不同。所以复用
    `explain_misses.analyse`:补上显式氢,再比标记加真实邻居序的置换宇称。
    """
    import csv
    import sys

    sys.path.insert(0, HERE)
    from explain_misses import analyse

    miss = [r for r in rows if r.get("omgkit_hit") is False]
    if not miss:
        return
    want_rows = {r["row"] for r in miss}
    raw = {}
    with open(raw_csv) as fh:
        # 按 enumerate 的下标取,不按 CSV 的第一列 —— 那一列不是行号
        for i, row in enumerate(csv.DictReader(fh)):
            if i in want_rows:
                raw[i] = row["rxn_smiles"]

    print(f"== 没命中的 {len(miss)} 条逐个立体中心定责 ==")
    tally = Counter()
    for r in miss:
        # `closest`:去掉立体之后与真值相同的那一个预测
        import omgkit

        target = canon(r["truth"])
        tm = Chem.MolFromSmiles(target)
        Chem.RemoveStereochemistry(tm)
        flat_truth = Chem.MolToSmiles(tm)
        closest = None
        rxn = omgkit.parse_reaction(r["smarts"])
        mols = []
        for s in r["inputs"]:
            m = omgkit.parse_smiles(s)
            m.sanitize()
            mols.append(m)
        for oc in rxn.run_on_substrate(mols, max_products=2000):
            parts = []
            for p in oc.products:
                q = p.copy()
                try:
                    q.sanitize()
                    parts.append(canon(q.to_smiles()))
                except Exception:
                    parts = None
                    break
            if not parts or any(x is None for x in parts):
                continue
            s = ".".join(sorted(parts))
            pm = Chem.MolFromSmiles(s)
            if pm is None:
                continue
            Chem.RemoveStereochemistry(pm)
            if Chem.MolToSmiles(pm) == flat_truth:
                closest = s
                break
        if closest is None or r["row"] not in raw:
            tally["拿不到可比的预测"] += 1
            r["verdicts"] = ["拿不到可比的预测"]
            continue
        rec = {
            "direction": r["dir"],
            "closest": closest,
            "retro": r["smarts"],
            "fwd": r["smarts"],
        }
        err, verdicts = analyse(rec, raw[r["row"]])
        r["verdicts"] = verdicts if not err else [err]
        for v in r["verdicts"]:
            tally[v[0] if isinstance(v, (list, tuple)) else v] += 1
    for k, v in tally.most_common():
        print(f"  {k:36s} {v}")
    print()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tpl", default=os.path.join(ROOT, "data", "templates.jsonl"))
    ap.add_argument("--raw", default=os.path.join(ROOT, "data", "uspto50k_raw.csv"))
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "intramolecular.json"))
    args = ap.parse_args()

    tally = Counter()
    rows = []
    with open(args.tpl) as fh:
        for line in fh:
            rec = json.loads(line)
            for d in ("fwd", "retro"):
                inputs = rec["reactants"] if d == "fwd" else [rec["prod"]]
                nf = n_fragments(rec[d], d)
                if nf is None:
                    tally[f"{d}_unparsable"] += 1
                    continue
                tally[f"{d}_records"] += 1
                if nf <= len(inputs):
                    continue  # 契约跑得了
                tally[f"{d}_contract_gap"] += 1
                tally[f"{d}_gap_{nf}tpl_{len(inputs)}mol"] += 1
                rows.append(
                    {
                        "row": rec["row"],
                        "id": rec["id"],
                        "dir": d,
                        "n_tpl": nf,
                        "n_mol": len(inputs),
                        "smarts": rec[d],
                        "inputs": inputs,
                        "truth": truth_of(rec, d),
                    }
                )

    print("== 契约跑不了的分布 ==")
    for k, v in sorted(tally.items()):
        if "gap" in k:
            print(f"  {k:34s} {v}")
    print(f"\n合计 {len(rows)} 组\n")

    n = len(rows)
    for tag, fn, key in (
        ("rdchiral", rdchiral_recovers, "rdchiral"),
        ("omgkit run_on_substrate", omgkit_recovers, "omgkit"),
    ):
        print(f"== 换成'合成一张图'之后能拿回多少({tag})==")
        ok = fail = err = 0
        for r in rows:
            hit, note = fn(r["smarts"], r["inputs"], r["truth"])
            r[f"{key}_hit"], r[f"{key}_note"] = hit, note
            if hit is None:
                err += 1
            elif hit:
                ok += 1
            else:
                fail += 1
        print(f"  命中真值   {ok:4d}  ({100 * ok / n:.1f}%)")
        print(f"  没命中     {fail:4d}  ({100 * fail / n:.1f}%)")
        print(f"  跑不起来   {err:4d}  ({100 * err / n:.1f}%)\n")
    attribute_misses(rows, args.raw)

    print("== omgkit 命中的例子 ==")
    shown = 0
    for r in rows:
        if not r["omgkit_hit"] or shown >= 4:
            continue
        shown += 1
        print(f"  行 {r['row']} {r['id']} [{r['dir']}]  模板 {r['n_tpl']} 片段 / 输入 {r['n_mol']} 分子")
        print(f"     模板 {r['smarts'][:118]}")
        print(f"     输入 {r['inputs']}")
        print(f"     真值 {r['truth'][:100]}")

    with open(args.out, "w") as fh:
        json.dump({"tally": dict(tally), "rows": rows}, fh, ensure_ascii=False, indent=1)
    print(f"\n已写出 {args.out}")


if __name__ == "__main__":
    main()
