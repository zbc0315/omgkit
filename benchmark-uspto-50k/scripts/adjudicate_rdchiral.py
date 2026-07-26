"""拿 rdchiral 当第三方裁判,裁定"疑似引擎缺陷"的那几条。

# 为什么它有资格当裁判

这批模板就是 rdchiral 从这些反应里抽出来的,模板里的手性标记也是它写的。
所以"这条模板该怎么读"这个问题,它自己的应用器(`rdchiralRun`)给出的就是
**作者本意**。

于是:

- rdchiral 也还原不出原记录 → 模板本身没把这个构型编码进去,任何忠实执行
  模板的实现都做不到。不是 omgkit 的缺陷。
- rdchiral 还原得出而 omgkit 还原不出 → omgkit 没做到模板说的事,是缺陷。

这一档不能只看 RDKit 的 `RunReactants`:它对反应里的手性处理本来就弱
(同一批模板换个书写次序,它给出的产物会变),两边都错时说明不了问题。
"""

import argparse
import json
import os

from rdkit import Chem, RDLogger

from rdchiral.main import rdchiralReactants, rdchiralReaction, rdchiralRun

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def canon(s):
    m = Chem.MolFromSmiles(s)
    return Chem.MolToSmiles(m) if m else None


def run_rdchiral(tpl, inputs):
    """把模板作用在输入上;输入多于一个时拼成一个多片段分子。"""
    try:
        rxn = rdchiralReaction(tpl)
        rct = rdchiralReactants(".".join(inputs))
        return {canon(o) for o in rdchiralRun(rxn, rct) if canon(o)}
    except Exception as e:
        return {f"<异常:{type(e).__name__}:{e}>"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default=os.path.join(ROOT, "results", "miss_verdicts.jsonl"))
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "rdchiral_verdicts.jsonl"))
    args = ap.parse_args()

    tpls = {}
    for line in open(os.path.join(ROOT, "data", "templates.jsonl")):
        r = json.loads(line)
        if "row" in r and "retro" in r:
            tpls[r["row"]] = r

    n_tpl_limit = n_engine = 0
    with open(args.src) as fh, open(args.out, "w") as out:
        for line in fh:
            rec = json.loads(line)
            if not rec["verdict"].startswith("引擎"):
                continue
            t = tpls.get(rec["row"])
            if t is None:
                continue
            smarts = t[rec["direction"]]
            inputs = t["reactants"] if rec["direction"] == "fwd" else [t["prod"]]
            truth = canon(rec["truth"])
            got = run_rdchiral(smarts, inputs)
            hit = truth in got
            verdict = "引擎:rdchiral 做得到而 omgkit 没做到" if hit else "模板:rdchiral 自己也还原不出"
            if hit:
                n_engine += 1
            else:
                n_tpl_limit += 1
            print(f"row={rec['row']:6d} {rec['id']:18s} {rec['direction']:6s} {verdict}")
            print(f"    真值     {truth}")
            print(f"    omgkit   {rec['closest']}")
            for g in sorted(got)[:2]:
                print(f"    rdchiral {g}")
            out.write(json.dumps({**rec, "rdchiral": sorted(got), "final": verdict}) + "\n")

    print(f"\n模板极限 {n_tpl_limit} 条,omgkit 缺陷 {n_engine} 条")


if __name__ == "__main__":
    main()
