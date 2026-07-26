"""按记录号复现一条:把输入、模板、两个引擎的输出、真值摆在一起。

用法:python repro.py US07994164B2 [fwd|retro]
"""

import json
import os
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import rdChemReactions

import omgkit

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def rd_canon(s):
    m = Chem.MolFromSmiles(s)
    return Chem.MolToSmiles(m) if m else f"<读不回:{s}>"


def main():
    # 语料里同一个专利号会对应几十条反应(最多 79 条),所以**按行号**取,
    # 按 id 取会静默拿到另一条。行号就是 templates.jsonl 里的 row。
    key = sys.argv[1]
    want_dir = sys.argv[2] if len(sys.argv) > 2 else None

    tpl = None
    for line in open(os.path.join(ROOT, "data", "templates.jsonl")):
        r = json.loads(line)
        if (str(r.get("row")) == key) or (key.startswith("US") and r.get("id") == key):
            tpl = r
            break
    if tpl is None:
        print("找不到记录")
        return

    print(f"记录 {tpl['id']}  row={tpl['row']}  class={tpl['cls']}")
    print(f"  反应物 {tpl['reactants']}")
    print(f"  产物   {tpl['prod']}")
    print(f"  旁观者 {tpl['spectators']}")
    print(f"  逆向模板 {tpl['retro']}")
    print(f"  正向模板 {tpl['fwd']}")

    for direction in ("fwd", "retro"):
        if want_dir and direction != want_dir:
            continue
        smarts = tpl[direction]
        inputs = tpl["reactants"] if direction == "fwd" else [tpl["prod"]]
        truth = tpl["prod"] if direction == "fwd" else ".".join(tpl["reactants"])
        print(f"\n=== {direction} ===")
        print(f"  输入 {inputs}")
        print(f"  真值 {rd_canon(truth)}")

        r = omgkit.parse_reaction(smarts)
        n = r.num_reactant_templates
        import itertools

        og = set()
        for perm in itertools.permutations(
            [_prep_og(s) for s in inputs], n
        ):
            for oc in r.run(list(perm)):
                q = [p.copy() for p in oc.products]
                for x in q:
                    x.sanitize()
                og.add(rd_canon(".".join(x.to_smiles() for x in q)))
        print(f"  omgkit {sorted(og)}")

        rx = rdChemReactions.ReactionFromSmarts(smarts)
        rx.Initialize()
        rd = set()
        for perm in itertools.permutations([Chem.MolFromSmiles(s) for s in inputs], n):
            for ps in rx.RunReactants(perm):
                out = []
                ok = True
                for p in ps:
                    q = Chem.Mol(p)
                    for a in q.GetAtoms():
                        a.SetAtomMapNum(0)
                    try:
                        Chem.SanitizeMol(q)
                        out.append(Chem.MolToSmiles(q))
                    except Exception:
                        ok = False
                if ok:
                    rd.add(rd_canon(".".join(out)))
        print(f"  rdkit  {sorted(rd)}")


def _prep_og(s):
    m = omgkit.parse_smiles(s)
    m.sanitize()
    return m


if __name__ == "__main__":
    main()
