"""单条深挖:把模板、输入、带映射号的产物、真值摆在一起,按映射号对齐立体。

用法:python diag_case.py <row> <fwd|retro>
"""

import json
import os
import sys

from rdkit import Chem, RDLogger

import omgkit

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
TETRA = (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW)


def tags_by_map(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return {}
    m = Chem.AddHs(m)
    out = {}
    for a in m.GetAtoms():
        if a.GetChiralTag() in TETRA and a.GetAtomMapNum():
            out[a.GetAtomMapNum()] = (
                str(a.GetChiralTag()),
                [n.GetAtomMapNum() for n in a.GetNeighbors()],
            )
    return out


def main():
    row, direction = sys.argv[1], sys.argv[2]
    tpl = None
    for line in open(os.path.join(ROOT, "data", "templates.jsonl")):
        r = json.loads(line)
        if str(r.get("row")) == row:
            tpl = r
            break
    raw = None
    import csv

    with open(os.path.join(ROOT, "data", "uspto50k_raw.csv")) as fh:
        for i, r in enumerate(csv.DictReader(fh)):
            if i == int(row):
                raw = r["rxn_smiles"]
                break

    smarts = tpl[direction]
    inputs = tpl["reactants"] if direction == "fwd" else [tpl["prod"]]
    print(f"记录 {tpl['id']} row={row} {direction}")
    print(f"  模板 {smarts}")
    print(f"  输入 {inputs}")
    lhs, _, rhs = raw.split(">")
    src_s, dst_s = (lhs, rhs) if direction == "fwd" else (rhs, lhs)
    st, dt = tags_by_map(src_s), tags_by_map(dst_s)
    print(f"  输入侧标记(按映射号) {st}")
    print(f"  目标侧标记(按映射号) {dt}")

    r = omgkit.parse_reaction(smarts)
    mols = []
    for s in inputs:
        m = omgkit.parse_smiles(s)
        m.sanitize()
        mols.append(m)
    import itertools

    n = r.num_reactant_templates
    seen = set()
    for perm in itertools.permutations(mols, n):
        for oc in r.run(list(perm), atom_mapping=True):
            ps = []
            for p in oc.products:
                q = p.copy()
                q.sanitize()
                ps.append(q.to_smiles())
            s = ".".join(ps)
            if s in seen:
                continue
            seen.add(s)
            print(f"  omgkit(带映射号) {s}")
            print(f"    标记 {tags_by_map(s)}")


if __name__ == "__main__":
    main()
