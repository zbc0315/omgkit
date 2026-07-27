"""把"omgkit 命中而 RDKit 未命中"的每一条,按 RDKit **错在哪一层**分档。

分档要从粗到细,先问拓扑再问立体 —— 反过来的话,一条原子数都不对的输出会因为
"骨架比不上"被丢进"其他",于是最严重的一档反而不见了。

  无输出        RDKit 一个能净化的产物都没给出
  原子被复制    产物重原子数 > 输入(共享的部分被搬进了每一片)
  丢原子        产物重原子数 < 输入
  切分不同      原子数对,但分成几片/怎么分不同
  连接不同      原子数对、片数对,去掉立体之后骨架仍不同
  四面体手性    骨架相同,差在手性中心
  双键顺反      骨架相同,差在双键几何
  两者都有      骨架相同,手性与顺反都差

另报一档与命中无关的:**RDKit 的产物净化不过**(`n_bad`),那是拿不出手的输出。
"""

import argparse
import itertools
import json
import os
from collections import Counter

from rdkit import Chem, RDLogger
from rdkit.Chem import rdChemReactions, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
TETRA = (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW)


def canon(s):
    m = Chem.MolFromSmiles(s)
    return Chem.MolToSmiles(m) if m else None


def flat(s):
    m = Chem.MolFromSmiles(s)
    return None if m is None else Chem.MolToSmiles(m, isomericSmiles=False)


def heavy(s):
    m = Chem.MolFromSmiles(s)
    return None if m is None else m.GetNumHeavyAtoms()


def nfrag(s):
    m = Chem.MolFromSmiles(s)
    return None if m is None else len(Chem.GetMolFrags(m))


def formula(s):
    m = Chem.MolFromSmiles(s)
    return None if m is None else rdMolDescriptors.CalcMolFormula(m)


def stereo_kind(truth, pred):
    """骨架相同时,差异落在手性还是顺反。返回集合。"""
    t, p = Chem.MolFromSmiles(truth), Chem.MolFromSmiles(pred)
    if t is None or p is None:
        return {"读不了"}
    t, p = Chem.AddHs(t), Chem.AddHs(p)
    match = p.GetSubstructMatch(t, useChirality=False)
    if not match or len(match) != t.GetNumAtoms():
        return {"配不上"}
    kinds = set()
    for i in range(t.GetNumAtoms()):
        a, b = t.GetAtomWithIdx(i), p.GetAtomWithIdx(match[i])
        ta, tb = a.GetChiralTag(), b.GetChiralTag()
        if (ta in TETRA) != (tb in TETRA):
            kinds.add("手性")
            continue
        if ta in TETRA and tb in TETRA:
            na = [match[n.GetIdx()] for n in a.GetNeighbors()]
            nb = [n.GetIdx() for n in b.GetNeighbors()]
            cur, sw = list(na), 0
            ok = True
            for k, want in enumerate(nb):
                if cur[k] == want:
                    continue
                j = next((j for j in range(k + 1, len(cur)) if cur[j] == want), None)
                if j is None:
                    ok = False
                    break
                cur[k], cur[j] = cur[j], cur[k]
                sw += 1
            if ok and ((ta == tb) != (sw % 2 == 0)):
                kinds.add("手性")
    inv = {match[i]: i for i in range(t.GetNumAtoms())}
    tb_ = {
        frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx())): str(b.GetStereo())
        for b in t.GetBonds()
        if b.GetStereo() != Chem.BondStereo.STEREONONE
    }
    pb = {}
    for b in p.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if b.GetStereo() != Chem.BondStereo.STEREONONE and i in inv and j in inv:
            pb[frozenset((inv[i], inv[j]))] = str(b.GetStereo())
    for k in set(tb_) | set(pb):
        if tb_.get(k) != pb.get(k):
            kinds.add("顺反")
    return kinds or {"查不出"}


def run_rdkit(tpl, inputs):
    rx = rdChemReactions.ReactionFromSmarts(tpl)
    rx.Initialize()
    n = rx.GetNumReactantTemplates()
    mols = [Chem.MolFromSmiles(s) for s in inputs]
    if len(mols) < n:
        return set(), 0
    out, bad = set(), 0
    for perm in itertools.permutations(mols, n):
        for ps in rx.RunReactants(perm):
            frags, ok = [], True
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
            else:
                bad += 1
    return out, bad


def classify(truth, preds, n_in, n_frag_truth):
    if not preds:
        return "无输出", None
    # 先按拓扑挑"最接近"的那个预测
    tf = flat(truth)
    same_skel = [p for p in preds if flat(p) == tf]
    if same_skel:
        kinds = set()
        for p in same_skel:
            kinds |= stereo_kind(truth, p)
        kinds -= {"查不出"}
        if kinds == {"手性"}:
            return "四面体手性", same_skel[0]
        if kinds == {"顺反"}:
            return "双键顺反", same_skel[0]
        if {"手性", "顺反"} <= kinds:
            return "手性+顺反", same_skel[0]
        return "骨架相同但查不出差异", same_skel[0]
    over = [p for p in preds if (heavy(p) or 0) > n_in]
    if over:
        return "原子被复制", max(over, key=lambda s: heavy(s) or 0)
    under = [p for p in preds if (heavy(p) or 0) < n_in]
    if under:
        return "丢原子", min(under, key=lambda s: heavy(s) or 0)
    diff_frag = [p for p in preds if nfrag(p) != n_frag_truth]
    if diff_frag:
        return "切分不同", diff_frag[0]
    ft = formula(truth)
    if any(formula(p) == ft for p in preds):
        return "连接不同(同分异构)", preds[0]
    return "骨架不同", preds[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench", default=os.path.join(ROOT, "results", "bench.jsonl"))
    ap.add_argument("--tpl", default=os.path.join(ROOT, "data", "templates.jsonl"))
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "rdkit_failures.jsonl"))
    args = ap.parse_args()

    tpls = {}
    for line in open(args.tpl):
        r = json.loads(line)
        if "row" in r and "retro" in r:
            tpls[r["row"]] = r

    wins, sani_bad = [], Counter()
    for line in open(args.bench):
        r = json.loads(line)
        if "err" in r:
            continue
        for d in ("fwd", "retro"):
            a, b = r.get(f"omgkit_{d}"), r.get(f"rdkit_{d}")
            if not a or not b or "hit" not in a or "hit" not in b:
                continue
            if b.get("n_bad"):
                sani_bad[("rdkit", d)] += 1
            if a.get("n_bad"):
                sani_bad[("omgkit", d)] += 1
            if a["hit"] and not b["hit"]:
                wins.append((r["row"], d))

    tally, per_dir = Counter(), Counter()
    with open(args.out, "w") as out:
        for row, d in wins:
            t = tpls.get(row)
            if t is None:
                continue
            inputs = t["reactants"] if d == "fwd" else [t["prod"]]
            truth = canon(t["prod"] if d == "fwd" else ".".join(t["reactants"]))
            n_in = sum(heavy(s) or 0 for s in inputs)
            try:
                preds, bad = run_rdkit(t[d], inputs)
            except Exception:
                preds, bad = set(), 0
            kind, closest = classify(truth, sorted(preds), n_in, nfrag(truth))
            tally[kind] += 1
            per_dir[(d, kind)] += 1
            out.write(
                json.dumps(
                    {
                        "row": row,
                        "id": t["id"],
                        "direction": d,
                        "kind": kind,
                        "n_in": n_in,
                        "truth": truth,
                        "rdkit_closest": closest,
                        "rdkit_n_bad_outcomes": bad,
                        "tpl": t[d],
                    }
                )
                + "\n"
            )

    print(f"omgkit 命中而 RDKit 未命中:{len(wins)} 条\n")
    print("| RDKit 错在哪一层 | 条数 | 占比 |")
    print("|---|---|---|")
    for k, v in tally.most_common():
        print(f"| {k} | {v} | {100 * v / len(wins):.1f}% |")
    print("\n按方向:")
    for k, v in sorted(per_dir.items()):
        print(f"  {k[0]:6s} {k[1]:22s} {v}")
    print("\n另计:产物**净化不过**的反应数(与命中无关)")
    for k, v in sorted(sani_bad.items()):
        print(f"  {k[0]:7s} {k[1]:6s} {v}")
    print(f"\n已写出 {args.out}")


if __name__ == "__main__":
    main()
