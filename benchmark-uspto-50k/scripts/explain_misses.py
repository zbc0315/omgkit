"""对每一条 omgkit 未命中给出**逐个立体中心**的定责。

# 为什么必须逐中心定责

一条反应里可以同时有好几处立体分歧,成因还各不相同:一个中心是记录自己前后
矛盾,另一个中心是模板写了而实现没落实。按**整条**归责,只要有一处沾了"模板
写了",整条就被算到引擎头上 —— 实测这样会把 126 条记录问题误判成引擎问题。

所以这里把 omgkit 的输出与**带映射号的记录**逐原子配起来,每个有分歧的中心
单独定责,最后再汇总。

# 每个中心怎么定责

设输入侧(正向是反应物、逆向是产物)为 src,目标侧为 dst。

  记录·造不出   dst 有、src 没有、模板也没写    → 谁也造不出这个构型
  记录·漏写     src 有、dst 没有                → 记录漏写,保住它反而更全
  记录·自相矛盾 两侧都有但构型不同(按映射号邻居序的宇称判)
  引擎·丢了     src 有、dst 有、构型一致,而 omgkit 的输出没写
  引擎·翻了     同上,但 omgkit 写反了
  引擎·模板没落实 模板在这个映射号上写了手性,omgkit 的输出与 dst 不符

构型比对一律不用 CIP 跨反应比 —— 取代基一换,同一个空间构型的 CIP 字母就可能
变。跨反应用**映射号邻居序的置换宇称**;而 omgkit 的输出与 dst **构成相同**,
那一侧才可以放心用 CIP。
"""

import argparse
import csv
import json
import os
import re
from collections import Counter

from rdkit import Chem, RDLogger
from rdkit.Chem import rdChemReactions

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

TETRA = (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW)


# ---------------------------------------------------------------- 基本工具


def with_hs(mol):
    """把氢补成真原子。

    四面体标记相对**邻居的存储顺序**,而隐式氢在这个顺序里的位置取决于它在
    原串里写在哪 —— 靠"氢一律排最前"去猜,同一个构型换种写法就会判反。补成
    真原子之后,顺序里再没有猜的成分。
    """
    return Chem.AddHs(mol) if mol is not None else None


def neighbour_maps(atom):
    """邻居的映射号序列,按键的存储顺序。补过氢的分子上,氢的映射号是 0。

    四面体中心至多带一个氢,所以序列里至多一个 0,不会有歧义。
    """
    return [nb.GetAtomMapNum() for nb in atom.GetNeighbors()]


def permutation_is_odd(a, b):
    if sorted(a) != sorted(b) or len(set(a)) != len(a):
        return None
    cur, swaps = list(a), 0
    for i, want in enumerate(b):
        if cur[i] == want:
            continue
        j = next((j for j in range(i + 1, len(cur)) if cur[j] == want), None)
        if j is None:
            return None
        cur[i], cur[j] = cur[j], cur[i]
        swaps += 1
    return swaps % 2 == 1


def centres(mol):
    return {
        a.GetAtomMapNum(): (a.GetChiralTag(), neighbour_maps(a))
        for a in mol.GetAtoms()
        if a.GetChiralTag() in TETRA and a.GetAtomMapNum()
    }


def stereo_bonds(mol):
    out = {}
    for b in mol.GetBonds():
        if b.GetStereo() == Chem.BondStereo.STEREONONE:
            continue
        i, j = b.GetBeginAtom(), b.GetEndAtom()
        if not (i.GetAtomMapNum() and j.GetAtomMapNum()):
            continue
        refs = tuple(mol.GetAtomWithIdx(x).GetAtomMapNum() for x in b.GetStereoAtoms())
        out[frozenset((i.GetAtomMapNum(), j.GetAtomMapNum()))] = (str(b.GetStereo()), refs)
    return out


def template_chiral_maps(smarts):
    out = set()
    for tok in re.findall(r"\[[^\]]*\]", smarts):
        if "@" in tok:
            m = re.search(r":(\d+)\]", tok)
            if m:
                out.add(int(m.group(1)))
    return out


def product_determined_bonds(smarts):
    """产物模板里几何被**写全**了的双键(两端各有一根方向键),按映射号对给出。

    判据与 Rust 侧的 `honoured_directions` 一致。只看"串里有没有斜杠"是错的 ——
    孤零零一根 `/` 定不了任何几何(`F/C=CF` 就是)。
    """
    try:
        rxn = rdChemReactions.ReactionFromSmarts(smarts)
    except Exception:
        return set()
    out = set()
    for m in rxn.GetProducts():
        dirs = {b.GetIdx() for b in m.GetBonds() if b.GetBondDir() != Chem.BondDir.NONE}
        for b in m.GetBonds():
            if b.GetBondType() != Chem.BondType.DOUBLE:
                continue
            flanked = all(
                any(x.GetIdx() in dirs for x in end.GetBonds() if x.GetIdx() != b.GetIdx())
                for end in (b.GetBeginAtom(), b.GetEndAtom())
            )
            ma, mb = b.GetBeginAtom().GetAtomMapNum(), b.GetEndAtom().GetAtomMapNum()
            if flanked and ma and mb:
                out.add(frozenset((ma, mb)))
    return out


def contributing(lhs, rhs):
    """反应物侧里**参与反应**的片段(有映射号出现在产物侧),拼回一个串。"""
    pm = Chem.MolFromSmiles(rhs)
    if pm is None:
        return None
    pmaps = {a.GetAtomMapNum() for a in pm.GetAtoms() if a.GetAtomMapNum()}
    keep = []
    for frag in lhs.split("."):
        fm = Chem.MolFromSmiles(frag)
        if fm is not None and any(a.GetAtomMapNum() in pmaps for a in fm.GetAtoms()):
            keep.append(frag)
    return ".".join(keep) if keep else None


# ---------------------------------------------------------------- 逐中心定责


def analyse(rec, rxn_smiles):
    lhs, _, rhs = rxn_smiles.split(">")
    contrib = contributing(lhs, rhs)
    if contrib is None:
        return "算不出参与反应的片段", []
    src_s, dst_s = (contrib, rhs) if rec["direction"] == "fwd" else (rhs, contrib)
    src = with_hs(Chem.MolFromSmiles(src_s))
    dst = with_hs(Chem.MolFromSmiles(dst_s))
    pred = with_hs(Chem.MolFromSmiles(rec["closest"])) if rec.get("closest") else None
    if src is None or dst is None or pred is None:
        return "读不了", []

    # dst(带映射号)与 pred(无映射号)逐原子配对;构成相同才配得上
    match = pred.GetSubstructMatch(dst, useChirality=False)
    if not match or len(match) != dst.GetNumAtoms():
        return "配不上", []

    tpl = rec["retro"] if rec["direction"] == "retro" else rec["fwd"]
    tpl_chiral = template_chiral_maps(tpl)
    tpl_bonds = product_determined_bonds(tpl)
    sc, dc = centres(src), centres(dst)
    sb, db = stereo_bonds(src), stereo_bonds(dst)

    verdicts = []

    def same_config(i):
        """dst 第 i 个原子与 pred 里对应原子的构型:True/False,标记缺一边时 None。

        不用 CIP —— CIP 按取代基优先级排,同一分子里相邻中心一起变时会连带
        翻个个儿,于是一处真差异会虚报成两处。这里直接比标记加邻居序的宇称,
        与参照系无关。
        """
        a, b = dst.GetAtomWithIdx(i), pred.GetAtomWithIdx(match[i])
        ta, tb = a.GetChiralTag(), b.GetChiralTag()
        if ta not in TETRA and tb not in TETRA:
            return True
        if (ta in TETRA) != (tb in TETRA):
            return None
        na = [match[n.GetIdx()] for n in a.GetNeighbors()]
        nb = [n.GetIdx() for n in b.GetNeighbors()]
        odd = permutation_is_odd(na, nb)
        if odd is None:
            return None
        return (ta == tb) == (not odd)

    # ---- 四面体中心 ----
    for i in range(dst.GetNumAtoms()):
        m = dst.GetAtomWithIdx(i).GetAtomMapNum()
        eq = same_config(i)
        if eq is True:
            continue
        tp = pred.GetAtomWithIdx(match[i]).GetChiralTag() in TETRA or None
        if not m:
            verdicts.append(("记录·无映射号的中心对不上", f"atom{i}"))
            continue
        if m in dc and m not in sc:
            verdicts.append(
                ("引擎·模板写了却没落实", m)
                if m in tpl_chiral
                else ("记录·目标要的构型输入侧没有", m)
            )
        elif m in sc and m not in dc:
            verdicts.append(("记录·目标侧漏写", m))
        elif m in sc and m in dc:
            odd = permutation_is_odd(sc[m][1], dc[m][1])
            if odd is None:
                verdicts.append(("记录·取代基换了,构型比不了", m))
            elif (sc[m][0] == dc[m][0]) != (not odd):
                verdicts.append(("记录·自相矛盾", m))
            elif tp is None:
                verdicts.append(("引擎·丢了", m))
            else:
                verdicts.append(("引擎·翻了", m))
        else:
            verdicts.append(("记录·两侧都没标这个中心", m))

    # ---- 双键顺反 ----
    pred_bs = {}
    inv = {match[i]: i for i in range(dst.GetNumAtoms())}
    for b in pred.GetBonds():
        if b.GetStereo() == Chem.BondStereo.STEREONONE:
            continue
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in inv and j in inv:
            pred_bs[frozenset((inv[i], inv[j]))] = str(b.GetStereo())
    dst_bs = {
        frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx())): str(b.GetStereo())
        for b in dst.GetBonds()
        if b.GetStereo() != Chem.BondStereo.STEREONONE
    }
    for k in set(dst_bs) | set(pred_bs):
        if dst_bs.get(k) == pred_bs.get(k):
            continue
        maps = frozenset(dst.GetAtomWithIdx(i).GetAtomMapNum() for i in k)
        if 0 in maps:
            verdicts.append(("记录·无映射号的双键对不上", sorted(maps)))
            continue
        if maps in db and maps not in sb:
            verdicts.append(
                ("引擎·模板写了却没落实", sorted(maps))
                if maps in tpl_bonds
                else ("记录·目标要的几何输入侧没有", sorted(maps))
            )
        elif maps in sb and maps not in db:
            verdicts.append(("记录·目标侧漏写", sorted(maps)))
        elif maps in sb and maps in db:
            if set(sb[maps][1]) != set(db[maps][1]):
                verdicts.append(("记录·参照原子换了,几何比不了", sorted(maps)))
            elif sb[maps][0] != db[maps][0]:
                verdicts.append(("记录·自相矛盾", sorted(maps)))
            elif k not in pred_bs:
                verdicts.append(("引擎·丢了", sorted(maps)))
            else:
                verdicts.append(("引擎·翻了", sorted(maps)))
        else:
            verdicts.append(("记录·两侧都没标这根双键", sorted(maps)))
    return None, verdicts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default=os.path.join(ROOT, "results", "stereo_diff.jsonl"))
    ap.add_argument("--raw", default=os.path.join(ROOT, "data", "uspto50k_raw.csv"))
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "miss_verdicts.jsonl"))
    ap.add_argument("--show", type=int, default=3)
    args = ap.parse_args()

    # 按**枚举下标**索引,与 templates.jsonl 的 row 同源。CSV 第一列不是行号
    # (50016 行里 49532 行对不上),拿它当键会静默取到另一条反应。
    raw = {}
    with open(args.raw) as fh:
        for i, row in enumerate(csv.DictReader(fh)):
            raw[i] = row["rxn_smiles"]

    by_centre = Counter()
    by_record = Counter()
    samples = {}
    with open(args.src) as fh, open(args.out, "w") as out:
        for line in fh:
            rec = json.loads(line)
            rxn = raw.get(rec["row"])
            if rxn is None:
                by_record["无原始记录"] += 1
                continue
            err, vs = analyse(rec, rxn)
            if err:
                by_record[err] += 1
                continue
            if not vs:
                by_record["查不出差异"] += 1
                continue
            for kind, _ in vs:
                by_centre[kind] += 1
            # 整条的结论:只要有一处是引擎的问题,这条就得修
            eng = [k for k, _ in vs if k.startswith("引擎")]
            label = eng[0] if eng else vs[0][0]
            by_record[label] += 1
            samples.setdefault(label, []).append((rec, vs))
            out.write(json.dumps({**rec, "verdict": label, "centres": vs}) + "\n")

    print("== 按**立体中心**计 ==")
    for k, v in by_centre.most_common():
        print(f"  {k:30s} {v}")
    print(f"  合计 {sum(by_centre.values())} 个中心")
    print("\n== 按**反应**计(有一处归引擎就算引擎)==")
    for k, v in by_record.most_common():
        print(f"  {k:30s} {v}")
    print(f"  合计 {sum(by_record.values())} 条")

    for label, rows in samples.items():
        if not label.startswith("引擎"):
            continue
        print(f"\n---- {label}({len(rows)} 条)----")
        for rec, vs in rows[: args.show]:
            print(f"  row={rec['row']} {rec['id']} [{rec['direction']}] {vs}")
            print(f"    真值   {rec['truth']}")
            print(f"    omgkit {rec['closest']}")
    print(f"\n已写入 {args.out}")


if __name__ == "__main__":
    main()
