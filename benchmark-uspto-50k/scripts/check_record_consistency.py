"""记录自己前后一致吗:反应物侧与产物侧在同一个映射号上说的是同一个构型吗。

# 为什么不能用 CIP 判

CIP 编号是按取代基**优先级**排的,而反应会换掉取代基 —— 酰胺氮换成胺氮,优先
级次序可能整个翻个个儿。同一个空间构型,反应前后的 CIP 字母完全可以不同。拿
CIP 跨反应比对,报出来的"翻转"一大半是假的。

# 用什么判

四面体标记本身相对**邻居的存储顺序**。把两侧的邻居都换算成映射号,就得到两个
可比的序列;两个序列之间的置换宇称,加上两侧的标记,合起来才是"空间构型变没变":

    构型相同  ⟺  (标记相同) == (置换为偶)

这与判断模板两侧次序是否对调用的是同一条式子。邻居里的隐式氢按"排在最前"处理
(SMILES 的约定:方括号里的氢算作第一个邻居),两侧口径一致即可。

# 输出

对每条反应,列出所有**两侧都带标记且邻居映射号集合相同**的中心,报它们构型是否
一致。邻居集合不同的中心(取代基被换掉了)单列一档 —— 那种情形靠这条式子判不了,
也不该硬判。
"""

import argparse
import json
import os
import re
from collections import Counter

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def neighbour_maps(atom):
    """邻居的映射号序列,按键的存储顺序;隐式氢记作 0 并排在最前。

    方括号里的隐式氢在 SMILES 里算第一个邻居,两侧同一套口径即可对比。
    """
    seq = []
    if atom.GetTotalNumHs() == 1:
        seq.append(0)
    seq.extend(nb.GetAtomMapNum() for nb in atom.GetNeighbors())
    return seq


def permutation_is_odd(a, b):
    """把 a 换成 b 需要几次对换;元素对不上时返回 None。"""
    if sorted(a) != sorted(b) or len(set(a)) != len(a):
        return None
    cur = list(a)
    swaps = 0
    for i, want in enumerate(b):
        if cur[i] == want:
            continue
        j = next((j for j in range(i + 1, len(cur)) if cur[j] == want), None)
        if j is None:
            return None
        cur[i], cur[j] = cur[j], cur[i]
        swaps += 1
    return swaps % 2 == 1


def tagged_centres(mol):
    out = {}
    for a in mol.GetAtoms():
        if a.GetChiralTag() in (
            Chem.ChiralType.CHI_TETRAHEDRAL_CW,
            Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
        ) and a.GetAtomMapNum():
            out[a.GetAtomMapNum()] = (a.GetChiralTag(), neighbour_maps(a))
    return out


def compare(rxn_smiles):
    """返回 (一致数, 翻转数, 判不了数, 翻转的映射号列表)。"""
    lhs, _, rhs = rxn_smiles.split(">")
    rm, pm = Chem.MolFromSmiles(lhs), Chem.MolFromSmiles(rhs)
    if rm is None or pm is None:
        return None
    r, p = tagged_centres(rm), tagged_centres(pm)
    same = flip = unknown = 0
    flipped = []
    for m in sorted(set(r) & set(p)):
        (tr, nr), (tp, np_) = r[m], p[m]
        odd = permutation_is_odd(nr, np_)
        if odd is None:
            unknown += 1
            continue
        if (tr == tp) == (not odd):
            same += 1
        else:
            flip += 1
            flipped.append(m)
    return same, flip, unknown, flipped


def template_maps_with_chirality(smarts):
    """模板里**写了手性**的映射号集合。"""
    out = set()
    for tok in re.findall(r"\[[^\]]*\]", smarts):
        if "@" not in tok:
            continue
        m = re.search(r":(\d+)\]", tok)
        if m:
            out.add(int(m.group(1)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--misses", default=os.path.join(ROOT, "results", "misses_attributed.jsonl"))
    ap.add_argument("--raw", default=os.path.join(ROOT, "data", "uspto50k_raw.csv"))
    ap.add_argument("--bucket", default="record-stereo")
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "record_consistency.jsonl"))
    args = ap.parse_args()

    import csv

    # 按**枚举下标**索引,与 extract_templates.py 写进 templates.jsonl 的 row 同源。
    # 不能用 CSV 第一列 —— 它不是行号(50016 行里有 49532 行对不上),
    # 拿它当键会静默取到另一条反应,而后面每一步都照样跑得通。
    raw = {}
    with open(args.raw) as fh:
        for i, row in enumerate(csv.DictReader(fh)):
            raw[i] = row["rxn_smiles"]

    tally = Counter()
    with open(args.misses) as fh, open(args.out, "w") as out:
        for line in fh:
            rec = json.loads(line)
            if args.bucket != "all" and rec.get("bucket") != args.bucket:
                continue
            rxn = raw.get(rec["row"])
            if rxn is None:
                tally["无原始记录"] += 1
                continue
            res = compare(rxn)
            if res is None:
                tally["读不了"] += 1
                continue
            same, flip, unknown, flipped = res
            tpl = rec["retro"] if rec["direction"] == "retro" else rec["fwd"]
            chiral_in_tpl = template_maps_with_chirality(tpl)
            covered = [m for m in flipped if m in chiral_in_tpl]
            if flip == 0:
                key = "记录前后一致" if (same or unknown) else "两侧都没有可比中心"
            elif covered:
                key = "记录说翻了,模板也写了手性"
            else:
                key = "记录说翻了,模板没写手性"
            tally[key] += 1
            out.write(
                json.dumps(
                    {
                        "id": rec["id"],
                        "row": rec["row"],
                        "direction": rec["direction"],
                        "verdict": key,
                        "same": same,
                        "flip": flip,
                        "unknown": unknown,
                        "flipped_maps": flipped,
                        "template_chiral_maps": sorted(chiral_in_tpl),
                        "truth": rec["truth"],
                        "closest": rec.get("closest"),
                    }
                )
                + "\n"
            )

    print(f"== 档次 {args.bucket} 的记录自洽性 ==")
    for k, v in tally.most_common():
        print(f"  {k:26s} {v}")
    print(f"\n合计 {sum(tally.values())},已写入 {args.out}")


if __name__ == "__main__":
    main()
