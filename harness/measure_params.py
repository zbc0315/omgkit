#!/usr/bin/env python3
"""从 MMFF94 优化后的几何里**量**键长与键角,而不是凭记忆写一张表。

    python3 measure_params.py <语料> <输出前缀>

口径(要能重跑,所以每一条都钉死):
- 构象:ETKDGv3,种子 0xf00d,再跑 MMFF94 最小化,**`maxIters=500`,只取返回 0
  (收敛)的分子**。`maxIters` 是承重的:用 RDKit 默认的 200 只收敛 7417 个
  分子,用 2000 收敛 8580 个,用 500 是 8526 个 —— 换个值整张表就对不上。
- 键长键:(较小元素, 较大元素, 键级) + 所在最小环尺寸(不在环里记 0)。
- 键角键:(中心元素, 中心配位数, 是否芳香) + **中心原子自己的最小环尺寸**
  + 三个原子共处的最小环尺寸。
  **中心那一维是后加的**,理由是量出来的:环丙烷的 H–C–H 真值 114.35°,
  而三个原子并不共处一环,于是它落进"共处环 0"那个桶,与无环 sp³ 碳混在一起 ——
  那个桶的 band 是 [106.3, 113.3],**罩不住真值**。judge 于是把"摆对了"读成红、
  把摆成 109.47(错 5°)读成绿。分开之后小环上的环外角才有自己的靶。
- 值按 0.001 Å / 0.1° 舍入后计数 —— 计数表能给出**精确的**中位数与分位数,
  又不用把几十万个浮点写到盘上。

输出两张表:<前缀>.bonds.tsv、<前缀>.angles.tsv,列是 键 计数 中位 p05 p95 均值。
"""

import math
import sys
from collections import Counter, defaultdict
from multiprocessing import Pool

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

ORDER = {Chem.BondType.SINGLE: "1", Chem.BondType.DOUBLE: "2",
         Chem.BondType.TRIPLE: "3", Chem.BondType.AROMATIC: "ar",
         Chem.BondType.DATIVE: "->"}


def min_ring(ri, *idxs):
    """含这些原子的最小环尺寸;不在同一个环里记 0。"""
    best = 0
    for r in ri.AtomRings():
        if all(i in r for i in idxs):
            if best == 0 or len(r) < best:
                best = len(r)
    return best


def one(job):
    _, smi = job
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    mh = Chem.AddHs(m)
    p = AllChem.ETKDGv3()
    p.randomSeed = 0xF00D
    try:
        if AllChem.EmbedMolecule(mh, p) < 0:
            return None
        # 收敛才要:没收敛的几何不是 MMFF 的平衡几何,当参数会把噪声带进来
        if AllChem.MMFFOptimizeMolecule(mh, maxIters=500) != 0:
            return None
    except Exception:  # noqa: BLE001
        return None

    c = mh.GetConformer().GetPositions()
    ri = mh.GetRingInfo()
    bonds, angles = Counter(), Counter()

    for b in mh.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        si, sj = mh.GetAtomWithIdx(i).GetSymbol(), mh.GetAtomWithIdx(j).GetSymbol()
        if si > sj:
            si, sj, i, j = sj, si, j, i
        o = ORDER.get(b.GetBondType(), "?")
        d = float(np.linalg.norm(c[i] - c[j]))
        key = (si, sj, o, min_ring(ri, i, j), round(d, 3))
        bonds[key] = bonds.get(key, 0) + 1

    for a in mh.GetAtoms():
        nb = [x.GetIdx() for x in a.GetNeighbors()]
        if len(nb) < 2:
            continue
        k = a.GetIdx()
        sym, deg, ar = a.GetSymbol(), len(nb), int(a.GetIsAromatic())
        for x in range(len(nb)):
            for y in range(x + 1, len(nb)):
                u, v = c[nb[x]] - c[k], c[nb[y]] - c[k]
                cs = float(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))
                ang = math.degrees(math.acos(max(-1.0, min(1.0, cs))))
                key = (
                    sym,
                    deg,
                    ar,
                    min_ring(ri, k),
                    min_ring(ri, k, nb[x], nb[y]),
                    round(ang, 1),
                )
                angles[key] = angles.get(key, 0) + 1
    return bonds, angles


def summarize(counter, path, valname):
    """把 (键..., 值) -> 次数 的表按键聚合,给出精确分位数。"""
    grp = defaultdict(Counter)
    for key, n in counter.items():
        grp[key[:-1]][key[-1]] += n
    rows = []
    for key, vals in grp.items():
        tot = sum(vals.values())
        xs = sorted(vals)
        run, q = 0, {}
        for x in xs:
            run += vals[x]
            for t in (0.05, 0.5, 0.95):
                if t not in q and run >= tot * t:
                    q[t] = x
        mean = sum(x * n for x, n in vals.items()) / tot
        rows.append((tot, key, q.get(0.5), q.get(0.05), q.get(0.95), mean))
    rows.sort(key=lambda r: -r[0])
    with open(path, "w") as f:
        f.write(f"#键\t计数\t中位{valname}\tp05\tp95\t均值\n")
        for tot, key, med, lo, hi, mean in rows:
            f.write("\t".join(str(x) for x in key) + f"\t{tot}\t{med}\t{lo}\t{hi}\t{mean:.4f}\n")
    return len(rows), sum(r[0] for r in rows)


def main():
    src, pre = sys.argv[1], sys.argv[2]
    jobs = []
    for i, line in enumerate(open(src), 1):
        line = line.strip()
        if line:
            jobs.append((i, line.split("\t")[0]))
    with Pool(6) as pool:
        res = [r for r in pool.map(one, jobs, chunksize=8) if r]
    B, A = Counter(), Counter()
    for b, a in res:
        B.update(b)
        A.update(a)
    nb, tb = summarize(B, pre + ".bonds.tsv", "键长")
    na, ta = summarize(A, pre + ".angles.tsv", "键角")
    print(f"用上 {len(res)}/{len(jobs)} 个分子(嵌入+MMFF 都收敛的)")
    print(f"键:{nb} 种组合,{tb} 条 → {pre}.bonds.tsv")
    print(f"角:{na} 种组合,{ta} 个 → {pre}.angles.tsv")


if __name__ == "__main__":
    main()
