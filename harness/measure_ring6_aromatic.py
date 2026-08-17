#!/usr/bin/env python3
"""按 §42 v2 量孤立芳香六元环的**逐位置角靶**,写成 TSV 给 Rust 那边的生成器用。

口径与 harness/measure_params.py 逐字一致(ETKDGv3 种子 0xf00d、MMFF94
maxIters=500、只取返回 0 的分子),否则与生产那张表不可比。

与 §41.1 那版探针的三处不同,都是送审要求的:
  1. 分档键是**环元素序列对旋转+反射的规范形**(12 种摆法取 min),
     不是排序后的组成 —— 按组成分会把哒嗪/嘧啶/吡嗪混成一桶,跨 6.1°。
  2. 靶取**均值**不是中位数。均值精确可加:Σ(重数×均值) ≡ 均值(逐环内角和),
     所以逐位置靶之和自动满足闭合要求的 720°(差的只是环自身的非平面性)。
  3. 同时输出 p05/p95 —— 权重 w ∝ 1/(p95−p05)² 用它,理由是"对齐 band 判据",
     不是正态假设(实测峰度 1.20~5.74,不正态)。

用法(**必须用 pin 了 numpy<2 的 venv**,系统 python 跑 rdkit 会 SIGSEGV):
    scratchpad/venv/bin/python scratchpad/ring6_targets.py > ring6.targets.tsv
"""

import math
import sys
from collections import defaultdict

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

OK_EL = set("H B C N O F Si P S Cl Se Br I".split())


def canon_ring(seq):
    """环序列对旋转+反射的规范形。返回 (规范串, 原位置 -> 规范位置 的映射)。"""
    n = len(seq)
    best, best_perm = None, None
    for refl in (False, True):
        base = list(reversed(seq)) if refl else list(seq)
        # 原下标在 base 里的位置
        idx = list(reversed(range(n))) if refl else list(range(n))
        for r in range(n):
            cand = tuple(base[(r + k) % n] for k in range(n))
            if best is None or cand < best:
                best = cand
                # base 里第 (r+k)%n 个 -> 规范位置 k
                perm = [0] * n
                for k in range(n):
                    perm[idx[(r + k) % n]] = k
                best_perm = perm
    return "".join(best), best_perm


def one(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None or any(a.GetSymbol() not in OK_EL for a in m.GetAtoms()):
        return None
    mh = Chem.AddHs(m)
    p = AllChem.ETKDGv3()
    p.randomSeed = 0xF00D
    try:
        if AllChem.EmbedMolecule(mh, p) != 0:
            return None
        if AllChem.MMFFOptimizeMolecule(mh, maxIters=500) != 0:
            return None
    except Exception:
        return None
    c = mh.GetConformer().GetPositions()
    ri = mh.GetRingInfo()
    cnt = defaultdict(int)
    for r in ri.AtomRings():
        for x in r:
            cnt[x] += 1
    out = []
    for r in ri.AtomRings():
        if len(r) != 6:
            continue
        if not all(mh.GetAtomWithIdx(x).GetIsAromatic() for x in r):
            continue
        if any(cnt[x] > 1 for x in r):
            continue  # 只要孤立环(一期范围)
        seq = [mh.GetAtomWithIdx(x).GetSymbol() for x in r]
        pat, perm = canon_ring(seq)
        angs = [0.0] * 6
        for pos in range(6):
            k, i, j = r[pos], r[(pos - 1) % 6], r[(pos + 1) % 6]
            pk, p0, p1 = c[k], c[i], c[j]
            u = [p0[t] - pk[t] for t in range(3)]
            v = [p1[t] - pk[t] for t in range(3)]
            du = math.sqrt(sum(t * t for t in u))
            dv = math.sqrt(sum(t * t for t in v))
            cs = sum(u[t] * v[t] for t in range(3)) / (du * dv)
            angs[perm[pos]] = math.degrees(math.acos(max(-1.0, min(1.0, cs))))
        out.append((pat, angs))
    return out


def main():
    corpus = sys.argv[1] if len(sys.argv) > 1 else "harness/corpus/large.smi"
    smis = [l.split("\t")[0].strip() for l in open(corpus) if l.strip()]
    per = defaultdict(lambda: [[] for _ in range(6)])
    rings = defaultdict(int)
    sums = defaultdict(list)
    for i, smi in enumerate(smis, 1):
        r = one(smi)
        if r:
            for pat, angs in r:
                rings[pat] += 1
                sums[pat].append(sum(angs))
                for j in range(6):
                    per[pat][j].append(angs[j])
        if i % 1000 == 0:
            print(f"  ...{i}/{len(smis)}", file=sys.stderr, flush=True)

    # **按 pattern 自身的自同构群把位置并轨。**
    # 自验发现:像吡嗪 CCNCCN 这种对称环,规范形由多个(旋转,反射)同时达到,
    # 于是"原位置 -> 规范位置"本来就不唯一(实测 4 种)。这不是 bug ——
    # 那些位置在化学上**就是等价的**,必须共用一个靶。所以不去纠结 tie-break,
    # 直接把等价位置并轨平均:表因此与 tie-break 怎么选无关,C1 那一侧也就不欠账。
    def orbits(pat):
        n = len(pat)
        auts = []
        for refl in (False, True):
            base = list(reversed(pat)) if refl else list(pat)
            idx = list(reversed(range(n))) if refl else list(range(n))
            for r in range(n):
                if "".join(base[(r + k) % n] for k in range(n)) != pat:
                    continue
                perm = [0] * n
                for k in range(n):
                    perm[idx[(r + k) % n]] = k
                auts.append(perm)
        parent = list(range(n))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for perm in auts:
            for a, b in enumerate(perm):
                ra, rb = find(a), find(b)
                if ra != rb:
                    parent[max(ra, rb)] = min(ra, rb)
        grp = defaultdict(list)
        for x in range(n):
            grp[find(x)].append(x)
        return list(grp.values())

    print("#pattern\t环数\t位置\t均值\tp05\tp95\t中位\t同轨位置")
    for pat in sorted(rings, key=lambda k: -rings[k]):
        n = rings[pat]
        smean = sum(sums[pat]) / n
        print(
            f"# {pat}\tn={n}\tΣ均值={sum(sum(v)/len(v) for v in per[pat]):.4f}"
            f"\t均值(逐环内角和)={smean:.4f}",
        )
        for orb in orbits(pat):
            v = sorted(x for j in orb for x in per[pat][j])
            mean = sum(v) / len(v)
            # 并轨前各位置的均值,拿来自验"它们真的等价"
            spread = max(sum(per[pat][j]) / len(per[pat][j]) for j in orb) - min(
                sum(per[pat][j]) / len(per[pat][j]) for j in orb
            )
            for j in orb:
                print(
                    "%s\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%s(轨内极差 %.4f)"
                    % (pat, n, j, mean, v[int(0.05 * len(v))],
                       v[min(int(0.95 * len(v)), len(v) - 1)], v[len(v) // 2],
                       ",".join(map(str, orb)), spread)
                )


if __name__ == "__main__":
    main()
