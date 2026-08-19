"""量**环的内扭转角**分布,给 `omgkit-conf` 的界矩阵定 1-4 那一档。

# 为什么要量

界矩阵里,环上一条 1-4 路径的扭转角范围直接决定那一对原子的界有多宽。
拍脑袋给"六元环 0–65°"是能过判据一(真实构象落在界内),但界会松 ——
而松的界正是判据二要拦的。**这两条只能靠实测的分位数同时满足。**

# 口径

- 键长键角表(`measure_params.py`)是在 ETKDGv3 + MMFF94 收敛的结构上量的,
  这里用同一套口径,数才能互相对得上。
- 逐个环,沿环走一圈取连续四元组 `(a,b,c,d)`,记 `|扭转角|`(0~180)。
- 分桶:**环尺寸 + 是否芳环 + 是否全 sp³**。

# 为什么必须加"全 sp³"这一维

只按"尺寸 + 芳香"分桶时,六元非芳那一桶的中位是 **20.9°**,而 p05=0.0、p95=61.2 ——
这个跨度是**两个总体被揉在一起**的signature:共轭近平面环(环己二烯、吡喃酮…)
在 0° 附近,全 sp³ 的椅式在 55° 附近。中位 20.9° 描述的是**两者都不是**的东西。

后果不是"界松一点"。理想椅的闭环条件是 `cos τ = −cos θ / (1 + cos θ)`,
键角 111.6°(同一张表的实测中位)解出 **|τ| = 54.4°**;反过来要让 τ = 20.9°,
键角得是 118.9°。**把这两个中位一起写进界矩阵,得到的是摆不出来的构型** ——
实测端到端:无环分子精修到恰好零,环己烷 1-4 越界 94.1%。

用法:

    python3 harness/measure_ring_torsion.py harness/corpus/large.smi harness/params/mmff.ringtorsion.tsv
"""

import math
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdDistGeom, rdMolTransforms  # noqa: F401

RDLogger.DisableLog("rdApp.*")

SEED = 0xF00D


def quantile(v, f):
    if not v:
        return float("nan")
    i = min(len(v) - 1, max(0, int(round((len(v) - 1) * f))))
    return v[i]


def main() -> int:
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    corpus, out = sys.argv[1], sys.argv[2]
    limit = int(sys.argv[3]) if len(sys.argv) > 3 else 10**9

    buckets: dict[tuple[int, int], list[float]] = {}
    n_ok = n_skip = 0
    for line in open(corpus, encoding="utf-8"):
        smi = line.split("\t")[0].strip()
        if not smi or n_ok >= limit:
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            n_skip += 1
            continue
        mol = Chem.AddHs(mol)
        p = rdDistGeom.ETKDGv3()
        p.randomSeed = SEED
        try:
            if rdDistGeom.EmbedMolecule(mol, p) < 0:
                n_skip += 1
                continue
            # 只收**收敛**的,与 measure_params.py 同口径
            if AllChem.MMFFOptimizeMolecule(mol, maxIters=500) != 0:
                n_skip += 1
                continue
        except Exception:  # noqa: BLE001
            n_skip += 1
            continue
        n_ok += 1
        conf = mol.GetConformer()
        for ring in mol.GetRingInfo().AtomRings():
            k = len(ring)
            arom = int(
                all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring)
            )
            # 全 sp³:环上每个原子都是 sp³ 杂化。这一维把"共轭近平面"与
            # "椅/船"两个总体分开 —— 见文件头。
            sp3 = int(
                all(
                    mol.GetAtomWithIdx(a).GetHybridization()
                    == Chem.HybridizationType.SP3
                    for a in ring
                )
            )
            for t in range(k):
                a, b, c, d = (
                    ring[t],
                    ring[(t + 1) % k],
                    ring[(t + 2) % k],
                    ring[(t + 3) % k],
                )
                try:
                    ang = rdMolTransforms.GetDihedralDeg(conf, a, b, c, d)
                except Exception:  # noqa: BLE001
                    continue
                if math.isfinite(ang):
                    buckets.setdefault((k, arom, sp3), []).append(abs(ang))

    with open(out, "w", encoding="utf-8") as fh:
        fh.write("#环尺寸\t芳香\t全sp3\t计数\t中位\tp05\tp95\t均值\n")
        for (k, arom, sp3), v in sorted(buckets.items()):
            v.sort()
            fh.write(
                f"{k}\t{arom}\t{sp3}\t{len(v)}\t{quantile(v, 0.5):.1f}\t"
                f"{quantile(v, 0.05):.1f}\t{quantile(v, 0.95):.1f}\t"
                f"{sum(v) / len(v):.4f}\n"
            )
    print(f"收敛 {n_ok} 个分子(跳过 {n_skip}),{len(buckets)} 个桶 → {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
