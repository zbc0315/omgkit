"""**这个脚本定的是要打的目标线。** 全语料 8833 个分子的实测(2026-08-18):

| 配置 | 失败 | 总耗时 | 中位 | 最慢 |
|---|---|---|---|---|
| A 嵌入 + 知识项(ETKDGv3 默认) | 46(0.52%) | 93.1 s | 3.6 ms | 8.62 s |
| B 跳过嵌入(随机坐标)+ 知识项 | 33(0.37%) | 251.5 s | 4.7 ms | 35.5 s |
| **C 嵌入 + 纯 DG(关掉 ET/K)** | **40(0.45%)** | **37.4 s** | **1.0 ms** | 5.67 s |
| D 跳过嵌入 + 纯 DG | 27(0.31%) | 135.0 s | 2.2 ms | 23.2 s |

读出来两件事:

1. **跳过嵌入慢 2.7 倍,不是快。** 度量矩阵嵌入不是成本 —— 它给优化器一个好起点,
   省下的迭代远超那次特征分解。"用别的办法造初值来砍掉 O(N³)"这条路
   (本仓一度打算走的)**就是被这一行否掉的**。
2. **ET/K 那层买的是"像真的",不是"成功率"。** 关掉它快 2.5 倍,失败反而少 6 个 ——
   实测扭转分布、平面环、反式酰胺正是后续力场自己会做的事。
   本仓要的是"交给力场优化的起点",所以目标线取 **C**。

量三件事:嵌入这一步有多吃劲、ET/K 那些知识项有多吃劲、代价各是多少。

A  ETKDGv3 默认(特征值嵌入 + ET + K)      —— 基线
B  ETKDGv3 + useRandomCoords=True          —— **跳过嵌入**,从随机坐标精修
C  纯 DG(关掉 ET 与 K)                     —— 知识项买了多少
D  纯 DG + useRandomCoords                  —— 两样都不要
全部固定 randomSeed=0xf00d,单进程。
"""
import sys, time, statistics
from rdkit import Chem, RDLogger
from rdkit.Chem import rdDistGeom as dg
RDLogger.DisableLog('rdApp.*')

corpus = sys.argv[1]
smis = []
for line in open(corpus, encoding='utf-8'):
    s = line.split('\t')[0].strip()
    if s: smis.append(s)

def mkparams(random_coords, knowledge):
    p = dg.ETKDGv3() if knowledge else dg.EmbedParameters()
    p.randomSeed = 0xf00d
    p.useRandomCoords = random_coords
    if not knowledge:
        p.useExpTorsionAnglePrefs = False
        p.useBasicKnowledge = False
        p.enforceChirality = True
    return p

for name, rc, kn in [("A 嵌入+知识(ETKDGv3)", False, True),
                     ("B 随机坐标+知识",       True,  True),
                     ("C 嵌入+纯DG",           False, False),
                     ("D 随机坐标+纯DG",       True,  False)]:
    ts, fail, n, worst, worst_smi = [], 0, 0, 0.0, ""
    t_all = time.perf_counter()
    for s in smis:
        m = Chem.MolFromSmiles(s)
        if m is None: continue
        m = Chem.AddHs(m)
        n += 1
        t0 = time.perf_counter()
        try:
            cid = dg.EmbedMolecule(m, mkparams(rc, kn))
        except Exception:
            cid = -1
        dt = time.perf_counter() - t0
        ts.append(dt)
        if dt > worst: worst, worst_smi = dt, s
        if cid < 0: fail += 1
    tot = time.perf_counter() - t_all
    ts.sort()
    print(f"{name}: 分子 {n}  失败 {fail}({100*fail/n:.2f}%)  "
          f"总 {tot:.1f}s  均 {1000*statistics.mean(ts):.1f}ms  "
          f"中位 {1000*statistics.median(ts):.1f}ms  最慢 {worst:.2f}s")
    print(f"    最慢那个:{worst_smi[:70]}")
    sys.stdout.flush()
