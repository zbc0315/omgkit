"""逐键算键长偏差 —— 与 `measure_params.py` 的"逐组合看中位"是**不同的量**。

    python3 harness/measure_bond_error.py

# 为什么要单独有这一个

`measure_params.py` 按"元素对 × 键级 × 最小环尺寸"分档,每档给中位数。拿那些
中位数算模型偏差、再按键数加权,得到的是**"落在偏差大的档里的键占多少"**,它
把档**内部**的离散度整个忽略掉了,系统性偏小。

实测同一批数据(8526 个分子、240703 条键):

| | 逐档看中位 | **逐键** |
|---|---:|---:|
| × 键级系数,平均误差 | 1.43% | **1.54%** |
| × 键级系数,>5% 的键 | 2.8% | **5.24%** |

差 1.9 倍。`omgkit-conformer::params` 的模块文档第一版登的是前者,是错的。
判据要的是后者 —— 一根键摆得对不对,与它属于哪个统计档无关。
"""
import sys
from collections import Counter
from multiprocessing import Pool
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog("rdApp.*")
ORDER = {Chem.BondType.SINGLE:"1", Chem.BondType.DOUBLE:"2", Chem.BondType.TRIPLE:"3",
         Chem.BondType.AROMATIC:"ar", Chem.BondType.DATIVE:"->"}
F = {"1":1.0070, "ar":0.9206, "2":0.8716, "3":0.7892}
pt = Chem.GetPeriodicTable()

def one(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None: return None
    mh = Chem.AddHs(m)
    p = AllChem.ETKDGv3(); p.randomSeed = 0xF00D
    try:
        if AllChem.EmbedMolecule(mh, p) < 0: return None
        if AllChem.MMFFOptimizeMolecule(mh, maxIters=500) != 0: return None
    except Exception: return None
    c = mh.GetConformer().GetPositions()
    out = []
    for b in mh.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        si = mh.GetAtomWithIdx(i).GetSymbol(); sj = mh.GetAtomWithIdx(j).GetSymbol()
        o = ORDER.get(b.GetBondType(), "?")
        if o == "?": continue
        d = float(np.linalg.norm(c[i]-c[j]))
        raw = pt.GetRcovalent(si) + pt.GetRcovalent(sj)
        out.append((abs(raw-d)/d, abs(raw*F.get(o,1.0)-d)/d))
    return out

if __name__ == "__main__":
    smis = [l.strip().split("\t")[0] for l in open("/Users/tom/Projects/momega/omgkit/harness/corpus/large.smi") if l.strip()]
    with Pool(6) as pool:
        res = [r for r in pool.map(one, smis, chunksize=8) if r]
    raw = np.array([x[0] for r in res for x in r]); fix = np.array([x[1] for r in res for x in r])
    print(f"用上 {len(res)} 个分子,{len(raw)} 条键(maxIters=500)")
    for name, e in (("rcov 和", raw), ("× 键级系数", fix)):
        print(f"  {name:<12} 加权平均 {100*e.mean():5.2f}%   >5% {100*(e>0.05).mean():5.2f}% ({int((e>0.05).sum())})   >10% {100*(e>0.10).mean():5.2f}% ({int((e>0.10).sum())})")
