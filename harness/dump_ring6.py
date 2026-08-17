#!/usr/bin/env python3
"""把参考构象里**每一个六元环自己的**内坐标与褶皱逐行导出。

`measure_pucker.py` 出的是**统计表**(逐桶的中位/p05/p95)。统计表判不了一件事:
**求解器准不准**。判据看到 θ 偏了,分不清是 [`close_six`] 解错了,还是喂给它的
目标键角本来就不对 —— 这两件事得分开判,不然接线时在 band 上让出去的每一度
都说不清是谁欠的。

这份逐环参考就是拿来分开判的:把**这个分子自己的**六根环内键长与六个环内键角
喂回求解器,要求它复现**这个分子自己的** Q 与 θ。复现得了,说明求解器是准的,
剩下的偏差全记在角模型头上;复现不了,那是求解器的账。

口径**完全继承 `measure_pucker`**(同一批构象:ETKDGv3 种子 0xf00d、嵌 100 个
构象各自 MMFF 极小取最低能;芳香标志在碰 MMFF 之前取)。共用的那几段是 import
过来的,不是抄过来的 —— 抄一份就会分家,分家之后逐环参考与统计表对不上。

    python3 dump_ring.py <语料> <输出.tsv> [种子] [环大小]

`种子` 默认 0xf00d,**入库的表一律用默认值**。给别的值只有一个用途:
量这套协议自己的噪声(同一个分子换个种子,它自己的 Q/θ 会差多少)——
那个差就是逐分子判据容差的下限。种子会写进文件头,换过种子的表一眼认得出来。

`环大小` 默认 6;给 5 就出五元环那一份。收环的条件与统计表一字不差:
非全芳香、杂化全认得出来。

**五元的 `t0..t2` 只有前两个有意义**(五元只有 5 个环内二面角,而这里按
"沿环走"的构造只需要 n−3 = 2 个);第三列对五元恒为 0,留着是为了两种尺寸
共用同一套列。

列(制表符分隔,一行一个环):

| 列 | 意思 |
|---|---|
| `花样` | 旋转+反射规范化的杂化串 + 稠合度,与统计表同一套标签 |
| `杂化` | **按本行的原子次序**的杂化串(没规范化 —— 规范化之后就找不到 sp² 在第几位了) |
| `元素` | 按本行的原子次序的元素串 |
| `双键伙伴` | 逐原子:环**外**双键伙伴的元素,没有记 `-`(逗号分隔) |
| `l0..l5` | 键长 Å,`li` = 第 i 与第 i+1 个原子之间(`l5` 闭合那根) |
| `a0..a5` | 键角度,`ai` = **第 i 个原子处**的环内键角 |
| `t0..t2` | 二面角度,`t0` = (0,1,2,3)、`t1` = (1,2,3,4)、`t2` = (2,3,4,5) |
| `Q` | Cremer–Pople 幅度(Å) |
| `theta` / `phi` | **列名随环大小变**:六元是 θ(折到 [0,90]),五元是相位 φ(折到 [0,18]) |
| `SMILES` | 出处 —— **不记语料行号**,行号会随语料变动静默指错分子 |

`l`/`a`/`t` 的下标约定与 [`omgkit_conformer::ring::close_six`] 的入参**逐位对齐**,
判据那边可以原样喂进去。
"""
import math, sys
from multiprocessing import Pool
import numpy as np
from rdkit import Chem, RDLogger
from measure_pucker import (
    HYB,
    NCONF,
    best_conformer,
    canon_pattern,
    cremer_pople,
    fused_tag,
    judgeable_ring,
)
RDLogger.DisableLog("rdApp.*")

def angle_deg(a, b, c):
    """b 处的键角(度)"""
    u, v = a - b, c - b
    cos = float(u @ v / (np.linalg.norm(u) * np.linalg.norm(v)))
    return math.degrees(math.acos(max(-1.0, min(1.0, cos))))

def dihedral_deg(a, b, c, d):
    """二面角 (a,b,c,d)(度,IUPAC 符号)"""
    b1, b2, b3 = b - a, c - b, d - c
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    m = np.cross(n1, b2 / np.linalg.norm(b2))
    return math.degrees(math.atan2(float(m @ n2), float(n1 @ n2)))

def dbl_partner(mh, idx, ring):
    """`idx` 的**环外**双键伙伴元素;没有返回 `-`"""
    for b in mh.GetAtomWithIdx(idx).GetBonds():
        if b.GetBondType() != Chem.BondType.DOUBLE:
            continue
        o = b.GetOtherAtomIdx(idx)
        if o not in ring:
            return mh.GetAtomWithIdx(o).GetSymbol()
    return "-"

def one(job):
    # **种子随 job 传,不走模块级全局。** macOS 上 `Pool` 默认是 spawn:
    # 子进程重新 import 本模块,`__main__` 里对全局的赋值**到不了子进程** ——
    # 那样跑出来的表会写着 0xbeef 的文件头、装着 0xf00d 的数据,
    # 而这种表恰恰是拿来量"换种子差多少"的,静默同种子等于量了个寂寞。
    smi, seed, size = job
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return []
    mh = Chem.AddHs(m)
    if not judgeable_ring(mh):
        return []
    got = best_conformer(mh, seed)
    if got is None:
        return []
    c, arom = got
    rows = []
    allr = mh.GetRingInfo().AtomRings()
    for ring in allr:
        if len(ring) != size:
            continue
        sym = "".join(HYB.get(mh.GetAtomWithIdx(i).GetHybridization(), "0") for i in ring)
        if "0" in sym:
            continue
        if all(arom[i] for i in ring):
            continue
        # **环的原子次序要自己验一遍。** 下面所有的键长/键角/二面角都压在
        # "相邻两个原子成键"上;RDKit 的 `AtomRings()` 一向是按环序给的,
        # 但这是别人家的实现细节,验一次不花钱,验不过就如实丢掉这个环。
        if any(
            mh.GetBondBetweenAtoms(ring[k], ring[(k + 1) % size]) is None
            for k in range(size)
        ):
            continue
        cp = cremer_pople(np.array([c[i] for i in ring]))
        if cp is None:
            continue
        Q, theta = cp
        p = [c[i] for i in ring]
        lens = [float(np.linalg.norm(p[(k + 1) % size] - p[k])) for k in range(size)]
        angs = [
            angle_deg(p[(k + size - 1) % size], p[k], p[(k + 1) % size])
            for k in range(size)
        ]
        # 沿环走需要 n−3 个二面角;为了两种尺寸共用列数,不足的补 0
        tors = [dihedral_deg(p[k], p[k + 1], p[k + 2], p[k + 3]) for k in range(size - 3)]
        tors += [0.0] * (3 - len(tors))
        rows.append(
            "\t".join(
                [
                    canon_pattern(sym) + fused_tag(ring, allr),
                    sym,
                    "".join(mh.GetAtomWithIdx(i).GetSymbol() for i in ring),
                    ",".join(dbl_partner(mh, i, ring) for i in ring),
                    # **六位小数不是为了好看。** 判据要把这些值原样喂回求解器,
                    # 再要求复现同一个 θ —— 喂进去的量一旦被舍入,解跟着挪,
                    # 挪多少取决于这个环有多病态,而那是环的性质、不是求解器的账。
                    # θ 对键角的敏感度实测约 1°/1°,舍到 1e-6 度就把这项压到
                    # 判据容差(0.1°)的十万分之一以下,红了就一定是求解器的事。
                    *(f"{v:.6f}" for v in lens),
                    *(f"{v:.6f}" for v in angs),
                    *(f"{v:.6f}" for v in tors),
                    f"{Q:.6f}",
                    f"{theta:.4f}",
                    smi,
                ]
            )
        )
    return rows

if __name__ == "__main__":
    corpus, out = sys.argv[1], sys.argv[2]
    seed = int(sys.argv[3], 0) if len(sys.argv) > 3 else 0xF00D
    size = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    # **同一个 SMILES 只跑一遍。** 语料 8863 行里只有 8749 个不同的 SMILES;
    # 逐行跑会让同一个分子的同一个环在表里出现两遍,而判据是按 SMILES 查表的 ——
    # 它会看到两个候选、判成"签名撞车"、退回池化的桶 band。实测这样白白放过了
    # 一个 0.457 Å 的真偏差。参考表的键是 SMILES,一个键就该只有一份答案。
    seen = set()
    smis = []
    for l in open(corpus):
        smi = l.split("\t")[0].strip()
        if smi and smi not in seen:
            seen.add(smi)
            smis.append(smi)
    rows = []
    with Pool() as pool:
        for r in pool.imap_unordered(one, [(s, seed, size) for s in smis], chunksize=8):
            rows.extend(r)
    # 排序后再写:多进程的完成次序是不定的,不排序每次重跑都是一份新 diff
    rows.sort()
    head = (
        ["#花样", "杂化", "元素", "双键伙伴"]
        + [f"l{k}" for k in range(size)]
        + [f"a{k}" for k in range(size)]
        + [f"t{k}" for k in range(3)]
        # **列名必须随环大小变。** 这一列六元装 θ、五元装 φ,是两个不同的量
        # (值域 [0,90] 与 [0,18],一个是折叠角、一个是相位)。先前无条件写
        # `theta`,判据那边跟着把五元的 φ 当 θ 用,混进了折成 Å 的褶皱距离里 ——
        # 一条实打实的错判,查了两轮才逮到。`measure_pucker.py` 一直是分开写的。
        + ["Q", "theta" if size == 6 else "phi", "SMILES"]
    )
    with open(out, "w") as f:
        f.write(f"#种子 {seed:#x}\tNCONF {NCONF}\t环大小 {size}\n")
        f.write("\t".join(head) + "\n")
        for r in rows:
            f.write(r + "\n")
    print(f"{len(rows)} 个 {size} 元环 → {out}")
