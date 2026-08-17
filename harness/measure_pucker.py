#!/usr/bin/env python3
"""从 MMFF94 优化后的几何里量**环的褶皱**。

**这里的口径与 measure_params.py 有一处故意不同,理由是实测出来的。**

键长表与键角表用的是"ETKDGv3 单次嵌入(种子 0xf00d)+ MMFF94 局部极小"。
那对键长键角没问题(椅式与扭船的键长几乎一样、键角差不到 1°),
**但对褶皱是错的** —— MMFF 是局部极小,跨不过椅式↔扭船那道势垒,
于是拿到的是 ETKDG 随机落点的抽签结果。实测:

    C1(N2CCCCC2)=C3C(C=CC=C3)=CC=N1   单次 87.1°(扭船) → 30 次取最低能 4.5°(椅式)
    C1(=C(C(N)=NC(=N1)Cl)[N+]([O-])=O)N2CCCCC2  单次 87.3°(扭船) → 30 次取最低能 4.9°(椅式)
    O=C1CCCCC1(非语料分子,手写对照) 单次  6.5°       → 30 次取最低能 6.5°(本来就是椅式)

按单次口径量,孤立的 233333 六元环里 58.7% 落在 θ>60°(船样),折后中位 85.6° ——
而那不是化学,是抽签。用这种分布做 [p05,p95] band,椅式与扭船会一起被放过,
判据等于没有。

所以褶皱表**多构象取最低能**:每个分子嵌 N 个构象、各自 MMFF 极小、取能量最低的
那一个。种子仍钉 0xf00d,构象数钉死,可重跑。


    python3 measure_pucker.py <语料> <输出前缀>

口径(**故意与键长/键角表不同**,理由见上;其余各条钉死以便重跑):
- 构象:ETKDGv3,种子 0xf00d,再跑 MMFF94 最小化,maxIters=500,只取收敛的分子。
- 环:RDKit 的 SSSR(GetRingInfo().AtomRings()),只取 5 元与 6 元。
- 花样:环上每个原子的杂化 2/3(SP2 记 2、SP3 记 3、其余记 0),按环序取,
  再按**旋转 + 反射**规范化成字典序最小的那个串 —— 于是 233333 与 333323 同一桶。
- 量:Cremer–Pople 褶皱幅度 Q(Å)。6 元再记 θ(度,0=椅式、90=船/扭船);
  5 元记幅度即可(相位 φ 与花样的对齐方式有关,先不进键)。
- 值按 0.001 Å / 0.1° 舍入后计数,分位数由计数表精确给出。

输出:<前缀>.pucker.tsv,列是 环大小 花样 量 计数 中位 p05 p95 均值
"""
import math, sys
from collections import Counter, defaultdict
from multiprocessing import Pool
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog("rdApp.*")

HYB = {Chem.HybridizationType.SP2: "2", Chem.HybridizationType.SP3: "3"}

# 每个分子嵌多少个构象取最低能 —— 钉死,换值整张表就对不上。
#
# **30 不够,量出来的。** 中位数在 30 与 100 之间几乎不动(233333i 3.7 → 3.4、
# 223333f1 48.3 → 48.3),但**尾分位数没收敛**:333333i 的 θ p95 从 **88.7 塌到
# 11.3**、333333f1 从 87.2 到 26.0。尾巴是由少数几个分子定的,它们要更多次采样
# 才找得到全局最低。band 的上沿正是尾分位数,所以这个差别直接影响判据的松紧。
NCONF = 100

def canon_pattern(sym):
    """旋转 + 反射下字典序最小"""
    n = len(sym)
    cands = []
    for s in (sym, sym[::-1]):
        for k in range(n):
            cands.append("".join(s[(k + i) % n] for i in range(n)))
    return min(cands)

def cremer_pople(coords):
    """返回 (Q, theta_deg or None)。coords 按环序。"""
    n = len(coords)
    c = coords - coords.mean(axis=0)
    # 平均平面:R1 x R2
    j = np.arange(n)
    R1 = (c * np.sin(2 * math.pi * j / n)[:, None]).sum(axis=0)
    R2 = (c * np.cos(2 * math.pi * j / n)[:, None]).sum(axis=0)
    nv = np.cross(R1, R2)
    nn = np.linalg.norm(nv)
    if nn < 1e-12:
        return None
    nv = nv / nn
    z = c @ nv
    Q = float(np.sqrt((z ** 2).sum()))
    if n == 5:
        # **五元环补相位**:只判幅度分不开信封(E)与扭式(T)——两者 Q 都是 0.395。
        #
        # 裸的 φ 依赖"从哪个原子起数",没法跨分子比。但绕环转一个原子会把 φ 平移
        # 144°,而 **144 mod 36 = 0**,所以 `φ mod 36` 对旋转不变;反射把 φ 取负,
        # 再折到 [0,18] 就对反射也不变。E → 0,T → 18。(数值验过。)
        q2c5 = math.sqrt(2.0 / n) * float((z * np.cos(4 * math.pi * j / n)).sum())
        q2s5 = -math.sqrt(2.0 / n) * float((z * np.sin(4 * math.pi * j / n)).sum())
        if math.hypot(q2c5, q2s5) < 1e-9:
            return (Q, None)
        phi = math.degrees(math.atan2(q2s5, q2c5)) % 36.0
        return (Q, min(phi, 36.0 - phi))
    if n != 6:
        return (Q, None)
    q2c = math.sqrt(2.0 / n) * float((z * np.cos(4 * math.pi * j / n)).sum())
    q2s = -math.sqrt(2.0 / n) * float((z * np.sin(4 * math.pi * j / n)).sum())
    q3 = math.sqrt(1.0 / n) * float((z * ((-1) ** j)).sum())
    q2 = math.hypot(q2c, q2s)
    if Q < 1e-9:
        return (Q, None)
    theta = math.degrees(math.atan2(q2, q3))
    # **折到 [0,90]。** 环的走向是任意的,椅式因此同时出现在 0° 与 180° 两个峰上,
    # 中位数会落在两峰之间的 90° —— 那恰好是扭船所在,band 于是把什么都放过。
    # 折过之后:椅式 → 0,船/扭船 → 90,与走向无关。
    theta = min(theta, 180.0 - theta)
    return (Q, theta)

def judgeable_ring(mh):
    """这个(已补氢的)分子里有没有值得判褶皱的 5/6 元环。

    先筛用的,省掉九成的构象活。口径与下面收环时的一致:5/6 元、非全芳香、
    杂化全认得出来(没有 `0`)。
    """
    for ring in mh.GetRingInfo().AtomRings():
        if len(ring) not in (5, 6):
            continue
        if all(mh.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            continue  # 全芳香环平面,褶皱没什么可判的
        sym = "".join(HYB.get(mh.GetAtomWithIdx(i).GetHybridization(), "0") for i in ring)
        if "0" not in sym:
            return True
    return False

def best_conformer(mh, seed=0xF00D):
    """嵌 `NCONF` 个构象、各自 MMFF 极小、取能量最低的那一个。

    返回 `(坐标, 芳香标志)`;嵌不出来或者一个都不收敛就返回 `None`。

    **`seed` 只有一个正当用途:量这套协议自己的噪声。** 入库的表一律是默认的
    0xf00d;换种子跑出来的表是拿来跟它比"同一个分子换个种子会差多少"的,
    那个差就是逐分子判据的容差下限。别拿换过种子的表入库 ——
    `dump_ring6.py` 会把种子写进文件头,就是为了让这种表一眼认得出来。

    **芳香标志在碰 MMFF 之前取好。** `MMFFOptimizeMoleculeConfs` 会就地把
    芳香性改写成 MMFF 自己的模型 —— 实测 776 个分子里 23 个被改过,
    于是"全芳香就不进表"这条过滤器会用错模型,把 21 个芳香环收进表。
    (判据那边用的是 omgkit 的模型 = RDKit 默认,两边必须对齐。)

    **这一段是 `measure_pucker` 与 `dump_ring6` 共用的**,谁都不许自己抄一份 ——
    两边量的是同一批构象,口径一旦分家,逐环参考就跟统计表对不上了。
    """
    arom = {a.GetIdx(): a.GetIsAromatic() for a in mh.GetAtoms()}
    p = AllChem.ETKDGv3(); p.randomSeed = seed
    try:
        ids = AllChem.EmbedMultipleConfs(mh, numConfs=NCONF, params=p)
        if not ids:
            return None
        res = AllChem.MMFFOptimizeMoleculeConfs(mh, maxIters=1000)
        ok = [(e, cid) for cid, (conv, e) in zip(ids, res) if conv == 0]
        if not ok:
            return None
        best = min(ok)[1]
    except Exception:
        return None
    return mh.GetConformer(best).GetPositions(), arom

def fused_tag(ring, allr):
    """稠合度标签:孤立 `i`,稠合 `f1`/`f2`/`f3`(封顶 3)。

    **稠合度要进键。** 稠环被伙伴钉住,褶皱不自由,与孤立环不是一档;
    而且钉的个数有分别 —— 实测 223333 里稠合 3 个环的 θ 中位 88.7(船),
    而整桶中位是 49.9(扭)。
    """
    nfused = sum(1 for o in allr if o != ring and len(set(ring) & set(o)) >= 2)
    return "i" if nfused == 0 else f"f{min(nfused, 3)}"

def one(job):
    _, smi = job
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    mh = Chem.AddHs(m)
    if not judgeable_ring(mh):
        return Counter()
    got = best_conformer(mh)
    if got is None:
        return None
    c, arom = got
    out = Counter()
    allr = mh.GetRingInfo().AtomRings()
    for ring in allr:
        n = len(ring)
        if n not in (5, 6):
            continue
        sym = "".join(HYB.get(mh.GetAtomWithIdx(i).GetHybridization(), "0") for i in ring)
        if "0" in sym:
            continue
        if all(arom[i] for i in ring):
            continue          # 全芳香环平面,不进表(用 MMFF **之前**取的标志)
        # **全 sp³ 也要进表。** 只收混合杂化的话,判据今天一个环都判不到 ——
        # 现在摆得成的只有全 sp² 与全 sp³ 冠状,而全 sp² 是平面。
        # 判不到东西的判据是个空壳,本仓在 `frame.rs` 上刚为这种写法栽过一次。
        cp = cremer_pople(np.array([c[i] for i in ring]))
        if cp is None:
            continue
        Q, theta = cp
        # 稠合度进键(理由见 `fused_tag`)。**元素没进键**,那是量过之后的决定:
        # 加杂原子计数会把桶从 69 涨到 274,"≥20 个样本的桶"覆盖率从 93% 掉到 63%
        # —— 逐桶比中位数那道闸会对三分之一的环失明。为修 2 个假阳性
        # (`CCCSSS` 那两个)不值当。
        pat = canon_pattern(sym) + fused_tag(ring, allr)
        out[(n, pat, "Q", round(Q, 3))] += 1
        if theta is not None:
            out[(n, pat, "theta" if n == 6 else "phi", round(theta, 1))] += 1
    return out

def summarize(counter, path):
    grp = defaultdict(Counter)
    for key, k in counter.items():
        grp[key[:-1]][key[-1]] += k
    rows = []
    for key, vals in grp.items():
        tot = sum(vals.values())
        xs = sorted(vals)
        run, q = 0, {}
        for v in xs:
            run += vals[v]
            for tag, frac in (("p05", 0.05), ("med", 0.5), ("p95", 0.95)):
                if tag not in q and run >= tot * frac:
                    q[tag] = v
        mean = sum(v * k for v, k in vals.items()) / tot
        rows.append((key, tot, q.get("med"), q.get("p05"), q.get("p95"), mean))
    rows.sort(key=lambda r: (-r[1], r[0]))
    with open(path, "w") as f:
        f.write("#环大小\t花样\t量\t计数\t中位\tp05\tp95\t均值\n")
        for key, tot, med, p05, p95, mean in rows:
            f.write(f"{key[0]}\t{key[1]}\t{key[2]}\t{tot}\t{med}\t{p05}\t{p95}\t{mean:.4f}\n")
    return len(rows), sum(r[1] for r in rows)

if __name__ == "__main__":
    corpus, pre = sys.argv[1], sys.argv[2]
    smis = [l.split("\t")[0].strip() for l in open(corpus) if l.strip()]
    jobs = list(enumerate(smis))
    total = Counter()
    ok = 0
    with Pool() as pool:
        for r in pool.imap_unordered(one, jobs, chunksize=32):
            if r is not None:
                ok += 1
                total.update(r)
    nk, nv = summarize(total, pre + ".pucker.tsv")
    print(f"用上 {ok}/{len(smis)} 个分子(嵌入+MMFF 都收敛的)")
    print(f"褶皱:{nk} 种 (环大小,花样,量) 组合,{nv} 个 → {pre}.pucker.tsv")
