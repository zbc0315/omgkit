#!/usr/bin/env python3
"""规范化的两条自指不变量,不需要任何外部参照。

1. **不动点**:把规范串重新读进来再规范化,必须还是同一串。
   现有判据盖不到 —— 重排判据换的是分子对象的编号,不经过"写出再解析"这一趟,
   而那一趟会同时改掉原子次序、键序**和氢的表示**。
2. **净化幂等**:净化两次要与净化一次相同。不幂等的话,调用方多调一次就得到
   另一个分子,而且不报错。

# 已知残留

参照原子的**挑法**还带着输入的痕迹:`perceive_bond_stereo` 取的是存储顺序里第一个
带方向的邻居,所以方向写在哪一侧,参照就落在哪一侧。双键一端挂着两个取代基时,
两种挑法都对,写出来的方向符号却落在不同的键上 —— 几何一样,串不一样。

实测语料上这样的有 **2 条**,而且都只差一个不携带信息的方向符号(外部实现读三串
都是同一个分子),第二次起就稳定。`--max-bad` 默认放到 2,多一条就红。

用法:python3 harness/check_canonical_fixpoint.py <语料.smi> [--max-bad N] [--limit N]
"""
import argparse
import sys

import omgkit

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("corpus")
ap.add_argument("--max-bad", type=int, default=2, help="容忍几条不动点例外")
ap.add_argument("--limit", type=int, default=10**9)
args = ap.parse_args()
path = args.corpus
limit = args.limit

n = fix_bad = idem_bad = skipped = 0
fix_examples, idem_examples = [], []

for line in open(path):
    smi = line.split()[0] if line.split() else ""
    if not smi or smi.startswith("#"):
        continue
    if n >= limit:
        break
    try:
        m = omgkit.parse_smiles(smi)
        m.sanitize()
    except ValueError:
        skipped += 1
        continue
    n += 1
    c1 = m.to_canonical_smiles()

    # 1. 不动点
    try:
        m2 = omgkit.parse_smiles(c1)
        m2.sanitize()
        c2 = m2.to_canonical_smiles()
    except ValueError as e:
        c2 = f"<读不回/净化不了: {e}>"
    if c1 != c2:
        fix_bad += 1
        if len(fix_examples) < 6:
            fix_examples.append((smi, c1, c2))

    # 2. 净化幂等
    m3 = omgkit.parse_smiles(smi)
    m3.sanitize()
    s_once = m3.to_canonical_smiles()
    try:
        m3.sanitize()
        s_twice = m3.to_canonical_smiles()
    except ValueError as e:
        s_twice = f"<第二次净化失败: {e}>"
    if s_once != s_twice:
        idem_bad += 1
        if len(idem_examples) < 6:
            idem_examples.append((smi, s_once, s_twice))

print(f"#mols\t{n}(另有 {skipped} 条本就净化不了,不计)")
print(f"不动点不成立 {fix_bad} 条")
for smi, a, b in fix_examples:
    print(f"  {smi}\n    第一次 {a}\n    第二次 {b}")
print(f"净化不幂等 {idem_bad} 条")
for smi, a, b in idem_examples:
    print(f"  {smi}\n    一次 {a}\n    两次 {b}")

bad = 0
if fix_bad > args.max_bad:
    print(f"\n不动点例外 {fix_bad} 条,超过上限 {args.max_bad}")
    bad = 1
if idem_bad:
    print(f"\n净化不幂等 {idem_bad} 条 —— 这一条不设上限")
    bad = 1
sys.exit(bad)
