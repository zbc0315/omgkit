#!/usr/bin/env python3
"""用外部实现裁判:图神经网络特征化要读的那一组原子/键描述符。

走的是**产品那条路** —— `omgkit.parse_smiles(...).sanitize()` 之后调
`atom_descriptors()` / `bond_descriptors()`,与写特征化的人调的是同一串。
自己拼输入就会漏掉净化与顺反折算,量到的是另一个分子。

# 比哪些,以及哪一档这条判据够不着

原子 12 项里有 **11 项**在这里逐原子与 RDKit 比。够不着的是 **Pauling 电负性
的值** —— RDKit 没有公开接口能读它(表埋在 `OxidationNumbers.cpp` 里),而拿
生成本表的同一份源码当参照是自证。

这条判据对电负性能看见什么、看不见什么,变异标定过,写清楚:

| 变异 | 这条判据 | `omgkit-core` 的 `pauling_electronegativity_table` |
|---|---|---|
| 把没有公认值的元素补成默认 2.0 | **红**(下限那一档被清空) | 红 |
| 把氟的值从 3.98 改成 3.99 | **绿 —— 看不见** | 红 |

所以电负性的**值**只有那条单元测试守着(它在 CI 主 job 里跑,不要 RDKit);
这里守的是"有值/没值两档都真的出现过"。

键 4 项里,"哪根键带构型"在这里与 RDKit 逐键比;**构型的值与参照原子**由
`check_bond_stereo.py` 比(它已经把参照归一、把小环双键过筛那一套写全了)。
两条判据分工,不各写一遍 —— 各抄一遍必然静默分岔。

# 顺反的口径

本实现给的是 `cis`/`trans`(相对记录的参照原子),不是 `Z`/`E`(要 CIP 优先级,
本仓库没有实现)。所以这里比的是"哪根键带构型"这个**集合**,拿 RDKit 高层解析
认的立体键当真值。值本身留给上面那条判据。

# 覆盖闸

每一项都配一条下限:一列恒为同一个值时,逐位比对当然全绿。分类量要求见过至少
两种取值;稀有的那几档(手性、顺反、非有限电荷、无电负性)各有一条计数下限。

用法:

    python3 harness/check_descriptors.py harness/corpus/large.smi \\
        --extra harness/corpus/hard.smi \\
        --extra harness/corpus/smoke.smi \\
        --extra harness/corpus/descriptors.smi

四份语料各有分工:`large` 是规模与真实分布,`hard` 是构型生成的难例,
`smoke` 里有非四面体立体(`@SP`/`@TB`/`@OH`)这类别处见不到的取值,
`descriptors` 专喂三档前三份都走不到的边界(无电负性的元素、Gasteiger 表外的
金属、同位素标注)。少喂一份,对应那几列就恒为同一个值 —— 逐位比对照样全绿。
"""

from __future__ import annotations

import argparse
import collections
import math
import pathlib
import sys

import denominator
import omgkit
import rdkit
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

# 留着显式氢,下标才与本实现逐个对齐:默认的 `removeHs=True` 会删掉 `[H]`,
# 整串下标错位,比出来的全是假分歧。
_KEEP_HS = Chem.SmilesParserParams()
_KEEP_HS.removeHs = False

HYBRIDIZATION = {
    Chem.HybridizationType.UNSPECIFIED: "unspecified",
    Chem.HybridizationType.S: "s",
    Chem.HybridizationType.SP: "sp",
    Chem.HybridizationType.SP2: "sp2",
    Chem.HybridizationType.SP3: "sp3",
    Chem.HybridizationType.SP2D: "sp2d",
    Chem.HybridizationType.SP3D: "sp3d",
    Chem.HybridizationType.SP3D2: "sp3d2",
}
CHIRAL_TAG = {
    Chem.ChiralType.CHI_UNSPECIFIED: "unspecified",
    Chem.ChiralType.CHI_TETRAHEDRAL_CW: "cw",
    Chem.ChiralType.CHI_TETRAHEDRAL_CCW: "ccw",
    Chem.ChiralType.CHI_ALLENE: "allene",
    Chem.ChiralType.CHI_SQUAREPLANAR: "square_planar",
    Chem.ChiralType.CHI_TRIGONALBIPYRAMIDAL: "trigonal_bipyramidal",
    Chem.ChiralType.CHI_OCTAHEDRAL: "octahedral",
}
BOND_ORDER = {
    Chem.BondType.UNSPECIFIED: "unspecified",
    Chem.BondType.SINGLE: "single",
    Chem.BondType.DOUBLE: "double",
    Chem.BondType.TRIPLE: "triple",
    Chem.BondType.QUADRUPLE: "quadruple",
    Chem.BondType.AROMATIC: "aromatic",
    Chem.BondType.DATIVE: "dative",
}

# 原子量:两边的表都是同一串十进制文本解出的 double,**要求逐位相同**。
# 头一版这里写着 1e-5,理由是"本实现存的是 f32,12.011 读回来是 12.0109996…"。
# 那是真的,但那是**可以修的**,不是该配容差的。容差一旦立在这儿,后面
# 抄错一位数字也照样绿。核心的元素表随后改成了 f64,这条就收成 0。
MASS_TOL = 0.0
# Gasteiger 电荷:同一套公式、同一个迭代次数、连"(a−b)+b 不化简成 a"都照抄了,
# 只剩邻居求和的**次序**不同 —— 实测差在 1e-17 量级(双精度的末一两位)。
# 收到 0 会红 48892 处,那些全不是分歧。这个数不许再往上放:它比实测差距
# 宽了五个数量级,已经足够;放宽等于给"抄错一位"留门。
CHARGE_TOL = 1e-12

# 刻意分歧,逐条钉死。键是 `(SMILES, 哪一处)`,值是 `(本实现, 参照, 理由)`。
#
# **这张表是双向的**:多一条新的分歧红,少一条也红。后者逼着改的人回来确认
# "消失"是对的,而不是把判据悄悄放松了一格。
#
# 三类的根因全在**参照那一侧**,而且都是本仓别处已经查明并记过的:
# 自由基中心的手性 RDKit 会清掉(见 `differential_l3.rs` 的 `NOT_CONVERGENT`)、
# 丙二烯轴手性 RDKit 完全不支持(见 `harness/requirements.lock` 里为什么要引
# Indigo)、超配位中心的四面体声称 RDKit 清掉。
PINNED = {
    ("[C@@H]1CCCCC1O", "原子 0 的 chiral_tag"): (
        "ccw", "unspecified",
        "两根重原子键 + 方括号里一个氢 ⇒ 参照把它读成自由基 [CH],"
        "`AssignStereochemistry` 随即清掉标记。differential_l3.rs 的 "
        "NOT_CONVERGENT 表里记的是同一条。",
    ),
    ("[C@H]1CCCCC1O", "原子 0 的 chiral_tag"): (
        "cw", "unspecified", "同上,是上一条的对映体。",
    ),
    ("[C@@H]1CCCC1", "原子 0 的 chiral_tag"): (
        "ccw", "unspecified", "同上,五元环。",
    ),
    ("[N@@H]1CCCC1", "原子 0 的 chiral_tag"): (
        "ccw", "unspecified", "同上,氮。",
    ),
    ("N[C@AL1]=C=C(O)F", "原子 1 的 chiral_tag"): (
        "allene", "unspecified",
        "**参照完全不支持丙二烯型轴手性** —— 逐条实测过六条路(读写、"
        "FindPotentialStereo、rdCIPLabeler、从三维反推、molblock 往返、"
        "带手性的子结构匹配)全把 @AL1 与 @AL2 当同一个东西。这一档的裁判"
        "是 Indigo,见 harness/check_allene.py。",
    ),
    ("[S@](C)(C)(C)C", "原子 0 的 chiral_tag"): (
        "ccw", "unspecified", "四配位中性硫(超价声称),参照清掉标记。",
    ),
    ("[Xe@](F)(F)(F)F", "原子 0 的 chiral_tag"): (
        "ccw", "unspecified", "四配位氙,同上。",
    ),
}

# 覆盖下限。**不是宽容度**:这些数守的是"这一档真的被比到过"。
MIN_WITH_CHIRAL = 200  # 带四面体手性的原子
MIN_WITH_EZ = 100  # 带顺反构型的键
MIN_AROMATIC = 10000  # 芳香原子
MIN_IN_RING = 10000  # 环上原子
MIN_CHARGED = 500  # 带形式电荷的原子
MIN_NO_ELECTRONEGATIVITY = 5  # 没有公认 Pauling 值的原子(只有边界语料喂得到)
MIN_INVALID_GASTEIGER = 100  # 电荷算不出来的原子
MIN_ISOTOPE_MASS = 10  # 标了同位素、原子量走精确质量那条路的原子


def reference_atoms(mol) -> list[dict]:
    return [
        {
            "atomic_num": a.GetAtomicNum(),
            "total_degree": a.GetTotalDegree(),
            "formal_charge": a.GetFormalCharge(),
            "chiral_tag": CHIRAL_TAG.get(a.GetChiralTag(), str(a.GetChiralTag())),
            "total_num_hs": a.GetTotalNumHs(),
            "hybridization": HYBRIDIZATION.get(
                a.GetHybridization(), str(a.GetHybridization())
            ),
            "is_aromatic": a.GetIsAromatic(),
            "is_in_ring": a.IsInRing(),
            "mass": a.GetMass(),
            "gasteiger_charge": a.GetDoubleProp("_GasteigerCharge"),
        }
        for a in mol.GetAtoms()
    ]


def reference_bonds(mol) -> dict:
    """{(小下标, 大下标): {...}}。按端点索引,不按遍历序 —— 遍历序是实现痕迹。"""
    out = {}
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        out[(min(i, j), max(i, j))] = {
            "order": BOND_ORDER.get(b.GetBondType(), str(b.GetBondType())),
            "is_conjugated": b.GetIsConjugated(),
            "is_in_ring": b.IsInRing(),
            "has_stereo": b.GetStereo() != Chem.BondStereo.STEREONONE,
        }
    return out


def record(bad, pinned_hit, smi, where, ours, theirs):
    """记一处分歧;在钉死表里、而且值也对得上,就归到那边,不算红。"""
    want = PINNED.get((smi, where))
    if want is not None and (str(want[0]), str(want[1])) == (str(ours), str(theirs)):
        pinned_hit.add((smi, where))
    else:
        bad.append((smi, where, ours, theirs))


def compare_atoms(smi, ours, theirs, bad, pinned_hit, seen):
    for i, (o, t) in enumerate(zip(ours, theirs)):
        for key in (
            "atomic_num",
            "total_degree",
            "formal_charge",
            "chiral_tag",
            "total_num_hs",
            "hybridization",
            "is_aromatic",
            "is_in_ring",
        ):
            seen[key].add(o[key])
            if o[key] != t[key]:
                record(bad, pinned_hit, smi, f"原子 {i} 的 {key}", o[key], t[key])
        if abs(o["mass"] - t["mass"]) > MASS_TOL:
            record(bad, pinned_hit, smi, f"原子 {i} 的 mass", o["mass"], t["mass"])
        seen["electronegativity"].add(o["electronegativity"] is None)

        # 电荷:先比"算不算得出来",两边都算得出来才比值。
        # 反过来做的话 nan != nan 会把每一个表外原子都报成分歧,而真正的
        # 问题(该失效的没失效)反倒淹在里面。
        ok_t = math.isfinite(t["gasteiger_charge"])
        seen["gasteiger_valid"].add(o["gasteiger_valid"])
        if o["gasteiger_valid"] != ok_t:
            record(
                bad, pinned_hit, smi,
                f"原子 {i} 的 gasteiger_valid", o["gasteiger_valid"], ok_t,
            )
        elif ok_t and abs(o["gasteiger_charge"] - t["gasteiger_charge"]) > CHARGE_TOL:
            record(
                bad, pinned_hit, smi, f"原子 {i} 的 gasteiger_charge",
                o["gasteiger_charge"], t["gasteiger_charge"],
            )


def compare_bonds(smi, ours, theirs, bad, pinned_hit, seen):
    keyed = {(min(b["begin"], b["end"]), max(b["begin"], b["end"])): b for b in ours}
    if set(keyed) != set(theirs):
        bad.append((smi, "键集合", sorted(keyed), sorted(theirs)))
        return
    for k, o in keyed.items():
        t = theirs[k]
        for key in ("order", "is_conjugated", "is_in_ring"):
            seen[f"bond.{key}"].add(o[key])
            if o[key] != t[key]:
                record(bad, pinned_hit, smi, f"键 {k} 的 {key}", o[key], t[key])
        has = o["stereo"] != "none"
        seen["bond.stereo"].add(o["stereo"])
        if has != t["has_stereo"]:
            record(bad, pinned_hit, smi, f"键 {k} 带不带构型", has, t["has_stereo"])
        if has and o["stereo_atoms"] is None:
            bad.append((smi, f"键 {k} 带构型却没有参照原子", None, "两个下标"))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("corpus", type=pathlib.Path)
    ap.add_argument("--extra", type=pathlib.Path, action="append", default=[])
    ap.add_argument(
        "--cap",
        type=int,
        default=0,
        help="允许有多少条没进比对(分母闸,不是宽容度)",
    )
    args = ap.parse_args()

    print(f"  RDKit {rdkit.__version__}")
    print(f"  omgkit {omgkit.__version__} @ {omgkit.__file__}")

    paths = [args.corpus, *args.extra]
    smiles = []
    for p in paths:
        for line in p.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                smiles.append(line.split("\t")[0])
    n_corpus = sum(denominator.corpus_size(p) for p in paths)

    bad = []
    skipped = collections.Counter()
    seen = collections.defaultdict(set)
    n_chiral = n_ez = n_aromatic = n_in_ring = n_charged = 0
    n_no_en = n_invalid_q = n_isotope = 0
    pinned_hit = set()
    compared = 0

    both_reject = 0
    for smi in smiles:
        ref = Chem.MolFromSmiles(smi, _KEEP_HS)
        try:
            mol = omgkit.parse_smiles(smi)
            mol.sanitize()
        except ValueError:
            # **两边都拒绝**不是覆盖漏洞,是两边对同一个超价分子的一致判断;
            # **只有本实现拒绝**才是 —— 那样这条分子的描述符从来没被比过,
            # 而判据会因为少了一个分歧源而变好看。两者必须分开数。
            if ref is None:
                both_reject += 1
            else:
                skipped["只有本实现净化不了(RDKit 认)"] += 1
            continue
        if ref is None:
            skipped["只有 RDKit 解析不了"] += 1
            continue
        ours_a = mol.atom_descriptors()
        if len(ours_a) != ref.GetNumAtoms():
            skipped["原子数对不上"] += 1
            continue
        from rdkit.Chem import rdPartialCharges

        rdPartialCharges.ComputeGasteigerCharges(ref)
        compare_atoms(smi, ours_a, reference_atoms(ref), bad, pinned_hit, seen)
        compare_bonds(
            smi, mol.bond_descriptors(), reference_bonds(ref), bad, pinned_hit, seen
        )
        compared += 1

        n_chiral += sum(1 for d in ours_a if d["chiral_tag"] in ("cw", "ccw"))
        n_aromatic += sum(1 for d in ours_a if d["is_aromatic"])
        n_in_ring += sum(1 for d in ours_a if d["is_in_ring"])
        n_charged += sum(1 for d in ours_a if d["formal_charge"] != 0)
        n_no_en += sum(1 for d in ours_a if d["electronegativity"] is None)
        n_invalid_q += sum(1 for d in ours_a if not d["gasteiger_valid"])
        n_isotope += sum(1 for a in ref.GetAtoms() if a.GetIsotope())
        n_ez += sum(1 for b in mol.bond_descriptors() if b["stereo"] != "none")

    print(f"逐原子/逐键比对 {compared} 条分子;分歧 {len(bad)} 处")
    print(f"  两边都拒绝(超价等,不计入分母):{both_reject} 条")
    for why, k in sorted(skipped.items()):
        print(f"  未进比对:{why} {k} 条")
    for f in bad[:15]:
        print(f"  ✗ {f[0]}:{f[1]} —— 本实现 {f[2]!r},RDKit {f[3]!r}")
    if len(bad) > 15:
        print(f"  ...(另有 {len(bad) - 15} 处)")
    if bad:
        return 1

    # **钉死表是双向的**:少一条同样红
    missing = sorted(set(PINNED) - pinned_hit)
    if missing:
        print(f"\n钉死的刻意分歧有 {len(missing)} 条这次没出现:")
        for smi, where in missing:
            print(f"  · {smi}:{where} —— {PINNED[(smi, where)][2]}")
        print(
            "分歧消失可能是好事(参照升级了、或本实现改了口径),但**要有人确认**:\n"
            "确认之后把这一条从 PINNED 里删掉。留着一条永不命中的例外,"
            "等于给判据挖了一个没人看的洞。"
        )
        return 1
    print(f"  刻意分歧(根因在参照侧,逐条钉死):{len(pinned_hit)} 条")

    why = denominator.verdict(
        n_corpus - both_reject, len(smiles) - both_reject, compared, args.cap
    )
    if why:
        print(f"\n{why}")
        return 1

    # ---- 覆盖闸:一列恒为同一个值时,逐位比对当然全绿 ----
    print(
        f"  取值覆盖:杂化 {sorted(seen['hybridization'])}\n"
        f"            手性 {sorted(seen['chiral_tag'])}\n"
        f"            键级 {sorted(seen['bond.order'])}\n"
        f"            顺反 {sorted(seen['bond.stereo'])}"
    )
    print(
        f"  计数:手性原子 {n_chiral}、顺反键 {n_ez}、芳香原子 {n_aromatic}、"
        f"环上原子 {n_in_ring}、带电原子 {n_charged}、\n"
        f"        无电负性原子 {n_no_en}、电荷算不出的原子 {n_invalid_q}、"
        f"同位素标注原子 {n_isotope}"
    )
    floors = [
        ("带四面体手性的原子", n_chiral, MIN_WITH_CHIRAL),
        ("带顺反构型的键", n_ez, MIN_WITH_EZ),
        ("芳香原子", n_aromatic, MIN_AROMATIC),
        ("环上原子", n_in_ring, MIN_IN_RING),
        ("带形式电荷的原子", n_charged, MIN_CHARGED),
        ("没有公认电负性的原子", n_no_en, MIN_NO_ELECTRONEGATIVITY),
        ("电荷算不出来的原子", n_invalid_q, MIN_INVALID_GASTEIGER),
        ("标了同位素的原子", n_isotope, MIN_ISOTOPE_MASS),
    ]
    for what, got, floor in floors:
        if got < floor:
            print(f"\n{what}只有 {got} 个,低于下限 {floor} —— 这一档被喂空了")
            return 1
    for key in (
        "hybridization",
        "chiral_tag",
        "is_aromatic",
        "is_in_ring",
        "electronegativity",
        "gasteiger_valid",
        "bond.order",
        "bond.is_conjugated",
        "bond.is_in_ring",
        "bond.stereo",
    ):
        if len(seen[key]) < 2:
            print(f"\n{key} 全语料只出现过一种取值 {seen[key]} —— 这一列恒定,比不出东西")
            return 1

    print("\n16 个描述符与外部实现逐位相同(电负性的值由单元测试钉住,见文件头)。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
