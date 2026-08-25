#!/usr/bin/env python3
"""用 RDKit 生成差分测试基准。

这是整个项目最重要的一块基础设施。omgkit 的每一层都必须能对 ChEMBL 全量
跑一遍逐字段比对,"分歧数"是看板上的头号 KPI。**先有尺子,再造东西。**

分层基准(与 omgkit-phase1-scope.md 的 L1..L3 对应):

  l1  纯解析,不净化(`MolFromSmiles(sanitize=False)`)
      比对:原子数/键数、邻接表、元素、电荷、同位素、显式氢、映射号、手性标记
      —— 这一层刻意不含任何化学语义,能把解析器的 bug 与净化的 bug 分开。

  l2  完整净化(`MolFromSmiles()` 默认)
      比对:在 l1 基础上增加 芳香标志、隐式氢数、显式价、杂化、环成员、最小环

  l3  规范化输出
      比对:canonical SMILES 逐字节 —— 王牌测试,一次覆盖 L1+L2+L3 的全部正确性

注意:**`removeHs` 不是净化的一部分。** `Chem.MolFromSmiles(smi)` 实际做了三件事:
解析 + 净化 + `removeHs`。后者会删掉显式 `[H]` 原子,**改变原子数**
(实测 `O1C[C@]1(CCCCCCCC)[H]`:12 原子 → 11 原子)。

omgkit 把 `removeHs` 划为独立操作,不属于 L2。因此 l2 基准默认
`--remove-hs` 关闭,保证 l1 与 l2 的图可逐项对齐。生成 l3 基准时应显式
打开 `--remove-hs`,以对齐用户实际调用的
`Chem.MolToSmiles(Chem.MolFromSmiles(x))`。

另一个已知的图变更:净化第 2 步 `cleanUpOrganometallics` 会把某些键改成配位键
并**交换端点**(实测 `CN(C)C[C-]12...[Fe++]...`:`(9,4) SINGLE` → `(4,9) DATIVE`)。
无向边集不变,故环感知等纯图算法不受影响。

输出为 JSONL,每行一条记录,字段用紧凑数组以便跑到 ChEMBL 规模(~240 万行)。

用法:
    python3 harness/oracle_pipeline.py --input corpus/smoke.smi --stage l1 \\
        --out baseline/smoke.l1.jsonl
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys

# RDKit 源码树里有个同名的 `rdkit/` 目录,cwd 落在里面时 import 会抓到源码
# 而非已安装的包,报 "partially initialized module"。提前给出可读的报错。
if (pathlib.Path.cwd() / "rdkit" / "__init__.py").exists() and not (
    pathlib.Path.cwd() / "rdkit" / "Chem" / "rdchem.so"
).exists():
    sys.exit(
        "当前目录下有 RDKit **源码**目录,它会遮蔽已安装的 rdkit 包。\n"
        "请换到 omgkit/ 或任何不含 rdkit 源码的目录再运行。"
    )

try:
    from rdkit import Chem, RDLogger
except ImportError:  # pragma: no cover
    sys.exit("需要 rdkit:pip install rdkit")

RDLogger.DisableLog("rdApp.*")  # 语料里必然有非法 SMILES,不刷屏

STAGES = ("l1", "l2", "l3")

# 净化的 12 个步骤,顺序取自 RDKit `MolOps.cpp:584 sanitizeMol()`。
# 有了这张表就能只跑管线的任意前缀,从而**逐步**验证 —— 否则单独实现某一步
# 时无从比对:l2 基准反映的是 12 步全跑完的最终状态,第 7 步芳香性感知会
# 改写芳香标志,最后还会重算一次价键。
SANITIZE_OPS = {
    "CLEANUP": Chem.SanitizeFlags.SANITIZE_CLEANUP,
    "CLEANUP_ORGANOMETALLICS": Chem.SanitizeFlags.SANITIZE_CLEANUP_ORGANOMETALLICS,
    "PROPERTIES": Chem.SanitizeFlags.SANITIZE_PROPERTIES,
    "SYMMRINGS": Chem.SanitizeFlags.SANITIZE_SYMMRINGS,
    "KEKULIZE": Chem.SanitizeFlags.SANITIZE_KEKULIZE,
    "FINDRADICALS": Chem.SanitizeFlags.SANITIZE_FINDRADICALS,
    "SETAROMATICITY": Chem.SanitizeFlags.SANITIZE_SETAROMATICITY,
    "SETCONJUGATION": Chem.SanitizeFlags.SANITIZE_SETCONJUGATION,
    "SETHYBRIDIZATION": Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION,
    "CLEANUPATROPISOMERS": Chem.SanitizeFlags.SANITIZE_CLEANUPATROPISOMERS,
    "CLEANUPCHIRALITY": Chem.SanitizeFlags.SANITIZE_CLEANUPCHIRALITY,
    "ADJUSTHS": Chem.SanitizeFlags.SANITIZE_ADJUSTHS,
}

# ---- 编码表:必须与 omgkit-core 的 #[repr(u8)] 判别值逐一对应 ----
# 任何一处对不上,差分测试会把编码错误报成化学错误,极难定位。

BOND_ORDER = {
    Chem.BondType.UNSPECIFIED: 0,
    Chem.BondType.SINGLE: 1,
    Chem.BondType.DOUBLE: 2,
    Chem.BondType.TRIPLE: 3,
    Chem.BondType.QUADRUPLE: 4,
    Chem.BondType.AROMATIC: 5,
    Chem.BondType.DATIVE: 6,
}

# 立体标记的**几何类别**。未列出的一律落到 3(其它),丙二烯轴手性即在此列。
CHIRAL_TAG = {
    Chem.ChiralType.CHI_UNSPECIFIED: 0,
    Chem.ChiralType.CHI_TETRAHEDRAL_CW: 1,
    Chem.ChiralType.CHI_TETRAHEDRAL_CCW: 2,
    Chem.ChiralType.CHI_SQUAREPLANAR: 4,
    Chem.ChiralType.CHI_TRIGONALBIPYRAMIDAL: 5,
    Chem.ChiralType.CHI_OCTAHEDRAL: 6,
}

BOND_DIR = {
    Chem.BondDir.NONE: 0,
    Chem.BondDir.ENDUPRIGHT: 1,    # '/'
    Chem.BondDir.ENDDOWNRIGHT: 2,  # '\'
}

HYBRIDIZATION = {
    Chem.HybridizationType.UNSPECIFIED: 0,
    Chem.HybridizationType.S: 1,
    Chem.HybridizationType.SP: 2,
    Chem.HybridizationType.SP2: 3,
    Chem.HybridizationType.SP3: 4,
    Chem.HybridizationType.SP2D: 5,
    Chem.HybridizationType.SP3D: 6,
    Chem.HybridizationType.SP3D2: 7,
}


def _atom_core(a) -> list[int]:
    """两个阶段共用的前 8 列。列号在此固定,新增列一律追加到行尾。"""
    return [
        a.GetAtomicNum(),
        a.GetFormalCharge(),
        a.GetIsotope(),
        a.GetNumExplicitHs(),
        a.GetAtomMapNum(),
        CHIRAL_TAG.get(a.GetChiralTag(), 3),
        int(a.GetIsAromatic()),
        int(a.GetNoImplicit()),
    ]


def _stereo_perm(a) -> int:
    """立体标记的类内排列序号。四面体的排列由标记本身表达,故为 0。"""
    return a.GetPropsAsDict().get("_chiralPermutation", 0)


def encode_atoms_l1(mol) -> list[list[int]]:
    """核心 8 列 + [排列序号]

    末三列在 l1 阶段就有意义:
      - 芳香 = 小写字母的字面声称(不是感知结果,那是 l2)
      - 不推断隐式氢 = 原子写在方括号中
      - 排列序号 = `@TB15` 里的 15;四面体为 0
    """
    return [_atom_core(a) + [_stereo_perm(a)] for a in mol.GetAtoms()]


def encode_atoms_l2(mol) -> list[list[int]]:
    """核心 8 列 + [隐式氢, 显式价, 杂化, 在环中, 最小环大小, 自由基电子数]
    + [排列序号]

    只跑部分净化步骤时,环信息可能尚未初始化。此时"最小环大小"列填 **-1**
    表示"该步骤未运行,数据不可用" —— 刻意区别于 0(表示"不在任何环中"),
    这样消费方拿它去比对会立刻炸出来,而不是悄悄比了个假值。

    注意:排列序号追加在**行尾**(l2 的第 14 列),不是跟在核心列后面。两个阶段
    的行长不同,若把新列插在中间,l2 独有的那 6 列会整体后移,所有硬编码
    列号的比对都会静默比错字段。**新列一律追加到行尾。**
    """
    out = []
    ri = mol.GetRingInfo()
    for a in mol.GetAtoms():
        idx = a.GetIdx()
        min_ring = 0
        if a.IsInRing():
            try:
                sizes = [s for s in range(3, 21) if ri.IsAtomInRingOfSize(idx, s)]
                min_ring = sizes[0] if sizes else 0
            except RuntimeError:
                min_ring = -1  # 环感知步骤未运行
        out.append(
            _atom_core(a)
            + [
                a.GetNumImplicitHs(),
                a.GetExplicitValence(),
                HYBRIDIZATION.get(a.GetHybridization(), 0),
                int(a.IsInRing()),
                min_ring,
                a.GetNumRadicalElectrons(),
            ]
            + [_stereo_perm(a)]
        )
    return out


def encode_rings(mol) -> list[list[int]] | None:
    """环集:每个环一个**已排序的原子下标列表**。

    环感知步骤未运行时返回 None(区别于 [] = 确实没有环),
    让误比对立刻炸出来。

    注意:刻意排序:RDKit 给的是遍历顺序,那是实现痕迹。实测 7553 条含环分子、
    每条随机重排原子编号 3 次,**环的原子集合无一改变** —— 集合是规范的,
    顺序不是。比对必须按集合。
    """
    try:
        ri = mol.GetRingInfo()
        return [sorted(r) for r in ri.AtomRings()]
    except RuntimeError:
        return None


def encode_bonds(mol) -> list[list[int]]:
    """[起点, 终点, 键级, 方向, 芳香, 在环中, 共轭]

    注意:**第 5 列(在环中)在 l1 阶段不可用于比对。** 实测:
    `MolFromSmiles(sanitize=False)` 之后 `RingInfo.NumRings()` 抛 RuntimeError,
    但 `bond.IsInRing()` 仍返回**真实**的环成员信息(对非芳香环也正确)——
    因为 RDKit 的 SMILES 解析器对未净化分子会调 `fastFindRings`。

    omgkit 把环感知划在 L2,所以 L1 的比对必须**跳过这一列**,到 L2 再比。
    详见 harness/README.md 的列规范。
    """
    return [
        [
            b.GetBeginAtomIdx(),
            b.GetEndAtomIdx(),
            BOND_ORDER.get(b.GetBondType(), 0),
            BOND_DIR.get(b.GetBondDir(), 0),
            int(b.GetIsAromatic()),
            int(b.IsInRing()),
            int(b.GetIsConjugated()),
        ]
        for b in mol.GetBonds()
    ]


def record(
    idx: int,
    smiles: str,
    name: str,
    stage: str,
    remove_hs: bool,
    ops: int | None = None,
) -> dict:
    """`ops` 非 None 时:先不净化地解析,再只跑指定的净化步骤位掩码。

    # 为什么 `removeHs` 不能留在解析这一步

    `MolFromSmiles(params)` 里 `params.sanitize=False` 且 `params.removeHs=True`
    时,RDKit 会把**方括号里的氢数、`noImplicit` 标志和手性标记一起抹掉**。
    实测(RDKit 2025.09.2,`[C@@H]1CCCCC1O`):

    | 解析参数 | 原子 0 的手性标记 | 显式氢 | `noImplicit` | 规范串 |
    |---|---|---|---|---|
    | `removeHs=False` | `CHI_TETRAHEDRAL_CCW` | 1 | True | `OC1[CH]CCCC1` |
    | `removeHs=True` | `CHI_UNSPECIFIED` | **0** | **False** | `OC1CCCCC1` |

    第二行不是"少了几个氢原子",是**另一个分子**:2-羟基环己基自由基
    变成了环己醇。语料里三条(`chiral-ring-open-cw` / `-ccw` /
    `chiral-cyclopentane`)都中招,而它们正是拿来区分对映体的用例。

    改成净化之后再 `Chem.RemoveHs` 就没有这个问题(上表第一行的结果),
    而且与 omgkit 那侧的调用顺序一致 —— 本文件开头就写着"omgkit 把
    `removeHs` 划为独立操作,不属于 L2"。
    """
    rec: dict = {"i": idx, "smi": smiles}
    if name:
        rec["name"] = name

    params = Chem.SmilesParserParams()
    # 净化一律自己控制,解析阶段不净化。
    #
    # 注意:关键:`MolFromSmiles(sanitize=True)` 做的**不止**那 12 步 —— 它之后还跑
    # `assignStereochemistry(cleanIt=True)`,那一步会清掉非真手性中心的标记,
    # 并顺带把显式氢改回隐式。那属于立体化学感知(L6),不属于净化。
    #
    # 之前只有"不指定 --sanitize-ops"时走 sanitize=True,于是"全 12 步"这份基准
    # 与逐步骤的基准量的**不是同一件事**,而这个差别只在带手性的分子上显形。
    params.sanitize = False
    # 显式关掉,不吃默认值(默认是 True)。**`removeHs` 挪到净化之后**,
    # 见下方那处调用与函数文档"为什么 removeHs 不能留在解析这一步"。
    params.removeHs = False
    try:
        mol = Chem.MolFromSmiles(smiles, params)
    except Exception as e:  # RDKit 偶尔抛而非返回 None
        rec["ok"] = False
        rec["err"] = f"{type(e).__name__}: {e}"
        return rec

    if mol is None:
        rec["ok"] = False
        rec["err"] = "MolFromSmiles returned None"
        return rec

    if stage != "l1":
        ops = ops if ops is not None else int(Chem.SanitizeFlags.SANITIZE_ALL)
        try:
            failed = Chem.SanitizeMol(mol, sanitizeOps=ops, catchErrors=True)
        except Exception as e:
            rec["ok"] = False
            rec["err"] = f"SanitizeMol {type(e).__name__}: {e}"
            return rec
        if int(failed) != 0:
            rec["ok"] = False
            rec["err"] = f"净化步骤失败: {failed!s}"
            return rec

    # 去氢排在净化**之后** —— 留在解析那一步会改分子,见函数文档。
    if remove_hs:
        # `sanitize=False`:这里不能再净化一遍,否则 `--sanitize-ops` 那条
        # "只跑指定步骤"的约定会被悄悄破坏。
        mol = Chem.RemoveHs(mol, sanitize=False)

    rec["ok"] = True
    rec["na"] = mol.GetNumAtoms()
    rec["nb"] = mol.GetNumBonds()

    if stage == "l1":
        rec["atoms"] = encode_atoms_l1(mol)
        rec["bonds"] = encode_bonds(mol)
    elif stage == "l2":
        rec["atoms"] = encode_atoms_l2(mol)
        rec["bonds"] = encode_bonds(mol)
        rec["rings"] = encode_rings(mol)
    else:  # l3
        rec["can"] = Chem.MolToSmiles(mol)
        # 打乱原子顺序后规范形式必须不变 —— 这条性质 omgkit 可以完全
        # 自测(不依赖 RDKit),此处顺带记录 RDKit 自己的表现作为参照。
        rec["can_renumbered"] = Chem.MolToSmiles(
            Chem.RenumberAtoms(mol, list(reversed(range(mol.GetNumAtoms()))))
        )
    return rec


def read_smi(path: pathlib.Path):
    """读 .smi:每行 `SMILES[<空白>名字]`,忽略空行与 # 注释。"""
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(None, 1)
            yield parts[0], (parts[1].strip() if len(parts) > 1 else "")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, type=pathlib.Path, help=".smi 语料")
    ap.add_argument("--out", required=True, type=pathlib.Path, help="输出 JSONL")
    ap.add_argument("--stage", choices=STAGES, default="l1")
    ap.add_argument("--limit", type=int, default=0, help="只处理前 N 条(0 = 全部)")
    ap.add_argument(
        "--remove-hs",
        action="store_true",
        help="额外执行 removeHs(会删掉显式 [H],改变原子数)。"
        "l1/l2 应保持关闭;生成 l3 基准时打开以对齐用户实际调用。",
    )
    ap.add_argument(
        "--sanitize-ops",
        default=None,
        help="逗号分隔的净化步骤名,只跑这些步骤(用于逐步验证)。"
        f"可选:{','.join(SANITIZE_OPS)}",
    )
    args = ap.parse_args()

    if args.remove_hs and args.stage == "l1":
        sys.exit("--remove-hs 对 l1 无意义:去氢要在净化之后跑,而 l1 不净化")

    ops = None
    if args.sanitize_ops is not None:
        if args.stage == "l1":
            sys.exit("--sanitize-ops 对 l1 无意义(l1 本就不净化)")
        names = [s.strip().upper() for s in args.sanitize_ops.split(",") if s.strip()]
        unknown = [n for n in names if n not in SANITIZE_OPS]
        if unknown:
            sys.exit(f"未知的净化步骤:{unknown};可选:{list(SANITIZE_OPS)}")
        ops = 0
        for n in names:
            ops |= int(SANITIZE_OPS[n])

    if not args.input.is_file():
        sys.exit(f"找不到语料 {args.input}")

    args.out.parent.mkdir(parents=True, exist_ok=True)

    n = n_ok = 0
    with args.out.open("w", encoding="utf-8") as fh:
        for i, (smi, name) in enumerate(read_smi(args.input)):
            if args.limit and n >= args.limit:
                break
            rec = record(i, smi, name, args.stage, args.remove_hs, ops)
            fh.write(json.dumps(rec, ensure_ascii=False, separators=(",", ":")) + "\n")
            n += 1
            n_ok += bool(rec.get("ok"))

    import rdkit

    print(f"已写出 {args.out}")
    print(f"  阶段: {args.stage}(removeHs={args.remove_hs}, "
          f"净化步骤={args.sanitize_ops or '全部'})")
    print(f"  记录: {n}(解析成功 {n_ok},失败 {n - n_ok})")
    print(f"  RDKit: {rdkit.__version__}")


if __name__ == "__main__":
    main()
