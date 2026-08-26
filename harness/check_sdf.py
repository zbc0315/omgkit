#!/usr/bin/env python3
"""读 SDF(多条记录):外部实现写的文件,我们逐条读出来是不是同一批分子。

# 与"读单条 molblock"是两件事

`check_molblock_read.py` 那条判据每次只喂一条,记录边界整个绕开了。SDF 多出来的
是**切分**:`$$$$` 在哪、`M  END` 在哪、数据字段属于哪一条、一条读不了之后
后面还读不读得下去。这些全是单条判据碰不到的。

# 三样都要比,而且条数要先对上

1. **条数**。我方每条都要占一行(读不了的也占),两侧条数不同就直接失败 ——
   静默跳过一条坏记录会让分母悄悄变小,而"零分歧"在那时依然成立。
2. **分子**。两侧各写回 SMILES,统一交给外部实现规范化再比(跨实现不能直接
   比规范串)。
3. **数据字段**。名字、顺序、多行值、同名重复,一样都不能丢。

# 一条坏记录夹在中间

语料里天然有:金属茂类配合物的键数超出 V2000 的表达能力,写出方自己就换成了
V3000,而我方明确拒收 V3000。它正好落在文件中间 —— 后面那几千条照读不误,
才说明切分是对的。这一档配上限,免得它悄悄长大。

用法:

    python3 harness/check_sdf.py --write  <out.sdf> <语料.smi>
    cargo run -q -p omgkit-io --release --example read_sdf -- <out.sdf> > <ours.txt>
    python3 harness/check_sdf.py --compare <out.sdf> <ours.txt>
"""
import argparse
import json
import pathlib
import sys

import rdkit
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

# 外部实现写成 V3000、我方拒收的条数上限。与 `check_molblock_read.py` 同一档
# (同一份语料、同一个写出方),实测 1 条:二茂铁。
MAX_V3000 = 3

# 带数据字段的条数下限。数据字段那一档一旦被喂空,"零分歧"照样成立 ——
# 而这份文件是判官自己写的,每条都该有字段。贴着现值(全部)留余量。
MIN_WITH_DATA = 8000

# 多行值那一档的条数下限。多行值是最容易在切分时被截断的一种,
# 它单独有个下限,免得写出方哪天不写了而判据一声不响。
MIN_MULTILINE = 8000

# **"骨架相同、立体不同"这一档的上限从别处借来,不在这里另立一个。**
#
# 这条判据管的是**切分**(`$$$$`、`M  END`、数据字段归属),立体感知的边界
# 是另一件事,已经由 `check_molblock_read.py` 在同一份语料上守着。两处各写一个
# 上限的话,那边收紧了这边不会跟着紧 —— 而"这边松着"没有任何理由,只是没人
# 想起来改。所以直接引用那两个常数。
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from check_molblock_read import MAX_CHIRAL_DIFF, MAX_EZ_ONLY  # noqa: E402

MAX_STEREO_DIFF = MAX_CHIRAL_DIFF + MAX_EZ_ONLY


def write_sdf(out_path, corpus, limit):
    """写一份带数据字段的 SDF。

    字段有意挑了三种:一行的、**多行的**、以及**某一行以 `>` 开头**的。

    最后一种是切分的经典陷阱:字段头也以 `>` 开头,按"见 `>` 就开新字段"去切
    的话,值会被拦腰截断、后半段变成一个名字读不出来的字段。放在行中间不算数
    —— 头一版就是那么写的,陷阱一次都没踩到。
    """
    n = 0
    with Chem.SDWriter(out_path) as w:
        for lineno, line in enumerate(open(corpus, encoding="utf-8")):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            smi = line.split("\t")[0].strip()
            m = Chem.MolFromSmiles(smi)
            if m is None:
                continue
            AllChem.Compute2DCoords(m)
            m.SetProp("_Name", f"第{lineno}条")
            m.SetProp("行号", str(lineno))
            m.SetProp("原串", smi)
            m.SetProp("备注", "第一行\n> 这一行以大于号开头,不是字段头\n第三行")
            w.write(m)
            n += 1
            if limit and n >= limit:
                break
    print(f"写了 {n} 条记录到 {out_path}")
    return n


def flat(smiles):
    """抹掉立体之后的规范串。读不了时给 None。"""
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    Chem.RemoveStereochemistry(m)
    try:
        m = Chem.RemoveHs(m)
    except Exception:  # noqa: BLE001
        return None
    return Chem.MolToSmiles(m)


def canon(smiles):
    """规范串。读不了时给 None。

    去显式氢那一步不能省,理由与 `check_molblock_read.py` 里记的一样:外部实现
    的 SMILES 解析器会保留承载方向键的显式氢,而两侧未必都写方向键。
    """
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    try:
        m = Chem.RemoveHs(m)
    except Exception:  # noqa: BLE001
        return None
    return Chem.MolToSmiles(m)


def their_records(path):
    """外部实现读出来的每条:`(分子或 None, 数据字段列表)`。

    用 `ForwardSDMolSupplier` 而不是 `SDMolSupplier`:后者会**跳过**读不了的
    记录,条数就对不上了 —— 而条数正是这条判据要比的东西之一。
    """
    with open(path, "rb") as f:
        sup = Chem.ForwardSDMolSupplier(f, sanitize=True, removeHs=False)
        for m in sup:
            if m is None:
                yield None, []
            else:
                yield m, [(k, m.GetProp(k)) for k in m.GetPropNames()]


def ours(path):
    """我方每条:`(第几条, SMILES 或 <…>, 数据字段)`。

    **字段走 JSON,不走自己拼的转义。** 头一版把换行编成 `\\n`,而语料里有条
    SMILES 是 `[H]/N=c/1\\nc[nH]s1` —— 那两个字符本来就是反斜杠加 n。解码时
    当成换行,判据就报了两条"数据字段不同",而读取器一点毛病没有。
    """
    for line in open(path, encoding="utf-8"):
        idx, smi, data = line.rstrip("\n").split("\t", 2)
        fields = [(k, v) for k, v in json.loads(data)] if data else []
        yield int(idx), smi, fields


def compare(sdf_path, ours_path):
    print(f"外部实现:RDKit {rdkit.__version__}")
    mine = list(ours(ours_path))
    theirs = list(their_records(sdf_path))
    if len(mine) != len(theirs):
        print(f"条数不同:我方 {len(mine)},外部实现 {len(theirs)} —— "
              "有一侧把读不了的记录跳过了,分母对不上")
        return 1

    same = v3000 = with_data = multiline = 0
    both_refused = stereo_diff = 0
    failures = []
    for (idx, smi, fields), (ref, ref_fields) in zip(mine, theirs):
        if smi.startswith("<"):
            if "V3000" in smi:
                v3000 += 1
                continue
            if ref is None:
                both_refused += 1
                continue
            failures.append(f"第 {idx} 条:外部实现读得出,我方 {smi}")
            continue
        if ref is None:
            failures.append(f"第 {idx} 条:外部实现读不了,我方读成了 `{smi}`")
            continue

        want = canon(Chem.MolToSmiles(ref))
        got = canon(smi)
        if got is None:
            failures.append(f"第 {idx} 条:我方写出的 `{smi}` 外部实现读不了")
            continue
        if got != want:
            if flat(smi) is not None and flat(smi) == flat(Chem.MolToSmiles(ref)):
                # 骨架一样、只有立体不同 —— 那是立体感知边界的事,由
                # `check_molblock_read.py` 在同一份语料上守着,上限也在那边。
                stereo_diff += 1
                continue
            failures.append(f"第 {idx} 条:我方 {got},外部实现 {want}")
            continue

        # 数据字段:名字、顺序、值全都要一样。外部实现不把 `_Name` 当数据字段
        # (它是标题,写在第一行),所以两侧比的是同一批。
        if fields != ref_fields:
            failures.append(f"第 {idx} 条:数据字段不同 —— 我方 {fields},外部实现 {ref_fields}")
            continue

        same += 1
        if fields:
            with_data += 1
        if any("\n" in v for _, v in fields):
            multiline += 1

    print(f"逐条一致 {same};不一致 {len(failures)};"
          f"两侧都拒收 {both_refused};我方拒收的 V3000 {v3000} 条(上限 {MAX_V3000})")
    print(f"  带数据字段的 {with_data} 条(下限 {MIN_WITH_DATA});"
          f"其中值是多行的 {multiline} 条(下限 {MIN_MULTILINE})")
    print(f"  骨架相同、只有立体不同的 {stereo_diff} 条"
          f"(上限 {MAX_STEREO_DIFF},借自 check_molblock_read.py)")
    for f in failures[:8]:
        print(f"  ✗ {f}")
    if failures:
        print("\n读出来不是同一批记录。")
        return 1
    if stereo_diff > MAX_STEREO_DIFF:
        print(f"\n只有立体不同的涨到 {stereo_diff} 条,超过上限 {MAX_STEREO_DIFF}")
        return 1
    if v3000 > MAX_V3000:
        print(f"\nV3000 那一档涨到 {v3000} 条,超过上限 {MAX_V3000}")
        return 1
    if with_data < MIN_WITH_DATA:
        print(f"\n带数据字段的只有 {with_data} 条,低于下限 {MIN_WITH_DATA} —— "
              "字段那一档被喂空了")
        return 1
    if multiline < MIN_MULTILINE:
        print(f"\n多行值只有 {multiline} 条,低于下限 {MIN_MULTILINE} —— "
              "最容易被切断的那一档没人验")
        return 1
    print("\n逐条一致(分子、数据字段、条数)。")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", nargs=2, metavar=("OUT", "CORPUS"))
    ap.add_argument("--compare", nargs=2, metavar=("SDF", "OURS"))
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()
    if args.write:
        write_sdf(args.write[0], args.write[1], args.limit)
        return 0
    if args.compare:
        return compare(args.compare[0], args.compare[1])
    ap.error("要么 --write 要么 --compare")
    return 2


if __name__ == "__main__":
    sys.exit(main())
