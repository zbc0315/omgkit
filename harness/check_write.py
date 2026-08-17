#!/usr/bin/env python3
"""用外部实现裁判 omgkit 的 SMILES 写出。

往返测试(`cargo test --test roundtrip_smiles`)只用自家解析器验证写出,
原理上留了一个漏网可能:**解析与写出共享同一个误解**,互为逆运算却都偏离
SMILES 的语义。手性尤其危险 —— 标记写反了,原子数、键集合、连通性全都对,
只有分子是镜像的,纯拓扑比对永远发现不了。

本脚本把两边各自规范化再比字符串,判官是外部实现,与 omgkit 无共谋。

用法:

    cargo run --release -p omgkit-io --example write_smiles -- \\
        harness/corpus/large.smi > /tmp/written.tsv
    python3 harness/check_write.py /tmp/written.tsv

输入是两列 TSV:原始 SMILES、omgkit 写出的 SMILES。

尚未写出的立体信息会被分桶而不是算作失败 —— 见下面 `--strict` 的说明。
"""

import argparse
import pathlib
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

# 关掉 removeHs:它会删掉显式 [H],而**带立体信息的 [H] 会被保留**。
# 于是"没写出双键立体"就间接表现为原子数不同,把一个已知缺口伪装成
# 结构错误。关掉之后两边的图才对得齐。
PARAMS = Chem.SmilesParserParams()
PARAMS.removeHs = False


# 尚未写出的两类立体信息。差别**仅限于**其中之一时不算失败,而是分桶计数:
# 它们是已登记的缺口,不是结构错误。做完之后用 --strict 收紧。
NOT_YET_WRITTEN = {
    "双键立体": lambda m: [
        (b.SetStereo(Chem.BondStereo.STEREONONE), b.SetBondDir(Chem.BondDir.NONE))
        for b in m.GetBonds()
    ],
    "配位几何立体": lambda m: [
        a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        for a in m.GetAtoms()
        if a.GetChiralTag()
        in (
            Chem.ChiralType.CHI_SQUAREPLANAR,
            Chem.ChiralType.CHI_TRIGONALBIPYRAMIDAL,
            Chem.ChiralType.CHI_OCTAHEDRAL,
        )
    ],
}


def canonical_without(mol, drop):
    """抹掉若干类立体信息后的规范 SMILES。`mol` 会被就地修改,故传副本。"""
    for name in drop:
        NOT_YET_WRITTEN[name](mol)
    return Chem.MolToSmiles(mol)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tsv", type=pathlib.Path, help="write_smiles 的输出")
    ap.add_argument("--limit", type=int, default=10, help="最多打印几条分歧")
    ap.add_argument(
        "--strict",
        action="store_true",
        help="不给尚未实现的立体信息开豁免。相应写出做完之后应当打开",
    )
    args = ap.parse_args()

    exact = 0
    chiral = 0
    excused = {name: 0 for name in NOT_YET_WRITTEN}
    bad = []

    for line in args.tsv.read_text().splitlines():
        parts = line.split("\t")
        if len(parts) != 2:
            continue
        orig, got = parts
        a = Chem.MolFromSmiles(orig, PARAMS)
        b = Chem.MolFromSmiles(got, PARAMS)
        if a is None:
            continue  # 原始就读不了,不是写出的问题
        if b is None:
            bad.append((orig, got, "写出结果无法解析"))
            continue

        if any(at.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for at in a.GetAtoms()):
            chiral += 1

        if Chem.MolToSmiles(a) == Chem.MolToSmiles(b):
            exact += 1
            continue

        # 逐类豁免:抹掉某一类尚未写出的立体信息后若一致,说明差别**仅限于**
        # 那一类,不是结构错误。归到对应的桶里而不是算失败。
        matched = None
        if not args.strict:
            for name in NOT_YET_WRITTEN:
                if canonical_without(Chem.Mol(a), [name]) == canonical_without(
                    Chem.Mol(b), [name]
                ):
                    matched = name
                    break
        if matched:
            excused[matched] += 1
        else:
            bad.append((orig, got, f"{Chem.MolToSmiles(a)}\n       != {Chem.MolToSmiles(b)}"))

    print(f"逐条规范形式相同 {exact} 条;含立体中心 {chiral} 条")
    for name, count in excused.items():
        if count:
            print(f"仅 {name} 不同 {count} 条({name}写出尚未实现,已登记)")
    if bad:
        print(f"\n分歧 {len(bad)} 条(最多列 {args.limit} 条):")
        for orig, got, why in bad[: args.limit]:
            print(f"  原: {orig}\n  写: {got}\n  {why}\n")
        return 1
    print("零分歧" + ("(严格模式)" if args.strict else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
