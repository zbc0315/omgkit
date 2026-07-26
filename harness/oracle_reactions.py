#!/usr/bin/env python3
"""生成产物生成的差分基准:「反应语料 × 分子语料」的全部产物。

每行 `反应序号<TAB>分子序号<TAB>产物多重集`,产物形如
`CCCl|CCO.CC`(组内用 `.` 连接同一组的多个产物分子,组间用 `|`,组已排序)。
只写有产物的组合。

# 判据是**规范 SMILES 的多重集**

产物原子的编号是构建顺序留下的痕迹,不是语义量。同一条反应在一个底物上可能有
多处反应位点,产物之间可能重复 —— 所以是多重集,不是集合。

# 净化失败的产物也要记

反应模板本来就能写出价键不合法的产物(`[c:1][H:2]>>[c:1][Br:2]` 这类)。
把它们悄悄丢掉会让两边"都少了同一条",看起来一致,其实谁也没验。所以记成
`<invalid>` 占位,一样参与比对。

用法:

    python3 harness/oracle_reactions.py \\
        --rxns harness/corpus/reactions.txt \\
        --mols harness/corpus/large.smi \\
        --limit-mols 300 \\
        --out harness/baseline/reactions.tsv
"""

import argparse
import pathlib
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

# 与 Rust 侧一致:一个组合最多取这么多组产物。
# 对称底物上产物组数会爆炸,两边必须用同一个上限,否则会出现"只是截断点不同"
# 的假分歧。
MAX_PRODUCT_SETS = 100

PARAMS = Chem.SmilesParserParams()
PARAMS.removeHs = False

INVALID = "<invalid>"


def read_lines(path: pathlib.Path, limit: int) -> list[str]:
    out = []
    for line in path.read_text(errors="replace").splitlines():
        tok = line.split()[0] if line.split() else ""
        if not tok or tok.startswith("#"):
            continue
        out.append(tok)
        if len(out) >= limit:
            break
    return out


def encode(product_sets) -> str:
    groups = []
    for tup in product_sets[:MAX_PRODUCT_SETS]:
        mols = []
        for p in tup:
            try:
                q = Chem.Mol(p)
                Chem.SanitizeMol(q)
                mols.append(Chem.MolToSmiles(q))
            except Exception:
                mols.append(INVALID)
        groups.append(".".join(mols))
    groups.sort()
    return "|".join(groups)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rxns", required=True, type=pathlib.Path)
    ap.add_argument("--mols", required=True, type=pathlib.Path)
    ap.add_argument("--out", required=True, type=pathlib.Path)
    ap.add_argument("--limit-mols", type=int, default=1 << 30)
    args = ap.parse_args()

    rxn_src = read_lines(args.rxns, 1 << 30)
    smis = read_lines(args.mols, args.limit_mols)

    rxns = []
    for r in rxn_src:
        try:
            rxns.append(AllChem.ReactionFromSmarts(r))
        except Exception:
            rxns.append(None)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    pairs = 0
    with_products = 0
    with args.out.open("w") as f:
        for mi, s in enumerate(smis):
            m = Chem.MolFromSmiles(s, PARAMS)
            if m is None:
                continue
            for ri, rxn in enumerate(rxns):
                # 只跑单反应物的模板 —— 多反应物要枚举底物组合,那是另一档测试
                if rxn is None or rxn.GetNumReactantTemplates() != 1:
                    continue
                pairs += 1
                try:
                    sets = rxn.RunReactants((m,))
                except Exception:
                    continue
                if not sets:
                    continue
                with_products += 1
                f.write(f"{ri}\t{mi}\t{encode(sets)}\n")

    print(f"已写出 {args.out}")
    print(f"  反应 {len(rxn_src)},分子 {len(smis)}")
    print(f"  组合 {pairs},有产物 {with_products}")
    print(f"  RDKit: {Chem.rdBase.rdkitVersion}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
