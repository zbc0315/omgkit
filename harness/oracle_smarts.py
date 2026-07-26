#!/usr/bin/env python3
"""生成 SMARTS 解析的差分基准。

每行一条记录:`{"smarts": ..., "ok": bool, "na": int, "nb": int}`。
`ok` 为假时 `na`/`nb` 无意义。

语料 `corpus/smarts.txt` 由 `--build-corpus` 从 RDKit 自带的数据文件抽取
(PAINS 过滤器、官能团层级、RLewis 库),不手写 —— 手挑的 SMARTS 会不自觉地
只挑自己已经实现的写法。

用法:

    python3 harness/oracle_smarts.py --build-corpus     # 重建语料
    python3 harness/oracle_smarts.py --out harness/baseline/smarts.jsonl
"""

import argparse
import csv
import io
import json
import pathlib
import re
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = pathlib.Path(__file__).parent
CORPUS = HERE / "corpus" / "smarts.txt"


def build_corpus(rdkit_data: pathlib.Path) -> int:
    """从 RDKit 的数据文件里抽 SMARTS。"""
    out = []
    # 制表符分隔,取第一个"看起来像 SMARTS"的字段
    for name in [
        "Functional_Group_Hierarchy.txt",
        "FunctionalGroups.txt",
        "SmartsLib/RLewis_smarts.txt",
    ]:
        p = rdkit_data / name
        if not p.exists():
            print(f"跳过(不存在): {name}", file=sys.stderr)
            continue
        for line in p.read_text(errors="replace").splitlines():
            line = line.strip()
            if not line or line.startswith(("//", "#")):
                continue
            for tok in re.split(r"\t+", line):
                tok = tok.strip()
                if len(tok) > 1 and any(c in tok for c in "[]$#=~@") and " " not in tok:
                    out.append(tok)
                    break

    pains = rdkit_data / "Pains" / "wehi_pains.csv"
    if pains.exists():
        for row in csv.reader(io.StringIO(pains.read_text(errors="replace"))):
            if row and row[0].strip():
                out.append(row[0].strip())

    uniq = sorted({s for s in out if s and not s.startswith('"')})
    CORPUS.write_text("\n".join(uniq) + "\n")
    return len(uniq)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=pathlib.Path, help="输出 JSONL")
    ap.add_argument(
        "--build-corpus",
        action="store_true",
        help="从 RDKit 数据文件重建 corpus/smarts.txt",
    )
    ap.add_argument(
        "--rdkit-data",
        type=pathlib.Path,
        default=pathlib.Path("../rdkit/Data"),
        help="RDKit 源码的 Data 目录",
    )
    args = ap.parse_args()

    if args.build_corpus:
        n = build_corpus(args.rdkit_data)
        print(f"语料写出 {n} 条 → {CORPUS}")
        if not args.out:
            return 0

    if not args.out:
        ap.error("需要 --out 或 --build-corpus")

    lines = [s.strip() for s in CORPUS.read_text().splitlines()]
    records = []
    ok = 0
    for s in lines:
        if not s or s.startswith("#"):
            continue
        m = Chem.MolFromSmarts(s)
        if m is None:
            records.append({"smarts": s, "ok": False, "na": 0, "nb": 0})
        else:
            ok += 1
            records.append(
                {"smarts": s, "ok": True, "na": m.GetNumAtoms(), "nb": m.GetNumBonds()}
            )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as f:
        for r in records:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")
    print(f"已写出 {args.out}")
    print(f"  记录: {len(records)}(可解析 {ok},不可解析 {len(records) - ok})")
    print(f"  RDKit: {Chem.rdBase.rdkitVersion}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
