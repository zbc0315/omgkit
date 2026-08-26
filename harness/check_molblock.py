#!/usr/bin/env python3
"""我们写的 V2000 molblock,外部实现读回来是不是同一个分子。

# 这条判据覆盖的是**文件格式本身**

`verify_stereo.py` 走的是 JSONL:原子表与坐标以数组交出去,判官在 Python 里
自己拼一个 RDKit 分子。那条路把文件格式整个绕开了 —— 计数行、原子块的价键
字段、`M CHG` / `M ISO` / `M RAD`、键块的字段宽度,全都没有人读过。

而 `.mol` / `.sdf` 交出去之后,别人读的正是那些字段。写错一个,对方拿到的就是
另一个分子:少写价键字段,`[CH]` 读回来成 `[CH3]`;漏 `M RAD`,自由基读成甲基;
计数行挤出格,整条记录直接读不了。

# 判什么

读回来的分子必须**满足输入 SMILES 指定的每一处立体**。用带手性的子结构匹配,
不是比规范串 —— 三维结构必然给桥头碳之类定一个构型,而输入 SMILES 往往不写。
多出来的立体信息不是错。这一条与 `verify_stereo.py` 同一个口径。

外部实现自己都读不回的那一档(实测是三配位磷)单独计数,不算我们摆错 ——
判官先证明自己看得见,它报的 0 才有意义。

用法:python3 harness/check_molblock.py <blocks.sdf> [--min-checked N]
"""
import argparse
import sys

import rdkit
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


def records(path):
    """切成 `(语料行号, SMILES, molblock)`。"""
    smi = None
    lineno = None
    buf = []
    for line in open(path, encoding="utf-8"):
        if line.startswith(">>> "):
            lineno, smi = line[4:].rstrip("\n").split("\t", 1)
            buf = []
        elif line.rstrip("\n") == "$$$$":
            yield lineno, smi, "".join(buf)
            buf = []
        else:
            buf.append(line)


def rdkit_can_read_itself(smi):
    """外部实现拿**它自己**嵌的构象走同一条检验,读得回来吗。"""
    from rdkit.Chem import AllChem

    ref = Chem.MolFromSmiles(smi)
    if ref is None:
        return False
    mh = Chem.AddHs(ref)
    if AllChem.EmbedMolecule(mh, randomSeed=0xF00D) != 0:
        return False
    Chem.AssignStereochemistryFrom3D(mh)
    return Chem.RemoveHs(mh).HasSubstructMatch(Chem.RemoveHs(ref), useChirality=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("blocks")
    ap.add_argument("--min-checked", type=int, default=500)
    args = ap.parse_args()
    print(f"外部实现:RDKit {rdkit.__version__}")

    ok = bad = unreadable = blind = skipped = 0
    failures = []
    for lineno, smi, block in records(args.blocks):
        ref = Chem.MolFromSmiles(smi)
        if ref is None:
            skipped += 1
            continue
        # **removeHs=False**:删氢会改原子数,而立体正是靠坐标读的
        m = Chem.MolFromMolBlock(block, removeHs=False)
        if m is None:
            unreadable += 1
            failures.append(f"第 {lineno} 行 {smi}:外部实现读不了这条 molblock")
            continue
        try:
            Chem.AssignStereochemistryFrom3D(m)
            got = Chem.RemoveHs(m)
        except Exception as e:  # noqa: BLE001
            unreadable += 1
            failures.append(f"第 {lineno} 行 {smi}:读回来净化不了({e})")
            continue
        if got.HasSubstructMatch(Chem.RemoveHs(ref), useChirality=True):
            ok += 1
        elif not rdkit_can_read_itself(smi):
            blind += 1
        else:
            bad += 1
            failures.append(
                f"第 {lineno} 行 {smi}:读回来是 {Chem.MolToSmiles(got)},"
                "与输入指定的立体不符"
            )

    checked = ok + bad
    print(f"读回来一致 {ok};不符 {bad};读不了 {unreadable};"
          f"判官够不着 {blind};输入本身跳过 {skipped}")
    if failures:
        for f in failures[:8]:
            print(f"  ✗ {f}")
    if bad or unreadable:
        print("\n我们写的 molblock 交出去之后成了另一个分子。")
        return 1
    if checked < args.min_checked:
        print(f"\n只比过 {checked} 条,低于下限 {args.min_checked} —— 判据被喂空了")
        return 1
    print("\n写出的 molblock 读回来逐条一致。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
