#!/usr/bin/env python3
"""我们写的 V2000 molblock,外部实现读回来是不是同一个分子。

# 这条判据覆盖的是**文件格式本身**

`verify_stereo.py` 走的是 JSONL:原子表与坐标以数组交出去,判官在 Python 里
自己拼一个 RDKit 分子。那条路把文件格式整个绕开了 —— 计数行、原子块的价键
字段、`M CHG` / `M ISO` / `M RAD`、键块的字段宽度,全都没有人读过。

而 `.mol` / `.sdf` 交出去之后,别人读的正是那些字段。写错一个,对方拿到的就是
另一个分子:少写价键字段,`[CH]` 读回来成 `[CH3]`;漏 `M RAD`,自由基读成甲基;
计数行挤出格,整条记录直接读不了。

# 读的是**真 SDF**,记录边界与数据字段一并验了

先前我方每条前面加一行 `>>> 行号\t原串`、后面手写 `$$$$`,判官照那个自造的
包装切。于是 `$$$$` 摆在哪、数据字段怎么写,**外部实现一次都没读过** ——
而那正是 `.sdf` 交出去之后别人要读的东西。

现在我方走 `write_sdf_record` 导真 SDF(行号与原串是两个数据字段),判官用
`ForwardSDMolSupplier` 当普通 SDF 读。**用 Forward 而不是 `SDMolSupplier`**:
后者会跳过读不了的记录,条数就对不上了 —— 而条数正是要比的东西之一。

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


def raw_blocks(path):
    """按 `$$$$` 切出每条记录的原文。只用来数条数、以及给失败的那条留个说法。"""
    buf = []
    for line in open(path, encoding="utf-8"):
        if line.rstrip("\n") == "$$$$":
            yield "".join(buf)
            buf = []
        else:
            buf.append(line)


def field(block, name):
    """从记录原文里抠一个数据字段的值。读不了的记录也拿得到,只为报错好看。"""
    lines = block.split("\n")
    for i, ln in enumerate(lines):
        if ln.startswith(">") and f"<{name}>" in ln and i + 1 < len(lines):
            return lines[i + 1]
    return "?"


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

    blocks = list(raw_blocks(args.blocks))
    # **removeHs=False**:删氢会改原子数,而立体正是靠坐标读的
    with open(args.blocks, "rb") as f:
        mols = list(Chem.ForwardSDMolSupplier(f, sanitize=True, removeHs=False))
    if len(mols) != len(blocks):
        print(f"条数不同:外部实现读出 {len(mols)} 条,文件里有 {len(blocks)} 条 —— "
              "我方写的记录边界(`$$$$`)与它切出来的对不上")
        return 1

    ok = bad = unreadable = blind = skipped = missing_field = 0
    failures = []
    for block, m in zip(blocks, mols):
        # 身份从**数据字段**里取:那一步也在验我方写的字段读不读得回来
        smi = m.GetProp("原串") if m is not None and m.HasProp("原串") else field(block, "原串")
        lineno = m.GetProp("行号") if m is not None and m.HasProp("行号") else field(block, "行号")
        if m is not None and not (m.HasProp("原串") and m.HasProp("行号")):
            missing_field += 1
            failures.append(f"第 {lineno} 行 {smi}:数据字段没读回来")
            continue
        ref = Chem.MolFromSmiles(smi)
        if ref is None:
            skipped += 1
            continue
        if m is None:
            unreadable += 1
            failures.append(f"第 {lineno} 行 {smi}:外部实现读不了这条记录")
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
          f"判官够不着 {blind};输入本身跳过 {skipped};数据字段丢了 {missing_field}")
    if failures:
        for f in failures[:8]:
            print(f"  ✗ {f}")
    if bad or unreadable or missing_field:
        print("\n我们写的 molblock 交出去之后成了另一个分子。")
        return 1
    if checked < args.min_checked:
        print(f"\n只比过 {checked} 条,低于下限 {args.min_checked} —— 判据被喂空了")
        return 1
    print("\n写出的 molblock 读回来逐条一致。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
