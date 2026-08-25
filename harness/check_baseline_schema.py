#!/usr/bin/env python3
"""基准文件与生成它的脚本**脱钩**了没有 —— 只比结构,不比数值。

# 为什么需要这条

提交 `61b8d58` 教会了 `dump_chirality.py` 收**三配位立体中心**(亚砜的 S、
膦的 P),却**没有重导** `harness/baseline/smoke.chirality.jsonl`。于是那个
提交声称落地的那一档,在主手性判官眼里根本不存在:committed 的基准里 247 个
中心全是四配位、连 `three_coordinate` 这个字段都没有,而当时的脚本导出来是
248 个、其中 8 个三配位。

**四个月里没有任何判据看得见这件事。** 判官照跑照绿 —— 它判的是基准里写着的
东西,而基准里那一档已经不写了。分母悄悄变小,没有任何红。

同一个坑第二次:`dump_chirality.py` 一直不导键的顺反列,于是这 150 个分子里
23 个带 E/Z 的分子,在两条手性判官那边是**没有顺反**的分子。

# 为什么只比结构,不比数值

这些基准的**数值**钉在各自的 RDKit 版本上(`smoke.chirality.jsonl` 是
2022.09.5,`smoke.bounds.jsonl` 是 2025.09.2 —— 各自的文件头写着为什么),
而 CI 只装得下一个版本。逐字节比会天天红,那种判据活不过一周。

**结构与版本无关**:记录有哪些键、`centers` 元素有哪些键、`bonds` 元组几列。
脚本长出新字段(或者删掉一个),这里当场看得见,而这正是上面两次漏掉的东西。

# 覆盖不到的:`dump_gram.py`

它读的是 `harness/baseline/rdkit_bounds.jsonl`(7.7 MB,没入库),CI 里跑不了。

用法:

    python3 harness/check_baseline_schema.py
"""

import json
import pathlib
import subprocess
import sys
import tempfile

import rdkit

REPO = pathlib.Path(__file__).resolve().parent.parent

# (基准, 生成它的脚本, 对照用的语料, 对照用的 limit, 最少要导出几条, 重导命令)
#
# 倒数第二个是**分母闸**:生成器一条都没导出来的话,两边的结构都是空的,
# "结构一致"就成了空过。喂空的判据不许打印"全部通过"。
#
# 最后那个是**红了之后照着敲的那一行**。判据报"该重导了"却不说怎么重导,
# 下一个人就得去翻三个文件 —— 而重导命令散在各个 oracle 的 `//!` 注释里。
# 注意重导 `smoke.chirality` / `smoke.lonepair` 要用 **2022.09.5**
# (`dump_chirality.py` 的文件头写着为什么与界基准的 2025.09.2 分开)。
PAIRS = [
    (
        "smoke.chirality.jsonl",
        "dump_chirality.py",
        "large.smi",
        "20",
        3,
        ".venv/bin/python harness/dump_chirality.py harness/corpus/large.smi"
        " harness/baseline/smoke.chirality.jsonl 150",
    ),
    (
        "smoke.lonepair.jsonl",
        "dump_chirality.py",
        "lonepair.smi",
        None,
        5,
        ".venv/bin/python harness/dump_chirality.py harness/corpus/lonepair.smi"
        " harness/baseline/smoke.lonepair.jsonl",
    ),
    (
        "smoke.bounds.jsonl",
        "dump_bounds.py",
        "large.smi",
        "20",
        3,
        # 冒烟档是**从全量档直接切**的(每 15 个取 1),不重新跑一遍 RDKit ——
        # 免得冒烟档与全量档撞上不同版本(见 `a6235ed`)。全量档 7.7 M,不入库。
        "python3 harness/dump_bounds.py harness/corpus/large.smi"
        " harness/baseline/rdkit_bounds.jsonl 400 && "
        "awk 'NR % 15 == 1' harness/baseline/rdkit_bounds.jsonl"
        " > harness/baseline/smoke.bounds.jsonl",
    ),
]


def schema(path):
    """一份 jsonl 的结构:记录键、`centers` 元素键、`bonds` 元组的列数。

    键取**并集与交集两样**,不取第一条。只比并集的话,"字段在一部分记录上丢了"
    是看不见的 —— 变异实测:把某**一个**中心的 `three_coordinate` 删掉,
    只比并集的版本退 0。交集把这一档收进来:并集说"哪里出现过",
    交集说"是不是处处都有",两个都对上才算结构一致。

    `bonds` 的列数本来就是集合,少一列当场就是另一个集合。
    """
    rec_u, cen_u, cols, n = set(), set(), set(), 0
    rec_i, cen_i = None, None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            n += 1
            rec_u |= set(r)
            rec_i = set(r) if rec_i is None else rec_i & set(r)
            for c in r.get("centers", []):
                cen_u |= set(c)
                cen_i = set(c) if cen_i is None else cen_i & set(c)
            for b in r.get("bonds", []):
                cols.add(len(b))
    return n, (
        frozenset(rec_u),
        frozenset(rec_i or ()),
        frozenset(cen_u),
        frozenset(cen_i or ()),
        frozenset(cols),
    )


def main() -> int:
    # 版本要打出来。结构与版本无关是这条判据成立的前提,不是可以默认的事。
    print(f"外部实现:RDKit {rdkit.__version__}(这条判据只比结构,与版本无关)")
    bad = 0
    with tempfile.TemporaryDirectory() as tmp:
        for name, script, corpus, limit, need, redo in PAIRS:
            base = REPO / "harness" / "baseline" / name
            fresh = pathlib.Path(tmp) / name
            cmd = [
                sys.executable,
                str(REPO / "harness" / script),
                str(REPO / "harness" / "corpus" / corpus),
                str(fresh),
            ]
            if limit is not None:
                cmd.append(limit)
            r = subprocess.run(cmd, capture_output=True, text=True, check=False)
            if r.returncode != 0 or not fresh.exists():
                print(f"✗ {name}:{script} 跑不起来(退出码 {r.returncode})")
                print(r.stdout[-800:] + r.stderr[-800:])
                bad += 1
                continue
            n_new, s_new = schema(fresh)
            n_old, s_old = schema(base)
            if n_new < need:
                print(f"✗ {name}:{script} 只导出 {n_new} 条(至少要 {need})—— 判据没东西可比")
                bad += 1
                print(f"    重导:{redo}")
                continue
            if s_new != s_old:
                labels = (
                    "记录字段(并集)",
                    "记录字段(处处都有的)",
                    "centers 字段(并集)",
                    "centers 字段(处处都有的)",
                    "bonds 列数",
                )
                print(f"✗ {name}:结构与 {script} 现在导出来的**不一样** —— 基准该重导了")
                for lab, a, b in zip(labels, s_old, s_new):
                    if a != b:
                        print(f"    {lab}:基准 {sorted(a)} vs 脚本 {sorted(b)}")
                print(f"    重导:{redo}")
                bad += 1
                continue
            print(f"✓ {name}:结构与 {script} 一致({n_old} 条基准,现导 {n_new} 条对照)")
    if bad:
        print(f"\n{bad} 份基准与生成它的脚本脱钩了 —— 重导命令在上面每一条的后面。")
        return 1
    print("\n全部基准与生成它的脚本结构一致。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
