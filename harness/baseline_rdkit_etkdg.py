#!/usr/bin/env python3
"""参照实现(RDKit ETKDGv3)在同一份语料上的失败率与耗时。

# 这个脚本为什么必须能跑

参照实现的失败率是这个项目的头号参照:`omgkit-conf/examples/feasibility.rs`
里那条硬闸(`MAX_INFEASIBLE_FRAC`)就是照着它定的。一个**量不出来的参照**
等于一句传说 —— 先前这里的语料路径写死成一个早已删掉的 worktree 的绝对路径
(`.claude/worktrees/agent-…/harness/corpus/large.smi`),脚本一行都跑不了,
那个数从此没法复核。

**当前值:0.41%**(large.smi,ETKDGv3,RDKit 2025.09.2)。更早的版本上量到的是
0.52%,仓库里若还有那个数,那是历史值 —— 见下面这一条。

# 这个数**跟 RDKit 版本走**

ETKDG 每个版本都在改。所以这里把版本打在最前面,引用这个数的地方也要连版本
一起写。仓库钉的是 `harness/requirements.lock` 里那一个。

用法:

    python3 harness/baseline_rdkit_etkdg.py [语料.smi]

口径与 `measure_params.py` 一致:单进程、`AddHs`、ETKDGv3、种子 0xf00d。
"""
import pathlib
import sys
import time

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

DEFAULT_CORPUS = pathlib.Path(__file__).resolve().parent / "corpus" / "large.smi"


def main() -> int:
    corpus = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_CORPUS
    if not corpus.is_file():
        print(f"语料不在:{corpus}", file=sys.stderr)
        return 1
    print(f"外部实现:RDKit {Chem.rdBase.rdkitVersion}")
    print(f"语料:{corpus}")

    smis = []
    for line in corpus.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            smis.append(line.split("\t")[0])
    print(f"语料行(去注释):{len(smis)}", flush=True)

    fail = exc = parse_fail = 0
    times = []
    worst = []
    t0 = time.time()
    for i, smi in enumerate(smis):
        m = Chem.MolFromSmiles(smi)
        if m is None:
            parse_fail += 1
            continue
        mh = Chem.AddHs(m)
        p = AllChem.ETKDGv3()
        p.randomSeed = 0xF00D
        t = time.time()
        try:
            if AllChem.EmbedMolecule(mh, p) < 0:
                fail += 1
        except Exception as e:  # noqa: BLE001
            exc += 1
            fail += 1
            print(f"  第 {i} 条抛异常:{smi[:70]}\n     {str(e)[:200]}", flush=True)
        times.append(time.time() - t)
        worst.append((times[-1], smi))
    tot = time.time() - t0

    if not times:
        print("一个分子都没嵌 —— 语料是空的?", file=sys.stderr)
        return 1
    times.sort()
    mean_ms = 1000 * sum(times) / len(times)
    print(
        f"墙钟合计 {tot:.1f} s;嵌入 {len(times)} 个;"
        f"平均 {mean_ms:.1f} ms;中位 {1000 * times[len(times) // 2]:.1f} ms"
    )
    print(f"p99 {1000 * times[int(0.99 * len(times))]:.0f} ms;最慢 {times[-1]:.2f} s")
    print(
        f"解析失败 {parse_fail};嵌入失败 {fail}"
        f"({100.0 * fail / len(times):.2f}%);C++ 异常 {exc}"
    )
    worst.sort(key=lambda w: -w[0])
    print("最慢的 5 个:")
    for dt, smi in worst[:5]:
        print(f"   {dt:.2f} s  {smi[:80]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
