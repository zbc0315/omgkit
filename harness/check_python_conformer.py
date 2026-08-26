#!/usr/bin/env python3
"""构型生成的 **Python 绑定**:与 Rust 侧逐位相同。

# 这条判据守什么

绑定那一层的规矩是"只做翻译,不做化学"。可"翻译"本身也会错:原子表对不上、
坐标顺序错位、补氢那一步漏掉 —— 每一种都让 Python 用户拿到一组**看着正常
其实是另一个分子**的坐标,而 Rust 侧的全套判据一概盖不到。

所以这里不另造真值,直接比:同一条 SMILES,`examples/dump_conformers` 走
Rust、`Mol.conformer()` 走 Python,**原子序数、形式电荷、键表、坐标必须逐位
相同**。两侧调的是同一个 `pipeline::conformer_for`,全程无随机数,所以"逐位
相同"是可以要求的 —— 差一位就说明翻译丢了东西。

顺带也就把 `verify_stereo.py` 那条外部判据继承过来了:Rust 那侧的坐标已经交给
RDKit 从三维读回过立体化学,Python 这侧与它逐位相同,自然同样成立。

# 分母闸

真正比过的条数低于 `--min-checked` 就直接失败。语料换了、上游筛选变了都可能
让这一档悄悄变空,而"零分歧"在那时依然成立 —— 那是最会骗人的一种绿。

用法:

    cargo run -q -p omgkit-conf --release --example dump_conformers -- \\
        harness/corpus/large.smi > /tmp/ours.jsonl
    python3 harness/check_python_conformer.py /tmp/ours.jsonl
"""
import argparse
import json
import sys

import omgkit


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("jsonl", help="dump_conformers 的输出")
    ap.add_argument("--min-checked", type=int, default=500)
    args = ap.parse_args()

    print(f"  omgkit wheel:{omgkit.__file__}")
    checked = 0
    failed = []
    unreadable = 0
    for line in open(args.jsonl, encoding="utf-8"):
        r = json.loads(line)
        smi = r["smiles"]
        try:
            conf = omgkit.parse_smiles(smi).conformer()
        except ValueError as e:
            unreadable += 1
            failed.append(f"{smi}:Rust 侧出了构型,Python 侧抛了 {e}")
            continue
        got = {
            "z": conf.mol.atomic_nums,
            "charge": conf.mol.formal_charges,
            "bonds": [[i, j] for i, j, _ in conf.mol.bonds],
            "xyz": [list(p) for p in conf.coords],
        }
        want = {
            "z": r["z"],
            "charge": r["charge"],
            "bonds": [[i, j] for i, j, _ in r["bonds"]],
            "xyz": r["xyz"],
        }
        checked += 1
        for key in ("z", "charge", "bonds", "xyz"):
            if got[key] != want[key]:
                failed.append(
                    f"{smi}:`{key}` 两侧不同"
                    f"(Rust {len(want[key])} 项 / Python {len(got[key])} 项)"
                )
                break

    print(f"逐位比过 {checked} 条;两侧不同 {len(failed)} 条;Python 侧生不出来 {unreadable} 条")
    if failed:
        for f in failed[:8]:
            print(f"  ✗ {f}")
        print(f"\n绑定与 Rust 侧对不上 —— 翻译丢了东西。")
        return 1
    if checked < args.min_checked:
        print(f"\n只比过 {checked} 条,低于下限 {args.min_checked} —— 判据被喂空了")
        return 1
    print("绑定与 Rust 侧逐位相同。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
