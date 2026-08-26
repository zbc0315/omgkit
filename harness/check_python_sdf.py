#!/usr/bin/env python3
"""读 SDF 的 **Python 绑定**:与 Rust 侧逐条相同。

# 这条判据守什么

绑定那一层"只做翻译",可这条路上翻译要多做的事不少:切分的结果要一条不落地
搬过来、数据字段要保序、**读不了的那条要留在它自己的位置上**。最后一条最容易
出错,而错法特别安静 —— 静默跳过之后条数变小,调用方数出来的与文件里的不符,
没有任何地方报错。

所以不另造真值,直接比:同一份 SDF,`examples/read_sdf` 走 Rust、
`omgkit.read_sdf` 走 Python,**条数、每条的规范 SMILES、每条的数据字段**都要
一样。Rust 那侧已经由 `check_sdf.py` 与外部实现比过,这一侧与它相同,那条
外部判据就继承了过来。

# 读不了的那条也要对上

不只是"两侧都失败"就算数:失败**落在第几条**必须一样。整个文件错位一条的话,
每条都还是"读得出",只是配错了对象 —— 而那种错在只比集合的判据下完全看不见。

用法:

    python3 harness/check_sdf.py --write /tmp/data.sdf harness/corpus/large.smi
    cargo run -q -p omgkit-io --release --example read_sdf -- /tmp/data.sdf > /tmp/ours.txt
    python3 harness/check_python_sdf.py /tmp/data.sdf /tmp/ours.txt
"""
import argparse
import json
import sys

import omgkit

# 带数据字段的条数下限。字段那一档一旦被喂空,"逐条相同"照样成立。
MIN_WITH_DATA = 8000

# 读不了的条数下限。**这一档是这条判据的重点** —— 语料里有一条(二茂铁,
# 外部实现写成了 V3000)。一条都没有的话,"坏记录留在原位"就没人验过了。
MIN_FAILED = 1


def rust_records(path):
    """Rust 侧每条:`(第几条, SMILES 或 <…>, 数据字段)`。"""
    for line in open(path, encoding="utf-8"):
        idx, smi, data = line.rstrip("\n").split("\t", 2)
        yield int(idx), smi, [(k, v) for k, v in json.loads(data)] if data else []


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("sdf")
    ap.add_argument("ours", help="examples/read_sdf 的输出")
    ap.add_argument("--min-checked", type=int, default=8000)
    args = ap.parse_args()

    print(f"  omgkit wheel:{omgkit.__file__}")
    rust = list(rust_records(args.ours))
    mine = omgkit.read_sdf(open(args.sdf, encoding="utf-8").read())
    if len(mine) != len(rust):
        print(f"条数不同:Python {len(mine)},Rust {len(rust)} —— 切分对不上")
        return 1

    checked = with_data = failed = 0
    failures = []
    for (idx, want_smi, want_data), rec in zip(rust, mine):
        rust_failed = want_smi.startswith("<")
        if (rec.error is not None) != rust_failed:
            failures.append(
                f"第 {idx} 条:Rust {'读不了' if rust_failed else '读得出'},"
                f"Python {'读不了' if rec.error else '读得出'}({rec.error or want_smi})"
            )
            continue
        if rust_failed:
            failed += 1
            continue
        got = rec.block.mol.to_canonical_smiles()
        if got != want_smi:
            failures.append(f"第 {idx} 条:Python `{got}` ≠ Rust `{want_smi}`")
            continue
        got_data = [tuple(kv) for kv in rec.data]
        if got_data != want_data:
            failures.append(f"第 {idx} 条:数据字段不同 —— Python {got_data},Rust {want_data}")
            continue
        # 坐标要与净化之后的原子表对齐,与单条那条判据同一条不变式
        if len(rec.block.coords) != rec.block.mol.num_atoms:
            failures.append(
                f"第 {idx} 条:坐标 {len(rec.block.coords)} 条,原子 "
                f"{rec.block.mol.num_atoms} 个"
            )
            continue
        checked += 1
        if got_data:
            with_data += 1

    print(f"逐条相同 {checked};两侧都读不了 {failed}(下限 {MIN_FAILED});"
          f"不一致 {len(failures)}")
    print(f"  带数据字段的 {with_data} 条(下限 {MIN_WITH_DATA})")
    for f in failures[:8]:
        print(f"  ✗ {f}")
    if failures:
        print("\nPython 绑定读出来的不是 Rust 读出来的那一批。")
        return 1
    if failed < MIN_FAILED:
        print(f"\n读不了的只有 {failed} 条,低于下限 {MIN_FAILED} —— "
              "\"坏记录留在原位\"这一档没人验过")
        return 1
    if with_data < MIN_WITH_DATA:
        print(f"\n带数据字段的只有 {with_data} 条,低于下限 {MIN_WITH_DATA} —— 字段那一档被喂空了")
        return 1
    if checked < args.min_checked:
        print(f"\n只比过 {checked} 条,低于下限 {args.min_checked} —— 判据被喂空了")
        return 1
    print("\n与 Rust 侧逐条相同。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
