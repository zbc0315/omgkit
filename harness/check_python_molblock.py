#!/usr/bin/env python3
"""读 molblock 的 **Python 绑定**:与 Rust 侧逐字符相同。

# 这条判据守什么

绑定那一层的规矩是"只做翻译,不做化学"。可读 molblock 这条路上,翻译要多做
两步 —— 净化,以及净化之后回来打立体标记。漏掉后一步不会报错,只会**静默地
把整个文件的立体丢掉**:分子照样合法、原子数照样对,只是 `@` 和 `/` 全没了。

所以这里不另造真值,直接比:同一批 molblock,`omgkit-io/examples/read_molblock`
走 Rust、`omgkit.parse_molblock` 走 Python,**写回来的规范 SMILES 必须逐字符
相同**。Rust 那一侧已经由 `check_molblock_read.py` 与外部实现比过,这一侧与它
相同,那条外部判据就继承了过来。

# 立体那一档要单独有个下限

"逐字符相同"在**两边都把立体丢光**时同样成立 —— 那正是这条判据要抓的故障,
而它会以全绿的样子出现。所以比过的串里带立体符号(`@` 或 `/`)的条数配一条
下限:低于它就说明这一档被喂空了,零分歧说明不了任何事。

用法:

    python3 harness/check_molblock_read.py --write  /tmp/in.sdf harness/corpus/large.smi
    cargo run -q -p omgkit-io --release --example read_molblock -- /tmp/in.sdf > /tmp/ours.txt
    python3 harness/check_python_molblock.py /tmp/in.sdf /tmp/ours.txt
"""
import argparse
import sys

import omgkit

# 比过的串里**带立体符号**的条数下限。实测 641 条 —— 贴着现值留了余量。
MIN_WITH_STEREO = 600


def blocks(path):
    """`>>> 行号\\t原始 SMILES` + molblock + `$$$$` 的装法,与 Rust 那侧同一个。"""
    smi = None
    lineno = None
    buf = []
    for line in open(path, encoding="utf-8"):
        if line.startswith(">>> "):
            lineno, smi = line[4:].rstrip("\n").split("\t", 1)
            buf = []
        elif line.rstrip("\n") == "$$$$":
            yield lineno, smi, "".join(buf)
        else:
            buf.append(line)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("sdf", help="molblock 输入,与 Rust 侧读的是同一份字节")
    ap.add_argument("ours", help="examples/read_molblock 的输出")
    ap.add_argument("--min-checked", type=int, default=8000)
    args = ap.parse_args()

    print(f"  omgkit wheel:{omgkit.__file__}")
    rust = {}
    for line in open(args.ours, encoding="utf-8"):
        lineno, _smi, got = line.rstrip("\n").split("\t", 2)
        rust[lineno] = got

    checked = with_stereo = refused_both = 0
    failures = []
    for lineno, smi, block in blocks(args.sdf):
        want = rust.get(lineno)
        if want is None:
            failures.append(f"第 {lineno} 行 {smi}:Rust 侧一行输出都没有")
            continue
        try:
            rec = omgkit.parse_molblock(block)
            # **坐标必须与净化之后的原子表一一对应。** 读进来之后要先净化才谈得上
            # 打立体标记,而净化那一步万一动了原子表(加原子、重排下标),坐标就
            # 与分子错位了 —— 分子照样合法、条数照样对,只有几何整个搬了家,
            # 而且一声不响。这一条在全语料上钉住,不只在一两个分子上试过。
            if len(rec.coords) != rec.mol.num_atoms:
                failures.append(
                    f"第 {lineno} 行 {smi}:坐标 {len(rec.coords)} 条,原子 "
                    f"{rec.mol.num_atoms} 个 —— 净化动过原子表"
                )
                continue
            got = rec.mol.to_canonical_smiles()
        except ValueError as e:
            # Rust 那侧把失败打印成 `<读不了:…>` / `<净化不了:…>`。两侧都拒收
            # 才算一致 —— 一侧读得出一侧读不出,正是翻译丢了东西的样子。
            if want.startswith("<"):
                refused_both += 1
            else:
                failures.append(f"第 {lineno} 行 {smi}:Rust 读得出 `{want}`,Python 抛了 {e}")
            continue
        if want.startswith("<"):
            failures.append(f"第 {lineno} 行 {smi}:Rust 拒收({want}),Python 读成了 `{got}`")
            continue
        if got != want:
            failures.append(f"第 {lineno} 行 {smi}:Python `{got}` ≠ Rust `{want}`")
            continue
        checked += 1
        if "@" in got or "/" in got or "\\" in got:
            with_stereo += 1

    print(f"逐字符相同 {checked} 条;两侧都拒收 {refused_both} 条;不一致 {len(failures)} 条")
    print(f"  其中带立体符号的 {with_stereo} 条(下限 {MIN_WITH_STEREO})")
    for f in failures[:8]:
        print(f"  ✗ {f}")
    if failures:
        print("\nPython 绑定读出来的不是 Rust 读出来的那个分子。")
        return 1
    if with_stereo < MIN_WITH_STEREO:
        print(f"\n带立体的只有 {with_stereo} 条,低于下限 {MIN_WITH_STEREO} —— "
              "立体那一档被喂空了,两侧一起丢光也是这个样子")
        return 1
    if checked < args.min_checked:
        print(f"\n只比过 {checked} 条,低于下限 {args.min_checked} —— 判据被喂空了")
        return 1
    print("\n与 Rust 侧逐字符相同。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
