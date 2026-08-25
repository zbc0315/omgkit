"""**分母闸** —— 判官共用这一份。

# 判官不可能自己发现自己少比了几个分子

`check_*.py` 都是"数分歧,分歧为 0 就退 0"。分歧数是**分子**,而分母是上游
喂的:`dump_*` 少喂几个分子,每一档都跟着变好看,判据照样退 0。喂**空文件**
进去最彻底 —— 一条分歧都没有,打印一片空白然后"全部通过"。

这不是假想:`check_write.py` 的注释里记着,先前那版拿空文件进去就是这么过的。

# 为什么只留一份

先前 `check_write.py` 自己写了一套,别的判官一套都没有。四份判官各抄一遍必然
分岔,而且是静默分岔(`omgkit-conf` 那边四个判官各抄一遍连接表重建,三份漏了
净化、一份漏了顺反列,四个月没人发现)。所以这里只留一份。

# 两个方向都要闸

`uncompared = n_corpus - compared` 传错语料时会是**负数**,而负数当然 <= 上限
—— 分母闸就静默失效了。实测:拿大语料的 TSV 配冒烟语料,先前那版打印
"没写出 -8690 行"然后退 0。所以 `rows > n_corpus` 单独判一次。
"""

from __future__ import annotations

import pathlib


def corpus_size(path: pathlib.Path) -> int:
    """语料的有效行数。**`#` 开头是注释、空行忽略** —— 不认注释的话,
    注释行会被算进分母,闸就永远差那么几行。"""
    return sum(
        1
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    )


def verdict(n_corpus: int, rows: int, compared: int, cap: int) -> str | None:
    """分母核得上返回 `None`,核不上返回该打印的理由。

    - `n_corpus`:喂给 dump 的语料有多少行
    - `rows`:dump 出来多少行
    - `compared`:真正进了比对的有多少条
    - `cap`:允许有多少条没进比对(**分母闸,不是宽容度**)
    """
    if rows > n_corpus:
        return (
            f"dump 出的行数({rows})比语料还多({n_corpus})—— 十有八九是 TSV 与语料"
            "对不上(传错文件了)。分母核不了,这条判据算出来的数没有意义"
        )
    uncompared = n_corpus - compared
    if uncompared > cap:
        return (
            f"语料里有 {uncompared} 条没真正被比对,超过上限 {cap} —— "
            "分歧数是分子,覆盖面是分母。别调大这个数,先查是哪一类分子进不了比对"
        )
    return None


def line(n_corpus: int, rows: int, compared: int, cap: int) -> str:
    """判据末尾那行分母账,四条判官打印同一个格式。"""
    return (
        f"  语料 {n_corpus} 行,dump 出 {rows} 行,真正比对 {compared} 条,"
        f"没比到 {n_corpus - compared} 条(上限 {cap})"
    )
