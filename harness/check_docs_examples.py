#!/usr/bin/env python3
"""判官:文档站里的 Python 示例真的跑得出文档写的结果。

# 为什么需要它

`gates.sh` 的四十来道闸**一条都不读 `docs/`**,`mkdocs build --strict` 只查链接
不查内容。于是代码往前走、文档留在原地:2026-08-30 的文档审查一次找出四处
"示例输出与实测不符",其中两处是**凭空编的报错回显**(文档写英文
`ValueError: unclosed branch`,而实现报的是中文"括号不匹配",整个仓库都没有
`unclosed branch` 这个串)。

这条闸把 `docs/**/*.md`、两份 README 里每一个 ```pycon 块当 doctest 跑。

# 判据的分母也要守

只跑到几个块就"全部通过"是最坏的一种绿。下限写死在 `MIN_BLOCKS` /
`MIN_STATEMENTS` 里,少于它当场红 —— 文档被搬走、glob 写错、解析器悄悄
返回空列表,都会先撞上这一条。

# 省略号

浮点数、路径、内存地址这类不该逐位钉的东西,文档里写 `...`,这里开
`doctest.ELLIPSIS`。`IGNORE_EXCEPTION_DETAIL` **不开** —— 报错消息的正文正是
这条闸最该守的东西。

用法:
    python3 harness/check_docs_examples.py            # 仓库根下全找
    python3 harness/check_docs_examples.py docs/guide # 只跑一个子树
"""

from __future__ import annotations

import doctest
import os
import pathlib
import re
import sys
import tempfile

# 少于这个数就说明"文档被搬走了/glob 写错了",不是"文档里本来就没例子"
MIN_BLOCKS = 12
MIN_STATEMENTS = 40

# ```python 与 ```pycon 都要按**出现顺序**过一遍:前者常常是 `import omgkit`
# 这类铺垫,只在后者里跑的话,除了第一块以外全会 `NameError`。
#
# **缩进的围栏也要收。** admonition(`!!! tip`)里的块要缩进四格,只匹配顶格的
# 围栏会把它们整批漏掉 —— 而那里面躺着好几段真正的示例输出。收进来之后按围栏
# 自己的缩进量去掉前缀,再交给 doctest。
FENCE = re.compile(
    r"^([ \t]*)```(python|pycon)[ \t]*$(.*?)^\1```[ \t]*$", re.MULTILINE | re.DOTALL
)

# 这些文件里的块是**写给人看的片段**,不是可执行示例(缺上下文、故意省略)。
# 名单必须逐条给理由 —— 不给理由的豁免是判据上的一个洞。
SKIP: dict[str, str] = {}


def mark_blank_lines(body: str) -> str:
    """把**预期输出内部**的空行换成 doctest 的 `<BLANKLINE>`。

    doctest 拿空行当"这一段预期输出到此为止"。而 molblock 的第三行(注释行)
    按规范就是空的,文档照实贴出来,doctest 就只比到第二行为止 —— 底下十几行
    连同那个改错了的键行**一起不在判据视野里**。

    规则:块里的一个空行,只有当它后面第一个非空行**不是** `>>>` 时,才算落在
    预期输出内部。真正用来分隔两个例子的空行后面跟的一定是下一个 `>>>`。
    这样文档里不必出现 `<BLANKLINE>` 这种只对判据有意义的记号。
    """
    lines = body.split("\n")
    out = []
    for i, ln in enumerate(lines):
        if ln.strip() == "":
            nxt = next((x for x in lines[i + 1 :] if x.strip() != ""), None)
            # 块尾的空行后面什么都没有 —— 它是围栏前的换行,不是预期输出的一部分。
            # 一并换成 `<BLANKLINE>` 的话,每个块的预期输出末尾都会凭空多一个空行。
            inside = nxt is not None and not nxt.lstrip().startswith(">>>")
            out.append("<BLANKLINE>" if inside else ln)
        else:
            out.append(ln)
    return "\n".join(out)


def blocks(path: pathlib.Path) -> list[tuple[int, str, str]]:
    """`(起始行号, 语言, 块内容)`,按出现顺序。行号是为了报错时指得回原文。"""
    text = path.read_text(encoding="utf-8")
    out = []
    for m in FENCE.finditer(text):
        line = text.count("\n", 0, m.start(3)) + 1
        indent = m.group(1)
        body = m.group(3)
        if indent:
            body = "\n".join(
                ln[len(indent) :] if ln.startswith(indent) else ln
                for ln in body.split("\n")
            )
        out.append((line, m.group(2), body))
    return out


def main() -> int:
    root = pathlib.Path(__file__).resolve().parent.parent
    targets = [pathlib.Path(a) for a in sys.argv[1:]] or [
        root / "docs",
        root / "README.md",
        root / "README.zh-CN.md",
    ]

    files: list[pathlib.Path] = []
    for t in targets:
        if t.is_dir():
            files.extend(sorted(t.rglob("*.md")))
        elif t.is_file():
            files.append(t)

    # **在临时目录里跑。** 文档里有好几个块往当前目录写文件(`out.sdf`、
    # `alanine.mol`)—— 那是给读者看的正常写法,但判据在仓库根跑就会把它们
    # 落进工作区。判据不许改它所在的仓库。
    #
    # 文档路径在 chdir **之前**就解析成绝对路径了(`files` 是从 `root` 拼的),
    # 所以换目录不影响读取。
    parser = doctest.DocTestParser()
    runner = doctest.DocTestRunner(
        optionflags=doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE,
        verbose=False,
    )

    import io

    n_blocks = 0
    n_statements = 0
    failed: list[str] = []
    tmp = tempfile.TemporaryDirectory()
    here = os.getcwd()
    os.chdir(tmp.name)
    for f in files:
        rel = f.relative_to(root)
        if str(rel) in SKIP:
            continue
        # **同一份文档里的块是接力的** —— 第一块 `import omgkit`,后面的块直接用。
        # 按块清空命名空间的话,除了第一块以外全会 `NameError`,而那不是文档的错。
        # **`import omgkit` 预先放进命名空间。** 指南各页的惯例是"入门那页展示
        # 一次导入,后面各页直接用",不是每页都重写一遍。判据按那个惯例来 ——
        # 反过来要求每页自带导入,红的会是文档的体例而不是文档的内容。
        globs: dict = {"omgkit": __import__("omgkit")}
        for line, lang, body in blocks(f):
            if lang == "python":
                # 铺垫块:跑得通就把名字留给后面的块用;跑不通不算这条闸的事
                # (它没有"文档写的结果"可比),但也不许它把命名空间弄脏。
                try:
                    exec(body, globs)  # noqa: S102 —— 判据自己就是要跑文档里的代码
                except Exception:  # noqa: BLE001
                    pass
                continue
            n_blocks += 1
            test = parser.get_doctest(
                mark_blank_lines(body), globs, str(rel), str(rel), line - 1
            )
            n_statements += len(test.examples)
            if not test.examples:
                continue
            buf = io.StringIO()
            res = runner.run(test, out=buf.write, clear_globs=False)
            # `get_doctest` 把 globs 复制了一份,跑完要合回来,下一块才接得上
            globs.update(test.globs)
            if res.failed:
                failed.append(f"{rel}:{line}\n{buf.getvalue().rstrip()}")

    os.chdir(here)
    tmp.cleanup()

    print(f"文档里的 pycon 块 {n_blocks} 个,可执行语句 {n_statements} 条")
    if n_blocks < MIN_BLOCKS or n_statements < MIN_STATEMENTS:
        print(
            f"只找到 {n_blocks} 个块 / {n_statements} 条语句,"
            f"低于下限 {MIN_BLOCKS} / {MIN_STATEMENTS} —— "
            "这条闸被喂空了,不是文档里没例子",
            file=sys.stderr,
        )
        return 1
    if failed:
        for f in failed:
            print(f"\n{f}", file=sys.stderr)
        print(f"\n{len(failed)} 个块跑不出文档写的结果。", file=sys.stderr)
        return 1
    print("每一个块都跑出了文档写的结果。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
