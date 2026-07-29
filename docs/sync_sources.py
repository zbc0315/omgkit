"""把 docs/ 之外的长文档复制进站点,构建前跑。

`docs/design.md` 与 `harness/README.md` 是仓库里的工程记录。GitHub 上直接浏览
它们的链接(README 里就有)不能断,所以**源文件留在原处不动**,这里只复制一份
到站点目录。复制出来的文件已 gitignore。

本来用的是 mkdocs-gen-files,但 i18n 插件处理不了它生成的虚拟文件
(`Unhandled file case`),页面会静默地从站点里消失 —— 构建不报错,导航里
却什么都没有。改成落成真实文件就没这个问题。

# 判据

源文件不在就**直接失败**,不生成空页。默认回落成空页的话,站点上会出现一个
看着正常、其实什么都没有的页面 —— 那比构建红了难发现得多。
"""

import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).parent.parent
BLOB = "https://github.com/zbc0315/omgkit/blob/main"

# `](目标)`,目标不是 URL、不是锚点、不是邮件
_LINK = re.compile(r"\]\((?!https?://|#|mailto:)([^)\s]+)\)")


def absolutise(body, src):
    """把指向**仓库文件**的相对链接改写成 GitHub 绝对地址。

    这两份全文原本躺在 `docs/` 和 `harness/`,里面的 `../README.md` 之类是相对
    仓库路径写的。复制到 `docs/dev/` 之后那些路径全都解析不到了 —— 而 MkDocs
    只会报 warning,非 strict 构建下就是一批**静默的死链**。

    它们指向的是仓库文件而不是站点页面,所以正确的去处是 GitHub,不是站内。
    """
    base = (ROOT / src).parent

    def repl(m):
        target, _, anchor = m.group(1).partition("#")
        if not target:
            return m.group(0)  # 纯锚点,站内解析
        resolved = (base / target).resolve()
        try:
            rel = resolved.relative_to(ROOT.resolve())
        except ValueError:
            return m.group(0)  # 指到仓库外,原样留着让 MkDocs 去报
        if not resolved.exists():
            return m.group(0)  # 指不到东西,原样留着 —— 别把死链藏起来
        return f"]({BLOB}/{rel}{'#' + anchor if anchor else ''})"

    return _LINK.sub(repl, body)

# (源文件, 站点路径, 页面标题, 这一页在讲什么)
SOURCES = [
    (
        "docs/design.md",
        "docs/dev/design-full.md",
        "设计(全文)",
        "逐层的设计取舍与验证方式。这是仓库里的工程记录原文,中文。",
    ),
    (
        "harness/README.md",
        "docs/dev/correctness-full.md",
        "差分测试基础设施(全文)",
        "每一条判据守的是什么、怎么证明它不会空过、覆盖多少。仓库里的工程记录原文,中文。",
    ),
]


def main():
    for src, dest, title, blurb in SOURCES:
        path = ROOT / src
        if not path.is_file():
            sys.exit(f"✗ 文档源文件不在:{src} —— 不生成空页,直接失败")
        lines = absolutise(path.read_text(encoding="utf-8"), src).splitlines()
        if lines and lines[0].startswith("# "):
            lines = lines[1:]          # 原文自带一级标题,去掉以免重复
        out = ROOT / dest
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(
            f"# {title}\n\n"
            f'!!! note ""\n'
            f"    {blurb}\n\n"
            f"    源文件:[`{src}`](https://github.com/zbc0315/omgkit/blob/main/{src})\n\n"
            + "\n".join(lines).lstrip("\n")
            + "\n",
            encoding="utf-8",
        )
        print(f"✓ {src} → {dest}")


if __name__ == "__main__":
    main()
