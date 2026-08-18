# Building and testing

## The five gates

一条命令跑全部:

```shell
bash harness/gates.sh
```

它就是下面这五条,顺序执行、任何一条非 0 就停:

```shell
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --release
cargo test --workspace                     # debug,让 debug_assert 真的跑到
cargo doc --workspace --no-deps --document-private-items
```

**语料判据要随算法一起回来。** 判据不进 CI,"违例不许涨"就只是文档里的一句话 ——
上一轮正是这么漏的:判据只进了本地脚本,推送时一条都不执行,而 CI 一直是绿的。

`harness/gates.sh` 里 `set -eo pipefail`,而且**判据后面不许接管道** ——
`set -e` 遇到管道只看最后一个命令的退出码,`tail` 永远成功,于是判据非 0 退出、
脚本照样打印"全部通过"。那是自己造的绿。要看少几行,把管道加在**外面**:
`bash harness/gates.sh 2>&1 | tail -40`。

These are exactly what
[CI](https://github.com/zbc0315/omgkit/actions/workflows/ci.yml) runs. Run them
before opening a pull request.

`cargo test --release` rather than debug: the differential tests run over a
large corpus and a debug build is too slow to be useful.

`cargo test --workspace` on top of it, in debug, for one reason: **`debug_assert!`
is compiled out of a release build.** The workspace has 37 of them — "bond
geometry was written but cis/trans was never perceived", "canonical ranks are not
injective", "a flip reached into another fragment" — and with only the release run
in CI, not one of them had ever been executed by a gate. A guard that never runs
is a guard that is not there. The debug run takes about 6 seconds for the whole
workspace, so it buys that back cheaply.

## Why `omgkit-py` is not in the default workspace members

```toml
default-members = [
    "crates/omgkit-core",
    "crates/omgkit-io",
    "crates/omgkit-chem",
    "crates/omgkit-match",
    "crates/omgkit-depict",
]
```

`omgkit-py` is a `cdylib` that leaves Python symbols unresolved at link time.
Building a test executable for it necessarily fails, so putting it in the
default set would turn every project-wide verification command red. It stays in
`members` so it shares the lockfile and target directory and `cargo -p` can
reach it.

## Do not add `panic = "abort"`

The Python extension relies on unwinding to catch Rust panics and turn them
into Python exceptions. With `abort`, a panic sends `SIGABRT` to the
interpreter — exit code 134, no exception, no traceback, nothing `try/except`
can catch, and the user loses their whole process along with any unsaved work.

Cargo does not allow overriding `panic` per package, so this is a workspace-wide
trade. The cost is the unwind tables and a small constant overhead, which does
not show up on the pipeline benchmark.

## The Python extension

```shell
pip install maturin
maturin build --release -m crates/omgkit-py/Cargo.toml --out dist
pip install --force-reinstall dist/omgkit-*.whl
python harness/test_python.py
```

## Differential tests

The tests come in two tiers.

**Smoke tier** — oracles are committed (about 680 KB), runs by default, green on
a fresh clone.

**Large-corpus tier** — marked `#[ignore]`, needs oracles you generate against
an external reference implementation:

```shell
cargo test --release -- --ignored
```

Generating the oracles, the column conventions for each layer, and what each
judge is guarding are documented in
[`harness/README.md`](https://github.com/zbc0315/omgkit/blob/main/harness/README.md)
— also readable here as
[the full text](correctness-full.md) (Chinese).

!!! warning "The gitignore rule for oracles is about tests, not filenames"

    Baselines that a non-`#[ignore]` test reaches must be committed, or a fresh
    clone fails immediately. Guessing by filename has gone wrong twice: once a
    rule written as `smoke.*.jsonl` excluded `smoke.matches.tsv` by suffix, once
    `smarts.jsonl` was missed because it is not named `smoke` but a non-ignored
    test hard-codes it.

    The way to verify a change to that rule is not to read it — it is to clone
    the repository somewhere else and run `cargo test`.

## Drawing: the gallery and the corpus audit

`omgkit-depict` has a third tier of its own — the properties a picture must
satisfy, run over the whole corpus rather than a hand-picked list:

```shell
# six decidable properties over 8831 molecules × 2 styles
cargo run -p omgkit-depict --release --example audit -- harness/corpus/large.smi

# eyeball it: 17 molecules × 2 styles × svg/png/jpg
cargo run -p omgkit-depict --release --features raster --example draw -- out/

# side by side with RDKit at the same bond length, then bind the lot
# into one PDF (byte-identical across runs; refuses to run if the manifest
# and the images on disk disagree in either direction)
python3 harness/compare_rdkit.py out/
python3 harness/make_gallery.py out/
```

!!! warning "Run the audit twice"

    One class of defect only shows up across processes: `HashMap` iteration
    order is seeded per run, so anything that sums positions in that order can
    silently give a different picture each time. Two runs that disagree is the
    only way to see it — no unit test can, because the seed is fixed within a
    process.

The `raster` feature is optional on purpose: without it the crate has **no
external dependencies** and emits SVG only.

## 初始构型:重新设计中

`omgkit-conf`(v2)已连同 v1 一起撤销。两版都栽在同一件事上:
**按分子的类别切分支** —— 无环走构造法、有环另说、超配位拒绝、累积双键再加一条。
每来一类分子就多一个分支,覆盖率停在 **14.3%**(1259 / 8831),
而 RDKit ETKDG 是 99.48%。

要的是**一个**对所有有机小分子都成立的算法,不是一堆特例的并集。
方案重写中(`dev-notes/`,不入库)。这一节等算法定下来再补。

代码在 git 历史里(`4eefccd` 及之前),需要参考时 `git show` 取。
`harness/params/` 里的实测参数表与 `harness/baseline_rdkit_etkdg.py` 留着 ——
那些是从语料量出来的**数据**与基线,与算法怎么写无关。

## Documentation

```shell
cargo doc --workspace --no-deps --open     # Rust API

pip install -r docs/requirements.txt        # the site
mkdocs serve                                # http://127.0.0.1:8000
```
