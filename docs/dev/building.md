# Building and testing

## The gates

一条命令跑全部:

```shell
bash harness/gates.sh
```

顺序执行、任何一条非 0 就停:

```shell
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --release
cargo test --workspace                     # debug,让 debug_assert 真的跑到
cargo doc --workspace --no-deps --document-private-items

# omgkit-conf 的三个外部判官,跑冒烟档(基准随仓库入库)
SMOKE=harness/baseline/smoke.bounds.jsonl
cargo run -p omgkit-conf --release --example smooth_oracle -- $SMOKE
cargo run -p omgkit-conf --release --example bounds_oracle -- $SMOKE
cargo run -p omgkit-conf --release --example eigen_oracle  -- $SMOKE \
    harness/baseline/smoke.gram_eigs.jsonl
```

(这一节的标题原先写着 "The five gates",而条数早已从四变五、现在是八 ——
数字写进标题就会掉队。现在不写数字:要知道有几道,数上面的命令。)

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
    "crates/omgkit-conf",
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

**Smoke tier** — oracles are committed (about 1.3 MB), runs by default, green on
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

## 初始构型:确定性距离几何(建设中)

`omgkit-conf` 的 v1 与 v2 都撤了。两版栽在同一件事上:**按分子的类别切分支** ——
无环走构造法、有环另说、超配位拒绝、累积双键再加一条。每来一类分子就多一个分支,
覆盖率停在 **14.3%**(1259 / 8831),而 RDKit ETKDG 是 99.48%。
旧代码在 git 历史里(`4eefccd` 及之前),需要参考时 `git show` 取。

现在这一版走的是距离几何,与 RDKit 同一条主干,**只在一处分岔**:
RDKit 在界矩阵里**逐对独立随机取**一组距离,取出来的表常常任何空间都摆不出来,
它的应对是作废整次尝试重掷;这里直接拿三角光滑化之后的上限矩阵 `U` 当参考距离表 ——
`U` 按构造满足三角不等式,而且**全程没有随机数**。

已落地的分块,每块都配了外部判官(判官不进 CI 就不是闸,所以三条都在上面的闸门里):

| 分块 | 判官 | 现状 |
|---|---|---|
| 三角光滑化 | RDKit 的 `GetMoleculeBoundsMatrix` 带/不带 smoothing | 逐位相同,最大偏差 5.3e-15 |
| 界矩阵 | 真实构象要落在界内 + 界宽不许比 RDKit 松 + `U` 要摆得进三维 | 越界 0.607%;宽度比 1.020;1-2/1-3/1-4 三档与 RDKit 逐位相同 |
| 特征分解 + 嵌入 | numpy `eigvalsh`(LAPACK)+ 真实构象精确回嵌 | 特征值偏差 5.96e-15;回嵌偏差 1.76e-11 Å |

还没做:误差函数与优化器、立体化学的有符号体积项、确定性的对称性破除与修复阶梯。

`harness/params/` 里的实测参数表与 `harness/baseline_rdkit_etkdg.py` 是从语料量出来的
**数据**与基线,与算法怎么写无关,一直留着。

## Documentation

```shell
cargo doc --workspace --no-deps --open     # Rust API

pip install -r docs/requirements.txt        # the site
mkdocs serve                                # http://127.0.0.1:8000
```
