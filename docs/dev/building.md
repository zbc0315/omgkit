# Building and testing

## The six gates

一条命令跑全部:

```shell
bash harness/gates.sh
```

它就是下面这六条,顺序执行、任何一条非 0 就停:

```shell
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --release
cargo test --workspace                     # debug,让 debug_assert 真的跑到
cargo run -p omgkit-conf --release --example conf_audit -- harness/corpus/large.smi
cargo doc --workspace --no-deps --document-private-items
```

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

## 初始构型:`omgkit-conf`

```shell
cargo run -p omgkit-conf --release --example conf_audit -- harness/corpus/large.smi
```

**这是第 5 道闸。** 它跑全语料的构造法 + 消撞,量的是"能不能当力场优化的起点",
**不是"像不像 MMFF 优化后的样子"**。

一期的范围是**无环分子**(1278 个)。有环的如实计数、不判 —— 是范围不是失败。

当前实测:

| | |
|---|---|
| 范围内原子 | 38609 / 38609 全摆好,退化 0、说不出原因 0 |
| 键长最大相对误差 | 3.044e-15(构造法这一条按定义是机器精度) |
| 键角超出容差 | 0.000° |
| 参数命中 | 键长 99.7% / 键角 99.8% 走实测表 |
| 非键地板(消撞前 → 后) | 0.011 → **1.489 Å** |
| 非键中位(消撞前 → 后) | 1.739 → **2.052 Å**(RDKit 是 2.07) |
| 消撞代价 | 74 μs/分子(RDKit ETKDG 中位 3.6 ms) |

### 判据自己会不会撒谎

这一条被独立审核逮到过四道**恒真**的闸,所以现在的每一条都做过变异验证。
记住这几个教训:

- **别拿"总数减去别的"当计数器。** 原来 `skipped_hypervalent` 是残差,
  于是守恒式是代数恒等式:往 BFS 里插一条静默 `continue` 丢掉 191 个原子,
  守恒闸、覆盖率闸**一个都不响**。现在每个桶在自己的分支里真的累加,
  残差单列 `unaccounted` 钉在 0。
- **单向放松的容差必须配一道上限闸。** 键角判据原来把兄弟角的解析容差
  减在**所有**原子对上,而容差自己没有上限(`sibling_skew(180°) = 180°`)——
  某些中心的角判据等于整个关掉。现在父–子容差 0、兄弟才有容差,而且封顶 35°
  并单独计数。
- **判据不该与被判的代码共用实现。** 判据原来用 `vsepr::arrangement()` 算期望值,
  而构造法用的也是它 —— 把 `Sp2 => Planar` 改成 `=> Tetrahedral`(每个羰基、
  酰胺、烯、硝基都变成锥形),**全套闸照绿**。现在片段分量判据自己算,
  另外补了三条对扭转敏感的棘轮(消撞前的中位/撞的对数、兄弟角推歪的中心数)。
- **量构造法的质量要在消撞之前。** 消撞会把大部分伤害盖掉:NeRF 标架写反之后
  消撞后的地板还有 1.3 Å,而消撞**前**的中位从 1.739 掉到 0.543。

### 消撞为什么只转扭转角

绕一根单键转动**整个子树**是刚体运动 —— 键长、键角**逐位不变**。
任何平移单个原子的消撞都会把已经精确到 3e-15 的那两项弄脏。
代价是它能做的事有上界:固定键角撑不开的拥挤,转扭转角也撑不开
(候选角从 24 加到 72,仍有 4 个分子的非键距离低于 1.6 Å)。
放宽键角是后面那一期的事。

接受准则**不带 epsilon**,所以"只会变好、不会变坏"是**推出来**的:
跨切口那部分的最小间距比不减、平手时罚和严格变小 ⟹ 全局按字典序单调不减。

### RDKit 基线

```shell
python3 harness/baseline_rdkit_etkdg.py harness/corpus/large.smi
```

crate 文档里引的那组数(8839 个分子 / 86.5 s / 中位 3.6 ms / 最慢单分子 8.11 s /
46 个失败)由它跑出来,单进程、种子 `0xf00d`、`AddHs`,与 `measure_params.py` 一致。

## Documentation

```shell
cargo doc --workspace --no-deps --open     # Rust API

pip install -r docs/requirements.txt        # the site
mkdocs serve                                # http://127.0.0.1:8000
```
