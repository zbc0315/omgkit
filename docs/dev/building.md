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
| 界矩阵 | 真实构象要落在界内 + 界宽不许比 RDKit 松 + `U` 要摆得进三维 | 越界 0.863%;宽度比 1.004;1-2/1-3/1-4 三档与 RDKit 逐位相同 |
| **界可行率(头号指标)** | 全语料 8831 个分子,界自相矛盾的比例 | **0.01%**(8831 个里 1 个),是 RDKit ETKDG 那 0.52% 的 1/52 |
| 手性中心 | 真值 = 真实构象上量出来的有符号体积 | 247 个中心,符号错 0、漏抽 0 |
| **端到端(产物)** | 精修前后各量一遍:越界、自穿、手性、耗时 | 1-2/1-3/1-4/长程 越界 **0.0 / 0.1 / 0.0 / 0.0%**;键交叉 1731 → **0**;手性 84.2% → **100%**;**0.97 ms/分子** |
| 自穿 | 先拿真实构象校准检测器(必须报 0),再量自己 | 真实构象 0 误报;我们嵌出来的环穿刺 1/400(0.2%),键交叉 1821 —— 后者正是精修要收拾的 |
| **通用性难例语料** | 68 个分子,照着算法的假设挑:笼状/张力、超配位、累积双键、超大环、少见元素、金属、自由基、两性离子 | 建界即空 0、界不可行 **0**(同一批分子 RDKit ETKDGv3 失败 2 个:SF₆ 与六氨合钴) |
| 特征分解 + 嵌入 | numpy `eigvalsh`(LAPACK)+ 真实构象精确回嵌 | 特征值偏差 5.96e-15;回嵌偏差 1.76e-11 Å |

整条流水线已经通了:**界矩阵 → 三角光滑化 → 取 `U` 当参考距离表 → 度量矩阵嵌入
→ 全局手性定向(离散一次)→ L-BFGS 精修**,全程无随机数。

还没做:确定性的重试阶梯、对称性破除、1-5 链式约束。**已知的一笔债**:饱和环的
1-4 扭转退回了全程(界宽比 7.68×),换来的是约束自洽 —— 正确修法是把
`ring_internal_torsion` 的分桶扩到"是否全 sp³"并用与键角自洽的值,要重跑参数表实测。

**整体手性不能指望三维精修去修**:翻转手性是反射(`det = −1`),不在 `SO(3)` 的
连通分支里,连续下降要走到镜像必须把分子压平,下降法不会付这个势垒。所以嵌入之后
**离散地**定一次全局定向 —— 这一步已经落地。实测(247 个中心):

| | 手性号正确的中心 |
|---|---|
| 嵌完直接看 | 53.0%(基本是掷硬币,与"定向任意"吻合) |
| 做一次全局反射后 | **86.2%** |

剩下的 13.8% 是个别中心相对多数错,全局反射按定义救不了,但**三维精修救得了** ——
翻一个中心只要它自己的体积过零(局部、有限势垒),而全局反射要求所有中心同时压平。
**所以四维先不做**,等精修落地再量剩下多少。(RDKit 一有手性中心就上四维,
为的正是让手性能连续翻:四维里 `(x₃, x₄)` 平面转 π 就把 `x₃` 送到 `−x₃`,
而四维两两距离精确不变。)

**"通用"的判定标准不是"在语料上都过",是"来了一类没见过的分子,要改的是不是只有
约束表"。** 代码里出现 `is_metal` / `is_macrocycle` 这类分子类别谓词就算违规 ——
前两版 `omgkit-conf` 正是死在这里。难例语料(`harness/corpus/hard.smi`)是这条标准
的探针:每一类精确攻击一条假设,红了就说明那条假设不成立,而修法只许是补表。

`harness/params/` 里的实测参数表与 `harness/baseline_rdkit_etkdg.py` 是从语料量出来的
**数据**与基线,与算法怎么写无关,一直留着。

## Documentation

```shell
cargo doc --workspace --no-deps --open     # Rust API

pip install -r docs/requirements.txt        # the site
mkdocs serve                                # http://127.0.0.1:8000
```
