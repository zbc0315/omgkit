# Building and testing

## The five gates

```shell
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --release
cargo test --workspace
cargo doc --workspace --no-deps --document-private-items

# 四条语料判据,各自带基线,超了就非 0 退出(合计约 10 秒)
cargo run -p omgkit-conformer --release --example audit3d -- harness/corpus/large.smi
cargo run -p omgkit-conformer --release --example angle_audit -- \
    harness/corpus/large.smi harness/params/mmff.angles.tsv
cargo run -p omgkit-conformer --release --example c1_audit -- harness/corpus/large.smi
cargo run -p omgkit-conformer --release --example pucker_audit -- \
    harness/corpus/large.smi harness/params/mmff.pucker.tsv
```

那四条判据**先前一条都不在闸门里**。方案里写着"补上的判据要么进 CI,要么至少进
闸门表",两处都没做 —— 于是"违例不许涨"只是文档里的一句话,没有执行机制。
它们合计约 10 秒、输入全部入库,进 CI 没有额外代价。

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
    "crates/omgkit-conformer",
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

## Conformers: 环褶皱判据

```shell
# 参考表(跑一次就好,约 25 分钟:每个分子嵌 100 个构象取最低能)
python3 harness/measure_pucker.py harness/corpus/large.smi harness/params/mmff

# 判据
cargo run -p omgkit-conformer --release --example pucker_audit -- \
    harness/corpus/large.smi harness/params/mmff.pucker.tsv
```

当前 **459 个判到 / 0 个违例**,超基线非 0 退出。判到的还只有全 sp³ 冠状环与
全 sp² 平面环 —— 混合杂化环一个都摆不成(1459 个环被跳过)。**这条判据真正开始
干活是在那一档落地之后**,现在先把线立住。

判据装了**三道防绕过的闸,每一道都做过变异验证**:逐桶比中位数(把 `crown` 的
±h 改成 0 摆平 → `6元333333i Q` 我们中位 0.000 vs 参考 0.564,红)、
`MIN_JUDGED = 400`(空参考表 → 判到 0,红)、`n_noref == 0`(空表 → 459,红)。
头一道承重:逐环判 band 堵不住"摆平",因为 **14 个桶的 Q p05 恰好是 0.0**,
而 `Q < 1e-9` 时 θ 根本算不出来、直接不判。

**基线是 0**:每一个违例都是真的,涨一个就得说清楚是哪个环。先前是 124,
全是"理想椅式 θ=0.00 vs band 下沿 0.40",我一度把它归成"撞理论极值的假阳性" ——
**独立审稿把这个理由推翻了**(参考里有 5 个真实环 θ 恰好 = 0.0)。
真正站得住的是尺度:0.40° 折成褶皱矢量弦长只有 **0.0039 Å**,而参考协议自己
换个 ETKDG 种子 p05 就在 0.6↔0.8 之间跳。所以 θ 加了一条**对称的、用 Å 定义的
地板**(弦长 < 0.02 Å 不算),Q 不给 —— 单侧化实测会让"把环全摆成椅子"那个
变异蒙混过关(逐环违例 480 → 3、约 60 个环静默)。

桶键是 (环大小, 杂化花样 + **稠合度**)。两处"审稿要求加、量完决定不加"的:
**元素没进键** —— 加杂原子计数会把桶从 69 涨到 274,"样本 ≥20 的桶"覆盖率从
93% 掉到 63%,逐桶比中位数那道主闸会对三分之一的环失明,为修 2 个假阳性不值当;
**五元环相位 φ 量了但不判** —— 参考分布在 φ 上几乎均匀(13 个大桶的 p05 全在
0~1.7、p95 全在 14.1~18.0,而值域就是 [0,18]),band 等于全域。那是化学不是缺陷:
环戊烷赝旋转势垒极低,参考自己没有偏好,判据就没有依据。

参考表有一处曾经写错、已修:`all(GetIsAromatic())` 那道过滤器原先跑在 MMFF
**之后**,而 `MMFFOptimizeMoleculeConfs` 会就地改写芳香标志(实测 776 个分子里
23 个被改过),于是 21 个芳香环被收进表。改成在碰 MMFF **之前**取标志,
表少了 122 个测量(`6元222222f` 151 → 102)。

混合杂化环那一档动手之前要立的三条判据里的第一条,先把**参考**量出来。
键是 (环大小, 杂化花样 + 稠合标记, 量),量的是 Cremer–Pople 的 Q 与 θ
(θ 折到 [0,90]:椅式 0、船/扭船 90)。

**这张表的口径与键长/键角表有一处故意不同,理由是实测出来的。** 那两张表用
"ETKDGv3 单次嵌入(种子 0xf00d)+ MMFF94 局部极小",对键长键角没问题;
**对褶皱是错的** —— MMFF 是局部极小,跨不过椅式↔扭船那道势垒,拿到的是 ETKDG
随机落点的抽签结果:

| | 单次(0xf00d) | 30 次取最低能 |
|---|---:|---:|
| `C1(N2CCCCC2)=C3C(C=CC=C3)=CC=N1` 的哌啶环 | **87.1°**(扭船) | **4.5°**(椅式) |
| 六元 `233333` 桶(一个 sp²+五个 sp³)中位 | **85.6°** | **4.9°** |

按单次口径,连**孤立**的 233333 环都有 57% 落在 θ>60°(船样)—— 那不是化学,
是抽签。用这种分布做 [p05, p95] band,椅式与扭船会一起被放过,判据等于没有。

所以褶皱表**多构象取最低能**(NCONF = 30,种子仍钉 0xf00d)。另外**稠合与否要
进键**:稠环被伙伴钉住、褶皱不自由,与孤立环不是一档。两处都做完之后:

| 桶 | 计数 | 中位 | p05 | p95 |
|---|---:|---:|---:|---:|
| `233333i`(孤立) | 143 | **3.7°** | 0.8 | **10.6** |
| `233333f`(稠合) | 49 | 10.8° | 2.5 | 89.3 |

孤立那一档的 band 才**排得掉扭船**(θ≈90),这正是这条判据存在的理由。

## Conformers: 头号契约(写法无关)的全语料判据

```shell
cargo run -p omgkit-conformer --release --example c1_audit -- harness/corpus/large.smi
```

同一分子用 `write_with_priority` 写成 4 种写法(原样 / 倒着 / 从中间起笔 / 按规范秩),
各自解析回来跑 `generate`,坐标按秩排好**逐位 `to_bits()`** 比,没有容差。
破约就非 0 退出。

**这条是补的,补之前契约在环那一档破着 24.2% 而三处写法不变性判据全绿** ——
它们挑的分子(乙酸乙酯、甲苯、异丁醇)一个都不含两个环系,而破约几乎全在有环那一档
(单键直连 705、靠链相连 1085、其余 39)。联苯 4 种写法坐标最大差 7.922 Å。

**它是抽样,数出来的是下界**:4 种写法加到 4+8 种,破约从 78 涨到 88;
语料里还有 8 个分子四种优先级写出同一个串。全语料判据与手写判据互相补,
不是谁替代谁。

当前 **10/8831 破约**(1 个多片段、9 个单片段有环:螺环、稠环,尚未定位),
超基线非 0 退出。

## Conformers: 键角判据

```shell
cargo run -p omgkit-conformer --release --example angle_audit -- \
    harness/corpus/large.smi harness/params/mmff.angles.tsv
```

参考表由 `harness/measure_params.py` 生成,与 `params.rs` 的键长表同一次跑出来
(ETKDGv3 种子 0xf00d + MMFF94 `maxIters=500`,只取收敛的 8526 个分子)。
键含**三原子共处的最小环尺寸** —— 环张力是真实几何,不是误差:同样是 sp³ 碳,
不在环里 109.5°、六元环内 111.5°、五元环内 **103.3°**。

这条判据是动"混合杂化环"之前必须先立的验收线:原有四条判据没有一条看得见
键角,把 sp³ 的环按平面多边形摆出来会让覆盖面大涨而一条判据都不红。

**超过基线就非 0 退出**(现为 6.2947%),抬基线必须逐桶说明。立线**当轮**就用它
逮住了一件大的:`attach_ring_system` 把环系接反了方向,第二个环折回来叠在第一个上。
按"有原子落在别的环里"量,**3180 → 902 个分子**;键角违例 7.66% → 6.35%,
非键重叠 65.51% → 52.19%。

它抓不到的同样要记住:六元环按平面摆能蒙混过去(现在的 109.47° 本来就在
`C 4 0 6` 的 band 里),键角在**镜像下不变**,也**完全不看扭转/褶皱**。
详见该 example 的模块文档。

## Conformers: the external stereo oracle

```shell
# 三维坐标 → V2000 molblock → RDKit 反读构型
cargo run -p omgkit-conformer --release --example dump_conformer -- \
    harness/corpus/large.smi > /tmp/conf.sdf
python3 harness/check_conformer_stereo.py /tmp/conf.sdf
```

`place` 里的手性是"算带号体积,符号不对就换两个儿子",而 crate 自己的判据用的是
**同一套符号约定** —— 两者共谋,那样的往返是空过的。RDKit 不知道这套约定,
它只看坐标。

当前实测(分子 29 个):**一致 29 / 读反 4 / 判官读不出 2**,判官结论
**不一致**。这一行原先写的是"5 一致 / 0 读反",早就陈了 —— 而它陈掉的正好是
**4 个真的读反**,例如 `CCCc1nc(on1)[C@@H](C)Cn2c(nc(n2)C)C`(原子 8:输入 S,
反读 R)。读不出的那 2 个是三配位的膦,RDKit 自己也从 3D 判不出 CIP,
那是判官的适用范围不是缺陷。手性反读另有任务,不在这一轮。


## Documentation

```shell
cargo doc --workspace --no-deps --open     # Rust API

pip install -r docs/requirements.txt        # the site
mkdocs serve                                # http://127.0.0.1:8000
```
