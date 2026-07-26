# omgkit

面向批处理与 GPU 的化学信息学工具箱。Rust 实现,附 Python 绑定。

> **状态:开发中。** 接口仍会变动,不建议用于生产。各项能力都有差分测试守着
> (见[验证](#验证)),但覆盖面还在扩,遇到问题请提 issue。

## 它解决什么

处理百万级分子时,瓶颈常常不在算法而在**内存布局**。把分子存成对象图,遍历一个
40 原子的分子就是 40 次随机访存;换成列式布局,同一属性的所有原子躺在一个连续
数组里,顺序访存、并行切分、零拷贝导出 numpy/Arrow、以及送进 GPU 都变成顺理成章
的事。omgkit 就是围绕这个决策搭起来的。

## 安装

### Rust

```toml
[dependencies]
omgkit-core = { git = "https://github.com/zbc0315/omgkit" }
omgkit-io   = { git = "https://github.com/zbc0315/omgkit" }
omgkit-chem = { git = "https://github.com/zbc0315/omgkit" }
omgkit-match = { git = "https://github.com/zbc0315/omgkit" }
```

按需取用:`core` 是数据结构,`io` 管 SMILES/SMARTS 的读写,`chem` 是净化管线,
`match` 是子结构匹配与反应。

### Python

需要 [maturin](https://github.com/PyO3/maturin)。abi3 编译,**一个 wheel 覆盖
Python 3.9 及以上**,不带系统依赖。

```bash
maturin build --release -m crates/omgkit-py/Cargo.toml --out dist
pip install dist/omgkit-*.whl
```

## 快速开始

### Python

```python
import omgkit

m = omgkit.parse_smiles("OC(=O)c1ccccc1N")
m.sanitize()
m.to_canonical_smiles()                  # 同一个分子,写法再怎么变都是同一串

q = omgkit.parse_smarts("[C](=[O])[OH]")
q.match(m)                               # [[1, 2, 0]] —— 按查询原子顺序给出分子原子下标

rxn = omgkit.parse_reaction("[C:1][OH:2]>>[C:1][Cl:2]")
for o in rxn.run([m], atom_mapping=True):
    o.products, o.reactants              # 两侧都带原子映射号
```

### Rust

```rust
use omgkit_io::{canon, smiles};

let a = smiles::parse("OC(=O)c1ccccc1N")?;
let b = smiles::parse("Nc1ccccc1C(O)=O")?;   // 同一个分子,写法不同
assert_eq!(canon::canonical_smiles(&a).smiles,
           canon::canonical_smiles(&b).smiles);
```

解析错误带精确到列的位置:

```rust
let err = smiles::parse("C1CC").unwrap_err();
println!("{}", err.render());
// C1CC
//  ^ 环闭合标号 1 未配对
```

批处理走 `MolBatch` —— 列是连续内存,单分子取出来是零拷贝视图:

```rust
use omgkit_core::{BondOrder, MolBatchBuilder, MolBuilder};

let mut ethanol = MolBuilder::new();
let c0 = ethanol.add_atom(6);
let c1 = ethanol.add_atom(6);
let o  = ethanol.add_atom(8);
ethanol.add_bond(c0, c1, BondOrder::Single)?;
ethanol.add_bond(c1, o,  BondOrder::Single)?;

let mut bb = MolBatchBuilder::new();
bb.push(&ethanol)?;
let batch = bb.finish();

assert_eq!(batch.atomic_nums(), &[6, 6, 8]);
assert_eq!(batch.mol(0).unwrap().degree(1), 2);
```

## 能做什么

| | Rust | Python |
|---|---|---|
| SMILES 解析(含立体、配位键、显式氢) | ✓ | ✓ |
| 净化:价键、隐式氢、环感知、kekulize、芳香性、共轭、杂化 | ✓ | ✓ |
| SMILES 写出(含四面体手性、双键顺反) | ✓ | ✓ |
| 规范 SMILES | ✓ | ✓ |
| SMARTS 解析 + 子结构匹配(VF2++,可选按手性判定) | ✓ | ✓ |
| SMARTS 写出(分子与反应) | ✓ | — |
| 反应模板与产物生成(可选原子映射号) | ✓ | ✓ |
| 显式氢并入氢计数(`removeHs`) | ✓ | ✓ |
| 列式批 / 零拷贝视图 | ✓ | — |

Python 侧刻意只暴露少数入口(`parse_smiles` / `parse_smarts` / `parse_reaction`
以及 `Mol`、`Query`、`Reaction` 上的方法),还在按需要扩。绑定层只做翻译,不放
任何判断分子的逻辑 —— 那样的逻辑只有 Python 用户碰得到,Rust 侧的整套差分测试
一概盖不到。

**尚未支持**:阻转异构(净化第 10 步)、轴手性、配位几何(`@SP`/`@TB`/`@OH`)
的写出。前两项与配位几何在 8839 条真实语料上的触发次数都是 **0**,所以在补到能
触发的语料之前不做 —— 没有用例守着的实现比没有实现更危险。

## 验证

每一项能力都对着一个**外部参照**逐条比对,而不是只测自己想得到的用例:

| 面 | 规模 | 结果 |
|---|---|---|
| SMILES 解析(逐字段) | 8839 条 | 零分歧 |
| 净化 12 步(逐步单独验) | 8839 条 × 每步 | 零分歧 |
| SMILES 写出(往返恒等) | 8839 条 | 全部往返成功 |
| 规范 SMILES(随机重排) | 8839 条 × 5 次 | 全部恒等 |
| SMARTS 解析 | 776 条真实模式 | 双向零分歧(756 条都成功、20 条都拒绝)|
| 子结构匹配 | 11591 组分子×模式 | 零分歧 |
| 反应产物 | 741 组 | 717 一致,24 组差在一处刻意的设计选择 |

判据怎么写才不会**空过**(测试通过只是因为压根没走到那条路)、覆盖面怎么算,
写在 [`harness/README.md`](harness/README.md)。各层的设计取舍与上表每个数字的
来龙去脉,写在 [`docs/design.md`](docs/design.md)。

整条管线还必须**线性于分子规模**。这类问题差分测试抓不到 —— 结果全对,只是慢,
而且在小分子上完全看不出来,所以另有一套盯增长曲线的测试。

## 开发

```bash
cargo test                       # 单元测试 + 差分测试(冒烟语料)
cargo test --release -- --ignored  # 大语料那一档
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --all --check
cargo doc --workspace --no-deps  # rustdoc 断链等已设为 deny
```

以上四道就是 CI 跑的闸门。冒烟基准随仓库入库,新克隆下来 `cargo test` 即为绿;
大语料那一档的基准要自己生成,方式见 `harness/README.md`。

## 许可

BSD-3-Clause。元素数据源自 Blue Obelisk Data Repository,经 RDKit
(BSD-3-Clause)转录。
