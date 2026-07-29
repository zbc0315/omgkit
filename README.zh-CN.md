# omgkit

[![CI](https://github.com/zbc0315/omgkit/actions/workflows/ci.yml/badge.svg)](https://github.com/zbc0315/omgkit/actions/workflows/ci.yml)

[English](README.md)

## 这是什么

omgkit 是一个用 Rust 写的化学信息学工具箱,附 Python 绑定。它围绕**列式分子表示**
搭建 —— 同一属性的所有原子躺在一个连续数组里,于是大批量分子的遍历、并行切分,
以及零拷贝交给 numpy / Arrow,都变成顺理成章的事。

  * [BSD-3-Clause 许可](LICENSE)
  * 核心数据结构与算法用 Rust 写,不含 `unsafe`
  * Python 3.9 及以上,由 [PyO3](https://pyo3.rs) 与
    [maturin](https://github.com/PyO3/maturin) 构建 —— abi3 编译,一个 wheel
    覆盖所有支持的版本,不带系统依赖
  * SMILES 解析与写出,含四面体手性、双键顺反、配位键、显式氢
  * 规范 SMILES
  * 净化管线:价键、隐式氢、环感知、kekulize、芳香性、共轭、杂化
  * SMARTS 解析、子结构匹配(VF2++,可选按立体判定),以及分子与反应的 SMARTS 写出
  * 反应模板与产物生成,可选同步给出原子映射号
  * 把模板丢弃的片段(酯化掉的那个水)收口成账平的副产物分子;记录本身不平时
    如实报"答不了",不编分子
  * 列式批(`MolBatch`)与零拷贝的单分子视图

**状态:开发中。** 接口仍在变。每一层都对着一个外部参照逐条比对过
(见[文档](#文档)),但接口还没稳到可以用于生产。遇到问题请提 issue。

## 安装

### Python

需要 [maturin](https://github.com/PyO3/maturin):

```shell-session
$ maturin build --release -m crates/omgkit-py/Cargo.toml --out dist
$ pip install dist/omgkit-*.whl
```

### Rust

```toml
[dependencies]
omgkit-core  = { git = "https://github.com/zbc0315/omgkit" }   # 数据结构
omgkit-io    = { git = "https://github.com/zbc0315/omgkit" }   # SMILES / SMARTS
omgkit-chem  = { git = "https://github.com/zbc0315/omgkit" }   # 净化
omgkit-match = { git = "https://github.com/zbc0315/omgkit" }   # 匹配与反应
```

按需取用,每一层只依赖它下面的层。

## 快速开始

```python
import omgkit

m = omgkit.parse_smiles("OC(=O)c1ccccc1N")
m.sanitize()
m.to_canonical_smiles()

q = omgkit.parse_smarts("[C](=[O])[OH]")
q.match(m)                      # 按查询原子顺序给出分子原子下标

rxn = omgkit.parse_reaction("[C:1][OH:2]>>[C:1][Cl:2]")
for outcome in rxn.run([m], atom_mapping=True):
    outcome.products, outcome.reactants
```

Rust 侧对应的入口在 `omgkit_io::smiles`、`omgkit_chem::sanitize` 与
`omgkit_match`,可直接运行的例子见各 crate 的文档。

## 文档

  * [`docs/design.md`](docs/design.md) —— 每一层做什么、为什么这么做,
    以及每个设计取舍是怎么验证的
  * [`harness/README.md`](harness/README.md) —— 差分测试基础设施:基准怎么生成、
    判据怎么写才不会空过
  * `cargo doc --workspace --no-deps --open` —— API 文档

## 参与

欢迎提 issue 与 PR。要过四道闸门,与 CI 跑的是同一套:

```shell-session
$ cargo fmt --all --check
$ cargo clippy --workspace --all-targets -- -D warnings
$ cargo test
$ cargo doc --workspace --no-deps
```

新克隆下来 `cargo test` 即为绿:冒烟基准随仓库入库。大语料那一档标了 `#[ignore]`,
它的基准要自己生成,方式见 [`harness/README.md`](harness/README.md)。

## 许可

代码按 [BSD-3-Clause](LICENSE) 发布。

测试语料与元素表是从别的项目再分发过来的,各有各的条款,逐文件的出处写在
[`THIRD-PARTY-NOTICES.md`](THIRD-PARTY-NOTICES.md)。
