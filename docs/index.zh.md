# omgkit

Rust 写的化学信息学工具箱,带 Python 绑定。

```python
import omgkit

m = omgkit.parse_smiles("OC(=O)c1ccccc1N")
m.sanitize()
m.to_canonical_smiles()          # 'Nc1ccccc1C(=O)O'
```

!!! warning "状态:开发中"

    接口在提交之间仍会变。每一层都对外部实现逐条比对过,但表面还没有稳定到
    可以用在生产上。欢迎提 issue。

## 能做什么

| | |
|---|---|
| **读写 SMILES** | 四面体手性、双键顺反、配位键、显式氢,以及规范 SMILES |
| **净化** | 价、隐式氢、环感知、凯库勒化、芳香性、共轭、杂化 |
| **子结构匹配** | SMARTS 解析、VF2++ 匹配、可选的立体敏感、分子与反应的 SMARTS 写出 |
| **应用反应模板** | 产物生成,原子映射号可选 |
| **推演副产物** | 酯化丢掉的那个水,重建成真正的分子;记录本身配不平时,**明说判不了**而不是猜 |
| **画结构式** | 2D 坐标与 SVG/PNG/JPEG 输出,两套绘图规范;画不好的地方**如实报出来** |
| **批处理** | 列式 `MolBatch`,逐分子零拷贝视图 |

## 和别的有什么不一样

**一个分子的性质由它自己决定,不由它被写成什么样决定。** 同一个结构的两种
SMILES 写法,必须净化成同一个东西、匹配同样的查询、按同样的方式反应。这话听着
理所当然,可大量边界情况恰恰藏在这里 —— 这也是每一层都要对着验的那条不变量。

从外面能看出来的三个后果:

**产物分子数由图决定,不由模板决定。** 一个反应模板重写的是**一张图**。出来
几个产物分子,取决于重写后的图有几个连通分量 —— 不取决于作者恰好写了几个产物
模板。这是与常见实现的一处
[**刻意分歧**](dev/correctness.md#a-deliberate-divergence),也正因如此,
应用一个切断环键的模板不会悄悄把原子复制一份。

**被丢弃的原子有交代。** 模板丢掉片段时,omgkit 会如实记下**具体是哪些原子**
被丢了,并能把它们收口成配平的副产物分子。记录本身配不平的时候 —— 比如还原剂
压根不在记录里 —— 它会明说判不了,而不是猜一个。

**每一条声明后面都有判据。** 正确性不是宣称出来的,是逐条对着外部实现比出来的,
而且每条判据都必须先证明自己**不会空过**。见[正确性](dev/correctness.md)。

## 安装

=== "Python"

    ```shell
    pip install maturin
    maturin build --release -m crates/omgkit-py/Cargo.toml --out dist
    pip install dist/omgkit-*.whl
    ```

=== "Rust"

    ```toml
    [dependencies]
    omgkit-core  = "0.0.1"
    omgkit-io    = "0.0.1"
    omgkit-chem  = "0.0.1"
    omgkit-match = "0.0.1"
    ```

按需取用,每一层只依赖它下面的层。完整说明见[安装](getting-started/install.md)。

## 接着看哪里

<div class="grid cards" markdown>

- **[五分钟上手](getting-started/quickstart.md)** —— 解析、净化、匹配、反应
- **[功能与用法](guide/index.md)** —— 一个能力一页,从任务出发
- **[Python API](api/python.md)** —— 每个可调用项,带签名
- **[开发者帮助](dev/index.md)** —— 构建、测试,以及推之前那一整套闸门

</div>

## 许可

[BSD-3-Clause](https://github.com/zbc0315/omgkit/blob/main/LICENSE)。测试语料与
元素表转自其他项目,各自带有自己的条款;每个文件的出处逐条记在
[`THIRD-PARTY-NOTICES.md`](https://github.com/zbc0315/omgkit/blob/main/THIRD-PARTY-NOTICES.md)。
