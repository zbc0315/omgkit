# 第三方材料与归属

本仓库自身的代码按 [BSD-3-Clause](LICENSE) 发布。除此之外,仓库里还**随附**了几份
来自他处的数据文件 —— 语料、元素表、SMARTS 模式集。它们各有各的出处与条款,
这份文件把"哪个文件来自哪里"逐条写清楚。

BSD-3-Clause 的第一条要求:以源码形式再分发时必须保留原始的版权声明。下面的
语料与元素表都取自 RDKit 仓库,所以那份声明照录在此:

```
BSD 3-Clause License

Copyright (c) 2006-2015, Rational Discovery LLC, Greg Landrum, and Julie Penzotti
and others
All rights reserved.
```

RDKit 的完整许可文本见 <https://github.com/rdkit/rdkit/blob/master/license.txt>。
下表的对应关系核对于 **RDKit Release_2025_09_2**。

## 逐文件对照

| 本仓库的文件 | 上游出处 | 上游许可 |
|---|---|---|
| `harness/corpus/large.smi`,命中 5285 条 | RDKit `Code/GraphMol/test_data/canonSmiles.smi` 与 `canonSmiles.long.smi` | BSD-3-Clause |
| 同上,命中 4900 条 | RDKit `Data/NCI/first_5K.smi`,内容取自 NCI/DTP 开放化合物集 | BSD-3-Clause(上游文件) |
| 同上,命中 979 条 | RDKit `Regress/Data/zinc.frags.500.q.smi` 与 `zinc.leads.500.q.smi`,内容取自 ZINC 数据库 | BSD-3-Clause(上游文件);ZINC 见下 |
| `harness/corpus/smarts.txt` | RDKit `Data/FunctionalGroups.txt`、`Data/SmartsLib/RLewis_smarts.txt`、`Data/Pains/wehi_pains.csv`,由 `harness/oracle_smarts.py --build-corpus` 合并去重 | BSD-3-Clause;PAINS 与 RLewis 见下 |
| `crates/omgkit-core/src/element_data.rs` | Blue Obelisk Data Repository 的元素数据,经 RDKit 转录,再由 `harness/gen_elements.py` 生成 | 见下 |
| `harness/baseline/*.jsonl`、`*.tsv` | 不是抄来的文件,是**跑 RDKit 生成的输出**(见 `harness/oracle_*.py`) | 派生数据 |

三个来源**彼此有重叠**(NCI 那批有一部分同时收在 `canonSmiles` 里),所以上面
三个数不该相加。要紧的是并集:去重后正好覆盖 `large.smi` 的全部条目,**余数为 0
—— 没有一条来路不明**。核对命令记在本文件末尾。

(`large.smi` 共 8863 行,生效 8839 条 —— 另外 24 行是**上游就注释掉的**用例,
连同 `#` 一起照搬了过来;去重后 8725 条互不相同的 SMILES。)

`harness/corpus/smoke.smi` 与 `harness/corpus/reactions.txt` 是本项目自己写的,
不涉及第三方。

## 需要一并引用的原始工作

**PAINS 模式**(`wehi_pains.csv` 那一批)出自:

> Baell, J. B.; Holloway, G. A. New Substructure Filters for Removal of Pan Assay
> Interference Compounds (PAINS) from Screening Libraries and for Their Exclusion
> in Bioassays. *J. Med. Chem.* **2010**, *53* (7), 2719–2740.

**RLewis 模式集**在上游文件头里注明由 Richard Lewis 收集并贡献。

**ZINC 子集**出自 ZINC 数据库:

> Irwin, J. J.; Shoichet, B. K. ZINC — A Free Database of Commercially Available
> Compounds for Virtual Screening. *J. Chem. Inf. Model.* **2005**, *45* (1), 177–182.

**NCI 子集**出自美国国家癌症研究所 DTP 的开放化合物集(NCI Open Database),
上游文件是 RDKit 收录的前 5000 条。

## 怎么核对上表

```bash
RD=/路径/到/rdkit          # RDKit 源码树
cut -f1 harness/corpus/large.smi | sort -u > /tmp/og.smi

cat $RD/Code/GraphMol/test_data/canonSmiles.smi \
    $RD/Code/GraphMol/test_data/canonSmiles.long.smi | awk '{print $1}' | sort -u > /tmp/a
cat $RD/Data/NCI/first_5K.smi | awk '{print $1}' | sort -u > /tmp/b
cat $RD/Regress/Data/zinc.{frags,leads}.500.q.smi | awk '{print $1}' | sort -u > /tmp/c

for f in /tmp/a /tmp/b /tmp/c; do comm -12 $f /tmp/og.smi | wc -l; done   # 5285 / 4900 / 979
comm -13 <(cat /tmp/a /tmp/b /tmp/c | sort -u) /tmp/og.smi | wc -l        # 0 —— 没有来路不明的
```

最后那一行是这张表的判据:**它必须是 0**。往语料里添东西时先跑一遍,不为 0 就
说明添进来的东西还没写进上表。前三个数只是各来源的命中量,来源之间有重叠,
**不要拿它们相加去凑总数**。

## 一处尚未核实到底的条款

元素数据的链路是 **BODR → RDKit → 本仓库**。RDKit 那一段是 BSD-3-Clause,已照录
在上;**BODR 自身的条款没有在本仓库里独立核实过**。这些数据是元素的物理常数
(符号、质量、价态、共价半径),事实性数据本身通常不构成受保护的表达,但真要把
仓库公开发布之前,这一条应当再确认一次,不要凭上面这句话就当已经清楚了。

## 语料为什么要随仓库带着

差分测试的价值全在"用**真实**语料、不手挑"。手挑的分子与 SMARTS 会不自觉地只挑
实现已经处理得了的写法,而真实语料里 `=!@`、`$(...)`、稠环上的方向键这类正是最
容易漏的。所以语料必须是外来的 —— 也因此才有了这份归属文件。
