#!/bin/bash
#
# 推之前在本地跑一遍全部闸门。**权威仍然是 `.github/workflows/ci.yml`** ——
# 这个脚本只是让本地跑一次不用手抄十几条命令,两边不一致时以 CI 为准。
#
# # 为什么这个文件必须存在
#
# 先前每一轮都是现敲一个临时脚本,而临时脚本里写的是:
#
#     set -e
#     cargo run ... --example some_audit ... | tail -3
#
# **`set -e` 遇到管道只看最后一个命令的退出码**,`tail` 永远成功。
# 于是判据非 0 退出、脚本照样跑到底,最后打印"全部十道闸通过"。
# 实测:拿一个会让三条闸变红的变异去跑,脚本退出码仍是 0。
# 那是**自己造的绿** —— 与 `docs/dev/` 里记过的 zsh `PIPESTATUS` 恒空是同一个病。
#
# 所以这里 `set -eo pipefail`,而且**不许在判据后面接管道**。要看少几行,
# 用 `bash harness/gates.sh 2>&1 | tail -40`,管道加在外面。
set -eo pipefail

cd "$(dirname "$0")/.."

# **步数只写一处。** 先前每一行都硬写着 `== 7/14 …`,加一道闸要改十几处,
# 漏一处就是个不会报错的假数。CI 的头注释里记过同一个坑(那里原先写着
# "四道闸门",而步骤早已加到八步)。这里由 `step` 计数,末尾自查:
# 改了步骤忘了改 `TOTAL`,脚本最后一行会红。
TOTAL=18
N=0
step() {
    N=$((N + 1))
    echo "== $N/$TOTAL $1"
}

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

# **最后四条判据要 RDKit,所以先在这里查一次。** 放在末尾的话,要等前面十四步
# (十来分钟 cargo)跑完才知道环境缺东西。没有就直接失败,**不跳过** ——
# 静默跳过的判据是最坏的一种,它让人以为跑过了。
PY=.venv/bin/python
if [ ! -x "$PY" ]; then
    echo "缺 $PY —— 最后四条判据要 RDKit。" >&2
    echo "  建法:python3 -m venv .venv &&" >&2
    echo "        .venv/bin/pip install --only-binary=:all: -r harness/requirements.lock" >&2
    exit 1
fi

step "fmt"
cargo fmt --all --check
step "clippy(警告即失败)"
cargo clippy -q --workspace --all-targets -- -D warnings
step "测试(release)"
cargo test -q --release
step "测试(debug —— 让 debug_assert 真的跑到)"
cargo test -q --workspace
step "文档"
cargo doc -q --workspace --no-deps --document-private-items

# ---- 拿预先烘好的基准比的判官 ----
#
# **判据不进这里就等于没有。** 先前 omgkit-conf 的三条判官全靠手动跑,
# 而 CI 与本脚本都只有上面那五步 —— 于是"界宽比 1.020"这类数是靠自觉维持的,
# 谁都可能在不知情的时候把它推回去而全程绿灯。
#
# 跑的是**冒烟档**(`smoke.bounds.jsonl`,27 个分子,随仓库入库),
# 因为全量基准 7.7 M 不入库。全量档在本地手动跑:
#
#   cargo run -p omgkit-conf --release --example bounds_oracle -- harness/baseline/rdkit_bounds.jsonl
#
# 冒烟档跑不出全量档的统计精度,但三条判官里有两条是**逐分子的硬判据**
# (光滑化要逐位相同、特征值要对上 LAPACK、真实构象要精确回嵌),
# 那两条在 27 个分子上照样抓得住错。
SMOKE=harness/baseline/smoke.bounds.jsonl
step "判官:三角光滑化 vs RDKit"
cargo run -q -p omgkit-conf --release --example smooth_oracle -- "$SMOKE"
step "判官:界矩阵(三条)"
cargo run -q -p omgkit-conf --release --example bounds_oracle -- "$SMOKE"
step "判官:特征分解 vs LAPACK + 精确回嵌"
cargo run -q -p omgkit-conf --release --example eigen_oracle -- "$SMOKE" harness/baseline/smoke.gram_eigs.jsonl

# 头号指标:界不可行的分子占比。跑**全语料 8831 个**,不跑冒烟档 ——
# 400 个样本上真实率 0.34% 只对应 1.4 个分子,泊松噪声足以让闸随机红绿。
# 语料随仓库入库(342 K),全程 0.7 秒。
step "判官:全语料界可行率(头号指标)"
cargo run -q -p omgkit-conf --release --example feasibility -- harness/corpus/large.smi

# **通用性难例语料。** large.smi 是药物样分子,在它上面全绿只说明"对药物样分子成立"。
# 这一份是照着算法的假设挑的:笼状/张力环、超配位、累积双键、超大环、少见元素、
# 金属、自由基、两性离子。一类分子在这里红了,答案必须是补一行约束表,不是加分支。
# 68 个分子,闸与全量档同一条(0.12%,对这个规模等于**一个都不许有**)。
# 这两步现在各跑**九条**:界可行(空区间 / 不可行)+ 硬不变量(原子重合 /
# 非有限数 / 该出构型没出)+ 几何(1-2 键 / 1-3 角 / 断键分子 / 键交叉分子)。
# 几何那四条先前只有端到端那条判官在看,而它跑的是 150 个药物样分子 ——
# 闸有、会让它红的数据也有,两者从没见过面。
# "该出没出"那条堵的是分母:几何四条的计数器都在构型生成成功之后才累加,
# 不给它配闸,任何让生成失败率上升的回归都会让几何闸变得更好看。
step "判官:难例语料(通用性 + 硬不变量)"
cargo run -q -p omgkit-conf --release --example feasibility -- harness/corpus/hard.smi
# 自穿:先拿真实构象校准检测器(必须报 0),再量我们自己的。
# 反过来做是自证 —— 检测器要是根本报不出东西,那个 0 只说明它没在看。
step "判官:自穿(先校准检测器,再量自己)"
cargo run -q -p omgkit-conf --release --example threading_oracle -- "$SMOKE"
step "判官:手性中心(真值取自真实构象)"
cargo run -q -p omgkit-conf --release --example chiral_oracle -- harness/baseline/smoke.chirality.jsonl

# **端到端。** 前面各条守一段,这一条守产物:分子进去、坐标出来,那组坐标满不满足化学。
# 精修前后各量一遍 —— 只报"之后"看不出精修有没有在干活。
step "判官:端到端构型(产物好不好)"
cargo run -q -p omgkit-conf --release --example conformer_oracle -- harness/baseline/smoke.chirality.jsonl

# **三配位立体中心**(亚砜/亚磺酰胺的 S、膦的 P:三根键 + 一对孤对)。
# 单独一条,因为上面那份基准里**一个这样的中心都没有** —— 于是这一档的
# 槽位约定在 CI 里从来没被验过:变异验证过,把三配位的槽位前两个对调
# (= 交付全部三配位中心的对映体),上面那条闸与全部单元测试**照样全绿**。
#
# 真值取自 RDKit 的**嵌入器**(它的 `AssignStereochemistryFrom3D` 读不回三配位 P,
# 但嵌入器认),号跨 seed 不稳的中心不进基准。
step "判官:三配位立体中心(孤对那一档)"
cargo run -q -p omgkit-conf --release --example conformer_oracle -- harness/baseline/smoke.lonepair.jsonl

# ---- 要外部实现(RDKit)的那四条 ----
#
# 上面几条都是拿预先烘好的基准比,所以不需要 RDKit。这四条不一样:它们把
# **当次**画出来 / 嵌出来的东西交给 RDKit 反读,基准没法预先烘。
#
# CI 里这四条在单独一个 job(`external`)里,版本钉在 `harness/requirements.lock`
# (RDKit 2025.09.2 —— 仓库里 `harness/baseline/` 那批基准就是它导的)。
# 开发机的 `.venv` 眼下是 2022.09.5,与 CI 不同:**这四条判据两边喂的是同一个
# RDKit**,版本变化会对消,两版都实测过退 0。判据自己会打印版本号,别靠记。
# 要跟 CI 完全对版就照 lock 重建 `.venv`。

# 楔形是"试 Up/Down、取反读回来对的那一个"构造出来的,而反读用的就是我们自己的
# `read_chirality` —— 拿它们往返是空过的。要问的是"别人照着这张图读,读出来是不是
# 同一个分子",那就必须把图交出去。
#
# 第二个参数是同一份语料:判据拿它核**分母**(dump 少喂几个分子,每一档都会
# 变好看 —— 实测空文件进去,先前那版打印"全部通过"并退 0)。
step "判官:楔形反读(别人照着图读构型)"
cargo run -q -p omgkit-depict --release --example dump_molblock -- harness/corpus/large.smi >"$WORK/blocks.txt"
"$PY" harness/check_wedge_readback.py "$WORK/blocks.txt" harness/corpus/large.smi

# 交付的三维坐标满不满足输入 SMILES 指定的每一处立体。完全绕开我们自己的任何公式。
# 第二个参数照例是语料,判据拿它核分母(实测空文件进去,先前那版打印
# "0/0 一致(0.00%)"并退 0)。
step "判官:交付坐标的立体化学(RDKit 从三维坐标读回)"
cargo run -q -p omgkit-conf --release --example dump_conformers -- harness/corpus/large.smi >"$WORK/ours.jsonl"
"$PY" harness/verify_stereo.py "$WORK/ours.jsonl" harness/corpus/large.smi

# 写出的外部裁判。**两个方向都跑** —— 规范那一条先前红了很久没人知道
# (规范写出丢掉超价原子的方括号,`Cl[I]Cl` → `ClICl`,外部读者补氢读成另一个分子),
# 而按存储顺序那一条一直是绿的:两个方向走不同分支,只跑一个等于只守一半。
#
# **`--strict` 是必须的。** 不加的话,判据会把"尚未写出的立体信息"分桶豁免,
# 而那个桶是**两侧一起抹掉**再比的 —— 于是"没写出 E/Z"和"E/Z 写反了"混成一档,
# 而且那个桶没有上限。独立审核实测:把写出器的单键方向符号一律写成 `/`
# (把全部顺式写成反式),判据打印"仅 双键立体 不同 149 条"然后**退 0**;
# 同样手法翻四面体手性则报 300 条分歧、退 1 —— 是这一档的洞,不是判官坏了。
# 大语料上两个豁免桶现值都是 0,所以 `--strict` 现在就能开(实测两个方向都退 0)。
step "判官:SMILES 写出(按存储顺序,严格)"
cargo run -q -p omgkit-io --release --example write_smiles -- harness/corpus/large.smi >"$WORK/written.tsv"
"$PY" harness/check_write.py "$WORK/written.tsv" harness/corpus/large.smi --strict
step "判官:SMILES 写出(规范,严格)"
cargo run -q -p omgkit-io --release --example write_smiles -- harness/corpus/large.smi --canonical >"$WORK/canon.tsv"
"$PY" harness/check_write.py "$WORK/canon.tsv" harness/corpus/large.smi --strict

# **自查。** 加了闸门忘了改 `TOTAL` 的话,这里红 —— 上面那些 `N/TOTAL`
# 就不会悄悄变成假数。
if [ "$N" -ne "$TOTAL" ]; then
    echo "闸门数对不上:实际跑了 $N 步,而 TOTAL 写着 $TOTAL" >&2
    exit 1
fi

echo
echo "$TOTAL 道闸全过。"
