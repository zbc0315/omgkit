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
TOTAL=32
N=0
step() {
    N=$((N + 1))
    echo "== $N/$TOTAL $1"
}

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

# **末尾那一批判据要 RDKit,所以先在这里查一次。** 放在末尾的话,要等前面十几步
# (十来分钟 cargo)跑完才知道环境缺东西。没有就直接失败,**不跳过** ——
# 静默跳过的判据是最坏的一种,它让人以为跑过了。
PY=${PY:-.venv/bin/python}
if [ ! -x "$PY" ]; then
    echo "缺 $PY —— 末尾那一批判据要 RDKit。" >&2
    echo "  建法:python3 -m venv .venv &&" >&2
    echo "        .venv/bin/pip install --only-binary=:all: -r harness/requirements.lock" >&2
    exit 1
fi

# **版本也要查,而且要在这里查。**
#
# 先前这里只查"有没有这个解释器",理由是"这批判据两边喂的是同一个 RDKit,
# 版本变化会对消"。那句话对**当时**那几条成立(参照侧与读回侧都是 RDKit),
# 对后来接进来的**不成立**:`check_smarts_chirality.py` 一侧是 RDKit、
# 另一侧是本实现,版本一换就没得对消。
#
# 实测:同一批 SMARTS 查询,RDKit 2022.09.5 与 2025.09.2 给出**相反**的匹配
# (`[C@@H](C)(N)O` 在 2025 上匹配 `C[C@H](N)O`,在 2022 上匹配它的对映体),
# 而两版对同一串**当 SMILES 读**完全一致 —— 2022 的 SMARTS 与它自己的 SMILES
# 读法自相矛盾,2025 修好了。拿 2022 跑,这条判据报 48 条"反了",一条都不是真的。
#
# 所以:版本对不上就**当场停**,不往下跑。读一个跑错版本的绿或红,比不跑更糟。
WANT_RDKIT=$(sed -n 's/^rdkit==//p' harness/requirements.lock)
HAVE_RDKIT=$("$PY" -c 'import rdkit; print(rdkit.__version__)' 2>/dev/null || echo "(装不上)")
# 锁里写 2025.9.2,而 `rdkit.__version__` 报 2025.09.2 —— 补零的差别,按数值比
norm() { echo "$1" | awk -F. '{printf "%d.%d.%d", $1, $2, $3}'; }
if [ "$(norm "$HAVE_RDKIT")" != "$(norm "$WANT_RDKIT")" ]; then
    echo "${PY} 的 RDKit 是 ${HAVE_RDKIT},而 harness/requirements.lock 钉的是 ${WANT_RDKIT}。" >&2
    echo "  仓库只认一个 RDKit 版本 —— 版本不对时,一侧是本实现的那几条判据" >&2
    echo "  (check_smarts_chirality)量出来的红绿都不作数。" >&2
    echo "  对版:$PY -m pip install --force-reinstall --only-binary=:all: -r harness/requirements.lock" >&2
    echo "  (若 pip list 里同时有 rdkit 与 rdkit-pypi,先卸掉 rdkit-pypi —— 两者装进同一个包目录,后装的赢)" >&2
    echo "  只想拿别的解释器跑一遍:PY=/path/to/python bash harness/gates.sh" >&2
    exit 1
fi

# **Indigo 也要钉。** 它只服务丙二烯轴手性那一条判据(RDKit 在那一档完全没有
# 能力,见 `harness/requirements.lock` 里的理由),而那条判据的一侧是本实现 ——
# 版本一换就没得对消,与上面 RDKit 同理。
WANT_INDIGO=$(sed -n 's/^epam\.indigo==//p' harness/requirements.lock)
HAVE_INDIGO=$("$PY" -c 'from indigo import Indigo; print(Indigo().version().split("-")[0])' 2>/dev/null || echo "(装不上)")
if [ "$(norm "$HAVE_INDIGO")" != "$(norm "$WANT_INDIGO")" ]; then
    echo "${PY} 的 Indigo 是 ${HAVE_INDIGO},而 harness/requirements.lock 钉的是 ${WANT_INDIGO}。" >&2
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

# ---- 要外部实现(RDKit)的那一批 ----
#
# 上面几条都是拿预先烘好的基准比,所以不需要 RDKit。下面这一批不一样:它们把
# **当次**画出来 / 嵌出来的东西交给 RDKit 反读,基准没法预先烘。
#
# CI 里这一批在单独一个 job(`external`)里,版本钉在 `harness/requirements.lock`
# (RDKit 2025.09.2 —— 仓库里 `harness/baseline/` 那批基准就是它导的)。
# 开发机的 `.venv` 眼下是 2022.09.5,与 CI 不同:**这一批判据两边喂的是同一个
# RDKit**,版本变化会对消,两版都实测过退 0。判据自己会打印版本号,别靠记。
# 要跟 CI 完全对版就照 lock 重建 `.venv`。

# **基准与生成它的脚本脱钩了没有。** 两档:
#
# - **结构档**(手性/孤对/界):数值各钉在各自的 RDKit 版本上(手性那份是
#   2022.09.5、界那份是 2025.09.2),逐字节比会天天红;而结构与版本无关
#   (实测:同一批脚本在两个版本下结构逐个相同,只是 lonepair 那份的分子数
#   15 vs 16 —— 值会变,结构不会)。
# - **逐字节档**(`smoke.l3.jsonl`):`oracle_pipeline.py` 的 l3 只吐字符串,
#   不含坐标与嵌入,同版本下跨平台逐字节相同(实测 macOS-arm64 与
#   Linux-x86_64 上 sha256 一致)。这一份实测有过两处脱钩,结构档一处也看不见:
#   失败记录的 `err` 字串变了(键没变),以及 `--remove-hs` 放在解析那一步
#   把三条手性用例换成了另一个分子(值变了,键没变)。
#
# **判据红了会把重导命令一起打出来**(命令实测能逐字节复现入库的那份)。
#
# 这一条守的是一个真发生过、四个月没人看见的洞:`61b8d58` 教会
# `dump_chirality.py` 收三配位立体中心,却没重导 `smoke.chirality.jsonl` ——
# 那个提交声称落地的那一档,在主手性判官眼里根本不存在。
step "判官:入库基准与生成它的脚本没脱钩(结构 + l3 逐字节)"
"$PY" harness/check_baseline_schema.py

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
# ---- 经 wheel 看 Rust 行为的那几条 ----
#
# 门槛比别的高一层:它们 `import omgkit`,看到的是**建出来的 wheel**,不是源码。
# 所以要先建再装。用户级 site-packages 里躺着旧的一份时,`import omgkit` 会
# 静默拿到它 —— 判据自己会把 `omgkit.__file__` 打出来。
#
# `check_byproducts.py` 进不来:它要 USPTO-50k 的 templates.jsonl,
# 那份语料不随仓库分发。
echo "== 建 wheel(下面几条判据经它看 Rust 行为)"
"$PY" -m maturin build --release -q -m crates/omgkit-py/Cargo.toml --out "$WORK/wheels"
"$PY" -m pip install -q --force-reinstall "$WORK"/wheels/omgkit-*.whl

step "判官:规范化的自指不变量(不动点 + 幂等)"
"$PY" harness/check_canonical_fixpoint.py harness/corpus/large.smi

# **配位几何的排列表每次重新量一遍。** `polyhedron.rs` 里那三张表是从 RDKit
# 穷举量出来的(72 / 2400 / 21600 种写法),进了源码之后就没有任何东西再核对它。
# "缺的那个顶点(方括号里的氢,或一个空的配位位置)插在列出顺序哪一位"也是
# 这么量出来的,一并在这里重量:9 族共 40224 种写法。
# 这条判据比的是**分组**(哪些写法落到同一个分子),不是字符串 —— 两个实现的
# 规范串本来就不一样,该一致的是"谁和谁是同一个分子"。
# 注意它有个盲区:分组看不见"给配体全局换个名字"。所以每一族里都放了几种
# 写法(起笔位置、环的书写方向),把被换的那两个配体的角色拆开 —— 变异实测
# 见 `harness/README.md`。
step "判官:配位几何的排列分组(与外部实现逐组比)"
"$PY" harness/check_stereo_perm.py
# 这一档**换 RDKit 版本会翻结论**:2022.09.5 与 2025.09.2 对同一批查询给出相反
# 的匹配,而两版对同一串当 SMILES 读完全一致 —— 2022 的 SMARTS 与它自己的
# SMILES 读法自相矛盾。仓库钉 2025.09.2;判据自己打印版本号。
# 开发机的 .venv 若是 2022.09.5,这一条会红 48 条,那不是回归。
# **丙二烯型轴手性的裁判是 Indigo,不是 RDKit。** 后者在这一档上完全没有能力
# (六条路实测都把 `@AL1` 与 `@AL2` 当成同一个东西)。这条判据比的同样是
# **分组**,而且每族都放了几种把配体角色拆开的写法 —— 变异实测见
# `harness/README.md`:少了那几种写法,"把端上的氢排到末尾"这条变异是全绿的。
step "判官:丙二烯轴手性的分组(与 Indigo 逐组比)"
"$PY" harness/check_allene.py
step "判官:SMARTS 手性的参照系(自带区分力闸)"
"$PY" harness/check_smarts_chirality.py
step "判官:产物侧手性的四种指令"
"$PY" harness/check_product_chirality.py
step "判官:Python 绑定"
"$PY" harness/test_python.py

# **构型生成的绑定:与 Rust 侧逐位比。** 绑定那一层"只做翻译",可翻译本身也会
# 错(原子表对不上、坐标错位、忘了生成时补过显式氢),而那种错只有 Python 用户
# 碰得到,Rust 侧的判据一概盖不到。两侧调的是同一个 `pipeline::conformer_for`、
# 全程无随机数,所以"逐位相同"是可以要求的。
# 变异实测:让绑定返回补氢**之前**那份分子(一个很自然的错),642 条里 641 条红。
step "判官:构型生成的 Python 绑定(与 Rust 逐位比)"
"$PY" harness/check_python_conformer.py "$WORK/ours.jsonl"

# ---- 先前只在本地手动跑的那几条 ----
#
# **闸不进 CI 就不是闸。** 下面这四条判官一直躺在 `harness/` 里,谁想起来谁跑一次
# —— 而接进来的那一刻,`check_bond_stereo.py` **就是红的**(小环双键那一条,
# 见下),`check_write.py --strict` 拿冒烟语料跑也是红的(非四面体立体写不出来,
# 现在按 SMILES 逐条钉死在 `NON_TETRAHEDRAL_GAP` 里)。
#
# 四条都补了**分母闸**(`harness/denominator.py`):它们原本只数分歧、
# 不数"该数到多少",喂个空文件进去打印一片空白然后退 0。
step "判官:双键顺反的感知(与外部实现逐根比)"
cargo run -q -p omgkit-io --release --example dump_bond_stereo -- harness/corpus/large.smi >"$WORK/bs.tsv"
"$PY" harness/check_bond_stereo.py "$WORK/bs.tsv" harness/corpus/large.smi
step "判官:写出时 E/Z 守不守恒"
cargo run -q -p omgkit-io --release --example dump_written -- harness/corpus/large.smi >"$WORK/dw.tsv"
"$PY" harness/check_ez.py "$WORK/dw.tsv" harness/corpus/large.smi
step "判官:净化之后写出忠不忠实"
cargo run -q -p omgkit-chem --release --example dump_sanitized -- harness/corpus/large.smi >"$WORK/san.tsv"
"$PY" harness/check_write_fidelity.py "$WORK/san.tsv" harness/corpus/large.smi
step "判官:SMARTS 写出(语义相同,自带区分力闸)"
cargo run -q -p omgkit-io --release --example dump_smarts_written -- harness/corpus/smarts.txt >"$WORK/sw.tsv"
"$PY" harness/check_smarts_write.py "$WORK/sw.tsv" --mols harness/corpus/large.smi --corpus harness/corpus/smarts.txt

step "判官:SMILES 写出(按存储顺序,严格)"
cargo run -q -p omgkit-io --release --example write_smiles -- harness/corpus/large.smi >"$WORK/written.tsv"
"$PY" harness/check_write.py "$WORK/written.tsv" harness/corpus/large.smi --strict
step "判官:SMILES 写出(规范,严格)"
cargo run -q -p omgkit-io --release --example write_smiles -- harness/corpus/large.smi --canonical >"$WORK/canon.tsv"
"$PY" harness/check_write.py "$WORK/canon.tsv" harness/corpus/large.smi --strict

# **产物生成。** 这条判官零容差,而它有 22 条**刻意分歧**(产物模板描述的是
# 反应中心的片段,不是"一个片段一个分子";环状底物 + 断环模板给出一个开环产物,
# 外部实现逐产物各搬一次,原子凭空变多)—— 于是它一直进不了 CI,谁想起来谁跑
# 一次。现在按(反应模板, 底物)**逐条钉死**,名单两个方向都红,才接得进来。
#
# 没接进来的这段时间里,文档写的 717 / 24 悄悄变成了 716 / 25,而多出来的那条
# 根本不是"刻意分歧":双键顺反的参照原子被反应删掉、又没人顶它的槽位时,
# 整根键的顺反被作废了。已修(见 reaction.rs 的
# `bond_stereo_rebases_to_the_other_side_when_nothing_fills_the_slot`),
# 现在是 719 / 22。
step "判官:产物生成(刻意分歧逐条钉死)"
cargo run -q -p omgkit-match --release --example dump_reactions -- harness/corpus/reactions.txt harness/corpus/large.smi 300 >"$WORK/rx.tsv"
"$PY" harness/check_reactions.py "$WORK/rx.tsv" --rxns harness/corpus/reactions.txt --mols harness/corpus/large.smi

# **冒烟语料也要跑写出。** CI 先前只拿 `large.smi` 跑这条判官,而那份语料里
# 一条非四面体立体(`@SP` / `@TB` / `@OH`)都没有 —— 于是"读得回来、写不出去"
# 这件事**一条判据都没守**。冒烟语料里有 6 条,现在按 SMILES 逐条钉死:
# 写出器补上任何一条这里当场红,逼着把它从名单里划掉。
step "判官:SMILES 写出(冒烟语料,严格 —— 非四面体立体那一档钉在这里)"
cargo run -q -p omgkit-io --release --example write_smiles -- harness/corpus/smoke.smi >"$WORK/wsmoke.tsv"
"$PY" harness/check_write.py "$WORK/wsmoke.tsv" harness/corpus/smoke.smi --strict

# **自查。** 加了闸门忘了改 `TOTAL` 的话,这里红 —— 上面那些 `N/TOTAL`
# 就不会悄悄变成假数。
if [ "$N" -ne "$TOTAL" ]; then
    echo "闸门数对不上:实际跑了 $N 步,而 TOTAL 写着 $TOTAL" >&2
    exit 1
fi

echo
echo "$TOTAL 道闸全过。"
