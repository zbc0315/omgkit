#!/bin/bash
#
# 推之前在本地跑一遍全部闸门。**权威仍然是 `.github/workflows/ci.yml`** ——
# 这个脚本只是让本地跑一次不用手抄十条命令,两边不一致时以 CI 为准。
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


echo "== 1/8 fmt"
cargo fmt --all --check
echo "== 2/8 clippy(警告即失败)"
cargo clippy -q --workspace --all-targets -- -D warnings
echo "== 3/8 测试(release)"
cargo test -q --release
echo "== 4/8 测试(debug —— 让 debug_assert 真的跑到)"
cargo test -q --workspace
echo "== 5/8 文档"
cargo doc -q --workspace --no-deps --document-private-items

# ---- 三个外部判官 ----
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
echo "== 6/8 判官:三角光滑化 vs RDKit"
cargo run -q -p omgkit-conf --release --example smooth_oracle -- "$SMOKE"
echo "== 7/8 判官:界矩阵(三条)"
cargo run -q -p omgkit-conf --release --example bounds_oracle -- "$SMOKE"
echo "== 8/8 判官:特征分解 vs LAPACK + 精确回嵌"
cargo run -q -p omgkit-conf --release --example eigen_oracle -- "$SMOKE" harness/baseline/smoke.gram_eigs.jsonl

echo
echo "八道闸全过。"
