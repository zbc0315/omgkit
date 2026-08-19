"""把真实分子的 Gram 矩阵特征值导出来,给 `omgkit-conf` 的特征分解当**外部判官**。

# 为什么需要外部判官

手写的 Jacobi 拿"重建残差 `A = ΣλᵥᵥT`"和"正交性 `VVᵀ = I`"自验是不够的:
那两条只说明**分解是自洽的**,不说明特征值是**对的**。
numpy 的 `linalg.eigh` 走的是 LAPACK(分治 / QR),与 Jacobi 是完全不同的算法,
它不知道我们怎么写的,它只认数。

# 导两套矩阵,因为它们病得不一样

- `eig_u`:拿**光滑化之后的上限矩阵 `U`** 当距离表建的 Gram。
  这是真实流水线上喂给嵌入的那张表,谱里带着负特征值。
- `eig_x`:拿 **MMFF 优化过的真实构象**的精确距离建的 Gram。
  这一套的正确答案还有第三条独立的验法 —— 它必须恰好有三个非零特征值。

Gram 矩阵的公式两边写的是同一个(见 `crates/omgkit-conf/src/embed.rs`),
所以这里比的是**特征分解**,不是公式;公式那一层由 embed.rs 里
"真实三维点集精确回嵌"那条单元测试单独钉死。

用法:

    python3 harness/dump_gram.py harness/baseline/rdkit_bounds.jsonl out.jsonl
"""

import json
import sys

import numpy as np


def gram(d: np.ndarray) -> np.ndarray:
    """由距离表算关于质心的 Gram 矩阵。与 embed.rs 的 `metric_matrix` 同一公式。"""
    n = d.shape[0]
    sq = d * d
    # Σ_{j<k} d_jk² / n² —— 只走上三角
    sum_sq = np.triu(sq, 1).sum() / (n * n)
    sq0 = sq.sum(axis=1) / n - sum_sq
    return 0.5 * (sq0[:, None] + sq0[None, :] - sq)


def main() -> int:
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    src, out = sys.argv[1], sys.argv[2]
    limit = int(sys.argv[3]) if len(sys.argv) > 3 else 10**9

    n_ok = 0
    with open(out, "w", encoding="utf-8") as fh:
        for line in open(src, encoding="utf-8"):
            rec = json.loads(line)
            if n_ok >= limit:
                break
            n = rec["n"]
            sm = np.array(rec["smoothed"], dtype=float).reshape(n, n)
            # 上三角是上限;补成对称的整张距离表
            u = np.triu(sm, 1)
            u = u + u.T
            eig_u = np.linalg.eigvalsh(gram(u))[::-1]

            eig_x = None
            if rec.get("coords"):
                x = np.array(rec["coords"], dtype=float)
                dx = np.linalg.norm(x[:, None, :] - x[None, :, :], axis=-1)
                eig_x = np.linalg.eigvalsh(gram(dx))[::-1]

            fh.write(
                json.dumps(
                    {
                        "smiles": rec["smiles"],
                        "n": n,
                        "eig_u": [float(v) for v in eig_u],
                        "eig_x": None if eig_x is None else [float(v) for v in eig_x],
                    }
                )
                + "\n"
            )
            n_ok += 1
    print(f"导出 {n_ok} 个分子的特征值(numpy {np.__version__} / LAPACK)→ {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
