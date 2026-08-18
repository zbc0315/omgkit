"""把 RDKit 的界矩阵导出成 JSONL,给 `omgkit-conf` 的三角光滑化当**外部判官**。

每个分子导两张矩阵:`raw`(未光滑)与 `smoothed`(RDKit 自己光滑过的)。
判据的用法是:把 `raw` 喂给我们的 `triangle_smooth`,结果与 `smoothed` 逐元素比。

**为什么要外部判官**:自己写的光滑化拿自己的性质去验(比如"再光滑一次不变"),
两边同错就同绿。RDKit 不知道我们怎么写的,它只认数。

约定与 RDKit 一致:**上三角是上限、下三角是下限**。

用法:

    python3 harness/dump_bounds.py harness/corpus/large.smi out.jsonl [分子数上限]

注意:装在系统 python 里的 rdkit 可能是旧版且与 numpy 2 冲突。
本仓用的是 2025.09.2(见 harness/embed_ablation.py 顶部那段版本说明)。
"""

import json
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import rdDistGeom

RDLogger.DisableLog("rdApp.*")


def main() -> int:
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    corpus, out = sys.argv[1], sys.argv[2]
    limit = int(sys.argv[3]) if len(sys.argv) > 3 else 10**9

    n_ok = 0
    n_skip = 0
    with open(out, "w", encoding="utf-8") as fh:
        for line in open(corpus, encoding="utf-8"):
            smi = line.split("\t")[0].strip()
            if not smi:
                continue
            if n_ok >= limit:
                break
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                n_skip += 1
                continue
            mol = Chem.AddHs(mol)
            try:
                raw = rdDistGeom.GetMoleculeBoundsMatrix(mol, doTriangleSmoothing=False)
                sm = rdDistGeom.GetMoleculeBoundsMatrix(mol, doTriangleSmoothing=True)
            except Exception:  # noqa: BLE001 —— RDKit 会抛 C++ 异常
                n_skip += 1
                continue
            fh.write(
                json.dumps(
                    {
                        "smiles": smi,
                        "n": int(raw.shape[0]),
                        "raw": [float(x) for x in raw.ravel()],
                        "smoothed": [float(x) for x in sm.ravel()],
                    }
                )
                + "\n"
            )
            n_ok += 1
    print(f"导出 {n_ok} 个分子,跳过 {n_skip} 个 → {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
