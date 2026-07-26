"""四组基准:omgkit 正向 / omgkit 逆向 / rdkit 正向 / rdkit 逆向。

每条反应记一行 JSON:四组各自的运行时间、是否命中、输出个数,加上重原子数。

# 怎么才算公平

**同一份模板。** 逆向模板由 rdchiral 从这条反应本身抽出,正向模板是它两侧对调。
两个引擎拿到的是**同一个字符串**,各自解析。

**同一套输入。** 输入分子的 SMILES 是同一串(RDKit 规范式、无映射号),
各自用自己的解析器读进来、各自净化 —— 这两步都在计时之外。

**同一个裁判。** 命中与否一律用 **RDKit 的规范 SMILES** 判。omgkit 的产物先写成
SMILES,再交给 RDKit 解析并规范化;RDKit 的产物直接净化后规范化。这样判据不偏
任何一方 —— 用 omgkit 自己的规范式判 omgkit 反而是自证。代价是 omgkit 要多过
一道"写出 + 被 RDKit 读回"的关,写错了就算未命中。这个偏差是**对 omgkit 不利**
的,宁可如此。

**同样多的调用。** 正向模板有 N 个反应物模板,记录里有 M 个参与反应的分子。
RunReactants 要求两者数目相等且**按位置对应**,所以枚举 M 取 N 的全部有序组合,
两个引擎枚举同一批、一个不落(不因为提前命中就收工 —— 那会让"运气好"看着像
"跑得快")。每条记录的组合数一并记下,想换算成单次调用的耗时随时可以除。

# 计时口径

计时只包住 `run` / `RunReactants` 本身。输入解析、净化、模板解析、结果转 SMILES
全在计时之外。先跑一次热身,再取 `--reps` 次的**最小值**(最小值受操作系统噪声
污染最小),中位数也一并记下。

有一处**必须挑明的不对称**:omgkit 的 Python 版 `run` 每次调用都会
`MolProps::compute` 重算一遍分子的查询性质(环集、最小环大小等),还会把输入
分子深拷贝一份;RDKit 这些是在 `MolFromSmiles` 净化时算好、跨调用复用的。
也就是说这里量到的 omgkit 耗时含着一段 RDKit 不付的固定开销。
`measure_overhead.py` 单独把这段量出来,分析时会从中扣除给出一个"去掉 API
开销"的对照值 —— 但**表里的主数字仍是用户真实付出的那个**。
"""

import argparse
import itertools
import json
import os
import sys
import time
from statistics import median

from rdkit import Chem, RDLogger
from rdkit.Chem import rdChemReactions

import omgkit

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

PERM_CAP = 24  # M 取 N 的有序组合上限;超了就截断并在行里标出来
MAX_PRODUCTS = 1000  # 两侧同一个上限,防止对称分子把输出撑爆


# ---------------------------------------------------------------- 结果规范化


def strip_maps(m):
    for a in m.GetAtoms():
        if a.GetAtomMapNum():
            a.SetAtomMapNum(0)
    return m


def rd_canon_of_product(p):
    """RDKit 反应产物 -> RDKit 规范 SMILES。产物是未净化的,得自己净化。"""
    try:
        q = Chem.Mol(p)
        strip_maps(q)
        Chem.SanitizeMol(q)
        return Chem.MolToSmiles(q)
    except Exception:
        pass
    try:  # 净化不过时退一步:先写出再读回,读回时 RDKit 会自己净化
        q = Chem.Mol(p)
        strip_maps(q)
        s = Chem.MolToSmiles(q)
        m = Chem.MolFromSmiles(s)
        return Chem.MolToSmiles(m) if m is not None else None
    except Exception:
        return None


def og_canon_of_product(p):
    """omgkit 反应产物 -> (交给 RDKit 判的)规范 SMILES。"""
    s = None
    try:
        q = p.copy()
        q.sanitize()
        s = q.to_smiles()
    except Exception:
        try:
            s = p.to_smiles()
        except Exception:
            return None
    if s is None:
        return None
    m = Chem.MolFromSmiles(s)
    return Chem.MolToSmiles(m) if m is not None else None


def joined(smis):
    """一组产物 -> 一个可比的规范串。任一处失败则整组作废。"""
    if any(x is None for x in smis):
        return None
    m = Chem.MolFromSmiles(".".join(smis))
    return Chem.MolToSmiles(m) if m is not None else None


# ---------------------------------------------------------------- 计时


def timed(fn, reps):
    """跑一次热身,再跑 reps 次。返回 (最小值, 中位数, 最后一次的结果)。"""
    out = fn()
    ts = []
    for _ in range(reps):
        t0 = time.perf_counter()
        out = fn()
        ts.append(time.perf_counter() - t0)
    return min(ts), median(ts), out


def assignments(mols, n):
    """M 个分子取 N 个的全部有序组合。返回 (组合列表, 是否被截断)。"""
    if n <= 0 or len(mols) < n:
        return [], False
    perms = []
    for i, p in enumerate(itertools.permutations(mols, n)):
        if i >= PERM_CAP:
            return perms, True
        perms.append(p)
    return perms, False


# ---------------------------------------------------------------- 一条反应


def bench_one(rec, reps, diff_sink=None):
    row = {
        "id": rec["id"],
        "row": rec["row"],
        "cls": rec.get("cls", ""),
        "n_heavy_r": rec["n_heavy_r"],
        "n_heavy_p": rec["n_heavy_p"],
        "n_reactants": len(rec["reactants"]),
    }

    # ---- 模板解析(计时之外)----
    tpl = {}
    for tag, smarts in (("fwd", rec["fwd"]), ("retro", rec["retro"])):
        try:
            r = rdChemReactions.ReactionFromSmarts(smarts)
            r.Initialize()
            tpl[("rdkit", tag)] = r
        except Exception:
            tpl[("rdkit", tag)] = None
        try:
            tpl[("omgkit", tag)] = omgkit.parse_reaction(smarts)
        except Exception:
            tpl[("omgkit", tag)] = None
    row["tpl_parse"] = {
        f"{e}_{d}": tpl[(e, d)] is not None
        for e in ("rdkit", "omgkit")
        for d in ("fwd", "retro")
    }

    # ---- 输入分子(计时之外)----
    rd_in = {"fwd": [], "retro": []}
    og_in = {"fwd": [], "retro": []}
    ok = True
    for smi in rec["reactants"]:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            ok = False
            break
        rd_in["fwd"].append(m)
        try:
            o = omgkit.parse_smiles(smi)
            o.sanitize()
            og_in["fwd"].append(o)
        except Exception:
            ok = False
            break
    pm = Chem.MolFromSmiles(rec["prod"])
    if pm is None:
        ok = False
    else:
        rd_in["retro"].append(pm)
        try:
            po = omgkit.parse_smiles(rec["prod"])
            po.sanitize()
            og_in["retro"].append(po)
        except Exception:
            ok = False
    if not ok:
        row["err"] = "input_prep"
        return row

    # ---- 真值 ----
    truth = {
        "fwd": joined([rec["prod"]]),
        "retro": joined(rec["reactants"]),
    }
    row["truth_ok"] = {k: v is not None for k, v in truth.items()}

    # ---- 四组 ----
    preds_by_key = {}
    for direction in ("fwd", "retro"):
        for engine in ("omgkit", "rdkit"):
            key = f"{engine}_{direction}"
            t = tpl[(engine, direction)]
            if t is None:
                row[key] = {"err": "template_parse"}
                continue
            n_tpl = (
                t.GetNumReactantTemplates()
                if engine == "rdkit"
                else t.num_reactant_templates
            )
            pool = (rd_in if engine == "rdkit" else og_in)[direction]
            perms, capped = assignments(pool, n_tpl)
            if not perms:
                row[key] = {
                    "err": "shape",
                    "n_tpl": n_tpl,
                    "n_pool": len(pool),
                }
                continue

            if engine == "rdkit":

                def call(t=t, perms=perms):
                    out = []
                    for p in perms:
                        out.extend(t.RunReactants(p, MAX_PRODUCTS))
                    return out

            else:
                # 先转成 list 再进计时区:RDKit 直接吃元组,omgkit 的绑定要 list。
                # 把这个转换留在计时里等于让 omgkit 白背一次列表构造。
                perms = [list(p) for p in perms]

                def call(t=t, perms=perms):
                    out = []
                    for p in perms:
                        out.extend(t.run(p, max_products=MAX_PRODUCTS))
                    return out

            try:
                tmin, tmed, outcomes = timed(call, reps)
            except Exception as e:
                row[key] = {"err": f"run:{type(e).__name__}"}
                continue

            # ---- 命中判定(计时之外)----
            conv = rd_canon_of_product if engine == "rdkit" else og_canon_of_product
            preds = set()
            n_bad = 0
            for oc in outcomes:
                mols = oc if engine == "rdkit" else oc.products
                s = joined([conv(x) for x in mols])
                if s is None:
                    n_bad += 1
                else:
                    preds.add(s)
            row[key] = {
                "t_min": tmin,
                "t_med": tmed,
                "n_perm": len(perms),
                "capped": capped,
                "n_tpl": n_tpl,
                "n_out": len(outcomes),
                "n_uniq": len(preds),
                "n_bad": n_bad,
                "hit": (truth[direction] is not None and truth[direction] in preds),
            }
            preds_by_key[key] = preds

    # 只要 omgkit 没命中就把双方的预测集原样留一份 —— 两边判得不一样的要挑
    # 典型案例,两边都没中的要逐条归因。不留就得为了几百条重跑整个基准。
    if diff_sink is not None:
        for direction in ("fwd", "retro"):
            a, b = f"omgkit_{direction}", f"rdkit_{direction}"
            if a not in preds_by_key or b not in preds_by_key:
                continue
            if row[a]["hit"] and row[b]["hit"]:
                continue
            diff_sink.write(
                json.dumps(
                    {
                        "id": rec["id"],
                        "row": rec["row"],
                        "direction": direction,
                        "winner": (
                            "omgkit"
                            if row[a]["hit"]
                            else ("rdkit" if row[b]["hit"] else "none")
                        ),
                        "truth": truth[direction],
                        "retro": rec["retro"],
                        "fwd": rec["fwd"],
                        "prod": rec["prod"],
                        "reactants": rec["reactants"],
                        "omgkit_preds": sorted(preds_by_key[a]),
                        "rdkit_preds": sorted(preds_by_key[b]),
                    }
                )
                + "\n"
            )
    return row


# ---------------------------------------------------------------- 主流程


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default=os.path.join(ROOT, "data", "templates.jsonl"))
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "bench.jsonl"))
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--start", type=int, default=0)
    args = ap.parse_args()

    recs = []
    with open(args.src) as fh:
        for line in fh:
            r = json.loads(line)
            if "err" in r:
                continue
            recs.append(r)
    if args.start:
        recs = recs[args.start :]
    if args.limit:
        recs = recs[: args.limit]
    print(f"待测 {len(recs)} 条,reps={args.reps}", file=sys.stderr, flush=True)

    diff_path = args.out.replace(".jsonl", "") + ".diff.jsonl"
    t0 = time.time()
    with open(args.out, "w") as out, open(diff_path, "w") as diff:
        for i, rec in enumerate(recs):
            try:
                row = bench_one(rec, args.reps, diff)
            except Exception as e:
                row = {"id": rec["id"], "row": rec["row"], "err": f"fatal:{e!r}"}
            out.write(json.dumps(row) + "\n")
            if (i + 1) % 500 == 0:
                out.flush()
                el = time.time() - t0
                print(
                    f"  {i + 1}/{len(recs)}  {el:.0f}s  ({(i + 1) / el:.1f} 条/s)",
                    file=sys.stderr,
                    flush=True,
                )
    print(f"完成,耗时 {time.time() - t0:.0f}s", file=sys.stderr)


if __name__ == "__main__":
    main()
