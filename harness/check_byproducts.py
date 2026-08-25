"""副产物收口的外部裁判:`capped` 那一档可以由外部实现**独立算出来**。

# 为什么这一档能被独立算出来

收口分两种。`bonded` 那一档要挑"成在哪两个原子之间",挑法是启发式的,外部实现
没有理由给出同一个答案;而 `capped` 那一档没有选择余地 —— 断口全是单键、不涉及
电荷、空价数正好等于要补的氢数,那就只有一件事可做:**在断口处切开、断头补氢**。
外部实现照样做得了,而且做的是同一件事。

于是这一档整个变成一次差分测试,判据与本实现无共谋。

# 它比手写用例强在哪

Rust 侧的 `tests/byproduct.rs` 靠十几条手写的规范 SMILES。手写用例有两个天花板:
数量,以及**作者想得到的形状**。这里拿真实语料跑上万条,而且**连立体化学一起比**
—— 手写用例里"一对对映体不该给出同一个副产物"这种事,不专门去想是写不出来的。

# 判据必须先证明自己不会空过

`--self-test` 往本实现的答案里**注入缺陷**,逐条确认判据抓得住:翻一个手性标记、
删一个原子、多补一个氢。抓不住就说明这一档比对是走过场。

另有 `--min-checked`:真正比过的条数低于它就直接失败。语料换了、口径改了都可能让
这一档悄悄变空,而"零分歧"在那时依然成立 —— 那是最会骗人的一种绿。

```bash
python3 harness/check_byproducts.py <templates.jsonl> --limit 3000
python3 harness/check_byproducts.py <templates.jsonl> --self-test
```

语料不随本仓库分发(见 `README.md` 的"USPTO-50k"一节),所以路径要自己给。
"""

import argparse
import collections
import json
import sys

import omgkit
import rdkit
from rdkit import Chem, RDLogger


def _banner() -> None:
    """**把版本与 wheel 的来路打出来。**

    这条判据是**经 wheel** 看 Rust 侧行为的,而 `import omgkit` 未必装的是
    刚建的那一份(用户级 site-packages 里可能躺着一个旧的)。改完 Rust 不重建
    就量,结论会稳稳指向错误的方向且毫无迹象 —— 这个坑本仓库栽过。

    RDKit 版本也要打:仓库钉 2025.09.2,开发机的 `.venv` 未必是同一个,
    而这一档**换版本会翻结论**(见 harness/README.md)。
    """
    print(f"外部实现:RDKit {rdkit.__version__}")
    print(f"  omgkit wheel:{omgkit.__file__}")


RDLogger.DisableLog("rdApp.*")

_PARAMS = Chem.SmilesParserParams()
_PARAMS.removeHs = False


def canon(smi):
    """交给外部实现规范化。读不了返回 None。"""
    m = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(m) if m is not None else None


def oracle(smi, discarded):
    """外部实现独立算出的副产物:在断口切开、断头补氢。

    判据不适用时返回 None —— **不适用与不一致要分开**,混起来的话判据会靠
    "大部分都不适用"显得很干净。
    """
    m = Chem.MolFromSmiles(smi, _PARAMS)
    if m is None:
        return None
    gone = set(discarded)
    cut = []
    for b in m.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if (i in gone) != (j in gone):
            if b.GetBondType() != Chem.BondType.SINGLE:
                return None                    # 只认单键断口,理由见模块文档
            cut.append(b.GetIdx())

    if not cut:
        if set(range(m.GetNumAtoms())) != gone:
            return None                        # 没有断口却也不是整分子被丢
        out = Chem.RWMol(m).GetMol()
    else:
        broken = Chem.FragmentOnBonds(m, cut, addDummies=True)
        n = m.GetNumAtoms()
        want = []
        for piece in Chem.GetMolFrags(broken):
            real = [a for a in piece if a < n]
            if not real:
                return None
            if set(real) <= gone:
                want.extend(piece)
            elif set(real) & gone:
                return None                    # 片段跨界,判据不适用
        if not want:
            return None
        rw = Chem.RWMol(broken)
        for a in rw.GetAtoms():
            if a.GetAtomicNum() == 0:          # 断头换成氢
                a.SetAtomicNum(1)
                a.SetIsotope(0)
                a.SetNoImplicit(False)
        for a in sorted(set(range(rw.GetNumAtoms())) - set(want), reverse=True):
            rw.RemoveAtom(a)
        out = rw.GetMol()
    try:
        Chem.SanitizeMol(out)
        return Chem.MolToSmiles(Chem.RemoveHs(out))
    except Exception:                          # noqa: BLE001
        return None


def mutate(smis, kind):
    """往本实现的答案里注入缺陷,用来证明判据抓得住。返回 None 表示这条注不进去。"""
    out = []
    touched = False
    for s in smis:
        m = Chem.MolFromSmiles(s)
        if m is None:
            return None
        if kind == "flip_stereo" and not touched:
            for a in m.GetAtoms():
                if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                    a.SetChiralTag(
                        Chem.ChiralType.CHI_TETRAHEDRAL_CW
                        if a.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CCW
                        else Chem.ChiralType.CHI_TETRAHEDRAL_CCW
                    )
                    touched = True
                    break
        elif kind == "drop_atom" and not touched and m.GetNumAtoms() > 1:
            rw = Chem.RWMol(m)
            rw.RemoveAtom(m.GetNumAtoms() - 1)
            m = rw.GetMol()
            touched = True
        elif kind == "extra_h" and not touched:
            for a in m.GetAtoms():
                if a.GetAtomicNum() == 1:
                    continue
                # 总氢必须在置 NO_IMPLICIT **之前**读:置上之后隐式那部分立刻
                # 不算数,读到 0 再加一就成了"把 CH3 改成 CH",而对 HCl 这类
                # 恰好还原成原样 —— 注入变成空操作,自证会误报"判据抓不住"。
                total = a.GetTotalNumHs()
                a.SetNoImplicit(True)
                a.SetNumExplicitHs(total + 1)
                touched = True
                break
        try:
            Chem.SanitizeMol(m)
            out.append(Chem.MolToSmiles(m))
        except Exception:                      # noqa: BLE001
            return None
    return out if touched else None


def main():
    _banner()
    ap = argparse.ArgumentParser()
    ap.add_argument("corpus", help="templates.jsonl(带 fwd 模板与 reactants)")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--min-checked", type=int, default=200,
                    help="真正比过的条数下限;低于它直接失败")
    ap.add_argument("--self-test", action="store_true",
                    help="注入缺陷,证明判据抓得住")
    args = ap.parse_args()

    stats = collections.Counter()
    caught = collections.Counter()
    injected = collections.Counter()
    bad = []

    with open(args.corpus) as f:
        for k, line in enumerate(f):
            if args.limit and k >= args.limit:
                break
            d = json.loads(line)
            try:
                rxn = omgkit.parse_reaction(d["fwd"])
                mols = []
                for s in d["reactants"]:
                    m = omgkit.parse_smiles(s)
                    m.sanitize()
                    mols.append(m)
            except Exception:                  # noqa: BLE001
                stats["读不了"] += 1
                continue
            if rxn.num_reactant_templates != len(mols):
                stats["位置式契约跑不了"] += 1
                continue
            try:
                outs = rxn.run(mols, byproducts=True)
            except Exception:                  # noqa: BLE001
                stats["run 抛异常"] += 1
                continue
            if not outs:
                stats["无输出"] += 1
                continue

            o = outs[0]
            if o.byproduct_verdict != "capped":
                stats["非 capped"] += 1
                continue
            b = o.byproduct_budget
            if b["charge_shift"] != 0 or b["need"] != b["open_valence"]:
                stats["涉及电荷/摘氢,判据不适用"] += 1
                continue

            want, ok = [], True
            for i, smi in enumerate(d["reactants"]):
                if not o.discarded[i]:
                    continue
                w = oracle(smi, o.discarded[i])
                if w is None:
                    ok = False
                    break
                want.extend(w.split("."))
            if not ok:
                stats["判据不适用"] += 1
                continue

            got = []
            for mol in o.byproducts:
                c = canon(mol.to_canonical_smiles())
                if c is None:
                    ok = False
                    break
                got.extend(c.split("."))
            if not ok:
                stats["副产物读不了"] += 1
                continue

            stats["比过"] += 1
            if sorted(got) == sorted(want):
                stats["一致"] += 1
            else:
                stats["分歧"] += 1
                if len(bad) < 10:
                    bad.append((d["id"], sorted(want), sorted(got)))

            if args.self_test:
                for kind in ("flip_stereo", "drop_atom", "extra_h"):
                    m = mutate(got, kind)
                    if m is None:
                        continue
                    injected[kind] += 1
                    if sorted(m) != sorted(want):
                        caught[kind] += 1

    for key, v in stats.most_common():
        print(f"  {v:7d}  {key}")
    if bad:
        print("\n分歧样例:")
        for rid, w, g in bad:
            print(f"  {rid}\n    判据 {w}\n    本库 {g}")

    failed = False
    if stats["比过"] < args.min_checked:
        print(f"\n判据空了:只比过 {stats['比过']} 条,下限是 {args.min_checked}")
        failed = True
    if stats["分歧"]:
        failed = True

    if args.self_test:
        print("\n注入缺陷 → 判据抓到:")
        for kind in ("flip_stereo", "drop_atom", "extra_h"):
            n, c = injected[kind], caught[kind]
            rate = f"{100 * c / n:.1f}%" if n else "—"
            print(f"  {kind:12s} 注入 {n:6d}  抓到 {c:6d}  {rate}")
            # 翻手性只在**本来带手性**的答案上注得进去,注不进的不算数;
            # 注得进却抓不住,说明这一档比对是走过场
            if n and c < n:
                print(f"    ✗ {kind} 有 {n - c} 条注入后判据仍说一致")
                failed = True

    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
