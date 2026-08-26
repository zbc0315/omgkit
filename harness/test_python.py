#!/usr/bin/env python3
"""Python 绑定的测试。

# 为什么测试在这一侧

`omgkit-py` 是 cdylib,链接时把 Python 符号欠着由宿主解释器提供 —— `cargo test`
给它建测试可执行文件必然失败,没人提供那些符号。所以它的测试只能从 Python 跑,
而这也是更有意义的一侧:要验的正是"从 Python 看过去对不对"。

# 测的是翻译,不是化学

绑定层只做"把参数翻过去、把结果翻回来"。化学正确性由 Rust 侧的全套差分测试
守着,在这里重测一遍既慢又盖不全。这里盯的是翻译层特有的失效:类型翻错、
错误信息丢失、对象语义(拷贝/可变性)不对、以及构建配置埋下的雷。

用法:

    maturin build --release -m crates/omgkit-py/Cargo.toml --out <目录>
    pip install <目录>/omgkit-*.whl
    python3 harness/test_python.py
"""

import pathlib
import re
import sys
import unittest

try:
    import omgkit
except ImportError:  # pragma: no cover
    sys.exit(
        "导入不到 omgkit。先构建并安装 wheel:\n"
        "  maturin build --release -m crates/omgkit-py/Cargo.toml --out /tmp/wheels\n"
        "  pip install --force-reinstall /tmp/wheels/omgkit-*.whl"
    )

REPO = pathlib.Path(__file__).resolve().parent.parent


def corpus(name: str, limit: int):
    """取语料前若干条 SMILES。"""
    path = REPO / "harness" / "corpus" / name
    out = []
    for line in path.read_text(errors="replace").splitlines():
        tok = line.split()[0] if line.split() else ""
        if tok and not tok.startswith("#"):
            out.append(tok)
        if len(out) >= limit:
            break
    return out


class Translation(unittest.TestCase):
    """翻译层本身。"""

    def test_parse_write_roundtrip(self):
        for smi in ["CCO", "OC(=O)c1ccccc1N", "F/C=C/F", "N[C@@H](C)C(=O)O"]:
            m = omgkit.parse_smiles(smi)
            self.assertEqual(m.to_smiles(), smi, f"{smi} 按存储顺序写出应当逐字节相同")

    def test_atomic_nums_is_a_list_of_int(self):
        """`Vec<u8>` 会被 PyO3 特判成 `bytes`,那是个很安静的类型错误。"""
        m = omgkit.parse_smiles("CCO")
        got = m.atomic_nums
        self.assertIsInstance(got, list, f"应当是 list,实际 {type(got).__name__}")
        self.assertTrue(all(isinstance(z, int) for z in got), "元素应当是 int")
        self.assertEqual(got, [6, 6, 8])

    def test_num_atoms_and_bonds(self):
        m = omgkit.parse_smiles("c1ccccc1")
        self.assertEqual(m.num_atoms, 6)
        self.assertEqual(m.num_bonds, 6)

    def test_copy_is_deep(self):
        """`copy()` 之后改副本不能动到原件 —— 净化是就地修改的。"""
        a = omgkit.parse_smiles("c1cc[nH]c1")
        b = a.copy()
        b.sanitize()
        self.assertEqual(a.to_smiles(), "c1cc[nH]c1", "原件被副本的净化改到了")

    def test_repr_is_informative(self):
        self.assertIn("Mol", repr(omgkit.parse_smiles("CC")))

    def test_version_is_exposed(self):
        self.assertRegex(omgkit.__version__, r"^\d+\.\d+\.\d+")


class Errors(unittest.TestCase):
    """失败要翻成异常,而且**不能丢信息**。"""

    def test_parse_error_carries_the_position(self):
        with self.assertRaises(ValueError) as cm:
            omgkit.parse_smiles("C1CC")
        text = str(cm.exception)
        self.assertIn("C1CC", text, "报错里应当带上原文")
        self.assertIn("^", text, "报错里应当有指出位置的插字号")

    def test_sanitize_error_is_a_value_error(self):
        # 五价碳,价键校验必然拒绝
        m = omgkit.parse_smiles("C(C)(C)(C)(C)C")
        with self.assertRaises(ValueError):
            m.sanitize()

    def test_parse_error_is_not_a_bare_message(self):
        """笼统的一句 "invalid" 等于把排查线索丢掉。"""
        with self.assertRaises(ValueError) as cm:
            omgkit.parse_smiles("C(((")
        self.assertGreater(len(str(cm.exception)), 10)


class AgainstTheCorpus(unittest.TestCase):
    """在真实分子上走一遍,确认 Python 这条路通到的是同一套实现。"""

    def test_canonical_is_independent_of_how_the_molecule_was_written(self):
        """规范 SMILES 的定义性质,从 Python 这一侧再验一遍。

        用 omgkit 自己写出的另一种写法当输入 —— 不需要外部参照就能自证。
        """
        checked = 0
        for smi in corpus("large.smi", 400):
            try:
                a = omgkit.parse_smiles(smi)
                a.sanitize()
            except ValueError:
                continue  # 语料里本就有净化不了的,那不是本档要管的事
            canon_a = a.to_canonical_smiles()
            # 换一种写法:按存储顺序写出再读回来
            b = omgkit.parse_smiles(a.to_smiles())
            b.sanitize()
            self.assertEqual(canon_a, b.to_canonical_smiles(), f"{smi} 的规范串不稳定")
            checked += 1
        self.assertGreater(checked, 300, f"只验到 {checked} 条,判据快空过了")


class BuildConfiguration(unittest.TestCase):
    """构建配置里那些会把隐患带回来的设置。"""

    def test_release_profile_does_not_abort_on_panic(self):
        """`panic = "abort"` 会让任何一处 panic 直接 SIGABRT 掉解释器。

        实测过:退出码 134,没有异常、没有回溯、`try/except` 拦不住。
        PyO3 靠展开把 panic 接住转成 `PanicException`,abort 把这条路断了。

        这一条守的是配置而不是行为 —— 要测行为就得在库里留一个故意 panic 的
        入口,那本身就是个不该有的东西。
        """
        text = (REPO / "Cargo.toml").read_text()
        # 只看没被注释掉的行
        live = [ln for ln in text.splitlines() if not ln.lstrip().startswith("#")]
        offending = [ln for ln in live if re.search(r'panic\s*=\s*"abort"', ln)]
        self.assertEqual(
            offending,
            [],
            "Cargo.toml 里又出现了 panic = \"abort\" —— 它会让 Python 进程在任何一处 "
            "Rust panic 上直接 SIGABRT。见该文件里 [profile.release] 上方的说明。",
        )



class RemoveHs(unittest.TestCase):
    """显式氢的合并 —— 化学正确性由 Rust 侧守,这里只验绑定接得对。"""

    def test_merging_drops_the_hydrogen_atoms(self):
        m = omgkit.parse_smiles("[H]OC([H])([H])C")
        n = m.remove_hs()
        self.assertEqual(n, 3, "三个显式氢都该合并")
        self.assertEqual(m.num_atoms, 3, "只剩 C、C、O")

    def test_information_carrying_hydrogens_are_kept(self):
        m = omgkit.parse_smiles("[2H]C")
        self.assertEqual(m.remove_hs(), 0, "氘是另一种核素,不该并成普通氢")

    def test_chirality_survives(self):
        """氢在邻居里的位置一变,手性标记就相对另一个参照系了。"""
        a = omgkit.parse_smiles("[H][C@](N)(O)F")
        a.remove_hs()
        a.sanitize()
        b = omgkit.parse_smiles("N[C@@H](O)F")
        b.sanitize()
        self.assertEqual(
            a.to_canonical_smiles(),
            b.to_canonical_smiles(),
            "合并氢之后成了镜像分子 —— 参照系没换对",
        )


class SanitizePerceivesBondStereo(unittest.TestCase):
    """`sanitize()` 要顺带把方向键换算成双键自己的顺反属性。

    方向是**写法**,顺反是**性质**。方向依附在某根单键上,反应把那根键删掉,
    几何就跟着没了 —— 哪怕双键本身根本没被碰过。产物照样合法、原子数照样对,
    只是顺反悄悄丢了,没有任何东西报错。

    Rust 侧把这一步留给调用方(净化那 12 步调不到上层的对称等价类)。到了
    Python 这边那就成了陷阱,所以绑定层并进 `sanitize()`,由这一档守着。
    """

    # 醚氧承载着 `/`,反应把它换成碘 —— 承载方向的那根键没了
    RXN = "[O;H0;D2:1]-[C;D3:2]>>[OH:1].I-[C;D3:2]"

    def _products(self, smiles):
        rxn = omgkit.parse_reaction(self.RXN)
        mol = omgkit.parse_smiles(smiles)
        mol.sanitize()
        out = []
        for outcome in rxn.run([mol], max_products=20):
            for p in outcome.products:
                p.sanitize()
                out.append(p.to_canonical_smiles())
        return sorted(out)

    def test_geometry_survives_losing_the_bond_that_carried_it(self):
        for smiles in ("CO/C(F)=C(/Cl)C", "CO/C(F)=C(\\Cl)C"):
            got = self._products(smiles)
            self.assertTrue(
                any("/" in p or "\\" in p for p in got),
                f"{smiles} 的产物 {got} 丢了双键几何 —— "
                "sanitize() 没把方向感知成双键自己的属性",
            )

    def test_the_two_geometries_stay_distinguishable(self):
        """上一条的防空过:两边都丢了几何时,它们同样都不含斜杠。"""
        self.assertNotEqual(
            self._products("CO/C(F)=C(/Cl)C"),
            self._products("CO/C(F)=C(\\Cl)C"),
            "一对顺反底物给出了同一批产物,几何被抹平了",
        )


class ByproductsAreOptOutAndLabelled(unittest.TestCase):
    """副产物收口在绑定层的三处约定。

    绑定层只做翻译,所以这里守的不是化学(那由 Rust 侧的
    `crates/omgkit-match/tests/byproduct.rs` 守),而是**翻译本身有没有丢东西**:
    开关默认状态、结论怎么翻成字符串、以及"未决时不给分子"这条约定有没有
    在翻译途中被抹掉。
    """

    ESTER = "[C:1](=[O:2])[OH:3].[OH:4][C:5]>>[C:1](=[O:2])[O:4][C:5]"

    def _run(self, smarts, smis, **kw):
        rxn = omgkit.parse_reaction(smarts)
        mols = []
        for s in smis:
            m = omgkit.parse_smiles(s)
            m.sanitize()
            mols.append(m)
        return rxn.run(mols, **kw)

    def test_default_is_off_and_says_so(self):
        """默认不算 —— 收口要在副本上净化产物,不是零开销。

        而且"没算"要与"算了但没有副产物"分得开:两者都给空列表,只有结论
        字符串区分得了。混成一个的话,一个坏掉的实现看起来会像"这批反应
        本来就没有副产物"。
        """
        (out,) = self._run(self.ESTER, ["CC(=O)O", "CCO"])
        self.assertEqual(out.byproduct_verdict, "off")
        self.assertEqual(out.byproducts, [])
        self.assertEqual(out.byproduct_budget, {})

    def test_water_comes_back_with_a_budget(self):
        (out,) = self._run(self.ESTER, ["CC(=O)O", "CCO"], byproducts=True)
        self.assertEqual(out.byproduct_verdict, "capped")
        self.assertEqual([b.to_canonical_smiles() for b in out.byproducts], ["O"])
        # 账要原样翻过来,不能只翻一个结论
        self.assertEqual(out.byproduct_budget["delta_h"], 2)
        self.assertEqual(out.byproduct_budget["remaining"], 0)

    def test_unresolved_never_ships_a_molecule(self):
        """收不了口就不给分子。

        编一个出来比不给更糟:它拓扑合法、能净化、看不出破绽,只是错的。
        这条约定在 Rust 侧成立,翻译途中也不能被绕过。
        """
        outs = self._run(
            "[N+:1](=[O:2])[O-:3]>>[NH2:1]",
            ["Cc1ccc(cc1)[N+](=O)[O-]"],
            byproducts=True,
        )
        self.assertTrue(outs)
        for o in outs:
            if o.byproduct_verdict.startswith("unresolved"):
                self.assertEqual(o.byproducts, [], o.byproduct_verdict)
        self.assertTrue(
            any(o.byproduct_verdict.startswith("unresolved") for o in outs),
            "这条反应应当有收不了口的 outcome,否则这一档判据是空过的",
        )

    def test_both_entry_points_agree(self):
        """`run` 与 `run_on_substrate` 对分子间底物必须给出同一个副产物。

        后者把输入**拼成一张图**跑,`discarded` 的下标基准因此不同。切不回去的
        后果不是报错,是**静默算错** —— 收口时拿拼接图的下标去索引原始分子,
        越界的被悄悄跳过。这一整个测试类原先一次都没碰过 `run_on_substrate`,
        缺口就是这么留下的。
        """
        smarts, smis = "[OH:3][C:4].[C:1][Cl]>>[C:1][O:3][C:4]", ["OCC", "CCCl"]
        a = self._run(smarts, smis, byproducts=True)[0]
        b = self._run_graph(smarts, smis, byproducts=True)[0]
        self.assertEqual(a.byproduct_verdict, b.byproduct_verdict)
        self.assertEqual(
            [m.to_canonical_smiles() for m in a.byproducts],
            [m.to_canonical_smiles() for m in b.byproducts],
        )
        self.assertEqual([m.to_canonical_smiles() for m in b.byproducts], ["Cl"])
        # 契约:两个入口都按**输入分子**给下标
        self.assertEqual(len(b.discarded), 2)

    def _run_graph(self, smarts, smis, **kw):
        rxn = omgkit.parse_reaction(smarts)
        mols = []
        for s in smis:
            m = omgkit.parse_smiles(s)
            m.sanitize()
            mols.append(m)
        return rxn.run_on_substrate(mols, **kw)

    def test_discarded_is_reported_even_with_the_switch_off(self):
        """`discarded` 是事实,不受开关影响 —— 它不是推断的产物。"""
        (out,) = self._run("[C:1](=[O:2])[OH:3].[N:4]>>[C:1](=[O:2])[N:4]",
                           ["CC(=O)O", "NC"])
        self.assertEqual(len(out.discarded), 2)
        self.assertEqual(len(out.discarded[0]), 1, "羧基那个羟基氧被丢掉了")


class ConformerTests(unittest.TestCase):
    """三维构型的翻译层。化学正确性由 Rust 侧那一套守着,这里盯翻译特有的失效。"""

    def test_conformer_does_not_touch_the_input_mol(self):
        """生成会补显式氢 —— 补在**副本**上,不能改到调用方手里那个分子。"""
        m = omgkit.parse_smiles("CCO")
        before = m.num_atoms
        conf = m.conformer()
        self.assertEqual(m.num_atoms, before, "原分子被改了")
        self.assertGreater(conf.mol.num_atoms, before, "构型里应当带上显式氢")

    def test_coords_line_up_with_the_conformers_own_mol(self):
        """坐标对应的是 `conf.mol` 的原子表,不是原分子的 —— 这是最容易翻错的一处。"""
        conf = omgkit.parse_smiles("CCO").conformer()
        self.assertEqual(len(conf.coords), conf.mol.num_atoms)
        for p in conf.coords:
            self.assertEqual(len(p), 3)
            for v in p:
                self.assertIsInstance(v, float)

    def test_same_molecule_gives_the_same_coordinates(self):
        """全程无随机数 —— 同一个分子每次都给同一组坐标。"""
        smi = "C[C@H](N)C(=O)O"
        a = omgkit.parse_smiles(smi).conformer()
        b = omgkit.parse_smiles(smi).conformer()
        self.assertEqual(a.coords, b.coords)

    def test_chirality_is_reported_and_correct(self):
        """手性中心的账要报出来,而且交付坐标上必须号对。"""
        conf = omgkit.parse_smiles("C[C@H](N)C(=O)O").conformer()
        self.assertEqual(conf.chiral_total, 1)
        self.assertEqual(conf.chiral_ok, conf.chiral_total)

    def test_failure_keeps_the_reason(self):
        """失败要抛 `ValueError`,而且**带上具体原因** —— 翻成一句笼统的话等于把
        排查线索丢掉。"""
        thorium = (
            "CC1=[O+][Th]234([O+]=C(C)C1)([O+]=C(C)CC(=[O+]2)C)"
            "([O+]=C(C)CC(=[O+]3)C)[O+]=C(C)CC(=[O+]4)C"
        )
        with self.assertRaises(ValueError) as cm:
            omgkit.parse_smiles(thorium).conformer()
        self.assertIn("界矩阵", str(cm.exception))

    def test_bonds_and_charges_are_plain_python(self):
        """`bonds` / `formal_charges` 要是普通 list —— `Vec<u8>` 会被 PyO3 特判成
        `bytes`,那种错很安静(索引出来仍是 int,长度也对,类型却错了)。"""
        conf = omgkit.parse_smiles("[NH3+]CC(=O)[O-]").conformer()
        self.assertIsInstance(conf.mol.bonds, list)
        self.assertIsInstance(conf.mol.formal_charges, list)
        self.assertEqual(len(conf.mol.bonds), conf.mol.num_bonds)
        self.assertEqual(len(conf.mol.formal_charges), conf.mol.num_atoms)
        self.assertIn(1, conf.mol.formal_charges)
        self.assertIn(-1, conf.mol.formal_charges)
        # 键级是数值:单键 1.0、芳香 1.5、双键 2.0
        orders = {o for _, _, o in conf.mol.bonds}
        self.assertTrue(orders <= {1.0, 1.5, 2.0, 3.0, 4.0, 0.0}, orders)


if __name__ == "__main__":
    unittest.main(verbosity=2)
