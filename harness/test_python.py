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


if __name__ == "__main__":
    unittest.main(verbosity=2)
