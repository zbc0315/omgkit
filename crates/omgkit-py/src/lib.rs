//! omgkit 的 Python 绑定。
//!
//! # 这一层只做翻译,不做化学
//!
//! 每个函数都应当是"把参数翻过去、把结果翻回来"。一旦在这里写了判断分子的
//! 逻辑,它就只有 Python 用户能碰到,Rust 侧的全套差分测试一概盖不到 ——
//! 那会成为整个项目里唯一没有验证覆盖的一块。
//!
//! # 错误一律翻成异常
//!
//! Rust 侧用 `Result` 表达失败,Python 侧用异常。翻译时**保留原始信息**:
//! 解析错误带着可打印的插字号视图,净化错误带着具体是哪个原子超价。
//! 翻成一句笼统的 "invalid molecule" 等于把排查线索丢掉。
//!
//! # 为什么这个 crate 允许 unsafe
//!
//! 工作区把 `unsafe_code` 设成 `deny`,而 PyO3 的过程宏会生成 unsafe 代码 ——
//! 跨语言边界本来就绕不开它。生成的代码没法逐处标注,只能整个 crate 放开。
//!
//! 放开的范围严格限于本 crate:这里除了宏生成的部分之外**不写**任何 unsafe,
//! 真正的实现全在下面四个纯安全的 crate 里。

#![allow(unsafe_code)]

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use omgkit_core::MolBuilder;

/// 杂化状态的名字。
///
/// 交出去的是**名字**,不是编号:编号只有这一层有,换一版就可能整体错位,
/// 而拿它去建特征列的人看不出来。名字改了则立刻炸在调用方脸上。
fn hybridization_name(h: omgkit_core::Hybridization) -> &'static str {
    use omgkit_core::Hybridization as H;
    match h {
        H::Unspecified => "unspecified",
        H::S => "s",
        H::Sp => "sp",
        H::Sp2 => "sp2",
        H::Sp3 => "sp3",
        H::Sp2d => "sp2d",
        H::Sp3d => "sp3d",
        H::Sp3d2 => "sp3d2",
    }
}

/// 手性标记的名字。记的是**几何类别**,不是 CIP 的 R/S。
fn chiral_tag_name(c: omgkit_core::ChiralTag) -> &'static str {
    use omgkit_core::ChiralTag as C;
    match c {
        C::Unspecified => "unspecified",
        C::Cw => "cw",
        C::Ccw => "ccw",
        C::Allene => "allene",
        C::SquarePlanar => "square_planar",
        C::TrigonalBipyramidal => "trigonal_bipyramidal",
        C::Octahedral => "octahedral",
    }
}

/// 键级的名字。与 `Mol.bonds` 给的数值键级是同一件事的两种写法 ——
/// 芳香在数值那边是 1.5,在这边是它自己一档。
fn bond_order_name(o: omgkit_core::BondOrder) -> &'static str {
    use omgkit_core::BondOrder as B;
    match o {
        B::Unspecified => "unspecified",
        B::Single => "single",
        B::Double => "double",
        B::Triple => "triple",
        B::Quadruple => "quadruple",
        B::Aromatic => "aromatic",
        B::Dative => "dative",
    }
}

/// 双键顺反的名字。`z`/`e` 依 CIP 优先级,`cis`/`trans` 依记录的参照原子。
fn bond_stereo_name(s: omgkit_core::BondStereo) -> &'static str {
    use omgkit_core::BondStereo as S;
    match s {
        S::None => "none",
        S::Z => "z",
        S::E => "e",
        S::Cis => "cis",
        S::Trans => "trans",
    }
}

/// 一个分子。
///
/// 对应 Rust 侧的 `MolBuilder` —— 可变的、逐分子的表示,适合建图与改写。
/// 列式的 `MolBatch` 是另一件事,等批处理接口再暴露。
#[pyclass(name = "Mol", module = "omgkit")]
#[derive(Clone)]
pub struct PyMol {
    inner: MolBuilder,
}

#[pymethods]
impl PyMol {
    /// 原子数。
    #[getter]
    fn num_atoms(&self) -> usize {
        self.inner.num_atoms()
    }

    /// 键数。
    #[getter]
    fn num_bonds(&self) -> usize {
        self.inner.num_bonds()
    }

    /// 逐原子的原子序数,按存储顺序。
    ///
    /// 返回 `list[int]`。**不能**直接返回 `Vec<u8>` —— PyO3 把它特判成
    /// `bytes`,于是 `mol.atomic_nums` 会得到 `b'\x06\x08'` 这种东西:
    /// 索引出来仍是 int,长度也对,但类型错了,而且错得很安静。
    #[getter]
    fn atomic_nums(&self) -> Vec<u16> {
        self.inner
            .atoms()
            .iter()
            .map(|a| u16::from(a.atomic_num))
            .collect()
    }

    /// 跑净化管线,再把双键的方向键换算成双键**自己的**顺反属性。
    ///
    /// **就地修改**,失败时抛 `ValueError`。
    ///
    /// 失败后分子可能已被部分修改 —— 需要"要么全成功要么不动"的调用方
    /// 应当先 `copy()`。这一条与 Rust 侧的语义一致,不在绑定层偷偷加保护:
    /// 悄悄多做一次深拷贝会让批处理的开销凭空翻倍,而调用方无从知情。
    ///
    /// # 为什么这里要多做一步顺反感知
    ///
    /// 净化那 12 步里没有它 —— 感知要用对称等价类,那在净化的**上一层**,
    /// 调不到(理由见 `omgkit_io::stereo` 的模块文档)。Rust 侧因此约定由调用方
    /// 在净化之后自己调一次。
    ///
    /// 可这条约定放到 Python 这边就成了陷阱:只有方向键的分子一旦被
    /// [`PyReaction::run`] 编辑,承载方向的那根单键可能被删掉,双键明明没被碰过,
    /// 几何却跟着没了 —— 产物照样合法、原子数照样对,只是顺反悄悄丢了。
    ///
    /// 方向是**写法**,顺反是**性质**。感知一次之后信息记在双键上,只要两个参照
    /// 原子还在就活得下来。绑定层是给人直接用的,把这一步并进来比留一条要人
    /// 记住的约定稳妥。
    fn sanitize(&mut self) -> PyResult<()> {
        omgkit_chem::sanitize(&mut self.inner).map_err(|e| PyValueError::new_err(e.to_string()))?;
        omgkit_io::stereo::perceive_bond_stereo(&mut self.inner);
        Ok(())
    }

    /// 写成 SMILES,按当前的原子存储顺序。
    fn to_smiles(&self) -> String {
        omgkit_io::smiles::write(&self.inner).smiles
    }

    /// 写成规范 SMILES —— 同一个分子无论原子怎么编号,都得到同一个字符串。
    fn to_canonical_smiles(&self) -> String {
        omgkit_io::canon::canonical_smiles(&self.inner).smiles
    }

    /// 把可以合并的显式氢并进邻居的氢计数,返回删掉的氢原子数。
    ///
    /// **就地修改,而且原子下标会全部改变** —— 删原子必然重排下标。带同位素、
    /// 电荷、映射号、自由基的氢,以及桥氢与承载双键方向的氢都会留着:多留一个
    /// 只是图里多个节点,删错会丢信息。
    fn remove_hs(&mut self) -> usize {
        omgkit_chem::remove_hs(&mut self.inner)
    }

    /// 逐原子的形式电荷,按存储顺序。返回 `list[int]`。
    #[getter]
    fn formal_charges(&self) -> Vec<i16> {
        self.inner
            .atoms()
            .iter()
            .map(|a| i16::from(a.formal_charge))
            .collect()
    }

    /// 逐键的 `(起点, 终点, 键级)`,按存储顺序。返回 `list[tuple[int, int, float]]`。
    ///
    /// 键级是**数值**:单键 1.0、芳香 1.5、双键 2.0、三键 3.0、四重 4.0、
    /// 未指定 0.0(配位键按 1.0 记)。这是 `BondOrder::as_double()` 的值 ——
    /// 不在绑定层另发明一套编号,那种编号只有 Python 这边有,Rust 侧的判据
    /// 一概盖不到。
    #[getter]
    fn bonds(&self) -> Vec<(u32, u32, f64)> {
        self.inner
            .bonds()
            .iter()
            .map(|b| (b.begin, b.end, f64::from(b.order.as_double())))
            .collect()
    }

    /// 生成一个三维构型。
    ///
    /// **不改动本分子**:内部先深拷贝一份,在那一份上净化、感知顺反、补显式氢,
    /// 再生成。所以返回的 `Conformer` 里那个 `mol` 的原子数通常
    /// 比这里多(多出来的是氢),而 `coords` 对应的是**它**的原子表,不是这个。
    ///
    /// 走的是 Rust 侧的 `omgkit_conf::pipeline::conformer_for` —— 那五步的顺序
    /// 与理由都在库里,绑定这一层一步化学都不做。
    ///
    /// 全程无随机数:同一个分子每次都给同一组坐标。
    ///
    /// 净化过不去、界矩阵自相矛盾时抛 `ValueError`。
    fn conformer(&self) -> PyResult<PyConformer> {
        let mut mol = self.inner.clone();
        let conf = omgkit_conf::pipeline::conformer_for(&mut mol)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(PyConformer {
            mol: PyMol { inner: mol },
            conf,
        })
    }

    /// 画一张二维结构式,返回 SVG 源码(字符串)。
    ///
    /// **不改动本分子**:内部先深拷贝一份,在那一份上净化、感知顺反、排布局。
    ///
    /// ```python
    /// open("aspirin.svg", "w").write(omgkit.parse_smiles("CC(=O)Oc1ccccc1C(=O)O").to_svg())
    /// ```
    ///
    /// | 参数 | 默认 | 说明 |
    /// |---|---|---|
    /// | `style` | `"ACS Document 1996"` | 绘图规范。另一套是 `"ChemDraw New Document"`(键长 30 pt,放到网页上更舒展) |
    /// | `aromatic_fill` | `False` | 芳香环底下要不要铺一层径向渐变 |
    /// | `fill_centre` | `"#ffffff"` | 渐变高光那一点的颜色 |
    /// | `fill_edge` | `"#add8e6"` | 渐变外缘的颜色(CSS 的具名色 `lightblue`) |
    ///
    /// 底色铺在**最底层**,不遮任何线条与原子标签。高光落在环心与**最靠左上
    /// 那个顶点**的中点上,与三维图的光源方向是同一个。
    ///
    /// 两个颜色只在 `aromatic_fill=True` 时起作用,但**写错了照样报错** ——
    /// 只在用得上时才校验的话,关着底色调好的一组颜色,开的那天才发现拼错了。
    ///
    /// 规范名或颜色不认识、净化过不去时抛 `ValueError`。
    #[pyo3(signature = (
        style = "ACS Document 1996",
        aromatic_fill = false,
        fill_centre = "#ffffff",
        fill_edge = "#add8e6",
    ))]
    fn to_svg(
        &self,
        style: &str,
        aromatic_fill: bool,
        fill_centre: &str,
        fill_edge: &str,
    ) -> PyResult<String> {
        let fill = omgkit_depict::style::AromaticFill {
            centre: colour("fill_centre", fill_centre)?,
            edge: colour("fill_edge", fill_edge)?,
        };
        let st = omgkit_depict::style::Style {
            aromatic_fill: aromatic_fill.then_some(fill),
            ..style_2d(style)?.clone()
        };
        let mut mol = self.inner.clone();
        omgkit_chem::sanitize(&mut mol).map_err(|e| PyValueError::new_err(e.to_string()))?;
        omgkit_io::stereo::perceive_bond_stereo(&mut mol);
        let d = omgkit_depict::generate(&mol, &st);
        Ok(omgkit_depict::svg::to_svg(
            &omgkit_depict::render::scene(&mol, &d, &st),
            &st,
        ))
    }

    /// 画一张二维结构图,写成 V2000 molblock(`.mol` 文件的内容)。
    ///
    /// **不改动本分子**:内部先深拷贝一份,在那一份上净化、感知顺反、排布局。
    ///
    /// # 立体靠**楔形**,不是坐标
    ///
    /// 二维图的手性写在键块第四列(1 实楔、6 虚楔)。为了把某个中心的构型画
    /// 出来,布局有时要**补一根显式 C–H** —— 楔形恰恰打在那根键上。所以写出去
    /// 的原子数可能比这个分子多,与
    /// `Conformer` 那边同理。
    ///
    /// 画不出构型的中心不会被硬画:那种中心在文件里就是"没写立体",而不是
    /// 随便给一个。
    ///
    /// # 与 `Conformer.to_molblock` 的分工
    ///
    /// 那个写**三维**:立体在坐标本身里,楔形是空的。这个写**二维**:所有 `z`
    /// 都是 0,立体全靠楔形。两种文件都合法,读的人按坐标是不是平的自己分。
    ///
    /// 芳香键会先凯库勒化,理由与三维那条一样。第二行是程序名,不写时间戳。
    ///
    /// 净化过不去、或者分子大到 V2000 装不下(原子或键超过 999)时抛 `ValueError`。
    #[pyo3(signature = (title = ""))]
    fn to_molblock_2d(&self, title: &str) -> PyResult<String> {
        let mut mol = self.inner.clone();
        omgkit_chem::sanitize(&mut mol).map_err(|e| PyValueError::new_err(e.to_string()))?;
        omgkit_io::stereo::perceive_bond_stereo(&mut mol);
        let d = omgkit_depict::generate(&mol, &omgkit_depict::style::Style::ACS_1996);
        // 画出来的那个分子才是要写的:补出来的显式氢也在里面,而楔形就打在
        // 那根 C–H 上。拿原分子写的话,那根键根本不存在,读的人看到的是
        // "没有立体信息"。
        let grown = d.drawn(&mol);
        let orders = omgkit_depict::render::drawn_orders(&grown);
        let coords: Vec<[f64; 3]> = d.coords.iter().map(|p| [p.x, p.y, 0.0]).collect();
        // **作者没写顺反的双键要标成交叉双键。** 布局总得把取代基摆在某一侧,
        // 于是图上每根双键都有一个确定的几何;不标的话,读的一方会把那个几何
        // 当成化学信息读走 —— 凭空多出一句作者没说过的话。
        let unknown = omgkit_io::stereo::unspecified_cis_trans(&grown);
        let rec = omgkit_io::molblock::Record {
            title,
            coords: &coords,
            wedges: &d.wedges,
            orders: &orders,
            unknown_stereo: &unknown,
        };
        omgkit_io::molblock::write_v2000(&grown, &rec)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// 画这个分子时**没能画好的地方**,返回 `dict[str, list]`。
    ///
    /// | 键 | 内容 |
    /// |---|---|
    /// | `degraded` | 布局不得不退化的地方(桥环等),每项一个说明串 |
    /// | `unresolved` | 消冲突之后仍然挤在一起的原子对 |
    /// | `crossings` | 仍然交叉的键对 |
    /// | `unwedged` | 没能画出构型的立体中心。配位几何(`@SP`/`@TB`/`@OH`)这一版画不出来,一律在这里 |
    /// | `misdrawn_stereo` | 画出来的几何与记录的顺反**不符**的双键。八元以上的环里的反式双键会落在这里 —— 环按凸多边形画,环内双键一律画成顺式 |
    ///
    /// 五个都是空的,这张图才把分子完整地表达出来了。
    ///
    /// 下标相对**被画的那个分子** —— 为了承载楔形可能补了显式氢,那时原子数比
    /// `num_atoms` 大;前 `num_atoms` 个与本分子逐项对应。
    ///
    /// 净化过不去时抛 `ValueError`。
    fn depiction_report<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let mut mol = self.inner.clone();
        omgkit_chem::sanitize(&mut mol).map_err(|e| PyValueError::new_err(e.to_string()))?;
        omgkit_io::stereo::perceive_bond_stereo(&mut mol);
        let d = omgkit_depict::generate(&mol, &omgkit_depict::style::Style::ACS_1996);
        let out = PyDict::new(py);
        out.set_item(
            "degraded",
            d.degraded
                .iter()
                .map(|g| format!("{g:?}"))
                .collect::<Vec<_>>(),
        )?;
        out.set_item("unresolved", d.unresolved.clone())?;
        out.set_item("crossings", d.crossings.clone())?;
        out.set_item("unwedged", d.unwedged.clone())?;
        out.set_item("misdrawn_stereo", d.misdrawn_stereo.clone())?;
        Ok(out)
    }

    /// 逐原子的描述符,按存储顺序。返回 `list[dict]`。
    ///
    /// 每个字典 12 个键:
    ///
    /// | 键 | 类型 | 说明 |
    /// |---|---|---|
    /// | `atomic_num` | int | 元素种类。0 是通配原子 `*` |
    /// | `total_degree` | int | 显式邻居数 + 总氢数 |
    /// | `formal_charge` | int | 形式电荷 |
    /// | `chiral_tag` | str | 手性标记的几何类别,不是 R/S |
    /// | `total_num_hs` | int | 显式声明 + 隐式推断,**不含**独立的 `[H]` 原子 |
    /// | `hybridization` | str | 杂化 |
    /// | `is_aromatic` | bool | 是否芳香 |
    /// | `is_in_ring` | bool | 是否在环上 |
    /// | `mass` | float | 标了同位素用该核素的精确质量,否则用标准原子量 |
    /// | `electronegativity` | float \| None | Pauling 电负性 |
    /// | `gasteiger_charge` | float | Gasteiger 部分电荷 |
    /// | `gasteiger_valid` | bool | 上一项算不算得出来 |
    ///
    /// # 交的是描述符,不是编码
    ///
    /// 分类量给的是**名字**(`"sp3"`、`"ccw"`),不是 one-hot,也不是整数编号。
    /// 词表该收哪些元素、留不留"其它"兜底桶,是特征化那一侧的决定 ——
    /// 在这里定死等于把某一个模型的词表焊进库里。
    ///
    /// # 两处"算不出"
    ///
    /// `electronegativity` 为 `None` 表示该元素**没有公认的 Pauling 值**
    /// (稀有气体、Pm/Eu/Tb/Yb/Fr 等);`gasteiger_valid` 为 `False` 表示该原子
    /// 落在 Gasteiger 参数表之外(多数金属),此时 `gasteiger_charge` 是
    /// `nan` 或 `inf`,并且**会沿着图扩散**——同一个分子里的碳也可能因此失效。
    /// 两者都如实交出,不拿 0 顶:那会让"不知道"和"恰好是 0"变成同一格。
    ///
    /// # 前置
    ///
    /// 要先 `sanitize()`。没净化的分子不会报错,只会让芳香、环、杂化、共轭、
    /// 隐式氢数全是解析时的占位值。
    fn atom_descriptors<'py>(&self, py: Python<'py>) -> PyResult<Vec<Bound<'py, PyDict>>> {
        omgkit_chem::atom_descriptors(&self.inner)
            .into_iter()
            .map(|d| {
                let e = PyDict::new(py);
                e.set_item("atomic_num", d.atomic_num)?;
                e.set_item("total_degree", d.total_degree)?;
                e.set_item("formal_charge", d.formal_charge)?;
                e.set_item("chiral_tag", chiral_tag_name(d.chiral_tag))?;
                e.set_item("total_num_hs", d.total_num_hs)?;
                e.set_item("hybridization", hybridization_name(d.hybridization))?;
                e.set_item("is_aromatic", d.is_aromatic)?;
                e.set_item("is_in_ring", d.is_in_ring)?;
                e.set_item("mass", d.mass)?;
                e.set_item("electronegativity", d.electronegativity)?;
                e.set_item("gasteiger_charge", d.gasteiger_charge)?;
                e.set_item("gasteiger_valid", d.gasteiger_is_valid())?;
                Ok(e)
            })
            .collect()
    }

    /// 逐键的描述符,按存储顺序。返回 `list[dict]`。
    ///
    /// 每个字典 7 个键:`begin`、`end`(两端原子下标)、`order`(键级的名字)、
    /// `is_conjugated`、`is_in_ring`、`stereo`(双键顺反)、`stereo_atoms`。
    ///
    /// # 顺反是 `cis`/`trans`,不是 `Z`/`E`
    ///
    /// `Z`/`E` 按 CIP 优先级定义,而 CIP 排序本库没有实现。这里给的是"相对
    /// `stereo_atoms` 那两个参照原子"的顺反。**两项必须一起看**:四取代双键上
    /// 参照挑得不同,同一个几何会得出相反的顺反值。带上参照之后,顺反与 Z/E
    /// 承载的几何信息相同,要 Z/E 的调用方自己排 CIP 换算。
    /// 没有顺反时 `stereo_atoms` 是 `None`。
    ///
    /// # 前置
    ///
    /// 要先 `sanitize()`。顺反尤其要注意:它由净化之后那一步方向键折算填写,
    /// 而 `sanitize()` 已经把两步并在一起了 —— 只跑 Rust 侧的净化不够,那样
    /// 每根双键的 `stereo` 都会是 `"none"`,而且不报错。
    fn bond_descriptors<'py>(&self, py: Python<'py>) -> PyResult<Vec<Bound<'py, PyDict>>> {
        omgkit_chem::bond_descriptors(&self.inner)
            .into_iter()
            .map(|d| {
                let e = PyDict::new(py);
                e.set_item("begin", d.begin)?;
                e.set_item("end", d.end)?;
                e.set_item("order", bond_order_name(d.order))?;
                e.set_item("is_conjugated", d.is_conjugated)?;
                e.set_item("is_in_ring", d.is_in_ring)?;
                e.set_item("stereo", bond_stereo_name(d.stereo))?;
                e.set_item("stereo_atoms", d.stereo_atoms.map(|[a, b]| (a, b)))?;
                Ok(e)
            })
            .collect()
    }

    /// 深拷贝。
    fn copy(&self) -> Self {
        self.clone()
    }

    fn __repr__(&self) -> String {
        format!(
            "<omgkit.Mol {} atoms, {} bonds>",
            self.inner.num_atoms(),
            self.inner.num_bonds()
        )
    }
}

/// 按名字取一套**二维**绘图规范。理由与 [`style_3d`] 那条一样:名字是公开
/// 契约,拼错了要报错,不能静默退回默认的那一套。
///
/// 名字用的就是 `Style::ALL` 里的全名(`"ACS Document 1996"` /
/// `"ChemDraw New Document"`)—— 不另起一套短名:两套拼法并存,迟早有人
/// 在一处加了新规范、忘了在另一处加短名。
fn style_2d(name: &str) -> PyResult<&'static omgkit_depict::style::Style> {
    omgkit_depict::style::Style::ALL
        .iter()
        .find(|s| s.name == name)
        .ok_or_else(|| {
            let all: Vec<&str> = omgkit_depict::style::Style::ALL
                .iter()
                .map(|s| s.name)
                .collect();
            PyValueError::new_err(format!(
                "不认识的绘图规范 {name:?};认识的是:{}",
                all.join("、")
            ))
        })
}

/// 读一个 `#rrggbb`,读不出来就报错并说清是哪个参数。
fn colour(what: &str, hex: &str) -> PyResult<[u8; 3]> {
    omgkit_depict::palette::parse_hex(hex)
        .ok_or_else(|| PyValueError::new_err(format!("{what} 不是 #rrggbb 形式的颜色:{hex:?}")))
}

/// 按名字取一套三维样式。**名字是公开契约的一部分**,不认识就报错并把
/// 认识的都列出来 —— 静默退回默认样式的话,拼错一个字母就会得到一张
/// "看着对但不是你要的那个样式"的图。
fn style_3d(name: &str) -> PyResult<&'static omgkit_depict::three::Style3D> {
    omgkit_depict::three::Style3D::ALL
        .iter()
        .find(|s| s.name == name)
        .ok_or_else(|| {
            let all: Vec<&str> = omgkit_depict::three::Style3D::ALL
                .iter()
                .map(|s| s.name)
                .collect();
            PyValueError::new_err(format!(
                "不认识的三维样式 {name:?};认识的是:{}",
                all.join("、")
            ))
        })
}

/// 一个三维构型。
///
/// 由 `Mol.conformer()` 产出。里面既有坐标,也有**坐标对应的
/// 那个分子** —— 生成时补了显式氢,原子表与原分子不是同一份。
#[pyclass(name = "Conformer", module = "omgkit")]
pub struct PyConformer {
    mol: PyMol,
    conf: omgkit_conf::pipeline::Conformer,
}

#[pymethods]
impl PyConformer {
    /// 坐标对应的那个分子(补过显式氢的那一份)。
    #[getter]
    fn mol(&self) -> PyMol {
        self.mol.clone()
    }

    /// 逐原子的 `(x, y, z)`,单位 Å,顺序与 [`mol`](Self::mol) 的原子表一致。
    #[getter]
    fn coords(&self) -> Vec<(f64, f64, f64)> {
        self.conf
            .coords
            .iter()
            .map(|p| (p[0], p[1], p[2]))
            .collect()
    }

    /// 精修之后的误差函数值。0 表示所有距离都落进了界内。
    #[getter]
    fn energy(&self) -> f64 {
        self.conf.energy
    }

    /// 精修**之前**的误差函数值 —— 与 [`energy`](Self::energy) 一起看才知道精修干了多少活。
    #[getter]
    fn energy_before(&self) -> f64 {
        self.conf.energy_before
    }

    /// 精修有没有收敛(梯度降到阈值以下)。
    #[getter]
    fn converged(&self) -> bool {
        self.conf.converged
    }

    /// 精修迭代了多少次。
    #[getter]
    fn iterations(&self) -> usize {
        self.conf.iterations
    }

    /// 手性中心总数。
    #[getter]
    fn chiral_total(&self) -> usize {
        self.conf.chiral_total
    }

    /// 其中在交付坐标上号正确的个数。**应当等于
    /// [`chiral_total`](Self::chiral_total)** —— 不等就是把某个中心摆成了对映体。
    #[getter]
    fn chiral_ok(&self) -> usize {
        self.conf.chiral_ok
    }

    /// 画成**三维分子图**,返回一段 SVG。
    ///
    /// `style` 取四套之一 —— 名字与半径都取自 Jmol 自己文档里的
    /// standard rendering styles:
    ///
    /// | `style` | 球半径 | 键(圆柱)半径 | 看什么 |
    /// |---|---|---|---|
    /// | `"space-filling"` | 100% 范德华半径 | 不画 | 分子占多大地方 |
    /// | `"ball-and-stick"`(默认) | 23% vdW | 0.15 Å | 键长键角、构型 |
    /// | `"stick"` | 与键同粗 | 0.30 Å | 骨架走向 |
    /// | `"wireframe"` | 不画 | 0.01 Å | 大体系、快速预览 |
    ///
    /// 按元素上 CPK(Jmol)色,**键的两半各随自己那一端的颜色**。
    ///
    /// ```python
    /// conf = omgkit.parse_smiles("CC(=O)Oc1ccccc1C(=O)O").conformer()
    /// open("aspirin.svg", "w").write(conf.to_svg())
    /// ```
    ///
    /// 画的是 [`mol`](Self::mol) 那一份(补过显式氢的)—— 三维图里氢是看得见的
    /// 实体,不画的话读图的人看到的是另一个分子。
    ///
    /// `style` 不认识时抛 `ValueError`,并把认识的四个列出来。
    #[pyo3(signature = (style = "ball-and-stick"))]
    fn to_svg(&self, style: &str) -> PyResult<String> {
        let st = style_3d(style)?;
        let d = omgkit_depict::three::depict(&self.mol.inner, &self.conf.coords, st)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(omgkit_depict::svg::to_svg(
            &d.scene,
            &omgkit_depict::style::Style::ACS_1996,
        ))
    }

    /// 三维图的诊断:视角定不定得下来,以及每个原子落在画布哪里。
    ///
    /// 返回的字典:
    ///
    /// | 键 | 内容 |
    /// |---|---|
    /// | `style` | 用的哪套样式 |
    /// | `width`、`height` | 画布尺寸(磅) |
    /// | `degenerate` | **主轴不唯一**。对称性强制两个主惯量相等时为真(甲烷、四氯化碳、氨、乙炔)。图不是错的,但它的取向没有承载任何信息 —— 别照着它比两个分子的姿态 |
    /// | `atoms` | 每个原子一项:`x`、`y`(画布坐标,磅)、`radius`(球半径,磅,不画球的样式是 0)、`depth`(深度,Å,越大越靠前) |
    ///
    /// 想在图上加标注就用 `atoms` —— SVG 里的圆没有原子号,从图形反推是猜。
    #[pyo3(signature = (style = "ball-and-stick"))]
    fn depiction_3d_report<'py>(
        &self,
        py: Python<'py>,
        style: &str,
    ) -> PyResult<Bound<'py, PyDict>> {
        let st = style_3d(style)?;
        let d = omgkit_depict::three::depict(&self.mol.inner, &self.conf.coords, st)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let out = PyDict::new(py);
        out.set_item("style", d.style_name)?;
        out.set_item("width", d.scene.width)?;
        out.set_item("height", d.scene.height)?;
        out.set_item("degenerate", d.view.degenerate)?;
        let atoms = pyo3::types::PyList::empty(py);
        for p in &d.placed {
            let one = PyDict::new(py);
            one.set_item("x", p.at.x)?;
            one.set_item("y", p.at.y)?;
            one.set_item("radius", p.radius)?;
            one.set_item("depth", p.depth)?;
            atoms.append(one)?;
        }
        out.set_item("atoms", atoms)?;
        Ok(out)
    }

    /// 写成 V2000 molblock(`.mol` 文件的内容),末尾带 `M  END`。
    ///
    /// 拼 `.sdf` 时每条后面接数据字段和 `$$$$`:
    ///
    /// ```python
    /// with open("out.sdf", "w") as f:
    ///     for smi in smiles_list:
    ///         conf = omgkit.parse_smiles(smi).conformer()
    ///         f.write(conf.to_molblock(title=smi))
    ///         f.write("$$$$\n")
    /// ```
    ///
    /// **芳香键会先凯库勒化** —— molblock 里没有"芳香键"这回事,留着它写出去
    /// 要么歧义、要么被读成饱和环。凯库勒化失败(比如芳香体系里有通配原子)时
    /// 抛 `ValueError`,不写一个读回来是另一个分子的文件。
    ///
    /// 第二行写的是程序名,**不写时间戳** —— 同一个分子每次写出都逐字节相同。
    #[pyo3(signature = (title = ""))]
    fn to_molblock(&self, title: &str) -> PyResult<String> {
        let mut kek = self.mol.inner.clone();
        omgkit_chem::kekulize(&mut kek)
            .map_err(|e| PyValueError::new_err(format!("凯库勒化失败:{e}")))?;
        let orders: Vec<_> = kek.bonds().iter().map(|b| b.order).collect();
        // 三维同理,而且一模一样地危险:嵌出来的构象给每根双键一个确定的二面角,
        // 作者没写顺反的那些键不标交叉,读回来就成了"作者说是顺式"。
        let unknown = omgkit_io::stereo::unspecified_cis_trans(&self.mol.inner);
        let rec = omgkit_io::molblock::Record {
            title,
            coords: &self.conf.coords,
            wedges: &[],
            orders: &orders,
            unknown_stereo: &unknown,
        };
        omgkit_io::molblock::write_v2000(&self.mol.inner, &rec)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn __repr__(&self) -> String {
        format!(
            "<omgkit.Conformer atoms={} energy={:.3e} converged={} chiral={}/{}>",
            self.conf.coords.len(),
            self.conf.energy,
            // Rust 的 bool 打出来是 true/false;这是给 Python 看的 repr
            if self.conf.converged { "True" } else { "False" },
            self.conf.chiral_ok,
            self.conf.chiral_total
        )
    }
}

/// 一个从 `.mol` / `.sdf` 文件读出来的记录:分子,加上它在文件里的坐标。
///
/// 由 `parse_molblock` 产出。分开成一个类而不是直接给 `Mol`,是因为坐标
/// **不在** `Mol` 里:molblock 的立体化学一半靠坐标表达,把坐标丢掉等于把这一半
/// 丢掉,而丢的时候一声不响。
#[pyclass(name = "Molblock", module = "omgkit")]
#[derive(Clone)]
pub struct PyMolblock {
    mol: PyMol,
    coords: Vec<[f64; 3]>,
    title: String,
    is_3d: bool,
}

#[pymethods]
impl PyMolblock {
    /// 分子。**已经净化过,立体也打上了** —— 二维靠楔形与平面投影,
    /// 三维靠有符号体积与二面角,两条路都接上了。
    #[getter]
    fn mol(&self) -> PyMol {
        self.mol.clone()
    }

    /// 逐原子的 `(x, y, z)`,顺序与 [`mol`](Self::mol) 的原子表一致。
    ///
    /// 二维图的 `z` 一律是 0。
    #[getter]
    fn coords(&self) -> Vec<(f64, f64, f64)> {
        self.coords.iter().map(|p| (p[0], p[1], p[2])).collect()
    }

    /// 文件第一行的标题。
    #[getter]
    fn title(&self) -> &str {
        &self.title
    }

    /// 坐标是不是三维的(有任何一个 `z` 不为 0)。
    ///
    /// 文件里没有哪个字段直说这件事,只能这么判 —— 与外部实现同法。
    #[getter]
    fn is_3d(&self) -> bool {
        self.is_3d
    }

    fn __repr__(&self) -> String {
        format!(
            "<omgkit.Molblock {} atoms {}D{}>",
            self.mol.inner.num_atoms(),
            if self.is_3d { 3 } else { 2 },
            if self.title.is_empty() {
                String::new()
            } else {
                format!(" {:?}", self.title)
            }
        )
    }
}

/// 读完之后共用的那三步:净化,再从坐标与楔形回来打立体标记。
///
/// 单条(`parse_molblock`)与整份 SDF(`read_sdf`)都走它。**只有这一处** ——
/// 两处各写一遍的话,迟早一边打了立体、另一边没打,而那种差别不报错。
fn finish(got: omgkit_io::molblock::Molblock) -> Result<PyMolblock, String> {
    let mut mol = got.mol;
    omgkit_chem::sanitize(&mut mol).map_err(|e| e.to_string())?;
    // 净化**之后**才打得上:手性要知道中心有几个隐式氢,顺反要用对称等价类。
    // 两个函数自己认出三维就整个不做,绑定这一层不加判断 —— 加了的话,
    // "什么时候读得出立体"就有了两个住处。
    let _ = omgkit_io::stereo::assign_chirality_2d(&mut mol, &got.coords, &got.wedges);
    let _ = omgkit_io::stereo::assign_bond_stereo_2d(&mut mol, &got.coords, &got.unknown_stereo);
    // 三维那两个。四个 `assign_*` 各自认出维数,不合的那一对返回 0 ——
    // 由它们自己判、而不是在这里 `if is_3d`,是为了让"什么时候读得出立体"
    // 只有一个住处。
    let _ = omgkit_io::stereo::assign_chirality_3d(&mut mol, &got.coords);
    let _ = omgkit_io::stereo::assign_bond_stereo_3d(&mut mol, &got.coords, &got.unknown_stereo);
    Ok(PyMolblock {
        mol: PyMol { inner: mol },
        coords: got.coords,
        title: got.title,
        is_3d: got.is_3d,
    })
}

/// 读一条 V2000 molblock(`.mol` 文件的内容,或 `.sdf` 里 `$$$$` 之前的一段)。
///
/// 读不出来时抛 `ValueError`,消息里说明是哪一行的什么字段 —— V3000 会被明确
/// 拒收,而不是当成 V2000 硬读出一个错分子。
///
/// # 它替你多做了两步,而且必须多做
///
/// 别的解析函数(如 `parse_smiles`)交回来的是**没净化**的分子,由调用方自己
/// 决定什么时候净化。这里不一样:
///
/// * SMILES 的立体写在串里,净化推迟不丢任何东西;
/// * molblock 的立体一半在**坐标与楔形**里,而那两样在 `Mol` 之外。给它们打上
///   标记要先知道每个原子有几个隐式氢、要用对称等价类 —— 两样都是净化算出来的。
///
/// 所以顺序只能是"读 → 净化 → 回来打立体标记",而中间那一步一旦交给调用方,
/// 漏了不会报错,只会**静默地把整个文件的立体丢掉**。与 `Mol.sanitize`
/// 把顺反感知并进来是同一个理由:绑定层是给人直接用的,把必须成对的两步拆开
/// 就是个陷阱。
///
/// # 二维三维都读立体
///
/// 二维靠楔形定手性、靠平面投影定顺反;三维靠有符号体积定手性、靠二面角定顺反。
/// 走哪条由坐标自己说了算(有任何一个 `z` 不为零就是三维),调用方不必分。
/// 读出来之后的那两步:净化,然后回来打立体标记。**只有这一处**。
///
/// 单条(`parse_molblock`)与整份 SDF(`read_sdf`)都走它。两处各写一遍的话,
/// 迟早一边打了立体、另一边没打,而那种差别不报错。
#[pyfunction]
fn parse_molblock(text: &str) -> PyResult<PyMolblock> {
    let got =
        omgkit_io::molblock::read_v2000(text).map_err(|e| PyValueError::new_err(e.to_string()))?;
    finish(got).map_err(PyValueError::new_err)
}

/// SDF 里的一条记录。
///
/// **读不了的那条也在这里**,`error` 不是 `None`,`block` 是 `None` ——
/// 见 `read_sdf` 的文档。
#[pyclass(name = "SdfRecord", module = "omgkit")]
pub struct PySdfRecord {
    block: Option<PyMolblock>,
    data: Vec<(String, String)>,
    error: Option<String>,
}

#[pymethods]
impl PySdfRecord {
    /// 分子那一段。这条读不了时是 `None`。
    #[getter]
    fn block(&self) -> Option<PyMolblock> {
        self.block.clone()
    }

    /// 数据字段,按文件里出现的顺序,`list[tuple[str, str]]`。
    ///
    /// **不是字典**:同名字段在真实文件里出现过(供应商把多次测量各写一行),
    /// 换成字典会静默地只留最后一条。名字重不重是调用方的判断。
    ///
    /// 这条读不了时是空的。
    #[getter]
    fn data(&self) -> Vec<(String, String)> {
        self.data.clone()
    }

    /// 这条读不了的原因;读得了时是 `None`。
    #[getter]
    fn error(&self) -> Option<String> {
        self.error.clone()
    }

    fn __repr__(&self) -> String {
        match (&self.error, &self.block) {
            (Some(e), _) => format!("<omgkit.SdfRecord 读不了:{e}>"),
            (None, Some(b)) => format!(
                "<omgkit.SdfRecord {} atoms, {} fields>",
                b.mol.inner.num_atoms(),
                self.data.len()
            ),
            (None, None) => "<omgkit.SdfRecord 空>".to_string(),
        }
    }
}

/// 逐条读一个 SDF(`.sdf` 文件的内容),返回 `list[SdfRecord]`。
///
/// # 读不了的那条**不抛异常,也不消失**
///
/// 抛异常会停在坏记录上,后面几千条一起丢掉;静默跳过会让**分母悄悄变小** ——
/// 调用方数出来的条数与文件里的不符,而没有任何地方报错。两种都不行。
///
/// 所以每条都在返回的列表里占一个位置:读不了的那条 `error` 是一句话、
/// `block` 是 `None`,后面的照读不误。怎么处理由调用方决定:
///
/// ```python
/// for i, rec in enumerate(omgkit.read_sdf(text)):
///     if rec.error:
///         print(f"第 {i} 条读不了:{rec.error}")
///         continue
///     print(rec.block.mol.to_canonical_smiles(), rec.data)
/// ```
///
/// 真实语料里这一档是有的:金属茂类配合物的键数超出 V2000 的表达能力,写出方
/// 自己就换成了 V3000,而 V3000 这里明确拒收。
///
/// # 整份读进内存
///
/// 文本本来就整份在内存里(参数就是个 `str`),这里再把每条都解析出来。
/// 超大文件(几十万条)的峰值内存要按这个估。
///
/// # 立体与 `parse_molblock` 同一条路
///
/// 每条都是"读 → 净化 → 回来打立体标记",与单条那个函数共用同一段代码。
/// 三维文件的立体同样读得出来(有符号体积定手性、二面角定顺反),
/// 走哪条由坐标自己说了算 —— 理由见 `parse_molblock`。
#[pyfunction]
fn read_sdf(text: &str) -> Vec<PySdfRecord> {
    omgkit_io::molblock::read_sdf(text)
        .map(|rec| match rec {
            Err(e) => PySdfRecord {
                block: None,
                data: Vec::new(),
                error: Some(e.to_string()),
            },
            Ok(rec) => match finish(rec.block) {
                Err(e) => PySdfRecord {
                    block: None,
                    data: rec.data,
                    error: Some(e),
                },
                Ok(block) => PySdfRecord {
                    block: Some(block),
                    data: rec.data,
                    error: None,
                },
            },
        })
        .collect()
}

/// 解析 SMILES,返回一个 `Mol`。
///
/// 解析失败时抛 `ValueError`,消息里带插字号视图指出出错在第几个字符。
///
/// **只解析,不净化。** 芳香标志、环信息、隐式氢数、杂化都还是空的 ——
/// 要用它们(或者要写规范 SMILES、要画图、要描述符)先调 `Mol.sanitize()`。
#[pyfunction]
fn parse_smiles(smiles: &str) -> PyResult<PyMol> {
    match omgkit_io::smiles::parse(smiles) {
        Ok(inner) => Ok(PyMol { inner }),
        // `render()` 给的是多行的插字号视图,比一句 "parse error" 有用得多
        Err(e) => Err(PyValueError::new_err(e.render())),
    }
}

/// 一个 SMARTS 查询。
#[pyclass(name = "Query", module = "omgkit")]
pub struct PyQuery {
    inner: omgkit_io::smarts::QueryMol,
}

#[pymethods]
impl PyQuery {
    /// 查询里的原子数。
    #[getter]
    fn num_atoms(&self) -> usize {
        self.inner.num_atoms()
    }

    /// 找出这个查询在分子里的全部匹配。
    ///
    /// 返回 `list[list[int]]`,每个内层列表按**查询原子顺序**给出对应的分子原子下标。
    ///
    /// 每次调用都重新算一遍分子的查询性质(环成员数、最小环大小等)。**不缓存**:
    /// `sanitize()` 之类的操作会就地改分子,缓存一旦失效就会静默给出错答案,
    /// 而那种错比多算一遍贵得多。要在同一个分子上匹配很多模式时,这一层
    /// 目前还没有复用入口。
    // Python 侧叫 `match`(软关键字,当方法名完全合法,`re.Pattern.match` 就是);
    // Rust 侧不能用这个名字,它是模式匹配的关键字
    #[pyo3(name = "match")]
    #[pyo3(signature = (mol, *, uniquify = true, max_matches = 0, use_chirality = true))]
    fn match_(
        &self,
        mol: &PyMol,
        uniquify: bool,
        max_matches: usize,
        use_chirality: bool,
    ) -> Vec<Vec<u32>> {
        let props = omgkit_match::MolProps::compute(&mol.inner);
        let opts = omgkit_match::MatchOptions {
            max_matches,
            uniquify,
            use_chirality,
        };
        omgkit_match::substructure_matches(&self.inner, &mol.inner, &props, opts)
    }

    fn __repr__(&self) -> String {
        format!("<omgkit.Query {} atoms>", self.inner.num_atoms())
    }
}

/// 解析 SMARTS。
#[pyfunction]
fn parse_smarts(smarts: &str) -> PyResult<PyQuery> {
    match omgkit_io::smarts::parse(smarts) {
        Ok(inner) => Ok(PyQuery { inner }),
        Err(e) => Err(PyValueError::new_err(e.render())),
    }
}

/// 一组产物,连同(可选的)带映射号的反应物副本。
#[pyclass(name = "Outcome", module = "omgkit")]
pub struct PyOutcome {
    /// 产物,每个产物模板一个分子。
    #[pyo3(get)]
    products: Vec<PyMol>,
    /// 带原子映射号的反应物副本。只在 `run(..., atom_mapping=True)` 时非空。
    #[pyo3(get)]
    reactants: Vec<PyMol>,
    /// 收口出来的副产物。只在 `run(..., byproducts=True)` 且账闭合时非空。
    #[pyo3(get)]
    byproducts: Vec<PyMol>,
    /// 收口的结论:`"nothing"` / `"capped"` / `"bonded(n)"` /
    /// `"unresolved(原因)"`;没开 `byproducts` 时为 `"off"`。
    ///
    /// **`"unresolved(...)"` 时 `byproducts` 必然为空** —— 收不了口就不给分子,
    /// 编一个出来比不给更糟:它拓扑合法、能净化、看不出破绽,只是错的。
    #[pyo3(get)]
    byproduct_verdict: String,
    /// 收口的原子账,键为 `open_valence` / `fragment_hydrogens` / `delta_h` /
    /// `need` / `remaining` / `delta_charge` / `fragment_charge` /
    /// `charge_shift`。用来自己复核结论。
    #[pyo3(get)]
    byproduct_budget: std::collections::BTreeMap<String, i64>,
    /// `discarded[i]` = 第 i 个输入分子里没有进入任何产物的原子下标。
    ///
    /// 这是**事实**,与收不收得了口无关 —— 收口失败时它照样有值,而那正是最
    /// 需要它的时候。
    #[pyo3(get)]
    discarded: Vec<Vec<u32>>,
}

#[pymethods]
impl PyOutcome {
    fn __repr__(&self) -> String {
        format!(
            "<omgkit.Outcome {} products, {} mapped reactants, byproducts {}>",
            self.products.len(),
            self.reactants.len(),
            self.byproduct_verdict
        )
    }
}

/// `Vec<Outcome>` → `Vec<PyOutcome>`。`run` 与 `run_on_substrate` 共用,免得
/// 两处翻译各写一遍、然后慢慢长歪。
fn translate_outcomes(
    outcomes: Vec<omgkit_match::Outcome>,
    inputs: &[MolBuilder],
    byproducts: bool,
) -> Vec<PyOutcome> {
    outcomes
        .into_iter()
        .map(|o| {
            let (by_mols, verdict, budget) = if byproducts {
                let by = omgkit_match::byproduct::reconstruct(inputs, &o);
                translate_byproducts(&by)
            } else {
                (Vec::new(), "off".to_string(), Default::default())
            };
            PyOutcome {
                products: o
                    .products
                    .into_iter()
                    .map(|inner| PyMol { inner })
                    .collect(),
                reactants: o
                    .reactants
                    .into_iter()
                    .map(|inner| PyMol { inner })
                    .collect(),
                byproducts: by_mols,
                byproduct_verdict: verdict,
                byproduct_budget: budget,
                discarded: o.discarded,
            }
        })
        .collect()
}

/// 把一次收口的结论翻成 Python 侧的两样东西:一句话与一张账。
fn translate_byproducts(
    by: &omgkit_match::byproduct::Byproducts,
) -> (Vec<PyMol>, String, std::collections::BTreeMap<String, i64>) {
    use omgkit_match::byproduct::{Unresolved, Verdict};
    let verdict = match by.verdict {
        Verdict::Nothing => "nothing".to_string(),
        Verdict::Capped => "capped".to_string(),
        Verdict::Bonded { bonds } => format!("bonded({bonds})"),
        Verdict::Unresolved(why) => {
            let name = match why {
                Unresolved::OddValence => "odd_valence",
                Unresolved::BudgetExceedsValence => "budget_exceeds_valence",
                Unresolved::HydrogenBudgetNegative => "hydrogen_budget_negative",
                Unresolved::TooManyBonds => "too_many_bonds",
                Unresolved::ProductsUnsanitizable => "products_unsanitizable",
                Unresolved::NoPairing => "no_pairing",
                Unresolved::FragmentUnsanitizable => "fragment_unsanitizable",
                Unresolved::SubstrateUnkekulizable => "substrate_unkekulizable",
                Unresolved::StrainedClosure => "strained_closure",
                Unresolved::BudgetMismatch => "budget_mismatch",
            };
            format!("unresolved({name})")
        }
    };
    let b = by.budget;
    let budget = [
        ("open_valence", i64::from(b.open_valence)),
        ("fragment_hydrogens", i64::from(b.fragment_hydrogens)),
        ("delta_h", i64::from(b.delta_h)),
        ("need", i64::from(b.need)),
        ("remaining", i64::from(b.remaining)),
        ("delta_charge", i64::from(b.delta_charge)),
        ("fragment_charge", i64::from(b.fragment_charge)),
        ("charge_shift", i64::from(b.charge_shift)),
    ]
    .into_iter()
    .map(|(k, v)| (k.to_string(), v))
    .collect();
    let mols = by
        .molecules
        .iter()
        .map(|m| PyMol { inner: m.clone() })
        .collect();
    (mols, verdict, budget)
}

/// 一条反应模板。
#[pyclass(name = "Reaction", module = "omgkit")]
pub struct PyReaction {
    inner: omgkit_io::smarts::Reaction,
}

#[pymethods]
impl PyReaction {
    /// 反应物模板的个数。调 `run` 时要给同样多的分子。
    #[getter]
    fn num_reactant_templates(&self) -> usize {
        self.inner.reactants.len()
    }

    /// 产物模板的个数。
    #[getter]
    fn num_product_templates(&self) -> usize {
        self.inner.products.len()
    }

    /// 对一组反应物跑这条反应。
    ///
    /// 每个反应物模板配一个**互不相同**的输入分子;分子数与模板数不等时返回
    /// 空列表 —— 那一档是 `run_on_substrate` 的形状(多个片段落在同一个分子上)。
    ///
    /// # 递入顺序不影响出不出产物
    ///
    /// **位置不是化学。** 谁先谁后是你敲键盘的顺序,不该决定这条反应跑不跑得
    /// 起来。所以顺序对不上时不会交白卷:引擎先试"第 i 个配第 i 个",给不出
    /// 产物才去找别的一一对应。顺序本来就对得上时,开销与只试那一种完全相同。
    ///
    /// 于是**返回空只剩一个意思:这批分子上没有反应位点**。
    ///
    /// 这一条是量出来的:USPTO-50k 正向语料按记录自带的分子顺序直接调用,
    /// 约 689 条交白卷;抽样 4000 条逐条核过,其中 59 条**全部**只是顺序对不上,
    /// 没有一条是真匹配不上 —— 而调用方拿到的是同一个空列表。
    ///
    /// `atom_mapping` 为真时,每个结果的 `reactants` 填上带映射号的反应物副本,
    /// 产物侧对应原子打同一个号 —— 两侧合起来就是一条完整的原子映射反应。
    /// `byproducts` 为真时,把模板丢弃的原子收口成分子填进
    /// `Outcome.byproducts`。默认关闭:收口要在副本上
    /// 净化产物才算得出氢预算,不是零开销。
    #[pyo3(signature = (reactants, *, max_products = 0, atom_mapping = false, byproducts = false))]
    fn run(
        &self,
        reactants: Vec<PyMol>,
        max_products: usize,
        atom_mapping: bool,
        byproducts: bool,
    ) -> Vec<PyOutcome> {
        let originals: Vec<MolBuilder> = reactants.iter().map(|m| m.inner.clone()).collect();
        let inputs: Vec<(MolBuilder, omgkit_match::MolProps)> = reactants
            .into_iter()
            .map(|m| {
                let props = omgkit_match::MolProps::compute(&m.inner);
                (m.inner, props)
            })
            .collect();
        let outs = omgkit_match::run_reactants(&self.inner, &inputs, max_products, atom_mapping);
        translate_outcomes(outs, &originals, byproducts)
    }

    /// 把整个反应物侧当作**一张图**上的查询来跑,而不是按位置配对。
    ///
    /// `run` 要求"第 i 个分子配第 i 个模板片段",于是模板片段数
    /// 比分子数多时直接返回空 —— 而那正是**分子内反应**的形状:两个片段落在
    /// 同一个分子上。本方法把输入拼成一张图,让每个片段在整张图上自由找位置,
    /// 只要求各片段匹配到的原子两两不重叠。
    ///
    /// - 分子间:与 `run` 结果一致,且不必再枚举输入的排列
    /// - 分子内:`run` 表达不了的那一档
    /// - 盐:阳离子与阴离子是同一个分子的两个组分,模板可以同时碰到
    ///
    /// 代价是搜索空间变大,耗时不如 `run` 可预测;要稳定耗时就用 `run`。
    #[pyo3(signature = (substrate, *, max_products = 0, atom_mapping = false, byproducts = false))]
    fn run_on_substrate(
        &self,
        substrate: Vec<PyMol>,
        max_products: usize,
        atom_mapping: bool,
        byproducts: bool,
    ) -> Vec<PyOutcome> {
        let originals: Vec<MolBuilder> = substrate.iter().map(|m| m.inner.clone()).collect();
        let inputs: Vec<(MolBuilder, omgkit_match::MolProps)> = substrate
            .into_iter()
            .map(|m| {
                let props = omgkit_match::MolProps::compute(&m.inner);
                (m.inner, props)
            })
            .collect();
        let outs = omgkit_match::run_on_substrate(&self.inner, &inputs, max_products, atom_mapping);
        translate_outcomes(outs, &originals, byproducts)
    }

    fn __repr__(&self) -> String {
        format!(
            "<omgkit.Reaction {} -> {}>",
            self.inner.reactants.len(),
            self.inner.products.len()
        )
    }
}

/// 解析反应 SMARTS(`反应物>试剂>产物` 的三段式)。
#[pyfunction]
fn parse_reaction(smarts: &str) -> PyResult<PyReaction> {
    match omgkit_io::smarts::parse_reaction(smarts) {
        Ok(inner) => Ok(PyReaction { inner }),
        Err(e) => Err(PyValueError::new_err(e.render())),
    }
}

#[pymodule]
fn omgkit(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_class::<PyMol>()?;
    m.add_class::<PyConformer>()?;
    m.add_class::<PyQuery>()?;
    m.add_class::<PyReaction>()?;
    m.add_class::<PyOutcome>()?;
    m.add_class::<PyMolblock>()?;
    m.add_class::<PySdfRecord>()?;
    m.add_function(wrap_pyfunction!(parse_smiles, m)?)?;
    m.add_function(wrap_pyfunction!(parse_smarts, m)?)?;
    m.add_function(wrap_pyfunction!(parse_reaction, m)?)?;
    m.add_function(wrap_pyfunction!(parse_molblock, m)?)?;
    m.add_function(wrap_pyfunction!(read_sdf, m)?)?;
    Ok(())
}
