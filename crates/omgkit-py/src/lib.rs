//! omgkit 的 Python 绑定。
//!
//! # 这一层只做翻译,不做化学
//!
//! 每个函数都应当是"把参数翻过去、把结果翻回来"。一旦在这里写了判断分子的
//! 逻辑,它就只有 Python 用户能碰到,Rust 侧的全套差分测试一概盖不到 ——
//! 那是本项目最不想要的那种代码。
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

use omgkit_core::MolBuilder;

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

/// 解析 SMILES。解析失败时抛 `ValueError`,消息里带插字号视图指出出错位置。
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
}

#[pymethods]
impl PyOutcome {
    fn __repr__(&self) -> String {
        format!(
            "<omgkit.Outcome {} products, {} mapped reactants>",
            self.products.len(),
            self.reactants.len()
        )
    }
}

/// 一条反应模板。
#[pyclass(name = "Reaction", module = "omgkit")]
pub struct PyReaction {
    inner: omgkit_io::smarts::Reaction,
}

#[pymethods]
impl PyReaction {
    /// 反应物模板的个数。调 [`run`](Self::run) 时要给同样多的分子。
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
    /// `reactants[i]` 要与第 i 个反应物模板对应,数目不符时返回空列表。
    ///
    /// `atom_mapping` 为真时,每个结果的 `reactants` 填上带映射号的反应物副本,
    /// 产物侧对应原子打同一个号 —— 两侧合起来就是一条完整的原子映射反应。
    #[pyo3(signature = (reactants, *, max_products = 0, atom_mapping = false))]
    fn run(
        &self,
        reactants: Vec<PyMol>,
        max_products: usize,
        atom_mapping: bool,
    ) -> Vec<PyOutcome> {
        let inputs: Vec<(MolBuilder, omgkit_match::MolProps)> = reactants
            .into_iter()
            .map(|m| {
                let props = omgkit_match::MolProps::compute(&m.inner);
                (m.inner, props)
            })
            .collect();
        omgkit_match::run_reactants(&self.inner, &inputs, max_products, atom_mapping)
            .into_iter()
            .map(|o| PyOutcome {
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
            })
            .collect()
    }

    /// 把整个反应物侧当作**一张图**上的查询来跑,而不是按位置配对。
    ///
    /// [`run`](Self::run) 要求"第 i 个分子配第 i 个模板片段",于是模板片段数
    /// 比分子数多时直接返回空 —— 而那正是**分子内反应**的形状:两个片段落在
    /// 同一个分子上。本方法把输入拼成一张图,让每个片段在整张图上自由找位置,
    /// 只要求各片段匹配到的原子两两不重叠。
    ///
    /// - 分子间:与 `run` 结果一致,且不必再枚举输入的排列
    /// - 分子内:`run` 表达不了的那一档
    /// - 盐:阳离子与阴离子是同一个分子的两个组分,模板可以同时碰到
    ///
    /// 代价是搜索空间变大,耗时不如 `run` 可预测;要稳定耗时就用 `run`。
    #[pyo3(signature = (substrate, *, max_products = 0, atom_mapping = false))]
    fn run_on_substrate(
        &self,
        substrate: Vec<PyMol>,
        max_products: usize,
        atom_mapping: bool,
    ) -> Vec<PyOutcome> {
        let inputs: Vec<(MolBuilder, omgkit_match::MolProps)> = substrate
            .into_iter()
            .map(|m| {
                let props = omgkit_match::MolProps::compute(&m.inner);
                (m.inner, props)
            })
            .collect();
        omgkit_match::run_on_substrate(&self.inner, &inputs, max_products, atom_mapping)
            .into_iter()
            .map(|o| PyOutcome {
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
            })
            .collect()
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
    m.add_class::<PyQuery>()?;
    m.add_class::<PyReaction>()?;
    m.add_class::<PyOutcome>()?;
    m.add_function(wrap_pyfunction!(parse_smiles, m)?)?;
    m.add_function(wrap_pyfunction!(parse_smarts, m)?)?;
    m.add_function(wrap_pyfunction!(parse_reaction, m)?)?;
    Ok(())
}
