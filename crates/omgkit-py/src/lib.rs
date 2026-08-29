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
        let rec = omgkit_io::molblock::Record {
            title,
            coords: &coords,
            wedges: &d.wedges,
            orders: &orders,
        };
        omgkit_io::molblock::write_v2000(&grown, &rec)
            .map_err(|e| PyValueError::new_err(e.to_string()))
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
        let rec = omgkit_io::molblock::Record {
            title,
            coords: &self.conf.coords,
            wedges: &[],
            orders: &orders,
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
    /// 分子。**已经净化过,立体也打上了**(二维图的情形)。
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
/// 三维文件同样眼下不读立体,理由见 `parse_molblock`。
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
    /// `reactants[i]` 要与第 i 个反应物模板对应,数目不符时返回空列表。
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
