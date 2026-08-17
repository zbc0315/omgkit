//! 只做一件事:告诉链接器"Python 符号先欠着"。
//!
//! 扩展模块被解释器 `dlopen` 进来时,`Py_*` 那些符号由**宿主解释器**提供,
//! 所以构建期不该去链 `libpython` —— 真去链了,得到的模块会绑死在某一个
//! Python 安装上;而且很多系统 Python(包括 macOS 自带的)根本不提供
//! 可链接的 libpython,构建当场就失败。
//!
//! 用 `rustc-cdylib-link-arg` 而不是 `.cargo/config.toml` 里的 rustflags:
//! 前者只作用于**本 crate 的 cdylib 目标**,后者会泼到整个工作区 ——
//! 那等于让所有测试可执行文件也不再检查未定义符号,把一类链接期错误
//! 推迟到运行期才炸。

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    // 其它平台各有各的办法(Linux 默认就允许未定义符号,Windows 必须真链),
    // 这里只处理 macOS。构建 wheel 时 maturin 会另行处理。
    if std::env::var("CARGO_CFG_TARGET_OS").as_deref() == Ok("macos") {
        println!("cargo:rustc-cdylib-link-arg=-undefined");
        println!("cargo:rustc-cdylib-link-arg=dynamic_lookup");
    }
}
