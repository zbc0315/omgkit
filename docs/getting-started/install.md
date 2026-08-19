# Installation

## Python

There is no wheel on PyPI yet, so build one. You need a Rust toolchain
(1.75 or newer) and [maturin](https://github.com/PyO3/maturin):

```shell
git clone https://github.com/zbc0315/omgkit
cd omgkit

pip install maturin
maturin build --release -m crates/omgkit-py/Cargo.toml --out dist
pip install dist/omgkit-*.whl
```

Check it:

```pycon
>>> import omgkit
>>> omgkit.parse_smiles("CCO").num_atoms
3
```

**One wheel covers every supported Python.** The extension is built against the
stable ABI (abi3, Python 3.9+), so a wheel built on 3.9 runs on 3.13 as well.
There are no system dependencies — nothing to `apt install`, no shared library
to find at runtime.

## Rust

Add the layers you need. Each depends only on the ones below it, so taking
`omgkit-io` alone gets you SMILES parsing without pulling in the reaction
engine.

```toml
[dependencies]
omgkit-core  = "0.0.1"   # data structures
omgkit-io    = "0.0.1"   # SMILES / SMARTS
omgkit-chem  = "0.0.1"   # sanitization
omgkit-match = "0.0.1"   # matching, reactions
```

Or with `cargo add`:

```shell
cargo add omgkit-core omgkit-io omgkit-chem omgkit-match
```

| Crate | Depends on | Gives you |
|---|---|---|
| `omgkit-core` | — | `MolBuilder`, `MolBatch`, `MolView`, the scalar types |
| `omgkit-io` | core | SMILES and SMARTS parsing and writing |
| `omgkit-chem` | core | the sanitization pipeline, valence, aromaticity |
| `omgkit-match` | core, io, chem | VF2++ matching, reaction templates, byproducts |
| `omgkit-py` | all of the above | the Python extension module (not a library) |

## Requirements

| | |
|---|---|
| Rust | 1.75 or newer (`rust-version` in the manifest) |
| Python | 3.9 or newer |
| Platform | anything Rust and PyO3 support; no system libraries |

## Building from a clone

The gates below are what CI runs. They pass on a fresh clone — the smoke
oracles for the differential tests are committed, so nothing has to be
generated first.

```shell
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --release
cargo test --workspace
cargo doc --workspace --no-deps --document-private-items

# the external judges for omgkit-conf, on the committed smoke baseline
cargo run -p omgkit-conf --release --example smooth_oracle -- harness/baseline/smoke.bounds.jsonl
cargo run -p omgkit-conf --release --example bounds_oracle -- harness/baseline/smoke.bounds.jsonl
cargo run -p omgkit-conf --release --example eigen_oracle  -- harness/baseline/smoke.bounds.jsonl harness/baseline/smoke.gram_eigs.jsonl
```

The large-corpus tier of the differential tests is marked `#[ignore]` and needs
oracles you generate yourself. See [Building and testing](../dev/building.md).
