# omgkit

[![CI](https://github.com/zbc0315/omgkit/actions/workflows/ci.yml/badge.svg)](https://github.com/zbc0315/omgkit/actions/workflows/ci.yml)

[中文说明](README.zh-CN.md)

## What is it?

omgkit is a cheminformatics toolkit written in Rust, with Python bindings. It is
built around a columnar molecule representation so that large batches stay cheap
to traverse, split across threads, and hand off to numpy or Arrow without a copy.

  * [BSD-3-Clause license](LICENSE)
  * Core data structures and algorithms in Rust, no `unsafe`
  * Python 3.9+ wrapper built with [PyO3](https://pyo3.rs) and
    [maturin](https://github.com/PyO3/maturin) — one abi3 wheel covers every
    supported version, and there are no system dependencies
  * SMILES parsing and writing, including tetrahedral chirality, double-bond
    geometry, dative bonds and explicit hydrogens
  * Canonical SMILES
  * A sanitization pipeline: valence, implicit hydrogens, ring perception,
    kekulization, aromaticity, conjugation, hybridization
  * SMARTS parsing, substructure matching (VF2++, optionally stereo-aware), and
    SMARTS writing for both molecules and reactions
  * Reaction templates and product generation, with optional atom-atom mapping
  * Reconstruction of the fragments a template discards (the water an
    esterification drops) into balanced byproduct molecules — or an explicit
    "cannot tell" when the record itself does not balance
  * Columnar batches (`MolBatch`) with zero-copy per-molecule views

**Status: under development.** The API still changes between commits. Every
layer is checked against an external reference implementation record by record
(see [Documentation](#documentation)), but the surface is not yet stable enough
for production use. Bug reports are welcome.

## Installation

### Python

Requires [maturin](https://github.com/PyO3/maturin):

```shell-session
$ maturin build --release -m crates/omgkit-py/Cargo.toml --out dist
$ pip install dist/omgkit-*.whl
```

### Rust

```toml
[dependencies]
omgkit-core  = { git = "https://github.com/zbc0315/omgkit" }   # data structures
omgkit-io    = { git = "https://github.com/zbc0315/omgkit" }   # SMILES / SMARTS
omgkit-chem  = { git = "https://github.com/zbc0315/omgkit" }   # sanitization
omgkit-match = { git = "https://github.com/zbc0315/omgkit" }   # matching, reactions
```

Take only the layers you need; each depends only on the ones below it.

## Getting started

```python
import omgkit

m = omgkit.parse_smiles("OC(=O)c1ccccc1N")
m.sanitize()
m.to_canonical_smiles()

q = omgkit.parse_smarts("[C](=[O])[OH]")
q.match(m)                      # molecule atom indices, in query atom order

rxn = omgkit.parse_reaction("[C:1][OH:2]>>[C:1][Cl:2]")
for outcome in rxn.run([m], atom_mapping=True):
    outcome.products, outcome.reactants
```

The Rust equivalents live in `omgkit_io::smiles`, `omgkit_chem::sanitize` and
`omgkit_match`; see the crate documentation for runnable examples.

## Documentation

  * [`docs/design.md`](docs/design.md) — what each layer does, why it is built
    that way, and how each design choice was validated
  * [`harness/README.md`](harness/README.md) — the differential-testing setup:
    how the oracles are generated and how a test is kept from passing
    vacuously
  * `cargo doc --workspace --no-deps --open` — API documentation

## Contributing

Issues and pull requests are welcome. Four gates have to pass, and they are the
same four that CI runs:

```shell-session
$ cargo fmt --all --check
$ cargo clippy --workspace --all-targets -- -D warnings
$ cargo test
$ cargo doc --workspace --no-deps
```

`cargo test` is green on a fresh clone: the smoke oracles are committed. The
large-corpus tier is marked `#[ignore]` and needs oracles you generate yourself
— see [`harness/README.md`](harness/README.md).

## License

Code released under the [BSD-3-Clause license](LICENSE).

Test corpora and the element table are redistributed from other projects and
carry their own terms; each file is traced to its origin in
[`THIRD-PARTY-NOTICES.md`](THIRD-PARTY-NOTICES.md).
