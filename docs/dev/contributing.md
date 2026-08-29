# Contributing

Issues and pull requests are welcome. The project is under development and the
API still moves, so bug reports are especially useful.

## Before you open a pull request

```shell
bash harness/gates.sh
```

That is the whole suite, and it needs a Python environment with the pinned
RDKit (`harness/requirements.lock`) because most gates compare against it. The
Rust-only part, which needs none of that, is:

```shell
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --release
cargo test --workspace
cargo doc --workspace --no-deps --document-private-items
```

## If you touch the layout, regenerate the figures

Every structure in the documentation is drawn by `omgkit-depict` itself:

```shell
python3 docs/figures/make_figures.py
```

Two runs are byte-identical — the layout has no random component — so a diff in
`docs/assets/` after a run means the change moved something on the page.

!!! warning "No gate watches this"

    Regenerating the figures is **not** part of `harness/gates.sh` and not part
    of CI: `draw` needs the `raster` feature, which pulls in resvg. So a layout
    change that alters what the figures look like will pass every gate while the
    committed SVGs quietly go stale. Until that is gated, it is on the author to
    re-run the script and look at the diff.

## If you touch the element table, regenerate it — do not hand-edit

`crates/omgkit-core/src/element_data.rs` is generated, and it is large: 119
elements, 3111 isotope masses, 93 Pauling electronegativities, the early-atom
table. It carries a `@generated` header for that reason.

```shell
python3 harness/gen_elements.py --rdkit ../rdkit --out crates/omgkit-core/src/element_data.rs
```

It needs an RDKit **source checkout** at the version pinned in
`harness/requirements.lock`; the script refuses to run if the source and the
installed package disagree, because a table from one version compared against
baselines from another produces divergences that are not real.

!!! warning "No gate watches this either, and here is what covers it instead"

    Re-generating and diffing cannot be a CI gate — CI has no RDKit source tree.
    What guards the table is a unit test in `omgkit-core` pinning a few textbook
    values, **which elements have no value at all**, and the row counts. That
    catches a shifted or truncated table, which is the realistic failure: the
    parser once silently dropped H-1 and Ca-47 because two of the sixty-odd raw
    string blocks put their first data row on the same physical line as the
    opening `R"DAT(`.

    What it does not catch is a single digit changed in place. Nothing does —
    RDKit exposes no API for the electronegativity table, so there is no
    independent reader to compare against. Do not hand-edit the file.

## What a change needs

**A test that goes red when the change is reverted.** Please actually revert it
and watch the test fail. A test that has only ever been seen passing might be
checking nothing — that has happened here, and the fix was to write a different
test, not to trust the first one.

**A reason, in the code.** Where a decision could reasonably have gone the other
way, say why it went this way. The comments in this repository lean heavy on
*why* and light on *what*, on purpose: the *what* is right there in the code,
and the *why* is what gets lost.

**Nothing silently narrowed.** If a change caps coverage — a top-N, a skipped
retry, a sampling step — say so in the output. Silent truncation reads as full
coverage when it is not.

## Style

Comments and documentation in this repository are currently mostly Chinese on
the Rust side and English in the user-facing documentation. Either is fine in a
pull request; matching the surrounding file is better than switching mid-file.

`unsafe_code = "deny"` across the workspace. `missing_docs` is on, and rustdoc's
broken-link lints are `deny` — a broken intra-doc link fails the build.

## Reporting a bug

The most useful report has the SMILES or SMARTS that triggers it, what you
expected, and what you got. If the difference is against another toolkit,
saying which one and which version helps — some divergences here are
[deliberate](correctness.md#a-deliberate-divergence) and it is worth knowing
quickly which kind you have hit.
