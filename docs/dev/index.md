# Developing

<div class="grid cards" markdown>

- **[Building and testing](building.md)** — the gates, and how to run the
  differential tests
- **[Design](design.md)** — the layer stack and the invariants each layer keeps
- **[Correctness](correctness.md)** — how the claims are checked, and how a
  test is kept from passing vacuously
- **[Contributing](contributing.md)** — what a change needs before it lands

</div>

## The short version

```shell
git clone https://github.com/zbc0315/omgkit
cd omgkit
cargo test --release
```

That is green on a fresh clone. The smoke oracles for the differential tests
are committed, so nothing has to be generated first.

The whole suite is one command:

```shell
bash harness/gates.sh
```

It needs a Python environment with the pinned RDKit
(`harness/requirements.lock`), because most gates compare omgkit against it
record by record. The Rust-only part, which runs without any of that, is:

```shell
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --release
cargo test --workspace
cargo doc --workspace --no-deps --document-private-items
```

(How many gates there are is deliberately not written here. The count has gone
from four to five to forty as judges were added, and a number in prose falls
behind without anything reporting it — `TOTAL` in `harness/gates.sh` is the one
place it lives. See [Building and testing](building.md).)

## The house rules

**Every fix needs a test that goes red when the fix is reverted.** Not a test
that passes — a test that has been shown to fail without the change. A test
nobody has seen fail is a test that might be checking nothing.

**A judge must prove it does not pass vacuously.** Several of the checks here
carry a `--self-test` mode that injects known defects and confirms the judge
catches them. Zero divergences is only good news if the comparison actually
ran; "zero divergences over zero comparisons" is the most convincing-looking
green there is.

**Confirm an injected defect was actually injected.** A self-test that reports
"the judge missed it" is more often a self-test whose injection was a no-op.
This has happened more than once here and it is documented where it did.

**Numbers in comments are gated too.** Where a comment says "measured N cases",
a test re-measures N. Otherwise the comments drift away from the code and
nobody notices.
