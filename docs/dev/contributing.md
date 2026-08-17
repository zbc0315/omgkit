# Contributing

Issues and pull requests are welcome. The project is under development and the
API still moves, so bug reports are especially useful.

## Before you open a pull request

```shell
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --release
cargo test --workspace
cargo doc --workspace --no-deps --document-private-items
```

All five must pass. They are what CI runs.

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
