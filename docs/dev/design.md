# Design

The layer stack, and the invariant each layer has to keep.

!!! note "Full text"

    This page is an overview. The complete design record — every trade-off, why
    it went that way, and how it was validated — is
    [`docs/design.md`](design-full.md), in Chinese.

## The one invariant

**A molecule's properties are decided by the molecule, not by how someone chose
to write it down.**

Every layer is checked against that. It sounds obvious and it is where most of
the hard cases live: the order neighbours happen to be stored in, which atom a
template happened to be written around, whether an aromatic ring was written
lower-case or kekulized. Each of those is a way for the *writing* to leak into
the *answer*.

## Two representations

| Type | Mutability | Layout | For |
|---|---|---|---|
| `MolBuilder` | mutable | array of structs, one molecule | building: parsing, constructing products |
| `MolBatch` | immutable | struct of arrays + CSR, across molecules | algorithms, parallelism, zero-copy export |

Conversion goes both ways and the round trip is identity, checked by test.

## The layers

| Layer | What it does | The thing it must not get wrong |
|---|---|---|
| **L1** SMILES parsing | string → graph, verbatim | record what the author *said*, not what is true — aromaticity here is a claim, not a perception |
| **L2** Sanitization | valence → implicit H → rings → kekulization → aromaticity → conjugation → hybridization | the order; each stage consumes the previous one's output |
| **L3** SMILES writing | graph → string | stereo tags are relative to neighbour order, so any reordering must rebase them |
| **L3** Canonical ordering | a numbering that depends only on the structure | two writings of one molecule must land on one string |
| **L4** SMARTS parsing | pattern → query | a query is not a molecule; primitives that look alike behave differently |
| **L5** Substructure matching | VF2++ | chirality is a property of the **mapping**, not of an atom pair |
| **L7** Reactions | template application | product count comes from connected components of the rewritten graph |
| **L8** Python bindings | PyO3 | do not silently add safety the Rust side does not have |

## Three decisions worth knowing about

### Aromaticity at L1 is a claim, not a fact

A lower-case atom in the input says *the author thinks this is aromatic*. It is
recorded as such and left alone. Only L2 decides whether it actually is. Mixing
the two means an input's formatting can change the perceived chemistry.

### Stereo tags are relative, so every reordering rebases

A tetrahedral tag is defined against the order the neighbours are stored in.
Removing a hydrogen, cutting a bond, reordering for canonical output — each
changes that order, and each therefore has to rebase the tag. Getting it wrong
raises nothing; it produces the mirror image. Every operation that reorders
neighbours has a stereo-specific judge behind it for this reason.

### Product count comes from the graph

A reaction template rewrites one graph. All product templates are built into
that one graph, atoms outside the template are moved exactly once, and the
result is split by connected components at the end.

This diverges from the common implementation exactly when a template cuts a
bond in a ring that extends beyond the template — see
[Correctness](correctness.md#a-deliberate-divergence). The criterion that
settles it is mass conservation: the product must not contain more heavy atoms
than the substrate.

## Complexity is a first-class concern

Where an algorithm's cost matters, there is a test guarding the exponent, not
just a comment claiming it. Those guards measure interleaved and are calibrated
by injecting a real defect and confirming the guard goes red — a threshold
nobody has seen fail is a threshold that might be checking nothing.

## Dependencies

Deliberately few. The Python wheel has no system dependencies: one abi3 wheel,
nothing to install alongside it, no shared library to find at runtime.
