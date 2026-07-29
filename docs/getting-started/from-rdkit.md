# Coming from RDKit

If you already know RDKit, this page is the fastest way in. Every omgkit output
below was produced by running the code.

## The operations you use most

| Task | RDKit | omgkit |
|---|---|---|
| Parse SMILES | `Chem.MolFromSmiles(s)` | `omgkit.parse_smiles(s)` |
| Parse without sanitizing | `Chem.MolFromSmiles(s, sanitize=False)` | `omgkit.parse_smiles(s)` — parsing never sanitizes |
| Sanitize | `Chem.SanitizeMol(m)` | `m.sanitize()` |
| Write SMILES | `Chem.MolToSmiles(m, canonical=False)` | `m.to_smiles()` |
| Write canonical SMILES | `Chem.MolToSmiles(m)` | `m.to_canonical_smiles()` |
| Atom count | `m.GetNumAtoms()` | `m.num_atoms` |
| Bond count | `m.GetNumBonds()` | `m.num_bonds` |
| Atomic numbers | `[a.GetAtomicNum() for a in m.GetAtoms()]` | `m.atomic_nums` |
| Copy | `Chem.Mol(m)` | `m.copy()` |
| Remove explicit H | `Chem.RemoveHs(m)` | `m.remove_hs()` |
| Parse SMARTS | `Chem.MolFromSmarts(s)` | `omgkit.parse_smarts(s)` |
| Match | `m.GetSubstructMatches(q)` | `q.match(m)` |
| Parse reaction | `AllChem.ReactionFromSmarts(s)` | `omgkit.parse_reaction(s)` |
| Run reaction | `rxn.RunReactants((a, b))` | `rxn.run([a, b])` |

## Four differences that will bite you

### 1. Parsing never sanitizes

`Chem.MolFromSmiles` sanitizes by default and returns `None` when that fails.
`omgkit.parse_smiles` only parses, and raises `ValueError` only when the string
itself is malformed. Sanitization is a separate, explicit step.

```python
m = omgkit.parse_smiles("OC(=O)c1ccccc1N")   # parsed, not sanitized
m.sanitize()                                  # now aromaticity, rings, H counts
```

This means an unsanitized molecule is a normal, usable object — you can inspect
it, and you decide when the pipeline runs. It also means **you must remember to
call `sanitize()`** before anything that depends on perceived properties
(matching ring or aromaticity criteria, canonical output, reactions).

### 2. Errors are exceptions, not `None`

RDKit returns `None` and writes to a log. omgkit raises `ValueError`, with a
caret pointing at the problem:

```pycon
>>> omgkit.parse_smiles("CC(C")
Traceback (most recent call last):
  ...
ValueError: unclosed branch
    CC(C
      ^
```

There is no silent-failure path, so there is no `if m is None` to forget.

### 3. `sanitize()` modifies in place and can leave a partial result

RDKit's `SanitizeMol` also mutates, but the idiom of passing molecules around
by value hides it. Here it is explicit: if you need all-or-nothing, copy first.

```python
safe = m.copy()
try:
    m.sanitize()
except ValueError:
    m = safe
```

The binding deliberately does not copy for you — a hidden deep copy would
double the cost of every batch with no way for you to know.

### 4. Product count comes from the graph, not from the template

This is the one that actually changes answers.

When a reaction template is applied, RDKit produces one product molecule per
**product template**. omgkit rewrites one graph and then counts its **connected
components**. For most templates these agree. They diverge exactly when the
template cuts a bond in a ring that extends beyond the template — RDKit
duplicates the atoms outside the template into both products; omgkit returns
one molecule, because the rewritten graph is still connected.

The [Correctness](../dev/correctness.md#a-deliberate-divergence) page has the
measured extent of this on a real corpus. It is a deliberate divergence, not an
incompatibility to work around.

## What omgkit does that RDKit does not

**Byproduct reconstruction.** When a template discards atoms, omgkit tells you
exactly which ones (`Outcome.discarded`) and can rebuild them into balanced
molecules (`Outcome.byproducts`) — or say explicitly that the record does not
balance:

```pycon
>>> rxn = omgkit.parse_reaction("[C:1](=[O:2])[OH].[N:3]>>[C:1](=[O:2])[N:3]")
>>> out = rxn.run([acid, amine], byproducts=True)[0]
>>> [p.to_canonical_smiles() for p in out.byproducts], out.byproduct_verdict
(['O'], 'capped')
```

**Whole-substrate matching.** `run_on_substrate` treats the entire reactant
side as one graph, so a template whose two fragments land on the *same*
molecule — an intramolecular reaction — works without special handling.

## What RDKit has that omgkit does not

Quite a lot. omgkit is a focused toolkit, not a replacement. Missing here:
descriptors and fingerprints, conformers and 3D, file formats other than
SMILES/SMARTS, drawing, force fields, and a very long tail of chemistry
utilities.

The Python surface is also narrower than the Rust one — `MolBatch` is not
wrapped yet. See [Python API](../api/python.md#not-yet-exposed).
