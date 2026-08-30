# The sanitization pipeline

```python
m = omgkit.parse_smiles("OC(=O)c1ccccc1N")
m.sanitize()
```

One call runs the whole pipeline. It modifies in place and raises `ValueError`
if the molecule cannot be sanitized — an impossible valence, an aromatic ring
that cannot be kekulized.

## What it does, in order

| Stage | What it settles |
|---|---|
| Valence | how many bonds each atom can carry, given its element and charge |
| Implicit hydrogens | the hydrogen count that follows from valence and the bonds present |
| Ring perception | which atoms and bonds are in rings, ring sizes, ring membership counts |
| Kekulization | assigning alternating single/double bonds to aromatic rings |
| Aromaticity | which rings actually *are* aromatic, as opposed to written lower-case |
| Conjugation | which bonds are part of a conjugated system |
| Hybridization | per-atom hybridization state |

Then a final step converts the SMILES *directional* bonds (`/` and `\`) into
the double bonds' own cis/trans property. Direction is a property of the
single bonds around the double bond in the input syntax; geometry is a property
of the double bond itself. Keeping the two apart is what lets the writer put
the slashes back in a different but equivalent place.

## What runs on it afterwards

Sanitization only *fills in* the molecule; nothing here computes anything a
caller did not ask for. Two things read that output and are deliberately not
part of the pipeline:

| | Needs from sanitization |
|---|---|
| [Descriptors for ML](descriptors.md) | aromaticity, ring membership, hybridization, conjugation, implicit hydrogen counts, plus the cis/trans step above |
| [Drawing](depict.md) and [3D structures](conformers.md) | the same, plus kekulization |

Calling them on an unsanitized molecule does not raise. It hands back a full set
of values made of parse-time placeholders — every atom `unspecified`, every bond
geometry `none` — which is the failure mode worth knowing about, because nothing
about the result looks wrong.

## The order is not arbitrary

Each stage consumes what the ones before it produced. Aromaticity perception
needs rings; kekulization needs valences and hydrogen counts; hydrogen counts
need valences, which need charges. Running them out of order does not error —
it produces a molecule that is internally inconsistent, which is far harder to
notice.

## Charge changes valence, and how it changes depends on the element

A positive charge does not always add a valence slot and a negative charge does
not always remove one. Nitrogen gains one when positive; carbon loses one.
Oxygen loses one when negative; boron gains one when negative.

Getting the sign right but the element wrong is the kind of mistake that
produces a plausible molecule with the wrong hydrogen count. `omgkit_chem`
exposes this as a table rather than a rule.

## Failure leaves a partial result

```python
safe = m.copy()
try:
    m.sanitize()
except ValueError:
    m = safe
```

The binding does not add that copy for you. A hidden deep copy would double the
cost of every batch and the caller would have no way to know it was happening.
