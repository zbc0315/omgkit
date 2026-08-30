# Reading and writing `.mol` and `.sdf` files

V2000 molblocks and multi-record SDF, in both directions and in both
dimensionalities. Stereochemistry survives the round trip.

| | |
|---|---|
| [`parse_molblock(text)`](../api/python.md#omgkit.parse_molblock) | one molblock — the contents of a `.mol`, or one record of an `.sdf` |
| [`read_sdf(text)`](../api/python.md#omgkit.read_sdf) | every record of an SDF, in order |
| [`Mol.to_molblock_2d()`](../api/python.md#omgkit.Mol.to_molblock_2d) | write a **2D** molblock — layout and wedges are computed for you |
| [`Conformer.to_molblock()`](../api/python.md#omgkit.Conformer.to_molblock) | write a **3D** molblock from generated coordinates |

## Writing 2D

```pycon
>>> m = omgkit.parse_smiles("C[C@H](N)C(=O)O"); m.sanitize()
>>> print(m.to_molblock_2d(title="L-alanine"))
L-alanine
  omgkit

  6  5  0  0  0  0  0  0  0  0999 V2000
   -1.0000   -0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5000   -0.0000    0.0000 C   0  0  0  0  0  4  0  0  0  0  0  0
   -1.0000    0.8660    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000   -0.8660    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.8660    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  6  0  0  0
  2  3  1  0  0  0  0
  2  4  1  0  0  0  0
  4  5  2  3  0  0  0
  4  6  1  0  0  0  0
M  END
```

The `6` in the fourth column of the first bond line is a **down wedge** — that
is where the `@` went. `to_molblock_2d` runs the whole chain for you: sanitize,
perceive double-bond geometry, lay out coordinates, assign wedges, and write
**the molecule that was drawn**.

!!! info "A double bond with no configuration is written as **crossed**"

    Laying a molecule out gives **every** double bond a definite geometry — the
    substituents have to go somewhere. If the input never said whether a bond is
    cis or trans, writing that geometry without saying "unknown" hands the reader
    a configuration the author never gave.

    So a double bond that could have a configuration but does not gets stereo
    code `3` (a crossed double bond, the V2000 way of saying *either*). Bonds
    that already have a configuration are written with their real geometry, and
    bonds that cannot have one — inside a small ring, or with the same two
    substituents on one end — are not marked at all, since "unknown" would be
    just as untrue there.

    This was measured: before the marks were written, **551 of 8831 molecules
    (6.2%)** came back from a round trip carrying cis/trans that the input never
    had. The gate that watches it is `harness/check_molblock_roundtrip.py`.

    The same applies to 3D: an embedded conformer gives every double bond a
    definite torsion, so `Conformer.to_molblock()` marks them too.

!!! info "The molecule written may have one more atom than the one you passed"

    Some stereocentres have all three of their heavy bonds inside a ring; the
    only wedge that can express the configuration is on an explicit C–H, so the
    layout adds one. That hydrogen is in the file, and the wedge is on it. A
    centre whose configuration cannot be drawn is left **unmarked** rather than
    given an arbitrary wedge.

## Writing 3D

[`Conformer.to_molblock()`](conformers.md#writing-it-out) writes generated
coordinates. There the stereochemistry is in the coordinates themselves and the
wedge column is empty. Both kinds of file are valid; a reader tells them apart
by whether any `z` is non-zero.

## Reading

```pycon
>>> block = m.to_molblock_2d(title="L-alanine")   # or: open("L-alanine.mol").read()
>>> rec = omgkit.parse_molblock(block)
>>> rec
<omgkit.Molblock 6 atoms 2D "L-alanine">
>>> rec.is_3d, rec.title
(False, 'L-alanine')
>>> rec.mol.to_canonical_smiles()
'C[C@@H](C(O)=O)N'
>>> rec.coords[0]
(-1.0, -0.866, 0.0)
```

`Molblock` keeps the molecule and its coordinates **together** on purpose. Half
the stereochemistry of a molblock lives in the coordinates and wedges, which are
not part of `Mol`; handing back only the molecule would drop that half without a
sound.

### Stereochemistry, both dimensionalities

Reading stereochemistry off a file takes two completely separate routes, and
which one runs is decided by the coordinates rather than by the caller:

| | chirality | double-bond geometry |
|---|---|---|
| **2D** | the wedge bonds | which side of the double-bond axis each reference atom is on |
| **3D** | the signed volume of the four ligands | the sign of the torsion |

```pycon
>>> conf = omgkit.parse_smiles("C[C@H](N)C(=O)O").conformer()
>>> back = omgkit.parse_molblock(conf.to_molblock())
>>> back.is_3d
True
>>> back.mol.to_canonical_smiles() == conf.mol.to_canonical_smiles()
True
```

!!! tip "Comparing a round trip: call `remove_hs()` first"

    Writing may add an explicit hydrogen to carry a wedge, so the molecule that
    comes back can have one more atom than the one you wrote. Merge those back
    before comparing:

    ```python
    back = omgkit.parse_molblock(m.to_molblock_2d()).mol
    back.remove_hs()
    back.to_canonical_smiles() == m.to_canonical_smiles()
    ```

    Across the 8831-molecule corpus this holds for 8829; the two that differ are
    trivalent phosphorus, where omgkit reads a centre the external
    implementation does not.

Bonds the file marks as **explicitly unknown** — a crossed double bond, a wavy
single bond — are left unassigned instead of being read off the drawing. "The
author says the configuration is unknown" is information, and turning it into a
configuration by measuring the picture destroys it.

!!! info "It sanitizes for you — and it has to"

    The other parse functions hand back an **un-sanitized** molecule and leave
    the timing to you. This one cannot: reading wedges onto atoms needs
    implicit-hydrogen counts and symmetry classes, and both come out of
    sanitization. Leaving that step to the caller would mean one forgotten call
    silently drops the stereochemistry of a whole file — no error, same atom
    count, just no `@` and no `/`.

## Reading a whole SDF

```python
for i, rec in enumerate(omgkit.read_sdf(open("library.sdf").read())):
    if rec.error:
        print(f"record {i}: {rec.error}")
        continue
    print(rec.block.mol.to_canonical_smiles(), rec.data)
```

!!! warning "A record that cannot be read does not raise, and does not vanish"

    Raising would stop at the bad record and throw away everything after it.
    Skipping would make the count quietly smaller than the file's — the caller
    counts records, gets a different number, and nothing anywhere reports it.

    So every record keeps its place: a bad one has `error` set to a sentence and
    `block` set to `None`, and the records after it read fine.

    ```pycon
    >>> one = omgkit.parse_smiles("CCO"); one.sanitize()
    >>> two = omgkit.parse_smiles("c1ccccc1"); two.sanitize()
    >>> def record(mol, smi):
    ...     return mol.to_molblock_2d() + "> <SMILES>" + chr(10) + smi + chr(10) * 2 + "$$$$" + chr(10)
    >>> text = record(one, "CCO") + record(two, "c1ccccc1")
    >>> text += one.to_molblock_2d()[:40]     # third record: cut off mid-file
    >>> recs = omgkit.read_sdf(text)          # three records, the last one truncated
    >>> [(r.error, r.data) for r in recs]
    [(None, [('SMILES', 'CCO')]),
     (None, [('SMILES', 'c1ccccc1')]),
     ('molblock 截断了:缺M  END', [])]
    ```

    Real files have these. A ferrocene-type complex has more bonds on the metal
    than V2000 can express, so writers emit V3000 for it — and V3000 is refused
    here rather than misread as V2000.

### Data fields are a list of pairs, not a dict

```pycon
>>> recs[0].data
[('SMILES', 'CCO')]
```

Repeated field names do occur — a vendor writing one line per measurement — and
a dict would silently keep only the last. Whether the names in *your* files are
unique is your call to make: `dict(rec.data)` is one line away.

## Writing an SDF

An SDF record is a molblock, then data fields, then `$$$$`:

```python
with open("out.sdf", "w") as f:
    for smi in ["CCO", "c1ccccc1"]:
        m = omgkit.parse_smiles(smi); m.sanitize()
        f.write(m.to_molblock_2d(title=smi))
        f.write(f"> <SMILES>\n{smi}\n\n")
        f.write("$$$$\n")
```

The Rust API has `omgkit_io::molblock::write_sdf_record`, which does the same
and **refuses** a data value that would break the record boundary — a blank line
or a lone `$$$$` inside a value would silently split one record into two.

## What is not supported

**V3000 is refused, not misread.** A V3000 file raises with a message saying so.
Writing V3000 is not implemented, so molecules with more than 999 atoms or bonds
cannot be written as a molblock at all — that raises too.

## Determinism

Writing the same molecule twice gives **byte-identical** output. The second line
of the block is the program name with no timestamp, and the 2D layout has no
random component: any SMILES for the same structure gives the same coordinates.
