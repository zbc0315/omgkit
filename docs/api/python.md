# Python API

Everything the `omgkit` module exposes. Signatures and descriptions are
generated from the source, so they cannot drift away from the code.

!!! info "Import once"

    ```python
    import omgkit
    ```

    There are no submodules to reach into — the five reading functions and the
    seven classes below are the whole surface.

## At a glance

| Callable | What it gives you |
|---|---|
| [`parse_smiles`](#omgkit.parse_smiles) | a [`Mol`](#omgkit.Mol) from a SMILES string |
| [`parse_smarts`](#omgkit.parse_smarts) | a [`Query`](#omgkit.Query) from a SMARTS pattern |
| [`parse_reaction`](#omgkit.parse_reaction) | a [`Reaction`](#omgkit.Reaction) from a reaction SMARTS |
| [`parse_molblock`](#omgkit.parse_molblock) | a [`Molblock`](#omgkit.Molblock) from the contents of a `.mol` file |
| [`read_sdf`](#omgkit.read_sdf) | a list of [`SdfRecord`](#omgkit.SdfRecord) from the contents of an `.sdf` file |
| [`Mol`](#omgkit.Mol) | a molecule — sanitize it, write it back out |
| [`Query`](#omgkit.Query) | a substructure query — match it against a molecule |
| [`Reaction`](#omgkit.Reaction) | a reaction template — run it on reactants |
| [`Outcome`](#omgkit.Outcome) | one result of running a reaction |
| [`Conformer`](#omgkit.Conformer) | a 3D structure — coordinates plus the molecule they belong to |
| [`Molblock`](#omgkit.Molblock) | one record read from a `.mol`/`.sdf` file — molecule plus its coordinates |
| [`SdfRecord`](#omgkit.SdfRecord) | one entry of an `.sdf` — that molblock, its data fields, or why it could not be read |

## 3D coordinates

[`Mol.conformer()`](#omgkit.Mol.conformer) turns a molecule into one 3D
structure. It is deterministic — no random seed, no retries you have to
configure; the same molecule always comes back with the same coordinates.

```pycon
>>> import omgkit
>>> conf = omgkit.parse_smiles("C[C@H](N)C(=O)O").conformer()
>>> conf
<omgkit.Conformer atoms=13 energy=0.000e0 converged=True chiral=1/1>
>>> conf.mol.num_atoms          # 13, not 6 — generation adds explicit hydrogens
13
>>> conf.coords[0]
(-1.191..., -0.899..., -0.073...)
```

!!! warning "The coordinates belong to `conf.mol`, not to the molecule you called it on"

    Generation needs explicit hydrogens, so it works on a copy and adds them
    there. Your molecule is left untouched; `conf.coords` lines up with
    `conf.mol`, which has more atoms.

`ValueError` is raised when the molecule cannot be sanitized, or when its
distance bounds contradict each other (about 1 molecule in 8831 on a drug-like
corpus). The message says which.

### Saving it

[`Conformer.to_molblock()`](#omgkit.Conformer.to_molblock) gives you the
contents of a `.mol` file. An `.sdf` is those records separated by `$$$$`:

```python
with open("out.sdf", "w") as f:
    for smi in ["CCO", "C[C@H](N)C(=O)O"]:
        conf = omgkit.parse_smiles(smi).conformer()
        f.write(conf.to_molblock(title=smi))
        f.write("$$$$\n")
```

Aromatic bonds are kekulized on the way out — a molblock has no aromatic bond
type, and writing one as a single bond would turn thiophene into
tetrahydrothiophene without saying so. The second line of the block is the
program name with **no timestamp**, so writing the same molecule twice gives
byte-identical output.

## Reading `.mol` files

[`parse_molblock`](#omgkit.parse_molblock) reads one V2000 molblock — the
contents of a `.mol` file, or one record of an `.sdf` up to its `$$$$`:

```pycon
>>> rec = omgkit.parse_molblock(open("aminoethanol.mol").read())
>>> rec
<omgkit.Molblock 4 atoms 2D "C[C@H](N)O">
>>> rec.mol.to_canonical_smiles()   # the wedge bond in the file is what makes this @
'C[C@H](N)O'
>>> rec.is_3d, len(rec.coords), rec.coords[0]
(False, 4, (-1.299, -0.75, 0.0))
```

!!! info "It sanitizes for you — and it has to"

    The other parse functions hand back an **un-sanitized** molecule and leave
    the timing to you. This one cannot. A molblock keeps half its stereochemistry
    in the **coordinates and wedge bonds**, which live outside `Mol`; reading them
    onto the atoms needs implicit-hydrogen counts and symmetry classes, and both
    of those come out of sanitization. Leaving that step to the caller would mean
    a forgotten call silently drops the stereochemistry of the whole file — no
    error, same atom count, just no `@` and no `/`.

!!! warning "3D files: the molecule comes back without stereochemistry"

    Stereochemistry in a 3D file lives in the coordinates themselves, which is a
    different route and is not implemented yet. When `is_3d` is true the molecule
    carries **no stereo tags at all** — that is "not implemented", not "this
    molecule is flat". The 2D route (wedges for chirality, coordinates for
    double-bond geometry) is complete.

Bonds the file marks as *explicitly unknown* — a crossed double bond, a wavy
single bond — are left unassigned rather than being read off the drawing.

## Reading `.sdf` files

[`read_sdf`](#omgkit.read_sdf) reads every record of an SDF at once:

```python
for i, rec in enumerate(omgkit.read_sdf(open("library.sdf").read())):
    if rec.error:
        print(f"record {i}: {rec.error}")
        continue
    print(rec.block.mol.to_canonical_smiles(), dict(rec.data))
```

!!! warning "A record that cannot be read does not raise, and does not vanish"

    Raising would stop at the bad record and throw away everything after it.
    Skipping would make the count quietly smaller than the file's — the caller
    counts records and gets a different number, with nothing anywhere reporting
    it. So every record keeps its place in the list: a bad one has `error` set
    to a sentence and `block` set to `None`, and the ones after it read fine.

    Real files have these. A ferrocene-type complex has more bonds on the metal
    than V2000 can express, so writers emit V3000 for it — and V3000 is refused
    here rather than misread.

`data` is a **list of pairs, not a dict**: repeated field names do occur (a
vendor writing one line per measurement), and a dict would silently keep only
the last. Values spanning several lines are joined with `\n`, and a line that
*starts* with `>` inside a value is part of the value, not a new field.

The whole file is parsed at once — budget peak memory accordingly for very
large libraries.

Writing a 2D molblock back out is not available from Python: choosing where the
wedges go is part of the drawing layer, which the bindings do not carry.
Writing a **3D** one is — see [`Conformer.to_molblock`](#omgkit.Conformer.to_molblock).

## Errors

Every parse function raises `ValueError` on malformed input.
[`Mol.sanitize`](#omgkit.Mol.sanitize) raises `ValueError` when the molecule
cannot be sanitized — for instance an impossible valence.

`parse_smiles` puts a caret view in the message pointing at the offending
character:

```pycon
>>> omgkit.parse_smiles("CC(C")
Traceback (most recent call last):
  ...
ValueError: unclosed branch
    CC(C
      ^
```

!!! warning "`sanitize` is in place and may leave a partial result"

    If sanitization fails, the molecule may already have been modified. Callers
    that need all-or-nothing should `copy()` first. The binding deliberately
    does **not** add that copy for you: a hidden deep copy would double the cost
    of every batch, and the caller would have no way to know.

---

## Parsing

::: omgkit.parse_smiles

::: omgkit.parse_smarts

::: omgkit.parse_reaction

::: omgkit.parse_molblock

::: omgkit.read_sdf

---

## Mol

::: omgkit.Mol
    options:
      members: true

---

## Query

::: omgkit.Query
    options:
      members: true

---

## Reaction

::: omgkit.Reaction
    options:
      members: true

---

## Outcome

::: omgkit.Outcome
    options:
      members: true

---

## Conformer

::: omgkit.Conformer
    options:
      members: true

---

## Molblock

::: omgkit.Molblock
    options:
      members: true

---

## SdfRecord

::: omgkit.SdfRecord
    options:
      members: true

---

## Not yet exposed

The Rust side has more than the Python side. These are reachable from Rust
today and are not yet wrapped:

| Rust | What it is |
|---|---|
| `omgkit_core::MolBatch` | the columnar batch and its zero-copy per-molecule views |
| `omgkit_io::smarts` writing | SMARTS output for molecules and reactions |
| `omgkit_chem` individual stages | running one sanitization stage at a time |
| `omgkit_io::molblock::write_v2000` with wedges | writing a **2D** molblock (needs the drawing layer to place the wedges) |
| SDF **writing** | multi-record output with data fields — only single-record 3D output exists today |

See the [Rust API](rust.md) if you need them.
