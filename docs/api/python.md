# Python API

Everything the `omgkit` module exposes. Signatures and descriptions are
generated from the source, so they cannot drift away from the code.

!!! info "Import once"

    ```python
    import omgkit
    ```

    There are no submodules to reach into — the three parse functions and the
    five classes below are the whole surface.

## At a glance

| Callable | What it gives you |
|---|---|
| [`parse_smiles`](#omgkit.parse_smiles) | a [`Mol`](#omgkit.Mol) from a SMILES string |
| [`parse_smarts`](#omgkit.parse_smarts) | a [`Query`](#omgkit.Query) from a SMARTS pattern |
| [`parse_reaction`](#omgkit.parse_reaction) | a [`Reaction`](#omgkit.Reaction) from a reaction SMARTS |
| [`Mol`](#omgkit.Mol) | a molecule — sanitize it, write it back out |
| [`Query`](#omgkit.Query) | a substructure query — match it against a molecule |
| [`Reaction`](#omgkit.Reaction) | a reaction template — run it on reactants |
| [`Outcome`](#omgkit.Outcome) | one result of running a reaction |
| [`Conformer`](#omgkit.Conformer) | a 3D structure — coordinates plus the molecule they belong to |

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

## Not yet exposed

The Rust side has more than the Python side. These are reachable from Rust
today and are not yet wrapped:

| Rust | What it is |
|---|---|
| `omgkit_core::MolBatch` | the columnar batch and its zero-copy per-molecule views |
| `omgkit_io::smarts` writing | SMARTS output for molecules and reactions |
| `omgkit_chem` individual stages | running one sanitization stage at a time |

See the [Rust API](rust.md) if you need them.
