# Python API

Everything the `omgkit` module exposes. Signatures and descriptions are
generated from the source, so they cannot drift away from the code.

!!! info "Import once"

    ```python
    import omgkit
    ```

    There are no submodules to reach into — the three parse functions and the
    four classes below are the whole surface.

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

## Not yet exposed

The Rust side has more than the Python side. These are reachable from Rust
today and are not yet wrapped:

| Rust | What it is |
|---|---|
| `omgkit_core::MolBatch` | the columnar batch and its zero-copy per-molecule views |
| `omgkit_io::smarts` writing | SMARTS output for molecules and reactions |
| `omgkit_chem` individual stages | running one sanitization stage at a time |

See the [Rust API](rust.md) if you need them.
