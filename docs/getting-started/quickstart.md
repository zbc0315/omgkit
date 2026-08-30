# Quickstart

Every output on this page was produced by running the code. If yours differs,
that is a bug — please report it.

```python
import omgkit
```

## Parse and write

```pycon
>>> m = omgkit.parse_smiles("OC(=O)c1ccccc1N")
>>> m.num_atoms, m.num_bonds
(10, 10)
>>> m.to_smiles()
'OC(=O)c1ccccc1N'
```

`to_smiles` writes in the current atom storage order, so it gives back what you
put in. To get a form that depends only on the structure, sanitize first and
ask for the canonical one:

```pycon
>>> m.sanitize()
>>> m.to_canonical_smiles()
'c1cccc(c1C(O)=O)N'
```

**Canonical means canonical.** Two ways of writing the same molecule land on
the same string:

```pycon
>>> a = omgkit.parse_smiles("OC(=O)c1ccccc1N"); a.sanitize()
>>> b = omgkit.parse_smiles("Nc1ccccc1C(O)=O"); b.sanitize()
>>> a.to_canonical_smiles() == b.to_canonical_smiles()
True
```

That is the invariant the whole library is built to keep: a molecule's
properties are decided by the molecule, not by how someone chose to write it
down.

## Sanitize

`sanitize()` runs the whole pipeline — valence, implicit hydrogens, ring
perception, kekulization, aromaticity, conjugation, hybridization — and then
converts directional bonds into the double bonds' own cis/trans property.

```pycon
>>> m = omgkit.parse_smiles("c1ccccc1")
>>> m.sanitize()
>>> m.to_canonical_smiles()
'c1ccccc1'
```

!!! warning "It modifies in place, and a failure can leave a partial result"

    ```python
    safe = m.copy()
    try:
        m.sanitize()
    except ValueError:
        m = safe          # m may have been partly modified
    ```

    The binding does not make that copy for you on purpose — a hidden deep copy
    would double the cost of every batch, and you would have no way to know.

## Match a substructure

```pycon
>>> m = omgkit.parse_smiles("OC(=O)c1ccccc1N"); m.sanitize()
>>> q = omgkit.parse_smarts("[CX3](=O)[OX2H1]")
>>> q.match(m)
[[1, 2, 0]]
```

Each inner list is one match, giving molecule atom indices **in query atom
order**. Here query atom 0 (the carbonyl carbon) is molecule atom 1, query
atom 1 (the `=O`) is molecule atom 2, and query atom 2 (the `OH`) is molecule
atom 0.

`uniquify=False` keeps symmetry-equivalent duplicates. Benzene shows the
difference — six aromatic bonds, each matchable in two directions:

```pycon
>>> bz = omgkit.parse_smiles("c1ccccc1"); bz.sanitize()
>>> len(omgkit.parse_smarts("cc").match(bz))
6
>>> len(omgkit.parse_smarts("cc").match(bz, uniquify=False))
12
```

`max_matches=n` stops after *n* matches; `0` (the default) means no limit.

## Run a reaction template

```pycon
>>> rxn = omgkit.parse_reaction("[C:1](=[O:2])[OH].[N:3]>>[C:1](=[O:2])[N:3]")
>>> rxn.num_reactant_templates, rxn.num_product_templates
(2, 1)

>>> acid  = omgkit.parse_smiles("CC(=O)O");  acid.sanitize()
>>> amine = omgkit.parse_smiles("CCN");      amine.sanitize()

>>> out = rxn.run([acid, amine])[0]
>>> [p.to_canonical_smiles() for p in out.products]
['CC(NCC)=O']
```

`reactants[i]` goes with reactant template *i*. If the counts do not match you
get an empty list back — see [`run_on_substrate`](../guide/reactions.md#intramolecular)
for the case where two template fragments land on one molecule.

### Atom mapping

```pycon
>>> out = rxn.run([acid, amine], atom_mapping=True)[0]
>>> [r.to_smiles() for r in out.reactants]
['[CH3:1][C:2](=[O:3])O', '[CH3:4][CH2:5][NH2:6]']
>>> [p.to_smiles() for p in out.products]
['[C:2](=[O:3])([N:6][CH2:5][CH3:4])[CH3:1]']
```

The two sides together are a complete atom-mapped reaction: every product atom
carries the number of the reactant atom it came from.

### Byproducts

The template above throws away the acid's `OH`. Ask for it back:

```pycon
>>> out = rxn.run([acid, amine], byproducts=True)[0]
>>> [p.to_canonical_smiles() for p in out.byproducts]
['O']
>>> out.byproduct_verdict
'capped'
>>> out.discarded
[[3], []]
```

`discarded` is the **fact**: atom 3 of the first input went into no product.
`byproducts` is the **inference**: those atoms, closed into a real molecule.
`byproduct_verdict` says how it was closed, and it can say it could not be:

| Verdict | Meaning |
|---|---|
| `'off'` | you did not pass `byproducts=True` |
| `'nothing'` | no atoms were discarded |
| `'capped'` | closed by adding hydrogens only — no choices, so high confidence |
| `'bonded(n)'` | *n* extra bonds were formed; which atoms they join is a heuristic |
| `'unresolved(...)'` | could not be closed, with the reason |

**When the verdict is `unresolved`, `byproducts` is empty.** Making one up would
be worse than giving nothing: it would be a well-formed molecule with nothing
wrong on its face.

See [Byproduct reconstruction](../guide/byproducts.md) for what the budget
means and when the answer is trustworthy.

## Generate a 3D structure

```pycon
>>> conf = omgkit.parse_smiles("C[C@H](N)C(=O)O").conformer()
>>> conf
<omgkit.Conformer atoms=13 energy=0.000e0 converged=True chiral=1/1>
>>> conf.coords[0]
(-1.1906..., -0.8985..., -0.0731...)
```

No random seed and no retry loop — the same molecule always gives the same
coordinates. `chiral=1/1` says every stereocentre came out with the right sign.
The extra atoms are the explicit hydrogens generation adds; `conf.coords` lines
up with `conf.mol`, not with the molecule you called it on. See
[3D structures](../guide/conformers.md).

## Read and write `.mol` files

```pycon
>>> m = omgkit.parse_smiles("C[C@H](N)C(=O)O"); m.sanitize()
>>> block = m.to_molblock_2d(title="L-alanine")   # layout and wedges computed for you
>>> back = omgkit.parse_molblock(block)
>>> back.mol.to_canonical_smiles() == m.to_canonical_smiles()
True
```

`read_sdf` reads a whole multi-record file, and a record it cannot read keeps
its place in the list instead of raising or vanishing. See
[Reading and writing `.mol`/`.sdf`](../guide/molfiles.md).

## Featurize for a model

```pycon
>>> m = omgkit.parse_smiles("CC(=O)Oc1ccccc1C(=O)O"); m.sanitize()
>>> m.atom_descriptors()[1]
{'atomic_num': 6, 'total_degree': 3, 'formal_charge': 0,
 'chiral_tag': 'unspecified', 'total_num_hs': 0, 'hybridization': 'sp2',
 'is_aromatic': False, 'is_in_ring': False, 'mass': 12.011,
 'electronegativity': 2.55, 'gasteiger_charge': 0.3075...,
 'gasteiger_valid': True}
>>> m.bond_descriptors()[1]
{'begin': 1, 'end': 2, 'order': 'double', 'is_conjugated': True,
 'is_in_ring': False, 'stereo': 'none', 'stereo_atoms': None}
```

Twelve values per atom and seven per bond, including Gasteiger partial charges.
Categorical ones come back as names, not one-hot vectors — the vocabulary is
your featurizer's decision. `electronegativity` is `None` where the element has
no accepted Pauling value and `gasteiger_valid` is `False` where the charge
could not be computed; neither is filled in with a default. See
[Descriptors for ML](../guide/descriptors.md).

## Where to go next

- [Guides](../guide/index.md) — one page per capability
- [Python API](../api/python.md) — every callable, with signatures
