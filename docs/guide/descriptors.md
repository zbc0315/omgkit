# Graph descriptors for machine learning

The per-atom and per-bond quantities a graph neural network reads, in one call
each, computed the same way the reference implementation computes them.

| | |
|---|---|
| [`Mol.atom_descriptors()`](../api/python.md#omgkit.Mol.atom_descriptors) | twelve values per atom |
| [`Mol.bond_descriptors()`](../api/python.md#omgkit.Mol.bond_descriptors) | seven values per bond |

```pycon
>>> import omgkit
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

## What you get

### Per atom

| Key | Type | |
|---|---|---|
| `atomic_num` | `int` | element; `0` is the wildcard `*` |
| `total_degree` | `int` | explicit neighbours **plus** total hydrogens |
| `formal_charge` | `int` | |
| `chiral_tag` | `str` | `unspecified`, `cw`, `ccw`, `allene`, `square_planar`, `trigonal_bipyramidal`, `octahedral` |
| `total_num_hs` | `int` | declared + implied; **not** free-standing `[H]` atoms in the graph |
| `hybridization` | `str` | `unspecified`, `s`, `sp`, `sp2`, `sp3`, `sp2d`, `sp3d`, `sp3d2` |
| `is_aromatic` | `bool` | |
| `is_in_ring` | `bool` | |
| `mass` | `float` | the exact mass of the labelled nuclide if an isotope was given, otherwise the standard atomic weight |
| `electronegativity` | `float \| None` | Pauling |
| `gasteiger_charge` | `float` | Gasteiger–Marsili partial charge (PEOE) |
| `gasteiger_valid` | `bool` | whether the line above is a real number |

### Per bond

| Key | Type | |
|---|---|---|
| `begin`, `end` | `int` | atom indices |
| `order` | `str` | `unspecified`, `single`, `double`, `triple`, `quadruple`, `aromatic`, `dative` |
| `is_conjugated` | `bool` | |
| `is_in_ring` | `bool` | |
| `stereo` | `str` | `none`, `cis`, `trans` (see below) |
| `stereo_atoms` | `tuple[int, int] \| None` | the two atoms `stereo` is measured against |

## Descriptors, not an encoding

Categorical values come back as **names** — `"sp3"`, `"ccw"` — never one-hot
vectors and never integer codes. Which elements go in your vocabulary, whether
you keep an "other" bucket, how you scale the three continuous values: those are
decisions belonging to your featurizer, and freezing one model's answer into a
library would make it wrong for the next model. Integer codes would be worse
still: they are invisible when they shift.

Building a 45-dimensional atom vector on top is a dozen lines:

```python
ELEMENTS = ["C", "N", "O", "F", "P", "S", "Cl", "Br", "I"]
ATOMIC_NUMS = [6, 7, 8, 9, 15, 16, 17, 35, 53]

def one_hot(value, vocabulary):
    """Last slot is the 'other' bucket."""
    v = [0] * (len(vocabulary) + 1)
    v[vocabulary.index(value) if value in vocabulary else -1] = 1
    return v

def atom_vector(d):
    return (
        one_hot(d["atomic_num"], ATOMIC_NUMS)
        + one_hot(d["total_degree"], list(range(6)))
        + one_hot(d["formal_charge"], [-2, -1, 0, 1, 2])
        + one_hot(d["hybridization"], ["sp", "sp2", "sp3", "sp3d", "sp3d2"])
        + [d["is_aromatic"] * 1, d["is_in_ring"] * 1]
        + [d["mass"] * 0.01]
        + [(d["electronegativity"] or 2.0) / 4]
        + [d["gasteiger_charge"] if d["gasteiger_valid"] else 0.0]
        + [d["gasteiger_valid"] * 1]
    )
```

## Two ways a value can be missing

`electronegativity` is `None` when the element **has no accepted Pauling value**
— the noble gases, and Pm, Eu, Tb, Yb, Fr. That is different from "not measured
yet" and different from zero.

`gasteiger_valid` is `False` when the atom falls outside the Gasteiger parameter
set, which covers H, C, N, O, F, Si, P, S, Cl, Br, I, B, Be, Mg and Al. Most
metals do not, and the failure **spreads along the graph**: a carbon bonded to a
sodium ends up with a `nan` charge too. Across the four corpora the gate runs on
(9 051 molecules), 4 192 atoms came out invalid — and RDKit produces `nan` for
exactly the same ones.

Both are reported honestly rather than filled in with a default. A default would
merge "we do not know" with "the value happens to be that", and your featurizer
is exactly the code that needs to tell those apart — which is what the
`gasteiger_valid` flag is for.

## `cis`/`trans`, not `Z`/`E`

`Z` and `E` are defined by CIP priority, and omgkit does not implement CIP
ranking. What you get instead is the geometry relative to the two atoms named in
`stereo_atoms`, matching RDKit's `SetBondStereoFromDirections`
(`STEREOCIS`/`STEREOTRANS`) rather than `AssignStereochemistry`'s
`STEREOZ`/`STEREOE`.

The two carry the **same geometric information**: given the reference atoms you
can convert one to the other by ranking them yourself. That is why
`stereo_atoms` is part of the descriptor and not an extra — read on its own,
`stereo` means nothing, because on a tetrasubstituted double bond a different
choice of reference flips the label.

## Sanitize first

Aromaticity, ring membership, hybridization, conjugation and implicit hydrogen
counts are all filled in by sanitization; double-bond geometry is filled in by
the step right after it. `Mol.sanitize()` does both.

Skipping it does not raise — you get a full set of descriptors made of parse-time
placeholders.

## The Rust API

```rust
use omgkit_chem::{atom_descriptors, bond_descriptors, sanitize};

let mut m = omgkit_io::smiles::parse("C/C=C/C(=O)O")?;
sanitize(&mut m)?;
omgkit_io::stereo::perceive_bond_stereo(&mut m);   // see below

let a = atom_descriptors(&m);   // Vec<AtomDescriptors>
let b = bond_descriptors(&m);   // Vec<BondDescriptors>
assert!(a[0].gasteiger_is_valid());
```

The enums come back as enums here (`Hybridization::Sp3`,
`BondStereo::Trans`), not as strings — the strings exist only at the Python
boundary, where there is no enum to hand over.

!!! warning "Rust callers must make that second call themselves"

    `Mol.sanitize()` in Python runs **two** things: `omgkit_chem::sanitize`, and
    then `omgkit_io::stereo::perceive_bond_stereo`, which turns the `/` and `\`
    of the input into the double bond's own cis/trans property. The second lives
    in `omgkit-io`, a sibling of `omgkit-chem` rather than a dependency, so
    `sanitize` cannot call it.

    Skipping it does not raise. Every bond just reports `BondStereo::None`, and
    nothing about the result looks wrong.

## How this is checked

Every one of the sixteen descriptors is compared atom-by-atom and bond-by-bond
against RDKit over four corpora (9 051 molecules) by
`harness/check_descriptors.py`, which runs in CI. Seven deliberate divergences
are pinned by name, each one a documented blind spot on the reference side — for
instance RDKit cannot represent allene axial chirality at all.

The one thing that gate cannot see is the **value** of the Pauling
electronegativity, since RDKit exposes no API for it; that is pinned by a unit
test in `omgkit-core` instead. Both facts are written down in the judge's header,
along with the mutation runs that establish them.
