# Reading and writing SMILES

## Parsing

```python
m = omgkit.parse_smiles("O[C@@H](C)c1ccccc1")
```

Parsing **does not sanitize**. What you get back is exactly what the string
said: the atoms in the order they were written, the bonds as written, aromatic
lower-case atoms recorded as aromatic because the author claimed so — not
because anything has verified it yet.

Malformed input raises `ValueError` with a caret view:

```pycon
>>> omgkit.parse_smiles("CC(C")
Traceback (most recent call last):
  ...
ValueError: CC(C
    ^ 括号不匹配
```

## Two ways to write it back

| | Depends on | Use it for |
|---|---|---|
| `to_smiles()` | the current atom storage order | round-tripping, debugging, seeing what you built |
| `to_canonical_smiles()` | the structure only | keys, deduplication, comparison |

```pycon
>>> m = omgkit.parse_smiles("OC(=O)c1ccccc1N")
>>> m.to_smiles()
'OC(=O)c1ccccc1N'
>>> m.sanitize(); m.to_canonical_smiles()
'c1cccc(c1C(O)=O)N'
```

`to_canonical_smiles()` is the one that keeps the promise: any two writings of
the same molecule give the same string.

```pycon
>>> a = omgkit.parse_smiles("OC(=O)c1ccccc1N"); a.sanitize()
>>> b = omgkit.parse_smiles("Nc1ccccc1C(O)=O"); b.sanitize()
>>> a.to_canonical_smiles() == b.to_canonical_smiles()
True
```

## What survives a round trip

Tetrahedral chirality, double-bond geometry, dative bonds, isotopes, formal
charges, radical electrons, atom map numbers and explicit hydrogens all parse
and write back.

Chirality deserves a note. A tetrahedral tag (`@` / `@@`) is **relative to the
order the neighbours are stored in**. Any operation that changes that order —
removing a hydrogen, cutting a bond, reordering for canonical output — has to
rebase the tag, or the molecule silently becomes its mirror image. Nothing
raises; you just get the wrong enantiomer. This is why the writing layer has
its own differential judge for stereochemistry rather than trusting round-trip
equality alone.

## Explicit hydrogens

```pycon
>>> h = omgkit.parse_smiles("[H]C([H])([H])O"); h.sanitize()
>>> h.num_atoms
5
>>> h.remove_hs()
3
>>> h.num_atoms, h.to_smiles()
(2, 'CO')
```

`remove_hs()` folds hydrogens into their neighbour's hydrogen count and returns
how many it removed. **It modifies in place and every atom index changes** —
deleting atoms necessarily renumbers.

It deliberately keeps hydrogens that carry information: isotopes, charges, atom
map numbers, radicals, bridging hydrogens, and hydrogens that carry a double
bond's direction. Keeping one extra hydrogen only costs a node in the graph;
removing the wrong one loses information.
