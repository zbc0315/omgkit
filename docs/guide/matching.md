# Substructure matching

```pycon
>>> m = omgkit.parse_smiles("OC(=O)c1ccccc1N"); m.sanitize()
>>> q = omgkit.parse_smarts("[CX3](=O)[OX2H1]")
>>> q.match(m)
[[1, 2, 0]]
```

Each inner list is one match, giving **molecule atom indices in query atom
order**. Query atom 0 is the carbonyl carbon, which is molecule atom 1; query
atom 2 is the hydroxyl oxygen, which is molecule atom 0.

That ordering is the useful one: it lets you go straight from "the third atom
of my pattern" to "this atom of the molecule" without re-deriving the
correspondence.

## Options

```python
q.match(mol, uniquify=True, max_matches=0, use_chirality=True)
```

**`uniquify`** collapses matches that differ only by a symmetry of the query.
Benzene makes the difference visible — six aromatic bonds, each matchable from
either end:

```pycon
>>> bz = omgkit.parse_smiles("c1ccccc1"); bz.sanitize()
>>> len(omgkit.parse_smarts("cc").match(bz))
6
>>> len(omgkit.parse_smarts("cc").match(bz, uniquify=False))
12
```

**`max_matches`** stops the search early. `0` means no limit.

**`use_chirality`** makes matching stereo-aware.

## Chirality is a property of the *mapping*

This is the part that is easy to get wrong. Whether a match respects
stereochemistry cannot be decided by looking at the query atom and the molecule
atom in isolation — a tetrahedral tag is relative to neighbour order, so
whether the tags agree depends on **which query neighbour mapped to which
molecule neighbour**. The check therefore belongs to the mapping, not to the
atom pair.

An implementation that compares tags atom-by-atom will accept mirror images.
It will also pass every test that only uses queries whose neighbours happen to
be listed in the same order as the molecule's — which is most hand-written
tests. This is why the stereo judge in the test harness has to prove it can
tell a molecule from its enantiomer before its results count for anything.

## Query properties are recomputed every call

`match` recomputes the molecule's query-relevant properties (ring membership
counts, smallest ring size, and so on) on each call rather than caching them.

That is deliberate. `sanitize()` and friends modify a molecule in place; a
stale cache would give silently wrong answers, and a wrong answer is far more
expensive than recomputing. If you are matching many patterns against one
molecule and this shows up in a profile, that is a good bug report.
