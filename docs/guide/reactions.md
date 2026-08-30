# Reaction templates

A template describes the **reaction centre**. Everything else in the molecule
comes along on its own, which is why a template can be written small:
`[C:1][OH:2]>>[C:1][Cl:2]` says "a hydroxyl becomes a chloride" and says nothing
about the rest of the molecule.

![Esterification](../assets/esterification.svg)

*The template used here is
`[C:1](=[O:2])[OH:3].[OH:4][C:5]>>[C:1](=[O:2])[O:4][C:5]`. The benzene ring is
never mentioned in it; the water on the right is not in it either — that comes
from [byproduct reconstruction](byproducts.md).*

```pycon
>>> acid = omgkit.parse_smiles("CC(=O)O"); acid.sanitize()
>>> amine = omgkit.parse_smiles("CCN"); amine.sanitize()
>>> rxn = omgkit.parse_reaction("[C:1](=[O:2])[OH].[N:3]>>[C:1](=[O:2])[N:3]")
>>> rxn.num_reactant_templates, rxn.num_product_templates
(2, 1)
>>> out = rxn.run([acid, amine])[0]
>>> [p.to_canonical_smiles() for p in out.products]
['CC(NCC)=O']
```

Each reactant template gets a **different** input molecule. Mismatched counts
give an empty list — see [Intramolecular](#intramolecular) below for the case
that needs.

!!! tip "The order you hand the molecules in does not decide whether it runs"

    Position is not chemistry. `rxn.run([amine, acid])` gives the same products
    as `rxn.run([acid, amine])`: the engine tries "template *i* with molecule
    *i*" first, and only if that yields nothing does it look for another
    one-to-one assignment. When your order already lines up, the cost is exactly
    the same as trying that one assignment alone.

    So an empty list means one thing only: **there is no reaction site on these
    molecules.**

## Atom mapping

```pycon
>>> out = rxn.run([acid, amine], atom_mapping=True)[0]
>>> [r.to_smiles() for r in out.reactants]
['[CH3:1][C:2](=[O:3])O', '[CH3:4][CH2:5][NH2:6]']
>>> [p.to_smiles() for p in out.products]
['[C:2](=[O:3])([N:6][CH2:5][CH3:4])[CH3:1]']
```

`Outcome.reactants` is filled only when you ask. The two sides together are a
complete atom-mapped reaction — every product atom carries the number of the
reactant atom it came from, including atoms the template never mentioned.

Numbers are assigned to the **whole** reactant side, not just the template's
matched atoms, so the mapping is total rather than partial.

## Product count comes from the graph, not from the template

A template rewrites **one graph**. How many product molecules come out is
decided by how many connected components that rewritten graph has — not by how
many product templates were written.

For most templates this agrees with the common implementation. They diverge
exactly when the template cuts a bond in a ring that extends beyond the
template. The other implementation duplicates the atoms outside the template
into both products; omgkit gives one molecule, because the rewritten graph is
still connected.

This is a [deliberate divergence](../dev/correctness.md#a-deliberate-divergence),
with the extent measured on a real corpus.

## Intramolecular

`run` gives each template fragment a **different** molecule. When there are more
template fragments than molecules — two fragments landing on the *same*
molecule — that shape cannot be expressed, and `run` returns nothing.

`run_on_substrate` treats the whole reactant side as one graph and lets each
fragment find its own place, requiring only that the fragments' matched atoms
do not overlap:

```python
out = rxn.run_on_substrate([one_molecule])
```

| | `run` | `run_on_substrate` |
|---|---|---|
| Intermolecular | yes | yes |
| Sensitive to the order you pass molecules in | no | no |
| Intramolecular | no | yes |
| Salts (cation and anion as components of one molecule) | no | yes |
| Cost | predictable | larger search space, less predictable |

Use `run` when you want predictable timing.

## Limiting output

`max_products=n` caps how many outcomes are generated; `0` means no limit. A
template that matches in many places on a symmetric molecule can produce a lot
of equivalent outcomes.

## Byproducts

Pass `byproducts=True` to get the discarded fragments closed into molecules.
That has its own page: [Byproduct reconstruction](byproducts.md).
