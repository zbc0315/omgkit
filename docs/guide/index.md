# Guides

One page per capability, organised by what you are trying to do rather than by
which crate the code lives in.

<div class="grid cards" markdown>

- **[Reading and writing SMILES](smiles.md)**

    Parsing, storage order vs canonical output, chirality and double-bond
    geometry, explicit hydrogens.

- **[The sanitization pipeline](sanitize.md)**

    What each stage does, in what order, and why the order matters.

- **[Substructure matching](matching.md)**

    SMARTS queries, match ordering, uniquification, and stereo-aware matching.

- **[Reaction templates](reactions.md)**

    Applying templates, atom mapping, intramolecular reactions, and where the
    product count comes from.

- **[Byproduct reconstruction](byproducts.md)**

    Rebuilding the fragments a template discards — and reading the verdict that
    tells you how much to trust the answer.

- **[Reading and writing `.mol`/`.sdf`](molfiles.md)**

    V2000 molblocks and multi-record SDF, in 2D and 3D, and how
    stereochemistry survives the round trip.

- **[Drawing structures](depict.md)**

    2D coordinates and SVG/PNG/JPEG output, the two drawing styles, and what
    gets reported when a structure cannot be drawn well.

- **[3D structures](conformers.md)**

    One deterministic conformer per molecule — why there is no random seed,
    and what the numbers that come back with it mean.

- **[Descriptors for ML](descriptors.md)**

    The per-atom and per-bond values a graph neural network reads — including
    Gasteiger partial charges — and why they come back as names, not one-hot
    vectors.

- **[Batches](batches.md)**

    The columnar representation, and what it is for.

</div>

## The one idea behind all of them

**A molecule's properties are decided by the molecule, not by how someone chose
to write it down.**

Every layer is checked against that. Two SMILES for the same structure must
sanitize to the same thing, match the same queries, and react the same way.
Where an implementation choice could break it — the order neighbours happen to
be stored in, which atom a template happened to be written around — the design
notes say so explicitly, because those are exactly the places where bugs are
invisible from the outside.
