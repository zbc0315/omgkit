# Correctness

Chemical semantics is full of undocumented edge cases. **You cannot get them
right by reading the specification.** So every layer here is compared, record by
record, against an external reference implementation, over a large corpus.

!!! note "Full text"

    This page is an overview. The complete testing record — every judge, what it
    guards, and the traps found while building it — is
    [`harness/README.md`](correctness-full.md), in Chinese. It is the longest
    document in the repository for a reason.

## A judge must prove it does not pass vacuously

This is the rule everything else follows from.

Zero divergences is only good news **if the comparison actually ran**. "Zero
divergences over zero comparisons" is the most convincing-looking green there
is, and a corpus change or a convention change can produce it silently.

So the judges here carry two things:

**A `--self-test` mode** that injects known defects and confirms they are
caught. Flip a stereo tag, drop an atom, add a hydrogen — if the judge still
says *identical*, that tier of the comparison is theatre.

**A `--min-checked` floor.** If the number of records actually compared falls
below it, the run fails, even with zero divergences.

### And confirm the injected defect was injected

A self-test that reports *the judge missed it* is, more often than not, a
self-test whose injection did nothing.

A worked example from this repository: a byproduct judge was checked by adding
one hydrogen to the byproduct. The balance rate only fell from 100% to 95.9%,
which looked like a judge passing vacuously. It was not. The most common
byproducts are HCl and H₂O; adding a hydrogen gives ClH₂ and OH₃, sanitization
fails, and the injection function returned the input unchanged. The injection
was a no-op on the bulk of the corpus. The fix was to count only the records
where the injection actually changed something, and to add a second injection
that works on every molecule.

## Coverage is a table, not a feeling

The hardest defects are not "a judge got it wrong" — they are the ones **in the
gaps between judges**. One judge tests writing before sanitization, another
tests the object after sanitization, and nothing walks the path *writing after
sanitization*.

So the coverage table matters more than any single judge:

| Path | Guarded by |
|---|---|
| Parsing (L1) | `differential_l1*` |
| Each sanitization stage / whole pipeline (L2) | `differential_l2_*` (compares molecule **object fields**) |
| Writing, unsanitized | `roundtrip_smiles` + `check_ez.py` |
| Writing, **after sanitization** | `check_write_fidelity.py` |
| **Canonical** SMILES fidelity | `check_write_fidelity.py` + `dump_canonical` |
| Uniqueness of the canonical ordering | `canonical_invariance` (needs no external reference) |
| Double-bond geometry perception | `check_bond_stereo.py` |
| SMARTS parsing (L4) | `oracle_smarts.py` |
| SMARTS **writing** (L3) | `roundtrip_smarts` (idempotence) + `check_smarts_write.py` (semantics) |
| Substructure matching (L5) | `oracle_matches.py` |
| Product generation (L7) | `check_reactions.py` |
| Closing discarded atoms (L7) | `tests/byproduct.rs` — the criterion is mass conservation, independent of the record |
| Closing, **real corpus + external judge** | `check_byproducts.py` (with `--self-test`) |
| Product generation, **real corpus, both directions** | `bench_reactions.py` |
| Python bindings (L8) | `test_python.py` |
| SMARTS **chirality** reference frame | `check_smarts_chirality.py` (with a discriminating-power check) |
| Product-side **chirality**, four instruction kinds | `check_product_chirality.py` (with a discriminating-power check) |
| Drawing (L9) | `omgkit-depict --example audit` over the whole corpus — see below |

**When adding a path, ask this table first.** Which cell does it fall in? No
cell means a new gap.

## A deliberate divergence

`check_reactions.py` reports **717 identical / 24 different**. Those 24 are not
defects — they are a chosen semantics, and they all have the same shape:
a template like `[C:1][O:2][C:3]>>[C:1][O:2].[C:3]` applied to a **cyclic**
ether or lactone.

```text
substrate  [C@]12(C(OCC1C(=C)CC[C@@H]2O)=O)C          a bicyclic lactone
external   C=C1CC[C@H](O)[C@](C)(C(=O)O)C1 . C=C1CC[C@H](O)[C](C)C1C
           ↑ two molecules — the carbocycle is duplicated into both,
             so atoms appear from nowhere
omgkit     C=C1CC[C@H](O)[C@](C)(C(=O)O)C1C
           ↑ one molecule, ring-opened, not one atom more or fewer
```

A product template describes a **fragment** of the reaction centre, not "one
fragment, one molecule". Whether the fragments are actually separate depends on
whether the atoms *outside* the template still connect them.

The criterion that settles it is **mass conservation**: a product must not
contain more heavy atoms than the substrate. Copying shared atoms into each
product leaves topology, valence and stereochemistry all self-consistent — only
the atom count is wrong, and no other judge was looking at that number. So a
dedicated test was added for it (`products_never_invent_atoms`).

**If that count of 24 changes, investigate.** Going up means a new shape has
appeared; going down means the behaviour was reverted by accident.

## Chasing a hit rate assumes the ground truth is right

Real-corpus "ground truth" is what someone recorded, which is not the same as
chemical fact. The most common gap when applying templates in reverse is an
**underdetermined record**: the recorded reactant has no configuration marked
while the substrate determines it completely.

In that case, giving the product a configuration is **correct**. Demanding an
exact match against the underdetermined record forces the implementation to
throw information away. Cases like that are recorded as corpus gaps, not
counted as misses.

## Performance has gates too

Where complexity matters there is a test guarding the exponent. Those guards:

- measure **interleaved**, so drift in machine state hits all sizes equally
- run at sizes large enough for the term of interest to dominate
- have thresholds **calibrated by injecting a real defect** and confirming the
  guard goes red

That last one has its own trap. Calibrating one guard, the injected quadratic
term produced no slowdown at all — because the shape injected was
loop-invariant and the compiler hoisted it out of the inner loop. A threshold
calibrated against an injection that never ran is worth nothing.

## Drawing: six decidable properties, run over the whole corpus

A picture cannot be judged by "does it look right". What *can* be decided is
whether it says something false about the molecule. Six properties, all run
over `harness/corpus/large.smi` — 8831 molecules × 2 styles:

```shell
cargo run -p omgkit-depict --release --example audit -- harness/corpus/large.smi
```

| Property | Precondition |
|---|---|
| Writing-independent: any SMILES for the same molecule draws the same primitives | none |
| No two atoms drawn on the same point | none |
| No bond angle below 90° at an atom of degree ≤3 | layout not degraded |
| Ring double bonds: both lines land inside a ring that contains the bond | layout not degraded |
| Wedges reach the canvas: as many wedge primitives as `Depiction` recorded | none |
| Wedges read back: every drawn stereocentre reads back as its recorded configuration | none |
| Bond lengths all equal | layout not degraded |
| Nothing drawn outside the canvas | none |

The two preconditions are not excuses. A degraded layout has no reliable shape
to begin with — and *that it degraded* is already reported in
`Depiction::degraded`, so the caller is not being told a comfortable lie.

### Hand-picked molecules only cover the model you already have

Two whole classes of defect got through unit judges before the corpus run found
them:

- "ring double bonds stay inside the ring" listed eight molecules, all
  single-ring or ortho-fused — so **bridged shared bonds**, where both rings
  sit on the *same* side, were never exercised.
- "a wedge starts at the stereocentre it describes" listed three molecules,
  none with a P(V) centre — so a wedge assigned to a **double bond**, which the
  renderer silently drops while `unwedged` stays empty, was never exercised.

Neither needed a new idea to find. They needed a corpus.

### The audit nearly passed vacuously itself

Its writing-independence check first verifies that the rewritten SMILES is
still the same molecule. That guard used the storage-order writer, so almost
every molecule compared unequal and was **skipped** — the check reported 3
violations out of 17662 while doing nearly nothing. With the canonical writer
it reports 132. The audit now counts how many comparisons actually happened and
fails loudly if the answer is zero.

### …and then a second time: the shuffler was the identity one time in nine

Writing-independence is checked by **rewriting the SMILES and drawing it
again**. The rewrite worked by handing the writer a priority order, and that
order was built from a multiplicative hash:

```rust
(0..n).map(|i| (i * 2_654_435_761 * k) % (n * 7 + 13))
```

**For many `n` that is not a permutation at all.** Atoms collide onto the same
priority, the sort degenerates to the identity, and the "rewritten" SMILES is
character-for-character the original. Measured over the whole corpus:

| | |
|---|---:|
| rewrite identical to the original | **2874 / 26493 = 10.85%** |
| fixed points in the priority order | 18.28% (a uniform permutation gives ~1/n) |

Benzene came out as `c1ccccc1` all three times; ethanol as `CCO`. **The crate's
headline contract was being verified by a shuffler that did nothing one time in
nine.**

It is now splitmix64 + Fisher–Yates, retried up to 8 seeds until the canonical
labelling actually changes — if `canonical_ranks` comes back identical, the
storage order never moved and the comparison would be free.

| | old hash | real permutation |
|---|---:|---:|
| real violations found | 137 | **257** |
| cases never compared at all | **1156 (6.5%)** | **2** |

**More than half the violations were invisible.**

This also exposed a reporting defect in the judge itself: "could not check" was
being counted in the violations column. Swapping in a *worse* shuffler
therefore pushed the violation count from 259 up to 1293 — and all 1156 of the
increase were cases it had failed to check. **A bigger number looked like more
diligence and meant the opposite.** It now reports three separate rows: cases
checked, cases that got the full five comparisons, and cases never checked.

The shuffler has its own guard (`check_the_shuffler`, run at the top of `main`
rather than under `#[cfg(test)]`, because `cargo test` does not run an
example's tests). Mutation-verified: putting the multiplicative hash back
panics immediately with `n=2 seed=0: that is not a permutation`.

!!! warning "The numbers below predate the fix"

    Every figure in the sections that follow — 129, 125, 134 — was measured
    with the defective shuffler and reflects only what it could see. They
    remain valid as **relative** comparisons between successive changes, which
    is what they were for. They are not the absolute level; see
    [Current standing](#current-standing).

### One class of defect no unit test can see

The violation count moved between runs of the **same binary on the same
corpus**: 141, 142, 141. `HashMap`'s hasher is seeded randomly per process, so
iteration order changes, and the layout summed positions in that order —
last-bit differences flipped branches.

For a drawing library that is a hard failure: regenerate a figure and it
silently changes. All `HashMap`/`HashSet` in `omgkit-depict` are now
`BTreeMap`/`BTreeSet`, and three consecutive runs give the same number.

**A unit test cannot catch this** — within one process the seed is fixed. It
takes running the audit twice.

### A bond that stops too far from its label

Bond lines stop short of an atom label so they do not run through the letters.
They used to stop at the label box's **circumscribed circle**. A circle
contains the box, so a line can never touch the text — at the cost of stopping
far too early **in whichever direction the box is narrow**. Approaching a wide
label vertically wastes exactly `hypot(w,h) − h`. Over 129330 labelled bond
ends in the corpus:

| | |
|---|---:|
| average over-trim | **0.075 bond lengths** |
| over-trimmed by more than 0.1 | 28259 (**21.9%**) |
| worst | **0.39 bond lengths** — `[NH2+]`, nearly 40% of a bond of white space |

It now intersects the axis-aligned box along the bond direction (the label is
set `text-anchor="middle"`, so the box is centred on the atom). Bonds whose two
labels together exceed the bond length dropped from **2.77% to 1.26%** of
283604 bonds.

**Two judges are needed here; either one alone can be fooled:**

| Judge | Guards | How it is fooled alone |
|---|---|---|
| `a_bond_stops_at_the_label_box_not_at_its_circumscribed_circle` | tight enough | make `trim` cut nothing — it looks better, and the letters get struck through |
| `no_drawn_line_runs_across_an_atom_label` | far enough | go back to the circumscribed circle — still green |

Mutation-verified and orthogonal: the circle only reddens the first; cutting
nothing reddens both.

The remaining 1.26% genuinely do not fit — an ACS label is 0.69 bond lengths,
and `O⁻—N⁺` needs 1.375 of clearance. A flip cannot help (the bond length is
fixed), so it is not a violation; it is counted separately as **1818 pictures
(10.3%) with a label that does not fit on its bond**. Reducing it further needs
per-glyph boxes, the way RDKit's `StringRect` does it.

### Two atoms on one point invent a ring that is not there

When two atoms land on the same point their bonds meet end to end, and the
picture gains a ring the molecule does not have. A reader has no way to tell it
is fake. One triterpene drew a three-membered ring whose three sides were each
exactly one bond length.

This was **1064 / 17662 (6.0%)**, and every distance was *exactly* zero — the
layout walks unit steps on a 30° lattice, so two branches hitting the same
lattice point is systematic, not floating-point noise. All of the sampled cases
were reported in `unresolved`, so "say what you could not draw" held; but
"unresolved collision" badly understates *a ring that is not there*.

The fix is to check whether a position is already taken before placing a
substituent there, and step round by 30° if it is. An angle off its ideal is
merely ugly; it does not make anyone misread the structure. **89 cases remain**
(0.5%), all of them where five steps were not enough to find room.

### Crossing bonds: the first classification was wrong

381 crossings were first put down as "320 the flip operator can reach but does
not fix". That was wrong. Bucketed by what the two crossing bonds actually are:
**330 (86.6%) are a fused ring system self-intersecting**, and that set is
exactly the set where `degraded` is non-empty.

Collision relief cannot reach them **geometrically**: a fused ring system is
2-connected and ring bonds are excluded from the flip candidates, so two bonds
inside one system keep their relative positions under every flip. Enumerating
all 2^k reachable configurations confirms it — only **19 of 381** can be
cleared by flips at all.

`rings::relax` is local descent from a single starting guess. Giving it five
starts, all derived from canonical rank, and picking by (self-intersections,
worst bond-length deviation, rank-ordered quantised coordinates):

| | before | multi-start |
|---|---:|---:|
| pictures with a crossing | 381 | **281** |
| of those, ring-system self-intersection | 330 | 230 |
| writing-independence violations | 129 | **125** |

The start that does the work is "**lay the largest ring out as a regular
polygon**, then spread the rest outward from placed neighbours". The other four
all put every atom on one circle — too alike topologically, and the spring
descent falls into the same bad minima.

Two writing dependences were introduced and caught while doing it: picking the
anchor by `neighbors` storage order (violations 129 → 349), and using
`BTreeMap`'s index iteration order as the tie-break key. Both now go by
canonical rank.

### The 51 that remain on a clean layout

These layouts are clean and still cross. Reported in `Depiction::crossings`.

Two attempts, both near-useless, recorded so nobody repeats them:

- **Avoid already-drawn bonds when placing a substituent** — fixed exactly one
  case (382 → 381). By the time the BFS places it, the bond it will cross
  usually is not drawn yet. The check earns its keep elsewhere: atom
  coincidences dropped 89 → 74.
- **Rank crossings above collision depth in the score** — they already were;
  the code was misread. Running it the other way round raised crossings from
  381 to 415 with no reduction in collisions.

Clearing them needs the operator generalised from "flip one bond" to "take an
articulation point and one component hanging off it, and place that component
in one of 24 lattice poses". That subsumes terminal redirection, angle opening
and rotating a whole system about a spiro atom, and it changes exactly one bond
angle — the one at the articulation point. An out-of-tree experiment clears all
51, at the price of adding a level to the score (**coincidences ahead of
crossings**, or it trades a phantom ring for a crossing) and giving up "bond
angle equals its ideal" — which no audit property guards today. **Not
implemented.**

### A 60° angle is not merely ugly

Stepping round by 30° to dodge a taken spot can land two steps out, giving a 60°
kink — which reads as a three-membered ring that is not there. Measured on a
nitrogen mustard's `N—CH₂—CH₂—Cl` arms: 60.1°.

**Rejecting** the pinched directions outright was tried: narrow angles 334 →
227, at a cost of 494 more unresolved collisions and 2.8 points of clean rate.
The rejected direction is often the only one that does not collide. Not worth
it.

Reordering instead — at the same step size, try the side that widens the angle
first — improves everything at once with no trade: narrow angles 334 → **287**,
collisions 1189 → 1167, crossings 281 → 278, clean 91.3% → **91.4%**.

### Current standing

```
property              checked   violations
inside the canvas       17662            0
writing-independent     17660          257
  … full 5 rewrites     17405            0
  … never checked           2            0
no atoms coincident     17662           76
wedges readable           602            0
wedges drawn              602            0
ring double bonds       12731            0
no pinched angle        16155          287
bond lengths equal      16155            0
```

Writing-independence is **257 / 17662 (1.5%)**. That is higher than the 134
recorded earlier — **not because the pictures got worse, but because replacing
the shuffler exposed the half that had been invisible** (see above).

The attribution below was done on the old set of 45 pairs. It has **not been
redone** on the new 257 and is pending recheck. Fitting each failing pair with
its best rigid transform:

| Lines still differing | Pairs |
|---|---|
| 0 | 2 — pose only |
| **1** | **35 — one substituent points a different way** |
| 2–4 | 4 |

Nine cases in ten come down to a single bond; the canonical orientation step
then flips the whole picture to minimise its key, which is why the primitive
count reads "all different". The remaining root cause is localised to
substituent placement (ring systems themselves come out point-for-point
identical) and is documented in `harness/README.md`, along with the hypotheses
already ruled out.
