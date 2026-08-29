# 3D structures

`Mol.conformer()` turns a molecule into one three-dimensional structure. It is
meant as the **starting point a force field refines from** — not a conformer
ensemble, not a search, and not something that should be compared against an
experimental geometry.

```pycon
>>> import omgkit
>>> conf = omgkit.parse_smiles("C[C@H](N)C(=O)O").conformer()
>>> conf
<omgkit.Conformer atoms=13 energy=0.000e0 converged=True chiral=1/1>
>>> conf.coords[0]
(-1.1906..., -0.8985..., -0.0732...)
```

## It is deterministic

**There is no random seed, and there is no retry loop.** The same molecule
always comes back with the same coordinates:

```pycon
>>> omgkit.parse_smiles("CCO").conformer().coords == omgkit.parse_smiles("CCO").conformer().coords
True
```

That is not a convenience — it is where the difference in failure rate comes
from. The usual pipeline is:

```text
bounds matrix → triangle smoothing → sample one distance per pair from its
interval → metric-matrix embedding → refine → if it fails, re-sample and repeat
```

The problem is the sampling step. Each pair is drawn **independently**, so the
resulting table often describes distances no 3D arrangement can satisfy — "A is
3 m from B, B is 3 m from C, A is 10 m from C" is easy to write down and
impossible to draw. The standard response is to discard the whole attempt and
re-draw, up to `10×N` times. When the cause is structural rather than unlucky,
**all `10×N` retries fail the same way**.

omgkit keeps the embedding and replaces only that step. After triangle
smoothing, the upper-bound matrix `U` satisfies the triangle inequality by
construction — it *is* a metric, a distance table that can be drawn — so it is
used directly as the reference table. On the same batch of bounds matrices:

| | spectrum: negative share | spectrum: top-3 share | geometry: out-of-bounds RMS |
|---|---:|---:|---:|
| sampled from the interval | 0.284 | 0.463 | 0.812 Å |
| **`U` used directly** | **0.042** | **0.889** | **0.323 Å** |

## Failure rate

On a 8831-molecule drug-like corpus:

| | failures | |
|---|---:|---|
| RDKit ETKDGv3 2025.09.2 | 36 (0.41%) | mostly metal complexes |
| **omgkit** | **1 (0.01%)** | the one case has contradictory distance bounds |

**Only contradictory bounds count as failure** — the case where not even a
self-consistent distance table exists. An embedding that will not fit into
three dimensions is not a failure (the extra dimension is collapsed and
refinement recovers it), and refinement that does not converge is not a failure
either (you get the best coordinates so far, plus an honest residual). The
product of this step is a *starting point*: one that is slightly off can be
fixed, and having none cannot.

### What "not a failure" costs

Those two non-failures are real, and they have numbers. On the same corpus, out
of 8830 structures produced:

| | | |
|---|---:|---|
| refinement did not converge | 443 (5.02%) | you get the best coordinates so far, and `converged` is `False` |
| at least one bond outside its bounds by >0.1 Å | 13 (0.15%) | mostly thiocyanate-bridged metal complexes |
| 1-2 bond distances outside bounds | 19 / 256955 pairs (0.007%) | |
| 1-3 angle distances outside bounds | 584 / 434872 (0.134%) | |
| bond crossings | 0 | |

`converged` and `energy` are on `Conformer` for exactly this reason: a structure
that did not converge is still worth handing to a force field, but you should be
able to tell.

### Metal complexes

The metal-complex gap comes from the bounds table rather than from the solver.
For a centre with more than four ligands, RDKit writes the 1-3 bound as
`[1.0, 1.2×(b₁+b₂)]` — the source comment calls the lower value "an arbitrary
min angle", and the upper one is wider than 180°. That leaves six ligands with
effectively no angular constraint, the metal lands on the centroid, and the
embedding check then rejects it. omgkit uses an angle envelope indexed by
coordination number, which is the same code path for every coordination number.

## What comes back

| | |
|---|---|
| `coords` | `(x, y, z)` per atom, in Å |
| `mol` | **the molecule the coordinates belong to** — see the warning below |
| `energy` | the error function after refinement; `0` means every distance landed inside its bounds |
| `energy_before` | the same before refinement, so the two together say how much work refinement did |
| `converged` | whether the gradient reached the tolerance |
| `iterations` | how many refinement steps it took |
| `chiral_total` / `chiral_ok` | how many stereocentres there are, and how many came out with the right sign |

!!! warning "The coordinates belong to `conf.mol`, not to the molecule you called it on"

    Generation needs explicit hydrogens, so it works on a **copy** and adds them
    there. Your molecule is left untouched, and `conf.coords` lines up with
    `conf.mol`, which has more atoms.

    ```pycon
    >>> m = omgkit.parse_smiles("C[C@H](N)C(=O)O")
    >>> conf = m.conformer()
    >>> m.num_atoms, conf.mol.num_atoms
    (6, 13)
    ```

**`chiral_ok` should always equal `chiral_total`.** Anything else means a centre
came out as its enantiomer, and it is worth checking rather than assuming:

```pycon
>>> conf.chiral_total, conf.chiral_ok
(1, 1)
```

Global chirality is fixed by one discrete step, not by refinement. Flipping the
handedness of a structure is a **reflection**, which is not in the connected
component of `SO(3)`: a continuous descent from a structure to its mirror image
would have to pass through a completely flat molecule, and no penalty weight
makes a descent method pay that. So both mirror images are scored once and the
better one is kept.

## Writing it out

`Conformer.to_molblock()` gives the contents of a `.mol` file. An `.sdf` is
those records separated by `$$$$`:

```python
with open("out.sdf", "w") as f:
    for smi in ["CCO", "C[C@H](N)C(=O)O"]:
        conf = omgkit.parse_smiles(smi).conformer()
        f.write(conf.to_molblock(title=smi))
        f.write("$$$$\n")
```

```text
L-alanine
  omgkit

 13 12  0  0  0  0  0  0  0  0999 V2000
   -1.1906   -0.8985   -0.0732 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2162    0.1825    0.3404 C   0  0  0  0  0  4  0  0  0  0  0  0
```

Aromatic bonds are kekulized on the way out — a molblock has no aromatic bond
type. The second line is the program name with **no timestamp**, so writing the
same molecule twice gives byte-identical output.

Reading such a file back gives the stereochemistry back too; see
[Reading and writing `.mol` and `.sdf` files](molfiles.md).

## Errors

`ValueError` is raised when the molecule cannot be sanitized, or when its
distance bounds contradict each other. The message says which, and for
contradictory bounds it names the pair of atoms whose interval came out empty.

## What this is not

- **Not a conformer ensemble.** One molecule in, one structure out.
- **Not a search.** No torsion scanning, no clustering, no energy window.
- **Not a force-field-minimised geometry.** The error function here is distance
  and chirality violation, not MMFF or UFF energy. Bond lengths and angles are
  reasonable and the stereochemistry is right; everything else is what a real
  force field is for.

## The Rust API

```rust
use omgkit_conf::pipeline::conformer_for;

let mut m = omgkit_io::smiles::parse("C[C@H](N)C(=O)O")?;
let conf = conformer_for(&mut m)?;   // m is sanitized and given explicit Hs in place
assert_eq!(conf.coords.len(), m.num_atoms());
```

`conformer_for` takes `&mut MolBuilder` and modifies it — that is where the
explicit hydrogens go. The Python binding copies first, which is why
`conf.mol` is a separate molecule there.
