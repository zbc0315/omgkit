# Batches

Molecules have two representations.

| Type | Mutability | Layout | For |
|---|---|---|---|
| `MolBuilder` | mutable | array of structs, one molecule | building: parsing, constructing reaction products |
| `MolBatch` | immutable | struct of arrays + CSR, across molecules | algorithms, parallelism, zero-copy export |

They convert both ways through `MolBatchBuilder::push` and
`MolView::to_builder`, and the round trip is identity — checked by test.

!!! warning "Python has `Mol` only"

    The Python `Mol` wraps `MolBuilder`. `MolBatch` is not wrapped yet; it is
    reachable from Rust today. See
    [Python API](../api/python.md#not-yet-exposed).

## Why columnar

The usual molecule graph is an adjacency list plus heap-allocated atom objects,
so vertex properties live behind pointers — walking the atom properties of a
40-atom molecule is 40 random accesses.

`MolBatch` puts every atom's value for one property into a single contiguous
array, with adjacency in CSR form. Three things follow directly from that
layout:

1. **Sequential access** — walking a property is one contiguous scan, not a
   pointer chase
2. **Zero copy** — a column can be exposed as a numpy or Arrow buffer as is
3. **Cheap parallelism** — split on molecule offsets, no shared state

These are consequences of the layout and hold without measuring anything. How
much end-to-end speed the layout buys is a separate question that needs
benchmarks, and is not claimed here.

## Index convention

Inside a batch, indices are **global** — continuous across molecules. That is
the natural form for CSR. User-facing **local** indices (0-based within a
molecule) come from `MolView`, which holds a reference to the batch plus a
molecule number and copies nothing.
