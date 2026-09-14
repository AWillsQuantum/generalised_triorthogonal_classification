# Protocol enumeration from compact supports

The entry points below operate on mathematical support records, without
external data paths or task-scheduling dependencies. Build `utsp-native`
using `build.py`. Python readers require the packages in the release's
`requirements.txt`.

## A small complete example

From the release root, with the executable at `build/utsp-native`:

```text
python code/compact_to_native.py data/spaces/c16/m04/spaces_000000000.utspace.zst --output work/base.utsp
build/utsp-native manifest-origin-reduce --input work/base.utsp --output work/pointed.utsp --maximum-protocol-length 54 --workers 1
build/utsp-native manifest-catalogue --input work/pointed.utsp --output work/q1.json --q 1 --marked-orbits --minimum-distance 3 --workers 1
```

The first command writes the affine support, not an already chosen protocol.
The second enumerates the five pointing cases in
`theory/support_correspondence.md`, reducing origins by exact support
automorphisms. The last enumerates logical subspaces and computes their
outputs, exact distances and Pareto witnesses. This example recovers the
15-input one-T protocol at distance three and footprint five.

`compact_to_native.py --start I --count K` selects a bounded interval within
one compact shard. The source identity includes `(c,m,index)`. Splitting
an input interval does not change the underlying mathematical problem.
All source intervals must be covered for a complete classification.

## Complete logical-subspace enumeration

For a fixed input manifest and `q=1,...,8`, `manifest-catalogue
--marked-orbits` uses `enumerate_complete_nondegenerate_subspace_orbits`.
It combines the hereditary nondegenerate branch with the exceptional
primitive seeds in dimensions 2, 3, 5 and 7. It does not assume every
nondegenerate tensor has a nondegenerate hyperplane. It deduplicates the
branch union by exact marked-code canonical keys.

`--final-nondegenerate-orbits` instead keeps all intermediate tensor types
and imposes nondegeneracy only in the final dimension. It is an independent
enumeration route, useful for cross-checks; it can be substantially slower.
The unflagged route enumerates raw isotropic subspaces, and is intended for
small validation cases. These routes cover the same intrinsic outputs.

The q5 target-chain and q4 hitting-set variants are alternative complete
finite covers when combined with their primitive branches. Their validity
depends on the eligible-parent rule proved in
`theory/canonical_augmentation.md`. Intermediate Pareto pruning is not
permitted. The fresh finite census and q3-chain checks can be reproduced by:

```text
python code/verify_tensor_authorities.py --native build/utsp-native
```

## Distance and the finite logical cutoff

`--minimum-distance D` imposes all lower-weight zero-sum label relations
before enumeration. The output records still distinguish exact distances;
they do not merge different distances into a threshold objective. A support
with quotient dimension below `q` is excluded by linear algebra.

Enumerating q through eight is enough in the present range only after the
numerical absence of all q8 protocols and both primitive q7 types is
established on every relevant support. The arbitrary-dimensional theorem
in `theory/logical_dimension_closure.md` then closes higher dimensions.
The finite tensor certificates alone do not establish those protocol
absence claims.

Full reproduction also needs the certified support lists and the global
distance/dominance partition. A successful run on a single interval is
evidence for that interval, not for the whole catalogue. These interfaces
are usable now; the release status records whether the complete coverage
evidence has been assembled.
