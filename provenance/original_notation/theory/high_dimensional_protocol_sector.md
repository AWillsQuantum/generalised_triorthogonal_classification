# High-Dimensional Stabiliser Supports

For protocol lengths 53 and 54 the even support correspondence uses only
unital supports of lengths 52 and 54. Separate the source spaces by their
affine dimension `m`. This note concerns the complete sector `m>=11`,
which contains 620,475 source classes. It is a source-dimension condition,
not a condition on the protocol's logical dimension `q`.

The stabiliser dimension is `h=m` or `m+1`, according to the pointing case.
All five cases are included. The 627 blocks in
`data/protocol_sectors/length54_high_dimensional.json` partition the source
indices without overlap or gaps. Each source interval is bound to the
ordered support digest of its compact catalogue, independently of the
chosen spelling of a representative matrix.

For each pointed support, let `W_d` be the logical-label quotient with all
zero-syndrome errors of weight less than `d` annihilated. The data record
the dimension distributions of `W_3` and `W_4`. A `q`-dimensional logical
space is enumerated precisely when it is common totally isotropic for the
mixed-overlap forms on `W_d`.

The finite census uses `d=4` for `q=1,...,4` and `d=3` for `q=5,...,8`.
These thresholds suffice for the new protocol lengths: all exact-distance-
three candidates with `q<=4` are strictly dominated, as established in
`output_profile_pruning.md`. Candidates with larger exact distance are not
discarded by that comparison.

The census has no nondegenerate outputs for `q=2,...,6` at the stated
thresholds. In fact there is no common isotropic subspace of dimension
seven, even with degenerate output. Since every higher-dimensional common
isotropic space contains a seven-dimensional one, this excludes every
`q>=7` directly; no assumption about eligible nondegenerate parents is
required for this sector.

The one-dimensional branches produce 464 finite-block Pareto witnesses.
For each, the verifier checks full rank, projectivity, all mixed overlaps,
the nonzero one-dimensional output tensor, exact distance and its leading
error coefficient. Each witness is weakly dominated by a distributed
frontier point at the same exact distance. Weak dominance permits equality
of both objectives: the purpose here is to show that no new metric point
is needed. It is not a claim that the finite candidate matrices are all
equivalent to the frontier matrices.

## Independent Checks and Reproduction

Run the finite-cover and witness checks with:

```text
python code/verify_high_dimensional_sector.py --output build/high_dimensional_sector.json
```

This recomputes the interval checks and all matrix parameters. It checks
the consistency of the supplied exhaustive enumeration counts; it does
not itself re-enumerate every logical subspace.

To re-enumerate a block directly from the compact source catalogue:

```text
python code/replay_high_dimensional_block.py --block 0 --native build/native/utsp-native --work-directory build/high_dimension_0 --output build/high_dimension_0.json
```

The function constructs the source manifest, computes the pointing cover,
enumerates all eight branches and compares the exact counts and Pareto
metric sets. The same function accepts every block index from 0 to 626.
Different matrix witnesses at an identical metric point are permitted;
their validity and exact distance are checked independently.

The complete protocol theorem additionally requires all lower-dimensional
source sectors and all shorter lengths. This finite sector is not a
standalone certificate of the full classification.
