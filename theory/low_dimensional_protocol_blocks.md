# Lower-Dimensional Finite Protocol Blocks

The source partition in
`data/protocol_sectors/length54_low_dimensional_blocks.json` contains
21,077 blocks. Together they cover every length-54 source of affine
dimensions eight through ten, and auxiliary length-52 sources of affine
dimensions seven through ten. There are 300,071,794 source entries in
these blocks. The higher-dimensional source partition is described in
`high_dimensional_protocol_sector.md`.

The auxiliary length-52 entries are not needed to produce a new Pareto
point at length 53: both possible constructions append a zero column,
and are therefore dominated by a length-52 protocol with the same output,
exact distance and row footprint. See the zero-column lemma in
`support_correspondence.md`. The finite block cover includes 186 of the
188 length-52 cubic classes; the lemma applies to the entire length-52
parent sector, including the other two entries.

## Pointing and Distance Profiles

All origins are considered, with exact duplicate pointed sets removed.
For a length-54 source in affine dimension `m`, this gives exactly
`2^m+55` pointed supports at protocol lengths 53 and 54: all translated
even supports, the 54 choices of removed point, and one affine-hyperplane
embedding. Such a source has no nonzero translation period. Otherwise,
choosing a linear functional nonzero on a period vector would give a
stabiliser row of weight 27, contradicting even overlap.

A length-52 source contributes only the zero-column extensions at length
53. If its translation group is trivial, there are `2^m-51` distinct
pointed supports; if the translation group has dimension one, there are
`(2^m-52)/2+1`. Dimension two is impossible: two linear coordinates on a
period plane would have overlap `52/4=13`, which is odd. The block counts
are consistent with exactly 21,500 sources having a one-dimensional
translation period.

For every block, the complete histograms of distance-three and
distance-four logical-label quotient dimensions sum to the pointing
count. There are 172,729,013,588 pointed supports in these finite domains.
The q1 through q4 branches enumerate every common isotropic logical space
in the distance-four quotient. Quotient dimensions above the raw threshold
are handled by the exact marked-code quotient; its orbit counts are not
misrepresented as counts of all raw subspaces.

The q5 branch enumerates the distance-three quotient dimensions five
through fourteen using raw subspaces. Dimensions below five cannot contain
a q5 logical space. There are exactly 47,848 pointed supports of quotient
dimension at least fifteen, which form a separate finite enumeration
obligation. They are not included in the small-quotient completeness claim.

## Witness and Successor Checks

`verify_low_dimensional_blocks.py` checks the source intervals, pointing
counts, both quotient histograms, every logical branch domain, and its
binding to the complete q5 output profiles. It independently checks all
7,859 recorded finite-block Pareto witnesses, including their output
equivalences, exact distances, error coefficients and frontier dominators.
Every comparison fixes exact distance. The q1 through q4 exact-distance-
three omissions are separately justified by `verify_low_q_pruning.py`.

The q5 output profiles feed the necessary-profile and refinement tests in
`output_profile_pruning.md`. Those tests cover successors within the
small-quotient domain. The larger-quotient sector and the complete shorter-
length classification are independent obligations in the global proof.

```text
python code/verify_low_dimensional_blocks.py --output build/low_dimensional_blocks.json
```

This command rechecks domain consistency and all supplied matrix witnesses.
It does not claim to have freshly repeated the entire subspace census.
The accelerated mathematical primitives for reproduction are
`manifest-origin-reduce`, `manifest-initial-label-caches`, and
`manifest-catalogue` in the native kernel. The paired label-space operation
can use its `--source-schur` option to share the linear-algebra work between
pointings of the same source.
