# Finite certificates for graph-lift sectors

If `2^m-1 > binomial(c,2)`, every c-point support has a zero-fibre direction.
The full-rank target support is then a graph over a full-rank parent of
affine dimension `m-1`. The graph functions form the kernel of the quadratic
evaluation matrix, modulo affine functions. The zero quotient element gives
a rank-deficient graph and is excluded. Thus each parent contributes exactly
`2^ell-1` graph-lift cosets, where

```text
ell = c - rank(E_2(parent)) - m.
```

All cosets, not just a count of them, can be enumerated independently with
binary elimination. This gives a compact finite proof of coverage when
paired with an explicit affine map from every lift to its catalogue class.
The maps are independently replayable without repeating affine search.

## Verified upper sectors at length 54

| Target affine dimension | Parent classes | Non-affine graph lifts | Target classes |
|---:|---:|---:|---:|
|13|21,585|34,816|938|
|14|938|942|47|
|15|47|23|2|
|16|2|0|0|
|17|0|0|0|

The 938 dimension-13 and 47 dimension-14 targets have distinct exact
invariant tuples consisting
of direction multiplicities, point-local difference profiles and affine
hyperplane intersection counts. Their inequivalence therefore needs no
negative backtracking decision. Every graph lift has an explicit positive
affine map into one of these classes. The two dimension-15 targets are
also independently inequivalent. Both have zero graph-lift quotient, closing
dimension 16; the empty predecessor then closes dimension 17.

These are complete induction steps, conditional on completeness of their
predecessor catalogues. They do not by themselves prove the predecessor
classification.

## Algorithms and reproduction

`code/verify_zero_fibre_spaces.py` is the dimension-generic reference
enumerator, exact quotient and witness checker. It works at multiple support
lengths, provided the stated zero-fibre inequality holds. Its Python affine
search is intended for small sectors or for checking supplied maps.

The accelerated length-54 route uses:

1. `code/prepare_zero_fibre_input.py` to expand any parent interval from the
   compressed catalogue and compute its quadratic-rank profile using the
   native source-profile kernel;
2. `mXX_fused_signature_kernel` to enumerate the lift quotient in Gray order
   and form affine invariant buckets;
3. `code/verify_separated_graph_cover.py` to enumerate the same affine-shear
   cosets independently, check their exact coverage, apply the native affine
   quotient and replay every positive map.

The third step rejects a sector if the target invariants do not separate
all classes. Such a sector requires the general exact quotient with negative
decisions as well; an invariant collision is never interpreted as equivalence.

For example, after building the native space kernels:

```text
python -B code/prepare_zero_fibre_input.py --dimension 14 --count 938 --native-directory build/spaces --work-directory work/m14 --output work/m14_input.json
build/spaces/m14_fused_signature_kernel --input work/m14/sources.bin --source-start 0 --source-count 938 --masks-output work/m14/candidates.bin --output work/m14/signatures.bin --threads 4
python -B code/verify_separated_graph_cover.py --dimension 14 --candidate-masks work/m14/candidates.bin --native-directory build/spaces --work-directory work/m14/quotient --output work/m14_cover.json
```

The adapter has the same interface in target dimensions 11 through 15.
Its neutral source-profile input wrapper contains no asserted affine
signature or class multiplicity; the rank kernel independently checks the
support's weight, full affine rank and moments.

The supplied dimension-13 and dimension-14 proofs contain the source
profiles and every affine map needed for a much shorter independent replay,
without native executables:

```text
python -B code/verify_separated_graph_cover.py --certificate certificates/zero_fibre_c54_m13.json --output work/m13_replay.json
python -B code/verify_separated_graph_cover.py --certificate certificates/zero_fibre_c54_m14.json --output work/m14_replay.json
python -B code/verify_zero_fibre_spaces.py --length 54 --dimension 15 --witnesses certificates/zero_fibre_c54_m15.json --output work/m15_replay.json
python -B code/verify_zero_fibre_spaces.py --length 54 --dimension 16 --output work/m16_replay.json
python -B code/verify_zero_fibre_spaces.py --length 54 --dimension 17 --output work/m17_replay.json
```

The checks bind the compact catalogue intervals and reject missing cosets,
repeated cosets, affine (rank-deficient) graphs and incorrect affine maps.

## Affine search refinement

The dimension-13 affine kernel refines the difference-labelled complete
graph on the support points. Initially a point is labelled by the multiset
of difference multiplicities to the other points. Each refinement replaces
its label by its preceding label together with the multiset of pairs
`(label of other point, difference multiplicity)`.

Every affine equivalence preserves the initial labels and, by induction,
each refined label. Requiring equal labels during affine backtracking
therefore cannot remove a valid map. Choosing an origin in a smallest label
class changes only the search order. The refined labels are not assumed to
be complete orbit invariants; all positive decisions still supply a full
invertible affine map. `verify_point_colour_refinement.py` compares the
native labels with an independent reference and checks affine covariance.

## Spanning-codeword graph backend

An alternative exact backend is `code/native/auxiliary/affine_codeword_graph.cpp`.
Let H be the binary row space generated by the unit and coordinate rows of
a support. Select all nonzero words of H of weight at most t, where t is
the smallest threshold at which those words span H. Form a bipartite graph
with one colour for coordinate vertices and a distinct colour for codewords
of each selected weight; join a word to the coordinates in its support.

This graph is a complete coordinate-equivalence representation of H.
Forward, a coordinate equivalence preserves weights, spanning ranks and
the minimal threshold, and therefore induces a colour-preserving graph
isomorphism. Conversely, a graph isomorphism induces a coordinate permutation
mapping the selected words, and hence their full row span H, to the other
row space. Since coordinate permutations fix the unit vector, this is
exactly affine equivalence of the unital supports.

Bliss canonicalises the graph. The resulting coordinate order determines
an affine frame and a canonical support. Every returned frame is checked
directly, so positive equivalences carry independently replayable witnesses.
The construction works for varying lengths and affine dimensions, with the
implemented fixed-width limits `n <= 63` and `m <= 17`.

`verify_affine_codeword_graph.py` compares the backend with all 1,344 affine
maps on all 149 full-affine-rank subsets of F_2^3, and with affine-transformed
length-54 examples in dimensions 13 through 15. It also independently
enumerates the codewords to check the chosen spanning threshold. The graph
backend is selected by adding `--graph-backend build/affine_codeword_graph`
to `verify_separated_graph_cover.py`. It produces the same finite affine-map
proof format and uses the same independent witness replay.
