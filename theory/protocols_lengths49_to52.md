# Protocol sectors at lengths 49 to 52

The source spaces, logical quotients and finite subspace censuses are
separate objects. A source catalogue interval identifies unital spaces;
its pointed supports identify stabiliser spaces; a common-isotropic
subspace of a distance-filtered logical quotient specifies logical rows.
The classification concerns the full projective matrix scope, including
matrices without a completion to the restricted puncturing construction.

The relevant mathematical reductions are in `support_correspondence.md` and the
other theory files linked from the main README. This file specifies the
finite domains at the two intermediate cutoffs. A sector verification
is conditional on the completeness of its parent-space catalogue and
on any explicitly separated logical sectors. It is not by itself a
certificate for the entire length-54 classification.

## Exact-distance and Pareto pruning

For distance at least three, distinct full columns have distinct
stabiliser syndromes: otherwise a weight-two vector would have zero
stabiliser syndrome and nonzero logical syndrome. Thus an h-dimensional
stabiliser has at most 2^h columns. At n >= 49 this gives h >= 6.

For each of the 22 intrinsic outputs on at most four qubits, the release
contains a distance-three protocol with n < 49 and S <= q+6. Consequently,
every new exact-distance-three candidate with q <= 4 is Pareto dominated.
`verify_low_q_pruning.py` checks both the complete set of these outputs
and their incumbents. This reduction does not exclude a higher exact
distance: those candidates are enumerated with distance floor four.

Logical branches with q >= 5 use distance floor three. Every intrinsic
candidate contributing a finite Pareto point has its exact distance and
leading error coefficient evaluated, not merely a lower bound on distance.

A zero stabiliser syndrome must also have zero logical entries at distance
at least three. Removing such a zero column preserves the output, row
rank and exact logical distance, and reduces n. Hence the extensions of
shorter even parent spaces by a zero column do not add Pareto points at
the next odd length. The new source lengths here are 50 and 52.

## Parent length 50

`data/protocol_sectors/length50_source_domain.json` partitions all 304,611
parent classes into 309 disjoint catalogue intervals. The source counts
by affine dimension m are:

| m | Classes |
|---|--------:|
| 8 | 95,657 |
| 9 | 180,054 |
| 10 | 26,288 |
| 11 | 2,415 |
| 12 | 186 |
| 13 | 11 |

The complete origin-reduced domain has 82,266,913 pointed supports of
protocol lengths 49 and 50. Each block records the full distance-three
and distance-four quotient histograms. The branches are q=1,...,4 with
distance floor four and q=5,...,8 with distance floor three. A shared
distance-filtered quotient is used without changing the enumerated
mathematical domain.

`enumerate_source_block.py` reconstructs any block from the compressed
space catalogue and calculates all eight raw common-isotropic censuses.
`verify_source_block_sector.py` verifies the interval cover, eligibility
counts, exact subspace and radical masses, all finite witness parameters,
and the union with the preceding frontier. A complete zero common-isotropic
q=8 census excludes every larger logical dimension by subspace containment.
The complete cumulative finite metrics are required to equal the 53
released frontier points of length at most 50.

All 309 blocks have passed a fresh complete replay. The aggregate
verification in `certificates/through50_protocol_sector.json` confirms
82,266,913 pointed supports, all exact census totals, the absence of
common-isotropic q=8 subspaces, 19 independently checked local witnesses,
and exact agreement with those 53 cumulative frontier points.

## Parent length 52

`data/protocol_sectors/length52_source_domain.json` identifies all 7,519,688
parent classes. Of these, 7,519,686 are in 7,711 main intervals; the two
remaining cubic sources have indices 186 and 187 in the m=7 catalogue.
All catalogue indices are zero based. The full main source data have been
compared record by record with their compressed supports.

The main intervals divide into:

- 7,536 blocks containing 7,519,511 sources, with explicit small-quotient
  censuses in `length52_small_quotient_blocks.json`.
- 175 individual cubic sources with a separately specified pointed domain.
  Their 8,923 origin representatives are in
  `length52_selected_cubic_domain.json` and its 14,261-byte compressed
  pointing file. A fresh complete origin-orbit calculation reproduced
  every pointed set exactly; its certificate is
  `certificates/length52_selected_cubic_origins.json`.

The first part contains 4,163,769,607 pointed supports and 69,741 finite
branch or theoretical-exclusion records. Its complete raw q=5 census for
quotient dimensions 5 through 14 contains 315,663,582,926 isotropic
subspaces, including 365,738 intrinsic ones. Every one of its 5,114
positive supports is bound to the complete output profile in
`data/protocol_profiles/length52_small_quotient.json`. The 5,068 pointed
supports of larger quotient dimension have the separate complete finite
cover described below; they are not implicitly declared empty.

Each finite branch specifies a quotient-dimension interval and an interval
of the selected pointed input list. For each relevant quotient dimension,
the input intervals must form an exact disjoint cover. Their quotient
histograms must sum to the complete selected-input histogram. This checks
coverage even when a quotient stratum has several finite pieces.

For raw origins, translation periods must be included in the cardinality
check. A length-52 unital support has translation period dimension p=0 or
1. Indeed, a period of dimension at least three forces divisibility of the
length by eight. A period of dimension two would give an odd number of
four-point cosets, contradicting the quadratic moment formed by dual
coordinates on that period. For a source of affine dimension m, the number
of distinct pointed supports at lengths 51 and 52 is therefore

    (2^m + 52) / 2^p + 1.

The summands count translations of the whole support, punctures at a
support point, and the affine-hyperplane embedding, respectively.
The complete raw source cover has 21,491 sources with p=1. Both the total
and the three nonzero parity-case counts are checked independently.

## Constant-stabiliser cubic cases

Eleven cubic sources have a nonzero distance-four quotient only for the
affine-hyperplane embedding. Independent calculation of all 1,991 pointed
supports confirms this, and confirms that each positive case contains the
constant stabiliser row.

If the stabiliser contains 1, the mixed-overlap equations imply that every
logical row is even and every pair of logical rows has even overlap.
The logical trilinear tensor is consequently alternating. In dimensions
one and two it is zero. Every alternating trilinear form in dimension four
has a nonzero radical: after choosing a volume form, it is the contraction
of that volume form by a vector, which lies in its radical. Thus none of
these three logical dimensions has an intrinsic output.

Dimension three is different: its nonzero alternating tensor is CCZ.
The complete direct q=3 subspace censuses are supplied in
`certificates/length52_cubic_blocks/`. Their intrinsic candidates all have
n=52, S=11 and exact distance four. They are dominated by the released
n=48, S=10, distance-four CCZ protocol. This is a Pareto exclusion, not
a claim that the intrinsic subspaces are absent.

`verify_cubic_distance_four.py` checks the constant-row premise, the
complete quotient distributions, the raw q=3 censuses, and the exact-distance
dominators. It also independently enumerates the 16 alternating tensors in
dimension four. `verify_length52_small_blocks.py` includes this check in
its finite source-block verification.

## Selected cubic logical sectors

The 175 selected cubic sources have 8,923 pointed supports. Their complete
five-dimensional output profiles and the 274 six-dimensional profiles
selected by the necessary-profile theorem are recorded in
`length52_selected_cubic_censuses.json`. The respective weighted numbers
of intrinsic subspaces are 113,524,604 and 130,894. These are marked-orbit
censuses with exact orbit weights, not counts of all raw isotropic subspaces.

Independent linear calculations give a nonzero distance-four quotient
only on 175 constant-stabiliser supports. Every distance-five quotient is
zero. The alternating-tensor argument above therefore closes q=1,2,4,
and every possible q=3 distance-four metric is dominated by the shorter
CCZ incumbent. Each positive q=5 or q=6 output profile belongs to a support
whose distance-four quotient is zero, so its exact distance is three.
Its complete output keys then determine all possible metrics: n is the
support length and S=h+q.

Every normalised support is explicitly bound to the primitive-target
domain. Its two primitive q=5 targets are absent. The complete q=6
profiles are contained in the five profiles covered by the finite
seven-dimensional anchored exclusion. The primitive restriction theorem
and descent through intrinsic hyperplanes exclude every q>=7.
`verify_selected_cubic_sector.py` checks these implications, the complete
profile-to-input correspondence and every metric comparison. A fresh
positive-case q=5 chain replay is also supplied. The verifier does not
describe the retained complete profiles as a fresh full-domain replay.

## Larger quotients

`length52_large_quotient_domain.json` specifies all 5,068 selected pointings
from 70 parent classes. Zero-column deletion and exact linear equivalence
reduce them to 441 stabiliser geometries. All 5,068 canonicalisations have
been independently repeated, and their per-block cardinalities agree with
the complete quotient histograms of the small-block cover.

There are 179 positive q=5 profiles and 654,868 weighted intrinsic
subspaces. For every positive profile, the distance-four quotient has
dimension less than five. This is enough to establish exact distance
three for every q=5 candidate; the quotient need not be zero. In particular,
case 385 has a one-dimensional distance-four quotient and a positive q=5
profile. A fresh chain replay reproduces its four marked orbits and
2,048 weighted intrinsic subspaces.

The primitive-target zero censuses cover all 441 geometries: 393 use the
shared primitive domain and 48 have separate complete target results.
Only profiles 61 and 109 pass the q=6 necessary-profile test. Both have
fresh complete q=6 censuses, explicitly bound to these same stabiliser
geometries. Their profiles satisfy the anchored q=7 exclusion; the
primitive restriction theorem closes every larger dimension. All q=5
and q=6 metrics are covered by the released frontier. The q<=4
higher-distance cases belong to the already bound full source-block
censuses, without imposing the q=5 small-quotient cutoff on them.

`verify_length52_large_quotients.py` verifies this composition. Its
reproduction companion, `replay_length52_large_quotient.py`, handles
either the primitive target census or the nonprimitive chain census
on any one of the 441 explicit inputs.

## Two separate cubic parents

The m=7 parent classes with indices 186 and 187 have a combined 362 raw
pointings at protocol lengths 51 and 52. Their finite domain contains
58 nonzero-syndrome h=7 representatives and two h=8 affine-hyperplane
representatives. Exact canonicalisation of every raw pointing, after
zero-column deletion, proves that these 60 inputs cover the full domain.

The h=7 cases have zero distance-four quotient. Their complete raw q=5
censuses contain 90,201,752 isotropic subspaces and 26,928 intrinsic ones;
51 of the 58 profiles are positive. One profile passes the q=6 test. Its
fresh q=6 census is the case on parent 186 with origin 104, and it is bound
by exact support equivalence. The same tensor-profile and primitive
restriction arguments exclude all q>=7.

Both h=8 inputs have equal distance-three and distance-four label spaces,
a constant stabiliser row, and zero distance-five quotient. For parent
186, the complete q=3 census contains 497,739 isotropic subspaces, 252 of
which are intrinsic CCZ outputs. Their n=52, S=11, distance-four metric
is dominated. Its q=5 census has 3,439 subspaces, all with zero tensor.
For parent 187, all 1,327,299 isotropic three-spaces have zero tensor;
an alternating tensor on any larger isotropic subspace is consequently
zero as well. All three higher-distance censuses have been freshly
repeated. These two parents add no Pareto metrics.

The data are in `length52_separate_cubic_censuses.json`; the full input,
distance, output and successor checks are reproduced by
`verify_cubic_parent_sector.py`. Its `--replay-case` option also repeats
a chosen retained q=5 raw census.

## Reproduction interfaces

After building the native mathematical executable as described in the
README, examples of finite-domain reproduction are:

```sh
python code/enumerate_source_block.py --block 0 --native build/utsp-native --work-directory work/c50 --output work/c50_block0.json
python code/replay_length52_block.py --block 604 --q 3 --native build/utsp-native --work-directory work/c52 --output work/c52_block604_q3.json
python code/verify_length52_small_blocks.py --output work/c52_small_quotient_check.json
python code/replay_selected_cubic_profile.py --q 5 --pointing 5653 --native build/utsp-native --work-directory work/c52_selected --output work/c52_selected_replay.json
python code/verify_length52_large_quotients.py --native build/utsp-native --work-directory work/c52_large --output work/c52_large_check.json
python code/verify_cubic_parent_sector.py --native build/utsp-native --work-directory work/c52_cubic --output work/c52_cubic_check.json
```

The length-52 reproduction interface constructs raw origins for affine
dimension at most ten and affine-automorphism orbit representatives above
that. It shares the two distance-filtered label spaces. The q=5 default
quotient cutoff is 14; increasing the cutoff is an explicit change to the
finite domain, not part of the small-quotient completeness claim.

The sector certificates, the preceding cutoff and the small-quotient
output-profile exclusions have passed the combined check in
`certificates/through52_protocol_cover.json`, giving exactly 73 frontier
points. `protocol_cover_composition.md` states the complete inference and
its finite-census premises. Independent parent-space completeness remains
a separate obligation in `AUDIT_STATUS.json`.
