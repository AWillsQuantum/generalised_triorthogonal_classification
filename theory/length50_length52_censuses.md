# Finite space censuses at lengths 50 and 52

The recurrence and equivalence are those of `space_classification.md`.
All spaces, not only indecomposable ones, occur in the source and target
domains. The following are affine dimensions m, so each unital generator
matrix has m+1 rows.

| Length | m | Complete candidate domain after the stated reductions | Affine classes |
|---:|---:|---:|---:|
|50|8|6,663,276|95,657|
|50|9|63,532,224|180,054|
|50|10|3,532,589|26,288|
|50|11|70,027|2,415|
|50|12|3,596|186|
|50|13|127|11|
|50|14|0|0|
|52|8|145,353,817|1,321,156|
|52|9|126,307,993|5,601,443|
|52|10|73,689,787|551,228|
|52|11|2,681,842|42,404|
|52|12|160,954|3,063|
|52|13|5,093|199|
|52|14|84|7|
|52|15,16|0|0|

There are also 188 length-52 classes at m=7, from the complete cubic
base. The total counts are 304,611 and 7,519,688 respectively.

## Domain and exact quotient

At m=8, pair averaging and the complete base leave core lengths 44 and
48 for target 50, and 44,48,52 for target 52. The N<=3 full-span lemma
excludes proper-span cores on those routes. For N>=4 the only potential
proper-span core lengths are 42 and 44, and the complete base has no such
core of dimension below seven.

The 131 primary length-50 cores and two alternative cores partition the
entire domain. At length 52 there are 317 primary and four alternative
cores. The alternative-direction witnesses are in
`alternative_contraction_cover.md`. The complete raw marked pair counts,
including alternative families, are 446,234,624 and 9,332,425,028.
The stabiliser action quotients the label/fibre choices; a subsequent
full affine quotient forgets the marked direction.

At m=9, all N=0,1,2 sources are included. For length 52 only lifts whose
chosen direction attains the minimum pair multiplicity are retained.
Every support has such a direction, so this does not remove an affine
class. The complete full-rank domain before this filter has 1,882,307,772
lifts. At length 50 there is no minimum-direction filter. At m=10 both
N=0,1 routes are retained. For m>=11 the pair bound forces N=0, reducing
the problem to all nonaffine solutions of the graph-lift equations.

Signatures only partition the candidate domain; equal signatures are
resolved by exact affine transporters. The retained censuses specify
the bucket and class distributions, complete interval counts, positive
affine-witness checks and negative exact-search checks. They bind the
representatives to the compact catalogue. Individual assignment streams
for these sectors are not included: independently checking every such
decision requires regeneration, not merely inspecting the census totals.

`code/verify_preceding_space_censuses.py` checks this composition,
reconstructs every marked-core profile from the complete affine base,
and checks the input partitions and aggregate quotient accounting.
Its certificate does not claim a fresh full exact-quotient computation.

## Accelerated reproduction

Build the common native algorithms and a required configuration, for
example length 52 and m=10:

```text
python -B code/space_native/build.py --build-dir build/spaces
python -B code/space_native/generic/build.py --length 52 --dimension 10 --include-lifts --build-dir build/spaces
```

The same second command supports lengths 50,52,54 and dimensions 9,10,11.
The configured routines retain the binary elimination, Gray-code lift
traversal, exact minimum-direction filter and fused signature algorithms.
The unfiltered length-50 m=9 configuration omits the unnecessary minimum
calculation. `certificates/preceding_space_native_lifts.json` compares 24
complete native families with independent binary equations modulo affine
shears, and compares fused signatures with a separate signature kernel.

Generate a bounded source interval directly from the compact catalogue:

```text
python -B code/enumerate_preceding_space_interval.py --length 52 --dimension 10 --multiplicity 0 --first 0 --count 100 --native-directory build/spaces --work-directory work/c52_m10_n0_000000000
```

The profile is computed before enumeration and gives the exact number
of raw lifts. The default one-million-candidate budget rejects larger
intervals before allocation; subdivide them or explicitly increase
`--maximum-candidates` when appropriate. For full reproduction, partition
each source domain into disjoint exhaustive intervals, including every
multiplicity in the census. The output is a sorted signature ledger of
32 signature bytes followed by a little-endian ambient bitmap. No matrix
JSON expansion or whole-parent collection in Python is required.

At m=8 reconstruct the marked inputs with:

```text
python -B code/prepare_preceding_marked_inputs.py --length 52 --output-directory work/c52_m08_inputs
```

Each listed task is accepted by `build/spaces/c52_m08_marked_orbits` with
`--input TASK --output RECORDS --summary SUMMARY`. The analogous executable
is provided for length 50. The input manifest identifies primary and
alternative families, which must both be accounted for. The primary
quotient uses the primary records; the separate explicit alternative
cover certifies the other families.

For the smaller m>=12 domains, `code/finite_space_census.py` provides a
complete dimension-configurable recurrence and exact quotient, with
independent replay of every positive affine map and canonical separation.
See `finite_contraction_censuses.md` for its interface. The empty m=14
length-50 sector closes m=15,16, and the empty m=16 length-52 sector
closes m=17. The independent dimension bound excludes all larger m.

Combining sorted interval ledgers requires a full merge by signature and
support, followed by an exact quotient of each whole signature bucket;
independent quotients of arbitrary source intervals are not sufficient.
`code/space_native/length54/merge_m09_signature_batches.cpp` and
`merge_m10_signature_batches.cpp` have the corresponding dimension-fixed
record layouts, independent of support length. The generic
`exact_quotient_c52_m10` executable accepts the merged ledger and bucket
index. `code/finite_space_census.py` supplies the reference bucket,
positive-witness replay and catalogue canonical-matching operations.
