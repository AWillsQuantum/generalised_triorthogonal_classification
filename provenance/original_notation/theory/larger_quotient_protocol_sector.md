# The larger-quotient protocol sector

For a pointed stabiliser support `X`, write `Q_d(X)=V_d(X)/L_X` for the
distance-filtered logical label quotient defined in
`support_correspondence.md`. The lower-dimensional source blocks divide
the five-logical-row problem into `5 <= dim Q_3 <= 14` and
`dim Q_3 >= 15`. This is a partition by linear algebra, not an assumption
about which outputs are likely to occur.

## Exact input cover

The latter domain has 47,848 pointed supports in 651 nonempty source
blocks, arising from 726 distinct affine source spaces. The complete
pointing coordinates, source indices and choices of origin are in
`data/protocol_sectors/large_quotient_domain.json`.

Every listed pointing is reconstructed from the compressed source
catalogue. Its quotient dimension is independently recomputed. The
selected cardinality in each block equals the corresponding tail of
the complete quotient histogram. Exact linear equivalence and coordinate
permutation reduce these pointings to 15,637 support classes. The classes
are not distinguished merely by a coarse invariant: their keys encode
exact canonical binary row spaces.

`verify_large_quotient_domain.py` checks these statements. With `--native`,
it also recomputes every canonical key and every support automorphism
order. The source-block histograms themselves can be reproduced by the
complete pointing and label-cache operations. Checking the selected
inputs does not constitute a fresh enumeration of all origins.

## Removing zero columns

At distance at least three, a zero stabiliser syndrome forces the whole
matrix column to be zero. Deleting it preserves the output, exact distance,
leading error coefficient and row count, and reduces protocol length.
Consequently zero-column matrices are Pareto-redundant. Deletion also
identifies their distance-three label quotients and common-isotropy
conditions with those of the shorter support.

Normalising the 15,637 larger-quotient classes gives 15,340 classes.
Combining these with the finite length-52, affine-dimension-seven input
family gives 24,701 inputs and 17,993 distinct nonzero supports. This
combined set is specified by
`data/protocol_sectors/primitive_support_union.json`:

| Length n | Stabiliser dimension h | Support classes |
|---:|---:|---:|
| 51 | 7 | 2,653 |
| 52 | 7 | 3,572 |
| 52 | 8 | 479 |
| 53 | 8 | 297 |
| 54 | 8 | 10,981 |
| 54 | 9 | 11 |

The input map is explicit. Support keys reconstruct the representative
rows and coordinates exactly. The independent verifier checks rank,
triorthogonality and the distance-three and distance-four label dimensions.
Its optional native check certifies all normalisation equivalences,
including the equivalent coordinates used for individual finite censuses.
It does not infer a distance-four exclusion just from distance-three
data; some of these supports have nonzero `Q_4`.

## The exceptional five-dimensional tensors

There are two primitive, nondegenerate tensor types in dimension five:
`ell wedge omega` and `ell wedge omega + ell^3`. Their finite target
signatures in the native monomial convention are `0x60000` and `0x60001`.
Their complete tensor classification and the admissible-parent proof are
in `logical_dimension_closure.md` and `canonical_augmentation.md`.

The target census uses a nondegenerate three-dimensional restriction as
its seed, allows the required degenerate intermediate restrictions, and
tests parent eligibility within exactly that reachable family. The
recorded result is zero for both target types on all 17,993 supports.
The data contain the per-support terminal counts, with exact input keys.

This result also excludes primitive seven-dimensional protocols on these
supports: either primitive seven-dimensional type has the corresponding
primitive five-dimensional restriction. Restriction of logical rows
cannot decrease distance. It does not, by itself, exclude nonprimitive
five-, six- or seven-dimensional tensors.

The normalised geometry checks and a configurable interval of primitive
censuses can be reproduced with:

```text
python -B code/verify_primitive_support_sector.py --output work/primitive_geometry.json
python -B code/verify_primitive_support_sector.py --native build/utsp-native --work-directory work/primitive --canonical-check --replay-start 0 --replay-count 2 --workers 4 --output work/primitive_sample.json
```

The second command checks only the stated primitive replay interval.
Taking `--replay-count 17993` recomputes the whole primitive census;
smaller disjoint intervals are equally valid when their union is checked.
The certificate distinguishes geometric verification from enumeration
replay and does not label a sample as a complete replay.

## A finite nonprimitive subsector

`data/protocol_sectors/nonprimitive_marked_subsector.json` gives 43
nonzero supports with 431 marked orbits of nonprimitive five-dimensional
logical subspaces. Every orbit has its exact canonical key, orbit size,
quotient representative and protocol matrix. The q3 chain cover generates
these orbits; the separate primitive target check supplies its two
exceptions. Together these are a complete five-dimensional census on the
43 supports.

All 431 matrices have exact distance three. Their complete output
profiles fail every necessary six-dimensional profile, and the primitive
seven-dimensional restriction is absent. The hyperplane and primitive
restriction theorems then exclude all logical dimensions at least six
on these 43 supports. Their matrices are covered by the released frontier
at the same output and exact distance, with no larger n or S.

```text
python -B code/verify_large_quotient_domain.py --native build/utsp-native --work-directory work/large_domain --workers 4 --output work/large_domain.json
python -B code/verify_nonprimitive_subsector.py --native build/utsp-native --work-directory work/nonprimitive --workers 4 --output work/nonprimitive.json
```

The latter replay compares all marked keys and orbit sizes, not just the
number of Pareto points.

## The chain-sector aggregate

The remaining finite domain is partitioned explicitly in
`data/protocol_sectors/nonprimitive_chain_sector.json`:

- 4,051 length-53 supports contain zero and are dominated at every logical
  dimension by the corresponding length-52 supports.
- The 43-support nonprimitive census is the finite subsector above.
- 11,543 supports form the remaining q3-chain domain. The aggregate census
  has 4,717 positive inputs and 6,826 zero inputs. Its 5,632 output
  candidates reduce to ten sector-frontier witnesses.

The successor tests retain 73 six-dimensional obligations. Their complete
aggregate has 177 marked orbits, representing 12,906 nondegenerate logical
subspaces, and one output key, `0x400a08a400`. Its sector frontier has one
point. All eleven q5 and q6 frontier witnesses are independently checked
and covered by the released Pareto frontier at their exact distance.
The q6 output profile excludes nonprimitive q7 extensions; the separate
primitive target census excludes primitive q7. The restriction theorem
then closes the higher logical dimensions on this finite domain.

A constant stabiliser row gives a useful additional exclusion. Its overlap
with one logical row forces every logical weight even, and its overlap
with two logical rows forces every logical pair overlap even. Thus the
logical tensor is alternating. The complete five-dimensional tensor census
has just one nondegenerate alternating orbit, the primitive type
`ell wedge omega`. On a support where that primitive type is absent, the
constant row therefore excludes all intrinsic q5 outputs without a
nonprimitive enumeration. The verifier records an explicit coefficient
vector for every such constant row in the chain domain.

These data are aggregate census results over an explicit finite input set,
not individual result files for every input. `verify_chain_sector.py`
checks their input partition, count consistency, successor-profile
conditions and all frontier witnesses. Its exclusion conclusions use the
stated completeness of those aggregate censuses. It does not present
consistency checks as a fresh enumeration of all logical subspaces.

The essential accelerated census can be reproduced on any input interval:

```text
python -B code/replay_chain_sector.py --q 5 --start 0 --count 1 --native build/utsp-native --work-directory work/q5_interval --workers 4 --output work/q5_interval.json
python -B code/replay_chain_sector.py --q 6 --start 0 --count 1 --native build/utsp-native --work-directory work/q6_interval --workers 4 --output work/q6_interval.json
```

The full index ranges are `[0,11543)` and `[0,73)`, respectively. Complete
reproduction must cover them without gaps and aggregate their exact
output, distance and footprint results. The primitive q5 target census
is a separate part of the complete q5 result. None of the input-cover,
primitive-only, finite-subsector or aggregate-consistency certificates
alone is a global classification certificate.
