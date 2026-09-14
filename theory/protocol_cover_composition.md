# Composition of the finite protocol cover

The cover concerns full-projective matrices, exact distance at least three,
arbitrary intrinsic logical dimension and CNOT+S output equivalence.
At each fixed output and exact distance, the frontier minimises `(n,S)`.
Distance is part of the key, not a quantity on which dominance is applied.

The terminal finite sectors give the following cumulative frontier sizes.

| Maximum protocol length | Frontier points |
|---:|---:|
|48|52|
|50|53|
|52|73|
|54|74|

The final frontier has 62 distinct outputs. Its largest protocol length
is 53: no length-54 candidate adds an undominated point.

## Source and successor partitions

Through length 48 the main source family, two separate cubic sources and
the finite shorter-source blocks are disjoint and exhaust the parent
catalogue. The last family contains 1,708 classes at lengths 44 and 46.
All five pointing cases give 2,192,216 supports. Its 323 complete raw
logical censuses enumerate 124,099,558 common-isotropic subspaces,
including 1,973,866 with nondegenerate output tensors. All 217 emitted
finite-frontier witnesses are strictly dominated by the distributed
frontier. This family therefore changes no cumulative Pareto count.

There is a direct independent check on its pointing counts. If `T` is
the translation-period subspace of a support `Y`, the numbers of distinct
point sets in the five cases are

```text
2^m/|T|, |Y|/|T|, (2^m-|Y|)/|T|, 1, 1.
```

Equality of two translated supports is exactly equality modulo `T`.
Removing or adjoining zero preserves this equivalence, and the two
hyperplane cases have a different ambient dimension. These counts apply
when the entire length interval `|Y|-1,...,|Y|+1` is in scope, as here.

The largest logical quotient dimension in these blocks is twelve. Every
dimension through the smaller of eight and that quotient dimension is
explicitly enumerated. Whenever the quotient dimension exceeds eight,
there is no common-isotropic eight-dimensional subspace, so containment
excludes every higher logical dimension, independently of tensor type.

At length 52 the parent domain has 7,519,688 classes. It is partitioned
into 7,519,511 small-quotient source classes, 175 separately resolved cubic
classes, and two further cubic classes. The 5,068 deferred larger-quotient
pointings are already in the first source family; they are not additional
parents or pointings. The complete pointing count is 4,163,778,892.

At length 54 all 293,172,583 parent classes appear exactly once in 21,124
blocks. Necessary tensor profiles close most small-quotient successors.
The remaining 197 selected pointings have 175 exact linear support classes.
Their complete six-dimensional censuses and seven-dimensional profile
exclusions close all successors. Each selected pointing is checked against
the precise input position in its source block, not merely against an
equivalent parent support.

The 47,848 larger-quotient pointings have 15,637 support classes, partitioned
into 4,051 zero-column dominated classes, 43 explicitly marked cases and
11,543 chain cases. The separate primitive-target union has 17,993 classes.
Exact canonicalisation binds these domains to their finite logical inputs.
The primitive restriction and necessary-profile lemmas close all higher
logical dimensions, including dimensions not explicitly enumerated.

Both composition verifiers collect the terminal finite metrics and recompute
the nondominated union with the preceding cutoff. Every finite source and
every logical successor must have an exclusion or a complete census.
The source, profile, canonicalisation and result bindings are recorded in
`certificates/through52_protocol_cover.json` and
`certificates/through54_protocol_cover.json`.

## Reproduction and certificate scope

After building the native protocol executable, the composition checks are:

```text
python -B code/verify_finite_protocol_sources.py --native build/utsp-native --work-directory work/finite-sources --workers 4 --independent-labels --output certificates/finite_shorter_protocol_sources.json
python -B code/verify_protocol_cover.py --native build/utsp-native --work-directory work/cover52 --workers 4 --output certificates/through52_protocol_cover.json
python -B code/verify_through54_cover.py --native build/utsp-native --work-directory work/cover54 --workers 4 --output certificates/through54_protocol_cover.json
```

Use initially empty work directories. The checks repeat geometry, source
positions, exact input equivalences, candidate validity, profile implications
and frontier composition. They do not repeat every logical census. Complete
re-enumeration interfaces for those finite domains are given in the sector
documents, including interval-based chain enumeration.

The conclusions are conditional on completeness of the stated finite
logical censuses and of the parent-space catalogues. In particular, checking
the aggregate chain results is not a fresh enumeration of all chain inputs.
The parent-space completeness proof is a separate part of the release.
This distinction prevents a valid frontier witness or a matching aggregate
count from being mistaken for a global completeness certificate.
