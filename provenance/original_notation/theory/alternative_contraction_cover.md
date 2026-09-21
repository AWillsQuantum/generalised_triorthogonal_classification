# Alternative-contraction covers

This reduction avoids a second affine quotient of children which are
already covered by another contraction direction. It applies to the
length-50, length-52 and length-54 sectors in affine dimension eight.

## Covering lemma

Let a complete finite family of contraction cores be partitioned into
P and E. Suppose the inverse-contraction construction is exhaustive for
every core in P. For every marked child of a core in E, suppose there
exists a nonzero direction whose odd-fibre contraction is affinely
equivalent to a core in P and satisfies that family's admissibility
conditions. Then the children of P alone cover every affine class in the
union of the two families.

Indeed, an invertible affine change of coordinates sends the certifying
direction to the marked last coordinate, and identifies its contracted
core with the representative in P. The resulting support is one of the
admissible lifts of that representative. Exhaustive marked-orbit
enumeration therefore includes an equivalent support. The argument does
not require a lookup in the final child catalogue.

The certificate must cover every marked child of E, not merely every
core in E. A core can have children with different alternative directions.
Checking one child or sampling directions cannot establish this lemma's
finite premise.

## The finite core partition

The complete seven-dimensional catalogues contain 35, 98 and 188 classes
at lengths 44, 48 and 52. Use zero-based indices throughout. Set

```text
E = {(48,7,96), (48,7,97), (52,7,186), (52,7,187)}.
```

All other seven-dimensional classes of these lengths belong to P.
Only cores whose length does not exceed the child length are admissible.
The complete inputs are reconstructed from the compact catalogue and
the complete affine stabilisers in `certificates/low_dimensional_spaces.json`.
At length 54 this partitions all 321 cores into 317 primary and four
alternative families.

For a contracted core C, let

```text
D_C(a) = |C intersect (C+a)|,   a != 0.
```

Both the multiset of D_C(a), and the multiset over x in C of the sorted
lists `(D_C(x+y): y in C, y != x)`, are affine invariants. A full-rank
triorthogonal core of an allowed length belongs to P whenever these
invariants separate it from every member of E of the same length.
This implication uses completeness of the seven-dimensional catalogue.
Equality of the invariants is never treated as proof of equivalence:
it leaves the candidate unresolved. A weight-44 core is immediately
in P, because E contains no weight-44 class.

## Exact finite checks

`code/verify_alternative_contractions.py` constructs all ten applicable
marked families, using the same affine-shear quotient and complete core
stabilisers as the primary enumeration. A separate native checker tests
every child's length, affine rank and degree-three moments, checks that
its marked contraction is the declared source, and records one
certifying direction. A second pass replays every recorded direction.
Python reference checks independently sample each stream.

For zero complete fibres, the affine graph is the unique rank-deficient
coset; it is removed. All other cases here have positive complete fibres
and full-rank cores, so every lift has full affine dimension eight.

The data in `data/space_contractions/alternative_cover/` contain native
inputs, compressed candidate streams and compressed direction streams.
Each candidate is four little-endian unsigned 64-bit support-mask words
followed by its unsigned 64-bit marked orbit size. The following direction
stream has exactly one byte per candidate. Direction zero would mean an
unresolved case and is forbidden in a complete cover certificate.

The row counts and orbit-size sums are independently compared with the
complete generation summaries and the exact pre-enumeration fibre counts.
These are marked-orbit masses, not affine orbit masses under AGL(8,2).

| Child length | Core families | Marked representatives | Raw marked pairs |
|---:|---:|---:|---:|
|50|2|9,360|327,680|
|52|4|334,840|13,008,894|
|54|4|8,471,112|341,508,096|

All ten families are fully covered: there are no unresolved children.
`certificates/alternative_contractions.json` records full regeneration
and replay. The compressed proof data occupy approximately 40 MB.

```sh
python code/space_native/build.py --build-dir build/space \
  --select c50_m08_marked_orbits c52_m08_marked_orbits \
  m08_marked_orbit_kernel alternative_contraction_cover
python code/verify_alternative_contractions.py \
  --native-directory build/space --work-directory audit/alternative_cover \
  --output audit/alternative_cover.json
```

The default mode replays the distributed proof streams. The `--create`
mode regenerates complete marked families and constructs all witnesses;
it requires an output root in which these proof data do not yet exist.
All code uses only public mathematical inputs. The primary child census
and its exact affine quotient remain separate proof obligations.
