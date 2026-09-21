# Necessary Output Profiles

Fix a pointed stabiliser support and a lower distance bound. Restricting
the logical row space preserves all mixed overlap equations and cannot
lower distance. Consequently every nondegenerate restriction of a realised
logical tensor must occur among the outputs on that same support and at
that same distance floor. This concerns the complete output profile, not
only the output of a Pareto-optimal intermediate protocol.

For a tensor `T` of dimension `q`, write `P(T)` for the set of intrinsic
output keys of its nondegenerate hyperplane restrictions. If `A` is the
complete output profile in dimension `q-1`, a necessary condition for
realising `T` is `P(T) subset A`.

## Six-Dimensional Necessary Profiles

A symmetric trilinear tensor satisfying `T(x,x,y)=T(x,y,y)` is determined
by a symmetric bilinear repeated-index form and an alternating part.
In dimension six the symmetric bilinear form has ten congruence types:
the zero form; one nonalternating type for each rank one through six;
and one alternating type for each positive even rank. Each fixed normal
form has `2^binomial(6,3)=2^20` alternating parts.

Thus ten sets of `2^20` tensors contain representatives of every tensor
orbit. They need not be disjoint as tensor orbits, which is immaterial for
a necessary-profile test. For every nondegenerate tensor compute all
63 hyperplane restrictions, discard degenerate restrictions, and identify
the remaining restrictions with the exact five-dimensional tensor census.
Keep the inclusion-minimal profiles. Every realisable nondegenerate
six-dimensional tensor must contain at least one of these profiles in
the support's available five-dimensional outputs.

The kernel is `q6_hyperplane_output_profile_census.cpp`. The computation
checks the entire ten-slice domain, not a sample. The absence of primitive
six-dimensional tensors, proved in `logical_dimension_closure.md`, ensures
the necessary profile is nonempty.

## Seven-Dimensional Anchored Exclusions

Suppose a nondegenerate seven-dimensional tensor has a nondegenerate
hyperplane whose output key is `K`. Fix that hyperplane to the canonical
six-dimensional tensor representing `K`. The tensor coordinate count is

```text
D(q) = q + binomial(q,2) + binomial(q,3),
D(7)-D(6) = 22.
```

All possibilities therefore occur among `2^22` anchored extensions.
For each extension, inspect its nondegenerate six-dimensional hyperplanes.
Every resulting output key must lie in the support's complete
six-dimensional output profile. Reject the extension otherwise. Finally
reject a degenerate seven-dimensional tensor itself.

The kernel `q7_q6_output_profile_gate.cpp` implements this test. It uses a
faster recogniser for certain six-dimensional outputs based on their
five-dimensional restriction profiles. Before that recogniser is used,
the ten-slice census checks every matching normal form with the direct
tensor canonicaliser and verifies that the profile determines the claimed
output key uniquely. `verify_output_profile_gates.py` performs that check
before enabling the optimisation.

An anchored exclusion does not exclude primitive seven-dimensional
tensors, because these have no nondegenerate hyperplane. They are checked
separately using the primitive tensor theorem and their lower-dimensional
restriction profiles. Nor does the finite tensor test alone certify that
a protocol enumeration has supplied every available output profile.
Both obligations must be bound in the global coverage certificate.

## Typed Keys and Exact Distance

An output key is a tuple of unsigned 64-bit words. Hexadecimal text is a
serialization, not a different mathematical type. Both sides of every
profile comparison are normalised to integer-word tuples before subset
tests. Empty keys and words outside the unsigned range are rejected.

The profiles use a distance floor. The final Pareto frontier instead fixes
the exact distance. When a branch is removed by Pareto dominance, the
incumbent must have the same intrinsic output and the same exact distance,
with no greater length or row footprint. Any possible larger exact
distance remains a separate obligation; it cannot be removed merely by
the existence of a shorter distance-three incumbent.

To recompute the finite gates:

```text
python code/verify_output_profile_gates.py --binary-directory build/native --output-directory build/output_profiles
```

The 175-input example in
`data/protocol_sectors/length54_q6_small_quotient.json` can be re-enumerated by
`verify_protocol_case_census.py`. Its certificate checks the complete
output profile of each listed support and the total isotropic and
nondegenerate counts. It is an explicitly bounded sector, not a certificate
for all pointed supports of length 54.

## Finite Support Refinement

The necessary-profile check can first be applied to the union of the
profiles in an input block. If the union contains no necessary profile,
every support in the block is excluded. If it does contain one, refine
to the individual support profiles: a match in a union need not occur
on any one support.

`verify_profile_refinement.py` checks this refinement for 21,077 blocks:
20,915 are excluded at the union level; the remaining 162 contain 33,758 positive
five-dimensional support profiles, of which 197 pass the six-dimensional
necessary test. Exact pointed-support canonicalisation identifies these
with the 175 cases in the finite q6 census. Their pointing origins and every
support point are independently recovered from the compact space catalogue
by `verify_pointed_sources.py`.

`verify_primitive_restrictions.py` exhibits primitive five-dimensional
restrictions of both primitive seven-dimensional tensors. Absence of both
required keys from a complete q5 profile excludes these exceptional q7
branches. `verify_finite_sector_closure.py` combines this with the anchored
q7 gate, complete q6 censuses, the distance-four label equations and explicit
catalogue dominators. It checks 178 finite supports, including three
larger-quotient cases at shorter lengths. Every dominance comparison fixes
exact distance three; the zero distance-four quotient excludes a missed
higher-distance metric on those supports.

These are coverage statements for the stated finite domains. The global
certificate must additionally show that its source partition covers all
pointed supports, including the blocks handled by other enumerations and
the theoretically excluded branches.

## Small Outputs at Exact Distance Three

There are 22 intrinsic nonzero output classes on at most four qubits:
one, two, six and thirteen in dimensions one through four, respectively.
They can be recovered without a separate list by taking the 23 degenerate
classes of the complete five-dimensional tensor census, discarding the zero
tensor, and restricting each tensor to a complement of its radical.
Adding zero directions preserves and reflects tensor equivalence, so this
procedure is complete and introduces no additional identifications.

Every projective protocol of distance at least three has distinct
stabiliser syndromes, hence `n<=2^h`. At lengths 49 through 54 this gives
`h>=6` and therefore `S=q+h>=q+6`. For each of the 22 outputs the distributed
frontier contains a distance-three protocol with `n<49` and `S<=q+6`.
These witnesses strictly dominate every exact-distance-three candidate in
that length interval. This does not exclude candidates of larger exact
distance, which require the corresponding filtered enumeration.

`verify_low_q_pruning.py` derives the radical complements, records explicit
equivalences with the frontier outputs, and checks all 22 metric
inequalities. The frontier matrices are independently verified by
`verify_protocols.py`.
