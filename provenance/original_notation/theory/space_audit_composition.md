# Composition of the parent-space classification

`code/verify_space_classification_cover.py` joins the finite space proofs
in increasing length and then increasing affine dimension. For every
permitted (c,m), it checks an explicit census or a proved empty successor;
it does not infer completeness from the number of files present.

The complete RM(3,7) orbit-mass base supplies all m<=7 sectors through
c=54. For m>=8 the pair-averaging bound excludes the empty contraction
core. Every possible positive core is a smaller-dimensional support of
length at most c, so all predecessors appear earlier in the induction.
The dimension bound

```text
m <= floor((c+floor(c/16))/3)-1
```

is independent of the finite censuses. Even retaining its conservative
maximum m=18 at c=54 causes no gap: the empty m=17 predecessor has no
nonaffine graph lifts. Odd lengths and positive lengths below 16 are
excluded by the moment and Reed-Muller arguments.

The proof components are:

1. The complete low-dimensional affine base and its full stabilisers.
2. The independently regenerated finite contraction domains through
   length 46, the complete marked c=48,m=8 domain, and all c=48,m>=9
   finite domains, with positive affine witnesses and independent
   separation of their class representatives.
3. The finite c=50,52 censuses, their complete source partitions and
   alternative-direction covers, described in `length50_length52_censuses.md`.
4. The finite c=54,m=8 and m=9,...,12 censuses, described in
   `length54_m08_census.md` and `length54_middle_dimensions.md`.
5. Every graph-lift coset and explicit affine map at c=54,m=13,...,17,
   and the empty-successor arguments at the remaining upper dimensions.

The resulting catalogue contains 301,029,259 affine classes. All distributed
supports have separately passed rank, degree-at-most-three moment, index
interval and lossless-decoding checks.

## What the composition certificate establishes

`certificates/space_classification_cover.json` checks that the declared
finite results close every mathematical sector and bind to the complete
distributed catalogue. No predecessor-completeness assumption remains
external to that induction. The correctness and exhaustiveness of the
listed large finite censuses remain explicit computational premises.

This distinction matters. The release contains the mathematical input
domains, accelerated enumeration and exact-equivalence algorithms, finite
census results and verification evidence. Some small and medium domains
also contain every compressed affine witness, allowing inexpensive full
replay. For the largest domains, an independent auditor must regenerate
the finite enumeration and exact quotient to check every computational
decision. Aggregate counts, source hashes and successful sample tests
are not a substitute for that regeneration, and are not presented as an
unconditional machine-checked proof.

The listed `finite_censuses_for_independent_regeneration` identify these
domains precisely. A fresh complete space-induction replay through length
48 is already included; a fresh rerun of every larger finite quotient is
not claimed. The protocol-completeness argument then consumes this
space classification as one of its mathematical inputs, separately from
the validity and Pareto properties of individual protocol witnesses.
