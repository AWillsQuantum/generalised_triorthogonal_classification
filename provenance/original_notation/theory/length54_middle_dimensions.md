# Length 54 in affine dimensions nine through twelve

## Complete contraction domains

For a 54-point support in affine dimension m, the 1,431 unordered pairs
are partitioned by their nonzero differences. Some direction therefore
has at most `floor(1431/(2^m-1))` complete fibres. This gives precisely
the following complete source cover.

| Target m | Complete fibres N | Core length | Core dimension | Core classes |
|---:|---:|---:|---:|---:|
|9|0|54|8|17,777,766|
|9|1|52|8|1,321,156|
|9|2|50|8|95,657|
|10|0|54|9|259,202,914|
|10|1|52|9|5,601,443|
|11|0|54|10|15,617,101|
|12|0|54|11|552,230|

All these cores have full affine rank. For N=0 this follows from
injectivity of projection on the support. More generally, suppose an
affine function vanishes on the singleton fibres. On the full fibres,
its values form an even word. A word of weight two would force the two
selected fibre points to be equal by the coordinate moment equations.
Thus any nonzero such word has weight at least four, which is impossible
for N<=2. No proper-span source family is omitted here.

These source classes are read from the final compact catalogues. Their
completeness is an inductive premise, separate from the lift census.

## Exact lift profiles

Let C span m-1 affine dimensions, and let E_2(C) evaluate square-free
monomials of degree at most two on its points. The dimension of the
homogeneous lift functions modulo affine shears is

```text
ell(C) = |C| - rank(E_2(C)) - m.
```

For each exterior point y, reduce its quadratic evaluation vector modulo
the span of the evaluation vectors of C. A full-fibre set is compatible
exactly when its reduced labels sum to zero. If e_N(C) is the number of
compatible sets of size N, the number of full-rank marked lifts is

```text
N=0: 2^ell(C)-1,
N>0: e_N(C) * 2^ell(C).
```

The subtracted zero coset is exactly the affine graph of insufficient
rank. Positive complete fibres supply the new affine direction. For N=1,
e_1 counts zero labels; for N=2 it counts pairs of equal labels. These
formulas compute the complete lift count before constructing any child.

In dimensions eleven and twelve all directions chosen by the bound have
N=0. Fresh complete source-profile checks are supplied, with each parent
read from the compressed catalogue and validated by the native kernel.
The independent Python reference checks quadratic ranks on a stratified
sample. The input preparation streams records rather than retaining the
whole parent catalogue in memory.

## Minimum-direction reduction at m=9

A marked occurrence may be discarded when another direction has fewer
complete fibres. Every affine class has a globally minimum direction,
and the preceding bound puts that minimum in {0,1,2}. Retaining all ties
at the minimum therefore preserves at least one occurrence of every
target class. It does not yet quotient the full affine group.

The complete finite m=9 census has 61,732,322,090 full-rank lifts before
this filter and 5,322,436,668 retained marked occurrences. All 95,657
two-fibre sources are included; their 36,589,823,888 marked lifts leave
no occurrence whose marked direction is globally minimum. This last
statement is a finite census result, not an assumed theorem or a reason
to skip those sources in a reproduction.

The retained candidates are first quotiented in two subsets and then
globally. A local representative carries its number of original candidate
members. This weight must be preserved in the merge and global quotient.
There are 516,920,989 local representatives and 259,202,914 final affine
classes. Different invariant buckets are inequivalent; within a bucket
the exact affine transporter remains mandatory. Both positive affine
witnesses and negative exact comparisons are independently checked.

## Remaining finite quotients

The m=10 census uses the complete N=0,1 lift domain without a
minimum-direction filter. Its 3,138,829,689 marked lifts yield 15,617,101
affine classes. There are 531 contiguous parent-profile intervals. The
N=1 rank histogram alone does not determine the lift total: the compatible
fibre count and lift dimension must be evaluated together for each source.

| Target m | Contributing cores | Graph lifts | Exact affine classes |
|---:|---:|---:|---:|
|11|9,952,735|39,726,605|552,230|
|12|394,799|1,275,117|21,585|

The zero-contribution cores are included in the source-profile domain,
even though they produce no candidate. Complete graph enumeration,
coarse invariant sorting and exact affine quotient then give the table.
The member counts of the final classes partition the graph lifts.

## Evidence and reproduction

`data/space_contractions/middle_sectors/` contains the source domains,
finite generation counts, exact quotient accounting, and bindings to
the final catalogues and accelerated kernels. These are aggregate finite
censuses, not the individual affine assignment streams. Checking their
arithmetic does not amount to freshly repeating all affine searches.
`code/verify_middle_space_censuses.py` makes that distinction explicit and
checks source coverage, rank-nullity totals, weighted-member preservation,
quotient intervals and output counts.

The mathematical input interfaces are:

```sh
python code/prepare_contraction_interval.py --dimension 9 --multiplicity 2 \
  --first 17 --count 256 --native-directory build/space \
  --work-directory audit/m09_n2 --output audit/m09_n2.json
python code/prepare_contraction_interval.py --dimension 10 --multiplicity 1 \
  --first 17 --count 256 --native-directory build/space \
  --work-directory audit/m10_n1 --output audit/m10_n1.json
python code/prepare_zero_fibre_input.py --dimension 12 --count 552230 \
  --native-directory build/space --work-directory audit/m12 \
  --output audit/m12_sources.json
python code/verify_middle_space_censuses.py --output audit/middle_censuses.json
```

`prepare_contraction_interval.py` accepts every admissible multiplicity
at m=9 and m=10 and any contiguous parent interval. The m=9 output is an
ambient bitmap followed by the `<4B4x3Q` source suffix; it is accepted
directly by `m09_minimum_filter_kernel` (and its zero-fibre specialisation).
The interval indices refer to the compact parent catalogue, not to an
ordering by profile difficulty. The m=10 outputs are a neutral weighted
parent wrapper and a four-byte-per-source profile file, consumed directly
by `m10_direct_fused_signature`. The neutral wrapper has zero signature
words and member count one: it is not an affine-invariant or mass ledger.

For m>=11, `prepare_zero_fibre_input.py` produces the source records for
the corresponding fused graph-lift kernels. Their candidate masks are
sorted by affine invariants and passed to exact quotient and independent
witness-replay kernels. All these fixed-width mathematical operations are
in `code/space_native/length54/`; `m09_profile_interval.cpp` reuses the
same validated rank-and-fibre calculation with a streaming input interface.
An independent reproduction may choose any subdivision of the finite
domain, provided it accounts for every source and does not split an
invariant bucket during the final exact quotient.

The complete m=11 and m=12 source profiles have also been recomputed
from the distributed compact parents. The m=12 candidate ledger can be
regenerated in full with:

```sh
python code/verify_zero_fibre_ledger.py --dimension 12 \
  --sources audit/m12/sources.bin --native-directory build/space \
  --work-directory audit/m12_ledger --output audit/m12_ledger.json
```

The complete regenerated sorted ledger and an independently reconstructed
signature-bucket index match their retained binary digests exactly. This
checks all 1,275,117 candidates, not a sample; the exact affine quotient is
a separate stage. The verifier has an explicit working-memory guard.
