# Complete finite contraction censuses

`code/finite_space_census.py` implements the inverse-contraction recurrence
of `space_classification.md` without requiring a separate algorithm for
each length or affine dimension. Its configured native kernels use the
same affine invariants and exact basis-image search throughout.

For a target `(c,m)`, the domain consists of every singleton core of length
`c-2N`, where `0<=N<=floor(binomial(c,2)/(2^m-1))`. All intrinsic core
dimensions allowed by the global dimension bound are included. There is
one useful exact reduction: when `N<=3`, a full-rank target must have a
full-rank quotient core. For `N=0` this follows from injectivity of the
projection. Otherwise an affine function vanishing on the core gives an
even binary word on the N complete fibres. A nonzero such word would
have weight two, and the linear-coordinate moments would force its two
distinct fibre points to coincide. Thus the affine function is zero.

For larger N the code retains proper-span cores and checks the affine
rank of each reconstructed candidate. A fixed embedding of each intrinsic
core is sufficient because all of its embeddings into the quotient
ambient space lie in one affine orbit. The full complement of that
embedding is available to the complete fibres.

For each core, binary elimination gives the quadratic evaluation rank,
the syndrome of every external point, and a basis of lift functions
modulo affine shears. Two independent fixed-cardinality procedures count
and enumerate the zero-XOR fibre subsets. Each compatible subset is
combined with every lift-basis coefficient. The zero coefficient is
discarded only for N=0, where it is exactly the rank-deficient affine
graph. No marked-orbit reduction or minimum-direction filter is needed
for these finite domains.

## Exact quotient and independent checks

Candidate supports are sorted by affine-invariant signatures. Signature
collisions are resolved by the exact affine transporter, never by assuming
that equal signatures imply equivalence. Each candidate is assigned to a
class with an explicit affine map from its representative.

The independent array verifier checks every candidate's affine rank and
all moments of degree at most three, then evaluates every affine map on
every support point. A map onto a full-affine-rank target is necessarily
invertible. It also checks the class member counts, source intervals,
representative positions and the absence of trailing records.

Finally, the codeword-incidence canonicaliser independently separates
every pair of output classes. Its chosen codewords span the complete
unital row space, so this is an exact affine-equivalence test. Thus the
coverage inference combines an exhaustive finite generation domain,
explicit positive witnesses and independent separation of the classes.
Completeness of the predecessor catalogues remains an explicit inductive
premise.

For example, a configured kernel pair is built with:

```sh
python code/space_native/generic/build.py --length 44 --dimension 8 \
  --build-dir build/configured
python code/finite_space_census.py --length 44 --dimension 8 \
  --native-directory build/configured \
  --graph-native build/native/affine_codeword_graph \
  --work-directory audit/c44_m08 --output audit/c44_m08.json
```

The explicit candidate budget guards against accidentally opening a much
larger domain. For a large catalogue one can subdivide candidate generation
and perform an external signature merge before the same exact quotient.

## Distributed finite proofs

`data/space_contractions/finite_sectors/index.json` specifies 54 finite
sectors: all required dimensions above seven through length 46, and the
length-48 successors in dimensions nine through fifteen. They contain
26,799,696 candidates and 22,339 affine classes. Length 44 has 863 classes
in total, with dimension counts `35,346,346,111,23,2` for `m=7,...,12`.
Length 46 has 2,403, with counts `980,1078,291,46,7,1` for `m=8,...,13`.

Each sector includes its exact predecessor domains, quadratic profiles,
class representatives, catalogue-index permutation and explicit affine
frames. The complete positive-assignment streams and signature bucket
indices are losslessly compressed. Candidate bitmaps and sorted ledgers
are deterministically regenerable and bound by their digests; they need
not be distributed as additional copies. These data are independent of
any particular job or task partition.

The public collection verifier checks the entire source cover, matches
all representatives to the compact catalogue and reads back every
compressed witness stream. With native executables supplied, it also
regenerates all candidates, replays every positive affine map and
independently recomputes every class's codeword-incidence canonical form:

```sh
python code/verify_finite_space_censuses.py \
  --native-directory build/configured \
  --graph-native build/native/affine_codeword_graph \
  --work-directory audit/finite-spaces \
  --output certificates/finite_space_censuses.json
```

Build a configured signature/exact-quotient pair for each listed `(c,m)`
using `code/space_native/generic/build.py`. The verifier's positive-witness
replay needs only the signature executable; `finite_space_census.py`
also reproduces the exact quotient from scratch.

The induction starts from the complete dimension-at-most-seven orbit-mass
base. The length-48 successors additionally depend on the separately
certified dimension-eight marked-contraction sector described in
`length48_marked_contractions.md`. Its dimension-sixteen sector is empty
without further enumeration: the pair-average bound is zero and its
dimension-fifteen predecessor is empty.
