# Length 54 in affine dimension eight

The inverse-contraction bound gives a direction with at most five full
fibres. Seven-dimensional cubic Reed-Muller weights force the admissible
multiplicities to be 1, 3 or 5, and therefore the core lengths to be 52,
48 or 44. The rank argument in `space_classification.md` excludes
proper-span cores in these cases.

## Finite domain and quotient

`data/space_contractions/c54_m08/domain.json` gives the complete 321-core
domain. The primary family P comprises all 35 length-44 cores, the first
96 length-48 cores and the first 186 length-52 cores. The other four
families E have the complete alternative-contraction proof described in
`alternative_contraction_cover.md`. This is a cover, not a disjoint
partition of the resulting affine classes: E contributes no new class.

| Quantity | Value |
|:---|---:|
|Primary cores|317|
|Alternative cores|4|
|Primary raw marked pairs|154,509,180,928|
|Alternative raw marked pairs|341,508,096|
|All raw marked pairs|154,850,689,024|
|Primary marked-orbit representatives|2,630,277,993|
|Primary coarse invariant buckets|17,777,737|
|Exact affine classes|17,777,766|

One primary core has no compatible fibre set and hence contributes no
child. It is explicitly included in the input domain and the zero census.

The coarse invariant has 17,777,708 buckets containing one affine class
and 29 containing two. Exact affine comparison, rather than invariant
equality, is essential. The finite quotient records 2,612,500,227 positive
and 1,175 negative transporter comparisons; every candidate has a checked
map to its assigned representative, including identity maps for the
representatives themselves. The count of marked-orbit representatives is
not the order of an AGL(8,2) orbit, nor its orbit mass.

## Evidence and its scope

`data/space_contractions/c54_m08/primary_census.json` contains the 317
primary marked-orbit counts, all 5,262 contiguous quotient intervals,
their exact comparison counts, complete bucket-size distributions, and
the binding to the ordered compact representative catalogue. Primary
core coordinates and lift parameters are bound to the same complete
public affine census used to reconstruct the native inputs.

This file is aggregate finite census evidence. It does not contain the
individual primary affine assignment streams. Its arithmetic and source
bindings are checked by `code/verify_space_m08_cover.py`; that check does
not pretend to repeat billions of transporter comparisons. A fresh full
check of the primary computational premise requires re-enumeration with
the supplied mathematical kernels. In contrast, the alternative-family
candidate streams and explicit direction witnesses are distributed and
have complete fresh generation and replay certificates.

The compact encoding certificate proves equality with the ordered
representatives and independently checks all weights, ranks and moments.
It does not by itself prove inequivalence or completeness. The induction
step combines the complete core census, exhaustive marked generation,
exact affine quotient, and alternative-contraction lemma.

```sh
python code/verify_space_m08_cover.py \
  --work-directory audit/c54_m08_inputs \
  --output audit/c54_m08_cover.json
```

This reconstructs all 321 native inputs, verifies the primary and
alternative partition, checks all aggregate quotient intervals, and
checks the compressed catalogue bindings. It leaves the two primary
finite enumeration premises explicitly labelled in its result.

## Essential accelerated operations

The primary census can be reproduced with the standalone C++20 units in
`code/space_native/length54/`, in the following mathematical order:

1. `m08_marked_orbit_kernel.cpp` enumerates the compatible full-fibre sets,
   their core-stabiliser orbits, and all lift-label orbits modulo affine
   shears. Inputs are reproduced by `prepare_space_contractions.py`.
   Its output record is four 64-bit support words and the marked orbit
   size. `verify_marked_contractions.py` gives a complete single-core
   reproduction and independent moment, rank and contraction check.
2. `m08_binary_signature_ledger.cpp` computes affine invariants and sorts
   finite candidate intervals. `merge_m08_signature_shards.cpp` merges
   these sorted streams, detecting literal duplicates. A record is four
   64-bit invariant words followed by four support words. Invariant
   hashing is used only for a coarse partition, never to assert equality
   of affine classes.
3. `m08_exact_affine_quotient.cpp` performs complete affine transporter
   tests inside each entire invariant bucket and emits assignments with
   explicit affine maps. Bucket boundaries must not be split.
4. `m08_verify_and_extract_quotient.cpp` independently replays every map,
   extracts representatives and records the number of assigned candidates.
   Its representative record is four invariant words, four support words,
   and an unsigned 64-bit member count.
5. `verify_m08_representative_rows.cpp` checks every representative's
   weight, rank and degree-three moments. The support words can then be
   encoded with `space_codec.py`; affine class ordering is not a
   mathematical invariant, so a reproduction using different core-group
   generators can also compare the resulting exact affine classes.

The quotient search and witness replay are distinct source units. To
audit a negative equivalence decision, rerun the exact search or use a
complete independent canonicalisation; replaying positive maps alone
would not establish inequivalence.
