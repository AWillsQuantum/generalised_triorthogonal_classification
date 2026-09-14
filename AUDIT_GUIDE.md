# Audit guide

## Claims and computational evidence

There are three separate obligations: validity of each representative,
inequivalence of representatives, and exhaustion of the permitted domain.
The mathematical reductions and complete finite input domains are supplied
with code that can regenerate the finite censuses. Some large censuses are
stored as aggregate results rather than individual assignment streams.
Their checksums and consistent totals do not independently prove every
enumeration decision. Checking those decisions requires regeneration.

`certificates/classification.json` joins the space induction and protocol
cover, binds the finite inputs, verifies the exact frontier witnesses and
lists the retained-census premises explicitly. It is not an unconditional
formally machine-checked proof. Nor does it claim that every large finite
census was freshly repeated by the composition verifier.

The complete space induction covers every dimension allowed by the
independent Kneser bound and closes every contraction predecessor.
The protocol source intervals partition precisely those catalogues.
All logical successor obligations, including unrestricted logical dimension,
have either a complete finite census or a stated exclusion. The proofs also
justify removing zero complete columns and Clifford-trivial logical factors
at fixed **exact distance**, not merely at a distance lower bound.

## Short audit

Run the four verification commands in the README. They check all distributed
file hashes, unit tests, every frontier matrix and the global composition.
The complete space-data read-back is separate and longer. Its retained
certificate is `certificates/all_spaces.json`; reproducing it reads every
space record but does not repeat the full affine quotient.

The source and census certificates are bound to mathematical data files,
not to execution environments. Build checks and sampled differential tests
certify their stated tests only. A sample is never labelled a complete
finite-domain replay.

## Native build and a complete small reproduction

From the release root:

```sh
python -B code/native/build.py --build-dir ../audit-build/protocol --test
python -B code/space_native/build.py --build-dir ../audit-build/spaces
python -B code/space_native/generic/build.py --length 52 --dimension 10 --include-lifts --build-dir ../audit-build/spaces
```

The protocol build supports `--gmp-include` and `--gmp-library` for a
nonstandard installation. The configured space build supports multiple
lengths and dimensions, as described in its README and the sector documents.
Builds require no other copy of the classification repository.

`code/native/PROTOCOL_ENUMERATION.md` gives a complete reproduction starting
from the 16-point affine support: reconstruct the support, enumerate its
five pointing cases and enumerate logical rows to recover the 15-input T
protocol. The raw and marked-orbit enumeration routes provide independent
cross-checks. Do not discard intermediate logical spaces merely because
their protocols are Pareto-dominated: their descendants may have new outputs.

A second complete example regenerates the three length-32 spaces of affine
dimension eight, including all ten candidate lifts, every positive affine
map and independent separation of the three classes:

```sh
python -B code/space_native/generic/build.py --length 32 --dimension 8 --build-dir ../audit-build/spaces
python -B code/finite_space_census.py --length 32 --dimension 8 --native-directory ../audit-build/spaces --graph-native ../audit-build/protocol/affine_codeword_graph --work-directory ../audit-c32-m08 --output ../audit-c32-m08.json
```

Both examples, the complete protocol-domain composition, 171 Python tests,
four native protocol regressions, 82 standalone space units and seven
configured space variants are covered by
`certificates/standalone_reproduction.json`. That certificate binds the
tested sources and explicitly does not claim a fresh replay of every large
finite census.

## Space reproduction map

| Domain | Mathematical specification and complete interfaces |
|---|---|
| Affine dimensions at most seven | `theory/cubic_seed_completeness.md` |
| Finite recurrence through length 46; length 48 above dimension eight | `theory/finite_contraction_censuses.md` |
| Length 48, dimension eight | `theory/length48_marked_contractions.md` |
| Lengths 50 and 52 | `theory/length50_length52_censuses.md` |
| Alternative marked-core families | `theory/alternative_contraction_cover.md` |
| Length 54, dimension eight | `theory/length54_m08_census.md` |
| Length 54, dimensions nine to twelve | `theory/length54_middle_dimensions.md` |
| Length 54, dimensions thirteen and above | `theory/zero_fibre_space_certificates.md` |
| Full inductive domain and bounds | `theory/space_audit_composition.md` |

The fresh finite proofs through length 48 include every generated candidate,
positive affine map and target separation. At lengths 50, 52 and 54, the
global certificate lists 21 retained finite censuses for full independent
regeneration. Further complete replays and differential checks certify
specified portions of these larger domains, without replacing their other
obligations. Signatures are bucket keys only; each shared bucket requires
an exact affine quotient, and weighted representative merges must preserve
the original candidate multiplicities.

## Protocol reproduction map

| Domain or reduction | Mathematical specification and complete interfaces |
|---|---|
| Parent support and all five pointing cases | `theory/support_correspondence.md` |
| Through length 48 | `theory/protocols_through48.md` |
| Finite shorter-source blocks | `theory/protocol_cover_composition.md` |
| Lengths 49 to 52 | `theory/protocols_lengths49_to52.md` |
| Length 54, low-dimensional source blocks | `theory/low_dimensional_protocol_blocks.md` |
| Length 54, high-dimensional source blocks | `theory/high_dimensional_protocol_sector.md` |
| Larger quotients, primitive targets and chains | `theory/larger_quotient_protocol_sector.md` |
| Necessary output profiles | `theory/output_profile_pruning.md` |
| Exceptional primitive outputs and arbitrary dimension | `theory/logical_dimension_closure.md` |
| Exact frontier composition | `theory/protocol_cover_composition.md` |

For full regeneration, cover every specified input interval without gaps
or overlaps. Merge exact output/distance/length/footprint results across
all intervals, together with the preceding frontier. Necessary-profile
tests must use complete profiles and include the primitive exceptions.
In particular, the 11,543 q5 chain inputs, their 73 q6 successors and the
17,993 primitive target inputs are different obligations. None alone is
the complete larger-quotient classification.

The complete enumeration may be subdivided arbitrarily, but subdivision
does not alter mathematical scope. The supplied interfaces express finite
mathematical tasks; no execution-service setup is needed to understand or
reproduce them.
