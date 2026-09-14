# Generalised triorthogonal protocols through length 54

This mathematical repository contains the Pareto classification of binary
generalised triorthogonal protocols of length `n <= 54` and exact distance
`d_Z >= 3`, together with the complete unital-space catalogues used in the
classification. `AUDIT_STATUS.json` records distribution verification.

## Scope and results

Matrices have full row rank and pairwise distinct complete columns; one
zero column is allowed. Outputs are identified under **CNOT+S equivalence**,
not full Clifford equivalence. Logical dimension is unrestricted. Only
nonzero intrinsic magic outputs are keys. At fixed output and **exact**
distance, the frontier minimises `(n,S)`, where `S` is the total number of
matrix rows, with all logical rows counted throughout. Compressed footprint
is not an objective. The reductions cover both completable and
non-completable matrices in this scope.

The frontier has **74 points for 62 outputs**: 67 at distance three, five
at distance four and two at distance five. The largest undominated protocol
has length 53. The counts are of Pareto points, not of all protocol matrices.

| Maximum protocol length | Frontier points |
|---:|---:|
|48|52|
|50|53|
|52|73|
|54|74|

`data/protocols/pareto_frontier.json` contains every frontier matrix, exact
distance, leading error coefficient and representative gate, with explicit
output-basis certificates. Gate representatives are heuristically simplified;
minimum gate-factor counts are not claimed.

The affine-equivalence space catalogue contains **301,029,259** representatives
through length 54, including decomposable spaces. Of these, 293,172,583 have
length 54. All compact space data together occupy 4.84 GB; length 54 occupies
4.73 GB. Matrices and indicator polynomials are reconstructed from the support
encoding, retaining the ambient dimension even when variables do not occur.

## Start here

Use Python 3.11 or later and the packages in `requirements.txt`. The C++20
build additionally requires a POSIX environment, `g++`, OpenMP, `ar` and GMP
development headers and library. Bliss source and notices are included;
see `THIRD_PARTY_NOTICES.md`.

From this directory:

```sh
python -B code/verify_release_integrity.py
python -B -m unittest discover -s tests -v
python -B code/verify_protocols.py
python -B code/verify_classification.py --output ../classification_check.json
```

The last command joins the finite space induction, all protocol source
partitions, logical successor closure and independently verified frontier.
It checks the retained census bindings; it does **not** freshly repeat every
large enumeration. `AUDIT_GUIDE.md` distinguishes these claims and gives a
map to the complete regeneration interfaces. A checksum or valid witness
alone is not a classification proof.

For a readable space example:

```sh
python -B code/read_spaces.py data/spaces/c16/m04/spaces_000000000.utspace.zst --limit 1 --matrix --polynomial --validate
```

To check every compact record, its rank, all unital overlap conditions and
every index interval and digest:

```sh
python -B code/verify_spaces.py --validity --require-complete
```

The included full read-back certificate covers all 301,029,259 records.
The batched checker uses bounded memory. `--reference` selects an independent
scalar implementation. Neither mode substitutes for the affine quotient
and completeness arguments.

## Organisation

- `theory/`: definitions, proofs, finite domains and reproduction specifications.
- `code/`: mathematical algorithms, independent readers and verifiers.
- `code/native/`: accelerated protocol enumeration and exact canonicalisation.
- `code/space_algorithms/`, `code/space_native/`: space recurrences and exact affine quotients.
- `data/`: compact representatives and finite mathematical censuses.
- `certificates/`: validity, equivalence, finite-cover and composition evidence.
- `tests/`: mathematical, format and corruption-rejection checks.
- `resource_estimates/`: the conditional length-56 CPU forecast and its inputs.

The theory starts with `theory/support_correspondence.md` and
`theory/space_classification.md`. Read `theory/canonical_augmentation.md`
and `theory/logical_dimension_closure.md` before changing logical-subspace
enumeration: nondegenerate outputs need not have nondegenerate hyperplanes.
The smallest complete native example is in
`code/native/PROTOCOL_ENUMERATION.md`.

The [length-56 resource estimate](resource_estimates/length56_core_hours.md)
gives a single planning estimate of **10^8 core-hours**, with reproducible
arithmetic and explicit assumptions for unmeasured work. It is not a
complexity bound or a further classification result.

The distribution is below the 25 GB target and 50 GB hard limit. Full
regeneration can require substantially more temporary storage than the
compressed distribution. Keep generated ledgers outside the release, or
use `verify_release_integrity.py --allow-extra` to permit local additions.
