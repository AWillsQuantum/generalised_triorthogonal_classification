# C++20 mathematical kernels

The source implements support pointings, label quotients, common isotropic
subspaces, exact coordinate-group actions, tensor equivalence, logical
extensions, distance and Pareto comparisons. It does not submit computations
to external services. The command-line program exposes the same mathematical
operations for finite input lists.

The dependency is Bliss 0.77, whose source archive and licence are included
in `vendor/bliss-0.77.zip`. It uses the system GNU Multiple Precision library
(GMP). A POSIX environment, Python 3.11+, a C++20 compiler, `ar`, and GMP
development headers and library are required.

```text
python build.py --test
build/utsp-native --help
```

`--build-dir` selects an alternative build directory. Nonstandard GMP
installations can be selected with `--gmp-include` and `--gmp-library`.
The latter accepts a library file path. The build verifies the vendored
Bliss archive's SHA-256 digest before extracting it.

## Mathematical interfaces

| Header | Operations |
|---|---|
| `native_core.hpp` | Binary row reduction, label quotients, mixed bilinear forms, distance |
| `origin_orbits.hpp` | Pointed supports modulo affine support automorphisms |
| `support_automorphisms.hpp` | Exact support automorphism actions |
| `marked_code_canonical.hpp` | Marked-code canonicalisation and complete admissible-parent extensions |
| `tensor_canonical.hpp` | Intrinsic output tensors and exact basis equivalence |
| `tensor_primitives.hpp` | Primitive exceptional tensor families |
| `tensor_extensions.hpp` | Hyperplane extensions of tensors |
| `tensor_unmarking.hpp` | Removing a distinguished tensor hyperplane |
| `tensor_decomposition.hpp` | Tensor direct sums |
| `protocol_catalogue.hpp` | Protocol witnesses and Pareto comparison |
| `sparse_completion.hpp` | Conditional unital completion |
| `puncture_embedding.hpp` | Parent-puncture embedding |

The support-first reduction, not existence of a projective completion,
governs the full matrix scope. The sparse-completion and puncture interfaces
are optional constructions and must not be imposed as completeness filters.

`tests/primitive_parent_regression.cpp` includes a positive 48-column matrix
and compares target-subspace enumeration with independent constructions.
The `--independent` option additionally checks the zero-hyperplane and
unseeded constructions. The other tests check stabiliser-sensitive caching,
nonfaithful quotient actions, canonical incidence selection and explicit
small-group actions against graph canonicalisation.

Passing these finite checks validates useful contracts, not the full
classification. Input-domain coverage and the corresponding exhaustive
outputs are separate evidence requirements.
