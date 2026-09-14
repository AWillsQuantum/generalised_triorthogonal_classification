# Accelerated Space Kernels

The `length54` directory contains standalone C++20 mathematical kernels.
Each source builds as a separate executable with `build.py`. These kernels
use packed fixed-width data, exact XOR elimination, Gray-code lift updates,
minimum-contraction filters, affine-invariant refinement and exact affine
transporters. They contain no external execution-service dependencies.

The dimension in a filename is the target affine dimension, not the
generator-matrix row count. An affine dimension `m` gives `m+1` rows.
Constants and fixed-width binary records are declared near the top of each
source. A record layout for one dimension must not be used for another.

| Family | Mathematical operation |
|---|---|
| `m08_marked_orbit_kernel` | Compatible fibre sets and lift classes under the marked-core stabiliser |
| `m09_minimum_filter_kernel*` | Inverse contraction with the least-intersection filter, including zero-fibre lifts |
| `*_source_profile*` | Quadratic ranks, compatible fibre counts, lift quotient dimensions |
| `*_fused_signature*`, `*_signature_kernel` | Lift enumeration, admissible-source filtering and affine invariant computation |
| `runtime_affine_signature` | The same invariant for a runtime ambient dimension |
| `*_exact*quotient*` | Exact affine classes within invariant buckets |
| `*weighted_representatives*`, `*representative_ledgers*` | Weighted unions of exact class representatives |
| `merge_*`, `project_*` | Exact sorted unions and projections of the mathematical ledgers |
| `*_decomposition_audit` | Independent direct-sum decomposition tests |
| `verify_*` | Independent row, count, invariant or affine-map verification |

The affine signatures are bucket keys, not equivalence certificates.
Hash collisions are harmless only because the exact transporter handles
every shared bucket. Likewise, the minimum-parent filter must compare
only sources admitted by the chosen source cover. Numerical signatures
alone do not establish that a support is present or absent.

Build on a POSIX system with a C++20 compiler and OpenMP:

```text
python code/space_native/build.py --build-dir build/spaces
```

Each executable validates its own argument set and reports its required
options on invalid input. Binary input adapters must also check record
count, dimension, support length, index intervals and source bindings.
`build_verification.json` certifies compilation only; it is not a
classification completeness certificate.
