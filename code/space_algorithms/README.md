# Space Algorithms

These modules supply the dimension-independent mathematical operations.
They depend only on the Python standard library.

| Module | Operation |
|---|---|
| `quadratic_syndromes.py` | Degree-two evaluation columns, quotient syndromes and quadratic closure |
| `inverse_contraction_extensions.py` | Exact lift equations, affine-shear quotient, zero-syndrome fibre sets, reconstruction |
| `kernel_patterns.py` | Small self-orthogonal kernel-pattern codes with minimum weight four |
| `exact_affine_transporter.py` | Exact affine equivalence with explicit affine maps |
| `affine_tools.py` | Affine frames and an independent transporter implementation |
| `fast_affine.py`, `support_invariants.py` | Affine-invariant refinement for equivalence testing |
| `decomposition_audit.py` | Direct-sum components, independently computed by circuits and the multiplicative stabiliser algebra |
| `triorthogonal_utils.py` | Binary rank, moment and polynomial utilities |

For a contraction core `C` in dimension `m-1`, construct
`CoreExtensionContext.build(m-1, C)`. The external labels are degree-two
evaluation columns reduced modulo the evaluation span of `C`. A fibre set
is compatible exactly when its labels sum to zero. The exact subset
iterators support fibre sizes zero through five, using sorted pair tables
for the larger sizes. The dynamic-programming and Walsh-character counters
provide independent checks on the number of compatible subsets.

For each compatible set `F`, compute its unreduced degree-two syndrome,
solve with `context.solver.solve`, and enumerate `context.lift_complement`
to quotient by affine shears. `context.reconstructed_support(F, f)` then
returns the lifted support. The new coordinate is the least significant
bit. Full affine rank must be checked on every retained lifted support;
rank-deficient cores themselves must not be discarded.

The separately written `../space_lifts.py` implements the same moment
equations in a small reference API, useful for differential testing. The
fixed-width C++ kernels in `../space_native` implement the larger traversals.
