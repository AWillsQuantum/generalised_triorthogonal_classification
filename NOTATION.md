# Resource notation

Current mathematical notation, adopted 20 September 2026:

| Quantity | Current | Historical |
|---|---|---|
| Protocol length / matrix columns | `n` | `n` |
| Spatial footprint / all independent matrix rows | `N` | `S` |
| Output qubits / logical rows | `k` | `q` |
| Stabiliser rows | `r` | `h` |

Thus `N = k + r`. The footprint counts all logical rows throughout; it is
not a scheduled or compressed footprint. The parent-space length `c` and
affine dimension `m` are unchanged. The phase gate `S` is also unchanged.

## Preserved evidence and compatibility

Frozen certificates, original source snapshots, historical computation logs,
hash-bound catalogues and compressed archives retain their exact bytes. Read
their resource symbols using the table above. A letter used for another
quantity is not renamed: for example a finite-field order, a unital-space
dimension, a polynomial degree or an index need not be a logical-row count.

Existing output IDs (`Q5_...`, `D5_...`), campaign names (`q5`, etc.), paths,
command-line options and legacy API/schema identifiers remain stable. In
particular, legacy `space_footprint_s`, `S`, `q`, `h` machine fields denote
`N`, `N`, `k`, `r` respectively when used as protocol resources. Descriptive
fields such as `logical_qubits` retain their names. This is an explicitly
documented compatibility boundary, not a different mathematical convention.

The current public witness catalogue uses schema
`triorthogonal-protocol-witnesses-v2`, with fields `n`, `N`, `k`, `r` and
Pareto objectives `n,N`. Its reversible notation adapter allows the original
mathematical checks to run unchanged. The preserved original release can
still be reconstructed byte for byte.

In the paper, columns are written `(lambda^(j), x^(j))`; a subscript indexes
a component, not a column. Named gates are upright. The double-pair count is
written `nu` and the generic Reed-Muller order `rho`, to avoid collisions
with the new footprint and stabiliser-row symbols.
