# Data formats and a worked example

## Protocol catalogue

`pareto_frontier.json` has schema `triorthogonal-protocol-witnesses-v1`.
Its `protocols` array contains the 74 frontier witnesses, with `output_id`,
`q`, `d_Z`, `n`, `S`, `generator_matrix_rows` and `error_coefficient`.
Each matrix row is a binary string of length `n`; the first `q` rows are
logical and the remaining rows are stabilisers. The `outputs` array resolves
each `output_id` to a representative gate and tensor/basis certificate.
The identifier alone is not the human-readable gate. Gate simplification
uses CNOT+S equivalence and does not certify a minimum number of factors.

The exact Z distance is the minimum Hamming weight of a vector `v` with
`G_0 v = 0` but `G_1 v != 0`, over the binary field. The error coefficient
counts such vectors at the minimum weight. `S` counts every matrix row.
Outputs and exact distances are separate Pareto comparisons in `(n,S)`.

## Space catalogue

Extracting `space_catalogues.zip` creates `data/spaces/index.json` and the
`.utspace.zst` shards. The index gives lengths `c`, affine dimensions `m`,
record intervals and checksums. The unital space dimension is `m+1`, not `m`.
An index is local to its `(c,m)` sector, not a global protocol identifier.

After streaming Zstandard decompression, each shard begins with the 32-byte
little-endian structure `<8sBBHQQI`: magic `UTSPACE1`, codec, ambient dimension
`m`, length `c`, first index, record count and record width. Records have that
fixed width and encode subsets of the `2^m` binary points, not matrix entries.
Codec 1 stores the colexicographic subset rank as a little-endian integer:
for increasing points `a_1 < ... < a_c`, the rank is
`sum(binomial(a_i,i), i=1,...,c)`. Codec 2 stores a support bitmap, with point
`a` represented by bit `a`. The independent implementation is GitHub's
`code/space_codec.py`.

An increasing support `a_1,...,a_c` expands to an all-one row followed by
`m` coordinate rows, most significant bit first. These are the columns
`(1,a_j)` of a generator matrix. The indicator polynomial equals one on
the support and zero elsewhere; its Boolean algebraic normal form is
recovered by a binary Mobius transform. **Keep `m` explicit**, even if a
variable does not occur in the polynomial.

For `c=16,m=4`, the support is all integers 0 through 15. Its generator is

```text
1111111111111111
0000000011111111
0000111100001111
0011001100110011
0101010101010101
```

The indicator polynomial is `1` in **four ambient variables**, not a
zero-variable polynomial for purposes of this catalogue. Using GitHub's
reader from an assembled release reproduces this example:

```sh
python -B code/read_spaces.py data/spaces/c16/m04/spaces_000000000.utspace.zst --limit 1 --matrix --polynomial --validate
```

## Evidence and preservation

`classification_evidence.zip` holds structured finite domains, contraction
assignments, retained censuses, mathematical certificates and profiles at
their original `data/` and `certificates/` paths. Inner `.bin.zst` streams
remain compressed. These are inputs and evidence for the mathematical
algorithms, not another list of frontier protocols.

`delivery_layout.json` maps **every** file in `ORIGINAL_MANIFEST.json` plus
the manifest itself to at least one delivery location. It explicitly lists
intentional copies, including the two frontier JSON files. The standard-library
`code/delivery.py` in GitHub checks all copies, all ZIP members and the exact
original byte hashes. It can reconstruct the original directory without
the original development repository. Current delivery manifests separately
cover the new documentation and packaging helpers.

The preserved certificates distinguish validity, inequivalence and exhaustion.
Their audit composition uses retained finite census premises; hashing or
reading all representatives does not independently repeat every classification
decision. The full theory and regeneration interfaces are in GitHub.
