# Lossless compact supports

Let a unital triorthogonal space have length `c` and dimension `m+1`.
Choose a generator matrix whose first row is all ones and write its columns
as `(1,x)`, where `x` belongs to a set `X` in `F_2^m`. Pairwise distinct
columns mean that `X` is a set rather than a multiset. Full row rank means
that `X` affinely spans `F_2^m`.

Recording `m`, `c` and `X` therefore records a complete representative of
the space. Its generator matrix, indicator function, polynomial and overlap
checks are recoverable; they need not be stored alongside it. A record's
position in a shard provides its index. The ambient dimension `m` remains
explicit, including when some variables do not occur in the polynomial.

## Coordinate convention

A point is an integer `x` in `[0,2^m)`. In its zero-padded binary expansion,
the leftmost bit is `x_1`, and the rightmost bit is `x_m`. Columns are ordered
by increasing point integer. The first generator row is all ones and the
remaining rows evaluate `x_1,...,x_m` in that order.

## Fixed-weight enumerative encoding

For `0 <= a_1 < ... < a_c < 2^m`, define

```text
R(X) = sum_{i=1}^c binomial(a_i,i).
```

This is the colexicographic rank, an integer from zero through
`binomial(2^m,c)-1`. It is a bijection: the largest point is the largest
`a_c` with `binomial(a_c,c) <= R`; subtract that term and repeat. Pascal's
identity gives the disjoint rank intervals for successive largest points,
and induction proves that the inverse is unique.

A fixed-width rank needs

```text
ceil(bit_length(binomial(2^m,c)-1)/8)
```

bytes, using unsigned little-endian integers. For the complete length-54
dimension counts, this is 9,081,437,867 bytes before compression or headers.
This is a bound for this encoding, not a claim of optimal compression.

## Bitmap encoding

Alternatively store the integer `sum(2^x for x in X)` in `ceil(2^m/8)`
little-endian bytes. Its fixed-weight redundancy can be recovered by a
general-purpose compressor, and structural similarity between nearby
records can make this representation smaller than compressed ranks.
Both encodings describe exactly the same support; neither changes its basis.

## Shards

An uncompressed shard consists of the following 32-byte header, followed by
fixed-size records:

| Field | Type | Meaning |
|---|---|---|
| magic | 8 bytes | ASCII `UTSPACE1` |
| codec | uint8 | 1 = colex rank; 2 = bitmap |
| m | uint8 | ambient affine dimension |
| c | uint16 | support cardinality |
| first_index | uint64 | zero-based index of the first record |
| count | uint64 | number of records |
| record_bytes | uint32 | byte width derived from the codec, m and c |

All integer fields are little-endian. No padding is inserted. Empty
catalogue sectors have no payload records. A reader rejects inconsistent
widths, invalid ranks, incorrect weights, truncation and trailing bytes.
Compression is outside this format and is declared by the catalogue index.

## Recovering the indicator polynomial

Start with a truth table that is one exactly on `X`. Apply the binary subset
Moebius transform: for each bit `b`, XOR the coefficient at `t xor b` into
the coefficient at `t` whenever `t` contains `b`. The result is the algebraic
normal form over `F_2`, with square-free monomials indexed by bit masks.
Repeating the transform recovers the truth table. The coordinate convention
above specifies which variable belongs to each monomial bit.

For a triorthogonal support of even weight, all moments of degree at most
three vanish, equivalently its indicator lies in `RM(m-4,m)`. Checking
these moments and affine rank proves validity, but does not prove either
inequivalence or completeness of a list; those are separate audit tasks.
