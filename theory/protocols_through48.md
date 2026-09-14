# Finite protocol classification through length 48

The parent-space catalogue through length 48 has 30,669 affine classes.
The finite protocol domain is the disjoint union of a 30,667-source sector
and two separately certified cubic sources, `(c,m,index)=(48,7,96)` and
`(48,7,97)`. Catalogue indices in this distribution are zero-based.
All five pointing constructions are those in `support_correspondence.md`.
The bound here is on protocol length, not just parent length.

## Source and coordinate data

`data/protocol_sectors/through48_source_domain.json` identifies every parent
by its compact-catalogue index. An optional coordinate frame is
`[t,Ae_0,...,Ae_(m-1)]`, defining the affine map `x -> t+A x` before the
pointing construction. Thus variable order is explicit and does not change
the implicit ambient dimension. All nonidentity frames in this sector
are variable-order reversals.

The main sector has exactly 4,194,716 pointed supports after quotienting
origins by the complete affine automorphism group of each parent. The
compact input file `through48_pointings.bin.zst` occupies 6,228,872 bytes.
Its decompressed header is the little-endian structure `<8sIQQI>`:
the magic `UTPDOM1` followed by a zero byte, version 1, the source count,
the record count, and the maximum protocol length. Each record is
`<IBi>`: source index, pointing-case byte, and signed origin. The five
case bytes 0 through 4 mean translate, remove the origin, add the origin,
embed in an affine hyperplane, and embed then adjoin zero. The last two
cases do not use the origin field. `code/protocol_domain.py` reconstructs
the full mathematical input, without expanded matrices in the data file.

Every pointed support has been independently reconstructed from the final
compact catalogue. A fresh complete origin-reduction replay also agrees
on every finite point set; its certificate is
`certificates/through48_origin_replay.json`. The replay compares complete
sets within each source, rather than assuming a particular internal
representative order.

## Logical census

The distance-three logical quotient has dimension at most 11 on 4,191,933
inputs. The remaining inputs are explicitly listed, with their individual
quotient dimensions, in `through48_logical_censuses.json`:

| Quotient dimension | Supports |
| ---: | ---: |
| 12 | 1,387 |
| 13 | 1,177 |
| 14 | 149 |
| 15 | 35 |
| 16 | 21 |
| 17 | 10 |
| 20 | 1 |
| 21 | 3 |

Each of these 2,783 quotient dimensions is independently checked using
the binary equations in `code/logical_spaces.py`. The finite data contain
36 census partitions, covering logical dimensions 1 through 7 and the
necessary dimension-eight successor domain. Quotient dimensions at most
14 use direct common-isotropic-subspace enumeration. Larger quotients use
complete marked-subspace orbits. The dimension-one census combines these
two constructions.

The complete 4,194,716-input quotient histogram has also been freshly
recomputed. `certificates/through48_quotient_partition.json` agrees on
every listed larger-quotient count and on the full smaller-quotient
complement, so the selected domain omits no larger quotient.

The data retain exact raw counts where available and explicitly distinguish
them from marked-orbit counts. They also contain 125 local frontier
witnesses. Independent checks of every witness's overlap conditions,
rank, projectivity, exact distance, leading error coefficient and output
equivalence agree. Merging their `(n,S)` metrics at fixed output and exact
distance gives precisely the 52 released frontier points with `n<=48`.

All dimension-seven outputs in the smaller-quotient sector are absent.
The larger-quotient sector has two output types, with 118 and 124
nondegenerate hyperplanes respectively. In particular neither is primitive.
The complete dimension-eight census on that sector is empty. The primitive
tensor theorem in `logical_dimension_closure.md` therefore excludes every
logical dimension at least eight on the main source sector. These finite
absence results are essential premises, not consequences of witness
validity alone.

## The two cubic sources

For each of the two separately listed sources, every origin is considered
directly. There are 177 distinct pointed supports of protocol length at
most 48 per source. Exact distance-three outputs with `q<=4` are strictly
dominated: every such intrinsic output has a witness with `n<=47` and
`S<=q+6`, while these source supports have `n>=47` and `S>=q+7`.
The complete small-output orbit authority supplies all 22 intrinsic
classes, not merely the classes observed on these two sources.

Fresh direct `q=5`, distance-three censuses contain respectively 204,713
and 18,384 common-isotropic subspaces, all with degenerate logical tensor.
Thus no intrinsic `q=5` output occurs. Every higher-dimensional
nondegenerate tensor either descends through a nondegenerate hyperplane
or is an odd-dimensional primitive tensor. The latter contains a primitive
five-dimensional restriction when its dimension exceeds five. Hence all
intrinsic `q>=5` outputs are excluded.

At distance at least four, all filtered quotients vanish except on the
single affine-hyperplane embedding of each source. These embeddings
contain the constant stabiliser row. Consequently their logical tensors
are alternating: singles and repeated-index pair coefficients vanish.
Nonzero intrinsic outputs in dimensions one and two are then impossible.
Fresh direct dimension-three censuses contain 3,347 and 811 isotropic
subspaces respectively, all degenerate. Any higher nondegenerate
alternating tensor descends through hyperplanes or has the primitive
form `ell wedge omega`, which contains a nonzero three-dimensional
restriction. Therefore all higher-distance outputs are excluded too.

`certificates/through48_small_sources.json` records these full-origin and
logical-census checks. Adding a zero column gives no further Pareto point:
at distance at least three the complete added column is zero, and its
deletion preserves the output, exact distance and footprint.

## Reproduction and evidence boundary

After building the native kernel as described in the main README:

```text
python -B code/verify_through48_sector.py --full-domain --output checks.json
python -B code/replay_protocol_origins.py --native build/utsp-native --work-directory origin_check --workers 4 --output origins.json
python -B code/verify_protocol_quotient_partition.py --native build/utsp-native --work-directory quotient_check --workers 4 --output quotients.json
python -B code/verify_small_source_closure.py --native build/utsp-native --work-directory cubic_check --workers 4 --output cubic.json
python -B code/replay_through48_branch.py --branch 31 --native build/utsp-native --work-directory branch_check --workers 4 --output branch.json
```

Branch indices 0 through 35 address every retained logical census.
The last command illustrates one complete partition, not the whole
classification. Its independent replay recovers 247 isotropic
seven-dimensional subspaces on the 1,387 dimension-12 quotients, all
degenerate. The release retains the remaining complete census results and
code for reproducing them; it does not describe all those large finite
censuses as freshly re-enumerated during verification.

The above closes the protocol source domain conditional on completeness
of its parent-space catalogues and the recorded complete finite censuses.
It does not by itself certify the higher-dimensional parent-space
classification or protocol lengths 49 through 54.
