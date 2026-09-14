# Generalised triorthogonal protocols through length 54

This is the **code and theory delivery** of the classification. The companion
[Zenodo dataset (DOI: 10.5281/zenodo.22740399)](https://doi.org/10.5281/zenodo.22740399)
supplies the compressed space catalogues, finite censuses and
certificates. Neither delivery needs the original development repository.

## Authorship and AI assistance

Developed by Adam Wills with substantial assistance from ChatGPT (OpenAI) in
implementation, testing and computational analysis. Adam Wills directed the
research and takes responsibility for the released results.

Please use [CITATION.cff](CITATION.cff) to cite the software. Adam Wills is the
citation author; ChatGPT is acknowledged here as an AI development assistant.

## Licence

The original data and documentation in both deliveries are licensed under
[Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).
See [LICENSE_DATA.md](LICENSE_DATA.md) for scope, including the data inside
the Zenodo archives and this repository's copy of the protocol catalogue.
This does not assign a licence to the project's software or change any
third-party terms.

## Mathematical scope

The protocol matrices have full row rank and pairwise distinct complete columns;
one zero column is permitted in the enumerated domain. The classification covers
`n <= 54` and exact distance `d_Z >= 3`, including non-completable matrices.
Outputs are identified under **CNOT+S equivalence**, not full Clifford equivalence.
Logical dimension is unrestricted; only nonzero intrinsic magic outputs are
keys. At fixed output and **exact distance**, the frontier minimises `(n,S)`,
where `S` is the total number of matrix rows, including every logical row.
Compressed footprint is not an objective.

There are **74 Pareto points for 62 output classes**: 67 at distance three,
five at distance four and two at distance five. The largest undominated length
is 53. These are frontier counts, not counts of all protocol matrices.

The complete catalogue is `data/protocols/pareto_frontier.json`. It contains
matrices, exact distances, leading error coefficients, representative gates and
output-basis certificates. Gate representatives are heuristically simplified,
not certified minimal. The Zenodo file `pareto_frontier.json` is byte-identical.

The space data comprise 301,029,259 affine-equivalence classes through length
54, including decomposable spaces; 293,172,583 have length 54. The small
`data/spaces/index.json` describes the complete dataset, **not** files all present
in this checkout. Only a length-16 example is included here. Obtain the other
shards from Zenodo.

## Use GitHub alone

Use Python 3.11 or later. The following checks use only the standard library:

```sh
python -B code/verify_release_integrity.py --allow-extra
python -B code/verify_protocols.py
python -B -m unittest discover -s delivery_tests -v
```

`--allow-extra` permits Git's `.git` directory and local build files; it still
checks every distributed file. Do not recreate the manifest to hide a mismatch.
The protocol check independently verifies all 74 witnesses and their output
equivalence separation; it does not establish exhaustion of the search domain.

For the small space example, first install `requirements.txt` in a virtual
environment (`python -m pip install -r requirements.txt`), then run:

```sh
python -B code/read_spaces.py data/spaces/c16/m04/spaces_000000000.utspace.zst --limit 1 --matrix --polynomial --validate
```

## Full audit using both deliveries

Download **all files** in the matching Zenodo record to an otherwise empty
directory, here called `../zenodo-download`. Keep both ZIP files intact. From
the GitHub checkout, run:

```sh
python -B code/delivery.py verify --zenodo ../zenodo-download
python -B code/delivery.py assemble --zenodo ../zenodo-download --destination ../classification-audit
cd ../classification-audit
python -B code/verify_release_integrity.py
python -B -m unittest discover -s tests -v
python -B code/verify_classification.py --output ../classification_check.json
```

The destination must not exist and must be outside both deliveries. Assembly
recovers the **exact original release**, including its original manifest and
README. It extracts the ZIP envelopes but leaves `.zst` streams compressed;
it does not expand millions of matrices. Allow about 5.6 GB of additional disk
space for assembly. Regenerating the full classification can require much more
temporary storage and computation.

The composition verifier joins the finite space induction, protocol source
partitions, logical successor closure and independently checked frontier.
It verifies the retained census bindings, **not** a fresh execution of every
large enumeration. Read `AUDIT_GUIDE.md` and the assembled README for the full
audit and native reproduction commands. Checksums and valid witnesses alone
are not completeness proofs. `AUDIT_STATUS.json` describes the original
classification audit; `DELIVERY_AUDIT.json` describes this packaging check.

## Organisation

- `theory/`: definitions, proofs and complete finite-domain specifications.
- `code/`: all original mathematical algorithms, readers and verifiers.
- `code/native/`, `code/space_native/`: accelerated C++20 implementations.
- `tests/`: original mathematical tests; run these after full assembly.
- `delivery_tests/`: standalone packaging and corruption-rejection tests.
- `resource_estimates/`: the conditional `10^8` core-hour length-56 forecast.
- `provenance/`: original README and original file manifest, retained unchanged.
- `delivery_layout.json`: every original path and its location in the deliveries.
- `packaging/`: reproducible delivery builder and Zenodo documentation templates.

Start the theory with `theory/support_correspondence.md` and
`theory/space_classification.md`. The native build requires a POSIX environment,
`g++` with C++20/OpenMP, `ar`, and GMP development files. Python algorithms use
NumPy and `zstandard`. Bundled Bliss source retains its own notices; see
`THIRD_PARTY_NOTICES.md`. A distribution licence for the project's original
software has not yet been assigned; data and documentation are covered by
`LICENSE_DATA.md`, and third-party components retain their own terms.

## Rebuild the delivery layout

After reconstructing the original release, the builder can reproduce its split
without modifying it. From this checkout:

```sh
python -B packaging/build_deliveries.py --source ../classification-audit --destination ../rebuilt-deliveries --templates .
```

The builder generates new package manifests. It does not copy the later audit
receipt, which must be regenerated by testing the new deliveries. Inner data
streams and original files are preserved byte-for-byte; ZIP bytes may depend
on the Python/zlib version. No command here uploads anything.
