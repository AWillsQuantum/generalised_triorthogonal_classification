# Classification data through length 54

This is the **data and evidence delivery** for the classification of binary
generalised triorthogonal magic-state distillation protocols. Its companion
GitHub repository contains all theory, source code, tests and independent
readers.
Neither delivery requires access to the private development repository.

The data delivery is prepared for Figshare. Its dataset DOI is pending.

## Authorship and AI assistance

Developed by Adam Wills with substantial assistance from ChatGPT (OpenAI) in
implementation, testing and computational analysis. Adam Wills directed the
research and takes responsibility for the released results.

Adam Wills is the citation author of this dataset; ChatGPT is acknowledged
here as an AI development assistant. The companion
[GitHub repository](https://github.com/AWillsQuantum/generalised_triorthogonal_classification)
provides the software and its `CITATION.cff`.

## Licence

The original data and documentation are licensed under
[Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).
See [LICENSE_DATA.md](LICENSE_DATA.md). This includes the contents of both
data archives and the identical data and documentation in the companion
GitHub delivery. Software and third-party components are excluded from this
licence grant and retain their own terms.

## Contents

- `pareto_frontier.json`: the complete, directly readable 74-point protocol
  catalogue, identical to GitHub's `data/protocols/pareto_frontier.json`.
- `space_catalogues.zip`: all compact unital-space shards and their index.
- `classification_evidence.zip`: finite mathematical censuses, contraction
  data, output profiles, seeds and all original certificates.
- `DATA_FORMAT.md`: definitions, encodings, example expansion and audit scope.
- `LICENSE_DATA.md`: the CC BY 4.0 data and documentation licence notice.
- `ORIGINAL_MANIFEST.json`: every original release file's size and SHA-256.
- `delivery_layout.json`: maps original files to this dataset and GitHub.
- `MANIFEST.json`: checksums of this delivery's own files (excluding itself).
- `DELIVERY_AUDIT.json`: retained results of the packaging and audit checks.

The ZIP archives use ZIP64. Already-compressed `.zst` members are stored
without recompression;
text evidence is losslessly DEFLATE-compressed. Extracting the envelopes
leaves the large mathematical streams compressed.

## Scope and results

Protocols have `n <= 54`, exact `d_Z >= 3`, full row rank and pairwise distinct
complete columns. One zero column is permitted in the enumerated scope.
The frontier is partitioned by nonzero intrinsic output and **exact distance**;
within each partition it minimises length `n` and total matrix rows `S`.
Every logical row counts towards `S`. Compressed footprint is not considered.
Output equivalence is **CNOT+S**, not full Clifford equivalence. Logical
dimension is unrestricted, and the reductions include non-completable matrices.

There are **74 points for 62 output classes**: 67 at distance three, five at
distance four and two at distance five. Each JSON witness includes its full
matrix, exact distance, leading error coefficient and output information;
the `outputs` array provides the representative gates and their basis
certificates. The largest undominated protocol has 53 inputs.

The spaces comprise **301,029,259 affine-equivalence classes** through length
54, including decomposable spaces, of which 293,172,583 have length 54.
These counts are different from the protocol frontier counts. The matrices
are recovered from compact supports with the ambient dimension explicit.

## Read or verify

The frontier can be read immediately with any JSON reader. The files in each
ZIP have their original relative paths (`data/...` or `certificates/...`).
For the space encoding and a small worked example, see `DATA_FORMAT.md`.

For complete verification, obtain the matching GitHub code delivery, place
all the files of this record in an otherwise empty directory, and run from
the GitHub checkout:

```sh
python -B code/delivery.py verify --data /path/to/figshare-download
python -B code/delivery.py assemble --data /path/to/figshare-download --destination /path/to/new-classification-audit
```

These two commands need only Python 3.11+ and no network access. Assembly
uses about 5.6 GB extra space and produces the byte-identical original release,
including its README and manifest. The destination must be new and outside
both deliveries. It does not expand all matrices or polynomials.

Install that release's `requirements.txt` in a virtual environment, then run
from the assembled directory:

```sh
python -B code/verify_release_integrity.py
python -B -m unittest discover -s tests -v
python -B code/verify_protocols.py
python -B code/verify_classification.py --output ../classification_check.json
```

Read the assembled `AUDIT_GUIDE.md` for native builds and regeneration.
Checksums establish preservation, not mathematical completeness. Witness
verification and composition of the retained finite censuses are distinct
from repeating every enumeration. Some large censuses are retained aggregate
results; full independent regeneration requires the supplied code and may
need substantial computation and temporary storage.

This dataset contains no source code. Obtain it from the companion GitHub
release; the layout's original manifest hash identifies the matching version.

The `delivery_layout.json` identifies these files as the `data` delivery;
verification is fully offline and does not depend on a hosting service.
