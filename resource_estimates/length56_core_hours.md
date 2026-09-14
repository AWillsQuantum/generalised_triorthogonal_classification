# Resource estimate for extending the protocol classification to length 56

## Estimate and scope

**The single order-of-magnitude planning estimate is 10^8 core-hours.**
This estimate was reassessed on 12 September 2026. It concerns completing the
Pareto classification from protocol length at most 54 to length at most 56,
using the present space-classification and compatible-logical-row methods.
The scope is the same as the release: full-row-rank matrices with distinct
complete columns, exact distance at least three, unrestricted intrinsic
logical outputs up to CNOT+S equivalence, and objectives `(n,S)` at fixed
output and exact distance. It is not a classification restricted to five
logical qubits. Already completed work through length 54 need not be repeated.

This is an **engineering forecast**, not a measured total, a complexity
bound, a completion guarantee or a statistical confidence interval. The
unmeasured parts are substantial. A representative new pilot could change
even the exponent. Core-hours measure CPU work, not queueing or transfer
time, and do not by themselves give a calendar ETA.

## Why the earlier extrapolation is not retained

The earlier feasibility calculation projected approximately 2.73 billion
core-hours for the dominant logical-row filter alone. Its expensive tail
used 16 completed timing cases: twelve with candidate-quotient dimension 15
and four with dimension 16. Their mean costs were respectively 29.95 and
48.07 core-hours. These are dimensions of the space of candidate logical
rows after quotienting by stabiliser rows, not logical output dimensions.
The expensive search in this measurement asks for five-dimensional logical
subspaces.

The observed hard cases were concentrated in seven sampled lifts of only
two length-48 cores. Those cores contribute respectively **199** and
**126,244** direction-marked extension classes, whereas the historical
`m=8, N=4` population used for calibration contains **14,453,126,733**.
Here `m` is the affine dimension of the length-56 parent support and `N`
is the number of doubly occupied pairs in its directional contraction.
The previous extrapolation averaged the selected core rates equally before
multiplying by the entire route population. That is not population weighting
when the selected cores have such different extension counts.

The revised calculation weights each sampled core by its own marked-class
count, and extrapolates within quadratic-nullity strata. It also distinguishes
raw lifts, direction-marked equivalence classes and origin orbits. The change
is an assessment of the forecast, not a correction to any completed
classification or to the measured timing results.

## Reproducible calculation

The retained input file contains all 256 sampled `N=4` parent supports,
their recomputed affine symmetry and quotient-dimension profiles, the
stratum populations and all 16 hard-key timings. The profiles were recomputed
locally; no new large enumeration or cluster benchmark was used.

Let `C` be a seven-dimensional contracted core and `M_C` its number of
direction-marked extension classes. For a sampled lift `(Y,a)`, let:

- `A(Y,a)` be the order of the affine automorphism group preserving both
  the support `Y` and its distinguished direction `a`;
- `d(Y)` be the number of affine-automorphism orbits of admissible directions;
- `p_k(Y)` be the number of origin-orbit constructions with candidate
  quotient dimension `k`, including the construction retaining the entire
  unital space as stabilisers when it has that dimension.

An admissible direction has at most six doubly occupied pairs. All such
directions in these retained samples have a full-affine-rank contracted core,
as required by the calibrated core census. Equivalent origins are counted
once. An origin in the support gives the length-55 construction with the
zero column deleted; retaining that redundant zero column cannot improve
the frontier at the same output, distance and footprint.

For a fixed core, a raw lift orbit has size proportional to `1/A(Y,a)`.
Consequently, the self-normalised, inverse-orbit-size weighted estimate of
its contribution is

```text
estimated contribution of C to dimension k
  = M_C * sum_i [ A(Y_i,a) * p_k(Y_i) / d(Y_i) ]
        / sum_i A(Y_i,a).
```

The denominator `d(Y)` prevents counting the same parent once for every
admissible direction orbit. In an exhaustive sum over marked classes,
`sum_(Y,a) p_k(Y)/d(Y) = sum_Y p_k(Y)`. The finite sampled ratio above is
only an estimator of that sum. In particular, the physical number of
directions is not interchangeable with their number of automorphism orbits.

The structural strata use `dim I_2(C)`, where `I_2(C)` is the vector space
of polynomials of degree at most two vanishing on the contracted support.
The historical populations used in the extrapolation are:

| `dim I_2(C)` | Core classes | Direction-marked extension classes | Sampled cores |
|---:|---:|---:|---:|
|0|59|14,294,858,292|12|
|1|22|157,332,063|2|
|2|7|660,046|1|
|3|3|248,111|0|
|4|3|26,725|0|
|5|2|1,496|1|

Within a sampled stratum, contributions are summed and divided by the
sampled marked-class mass before scaling to the full stratum mass.
For the two unsampled strata, the calculation uses the largest observed
hard-support rate per marked class, separately for dimensions 15 and 16.
This choice is an explicit planning assumption, not an upper bound.

The resulting model has approximately 1.03 million dimension-15 hard
supports and 76,800 dimension-16 hard supports. Multiplying by the measured
mean CPU costs gives the following accounting. Extra digits are retained
in the machine-readable arithmetic only for reproducibility.

| Component | Core-hours |
|---|---:|
| Reweighted hard logical-row filter | 34.5 million |
| Ordinary `m=8` and initial `m=9,10` filter allowance | 2.15 million |
| Additional work allowance | 36.7 million |
| Unrounded planning total | 73.4 million |
| **Single reported order of magnitude** | **10^8** |

The ordinary-filter allowance carries forward the earlier coarse projection;
it is not a newly certified source count. The additional allowance equals
the modelled filter subtotal. It is intended to cover source generation and
affine classification, other logical dimensions and higher distances,
primitive branches introduced by the completeness repair, validation and
CPU work lost to retries. This factor of two is an engineering assumption,
not measured overhead. Rounding is to the nearest power of ten.

## Limitations

Only sixteen raw lifts per selected core were sampled, deterministically.
Inverse-orbit-size weighting can give a large influence to a single highly
symmetric lift. One selected core represents each observed hard stratum;
their timings and frequencies are not known to represent every other core.
The source census used for this forecast predates the cubic-seed repair.
The two additional full-affine-rank length-48 seeds both have quadratic
nullity zero, but that does not establish that their extensions are cheap.
Historical calibration populations are not asserted to be the current
complete length-56 space catalogue.

The base model observes no hard tail in the sampled nullity-zero/one cores
or the other contraction routes. This is not a proof of absence and must
not be used to omit these routes from a classification. Unobserved hard
families, higher quotient dimensions, and expensive logical successors
could exceed the additional-work allowance. Conversely, further exact
deduplication or a better algorithm could reduce the work.

The proposed universal nonexistence result for the difficult five-logical-row
branch was disproved by explicit positive cases. No such theorem is assumed
here. Later faster profile-presence runs are not substituted for the older
full-enumeration timings: those operations need not have the same cost.

## Reproduction

From the release root, using only the Python standard library:

```sh
python -B code/estimate_length56_core_hours.py
```

This reproduces `length56_estimate.json` from `length56_inputs.json`.
The input records identify the historical measurement and sampling files
by path and SHA-256. Those larger historical campaign directories are not
needed for this arithmetic; this small package is not a replacement for
their raw timing logs or a new certificate of their population counts.

For an additional recheck of the support-profile calculations,
install `igraph` (tested with version 0.11.9) in an external environment and run:

```sh
python -B code/estimate_length56_core_hours.py --check-sample
```

This reconstructs the full binary row space on each support, computes its
affine automorphisms using a coloured incidence graph, enumerates direction
and origin orbits, and checks the stored quotient dimensions. Neither
command launches a protocol classification. The resource forecast is
separate from the completed classification certificates in this release.
