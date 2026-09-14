# Configured affine algorithms

`build.py` compiles the same signature and exact affine-quotient algorithms
for any even support length from 16 to 54 and affine dimension from 8 to 15.
`configuration.hpp` sets the bitmap width and record layout at compilation.
This range is a supported binary representation, not a claim that every
length/dimension pair contains spaces.

```sh
python -B code/space_native/generic/build.py --length 52 --dimension 10 --include-lifts --build-dir ../audit-build/spaces
```

The executables are `signature_c52_m10`, `exact_quotient_c52_m10` and,
with `--include-lifts`, `extensions_c52_m10`. Accelerated lift configurations
are available for lengths 50, 52, 54 and dimensions 9, 10, 11. Dimension
eleven also supplies `profile_c52_m11`. The minimum-direction filter is
enabled for the dimension-nine length-52 and length-54 lift routines and
disabled at length 50. Other dimensions use their stated full lift domains.

The complete generation and quotient interfaces, their binary records and
all source multiplicities are specified in
`theory/length50_length52_censuses.md`,
`theory/length54_middle_dimensions.md` and
`theory/finite_contraction_censuses.md`. Build reports certify compilation;
the separate differential and finite-domain checks certify mathematical
behaviour.
