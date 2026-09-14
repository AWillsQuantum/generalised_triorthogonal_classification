# Cubic Seed Completeness

The affine-dimension-seven seeds are represented by cubic Boolean
polynomials. `data/seeds/rm37.json` contains the polynomial representatives
from [Gillot and Langevin's affine classification data](https://langevin.univ-tln.fr/project/agl7/aglclass.html).
The following independent test establishes their completeness, without
assuming the number of representatives or any supplied stabiliser order.

For a support in `F_2^m`, construct the incidence graph between all `2^m`
points and all `2(2^m-1)` affine hyperplanes. Colour support points,
complement points and hyperplanes separately. For `m>=2`, an isomorphism
of these coloured graphs induces exactly an affine equivalence of supports.

Indeed, affine maps preserve hyperplanes. Conversely, preservation of
hyperplanes preserves their intersections and hence affine two-planes.
Four distinct points form an affine two-plane precisely when their sum
vanishes. Translating the image of zero to zero therefore makes the point
map additive, and thus linear. Hyperplanes are determined by their incident
point sets, so their permutation adds no kernel to the action on points.

The bundled Bliss library computes exact canonical forms and full graph
automorphism groups. GMP represents group orders exactly. Every returned
generator is independently checked on the complete affine point set. A
canonical affine frame supplies an explicit map from each canonical support
back to its input. `verify_affine_graph_small.py` compares this construction
with all 1,344 affine maps on all 256 three-variable Boolean functions.

For seven variables,

```text
|AGL(7,2)| = 128 product_{i=0}^6 (128-2^i)
          = 20,972,799,094,947,840.
```

Each polynomial is explicitly checked to have square-free degree at most
three. Distinct exact canonical supports give disjoint affine orbits. The
sum of their orbit sizes is

```text
sum |AGL(7,2)| / |Stab(f)| = 2^64
                           = |RM(3,7)|.
```

Thus the 3,486 orbits exhaust the cubic code. Restricting to full affine
rank gives, respectively, 35, 98 and 188 classes at weights 44, 48 and 52.
There is no weight-54 cubic function. The support correspondence identifies
these classes with the affine-dimension-seven unital spaces. Completeness
in larger affine dimensions is a separate inverse-contraction obligation.

To reproduce the finite certificate after building the native code:

```text
python code/verify_affine_graph_small.py --native build/native/affine_support_graph
python code/verify_rm37.py --native build/native/affine_support_graph --output-directory build/rm37
```

The resulting `orbits.json` records exact stabiliser orders, canonical
support masks and affine maps; `certificate.json` checks disjointness and
the full orbit mass. Neither a file hash nor a representative count is
used as a substitute for this mathematical exhaustion test.

## Lower Dimensions and Catalogue Binding

For `m=4`, the indicator is constant. For `m=5`, it is affine linear.
For `m=6`, it is quadratic. The affine quadratic normal forms are
`Q_r+c`, where `Q_r=x1*x2+...+x(2r-1)*x(2r)` and `c` is zero or one,
together with `Q_r+x(2r+1)` when `2r<m`. The zero polynomial and the affine
linear cases are included. This follows by putting the alternating polar
form in symplectic normal form. A linear term that vanishes on the radical can be
removed by translation. A nonzero restriction of the linear term to the
radical gives the displayed balanced form; translation removes its constant.

`verify_low_dimensional_spaces.py` checks these finite normal forms and the
complete cubic seed orbits against every released support of affine
dimension at most seven. It verifies exact affine maps, pairwise
inequivalence and equality of the two class sets, not just their counts.
Dimensions below four contain no nonempty unital projective support,
by the indicator-code correspondence. This gives a complete base for the
higher-dimensional induction.
