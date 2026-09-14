# Classification of the unital parent spaces

This section concerns positive-length projective unital triorthogonal binary
spaces, including decomposable spaces. Write `r=dim(H)` and `m=r-1`.
After choosing the unit row first, columns are `(1,y)` for distinct points
of an affinely spanning support `Y subset F_2^m`.

## 1. Indicators, equivalence and the minimum length

The overlap equations say that the indicator `1_Y` is orthogonal to every
square-free monomial of degree at most three. Reed-Muller duality gives
`1_Y in RM(m-4,m)`. A nonzero such word has weight at least 16. In
particular there are no positive-length spaces below length 16; every
length is even. For `m<4` the relevant dual is zero. The duality and minimum
weight facts are the standard Reed-Muller statements used in the
[Nezami-Haah correspondence](https://arxiv.org/abs/2107.09684).

A linear row-space isomorphism fixing the distinguished all-one row acts
on the other coordinates as an invertible affine map. Thus coordinate
equivalence of the spaces is exactly affine equivalence of the supports.
This includes decomposable row spaces: indecomposability is not an input
restriction.

## 2. A dimension bound independent of the census

Let `D=H^2` be the span of coordinatewise products of pairs of words of H.
Triorthogonality gives `D subset H^perp`. Define the multiplier algebra
`A={a in F_2^c: aD subset D}`. It contains the unit and is closed under
coordinatewise multiplication. Every such binary algebra has a basis of
disjoint nonempty block indicators `e_1,...,e_t`; consequently `dim(A)=t`.
One can see this by identifying coordinates on which every element of A
agrees, then taking products of separating elements and their complements
to recover the indicators of the resulting classes.

The projection of H to each block is itself projective, unital and
triorthogonal. Indeed, for u,v,w in H,

```text
sum_j (e_i u)_j(e_i v)_j(e_i w)_j = <e_i(uv),w> = 0,
```

because `uv in D`, `e_i D subset D` and `D` is orthogonal to H. The unit
and distinct-column properties persist under restriction to a coordinate
block. Each block therefore has at least 16 coordinates and
`t<=floor(c/16)`.

The coding Kneser inequality gives `dim(D)>=2r-t`; see Theorem 3.3 of
[Mirandola and Zemor, Critical pairs for the Product Singleton Bound](https://arxiv.org/pdf/1501.06419).
Combining it with `dim(D)<=c-r` proves

```text
ceil(log_2 c) <= m <= floor((c+floor(c/16))/3)-1.
```

This bound requires neither a shorter catalogue nor an assumption that H
is indecomposable. It gives `m<=18` at c=54. The sharper c=54 bound
`m<=17` follows as below: three multiplier blocks would each have length
at least16, and the possible shorter lengths force a total of either48
or at least56, not54. Indeed, lengths18,20 and22 have m<=6 by the general
bound and are absent from the affine-linear and quadratic Reed-Muller
weight spectra. Thus t<=2.

Multiplier blocks of D need not be direct summands of H. Only the
projections are used above. Direct decompositions of H itself are a
separate invariant and can be checked by column-matroid components or by
the multiplier algebra of H.

## 3. An exhaustive inductive construction

There is a dimension-generic recursion ordered first by support length and
then by affine dimension. A support other than the full affine cube has a
strictly lighter affine hyperplane section: otherwise every nonconstant
Fourier coefficient of its indicator would vanish, forcing a constant
indicator. The full cubes at lengths16 and32 are included directly.

In coordinates `(x,t)`, write the indicator as `g(x)+t h(x)`. The
sections have supports G and `G triangle K`, where `K=supp(h)`, and

```text
g in RM(m-4,m-1),  h in RM(m-5,m-1),
c = 2|G|+|K|-2|G intersect K|.
```

In the (m-1)-dimensional section ambient space, G is orthogonal to
degree-at-most-two evaluations, whereas K is orthogonal to degree-at-most-
three evaluations. These are different constraints. For example, an
affine three-cube is an admissible light-section support despite having a
nonzero cubic moment.

Choose G to be the light section, so `2|G|<c`. The cylinder
`G x F_2` is a unital support of this shorter length. Thus every possible
G is obtained from a shorter unital space that admits a translation
period, by quotienting along that period. The derivative core K has
length at most c; if its length is c its affine dimension is lower than
that of the target. All its possible affine embeddings, including ones
in a proper subspace of the quotient ambient space, must be retained.

Enumerate the relative affine embeddings of G and K, impose the displayed
intersection equation, reconstruct the two sections, and retain exactly
the full-affine-rank supports with the required moments. Quotient by the
full affine group, not merely the group preserving the selected section.
The light-section argument proves coverage, and the lexicographic order
proves that the recursion terminates in lower cases.

## 4. The inverse-contraction acceleration

For a nonzero direction a, let N_a count pairs `{y,y+a}` contained in Y.
Each unordered pair has one difference, so

```text
sum_(a!=0) N_a = binomial(c,2).
```

There is therefore a direction with
`N_a<=floor(binomial(c,2)/(2^m-1))`. Project along it. Write C for the
singleton fibres and F for the complete fibres; these sets are disjoint,
`|C|=c-2N_a`, and `|F|=N_a`. The support C inherits all moments of degree
at most three, since complete fibres contribute twice. It is a shorter
unital support, or a same-length one in lower dimension when N_a=0.

After choosing the last coordinate along a, every lift has the form

```text
{(x,f(x)):x in C} union {(z,0),(z,1):z in F}.
```

Let phi_2(x) contain the values of all square-free monomials of degree at
most two, including the constant. The remaining moment equations are
exactly the linear system

```text
sum_(x in C) f(x) phi_2(x) = sum_(z in F) phi_2(z).
```

Adding an affine function to f is an affine shear of the lift. Once a
particular solution is found, the homogeneous solutions can therefore be
quotiented by the affine evaluation space. This uses only exact binary
elimination. The admissible F are a fixed-cardinality zero-sum problem in
the quotient of the feature space by the span of `phi_2(C)`, allowing
dynamic programming or meet-in-the-middle enumeration.

Do not assume that C spans the whole quotient ambient space when N_a>0:
complete fibres may supply additional affine directions. A full-rank
check on each lift, or a proved deletion-kernel bound, is required.

At c=54 the complete direction cover is particularly small:

| target m | multiplicities | core lengths |
|---:|:---|:---|
|8|1,3,5|52,48,44|
|9|0,1,2|54,52,50|
|10|0,1|54,52|
|11 through17|0|54|

For m=8, RM(3,7) weight divisibility by four forces N_a to be odd.
The low-dimensional source counts are188,98 and35, respectively, for
lengths52,48 and44 in affine dimension seven.

These particular cores necessarily span seven dimensions. At lengths44
and52, the complete census in dimensions at most six is empty. At
length48 a proper-span core would require a nonzero affine function
vanishing on the singleton fibres. Its values on the three full fibres
would be an even binary word; a word of weight two contradicts the
linear-coordinate moment equations, since the two fibres are distinct.
Thus no such function exists. More generally this argument rules out
rank loss for N<=3. Once a core spans seven dimensions and N>0, every
lift has full affine dimension eight: an affine equation on the lift
vanishes in the last coordinate because a full fibre contains both values,
and then vanishes identically because the core spans the quotient.

`code/prepare_space_contractions.py` reconstructs all321native inputs from
the compact catalogues and their complete affine stabilisers. It checks
the compatible-fibre counts by two independent methods. This certifies
the input domain, not the subsequent child enumeration or affine quotient.

`alternative_contraction_cover.md` proves a finite cover of four core
families by alternative directions and gives their complete explicit
witness streams. `length54_m08_census.md` describes the primary census,
its exact affine quotient, and the distinction between aggregate census
checks and a fresh full re-enumeration.

## 5. Zero-fibre closure

When `2^m-1>binomial(c,2)`, some direction has no complete fibre.
Projection is then injective. If the target spans m dimensions, its core
spans m-1 and every lift is a Boolean graph. For a core C of affine
dimension s, valid graph functions modulo affine shears have dimension

```text
c-rank(E_2(C))-(s+1).
```

The zero coset is precisely the rank-deficient affine graph and is
discarded. Thus an all-zero lift-quotient profile proves the next level
empty. This, followed by empty-predecessor induction, closes the top
dimensions without an open-ended dimension search.

## 6. Exact affine quotient and evidence requirements

Moment ranks, hyperplane spectra and direction-multiplicity profiles are
invariants, not canonical labels. Unequal invariants prove inequivalence;
equal invariants only identify candidates for an exact test. A positive
affine equivalence is certified by a translation and an invertible binary
basis map taking one support exactly onto the other. A negative decision
requires an exhaustive exact transporter test or a complete canonical
labelling algorithm.

For each inductive sector, a completeness certificate must identify the
entire parent domain and its lift equations, account for every admissible
lift modulo the stated groups, and bind the exact affine quotient to the
distributed representatives. Counts or hashes alone do not prove that
domain coverage. The release separates these obligations from the
independent rank, moment and lossless-encoding checks of each space.
