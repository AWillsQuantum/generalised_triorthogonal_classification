# Complete augmentation through admissible parents

Fix a stabiliser support `X` and its finite coordinate automorphism group
`H`. Let `F_i` be an `H`-invariant family of admissible logical subspaces of
dimension `i`. The families can include isotropy, distance constraints and
output restrictions, but their definitions must be explicit.

## Minimum-admissible-parent lemma

Assume every member of `F_(i+1)` has a hyperplane in `F_i`, and one
representative of every `H`-orbit in `F_i` is available. For each parent `P`:

1. Enumerate all extensions in `F_(i+1)` modulo the exact stabiliser `H_P`.
2. Reject a child only when it has a hyperplane in `F_i` with strictly
   smaller invariant than `P`.
3. Canonicalise retained children under `H` and remove exact duplicates.

This retains every `H`-orbit in `F_(i+1)`. To prove it, choose a
minimum-invariant admissible hyperplane of any desired child. Its orbit
has a retained parent representative, and the extension-vector orbit
enumeration reaches the corresponding child. No smaller admissible
hyperplane can reject it. Equal invariants can retain more than one copy,
which exact canonicalisation subsequently removes.

The invariant need not distinguish all orbits. In contrast, the final
canonical key must be exact. Comparing with a hyperplane outside `F_i`
does not satisfy this lemma, even if that hyperplane is isotropic and has
the right dimension.

## Seeded target families

Let `R` be a downward-closed restriction family of target tensors of dimension
`t`. For a nondegenerate seed dimension `s`, the family at `s<=i<=t` consists
of isotropic `i`-spaces whose tensor belongs to `R` and which contain a
nondegenerate `s`-space. Below `s` there is no seed restriction.

Every generated flag preserves its seed. Conversely, any space in this
family contains a flag through a suitable seed; every member of the flag
passes `R` by downward closure. Applying the lemma at each level gives
completeness, provided competing parents satisfy the same reachability
condition. At dimension `s` this condition is nondegeneracy; above it the
condition is existence of the nondegenerate seed restriction.

The target family alone does not justify arbitrary extensions beyond `t`.
There one needs a new parent-family definition, a branch-union induction,
or an exact reachability predicate.

## Unions of target branches

Several seed branches can collectively supply the complete parent family.
It is sufficient that their union contains every orbit in `F_i` and that
every branch uses the same eligible-parent predicate. A child belongs to
the branch containing its minimum eligible parent; it need not be retained
in every seed branch from which it could be generated.

For example, a tensor-level cover of nonprimitive five-dimensional outputs
can use a family `F_3` of three-dimensional restrictions, `F_4` of
nondegenerate four-spaces containing a member of `F_3`, and `F_5` of
nondegenerate five-spaces containing a member of `F_4`. A finite tensor
cover and this lemma give separate, composable proofs. Primitive
five-dimensional tensors require their own complete construction.

Likewise, a four-dimensional hitting set defines `F_4` as the union of the
selected tensor orbits. A child comparison must be restricted to that
union, not to every four-dimensional hyperplane.

## Hereditary nondegenerate extension

Every nonprimitive nondegenerate tensor has a nondegenerate hyperplane.
Once the complete predecessor union is available, the lemma applies using
nondegenerate hyperplanes as the common parent family. Primitive tensors
are the exceptional branches described in `logical_dimension_closure.md`.
The predecessors here are all required subspace orbits, not merely the
predecessors that are Pareto optimal as protocols.

## Exact group actions and caching

When the action of `H` on a quotient has a kernel `K`, stabiliser orders
in the full group are `|K|` times those in the quotient image. The image can
be enumerated explicitly when small; its full closure, not only the given
generators, defines orbit membership. Nonfaithful actions must retain the
kernel factor.

An invariant cache is valid only for the full marked object determining the
eligible parent family and group action. Reusing a value across different
stabiliser or seed contexts requires a proof of invariance under that change.

The native tests check positive witness membership, independent primitive
constructions, canonical parent incidence, nonfaithful quotient actions and
explicit groups against graph canonicalisation. The lemma explains the
coverage condition that those implementation checks must preserve.
