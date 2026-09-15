# Independent factors of output magic states

This presentation update concerns only the 62 intrinsic output classes on the
74-point frontier for n <= 54 and d_Z >= 3. It changes no protocol matrix,
identifier, metric or coverage statement. Equivalence remains CNOT+S, not full
Clifford equivalence. A tensor product of output states does not imply a tensor
product of smaller distillation protocols.

## Exact decomposition

Let tau be the full symmetric trilinear form on V = F_2^q defining the output,
including its repeated-index values. Its radical is zero for intrinsic outputs.
The centroid consists of linear maps A satisfying

    tau(Ax,y,z) = tau(x,Ay,z) = tau(x,y,Az)

for every x,y,z. These are linear equations in the q^2 entries of A. A centroid
element P with P^2=P splits V as im(P) + ker(P); mixed tensor values vanish.
Conversely, every independent splitting yields such a projector.

The centroid is a commutative algebra. For A,B in it, moving the maps between
the three arguments gives

    tau(ABx,y,z) = tau(Bx,Ay,z) = tau(x,Ay,Bz)
                = tau(Ax,y,Bz) = tau(Ax,By,z) = tau(BAx,y,z).

The zero radical implies AB=BA. Therefore A -> A^2+A is F_2-linear on the
centroid. Its kernel is the finite Boolean algebra of idempotents. Its atoms,
the primitive idempotents, are unique, mutually orthogonal and sum to the
identity. Their images are precisely the indecomposable component spaces.
An output equivalence conjugates the centroid and permutes these components.
Thus the multiset of component isomorphism classes is intrinsic.

This argument uses the trilinear tensor, not tau(x,x,x), which loses information
in characteristic two. It makes no characteristic-zero polynomial assumption.
See also Belitskii and Sergeichuk, *Congruence of multilinear forms*, Theorem 9,
https://arxiv.org/abs/0710.0834.

## Registered normal form

All component classes needed here already occur among the frontier's
indecomposable outputs. We retain their previously simplified gate strings and
established output IDs. Factors are ordered by dimension, then ID. Equal
factors are collected into tensor powers. Nontrivial factors use local qubit
indices; the flat expanded gate instead uses consecutive disjoint blocks.

`code/output_factorisation.py` solves the centroid equations, computes the
primitive projectors and matches their restrictions exactly against this fixed
registry. An unknown component is rejected, not assigned a basis-dependent
name. Gate minimisation within factors is heuristic and frozen; the factorisation
and the normal form relative to this registry are exact. We do not claim that
factorisation-first notation minimises the global number of gate factors.

Each updated output stores the invertible change of basis from its previous
gate and the composed change from its unchanged protocol-tensor coordinates.
Checking all 2^q phases gives

    phase_previous(Bx) = phase_factored(x) + phase_Clifford(x) mod 8.

The correction contains only even linear coefficients and quadratic coefficients
divisible by four, hence only S powers and CZ gates, within CNOT+S equivalence.

## Audit replay

The original hash-bound catalogue is retained in
`provenance/pareto_frontier_before_factorisation.json`. The current catalogue is
`data/protocols/pareto_frontier.json`. `code/verify_output_factorisation.py` checks
the exact phase identities, indecomposability, component ordering and that
every other field is unchanged. `code/verify_protocols.py` independently checks
the current labels against all 74 original matrices, their exact distances
and error coefficients. Delivery assembly retains the original catalogue bytes
for replay of the original classification evidence; the presentation verifier
connects that result to the current labels without rewriting historical hashes.

`delivery_tests/test_output_factorisation.py` also checks random invertible
basis changes for every output, explicit T/CS/CCZ products, the previously
hidden six-T output, and rejection of invalid or altered certificates.
