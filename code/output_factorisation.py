"""Exact independent factors of intrinsic binary output tensors.

The centroid and its idempotents are computed over F_2, using the full
trilinear tensor (not the polynomial obtained by putting its arguments equal).
Gate labels are a presentation layer; protocol matrices and class IDs stay fixed.
"""

from collections import Counter
from copy import deepcopy
from itertools import combinations, groupby, product
import re

from protocol_checks import LogicalTensor, find_equivalence_basis, gate_tensor, verify_basis
from space_codec import binary_rank


def combine(vector, basis):
    result = 0
    for i, value in enumerate(basis):
        if vector >> i & 1:
            result ^= value
    return result


def reduced_basis(rows):
    pivots = {}
    for value in rows:
        for pivot in sorted(pivots, reverse=True):
            if value >> pivot & 1:
                value ^= pivots[pivot]
        if value:
            pivot = value.bit_length() - 1
            for other in pivots:
                if pivots[other] >> pivot & 1:
                    pivots[other] ^= value
            pivots[pivot] = value
    return [pivots[p] for p in sorted(pivots, reverse=True)]


def nullspace(rows, width):
    reduced = reduced_basis(rows)
    pivots = {row.bit_length() - 1 for row in reduced}
    result = []
    for free in range(width):
        if free in pivots:
            continue
        vector = 1 << free
        for row in reduced:
            if (row & vector).bit_count() & 1:
                vector ^= 1 << (row.bit_length() - 1)
        result.append(vector)
    return result


def matrix_columns(matrix, q):
    return [(matrix >> (i * q)) & ((1 << q) - 1) for i in range(q)]


def matrix_product(left, right, q):
    columns = matrix_columns(left, q)
    return sum(combine(value, columns) << (i * q)
               for i, value in enumerate(matrix_columns(right, q)))


def centroid_basis(tensor):
    q = tensor.q
    if not tensor.intrinsic():
        raise ValueError("Factorisation requires an intrinsic (radical-free) output")
    values = [[[tensor.value(1 << i, 1 << j, 1 << k)
                for k in range(q)] for j in range(q)] for i in range(q)]
    equations = []
    for i, j, k in product(range(q), repeat=3):
        first = sum(values[a][j][k] << (i*q+a) for a in range(q))
        second = sum(values[i][a][k] << (j*q+a) for a in range(q))
        third = sum(values[i][j][a] << (k*q+a) for a in range(q))
        equations.extend((first ^ second, first ^ third))
    return nullspace(equations, q*q)


def primitive_projectors(tensor):
    """Return every primitive centroid idempotent, in deterministic order."""
    q = tensor.q
    centroid = centroid_basis(tensor)
    for a, b in combinations(centroid, 2):
        if matrix_product(a, b, q) != matrix_product(b, a, q):
            raise ValueError("The intrinsic tensor's centroid is not commutative")
    # In a commutative F_2-algebra, A -> A^2 + A is linear. Its kernel is
    # the Boolean algebra of idempotents, with dimension equal to block count.
    images = [matrix_product(a, a, q) ^ a for a in centroid]
    equations = [sum(((a >> bit) & 1) << i for i, a in enumerate(images))
                 for bit in range(q*q)]
    boolean_basis = [combine(v, centroid) for v in nullspace(equations, len(centroid))]
    if not 1 <= len(boolean_basis) <= q:
        raise ValueError("Invalid idempotent algebra dimension")
    idempotents = [combine(v, boolean_basis) for v in range(1, 1 << len(boolean_basis))]
    primitive = []
    for a in idempotents:
        if matrix_product(a, a, q) != a:
            raise ValueError("Invalid centroid idempotent")
        if not any(b != a and matrix_product(a, b, q) == b for b in idempotents):
            primitive.append(a)
    identity = sum(1 << (i*q+i) for i in range(q))
    if combine((1 << len(primitive))-1, primitive) != identity:
        raise ValueError("Primitive projectors do not sum to the identity")
    if any(matrix_product(a, b, q) for a, b in combinations(primitive, 2)):
        raise ValueError("Primitive projectors are not orthogonal")
    return sorted(primitive, key=lambda a: (binary_rank(matrix_columns(a, q)), a))


def independent_blocks(tensor):
    projectors = primitive_projectors(tensor)
    if len(projectors) == 1:
        return [list(1 << i for i in range(tensor.q))]
    blocks = [reduced_basis(matrix_columns(a, tensor.q)) for a in projectors]
    flat = [v for block in blocks for v in block]
    if len(flat) != tensor.q or binary_rank(flat) != tensor.q:
        raise ValueError("The component spaces do not give a direct sum")
    for left, right in combinations(blocks, 2):
        if any(tensor.value(a, b, 1 << k) for a in left for b in right
               for k in range(tensor.q)):
            raise ValueError("Nonzero tensor overlap between components")
    return blocks


def restrict(tensor, basis):
    return LogicalTensor(tuple(tensor.words[v] for v in basis))


def compact_gate(tensor):
    terms = []
    for size, name in ((1, "T"), (2, "CS"), (3, "CCZ")):
        for term in combinations(range(tensor.q), size):
            arguments = term + (term[-1],) * (3-size)
            if tensor.value(*(1 << i for i in arguments)):
                terms.append(name + "".join(str(i+1) for i in term))
    return "".join(terms)


def gate_terms(gate):
    terms = re.findall(r"(CCZ|CS|T)([1-9]+)", gate)
    if not terms or "".join(name+indices for name, indices in terms) != gate:
        raise ValueError("Invalid compact gate expression")
    return terms


def expanded_gate(factors):
    terms, offset = [], 0
    for factor in factors:
        for name, indices in gate_terms(factor["representative_gate"]):
            terms.append((name, tuple(int(i)+offset for i in indices)))
        offset += factor["q"]
    return "".join(name + "".join(map(str, indices)) for name, indices in
                   sorted(terms, key=lambda t: (len(t[1]), t[1])))


def factor_label(factor, latex=False):
    gate = factor["representative_gate"]
    named = {(1, "T1"): "T", (2, "CS12"): "CS", (3, "CCZ123"): "CCZ"}
    if (factor["q"], gate) in named:
        name = named[factor["q"], gate]
        return r"\mathsf{" + name + "}" if latex else name
    if latex:
        return r"\allowbreak".join(r"\mathsf{" + name + "}_{" + indices + "}"
                                   for name, indices in gate_terms(gate))
    return gate


def factored_label(factors, latex=False):
    parts = []
    for _, group in groupby(factors, key=lambda f: (f["q"], f["representative_gate"])):
        copies = list(group)
        label = factor_label(copies[0], latex)
        composite = len(gate_terms(copies[0]["representative_gate"])) > 1
        if composite and (len(factors) > 1 or len(copies) > 1):
            label = (r"\bigl(" + label + r"\bigr)") if latex else "(" + label + ")"
        if len(copies) > 1:
            label += (r"^{\otimes " + str(len(copies)) + "}") if latex else "^(tensor " + str(len(copies)) + ")"
        parts.append(label)
    return (r"\otimes\allowbreak " if latex else " tensor ").join(parts)


def factor_registry(outputs):
    """Freeze the existing short representatives of indecomposable classes."""
    result = []
    for row in outputs:
        tensor = gate_tensor(row["q"], row["representative_gate"])
        if len(primitive_projectors(tensor)) == 1:
            result.append(({key: row[key] for key in ("q", "output_id", "representative_gate")}, tensor))
    return sorted(result, key=lambda pair: (pair[0]["q"], pair[0]["output_id"]))


def normal_form(tensor, registry):
    """Canonical factored labels relative to a fixed indecomposable registry.

    Unknown factors fail explicitly rather than receive basis-dependent labels.
    The returned basis is in the coordinates of the input tensor.
    """
    matched = []
    for block in independent_blocks(tensor):
        component = restrict(tensor, block)
        for label, target in registry:
            if target.q != component.q:
                continue
            basis = find_equivalence_basis(component, target)
            if basis is not None:
                matched.append((deepcopy(label), [combine(v, block) for v in basis]))
                break
        else:
            raise ValueError("Indecomposable factor absent from the representative registry")
    matched.sort(key=lambda pair: (pair[0]["q"], pair[0]["output_id"]))
    factors = [pair[0] for pair in matched]
    basis = [v for _, block in matched for v in block]
    gate = expanded_gate(factors)
    verify_basis(tensor, gate_tensor(tensor.q, gate), basis)
    return dict(q=tensor.q, representative_gate=gate, independent_factors=factors,
                factorised_representative_gate=factored_label(factors),
                factorised_representative_latex=factored_label(factors, True)), basis


def phase(q, gate, vector):
    return sum(1 << (len(indices)-1) for _, indices in gate_terms(gate)
               if all(vector >> (int(i)-1) & 1 for i in indices)) % 8


def clifford_correction(q, old_gate, new_gate, basis):
    verify_basis(gate_tensor(q, old_gate), gate_tensor(q, new_gate), basis)
    coefficients = [(phase(q, old_gate, combine(x, basis)) - phase(q, new_gate, x)) % 8
                    for x in range(1 << q)]
    for i in range(q):
        for x in range(1 << q):
            if x >> i & 1:
                coefficients[x] = (coefficients[x] - coefficients[x ^ (1 << i)]) % 8
    singles, pairs = [], []
    for mask, value in enumerate(coefficients):
        if not value:
            continue
        indices = [i+1 for i in range(q) if mask >> i & 1]
        if len(indices) == 1 and value % 2 == 0:
            singles.append([indices[0], value//2])
        elif len(indices) == 2 and value == 4:
            pairs.append(indices)
        else:
            raise ValueError("The full phase difference is not diagonal Clifford")
    return dict(S_powers=singles, CZ_pairs=pairs, truth_table_points_checked=1 << q)


def verify_factorisation(output):
    factors = output["independent_factors"]
    if not factors or sum(f["q"] for f in factors) != output["q"]:
        raise ValueError("Wrong factor dimensions")
    if factors != sorted(factors, key=lambda f: (f["q"], f["output_id"])):
        raise ValueError("Noncanonical factor ordering")
    for factor in factors:
        if len(primitive_projectors(gate_tensor(factor["q"], factor["representative_gate"]))) != 1:
            raise ValueError("A stated indecomposable factor splits further")
    if expanded_gate(factors) != output["representative_gate"]:
        raise ValueError("The factors do not expand to the displayed gate")
    if factored_label(factors) != output["factorised_representative_gate"]:
        raise ValueError("Incorrect factored label")
    if factored_label(factors, True) != output["factorised_representative_latex"]:
        raise ValueError("Incorrect factored LaTeX label")
    return True


def update_catalogue(original):
    result = deepcopy(original)
    registry = factor_registry(original["outputs"])
    for old, new in zip(original["outputs"], result["outputs"]):
        form, basis = normal_form(gate_tensor(old["q"], old["representative_gate"]), registry)
        new.update(form)
        new["gate_basis_in_tensor_coordinates"] = [combine(v, old["gate_basis_in_tensor_coordinates"])
                                                     for v in basis]
        new["factorisation_certificate"] = dict(basis_in_previous_gate_coordinates=basis,
            diagonal_clifford_correction=clifford_correction(old["q"], old["representative_gate"],
                                                             new["representative_gate"], basis))
    result["output_representative_presentation"] = dict(
        method="primitive-centroid-idempotents-v1",
        preferred_field="factorised_representative_gate",
        expanded_coordinate_field="representative_gate",
        factor_order="Increasing qubit count, then established output ID",
        factor_coordinates="Each tensor factor uses its own local qubit indices; the expanded gate uses consecutive disjoint blocks",
        scope="The 62 intrinsic output classes of the n<=54, d_Z>=3 frontier",
        minimum_gate_count_claim=False,
        description="Exact indecomposable output factors; fixed previously simplified representatives within factors. Protocol matrices and class IDs are unchanged.")
    verify_presentation_update(original, result)
    return result


PRESENTATION_FIELDS = ("factorised_representative_gate", "factorised_representative_latex",
                       "independent_factors", "factorisation_certificate")


def verify_presentation_update(original, updated):
    """Bind changed labels to the unchanged, hash-certified scientific catalogue."""
    stripped = deepcopy(updated)
    metadata = stripped.pop("output_representative_presentation")
    if metadata["method"] != "primitive-centroid-idempotents-v1":
        raise ValueError("Unknown output presentation method")
    if len(original["outputs"]) != len(updated["outputs"]):
        raise ValueError("Changed output count")
    for old, new, clean in zip(original["outputs"], updated["outputs"], stripped["outputs"]):
        if old["output_id"] != new["output_id"] or old["q"] != new["q"]:
            raise ValueError("Changed class identity or order")
        verify_factorisation(new)
        certificate = new["factorisation_certificate"]
        basis = certificate["basis_in_previous_gate_coordinates"]
        correction = clifford_correction(old["q"], old["representative_gate"],
                                         new["representative_gate"], basis)
        if correction != certificate["diagonal_clifford_correction"]:
            raise ValueError("Incorrect phase correction")
        composed = [combine(v, old["gate_basis_in_tensor_coordinates"]) for v in basis]
        if composed != new["gate_basis_in_tensor_coordinates"]:
            raise ValueError("Incorrect matrix-to-gate basis composition")
        for name in PRESENTATION_FIELDS:
            clean.pop(name)
        clean["representative_gate"] = old["representative_gate"]
        clean["gate_basis_in_tensor_coordinates"] = old["gate_basis_in_tensor_coordinates"]
    registered = {row["output_id"]: row for row in updated["outputs"]}
    for row in updated["outputs"]:
        for factor in row["independent_factors"]:
            authority = registered[factor["output_id"]]
            if (len(authority["independent_factors"]) != 1 or
                    any(authority[k] != factor[k] for k in ("q", "representative_gate"))):
                raise ValueError("Inconsistent registered factor representative")
    if stripped != original:
        raise ValueError("Scientific catalogue fields changed")
    return dict(status="pass", outputs=len(updated["outputs"]),
                protocols=len(updated["protocols"]), scientific_fields_unchanged=True,
                exact_indecomposable_factorisation=True, full_phase_equivalence=True,
                decomposition_counts=dict(sorted(Counter(len(r["independent_factors"])
                                                          for r in updated["outputs"]).items())))
