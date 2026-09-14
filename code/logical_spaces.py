"""Independent binary equations for the distance-filtered logical label space."""

from itertools import combinations

from space_codec import binary_rank, checked_support
from space_lifts import evaluation_rows, linear_solutions, reduce_vector


def zero_syndrome_constraints(points, maximum_weight):
    """Yield characteristic vectors of all small nonempty zero-sum subsets."""
    points = tuple(points)
    if maximum_weight < 0 or len(points) != len(set(points)):
        raise ValueError("Expected a set of syndromes and a nonnegative weight")
    positions = {point: index for index, point in enumerate(points)}
    if maximum_weight >= 1 and 0 in positions:
        yield 1 << positions[0]
    # A set of distinct binary syndromes has no weight-two zero-sum subset.
    if maximum_weight >= 3:
        for first, second in combinations(range(len(points)), 2):
            third = positions.get(points[first] ^ points[second])
            if third is not None and second < third:
                yield (1 << first) | (1 << second) | (1 << third)
    for weight in range(4, min(maximum_weight, len(points))+1):
        for indices in combinations(range(len(points)), weight):
            syndrome = 0
            mask = 0
            for index in indices:
                syndrome ^= points[index]
                mask |= 1 << index
            if not syndrome:
                yield mask


def logical_label_space(points, h, minimum_distance=3):
    """Return V_d/L and its common-isotropy forms, without enumerating subspaces."""
    points = checked_support(points, h)
    if not points or binary_rank(points) != h or not 3 <= minimum_distance <= len(points)+1:
        raise ValueError("Invalid spanning stabiliser support or distance threshold")
    if any(row.bit_count() & 1 for row in evaluation_rows(points, h, 3)[1:]):
        raise ValueError("Stabiliser support fails triorthogonality")
    equations = set(evaluation_rows(points, h, 2)[1:])
    equations.update(zero_syndrome_constraints(points, minimum_distance-1))
    stabilisers = evaluation_rows(points, h, 1)[1:]
    result = linear_solutions(sorted(equations), [0]*len(equations), len(points))
    if result is None:
        raise AssertionError("Homogeneous equations cannot be inconsistent")
    _, kernel = result
    pivots = {}
    for row in stabilisers:
        if any((row & equation).bit_count() & 1 for equation in equations):
            raise AssertionError("Stabiliser row violates a logical-label equation")
        reduced = reduce_vector(row, pivots)
        if not reduced:
            raise AssertionError("Stabiliser rows are dependent")
        pivots[reduced.bit_length()-1] = reduced
    basis = []
    for row in kernel:
        reduced = reduce_vector(row, pivots)
        if reduced:
            pivots[reduced.bit_length()-1] = reduced
            basis.append(reduced)
    forms = []
    for row in stabilisers:
        form = tuple(sum(((left & right & row).bit_count() & 1) << j
                         for j, right in enumerate(basis)) for left in basis)
        if any(form[i] >> i & 1 for i in range(len(basis))):
            raise AssertionError("A logical common-isotropy form is not alternating")
        forms.append(form)
    return dict(stabiliser_dimension=h, minimum_distance=minimum_distance,
                constraint_rank=len(points)-len(kernel),
                quotient_dimension=len(basis), basis=tuple(basis),
                common_isotropy_forms=tuple(forms))
