"""Dimension-independent inverse-contraction recurrence on positive cores."""

from itertools import combinations
from math import comb
from pathlib import Path
import sys

from space_codec import binary_rank, validate_unital_support

sys.path.insert(0, str(Path(__file__).resolve().parent / "space_algorithms"))
from inverse_contraction_extensions import CoreExtensionContext, xor_quadratic_evaluations


def multiplicity_bound(length, dimension):
    if type(length) is not int or length < 1 or type(dimension) is not int or dimension < 1:
        raise ValueError("Invalid support parameters")
    return comb(length, 2)//((1 << dimension)-1)


def iter_space_children(core, core_dimension, target_length, target_dimension):
    """Yield all admissible lifts of one core modulo affine shears.

    Affine equivalence between children is deliberately not imposed here.
    Taking all precursor classes, then an exact affine quotient, is complete
    whenever the averaging bound excludes the empty-core case. This includes
    every target of length at most 54 and affine dimension at least eight.
    """
    core = tuple(core)
    if not 16 <= target_length <= 54 or target_length % 2 or not 1 <= target_dimension <= 17:
        raise ValueError("Target is outside the supported finite classification range")
    if not 4 <= core_dimension < target_dimension:
        raise ValueError("The core must have smaller positive affine dimension")
    validate_unital_support(core, core_dimension)
    bound = multiplicity_bound(target_length, target_dimension)
    if 2*bound >= target_length:
        raise ValueError("This sector requires a separate empty-core or full-cube seed")
    difference = target_length-len(core)
    if difference < 0 or difference % 2 or difference//2 > bound:
        return
    multiplicity = difference//2
    # Every embedding of this intrinsic core into the quotient ambient space
    # is equivalent; append zero coordinates and retain the full complement.
    context = CoreExtensionContext.build(target_dimension-1, core, core_dimension+1)
    if multiplicity <= 5:
        fibres = context.iter_fiber_sets(multiplicity, 0)
    else:
        fibres = (selection for selection in combinations(context.external_points, multiplicity)
                  if context.solver.solve(xor_quadratic_evaluations(selection, target_dimension-1)) is not None)
    for full_fibres in fibres:
        labels = context.solver.solve(xor_quadratic_evaluations(full_fibres, target_dimension-1))
        if labels is None:
            raise AssertionError("The compatible-fibre iterator produced an inconsistent system")
        basis = context.lift_complement
        for index in range(1 << len(basis)):
            if index:
                labels ^= basis[(index & -index).bit_length()-1]
            support = context.reconstructed_support(full_fibres, labels)
            if binary_rank(x ^ support[0] for x in support) == target_dimension:
                yield support
