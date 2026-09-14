"""Exact witnessed affine transporters for supports in binary affine space."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from triorthogonal_utils import gf2_rank


CompactWitness = tuple[int, ...]


@dataclass
class SearchMetrics:
    target_origin_count: int = 0
    search_node_count: int = 0
    rejected_candidate_image_count: int = 0


@dataclass(frozen=True)
class PreparedSupport:
    dimension: int
    points: tuple[int, ...]
    difference_counts: tuple[int, ...]
    difference_multiset: tuple[int, ...]
    profiles: dict[int, tuple[int, ...]]
    profile_multiset: tuple[tuple[int, ...], ...]
    source_origin: int
    translated_mask: int
    source_basis: tuple[int, ...]
    source_spans: tuple[tuple[int, ...], ...]


def difference_counts(points: tuple[int, ...], dimension: int) -> tuple[int, ...]:
    counts = [0] * (1 << dimension)
    for index, left in enumerate(points):
        for right in points[index + 1 :]:
            counts[left ^ right] += 2
    return tuple(counts)


def local_profiles(
    points: tuple[int, ...], counts: tuple[int, ...]
) -> dict[int, tuple[int, ...]]:
    return {
        point: tuple(
            sorted(counts[point ^ other] for other in points if other != point)
        )
        for point in points
    }


def translated_mask(points: tuple[int, ...], origin: int) -> int:
    mask = 0
    for point in points:
        mask |= 1 << (point ^ origin)
    return mask


def prefix_spans(basis: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    spans: list[tuple[int, ...]] = [(0,)]
    for vector in basis:
        prior = spans[-1]
        spans.append(prior + tuple(vector ^ value for value in prior))
    return tuple(spans)


def choose_source_basis(
    translated_points: tuple[int, ...],
    source_origin: int,
    profiles: dict[int, tuple[int, ...]],
    counts: tuple[int, ...],
    dimension: int,
) -> tuple[int, ...]:
    features = {
        vector: (counts[vector], profiles[source_origin ^ vector])
        for vector in translated_points
        if vector
    }
    frequencies = Counter(features.values())
    ordered = sorted(
        features,
        key=lambda vector: (frequencies[features[vector]], features[vector], vector),
    )
    basis: list[int] = []
    for vector in ordered:
        if gf2_rank((*basis, vector)) > len(basis):
            basis.append(vector)
            if len(basis) == dimension:
                return tuple(basis)
    raise ValueError("Support does not affinely span the stated dimension")


def prepare_support(points: tuple[int, ...], dimension: int) -> PreparedSupport:
    points = tuple(sorted(int(point) for point in points))
    size = 1 << dimension
    if (
        not points
        or len(points) != len(set(points))
        or any(not 0 <= point < size for point in points)
    ):
        raise ValueError("Support must contain distinct ambient points")
    counts = difference_counts(points, dimension)
    profiles = local_profiles(points, counts)
    origin = points[0]
    vectors = tuple(sorted(point ^ origin for point in points))
    basis = choose_source_basis(vectors, origin, profiles, counts, dimension)
    return PreparedSupport(
        dimension=dimension,
        points=points,
        difference_counts=counts,
        difference_multiset=tuple(sorted(counts[1:])),
        profiles=profiles,
        profile_multiset=tuple(sorted(profiles.values())),
        source_origin=origin,
        translated_mask=translated_mask(points, origin),
        source_basis=basis,
        source_spans=prefix_spans(basis),
    )


def identity_witness(dimension: int) -> CompactWitness:
    return (0, *(1 << bit for bit in range(dimension)))


def apply_witness(point: int, witness: CompactWitness) -> int:
    image = witness[0]
    for bit, basis_image in enumerate(witness[1:]):
        if (point >> bit) & 1:
            image ^= basis_image
    return image


def witness_is_valid(
    source: PreparedSupport,
    target: PreparedSupport,
    witness: CompactWitness,
) -> bool:
    dimension = source.dimension
    return (
        target.dimension == dimension
        and len(witness) == dimension + 1
        and 0 <= witness[0] < (1 << dimension)
        and all(0 <= value < (1 << dimension) for value in witness[1:])
        and gf2_rank(witness[1:]) == dimension
        and {apply_witness(point, witness) for point in source.points}
        == set(target.points)
    )


def transporter(
    source: PreparedSupport,
    target: PreparedSupport,
    *,
    metrics: SearchMetrics | None = None,
) -> CompactWitness | None:
    """Return an exact affine source-to-target witness, or prove none exists."""

    if source.dimension != target.dimension or len(source.points) != len(target.points):
        return None
    if source.difference_multiset != target.difference_multiset:
        return None
    if source.profile_multiset != target.profile_multiset:
        return None
    if metrics is None:
        metrics = SearchMetrics()

    dimension = source.dimension
    source_origin = source.source_origin
    source_profile = source.profiles[source_origin]
    source_basis = source.source_basis
    source_spans = source.source_spans

    for target_origin in target.points:
        if target.profiles[target_origin] != source_profile:
            continue
        metrics.target_origin_count += 1
        target_mask = translated_mask(target.points, target_origin)
        target_vectors = tuple(
            sorted(
                point ^ target_origin
                for point in target.points
                if point != target_origin
            )
        )
        target_span_sets: list[set[int]] = []

        def search(
            depth: int,
            target_span: tuple[int, ...],
            images: tuple[int, ...],
        ) -> tuple[int, ...] | None:
            metrics.search_node_count += 1
            if depth == dimension:
                return images
            source_vector = source_basis[depth]
            source_span = source_spans[depth]
            if len(target_span_sets) <= depth:
                target_span_sets.append(set(target_span))
            else:
                target_span_sets[depth] = set(target_span)
            occupied = target_span_sets[depth]
            source_point_profile = source.profiles[source_origin ^ source_vector]
            for target_vector in target_vectors:
                if (
                    target_vector in occupied
                    or target.difference_counts[target_vector]
                    != source.difference_counts[source_vector]
                    or target.profiles[target_origin ^ target_vector]
                    != source_point_profile
                ):
                    continue
                valid = True
                new_half = []
                for source_old, target_old in zip(
                    source_span, target_span, strict=True
                ):
                    source_new = source_vector ^ source_old
                    target_new = target_vector ^ target_old
                    source_member = (source.translated_mask >> source_new) & 1
                    target_member = (target_mask >> target_new) & 1
                    if (
                        source_member != target_member
                        or source.difference_counts[source_new]
                        != target.difference_counts[target_new]
                        or (
                            source_member
                            and source.profiles[source_origin ^ source_new]
                            != target.profiles[target_origin ^ target_new]
                        )
                    ):
                        valid = False
                        metrics.rejected_candidate_image_count += 1
                        break
                    new_half.append(target_new)
                if not valid:
                    continue
                result = search(
                    depth + 1,
                    target_span + tuple(new_half),
                    images + (target_vector,),
                )
                if result is not None:
                    return result
            return None

        images = search(0, (0,), ())
        if images is None:
            continue
        full_source_span = source_spans[dimension]
        full_target_span = prefix_spans(images)[dimension]
        source_coordinates = {
            value: coefficient for coefficient, value in enumerate(full_source_span)
        }
        standard_basis_images = tuple(
            full_target_span[source_coordinates[1 << bit]]
            for bit in range(dimension)
        )
        linear_source_origin = full_target_span[source_coordinates[source_origin]]
        witness = (
            target_origin ^ linear_source_origin,
            *standard_basis_images,
        )
        if not witness_is_valid(source, target, witness):
            raise AssertionError("Basis-image search produced an invalid witness")
        return witness
    return None
