"""Dependency-free exact affine-equivalence tools for binary point supports."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class AffineMap:
    """Affine map represented by source-basis and target-image data."""

    m: int
    source_origin: int
    target_origin: int
    source_basis: tuple[int, ...]
    target_basis_images: tuple[int, ...]


def _rank(vectors: tuple[int, ...]) -> int:
    basis: dict[int, int] = {}
    for vector in vectors:
        value = int(vector)
        while value:
            pivot = value.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = value
                break
            value ^= basis[pivot]
    return len(basis)


def ordered_basis(vectors: tuple[int, ...], m: int) -> tuple[int, ...]:
    basis: list[int] = []
    for vector in sorted(set(vectors)):
        if vector and _rank((*basis, vector)) > len(basis):
            basis.append(vector)
            if len(basis) == m:
                break
    return tuple(basis)


def affine_basis_from_support(
    support: tuple[int, ...], m: int, origin: int
) -> tuple[int, ...]:
    basis = ordered_basis(tuple(point ^ origin for point in support), m)
    if len(basis) != m:
        raise ValueError("Support does not have full affine rank")
    return basis


def coordinate_map(basis: tuple[int, ...]) -> dict[int, int]:
    coordinates = {0: 0}
    for index, vector in enumerate(basis):
        coordinates.update(
            {
                value ^ vector: mask | (1 << index)
                for value, mask in list(coordinates.items())
            }
        )
    if len(coordinates) != 1 << len(basis):
        raise ValueError("Basis vectors are linearly dependent")
    return coordinates


def support_coordinate_masks(
    support: tuple[int, ...], origin: int, basis: tuple[int, ...]
) -> tuple[int, ...]:
    coordinates = coordinate_map(basis)
    return tuple(sorted(coordinates[point ^ origin] for point in support))


def apply_affine_map(point: int, affine_map: AffineMap) -> int:
    coordinates = coordinate_map(affine_map.source_basis)[
        point ^ affine_map.source_origin
    ]
    image = affine_map.target_origin
    for index, target_vector in enumerate(affine_map.target_basis_images):
        if (coordinates >> index) & 1:
            image ^= target_vector
    return image


def transform_support(
    support: tuple[int, ...], affine_map: AffineMap
) -> tuple[int, ...]:
    return tuple(sorted(apply_affine_map(point, affine_map) for point in support))


def affine_equivalence_map(
    source: tuple[int, ...], target: tuple[int, ...], m: int
) -> AffineMap | None:
    source = tuple(sorted(source))
    target = tuple(sorted(target))
    if (
        len(source) != len(target)
        or len(set(source)) != len(source)
        or len(set(target)) != len(target)
    ):
        return None
    target_set = set(target)
    source_origin = source[0]
    source_basis = affine_basis_from_support(source, m, source_origin)
    source_masks = support_coordinate_masks(source, source_origin, source_basis)
    for target_origin in target:
        target_vectors = tuple(
            point ^ target_origin for point in target if point != target_origin
        )
        if len(ordered_basis(target_vectors, m)) != m:
            continue
        result = _search_target_basis_images(
            m=m,
            source_masks=source_masks,
            target_set=target_set,
            target_origin=target_origin,
            candidates=tuple(sorted(set(target_vectors))),
            chosen=(),
        )
        if result is not None:
            affine_map = AffineMap(
                m=m,
                source_origin=source_origin,
                target_origin=target_origin,
                source_basis=source_basis,
                target_basis_images=result,
            )
            if transform_support(source, affine_map) == target:
                return affine_map
    return None


def affinely_equivalent(
    source: tuple[int, ...], target: tuple[int, ...], m: int
) -> bool:
    return affine_equivalence_map(source, target, m) is not None


def _search_target_basis_images(
    *,
    m: int,
    source_masks: tuple[int, ...],
    target_set: set[int],
    target_origin: int,
    candidates: tuple[int, ...],
    chosen: tuple[int, ...],
) -> tuple[int, ...] | None:
    depth = len(chosen)
    if depth == m:
        return (
            chosen
            if _images_of_source_masks(source_masks, target_origin, chosen)
            <= target_set
            else None
        )
    for candidate in candidates:
        if candidate in chosen or _rank((*chosen, candidate)) <= depth:
            continue
        next_chosen = (*chosen, candidate)
        if not _prefix_constraints_hold(
            source_masks, target_origin, next_chosen, target_set
        ):
            continue
        result = _search_target_basis_images(
            m=m,
            source_masks=source_masks,
            target_set=target_set,
            target_origin=target_origin,
            candidates=candidates,
            chosen=next_chosen,
        )
        if result is not None:
            return result
    return None


def _prefix_constraints_hold(
    source_masks: tuple[int, ...],
    target_origin: int,
    chosen: tuple[int, ...],
    target_set: set[int],
) -> bool:
    prefix_mask = (1 << len(chosen)) - 1
    for source_mask in source_masks:
        if source_mask & ~prefix_mask:
            continue
        image = target_origin
        for index, target_vector in enumerate(chosen):
            if (source_mask >> index) & 1:
                image ^= target_vector
        if image not in target_set:
            return False
    return True


def _images_of_source_masks(
    source_masks: tuple[int, ...],
    target_origin: int,
    target_basis_images: tuple[int, ...],
) -> set[int]:
    images: set[int] = set()
    for source_mask in source_masks:
        image = target_origin
        for index, target_vector in enumerate(target_basis_images):
            if (source_mask >> index) & 1:
                image ^= target_vector
        images.add(image)
    return images

