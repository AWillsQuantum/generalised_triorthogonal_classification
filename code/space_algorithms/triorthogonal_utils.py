"""Exact binary utilities shared by the length-54 classification scripts."""

from __future__ import annotations

import ast
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from typing import Iterable


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_write_json(path: Path, data: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def gf2_rank(vectors: Iterable[int]) -> int:
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


def affine_rank(points: Iterable[int]) -> int:
    values = tuple(int(point) for point in points)
    if not values:
        return -1
    base = values[0]
    return gf2_rank(point ^ base for point in values[1:])


def has_even_moments(m: int, support: Iterable[int], max_degree: int = 3) -> bool:
    points = tuple(int(point) for point in support)
    for degree in range(max_degree + 1):
        for monomial in combinations(range(m), degree):
            parity = 0
            for point in points:
                value = 1
                for variable in monomial:
                    value &= (point >> (m - 1 - variable)) & 1
                parity ^= value
            if parity:
                return False
    return True


def parse_yaml_scalar(raw_value: str) -> object:
    value = raw_value.strip()
    if value.startswith(("'", '"')):
        return ast.literal_eval(value)
    if value in {"true", "false"}:
        return value == "true"
    try:
        return int(value)
    except ValueError:
        return value


def load_generator_matrix_yaml(path: Path) -> list[dict[str, object]]:
    """Parse the deliberately simple generated YAML without a YAML dependency."""

    entries: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    list_mode: str | None = None
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped.startswith("- number:"):
            if current is not None:
                entries.append(current)
            current = {"number": int(stripped.split(":", 1)[1].strip())}
            list_mode = None
            continue
        if current is None:
            continue
        if stripped == "support_without_unital_row:":
            current["support"] = []
            list_mode = "support"
            continue
        if stripped == "generator_matrix_rows:":
            current["generator_matrix_rows"] = []
            list_mode = "generator_matrix_rows"
            continue
        if list_mode is not None and stripped.startswith("- "):
            values = current[list_mode]
            assert isinstance(values, list)
            values.append(parse_yaml_scalar(stripped[2:].strip()))
            continue
        if ":" not in stripped:
            continue
        list_mode = None
        key, raw_value = stripped.split(":", 1)
        if key in {
            "polynomial_latex",
            "r",
            "c",
            "rank",
            "triorthogonal",
            "no_repeated_columns",
        }:
            current[key] = parse_yaml_scalar(raw_value)
    if current is not None:
        entries.append(current)
    return entries


def support_strings_to_ints(support: Iterable[str]) -> tuple[int, ...]:
    return tuple(int(point, 2) if point else 0 for point in support)
