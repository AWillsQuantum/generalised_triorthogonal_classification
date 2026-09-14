"""Compare codeword graphs with an exhaustive affine census and exact codeword spans."""

import argparse
from collections import defaultdict
import hashlib
from itertools import permutations
import json
from pathlib import Path
import random

from affine_codeword_graph import canonicalise_many
from affine_graph import affine_images
from selected_spaces import read_selected
from space_codec import binary_rank
from verify_protocol_cover import require


def selected_words(points, m):
    rows = [(1 << len(points))-1]
    rows += [sum(((x >> i) & 1) << j for j, x in enumerate(points)) for i in range(m)]
    words, value = [], 0
    for i in range(1, 1 << len(rows)):
        value ^= rows[(i & -i).bit_length()-1]
        words.append(value)
    words.sort(key=lambda word: (word.bit_count(), word))
    pivots = {}
    for word in words:
        value = word
        while value:
            bit = value.bit_length()-1
            if bit not in pivots:
                pivots[bit] = value
                break
            value ^= pivots[bit]
        if len(pivots) == len(rows):
            threshold = word.bit_count()
            return threshold, sum(w.bit_count() <= threshold for w in words)
    raise ValueError("The codeword family is not full rank")


def apply(point, frame):
    result = frame[0]
    for i, column in enumerate(frame[1:]):
        if point >> i & 1:
            result ^= column
    return result


def verify(root, native):
    inputs = [(3, tuple(x for x in range(8) if mask >> x & 1)) for mask in range(1, 256)]
    inputs = [(m, points) for m, points in inputs if binary_rank(x ^ points[0] for x in points) == m]
    require(len(inputs) == 149, "The small affine census has the wrong input count")
    small_count = len(inputs)
    keys = {(54, 13, i) for i in (0, 376, 377, 383, 385, 937)}
    keys |= {(54, 14, i) for i in (0, 23, 46)} | {(54, 15, i) for i in (0, 1)}
    supports, bindings = read_selected(root, keys)
    rng = random.Random(5413)
    high_pairs = []
    for (_, m, _), points in sorted(supports.items()):
        basis = [1 << i for i in range(m)]
        for _ in range(10*m):
            a, b = rng.sample(range(m), 2)
            basis[a] ^= basis[b]
        frame = [rng.randrange(1 << m), *basis]
        transformed = tuple(sorted(apply(x, frame) for x in points))
        high_pairs.append(len(inputs))
        inputs.extend(((m, points), (m, transformed)))
    forms = canonicalise_many([native], inputs)
    for (m, points), form in zip(inputs, forms, strict=True):
        require(selected_words(points, m) == (form["maximum_codeword_weight"], form["selected_codewords"]),
                "The native graph does not contain exactly the selected spanning codewords")
    maps = [affine_images((origin, *basis), 3) for basis in permutations(range(1, 8), 3)
            if binary_rank(basis) == 3 for origin in range(8)]
    graph_to_exact, exact_to_graph = defaultdict(set), defaultdict(set)
    for (_, points), form in zip(inputs[:small_count], forms[:small_count], strict=True):
        exact = min(tuple(sorted(mapping[x] for x in points)) for mapping in maps)
        graph = tuple(form["canonical_points"])
        graph_to_exact[graph].add(exact)
        exact_to_graph[exact].add(graph)
    require(len(exact_to_graph) == len(graph_to_exact) == 5
            and all(len(v) == 1 for v in (*exact_to_graph.values(), *graph_to_exact.values())),
            "The graph and exhaustive affine partitions disagree")
    for i in high_pairs:
        require(forms[i]["canonical_points"] == forms[i+1]["canonical_points"], "An affine pair has different canonical keys")
    for name in ("code/native/auxiliary/affine_codeword_graph.cpp", "code/native/vendor/bliss-0.77.zip"):
        bindings[name] = hashlib.sha256((root / name).read_bytes()).hexdigest()
    return dict(schema="affine-codeword-graph-verification-v1", status="pass", exhaustive_small_supports=small_count,
        exact_affine_maps=len(maps), exhaustive_small_classes=5, larger_affine_pairs=len(high_pairs),
        every_affine_frame_replayed=True, every_selected_codeword_span_independently_checked=True,
        every_small_equivalence_and_inequivalence_checked=True, every_larger_affine_pair_agrees=True,
        dependencies=bindings)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.root.resolve(), args.native.resolve())
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps({k: v for k, v in result.items() if k != "dependencies"}, indent=2))


if __name__ == "__main__":
    main()
