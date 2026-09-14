"""Compare the graph backend with all AGL(3,2) maps on all 256 truth tables."""

import argparse
from collections import defaultdict
from itertools import permutations
import json
from pathlib import Path

from affine_graph import affine_images, canonicalise_many, rank


def verify(command):
    maps = [affine_images((origin, *basis), 3)
            for basis in permutations(range(1, 8), 3) if rank(basis) == 3
            for origin in range(8)]
    if len(maps) != 1344:
        raise AssertionError("Incorrect affine-map count")
    results = canonicalise_many(command, [(3, mask) for mask in range(256)])
    graph_by_reference = defaultdict(set)
    reference_by_graph = defaultdict(set)
    for mask, result in enumerate(results):
        orbit = [sum(1 << y for x, y in enumerate(images) if mask >> x & 1)
                 for images in maps]
        reference = min(orbit)
        graph = int(result["canonical_mask"], 16)
        graph_by_reference[reference].add(graph)
        reference_by_graph[graph].add(reference)
        if result["stabilizer_order"] != sum(image == mask for image in orbit):
            raise AssertionError("Incorrect graph stabiliser order")
    if len(graph_by_reference) != 10 or len(reference_by_graph) != 10:
        raise AssertionError("Incorrect graph orbit census")
    if not all(len(values) == 1 for values in (*graph_by_reference.values(), *reference_by_graph.values())):
        raise AssertionError("Graph/reference equivalence disagreement")
    return dict(status="pass", functions=256, affine_maps=1344, affine_classes=10,
                every_stabilizer_order_and_equivalence_class_checked=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = verify([args.native])
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
