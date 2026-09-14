"""Certify exact-distance-three Pareto pruning for all outputs on at most four qubits."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from protocol_checks import LogicalTensor, find_equivalence_basis
from space_codec import binary_rank
from verify_protocols import verify as verify_witnesses
from verify_tensor_authorities import restriction, standard_tensor


def verify(root):
    root = Path(root)
    census_path = root / "certificates/tensors/q5_orbits.json"
    frontier_path = root / "data/protocols/pareto_frontier.json"
    census_bytes, frontier_bytes = census_path.read_bytes(), frontier_path.read_bytes()
    census, frontier = json.loads(census_bytes), json.loads(frontier_bytes)
    if (census["status"] != "pass" or census["tensor_count"] != 1 << 25
            or len(census["orbits"]) != 88
            or sum(row["orbit_size"] for row in census["orbits"]) != 1 << 25):
        raise ValueError("The complete five-dimensional tensor census is missing")
    verify_witnesses(frontier_path)
    tensors = {}
    for row in frontier["protocols"]:
        if row["q"] <= 4:
            tensors[row["output_id"]] = LogicalTensor(tuple(int(r, 2) for r in row["generator_matrix_rows"][:row["q"]]))
    records = []
    for row in census["orbits"]:
        q = 5-row["radical_dimension"]
        if not 1 <= q <= 4:
            continue
        tensor = standard_tensor(5, int(row["canonical_standard_signature"], 16))
        radical = [v for v, (_, rank) in enumerate(tensor.contraction_invariants) if rank == 0]
        if binary_rank(radical) != 5-q:
            raise ValueError("Wrong radical dimension in the tensor census")
        span, complement = list(radical), []
        rank = binary_rank(span)
        for bit in range(5):
            value = 1 << bit
            if binary_rank([*span, value]) > rank:
                span.append(value)
                complement.append(value)
                rank += 1
        output = restriction(tensor, complement)
        if output.q != q or not output.intrinsic():
            raise ValueError("The radical complement does not induce the intrinsic tensor")
        matches = []
        for identifier, candidate in tensors.items():
            basis = find_equivalence_basis(output, candidate)
            if basis is not None:
                matches.append((identifier, basis))
        if len(matches) != 1:
            raise ValueError("An intrinsic small output is missing or repeated in the frontier")
        identifier, basis = matches[0]
        incumbents = [r for r in frontier["protocols"] if r["output_id"] == identifier
                      and r["d_Z"] == 3 and r["n"] < 49 and r["S"] <= q+6]
        if not incumbents:
            raise ValueError("No witness dominates every later exact-distance-three candidate")
        incumbent = min(incumbents, key=lambda r: (r["n"], r["S"], r["index"]))
        records.append(dict(q=q, output_id=identifier,
                            five_dimensional_standard_signature=row["canonical_standard_signature"],
                            radical_complement_basis=complement,
                            catalogue_basis_in_complement_coordinates=basis,
                            incumbent_protocol_index=incumbent["index"],
                            incumbent_n=incumbent["n"], incumbent_S=incumbent["S"],
                            candidate_minimum_S=q+6))
    if len(records) != 22 or len({r["output_id"] for r in records}) != 22:
        raise ValueError("The complete set of small intrinsic output orbits is not covered")
    return dict(schema="small-output-exact-distance-three-pruning-v1", status="pass",
                protocol_length_interval=[49, 54], exact_distance=3,
                logical_dimensions=[1, 2, 3, 4], intrinsic_outputs=22,
                outputs_by_dimension=dict(sorted(Counter(r["q"] for r in records).items())),
                stabiliser_dimension_lower_bound=6,
                lower_bound_reason="Distinct stabiliser syndromes imply n<=2^h; n>=49 gives h>=6",
                all_exact_distance_three_candidates_strictly_dominated=True,
                larger_exact_distances_excluded=False, witnesses=records,
                q5_census_sha256=hashlib.sha256(census_bytes).hexdigest(),
                frontier_sha256=hashlib.sha256(frontier_bytes).hexdigest(),
                is_global_completeness_certificate=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = verify(args.root)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="ascii")
    print(json.dumps({k: v for k, v in report.items() if k != "witnesses"}, indent=2))
