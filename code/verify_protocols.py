"""Verify all distributed frontier witnesses, not enumeration completeness."""

import argparse
from collections import Counter
import hashlib
from itertools import combinations
import json
from pathlib import Path

from protocol_checks import (LogicalTensor, check_matrix, exact_distance_and_coefficient,
                             find_equivalence_basis, gate_tensor, verify_basis)


def verify(path):
    data = json.loads(path.read_bytes())
    if data["schema"] != "triorthogonal-protocol-witnesses-v1":
        raise ValueError("Unsupported witness format")
    if (data["maximum_protocol_length"] != 54 or data["minimum_distance"] != 3
            or data["equivalence"] != "CNOT+S" or data["pareto_objectives"] != ["n", "S"]
            or data["pareto_partition"] != ["output_id", "exact d_Z"]):
        raise ValueError("Incorrect classification scope")
    outputs = {row["output_id"]: row for row in data["outputs"]}
    if len(outputs) != len(data["outputs"]):
        raise ValueError("Duplicate output identifiers")
    tensors = {}
    for index, protocol in enumerate(data["protocols"]):
        if protocol["index"] != index:
            raise ValueError("Noncontiguous protocol indices")
        output = outputs[protocol["output_id"]]
        q = output["q"]
        if protocol["q"] != q:
            raise ValueError("Conflicting logical dimensions")
        rows = protocol["generator_matrix_rows"]
        masks, columns = check_matrix(rows, q)
        tensor = LogicalTensor(masks[:q])
        if not tensor.intrinsic() or tensor.key_words() != output["tensor_key_words"]:
            raise ValueError("Incorrect output tensor")
        key = "Q" + str(q) + "_" + "_".join(word.removeprefix("0x") for word in tensor.key_words())
        if key != output["output_id"]:
            raise ValueError("Incorrect output identifier")
        verify_basis(tensor, gate_tensor(q, output["representative_gate"]),
                     output["gate_basis_in_tensor_coordinates"])
        distance, coefficient = exact_distance_and_coefficient(columns, q)
        if (protocol["n"], protocol["S"], protocol["d_Z"], protocol["error_coefficient"]) != (
                len(rows[0]), len(rows), distance, coefficient):
            raise ValueError("Incorrect protocol parameters")
        if not 1 <= len(rows[0]) <= 54 or distance is None or distance < 3:
            raise ValueError("Protocol outside scope")
        tensors[output["output_id"]] = tensor
    if set(tensors) != set(outputs):
        raise ValueError("Output without a witness")
    pair_checks = 0
    for a, b in combinations(tensors, 2):
        if tensors[a].q != tensors[b].q:
            continue
        if find_equivalence_basis(tensors[a], tensors[b]) is not None:
            raise ValueError("Equivalent outputs have different identifiers")
        pair_checks += 1
    metrics = set()
    for a in data["protocols"]:
        metric = (a["output_id"], a["d_Z"], a["n"], a["S"])
        if metric in metrics:
            raise ValueError("Duplicate Pareto metric point")
        metrics.add(metric)
        for b in data["protocols"]:
            if (a["output_id"], a["d_Z"]) == (b["output_id"], b["d_Z"]):
                if (b["n"] <= a["n"] and b["S"] <= a["S"]
                        and (b["n"] < a["n"] or b["S"] < a["S"])):
                    raise ValueError("Dominated frontier point")
    if len(data["protocols"]) != 74 or len(outputs) != 62:
        raise ValueError("Unexpected frontier size")
    return dict(schema="triorthogonal-witness-verification-v1", status="pass",
                catalogue_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                protocols=74, outputs=62, output_inequivalence_pairs_checked=pair_checks,
                distance_counts=dict(sorted(Counter(p["d_Z"] for p in data["protocols"]).items())),
                checks=dict(full_rank=True, distinct_columns=True, mixed_overlap_conditions=True,
                            intrinsic_outputs=True, gate_basis_certificates=True,
                            exact_distances=True, leading_error_coefficients=True,
                            pairwise_output_inequivalence=True, pareto_nondominance_within_catalogue=True),
                enumeration_completeness_established=False,
                note="This verifies witnesses, output separation and nondominance within the supplied list. Exhaustive coverage requires the separate enumeration evidence.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument("--catalogue", type=Path, default=root / "data/protocols/pareto_frontier.json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = verify(args.catalogue)
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
