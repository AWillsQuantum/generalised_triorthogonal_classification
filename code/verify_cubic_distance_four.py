"""Certify the higher-distance cubic exceptions by alternating output tensors."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from logical_spaces import logical_label_space
from protocol_checks import LogicalTensor, check_matrix, distance_up_to, find_equivalence_basis
from selected_spaces import read_selected
from space_codec import binary_rank
from space_lifts import evaluation_rows
from verify_length52_small_blocks import check_branch
from verify_small_source_closure import all_pointings
from verify_tensor_authorities import standard_tensor
from verify_protocol_case_census import normalized_keys
from verify_protocols import verify as verify_frontier


def verify(root, block_path, replay_directory):
    raw = block_path.read_bytes()
    data = json.loads(raw)
    if (data["schema"] != "finite-length52-small-quotient-blocks-v1"
            or data["matrix_scope"] != "full_projective" or data["protocol_length_interval"] != [51, 52]):
        raise ValueError("Wrong higher-distance cubic source domain")
    domain_path = root / data["source_domain"]
    domain_raw = domain_path.read_bytes()
    if hashlib.sha256(domain_raw).hexdigest() != data["source_domain_sha256"]:
        raise ValueError("Changed finite source domain")
    blocks = [block for block in data["blocks"] if block.get("constant_stabiliser_d4_closure")]
    if len(blocks) != 11:
        raise ValueError("Incomplete constant-stabiliser cubic source set")
    sources = []
    for block in blocks:
        first, end = block["source_interval"]
        if block["c"] != 52 or block["m"] != 7 or end != first+1:
            raise ValueError("The alternating cubic exception has a different source scope")
        sources.append((52, 7, first))
    points, bindings = read_selected(root, sources)
    # There are only 16 alternating trilinear forms in dimension four.
    # In dimensions one and two the only alternating trilinear form is zero.
    alternating = []
    for q in (1, 2, 4):
        triples = q*(q-1)*(q-2)//6
        offset = q+q*(q-1)//2
        tensors = [standard_tensor(q, value << offset) for value in range(1 << triples)]
        if any(tensor.intrinsic() for tensor in tensors):
            raise ValueError("An alternating output tensor is intrinsic in an excluded dimension")
        alternating.append(dict(q=q, tensors=len(tensors), intrinsic_tensors=0))
    frontier_path = root / "data/protocols/pareto_frontier.json"
    verify_frontier(frontier_path)
    frontier = json.loads(frontier_path.read_bytes())
    outputs = {(row["q"], tuple(int(word, 0) for word in row["tensor_key_words"])): row["output_id"]
               for row in frontier["outputs"]}
    tensors = {row["output_id"]: LogicalTensor(tuple(int(s, 2) for s in row["generator_matrix_rows"][:row["q"]]))
               for row in frontier["protocols"]}
    results, replay_bindings = [], {}
    for block, key in zip(blocks, sources, strict=True):
        cases = [row for row in all_pointings(points[key], 7, 52) if 51 <= len(row["points"]) <= 52]
        distribution, positive = Counter(), []
        for case in cases:
            support, h = case["points"], case["ambient_dimension"]
            quotient = logical_label_space(support, h, 4)["quotient_dimension"]
            distribution[quotient] += 1
            if quotient:
                rows = evaluation_rows(support, h, 1)[1:]
                if binary_rank([*rows, (1 << len(support))-1]) != h:
                    raise ValueError("A positive higher-distance quotient lacks the constant stabiliser")
                positive.append(dict(parity_case=case["parity_case"], origin=case["origin"],
                                     ambient_dimension=h, quotient_dimension=quotient))
        if len(cases) != block["pointing_supports"] or dict(distribution) != dict(block["distance_four_quotient_profile"]):
            raise ValueError("Independent distance-four quotient census differs")
        path = replay_directory / f"block{block['index']}_q3.json"
        replay_raw = path.read_bytes()
        replay = json.loads(replay_raw)
        if (replay["schema"] != "finite-length52-branch-replay-v1" or replay["status"] != "pass"
                or replay["source_domain_sha256"] != data["source_domain_sha256"]
                or any(replay[k] != block[k] for k in ("index", "c", "m", "source_interval", "protocol_length_interval"))
                or replay["pointing_supports"] != len(cases)
                or dict(replay["distance_four_quotient_profile"]) != dict(distribution)
                or replay["all_censuses_freshly_recomputed"] is not True):
            raise ValueError("Q3 census does not cover the verified cubic source domain")
        branch = replay["branch"]
        check_branch(branch, dict(distribution), 1)
        if (branch["q"] != 3 or branch["minimum_distance"] != 4 or branch["quotient_dimension_interval"] != [1, 64]
                or branch["statistics"]["marked_orbit_enumerated_supports"]):
            raise ValueError("Incomplete direct q3 census")
        retained, = [row for row in block["branches"] if row["q"] == 3]
        if retained != branch:
            raise ValueError("The finite block does not use the direct q3 census")
        witness_checks = []
        for witness in branch["witnesses"]:
            masks, columns = check_matrix(witness["generator_matrix_rows"], 3)
            tensor = LogicalTensor(masks[:3])
            key_words, = normalized_keys([witness["output_key_words"]])
            identifier = outputs.get((3, key_words))
            basis = find_equivalence_basis(tensor, tensors[identifier]) if identifier is not None else None
            distance, coefficient = distance_up_to(columns, 3)
            if (basis is None or not tensor.intrinsic() or distance is None or distance < 4
                    or (witness["n"], witness["S"], witness["d_Z"], witness["error_coefficient"])
                    != (len(columns), len(masks), distance, coefficient)):
                raise ValueError("Incorrect higher-distance cubic witness")
            dominators = [row["index"] for row in frontier["protocols"] if row["output_id"] == identifier
                          and row["d_Z"] == distance and row["n"] <= len(columns) and row["S"] <= len(masks)]
            if not dominators:
                raise ValueError("A cubic higher-distance witness is not covered by the frontier")
            witness_checks.append(dict(dominator=min(dominators), output_basis=list(basis)))
        replay_bindings[path.relative_to(root).as_posix()] = hashlib.sha256(replay_raw).hexdigest()
        results.append(dict(block=block["index"], source=dict(c=key[0], m=key[1], index=key[2]),
                            pointings=len(cases), distance_four_quotient_profile=sorted(distribution.items()),
                            positive_constant_stabiliser_cases=positive,
                            raw_q3_isotropic_subspaces=branch["statistics"]["isotropic_subspaces"],
                            raw_q3_nondegenerate_subspaces=branch["statistics"]["nondegenerate_subspaces"],
                            witness_checks=witness_checks, new_pareto_points=0))
    return dict(schema="alternating-cubic-higher-distance-verification-v1", status="pass",
        source_domain_sha256=data["source_domain_sha256"], finite_block_data_sha256=hashlib.sha256(raw).hexdigest(),
        all_origin_quotients_independently_recomputed=True, all_positive_quotients_have_constant_stabiliser=True,
        q3_complete_raw_censuses_pareto_dominated=True, q1_q2_q4_excluded_by_alternating_tensor_classification=True,
        all_constant_stabiliser_higher_distance_exceptions_closed=True, alternating_tensors=alternating,
        sources=results, source_shard_bindings=bindings, q3_census_bindings=replay_bindings,
        is_global_completeness_certificate=False)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=root / "data/protocol_sectors/length52_small_quotient_blocks.json")
    parser.add_argument("--censuses", type=Path, default=root / "certificates/length52_cubic_blocks")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify(root, args.input, args.censuses)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps(dict(status="pass", sources=len(result["sources"]),
                         pointings=sum(row["pointings"] for row in result["sources"])), indent=2))


if __name__ == "__main__":
    main()
