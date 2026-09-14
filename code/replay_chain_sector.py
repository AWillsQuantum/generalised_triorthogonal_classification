"""Reproduce a chosen interval of the finite q5 chain or complete q6 census."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

from protocol_checks import LogicalTensor, check_matrix, distance_up_to, find_equivalence_basis
from verify_protocol_case_census import normalized_keys, write_pointed_cases


def cases_for_interval(data, domain, q, start, count):
    if q not in (5, 6):
        raise ValueError("This finite census contains dimensions five and six")
    indices = data["q5_input_support_class_indices"] if q == 5 else [
        row["support_class_index"] for row in data["q6_candidates"]]
    if start < 0 or count < 1 or start+count > len(indices):
        raise ValueError("Invalid finite census interval")
    cases = []
    for position, index in enumerate(indices[start:start+count]):
        row = domain["classes"][index]
        h, n = row["ambient_dimension"], row["support_length"]
        words = [int(word, 16) for word in row["key_words"]]
        if len(words) != 1+2*h or words[0] != n | h << 8 or words[1:1+h] != words[1+h:]:
            raise ValueError("Invalid canonical support key")
        points = sorted(sum((words[1+i] >> j & 1) << i for i in range(h)) for j in range(n))
        cases.append(dict(index=position, input_index=start+position, support_class_index=index,
                          ambient_dimension=h, points=points))
    return cases, len(indices)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--q", type=int, choices=(5, 6), required=True)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--count", type=int, default=1)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("Workers must be positive")
    data_path = root / "data/protocol_sectors/nonprimitive_chain_sector.json"
    data_raw = data_path.read_bytes()
    data = json.loads(data_raw)
    domain_raw = (root / data["domain_dataset"]).read_bytes()
    domain = json.loads(domain_raw)
    cases, total = cases_for_interval(data, domain, args.q, args.start, args.count)
    work = args.work_directory.resolve()
    work.mkdir(parents=True, exist_ok=True)
    manifest = work / "inputs.utsp"
    with manifest.open("wb") as stream:
        write_pointed_cases(stream, dict(maximum_protocol_length=54, cases=cases))
    output = work / "census.json"
    mode = "--q5-q3-chain-cover-orbits" if args.q == 5 else "--marked-orbits"
    subprocess.run([str(args.native), "manifest-catalogue", "--input", str(manifest), "--output", str(output),
                    "--q", str(args.q), "--minimum-distance", "3", mode,
                    "--emit-positive-supports", "--workers", str(args.workers)], check=True)
    result = json.loads(output.read_bytes())
    expected_mode = "q5_q3_chain_cover_marked_code_orbits" if args.q == 5 else "complete_marked_code_orbits"
    if (result["status"] != "complete" or result["matrix_scope"] != "full_projective"
            or result["logical_qubits"] != args.q or result["minimum_distance"] != 3
            or result["enumeration_mode"] != expected_mode
            or result["support_range"] != dict(start=0, count=args.count, end_exclusive=args.count)
            or result["statistics"]["supports_processed"] != args.count):
        raise ValueError("Incomplete finite census")
    frontier_path = root / "data/protocols/pareto_frontier.json"
    frontier_raw = frontier_path.read_bytes()
    frontier = json.loads(frontier_raw)
    outputs = {(row["q"], tuple(int(word, 0) for word in row["tensor_key_words"])): row for row in frontier["outputs"]}
    witnesses = []
    for row in result["pareto_protocols"]:
        q = args.q
        masks, columns = check_matrix(row["generator_rows"], q)
        distance, coefficient = distance_up_to(columns, q)
        if (distance is None or distance < 3 or (distance, coefficient) != (row["d_Z"], row["error_coefficient"])
                or (row["protocol_length_n"], row["space_footprint_S"])
                != (len(row["generator_rows"][0]), len(row["generator_rows"]))):
            raise ValueError("Incorrect regenerated finite witness metrics")
        key, = normalized_keys([row["output"]["canonical_key_words"]])
        if (q, key) not in outputs:
            raise ValueError("Regenerated output is absent from the frontier")
        output_id = outputs[q, key]["output_id"]
        choices = [p for p in frontier["protocols"] if p["output_id"] == output_id and p["d_Z"] == distance
                   and p["n"] <= row["protocol_length_n"] and p["S"] <= row["space_footprint_S"]]
        if not choices:
            raise ValueError("Regenerated protocol is not covered by the frontier")
        chosen = min(choices, key=lambda p: (p["n"], p["S"], p["index"]))
        basis = find_equivalence_basis(LogicalTensor(masks[:q]),
            LogicalTensor(tuple(int(s, 2) for s in chosen["generator_matrix_rows"][:q])))
        if basis is None:
            raise ValueError("Regenerated output equivalence is incorrect")
        witnesses.append(dict(dominator=chosen["index"], output_basis=list(basis)))
    report = dict(schema="finite-chain-sector-interval-replay-v1", status="pass", q=args.q,
        start=args.start, count=args.count, total_input_classes=total,
        entire_dimension_replayed=args.start == 0 and args.count == total,
        complete_interval_enumeration=True, witness_checks=witnesses,
        input_support_class_indices=[row["support_class_index"] for row in cases],
        primitive_q5_requires_separate_target_census=args.q == 5,
        source_dataset_sha256=hashlib.sha256(data_raw).hexdigest(),
        domain_dataset_sha256=hashlib.sha256(domain_raw).hexdigest(),
        frontier_sha256=hashlib.sha256(frontier_raw).hexdigest(),
        is_global_completeness_certificate=False)
    args.output.write_text(json.dumps(report, indent=2)+"\n", encoding="ascii")
    print(json.dumps({k: v for k, v in report.items() if k != "witness_checks"}, indent=2))


if __name__ == "__main__":
    main()
