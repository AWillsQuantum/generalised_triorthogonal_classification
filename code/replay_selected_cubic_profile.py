"""Reproduce a selected cubic q5 chain census or complete q6 marked census."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

from protocol_checks import LogicalTensor, check_matrix, distance_up_to, find_equivalence_basis
from protocol_domain import open_domain, read_header, read_pointings, reconstruct, source_key
from selected_spaces import read_selected
from verify_protocol_case_census import normalized_keys, write_pointed_cases


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--q", type=int, choices=(5, 6), required=True)
    parser.add_argument("--pointing", type=int, required=True)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    path = root / "data/protocol_sectors/length52_selected_cubic_censuses.json"
    raw = path.read_bytes()
    data = json.loads(raw)
    if data["schema"] != "finite-selected-cubic-protocol-censuses-v1" or args.workers < 1:
        raise ValueError("Invalid finite census domain")
    profiles = {row["pointing_index"]: row for row in data[f"q{args.q}_profiles"]}
    if args.pointing not in profiles:
        raise ValueError("No recorded finite branch at this pointed input")
    expected = profiles[args.pointing]
    domain_raw = (root / data["source_domain"]).read_bytes()
    if hashlib.sha256(domain_raw).hexdigest() != data["source_domain_sha256"]:
        raise ValueError("Changed source domain")
    domain = json.loads(domain_raw)
    pointing_path = root / domain["pointing_data"]
    with pointing_path.open("rb") as stream:
        if hashlib.file_digest(stream, "sha256").hexdigest() != domain["pointing_sha256"]:
            raise ValueError("Changed finite pointing data")
    target = None
    with open_domain(pointing_path) as stream:
        header = read_header(stream)
        if header != dict(source_count=len(domain["sources"]), count=domain["pointing_count"], maximum_protocol_length=52):
            raise ValueError("Pointing header disagrees with the finite source domain")
        for i, record in enumerate(read_pointings(stream, header)):
            if i == args.pointing:
                target = record
    if target is None:
        raise ValueError("Pointed support is outside the finite domain")
    source = domain["sources"][target[0]]
    supports, bindings = read_selected(root, {source_key(source)})
    points, h = reconstruct(domain["sources"], supports, target)
    work = args.work_directory.resolve()
    work.mkdir(parents=True, exist_ok=True)
    manifest, output = work / "case.utsp", work / "census.json"
    with manifest.open("wb") as stream:
        write_pointed_cases(stream, dict(maximum_protocol_length=52,
            cases=[dict(index=0, points=points, ambient_dimension=h)]))
    mode = "--q5-q3-chain-cover-orbits" if args.q == 5 else "--marked-orbits"
    result = subprocess.run([str(args.native.resolve()), "manifest-catalogue", "--input", str(manifest),
        "--output", str(output), "--q", str(args.q), "--minimum-distance", "3", mode,
        "--workers", str(args.workers), "--emit-positive-supports"], capture_output=True, text=True)
    if result.returncode:
        raise ValueError("Native finite census failed: "+result.stderr.strip())
    value = json.loads(output.read_bytes())
    expected_mode = "q5_q3_chain_cover_marked_code_orbits" if args.q == 5 else "complete_marked_code_orbits"
    if (value["status"] != "complete" or value["matrix_scope"] != "full_projective"
            or value["logical_qubits"] != args.q or value["minimum_distance"] != 3
            or value["enumeration_mode"] != expected_mode
            or value["support_range"] != dict(start=0, count=1, end_exclusive=1)):
        raise ValueError("Incomplete finite case enumeration")
    keys = normalized_keys([row["output"]["canonical_key_words"] for row in value["pareto_protocols"]])
    if (keys != normalized_keys(expected["output_keys"])
            or any(value["statistics"][name] != expected["statistics"][name]
                   for name in ("nondegenerate_subspaces", "marked_code_orbits", "canonical_output_orbits"))):
        raise ValueError("Fresh complete output keys or marked subspace masses differ")
    frontier_path = root / "data/protocols/pareto_frontier.json"
    frontier_raw = frontier_path.read_bytes()
    frontier = json.loads(frontier_raw)
    outputs = {(row["q"], tuple(int(word, 0) for word in row["tensor_key_words"])): row["output_id"]
               for row in frontier["outputs"]}
    witnesses = []
    for row in value["pareto_protocols"]:
        q = args.q
        masks, columns = check_matrix(row["generator_rows"], q)
        distance, coefficient = distance_up_to(columns, q)
        tensor = LogicalTensor(masks[:q])
        key, = normalized_keys([row["output"]["canonical_key_words"]])
        identifier = outputs.get((q, key))
        if (not tensor.intrinsic() or identifier is None or distance != 3
                or (row["protocol_length_n"], row["space_footprint_S"], row["d_Z"], row["error_coefficient"])
                != (len(columns), len(masks), distance, coefficient)):
            raise ValueError("Incorrect finite replay witness")
        choices = [r for r in frontier["protocols"] if r["output_id"] == identifier and r["d_Z"] == distance
                   and r["n"] <= len(columns) and r["S"] <= len(masks)]
        if not choices:
            raise ValueError("Regenerated finite witness is not Pareto covered")
        chosen = min(choices, key=lambda r: (r["n"], r["S"], r["index"]))
        basis = find_equivalence_basis(tensor, LogicalTensor(tuple(int(s, 2) for s in chosen["generator_matrix_rows"][:q])))
        if basis is None:
            raise ValueError("Regenerated output key does not describe the witness")
        witnesses.append(dict(q=q, n=len(columns), S=len(masks), d_Z=distance, error_coefficient=coefficient,
                              output_id=identifier, generator_matrix_rows=row["generator_rows"],
                              dominator=chosen["index"], output_basis=list(basis)))
    certificate = dict(schema="selected-cubic-profile-replay-v1", status="pass", q=args.q,
        pointing_index=args.pointing, source_dataset_sha256=hashlib.sha256(raw).hexdigest(),
        source_domain_sha256=data["source_domain_sha256"], source_shard_bindings=bindings,
        complete_finite_case_recomputed=True, primitive_q5_target_census_separate=args.q == 5,
        output_keys=[list(key) for key in sorted(keys)],
        marked_code_orbits=value["statistics"]["marked_code_orbits"],
        nondegenerate_subspaces=value["statistics"]["nondegenerate_subspaces"],
        frontier_sha256=hashlib.sha256(frontier_raw).hexdigest(), witnesses=witnesses,
        is_global_completeness_certificate=False)
    args.output.write_text(json.dumps(certificate, indent=2)+"\n", encoding="ascii")
    print(json.dumps({k: v for k, v in certificate.items() if k not in ("witnesses", "source_shard_bindings")}, indent=2))


if __name__ == "__main__":
    main()
