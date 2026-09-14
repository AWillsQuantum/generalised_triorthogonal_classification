"""Replay the complete length-48 dimension-eight marked-contraction certificate."""

import argparse
import json
from pathlib import Path
import subprocess

import numpy as np

from affine_codeword_graph import canonicalise_many
from enumerate_c48_marked_spaces import check_contractions, generate
from finite_space_census import replay_quotient
from prepare_c48_marked_inputs import prepare
from prepare_zero_fibre_input import file_digest
from space_batch import decode_bitmaps, validate_batch
from space_codec import binary_rank
from verify_finite_space_censuses import check_sector_frames, unpack_stream
from verify_protocol_cover import Evidence, require
from verify_zero_fibre_ledger import build_buckets


def verify(root, work, native=None, graph=None, marked_native=None):
    evidence = Evidence(root)
    base_path = "data/space_contractions/c48_m08"
    data = evidence.read(base_path + "/census.json")
    require(data["schema"] == "marked-c48-space-sector-v1" and (data["c"], data["m"]) == (48, 8),
            "Different marked space sector")
    evidence.digest(data["base_certificate"]["path"], data["base_certificate"]["sha256"])
    work.mkdir(parents=True, exist_ok=True)
    inputs = work / "inputs"
    domain = prepare(evidence.root, inputs)
    require(data["source_count"] == domain["source_count"] == len(data["cases"]) == 152
            and data["raw_pairs"] == domain["raw_pairs"], "Incomplete source domain")
    for expected, actual in zip(data["cases"], domain["records"], strict=True):
        require(all(expected[k] == v for k, v in actual.items()), "Different complete marked core input")
        evidence.digest(base_path + "/" + actual["task"], actual["sha256"])
    if marked_native is not None:
        regenerated = generate(inputs, marked_native, work / "marked_generation")
        for actual, expected in zip(regenerated["records"], data["cases"], strict=True):
            require(actual["marked_binary_sha256"] == expected["marked_stream"]["decoded_sha256"]
                    and actual["trace_sha256"] == expected["schreier_trace"]["decoded_sha256"]
                    and all(actual[k] == v for k, v in expected["expected_output"].items()),
                    "Fresh marked enumeration differs")
    count = raw_mass = full_mass = deficient_mass = 0
    masks = work / "candidates.bin"
    with masks.open("xb") as output:
        for case in data["cases"]:
            binary = work / (Path(case["task"]).stem+".bin")
            unpack_stream(evidence, case["marked_stream"], binary)
            unpack_stream(evidence, case["schreier_trace"])
            core = tuple(map(int, (inputs / case["task"]).read_text(encoding="ascii").splitlines()[3].split()))
            full_count = deficient_count = case_full_mass = case_deficient_mass = 0
            with binary.open("rb") as stream:
                while raw := stream.read(4096*40):
                    require(len(raw) % 40 == 0, "Truncated marked support stream")
                    records = np.frombuffer(raw, dtype=np.uint8).reshape(-1, 40)
                    check_contractions(records[:, :32], core, case["multiplicity"])
                    points = decode_bitmaps(records[:, :32].copy().tobytes(), 8, 48)
                    weights = records[:, 32:].copy().view("<u8").reshape(-1)
                    require(np.all(weights > 0), "Empty marked orbit")
                    if case["proper_span_requires_output_rank_filter"]:
                        full = np.array([binary_rank(int(x) ^ int(p[0]) for x in p) == 8 for p in points], dtype=bool)
                    else:
                        full = np.ones(len(points), dtype=bool)
                    validate_batch(points[full], 8, 48)
                    output.write(records[full, :32].copy().tobytes())
                    full_count += int(np.sum(full))
                    deficient_count += int(np.sum(~full))
                    case_full_mass += int(np.sum(weights[full]))
                    case_deficient_mass += int(np.sum(weights[~full]))
            require(case["expected_output"] == dict(marked_orbits=full_count+deficient_count,
                    full_rank_orbits=full_count, full_rank_pairs=case_full_mass,
                    deficient_orbits=deficient_count, deficient_pairs=case_deficient_mass)
                    and case_full_mass+case_deficient_mass == case["raw_pairs"],
                    "Marked rank or mass accounting fails")
            count += full_count
            raw_mass += case["raw_pairs"]
            full_mass += case_full_mass
            deficient_mass += case_deficient_mass
    require(count == data["candidate_count"] and file_digest(masks) == data["candidate_bitmap_sha256"]
            and (raw_mass, full_mass, deficient_mass) == (data["raw_pairs"], data["full_rank_pairs"], data["deficient_pairs"]),
            "The complete marked union differs")
    bindings, keys = check_sector_frames(evidence.root, data)
    for path, sha in bindings.items():
        evidence.digest(path, sha)
    files = {}
    for name, row in data["affine_witness_streams"].items():
        files[name] = work / (name+".bin")
        unpack_stream(evidence, row, files[name])
    if native is not None:
        require(graph is not None, "Exact quotient replay requires an independent graph canonicaliser")
        ledger, buckets = work / "signatures.bin", work / "generated_buckets.bin"
        subprocess.run([str(native / "signature_c48_m08"), "--input", str(masks),
                        "--output", str(ledger), "--threads", "4"], capture_output=True, text=True, check=True)
        require(file_digest(ledger) == data["sorted_ledger_sha256"], "The complete signature ledger differs")
        build_buckets(ledger, buckets, 8, 48)
        require(file_digest(buckets) == data["affine_witness_streams"]["buckets"]["decoded_sha256"], "Different affine-invariant partition")
        replayed = replay_quotient(ledger, files["buckets"], files["classes"], files["assignments"], 48, 8)
        require(replayed["class_supports"] == [tuple(p) for p in data["class_supports"]]
                and replayed["member_counts"] == data["class_member_counts"]
                and replayed["candidate_count"] == count, "The exact affine quotient differs")
        forms = canonicalise_many([graph], [(8, p) for p in data["class_supports"]])
        require([tuple(r["canonical_points"]) for r in forms] == keys, "Independent canonical separation differs")
    return dict(schema="marked-c48-space-census-verification-v1", status="pass", c=48, m=8,
                core_classes=152, raw_pairs=raw_mass, full_rank_pairs=full_mass,
                deficient_pairs=deficient_mass, marked_representatives=count, affine_classes=data["class_count"],
                complete_core_domain_reconstructed=True, proper_span_core_included=True,
                all_marked_outputs_and_contractions_checked=True,
                marked_enumeration_freshly_recomputed=marked_native is not None,
                all_positive_affine_witnesses_freshly_replayed=native is not None,
                all_class_canonical_forms_freshly_recomputed=native is not None,
                complete_given_low_dimensional_base=True,
                conditional_premises=["The complete low-dimensional affine base and its full affine stabilisers",
                                      "The finite marked-orbit enumeration on the explicit source domain"],
                dependencies=evidence.bindings, is_global_completeness_certificate=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument("--native-directory", type=Path)
    parser.add_argument("--graph-native", type=Path)
    parser.add_argument("--marked-native", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.root.resolve(), args.work_directory.resolve(),
                    args.native_directory.resolve() if args.native_directory else None,
                    args.graph_native.resolve() if args.graph_native else None,
                    args.marked_native.resolve() if args.marked_native else None)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps({k:v for k,v in result.items() if k != "dependencies"}, indent=2))


if __name__ == "__main__":
    main()
