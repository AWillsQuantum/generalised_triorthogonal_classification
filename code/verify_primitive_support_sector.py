"""Check the finite primitive-target domain and optionally reproduce its census."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess

from logical_spaces import logical_label_space
from space_codec import binary_rank, checked_support
from space_lifts import evaluation_rows
from verify_protocol_case_census import write_pointed_cases


def support_key(points, h, words):
    n = len(points)
    values = [int(word, 16) for word in words]
    if (len(values) != 1+2*h or values[0] != n | h << 8
            or values[1:1+h] != values[1+h:] or any(row >> n for row in values[1:])):
        raise ValueError("Invalid unmarked support key")
    recovered = [sum((values[1+i] >> j & 1) << i for i in range(h)) for j in range(n)]
    if recovered != points or binary_rank(values[1:1+h]) != h:
        raise ValueError("Support key does not reconstruct its full-rank input")
    payload = json.dumps([n, h, words], separators=(",", ":")).encode("ascii")
    return hashlib.sha256(payload).hexdigest()


def validate_zero_result(result):
    expected = dict(logical_dimension=5, target_dimension=5, target_nondegenerate_seed=3,
                    target_signature_count=2, used_complete_primitive_target_recognizer=True,
                    sector="all", subspace_orbits=0, weighted_subspace_count=0)
    if any(result.get(key) != value for key, value in expected.items()):
        raise ValueError("Not a complete primitive-target zero result")
    for key in ("representatives", "canonical_keys", "orbit_sizes"):
        if result.get(key, []) != []:
            raise ValueError("A zero census has nonempty orbit data")


def read_data(path):
    raw = path.read_bytes()
    data = json.loads(raw)
    if (data["schema"] != "finite-primitive-support-sector-v1"
            or data["maximum_protocol_length"] != 54 or data["minimum_distance"] != 3
            or data["logical_qubits"] != 5 or data["target_nondegenerate_seed"] != 3
            or data["target_signatures"] != ["0x60000", "0x60001"]):
        raise ValueError("Wrong finite primitive-sector scope")
    return raw, data


def verify(data):
    cases = data["cases"]
    keys, dimensions, lengths, groups = set(), Counter(), Counter(), Counter()
    enumeration_variants = []
    for index, row in enumerate(cases):
        points, h = row["points"], row["ambient_dimension"]
        checked_support(sorted(points), h)
        if (row["index"] != index or not 0 < len(points) <= 54 or 0 in points
                or row["support_length"] != len(points)):
            raise ValueError("Invalid normalised support domain")
        key = support_key(points, h, row["key_words"])
        if key != row["canonical_key_sha256"] or key in keys:
            raise ValueError("Repeated or inconsistent support key")
        keys.add(key)
        order = row["automorphism_group_order"]
        if type(order) is not int or order <= 0:
            raise ValueError("Invalid support automorphism order")
        qdim = logical_label_space(sorted(points), h, 3)["quotient_dimension"]
        d4 = logical_label_space(sorted(points), h, 4)["quotient_dimension"]
        result = row["result"]
        if (result["quotient_dimension"] != qdim or result["subspace_orbits"] != 0
                or result["weighted_subspace_count"] != 0
                or any(type(value) is not int or value < 0 for value in result.values())):
            raise ValueError("Invalid primitive-target counts or label dimension")
        if "enumerated_points" in row:
            variant = row["enumerated_points"]
            checked_support(sorted(variant), h)
            if (support_key(variant, h, row["enumerated_key_words"])
                    != row["enumerated_canonical_key_sha256"]
                    or logical_label_space(sorted(variant), h, 3)["quotient_dimension"] != qdim):
                raise ValueError("Invalid equivalent enumeration input")
            enumeration_variants.append(dict(index=len(enumeration_variants), case_index=index,
                                             ambient_dimension=h, points=sorted(x for x in variant if x)))
        dimensions[qdim, d4] += 1
        lengths[len(points), h] += 1
        groups[order] += 1
    seen = set()
    for index, row in enumerate(data["normalisation_inputs"]):
        if row["index"] != index or not 0 <= row["case_index"] < len(cases):
            raise ValueError("Incomplete input normalisation map")
        points, h = row["points"], row["ambient_dimension"]
        checked_support(points, h)
        if (0 in points or binary_rank(points) != h
                or any(mask.bit_count() & 1 for mask in evaluation_rows(points, h, 3)[1:])):
            raise ValueError("Invalid normalised input geometry")
        target = cases[row["case_index"]]
        if (len(points), h) != (target["support_length"], target["ambient_dimension"]):
            raise ValueError("Normalisation changes length or stabiliser dimension")
        if row["family"] not in ("c52_m07_pointings", "larger_quotient_supports"):
            raise ValueError("Unknown finite input family")
        seen.add(row["case_index"])
    if seen != set(range(len(cases))):
        raise ValueError("A listed support is not covered by the input map")
    return dict(schema="finite-primitive-support-verification-v1", status="pass", cases=len(cases),
                normalisation_inputs=len(data["normalisation_inputs"]),
                enumeration_input_variants=len(enumeration_variants),
                distance_filtration=[dict(d3=a, d4=b, count=count) for (a, b), count in sorted(dimensions.items())],
                length_dimension_counts=[dict(n=n, h=h, count=count) for (n, h), count in sorted(lengths.items())],
                all_geometries_and_quotient_dimensions_independently_checked=True,
                normalisation_equivalences_freshly_recomputed=False,
                primitive_censuses_freshly_recomputed=False, is_global_completeness_certificate=False)


def canonical_check(data, native, work, workers):
    records = [dict(index=i, points=row["points"], ambient_dimension=row["ambient_dimension"],
                    case_index=row["case_index"]) for i, row in enumerate(data["normalisation_inputs"])]
    for row in data["cases"]:
        if "enumerated_points" in row:
            records.append(dict(index=len(records), case_index=row["index"],
                                points=sorted(x for x in row["enumerated_points"] if x),
                                ambient_dimension=row["ambient_dimension"]))
    manifest = work / "normalised_inputs.utsp"
    with manifest.open("wb") as stream:
        write_pointed_cases(stream, dict(maximum_protocol_length=54, cases=records))
    output = work / "canonical_keys.json"
    subprocess.run([str(native), "manifest-support-canonical-keys", "--input", str(manifest),
                    "--output", str(output), "--workers", str(workers)], check=True)
    result = json.loads(output.read_bytes())
    if (result["status"] != "complete" or len(result["records"]) != len(records)
            or result["support_range"] != dict(start=0, count=len(records), end_exclusive=len(records))):
        raise ValueError("Incomplete canonical support census")
    for row, actual in zip(records, result["records"], strict=True):
        expected = data["cases"][row["case_index"]]
        if (actual["support_index"] != row["index"] or actual["source_index"] != row["index"]
                or actual["original_record_index"] != row["index"]
                or actual["support_length"] != len(row["points"])
                or actual["ambient_dimension"] != row["ambient_dimension"]
                or actual["key_words"] != expected["key_words"]
                or actual["automorphism_group_order"] != expected["automorphism_group_order"]):
            raise ValueError("Input geometry is not equivalent to its assigned support")
    return dict(normalisation_equivalences_freshly_recomputed=True,
                native_canonical_inputs_checked=len(records))


def replay(data, native, start, count, workers, work):
    if start < 0 or count < 1 or start+count > len(data["cases"]):
        raise ValueError("Invalid finite replay interval")
    checked = []
    for row in data["cases"][start:start+count]:
        result = subprocess.run([str(native), "canonical-isotropic-orbits", "--ambient",
            str(row["ambient_dimension"]), "--points", ",".join(map(str, row["points"])),
            "--q", "5", "--target-signatures", "0x60000,0x60001", "--target-nondegenerate-seed", "3",
            "--workers", str(workers), "--emit", "100000", "--emit-orbit-data"],
            capture_output=True, text=True, check=True)
        raw = json.loads(result.stdout)
        validate_zero_result(raw)
        if raw["quotient_dimension"] != row["result"]["quotient_dimension"]:
            raise ValueError("Replayed quotient dimension differs")
        checked.append(row["index"])
        (work / f"primitive_case_{row['index']:05d}.json").write_text(json.dumps(dict(
            index=row["index"], status="pass", subspace_orbits=0, weighted_subspace_count=0,
            quotient_dimension=raw["quotient_dimension"]), indent=2)+"\n", encoding="ascii")
    return dict(primitive_censuses_freshly_recomputed=len(checked) == len(data["cases"]),
                primitive_replay_indices=checked)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=root / "data/protocol_sectors/primitive_support_union.json")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--native", type=Path)
    parser.add_argument("--work-directory", type=Path)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--canonical-check", action="store_true")
    parser.add_argument("--replay-start", type=int, default=0)
    parser.add_argument("--replay-count", type=int, default=0)
    args = parser.parse_args()
    if args.workers < 1 or args.replay_count < 0:
        parser.error("Invalid workers or replay count")
    if (args.canonical_check or args.replay_count) and (args.native is None or args.work_directory is None):
        parser.error("Native checks require --native and --work-directory")
    raw, data = read_data(args.input)
    result = verify(data)
    if args.work_directory:
        args.work_directory.mkdir(parents=True, exist_ok=True)
    if args.canonical_check:
        result.update(canonical_check(data, args.native, args.work_directory, args.workers))
    if args.replay_count:
        result.update(replay(data, args.native, args.replay_start, args.replay_count,
                             args.workers, args.work_directory))
    result["input_sha256"] = hashlib.sha256(raw).hexdigest()
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
