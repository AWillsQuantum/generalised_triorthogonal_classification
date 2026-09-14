"""Verify the larger-quotient source partition, geometries and exact equivalences."""

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import subprocess

from logical_spaces import logical_label_space
from selected_spaces import read_selected
from space_codec import checked_support
from verify_primitive_support_sector import support_key
from verify_protocol_case_census import write_pointed_cases


def pointed_support(points, m, parity, origin):
    points = checked_support(points, m)
    if parity in (0, 1, 2):
        if type(origin) is not int or not 0 <= origin < 1 << m:
            raise ValueError("Invalid source origin")
        result = sorted(x ^ origin for x in points)
        if parity == 1:
            if not result or result[0] != 0:
                raise ValueError("Removed origin does not belong to the source")
            result = result[1:]
        elif parity == 2:
            if 0 in result:
                raise ValueError("Added origin already belongs to the source")
            result = [0]+result
        return result, m
    if parity in (3, 4):
        result = [x | 1 << m for x in points]
        return ([0]+result if parity == 4 else result), m+1
    raise ValueError("Unknown pointing case")


def verify(root, data, primitive, blocks):
    if (data["schema"] != "finite-large-quotient-domain-v1" or data["maximum_protocol_length"] != 54
            or data["logical_qubits"] != 5 or data["minimum_distance"] != 3
            or data["minimum_quotient_dimension"] != 15):
        raise ValueError("Wrong larger-quotient scope")
    if (primitive["schema"] != "finite-primitive-support-sector-v1"
            or blocks["schema"] != "finite-low-dimensional-protocol-blocks-v1"):
        raise ValueError("Wrong supporting datasets")
    classes = data["classes"]
    class_keys, normalised_cases = set(), set()
    for index, row in enumerate(classes):
        h, n = row["ambient_dimension"], row["support_length"]
        values = [int(word, 16) for word in row["key_words"]]
        if len(values) != 1+2*h:
            raise ValueError("Invalid canonical key length")
        points = [sum((values[1+i] >> j & 1) << i for i in range(h)) for j in range(n)]
        checked_support(sorted(points), h)
        key = support_key(points, h, row["key_words"])
        if row["index"] != index or key != row["canonical_key_sha256"] or key in class_keys:
            raise ValueError("Invalid larger-quotient class partition")
        class_keys.add(key)
        mapped = primitive["normalisation_inputs"][row["normalisation_input_index"]]
        case = row["primitive_case_index"]
        if (mapped["case_index"] != case or mapped["family"] != "larger_quotient_supports"
                or mapped["ambient_dimension"] != h or mapped["points"] != sorted(x for x in points if x)
                or mapped["deleted_zero_column"] != (0 in points)):
            raise ValueError("Larger-quotient class missing from primitive input cover")
        normalised_cases.add(case)
    expected_by_block = {row["index"]: row["q5_large_quotient_supports"] for row in blocks["blocks"]}
    selected = defaultdict(list)
    multiplicities, dimensions = Counter(), Counter()
    keys = {(row["space"]["c"], row["space"]["m"], row["space"]["index"]) for row in data["pointings"]}
    supports, source_bindings = read_selected(root, keys)
    seen = set()
    for index, row in enumerate(data["pointings"]):
        block_index = row["block_index"]
        block = blocks["blocks"][block_index]
        source = row["space"]
        key = source["c"], source["m"], source["index"]
        if (row["index"] != index or block["index"] != block_index
                or (source["c"], source["m"]) != (block["c"], block["m"])
                or not block["source_interval"][0] <= source["index"] < block["source_interval"][1]):
            raise ValueError("Pointing does not belong to its claimed source block")
        expected, h = pointed_support(supports[key], source["m"], row["parity_case"], source["origin"])
        if row["points"] != expected or row["ambient_dimension"] != h:
            raise ValueError("Pointing does not reconstruct from compact source")
        identity = key, h, tuple(expected)
        if identity in seen:
            raise ValueError("Duplicate selected pointing")
        seen.add(identity)
        cls = classes[row["support_class_index"]]
        if (cls["support_length"], cls["ambient_dimension"]) != (len(expected), h):
            raise ValueError("Pointed class has wrong geometric parameters")
        dimension = logical_label_space(expected, h, 3)["quotient_dimension"]
        if dimension < 15 or dimension != primitive["cases"][cls["primitive_case_index"]]["result"]["quotient_dimension"]:
            raise ValueError("Incorrect larger-quotient selection")
        selected[block_index].append(row["block_selected_index"])
        multiplicities[row["support_class_index"]] += 1
        dimensions[dimension] += 1
    if set(multiplicities) != set(range(len(classes))):
        raise ValueError("Listed canonical class is not represented")
    for block_index, count in expected_by_block.items():
        if sorted(selected.get(block_index, [])) != list(range(count)):
            raise ValueError("The selected inputs do not cover their block census")
    if len(data["pointings"]) != sum(expected_by_block.values()):
        raise ValueError("Incorrect global larger-quotient cardinality")
    return dict(schema="finite-large-quotient-domain-verification-v1", status="pass",
                selected_pointings=len(data["pointings"]), support_classes=len(classes),
                normalised_classes=len(normalised_cases), nonempty_source_blocks=len(selected),
                distinct_source_spaces=len(keys), quotient_dimensions=dict(sorted(dimensions.items())),
                every_pointing_reconstructed_from_compact_catalogue=True,
                every_selected_quotient_dimension_independently_recomputed=True,
                block_cardinalities_match_complete_quotient_censuses=True,
                every_canonical_class_in_primitive_input_cover=True,
                exact_support_equivalences_freshly_recomputed=False,
                full_origin_censuses_freshly_recomputed=False,
                source_shard_bindings=source_bindings, is_global_completeness_certificate=False)


def canonical_check(data, native, work, workers):
    manifest = work / "large_quotient_inputs.utsp"
    with manifest.open("wb") as stream:
        write_pointed_cases(stream, dict(maximum_protocol_length=54, cases=data["pointings"]))
    output = work / "large_quotient_keys.json"
    subprocess.run([str(native), "manifest-support-canonical-keys", "--input", str(manifest),
                    "--output", str(output), "--workers", str(workers)], check=True)
    result = json.loads(output.read_bytes())
    count = len(data["pointings"])
    if (result["status"] != "complete" or len(result["records"]) != count
            or result["support_range"] != dict(start=0, count=count, end_exclusive=count)):
        raise ValueError("Incomplete support canonicalisation")
    for row, actual in zip(data["pointings"], result["records"], strict=True):
        expected = data["classes"][row["support_class_index"]]
        if (actual["support_index"] != row["index"] or actual["source_index"] != row["index"]
                or actual["original_record_index"] != row["index"]
                or actual["support_length"] != expected["support_length"]
                or actual["ambient_dimension"] != expected["ambient_dimension"]
                or actual["key_words"] != expected["key_words"]
                or actual["automorphism_group_order"] != expected["automorphism_group_order"]):
            raise ValueError("Selected input is not in its recorded exact support class")
    return dict(exact_support_equivalences_freshly_recomputed=True,
                canonical_inputs_checked=count)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=root / "data/protocol_sectors/large_quotient_domain.json")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--native", type=Path)
    parser.add_argument("--work-directory", type=Path)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    if args.workers < 1 or (args.native is None) != (args.work_directory is None):
        parser.error("Native canonicalisation requires an executable, work directory and positive workers")
    raw = args.input.read_bytes()
    data = json.loads(raw)
    paths = [root / data[key] for key in ("primitive_dataset", "block_dataset")]
    if any(not path.resolve().is_relative_to(root) for path in paths):
        raise ValueError("Dataset reference escapes the release")
    supporting = [path.read_bytes() for path in paths]
    result = verify(root, data, *map(json.loads, supporting))
    if args.native:
        args.work_directory.mkdir(parents=True, exist_ok=True)
        result.update(canonical_check(data, args.native, args.work_directory, args.workers))
    result["input_sha256"] = hashlib.sha256(raw).hexdigest()
    result["supporting_data_sha256"] = {path.relative_to(root).as_posix(): hashlib.sha256(content).hexdigest()
                                        for path, content in zip(paths, supporting, strict=True)}
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps({k: v for k, v in result.items() if k != "source_shard_bindings"}, indent=2))


if __name__ == "__main__":
    main()
