"""Re-enumerate a finite list of pointed supports with the native protocol kernel."""

import argparse
import hashlib
import json
from pathlib import Path
import struct
import subprocess

from compact_to_native import BASE_RECORD, MANIFEST_HEADER
from space_codec import binary_rank, checked_support
from space_lifts import evaluation_rows


def write_pointed_cases(stream, data):
    cases = data["cases"]
    if not 0 < len(cases) < 1 << 32:
        raise ValueError("Invalid finite case count")
    if [r["index"] for r in cases] != list(range(len(cases))):
        raise ValueError("Incomplete case indices")
    length_limit = data["maximum_protocol_length"]
    if not 0 < length_limit <= 64:
        raise ValueError("Invalid length limit")
    stream.write(MANIFEST_HEADER.pack(b"UTSPTS1\0", 1, len(cases), len(cases), length_limit, 1))
    for row in cases:
        for value in (f"case_{row['index']}".encode("ascii"), b"pointed_support"):
            stream.write(struct.pack("<H", len(value)))
            stream.write(value)
    for row in cases:
        m = row["ambient_dimension"]
        points = checked_support(row["points"], m)
        n = len(points)
        if not 0 < n <= length_limit or binary_rank(points) != m:
            raise ValueError("Pointed support lacks full linear rank or has invalid length")
        if any(row.bit_count() & 1 for row in evaluation_rows(points, m, 3)[1:]):
            raise ValueError("Pointed stabiliser support is not triorthogonal")
        # These are already pointed supports: no origin expansion is performed.
        stream.write(BASE_RECORD.pack(row["index"], n, m, m, 0, n, 0))
        stream.write(struct.pack(f"<{n}I", *points))


def normalized_keys(keys):
    result = set()
    for key in keys:
        if not key:
            raise ValueError("An output key must contain at least one word")
        words = tuple(int(word, 0) if isinstance(word, str) else word for word in key)
        if any(type(word) is not int or not 0 <= word < 1 << 64 for word in words):
            raise ValueError("Output words must be unsigned 64-bit integers")
        result.add(words)
    return result


def check_result(data, result):
    if (result["status"] != "complete" or result["matrix_scope"] != "full_projective"
            or result["logical_qubits"] != data["logical_qubits"]
            or result["minimum_distance"] != data["minimum_distance"]
            or result["enumeration_mode"] != "raw_isotropic_subspaces"
            or result["support_range"] != dict(start=0, count=len(data["cases"]), end_exclusive=len(data["cases"]))):
        raise ValueError("Incomplete or mismatched native census")
    stats = result["statistics"]
    if (stats["supports_processed"] != len(data["cases"])
            or not stats["isotropic_subspace_statistics_exact"] or not stats["radical_statistics_exact"]):
        raise ValueError("Native counts are not exact")
    for name in ("isotropic_subspaces", "nondegenerate_subspaces"):
        if stats[name] != sum(row["expected"][name] for row in data["cases"]):
            raise ValueError(f"Native census disagrees on {name}")
    actual = {}
    for row in result["positive_supports"]:
        index = row["manifest_support_index"]
        if (index in actual or not 0 <= index < len(data["cases"])
                or row["original_record_index"] != index
                or any(key["logical_qubits"] != data["logical_qubits"] for key in row["output_keys"])):
            raise ValueError("Invalid positive-support identity")
        actual[index] = normalized_keys(key["canonical_key_words"] for key in row["output_keys"])
    for row in data["cases"]:
        if actual.get(row["index"], set()) != normalized_keys(row["expected"]["output_keys"]):
            raise ValueError("Output profile differs from complete finite-sector evidence")
    return dict(status="pass", cases=len(data["cases"]), logical_qubits=data["logical_qubits"],
                minimum_distance=data["minimum_distance"],
                maximum_protocol_length=data["maximum_protocol_length"],
                isotropic_subspaces=stats["isotropic_subspaces"],
                nondegenerate_subspaces=stats["nondegenerate_subspaces"],
                every_support_output_profile_recomputed=True,
                is_global_completeness_certificate=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("Workers must be positive")
    raw = args.input.read_bytes()
    data = json.loads(raw)
    work = args.work_directory.resolve()
    work.mkdir(parents=True, exist_ok=True)
    manifest = work / "supports.utsp"
    with manifest.open("wb") as stream:
        write_pointed_cases(stream, data)
    output = work / "census.json"
    subprocess.run([str(args.native), "manifest-catalogue", "--input", str(manifest),
                    "--output", str(output), "--q", str(data["logical_qubits"]),
                    "--minimum-distance", str(data["minimum_distance"]),
                    "--emit-positive-supports", "--support-workers", str(args.workers)], check=True)
    result = check_result(data, json.loads(output.read_bytes()))
    result["input_sha256"] = hashlib.sha256(raw).hexdigest()
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
