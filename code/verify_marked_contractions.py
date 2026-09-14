"""Replay a complete marked contraction case and independently validate its outputs."""

import argparse
import hashlib
import json
from pathlib import Path
import struct
import subprocess

from space_codec import validate_unital_support


def verify(domain, case, native, work):
    domain = Path(domain).resolve()
    report = json.loads((domain / "domain.json").read_bytes())
    selected = [r for r in report["records"] if r["id"] == case]
    if len(selected) != 1:
        raise ValueError("Contraction case is absent or repeated")
    row, = selected
    path = (domain / row["task"]).resolve()
    if not path.is_relative_to(domain):
        raise ValueError("External contraction input")
    payload = path.read_bytes()
    if hashlib.sha256(payload).hexdigest() != row["sha256"]:
        raise ValueError("Changed contraction input")
    lines = payload.decode("ascii").splitlines()
    core = tuple(map(int, lines[3].split()))
    work = Path(work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    summary_path, binary_path = work / "summary.json", work / "children.bin"
    if summary_path.exists() or binary_path.exists():
        raise ValueError("The output case already exists")
    subprocess.run([str(native), "--input", str(path), "--output", str(binary_path),
                    "--summary", str(summary_path)], check=True, capture_output=True)
    summary = json.loads(summary_path.read_bytes())
    if (summary.get("all_checks_pass") is not True or summary.get("count_only") is not False
            or summary["parent_id"] != case
            or summary["raw_marked_pair_count"] != row["raw_marked_pairs"]
            or summary["compatible_fiber_set_count"] != row["compatible_fibre_sets"]
            or summary["excluded_rank_deficient_pair_count"] != 0):
        raise ValueError("Native marked enumeration accounting failed")
    record = struct.Struct("<5Q")
    count, mass = 0, 0
    with binary_path.open("rb") as stream:
        while payload := stream.read(record.size):
            if len(payload) != record.size:
                raise ValueError("Truncated marked support")
            *words, orbit_size = record.unpack(payload)
            mask = sum(word << (64*i) for i, word in enumerate(words))
            support = tuple(x for x in range(256) if mask >> x & 1)
            if len(support) != 54 or orbit_size <= 0:
                raise ValueError("Invalid marked output")
            validate_unital_support(support, 8)
            singleton = tuple(x for x in range(128) if ((mask >> (2*x)) & 3) in (1, 2))
            fibres = sum(((mask >> (2*x)) & 3) == 3 for x in range(128))
            if singleton != core or fibres != row["full_fibres"]:
                raise ValueError("A child contracts to the wrong core")
            count += 1
            mass += orbit_size
    if count != summary["direction_marked_orbit_count"] or mass != row["raw_marked_pairs"]:
        raise ValueError("Marked output stream does not account for the input domain")
    return dict(schema="marked-contraction-case-replay-v1", status="pass", case=case,
                input_sha256=row["sha256"], marked_representatives=count, raw_marked_pairs=mass,
                complete_fibre_orbits=summary["fiber_orbit_count"],
                every_output_rank_moments_and_contraction_checked=True,
                orbit_sizes_from_native_enumerator=True,
                is_global_completeness_certificate=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--domain", type=Path, default=Path(__file__).resolve().parents[1] / "data/space_contractions/c54_m08")
    parser.add_argument("--case", default="c52_m07_000000169")
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.domain, args.case, args.native, args.work_directory)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="ascii")
    print(json.dumps(result, indent=2))
