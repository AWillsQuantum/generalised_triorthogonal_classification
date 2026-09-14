"""Compare a complete accelerated zero-fibre lift family with reference equations."""

import argparse
import hashlib
import json
from pathlib import Path
import struct
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent / "space_algorithms"))
from inverse_contraction_extensions import CoreExtensionContext
from space_codec import validate_unital_support
from space_lifts import evaluation_rows, iter_lifts, reduce_vector


def verify(data, native, work):
    m, c, points = data["m"], data["c"], tuple(data["points"])
    if (m, c) != (8, 54):
        raise ValueError("This fixed-width kernel takes dimension-eight length-54 cores")
    validate_unital_support(points, m)
    context = CoreExtensionContext.build(m, points)
    lift_dimension = len(context.lift_complement)
    if lift_dimension > 16:
        raise ValueError("Reference comparison is bounded to 65,536 lift classes")
    expected_count = (1 << lift_dimension) - 1
    work.mkdir(parents=True, exist_ok=True)
    source = work / "source.bin"
    target = work / "lifts.bin"
    mask = sum(1 << x for x in points)
    words = tuple((mask >> (64*i)) & ((1 << 64)-1) for i in range(4))
    source.write_bytes(struct.pack("<4Q4BI3Q", *words, 0, 0, context.solver.rank,
                                   lift_dimension, 0, 0, 1, 1 << lift_dimension))
    if target.exists():
        raise ValueError("Lift output already exists")
    result = subprocess.run([str(native), "--input", str(source), "--output", str(target),
                             "--source-start", "0", "--source-count", "1",
                             "--expected-multiplicity", "0"], text=True, capture_output=True, check=True)
    summary = json.loads(result.stdout)
    if (summary["status"] != "complete_length54_m09_minimum_filter_batch_v1"
            or summary["sources_completed"] != 1 or summary["candidate_limit"] != 0
            or summary["processed_candidate_count"] != expected_count
            or summary["retained_distinct_batch_mask_count"] != expected_count):
        raise ValueError("Native lift traversal did not cover the exact family")
    affine_pivots = {}
    for row in evaluation_rows(points, m, 1):
        value = reduce_vector(row, affine_pivots)
        if value:
            affine_pivots[value.bit_length()-1] = value
    def class_key(support):
        if len(support) != c or sorted(x >> 1 for x in support) != list(points):
            raise ValueError("Lift does not have the prescribed singleton fibres")
        labels = sum((x & 1) << i for i, x in enumerate(sorted(support)))
        return reduce_vector(labels, affine_pivots)
    payload = target.read_bytes()
    if len(payload) != 64 * expected_count:
        raise ValueError("Incorrect native lift output length")
    actual = set()
    for offset in range(0, len(payload), 64):
        bits = int.from_bytes(payload[offset:offset+64], "little")
        support = tuple(x for x in range(512) if bits >> x & 1)
        validate_unital_support(support, 9)
        actual.add(class_key(support))
    expected = {class_key(support) for support in iter_lifts(points, (), m)}
    if len(expected) != expected_count or actual != expected:
        raise ValueError("Native and reference affine-shear classes disagree")
    return dict(status="pass", core_length=c, core_affine_dimension=m,
                target_affine_dimension=m+1, fibre_multiplicity=0,
                lift_quotient_dimension=lift_dimension, compared_full_rank_lifts=expected_count,
                all_native_outputs_independently_validated=True,
                equality_modulo_affine_shears=True,
                is_global_completeness_certificate=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path,
                        default=Path(__file__).resolve().parents[1] / "data/examples/inverse_contraction_core.json")
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    raw = args.input.read_bytes()
    report = verify(json.loads(raw), args.native, args.work_directory.resolve())
    report["input_sha256"] = hashlib.sha256(raw).hexdigest()
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
