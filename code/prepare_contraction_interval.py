"""Prepare accelerated m=9,10 inverse-contraction inputs from compact intervals."""

import argparse
from collections import Counter
import json
from pathlib import Path
import struct
import subprocess

from prepare_space_contractions import CoreExtensionContext, count_zero_xor_index_subsets
from prepare_zero_fibre_input import file_digest, PROFILE, SOURCE_SUFFIX
from selected_spaces import iter_interval
from verify_protocol_cover import Evidence, require


def reference_profile(points, dimension, multiplicity):
    context = CoreExtensionContext.build(dimension-1, points, dimension)
    return (context.solver.rank, len(context.lift_complement),
            count_zero_xor_index_subsets(context.external_labels, multiplicity))


def prepare(root, dimension, multiplicity, first, count, native, work, target_length=54):
    require(dimension in (9, 10) and multiplicity in range(12-dimension)
            and first >= 0 and count > 0 and target_length in (50, 52, 54), "Unsupported contraction interval")
    root, work, native = Path(root).resolve(), Path(work).resolve(), Path(native).resolve()
    require(not work.exists(), "Use a new contraction input directory")
    work.mkdir(parents=True)
    evidence = Evidence(root)
    length = target_length-2*multiplicity
    mask_bytes = (1 << (dimension-1))//8
    input_path, output_path = work / "parents.bin", work / ("sources.bin" if dimension == 9 else "profiles.bin")
    sample_indices = {first+i*(count-1)//127 for i in range(128)}
    samples, bindings = {}, {}
    with input_path.open("xb") as stream:
        for i, points in iter_interval(root, length, dimension-1, first, count, bindings):
            mask = sum(1 << x for x in points).to_bytes(mask_bytes, "little")
            stream.write(mask if dimension == 9 else bytes(32)+mask+struct.pack("<Q", 1))
            if i in sample_indices:
                samples[i] = reference_profile(points, dimension, multiplicity)
    for name, value in bindings.items():
        evidence.digest(name, value)
    if dimension == 9:
        stem = "m09_profile_interval" if target_length == 54 else "m09_configured_profile_interval"
        source = f"code/space_native/contractions/{stem}.cpp"
        evidence.digest("code/space_native/length54/m09_source_profile_kernel.cpp")
        command = [str(native / stem), str(input_path), str(output_path), str(multiplicity), str(first)]
        if target_length != 54:
            command.append(str(target_length))
        subprocess.run(command, check=True, capture_output=True)
        width = mask_bytes+SOURCE_SUFFIX.size
    else:
        source = "code/space_native/length54/m10_direct_source_profile.cpp"
        subprocess.run([str(native / "m10_direct_source_profile"), "--input", str(input_path), "--output", str(output_path),
                        "--summary", str(work / "summary.json"), "--first", "0", "--count", str(count),
                        "--multiplicity", str(multiplicity), "--expected-weight", str(length),
                        "--target-length", str(target_length)], check=True, capture_output=True)
        width = PROFILE.size
    require(output_path.stat().st_size == count*width, "The native profile interval is incomplete")
    histogram = Counter()
    compatible = full_rank = contributors = 0
    with output_path.open("rb") as stream, input_path.open("rb") as parents:
        for offset in range(count):
            i = first+offset
            payload = stream.read(width)
            if dimension == 9:
                kind, n, rank, lift, index, fibres, raw = SOURCE_SUFFIX.unpack(payload[mask_bytes:])
                require(payload[:mask_bytes] == parents.read(mask_bytes) and kind == n == multiplicity
                        and index == i and raw == fibres*(1 << lift), "Invalid native source record")
            else:
                rank, lift, fibres = PROFILE.unpack(payload)
            require(rank+lift+dimension == length and 0 <= lift < 54 and fibres >= 0,
                    "The exact profile violates rank-nullity")
            if i in samples:
                require((rank, lift, fibres) == samples[i], "The independent reference profile disagrees")
            require(multiplicity != 0 or fibres == 1, "A graph source has a different fibre count")
            candidates = fibres*(1 << lift)-int(multiplicity == 0)
            require(candidates >= 0, "Negative full-rank candidate count")
            compatible += fibres
            full_rank += candidates
            contributors += candidates > 0
            histogram[rank, lift] += 1
    return dict(schema="native-contraction-interval-input-v1", status="pass", c=target_length, m=dimension,
        multiplicity=multiplicity, parent_length=length, source_interval=[first, first+count], source_count=count,
        compatible_fibre_sets=compatible, full_rank_marked_lifts=full_rank, contributing_sources=contributors,
        profiles=[dict(quadratic_rank=rank, lift_dimension=lift, sources=n) for (rank, lift), n in sorted(histogram.items())],
        independent_reference_profiles=len(samples), bounded_memory_catalogue_stream=True,
        native_input_sha256=file_digest(input_path), native_output_sha256=file_digest(output_path),
        native_source_sha256=evidence.digest(source), dependencies=evidence.bindings)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--dimension", type=int, required=True)
    parser.add_argument("--length", type=int, default=54, choices=(50, 52, 54))
    parser.add_argument("--multiplicity", type=int, required=True)
    parser.add_argument("--first", type=int, default=0)
    parser.add_argument("--count", type=int, required=True)
    parser.add_argument("--native-directory", type=Path, required=True)
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.root, args.dimension, args.multiplicity, args.first, args.count,
                     args.native_directory, args.work_directory, args.length)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps({k: v for k, v in result.items() if k != "dependencies"}, indent=2))


if __name__ == "__main__":
    main()
