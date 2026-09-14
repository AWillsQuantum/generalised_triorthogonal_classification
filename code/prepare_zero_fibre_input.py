"""Convert compact parent intervals to accelerated zero-fibre kernel inputs."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import struct
import subprocess

from selected_spaces import iter_interval
from space_codec import binary_rank, checked_support
from space_lifts import evaluation_rows
from verify_protocol_cover import Evidence, require


PROFILE = struct.Struct("<BBH")
SOURCE_SUFFIX = struct.Struct("<4B4x3Q")


def bitmap(points, dimension, length=54):
    checked_support(points, dimension, length)
    return sum(1 << x for x in points).to_bytes((1 << dimension)//8, "little")


def source_record(points, dimension, source_index, quadratic_rank, lift_dimension, length=54):
    require(10 <= dimension <= 14 and 0 <= source_index < 1 << 64,
            "Invalid zero-fibre source parameters")
    require(length in (50, 52, 54) and 0 <= lift_dimension < length and quadratic_rank+lift_dimension+dimension+1 == length,
            "Invalid zero-fibre rank-nullity profile")
    return bitmap(points, dimension, length)+SOURCE_SUFFIX.pack(0, 0, quadratic_rank, lift_dimension,
                                                       source_index, 1, 1 << lift_dimension)


def prepare(root, dimension, first, count, native, work, length=54):
    require(11 <= dimension <= 15 and first >= 0 and count >= 1, "Invalid native source interval")
    require(length == 54 or (length in (50, 52) and dimension == 11), "Unsupported native graph profile configuration")
    evidence = Evidence(root)
    index = evidence.read("data/spaces/index.json")
    sectors = [s for row in index["lengths"] if row["c"] == length
               for s in row["sectors"] if s["m"] == dimension-1]
    require(len(sectors) == 1, "The parent sector is not uniquely declared")
    sector, = sectors
    require(sector["data_complete"] and sector["included_classes"] == sector["expected_classes"]
            and first+count <= sector["included_classes"], "The source interval is outside the complete parent data")
    for item in sector["encoding_certificates"]:
        evidence.digest(item["path"], item["sha256"])
    work.mkdir(parents=True, exist_ok=True)
    inputs, profiles, summary, sources = [work / name for name in
                                          ("profile_input.bin", "profiles.bin", "profile_summary.json", "sources.bin")]
    require(not any(p.exists() for p in (inputs, profiles, summary, sources)), "A native output path already exists")
    # The rank kernel ignores the signature and requires a positive member
    # count. This neutral wrapper is not a signature or class-mass ledger.
    bindings, samples = {}, {}
    sample_indices = {first + i*(count-1)//127 for i in range(128)}
    with inputs.open("xb") as stream:
        for i, points in iter_interval(evidence.root, length, dimension-1, first, count, bindings):
            stream.write(bytes(32)+bitmap(points, dimension-1, length)+struct.pack("<Q", 1))
            if i in sample_indices:
                samples[i] = binary_rank(evaluation_rows(points, dimension-1, 2))
    for name, digest in bindings.items():
        evidence.digest(name, digest)
    executable = native / (f"m{dimension:02d}_source_profile_kernel" if length == 54 else f"profile_c{length}_m{dimension:02d}")
    subprocess.run([str(executable), "--input", str(inputs), "--output", str(profiles), "--summary", str(summary),
                    "--first", "0", "--count", str(count)], check=True)
    raw = profiles.read_bytes()
    require(len(raw) == count*PROFILE.size, "The native source profile is incomplete")
    histogram = Counter()
    total_lifts = 0
    mask_bytes = (1 << (dimension-1))//8
    with sources.open("xb") as stream, inputs.open("rb") as incoming:
        for offset, (rank, lift_dimension, fibre_count) in enumerate(PROFILE.iter_unpack(raw)):
            require(fibre_count == 1, "A zero-fibre profile has a different fibre count")
            i = first+offset
            require(rank+lift_dimension+dimension == length, "Invalid native rank-nullity profile")
            if i in samples:
                require(rank == samples[i], "The independent Python quadratic rank disagrees")
            wrapped = incoming.read(32+mask_bytes+8)
            require(len(wrapped) == 32+mask_bytes+8 and wrapped[:32] == bytes(32)
                    and wrapped[-8:] == struct.pack("<Q", 1), "The neutral source wrapper has changed")
            stream.write(wrapped[32:-8]+SOURCE_SUFFIX.pack(0, 0, rank, lift_dimension, i, 1, 1 << lift_dimension))
            histogram[lift_dimension] += 1
            total_lifts += (1 << lift_dimension)-1
        require(not incoming.read(1), "Trailing neutral source records")
    report = json.loads(summary.read_bytes())
    require(report["source_count"] == count and int(report["full_rank_marked_extension_count"]) == total_lifts,
            "Native source-profile accounting does not close")
    source = (f"code/space_native/length54/m{dimension:02d}_source_profile_kernel.cpp" if length == 54
              else "code/space_native/generic/profile_m11.cpp")
    if length != 54:
        evidence.digest("code/space_native/generic/configuration.hpp")
    return dict(schema="native-zero-fibre-input-v1", status="pass", c=length, m=dimension,
        source_interval=[first, first+count], complete_parent_classes=sector["included_classes"],
        complete_parent_sector_included=first == 0 and count == sector["included_classes"],
        source_count=count, nonaffine_graph_lifts=total_lifts, lift_dimension_histogram=dict(sorted(histogram.items())),
        source_record_bytes=(1 << (dimension-1))//8+SOURCE_SUFFIX.size,
        source_record_format="little-endian ambient bitmap, then <4B4x3Q",
        source_suffix_fields=["source_kind=0", "multiplicity=0", "quadratic_rank", "lift_dimension",
                              "source_index", "compatible_fibre_count=1", "raw_lift_count"],
        source_sha256=file_digest(sources),
        native_profile_sha256=hashlib.sha256(raw).hexdigest(),
        native_source_sha256=evidence.digest(source),
        independent_python_rank_samples=len(samples), bounded_memory_catalogue_stream=True,
        dependencies=evidence.bindings)


def file_digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--dimension", type=int, required=True)
    parser.add_argument("--length", type=int, default=54, choices=(50, 52, 54))
    parser.add_argument("--first", type=int, default=0)
    parser.add_argument("--count", type=int, required=True)
    parser.add_argument("--native-directory", type=Path, required=True)
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.root, args.dimension, args.first, args.count,
                     args.native_directory.resolve(), args.work_directory.resolve(), args.length)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps({k: v for k, v in result.items() if k != "dependencies"}, indent=2))


if __name__ == "__main__":
    main()
