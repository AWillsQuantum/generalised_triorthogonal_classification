"""Regenerate a bounded interval of accelerated length-50 or length-52 space lifts."""

import argparse
import json
from pathlib import Path
import subprocess

from prepare_contraction_interval import prepare as prepare_contraction
from prepare_zero_fibre_input import file_digest, prepare as prepare_graph
from verify_protocol_cover import Evidence, require


def generate(root, length, dimension, multiplicity, first, count, native, work,
             maximum_candidates=1000000, threads=1):
    require(length in (50, 52) and dimension in (9, 10, 11)
            and multiplicity in range(max(1, 12-dimension)) and threads > 0
            and maximum_candidates > 0, "Unsupported accelerated lift interval")
    native, work = Path(native).resolve(), Path(work).resolve()
    require(not work.exists(), "Use a new interval directory")
    work.mkdir(parents=True)
    inputs = work / "input"
    if dimension == 11:
        require(multiplicity == 0, "Only graph lifts are possible at this dimension")
        profile = prepare_graph(root, dimension, first, count, native, inputs, length)
        total = profile["nonaffine_graph_lifts"]
    else:
        profile = prepare_contraction(root, dimension, multiplicity, first, count, native, inputs, length)
        total = profile["full_rank_marked_lifts"]
    (work / "input.json").write_text(json.dumps(profile, indent=2)+"\n", encoding="ascii")
    require(total <= maximum_candidates,
            "The exact lift profile exceeds the candidate budget; subdivide the source interval")
    evidence = Evidence(root)
    for path, digest in profile["dependencies"].items():
        evidence.digest(path, digest)
    evidence.digest("code/space_native/generic/configuration.hpp")
    evidence.digest(f"code/space_native/generic/extensions_m{dimension:02d}.cpp")
    masks, ledger = work / "candidates.bin", work / "ledger.bin"
    command = [str(native / f"extensions_c{length}_m{dimension:02d}")]
    if dimension == 9:
        command += ["--input", str(inputs / "sources.bin"), "--output", str(masks),
                    "--source-start", "0", "--source-count", str(count),
                    "--expected-multiplicity", str(multiplicity)]
    elif dimension == 10:
        command += ["--source-input", str(inputs / "parents.bin"), "--profile-input", str(inputs / "profiles.bin"),
                    "--output", str(ledger), "--source-start", "0", "--profile-start", "0",
                    "--source-count", str(count), "--multiplicity", str(multiplicity), "--threads", str(threads)]
    else:
        command += ["--input", str(inputs / "sources.bin"), "--masks-output", str(work / "native_masks.bin"),
                    "--output", str(ledger), "--source-start", "0", "--source-count", str(count),
                    "--threads", str(threads)]
    completed = subprocess.run(command, text=True, capture_output=True, check=True)
    report = json.loads(completed.stdout)
    width = (1 << dimension)//8
    if dimension == 9:
        require(report["status"].startswith("complete_") and report["sources_completed"] == count
                and report["candidate_limit"] == 0 and report["processed_candidate_count"] == total,
                "The native lift traversal stopped before exhausting its domain")
        retained = report["retained_occurrence_count"]
        candidates = report["retained_distinct_batch_mask_count"]
        require(masks.stat().st_size == candidates*width, "The candidate bitmap stream is incomplete")
        if candidates:
            evidence.digest("code/space_native/generic/signature.cpp")
            subprocess.run([str(native / f"signature_c{length}_m{dimension:02d}"),
                "--input", str(masks), "--output", str(ledger), "--threads", str(threads)],
                text=True, capture_output=True, check=True)
        else:
            ledger.write_bytes(b"")
        if length == 50:
            require(retained == total, "The unfiltered length-50 traversal lost a lift")
    else:
        require(report["status"].startswith("complete_") and report["source_count"] == count
                and report["full_rank_occurrence_count"] == total, "The native lift count is incomplete")
        retained, candidates = total, report["distinct_candidate_count"]
        with ledger.open("rb") as incoming, masks.open("xb") as outgoing:
            while raw := incoming.read(32+width):
                require(len(raw) == 32+width, "Truncated signature record")
                outgoing.write(raw[32:])
    require(0 <= candidates <= retained <= total
            and report["within_batch_duplicate_count"] == retained-candidates
            and ledger.stat().st_size == candidates*(32+width), "The finite output accounting does not close")
    result = dict(schema="accelerated-space-interval-v1", status="pass", c=length, m=dimension,
        multiplicity=multiplicity, source_interval=[first, first+count], source_count=count,
        full_rank_marked_lifts=total, retained_occurrences=retained, candidates=candidates,
        minimum_direction_filter=length == 52 and dimension == 9,
        complete_interval_enumerated=True, candidate_sha256=file_digest(masks), ledger_sha256=file_digest(ledger),
        input_profile=profile, dependencies=evidence.bindings, is_global_completeness_certificate=False)
    (work / "generation.json").write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--length", type=int, required=True)
    parser.add_argument("--dimension", type=int, required=True)
    parser.add_argument("--multiplicity", type=int, default=0)
    parser.add_argument("--first", type=int, default=0)
    parser.add_argument("--count", type=int, required=True)
    parser.add_argument("--maximum-candidates", type=int, default=1000000)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--native-directory", type=Path, required=True)
    parser.add_argument("--work-directory", type=Path, required=True)
    args = parser.parse_args()
    result = generate(args.root, args.length, args.dimension, args.multiplicity, args.first, args.count,
                      args.native_directory, args.work_directory, args.maximum_candidates, args.threads)
    print(json.dumps({k: v for k, v in result.items() if k not in ("dependencies", "input_profile")}, indent=2))


if __name__ == "__main__":
    main()
