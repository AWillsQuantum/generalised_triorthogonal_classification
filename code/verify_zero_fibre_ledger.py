"""Regenerate a complete graph-lift ledger and check its exact binary census binding."""

import argparse
from collections import Counter
import json
from pathlib import Path
import struct
import subprocess
from time import perf_counter

from prepare_zero_fibre_input import file_digest
from verify_protocol_cover import Evidence, require


SIGNATURE = struct.Struct("<4Q")
INTERVAL = struct.Struct("<2Q")


def build_buckets(ledger, output, dimension, length=54):
    """Independently check the sorted ledger and construct its signature index."""
    mask_bytes = (1 << dimension)//8
    size = SIGNATURE.size+mask_bytes
    current = prior_mask = None
    first = last = b""
    offset = count = records = buckets = 0
    histogram = Counter()
    with ledger.open("rb") as incoming, output.open("xb") as outgoing:
        def finish():
            outgoing.write(SIGNATURE.pack(*current)+INTERVAL.pack(offset, count)+first+last)
            histogram[count] += 1

        while raw := incoming.read(size):
            require(len(raw) == size, "Truncated signature record")
            signature = SIGNATURE.unpack(raw[:SIGNATURE.size])
            mask = raw[SIGNATURE.size:]
            value = int.from_bytes(mask, "little")
            require(value.bit_count() == length, "Candidate has a different support weight")
            if signature != current:
                if current is not None:
                    require(signature > current, "Signature ledger is not sorted")
                    finish()
                current, prior_mask = signature, None
                offset, count, first = records, 0, mask
                buckets += 1
            require(prior_mask is None or value > prior_mask, "Repeated or unsorted candidate in a signature bucket")
            prior_mask, last = value, mask
            count += 1
            records += 1
        if current is not None:
            finish()
    return dict(candidate_count=records, signature_bucket_count=buckets,
                bucket_size_distribution=dict(sorted(histogram.items())))


def verify(root, dimension, sources, native, work, threads, maximum_working_bytes):
    started = perf_counter()
    require(dimension in (11, 12) and threads >= 1, "Unsupported graph-ledger parameters")
    evidence = Evidence(root)
    census = evidence.read(f"data/space_contractions/middle_sectors/c54_m{dimension:02d}.json")
    profile = evidence.certificate(f"certificates/zero_fibre_native_c54_m{dimension:02d}_input.json",
        "native-zero-fibre-input-v1", ("complete_parent_sector_included", "bounded_memory_catalogue_stream"),
        dict(c=54, m=dimension))
    expected = census["generation"]
    binding = census["generation_binary_binding"]
    require(profile["source_count"] == census["source_domains"][0]["classes"]
            and profile["nonaffine_graph_lifts"] == expected["candidate_count"], "The complete source domain differs")
    source_bytes = profile["source_count"]*profile["source_record_bytes"]
    mask_bytes = (1 << dimension)//8
    estimated_bytes = source_bytes+profile["nonaffine_graph_lifts"]*(2*mask_bytes+32)
    require(estimated_bytes <= maximum_working_bytes, "The complete ledger exceeds the requested working-memory budget")
    require(sources.stat().st_size == source_bytes and file_digest(sources) == profile["source_sha256"],
            "The graph-lift source stream differs from its full-domain certificate")
    work.mkdir(parents=True, exist_ok=True)
    masks, ledger, buckets, summary = [work / name for name in
                                     ("masks.bin", "signatures.bin", "buckets.bin", "generation.json")]
    require(not any(p.exists() for p in (masks, ledger, buckets, summary)), "A generation output already exists")
    name = f"m{dimension:02d}_fused_signature_kernel"
    evidence.digest(f"code/space_native/length54/{name}.cpp", census["native_sources_sha256"][f"code/space_native/length54/{name}.cpp"])
    result = subprocess.run([str(native/name), "--input", str(sources), "--masks-output", str(masks),
        "--output", str(ledger), "--source-start", "0", "--source-count", str(profile["source_count"]),
        "--threads", str(threads)], check=True, capture_output=True, text=True)
    generated = json.loads(result.stdout)
    summary.write_text(json.dumps(generated, indent=2)+"\n", encoding="ascii")
    require(generated["source_count"] == profile["source_count"]
            and generated["full_rank_occurrence_count"] == expected["candidate_occurrence_count"]
            and generated["distinct_candidate_count"] == expected["candidate_count"]
            and generated["within_batch_duplicate_count"] == expected["within_batch_literal_duplicate_count"],
            "The regenerated candidate census differs")
    require(ledger.stat().st_size == binding["ledger_bytes"] and file_digest(ledger) == binding["ledger_sha256"],
            "The complete regenerated signature ledger differs from the retained digest")
    checked = build_buckets(ledger, buckets, dimension)
    require(checked["candidate_count"] == expected["candidate_count"]
            and checked["signature_bucket_count"] == expected["signature_bucket_count"]
            and {str(k): v for k, v in checked["bucket_size_distribution"].items()} == expected["bucket_size_distribution"],
            "The independently reconstructed bucket census differs")
    require(buckets.stat().st_size == binding["bucket_index_bytes"]
            and file_digest(buckets) == binding["bucket_index_sha256"], "The regenerated bucket index differs")
    return dict(schema="complete-zero-fibre-ledger-regeneration-v1", status="pass", c=54, m=dimension,
        source_count=profile["source_count"], candidate_count=checked["candidate_count"],
        signature_bucket_count=checked["signature_bucket_count"],
        complete_candidate_generation_freshly_recomputed=True, sorted_ledger_and_bucket_index_byte_identical=True,
        exact_affine_quotient_freshly_recomputed=False, estimated_peak_array_bytes=estimated_bytes,
        generation_seconds=generated["elapsed_seconds"], elapsed_seconds=perf_counter()-started,
        binary_binding=binding, dependencies=evidence.bindings)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--dimension", type=int, required=True)
    parser.add_argument("--sources", type=Path, required=True)
    parser.add_argument("--native-directory", type=Path, required=True)
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--maximum-working-bytes", type=int, default=2_000_000_000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.root, args.dimension, args.sources.resolve(), args.native_directory.resolve(),
                    args.work_directory.resolve(), args.threads, args.maximum_working_bytes)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps({k: v for k, v in result.items() if k != "dependencies"}, indent=2))


if __name__ == "__main__":
    main()
