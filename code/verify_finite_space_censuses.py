"""Check finite contraction proofs, optionally regenerating every affine witness input."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

from affine_codeword_graph import canonicalise_many
from finite_space_census import generate_domain, replay_quotient
from prepare_zero_fibre_input import file_digest
from selected_spaces import read_selected
from space_codec import binary_rank, checked_support, validate_unital_support
from space_recursion import multiplicity_bound
from verify_protocol_cover import Evidence, require
from verify_zero_fibre_ledger import build_buckets


def check_frame(points, canonical, frame, m):
    checked_support(points, m)
    checked_support(canonical, m, len(points))
    require(len(frame) == m + 1 and binary_rank(frame[1:]) == m
            and all(type(x) is int and 0 <= x < 1 << m for x in frame), "Invalid affine frame")
    images = []
    for point in canonical:
        value = frame[0]
        for bit, column in enumerate(frame[1:]):
            if point >> bit & 1:
                value ^= column
        images.append(value)
    require(sorted(images) == list(points), "An affine frame does not reproduce its support")


def unpack_stream(evidence, row, output=None):
    import zstandard
    path = evidence.path(row["path"])
    evidence.digest(row["path"], row["sha256"])
    require(path.stat().st_size == row["compressed_bytes"], "Incorrect compressed witness size")
    digest, size = hashlib.sha256(), 0
    with path.open("rb") as raw, zstandard.ZstdDecompressor().stream_reader(raw) as stream:
        target = output.open("wb") if output is not None else None
        try:
            while block := stream.read(1 << 20):
                digest.update(block)
                size += len(block)
                if target is not None:
                    target.write(block)
        finally:
            if target is not None:
                target.close()
    require(size == row["decoded_bytes"] and digest.hexdigest() == row["decoded_sha256"],
            "The lossless affine witness stream differs")


def check_sector_frames(root, data):
    c, m, count = data["c"], data["m"], data["class_count"]
    require(len(data["class_supports"]) == len(data["canonical_forms"]) == len(data["catalogue_frames"]) == count,
            "Incomplete finite class records")
    require(sorted(data["catalogue_indices"]) == list(range(count)), "The catalogue matching is not a permutation")
    targets, bindings = read_selected(root, {(c, m, i) for i in range(count)})
    keys = []
    for points, form, target_index, frame in zip(data["class_supports"], data["canonical_forms"],
                                               data["catalogue_indices"], data["catalogue_frames"], strict=True):
        validate_unital_support(points, m)
        key = tuple(form["canonical_points"])
        require((form["m"], form["n"]) == (m, c), "Wrong canonical-form parameters")
        check_frame(points, key, form["frame"], m)
        check_frame(targets[c, m, target_index], key, frame, m)
        keys.append(key)
    require(len(set(keys)) == count, "Repeated claimed affine canonical form")
    return bindings, keys


def replay_sector(root, data, evidence, native, graph, work):
    work.mkdir(parents=True, exist_ok=True)
    c, m = data["c"], data["m"]
    masks, ledger, buckets = [work / name for name in ("candidates.bin", "signatures.bin", "generated_buckets.bin")]
    sources = [evidence.path(row["support_file"]) for row in data["source_domains"] if "support_file" in row]
    domain = generate_domain(root, c, m, sources, masks, data["candidate_count"] + 1)
    require(domain["candidates"] == data["candidate_count"]
            and domain["candidate_sha256"] == data["candidate_bitmap_sha256"]
            and domain["source_profiles"] == data["source_profiles"], "The complete finite candidate domain differs")
    require(len(domain["source_domains"]) == len(data["source_domains"]), "A contraction core family is missing")
    for actual, expected in zip(domain["source_domains"], data["source_domains"], strict=True):
        require(all(actual[k] == expected[k] for k in ("c", "m", "classes", "multiplicity")), "A contraction family differs")
    if data["candidate_count"]:
        process = subprocess.run([str(native / f"signature_c{c}_m{m:02d}"), "--input", str(masks),
            "--output", str(ledger), "--threads", "4"], capture_output=True, text=True, check=True)
        (work / "signature.json").write_text(process.stdout, encoding="ascii")
        require(file_digest(ledger) == data["sorted_ledger_sha256"], "The complete sorted candidate ledger differs")
        build_buckets(ledger, buckets, m, c)
        require(file_digest(buckets) == data["affine_witness_streams"]["buckets"]["decoded_sha256"],
                "The independently generated signature partition differs")
        files = {}
        for name, row in data["affine_witness_streams"].items():
            files[name] = work / (name + ".bin")
            unpack_stream(evidence, row, files[name])
        replayed = replay_quotient(ledger, files["buckets"], files["classes"], files["assignments"], c, m)
        require(replayed["candidate_count"] == data["candidate_count"]
                and replayed["class_supports"] == [tuple(p) for p in data["class_supports"]]
                and replayed["member_counts"] == data["class_member_counts"], "Positive affine witness replay differs")
    else:
        require(data["class_count"] == 0, "Empty candidate domain has an output class")
    records = [(m, p) for p in data["class_supports"]]
    forms = canonicalise_many([graph], records) if records else []
    require([row["canonical_points"] for row in forms] == [row["canonical_points"] for row in data["canonical_forms"]],
            "Independent exact canonical separation differs")
    return dict(c=c, m=m, candidate_count=data["candidate_count"], class_count=data["class_count"],
                every_candidate_and_positive_affine_map_replayed=True, every_class_independently_separated=True)


def verify(root, native=None, graph=None, work=None):
    evidence = Evidence(root)
    root = evidence.root
    index = evidence.read("data/space_contractions/finite_sectors/index.json")
    catalogue = evidence.read("data/spaces/index.json")
    counts = {(row["c"], sector["m"]): sector["expected_classes"] for row in catalogue["lengths"] for sector in row["sectors"]}
    expected = {(c, m) for c in range(16, 47, 2) for m in range(8, (c+c//16)//3)}
    expected.update((48, m) for m in range(9, 16))
    require({(r["c"], r["m"]) for r in index["sectors"]} == expected
            and len(index["sectors"]) == len(expected), "The finite sector cover has a gap or overlap")
    seen, reports = {}, []
    total_candidates = total_classes = 0
    for entry in index["sectors"]:
        data = evidence.read(entry["path"], entry["sha256"])
        c, m = data["c"], data["m"]
        require(data["schema"] == "finite-affine-contraction-sector-v1"
                and data["class_count"] == counts[c, m] == entry["class_count"]
                and data["candidate_count"] == entry["candidate_count"], "Wrong finite sector cardinality")
        for source in data["source_domains"]:
            key = source["c"], source["m"]
            require(key[0] <= c and key[1] < m and source["classes"] == counts.get(key, 0),
                    "Invalid contraction predecessor or missing class")
            if "support_file" in source:
                require(key in seen and source["support_file"] == seen[key]["path"]
                        and source["support_file_sha256"] == seen[key]["sha256"], "Unproved or misbound supplied predecessor")
            else:
                require(source["support_source"] == "compact_catalogue"
                        and (key[1] <= 7 or key == (48, 8) or key in seen),
                        "An internal induction predecessor is unproved")
        bindings, _ = check_sector_frames(root, data)
        for path, sha in bindings.items():
            evidence.digest(path, sha)
        for row in data["affine_witness_streams"].values():
            unpack_stream(evidence, row)
        require(sum(data["class_member_counts"]) == data["candidate_count"]
                and len(data["class_member_counts"]) == data["class_count"], "Incomplete affine class mass")
        if native is not None:
            require(graph is not None and work is not None, "A complete replay needs native kernels, graph canonicaliser and work directory")
            report = replay_sector(root, data, evidence, Path(native).resolve(), Path(graph).resolve(),
                                   Path(work).resolve() / f"c{c}_m{m:02d}")
        else:
            report = dict(c=c, m=m, candidate_count=data["candidate_count"], class_count=data["class_count"])
        reports.append(report)
        seen[c, m] = entry
        total_candidates += data["candidate_count"]
        total_classes += data["class_count"]
        print(json.dumps(report), flush=True)
    require(multiplicity_bound(48, 16) == 0 and counts[48, 15] == counts[48, 16] == 0,
            "The upper empty zero-fibre successor remains open")
    return dict(schema="finite-affine-contraction-collection-verification-v1", status="pass",
        sectors=len(reports), candidates=total_candidates, classes=total_classes, reports=reports,
        complete_through_length=46, higher_length_conditional_sectors=dict(c=48, affine_dimensions=list(range(9, 17))),
        all_catalogue_records_matched_by_affine_frames=True, all_witness_streams_losslessly_verified=True,
        all_candidates_and_positive_maps_freshly_replayed=native is not None,
        all_class_canonical_forms_freshly_recomputed=native is not None,
        conditional_premises=["The complete RM(3,7) orbit-mass base through affine dimension seven",
                              "The separate complete length-48 affine-dimension-eight marked-contraction sector"],
        dependencies=evidence.bindings, is_global_completeness_certificate=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--native-directory", type=Path)
    parser.add_argument("--graph-native", type=Path)
    parser.add_argument("--work-directory", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.root, args.native_directory, args.graph_native, args.work_directory)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="ascii")
    print(json.dumps({k: v for k, v in result.items() if k not in ("reports", "dependencies")}, indent=2))


if __name__ == "__main__":
    main()
