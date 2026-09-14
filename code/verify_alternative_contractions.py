"""Construct or replay finite alternative-contraction covers in affine dimension 8."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import shutil
import struct
import subprocess
from time import perf_counter

from prepare_space_contractions import task_text
from selected_spaces import read_selected
from space_codec import binary_rank, validate_unital_support
from verify_protocol_cover import Evidence, require

CORES = ((48, 7, 96), (48, 7, 97), (52, 7, 186), (52, 7, 187))
DATA = "data/space_contractions/alternative_cover"
RECORD = struct.Struct("<5Q")
KERNELS = {50: "contractions/c50_m08_marked_orbits.cpp",
           52: "contractions/c52_m08_marked_orbits.cpp",
           54: "length54/m08_marked_orbit_kernel.cpp"}


def digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def contract(points, direction):
    """Odd fibres in the quotient by one nonzero vector."""
    require(type(direction) is int and 0 < direction < 256, "Invalid contraction direction")
    pivot = direction.bit_length()-1
    projected = set()
    for point in points:
        if point >> pivot & 1:
            point ^= direction
        value = (point & ((1 << pivot)-1)) | ((point >> (pivot+1)) << pivot)
        projected.symmetric_difference_update((value,))
    return tuple(sorted(projected))


def difference_signature(points):
    counts = Counter(x ^ y for x in points for y in points if x != y)
    histogram = tuple(sorted(Counter(counts[x] for x in range(1, 128)).items()))
    local = tuple(sorted(tuple(sorted(counts[x ^ y] for y in points if x != y)) for x in points))
    return len(points), histogram, local


def valid_cover(points, direction, excluded_signatures):
    core = contract(points, direction)
    if len(core) not in (44, 48, 52) or binary_rank(x ^ core[0] for x in core) != 7:
        return False
    validate_unital_support(core, 7)
    return difference_signature(core) not in excluded_signatures


def inputs(evidence):
    certificate = evidence.read("certificates/low_dimensional_spaces.json")
    require(certificate.get("status") == "pass" and certificate.get("complete_and_pairwise_inequivalent") is True
            and certificate.get("maximum_length") == 54 and certificate.get("maximum_affine_dimension") == 7,
            "The complete low-dimensional affine census is required")
    supports, bindings = read_selected(evidence.root, CORES)
    for name, value in bindings.items():
        evidence.digest(name, value)
        require(certificate["catalogue_sha256"].get(name) == value, "The affine census uses a different catalogue")
    forms = {(r["c"], r["m"], r["index"]): r for r in certificate["canonical_maps"]}
    excluded = ["M8_EXCLUDED_CORES_V1", "4"]
    for key in CORES:
        mask = sum(1 << x for x in supports[key])
        excluded.append(f"{mask & ((1 << 64)-1)} {mask >> 64}")
    cases = []
    for length in (50, 52, 54):
        for key in CORES:
            if key[0] > length:
                continue
            identifier = f"c{length}_from_c{key[0]}_m07_{key[2]:09d}"
            text, profile = task_text(identifier, supports[key], forms[key], length)
            cases.append(dict(id=identifier, target_length=length, core=list(key), text=text, profile=profile))
    return cases, "\n".join(excluded)+"\n", supports


def check_counts(summary, case):
    profile = case["profile"]
    require(summary.get("all_checks_pass") is True and summary.get("count_only") is False
            and summary["parent_id"] == case["id"]
            and summary["parent_weight"] == profile["core_length"]
            and summary["full_fiber_multiplicity"] == profile["full_fibres"]
            and summary["compatible_fiber_set_count"] == profile["compatible_fibre_sets"]
            and summary["lift_quotient_dimension"] == profile["lift_quotient_dimension"]
            and summary["core_stabilizer_order"] == profile["stabiliser_order"]
            and summary["raw_marked_pair_count"] == profile["raw_marked_pairs"]
            and summary["excluded_rank_deficient_pair_count"] == profile["rank_deficient_pairs"],
            "Marked enumeration does not account for the exact input domain")


def check_cover(summary, case, count):
    require(summary["status"] == "complete_alternative_contraction_cover"
            and summary["target_length"] == case["target_length"]
            and summary["candidate_count"] == summary["covered_count"] == count
            and summary["unresolved_count"] == 0
            and summary["raw_marked_mass"] == summary["covered_raw_marked_mass"] == case["profile"]["raw_marked_pairs"]
            and "0" not in summary["direction_counts"]
            and sum(summary["direction_counts"].values()) == count,
            "The alternative-contraction proof has uncovered cases")


def sample_proof(binary, proof, case, support, excluded_signatures, count):
    """Independent Python checks across each complete native proof stream."""
    require(binary.stat().st_size == count * RECORD.size and proof.stat().st_size == count,
            "The candidate and direction streams have inconsistent lengths")
    indices = sorted({i * (count-1) // 127 for i in range(128)})
    with binary.open("rb") as candidates, proof.open("rb") as directions:
        for index in indices:
            candidates.seek(index * RECORD.size)
            *words, mass = RECORD.unpack(candidates.read(RECORD.size))
            mask = sum(word << (64*i) for i, word in enumerate(words))
            points = tuple(x for x in range(256) if mask >> x & 1)
            directions.seek(index)
            direction, = directions.read(1)
            require(len(points) == case["target_length"] and mass > 0, "Invalid marked support")
            validate_unital_support(points, 8)
            require(contract(points, 1) == tuple(support), "The marked support has a different source")
            require(valid_cover(points, direction, excluded_signatures), "Invalid alternative-contraction direction")
    return len(indices)


def compress(source, target):
    import zstandard
    require(not target.exists(), "Refusing to replace compressed evidence")
    with source.open("rb") as incoming, target.open("xb") as outgoing:
        zstandard.ZstdCompressor(level=9, threads=2, write_checksum=True).copy_stream(incoming, outgoing)


def decompress(evidence, name, expected, destination):
    import zstandard
    evidence.digest(name, expected)
    with evidence.path(name).open("rb") as incoming, destination.open("xb") as outgoing:
        with zstandard.ZstdDecompressor().stream_reader(incoming) as decoded:
            shutil.copyfileobj(decoded, outgoing)


def verify(root, native, work, create=False):
    started = perf_counter()
    evidence = Evidence(root)
    native, work = Path(native).resolve(), Path(work).resolve()
    require(not work.exists(), "Use a new proof work directory")
    work.mkdir(parents=True)
    cases, excluded_text, supports = inputs(evidence)
    directory = evidence.path(DATA)
    directory.mkdir(parents=True, exist_ok=True)
    excluded_path = directory / "excluded_cores.txt"
    if create:
        require(not (directory / "domain.json").exists(), "The proof domain already exists")
        excluded_path.write_text(excluded_text, encoding="ascii")
        domain = None
    else:
        domain = evidence.read(DATA + "/domain.json")
        require(domain["schema"] == "alternative-contraction-domain-v1"
                and len(domain["cases"]) == len(cases), "Incomplete alternative-contraction domain")
    require(excluded_path.read_text() == excluded_text, "The excluded family differs from the complete catalogue")
    evidence.digest(DATA + "/excluded_cores.txt")
    signatures = {difference_signature(supports[key]) for key in CORES}
    cover_source = "code/space_native/contractions/alternative_contraction_cover.cpp"
    evidence.digest(cover_source)
    rows, samples = [], 0
    for case_index, case in enumerate(cases):
        identifier = case["id"]
        case_dir = work / identifier
        case_dir.mkdir()
        task_name = DATA + "/" + identifier + ".task"
        task = evidence.path(task_name)
        if create:
            task.write_text(case["text"], encoding="ascii")
        require(task.read_text() == case["text"], "The task does not reproduce from the complete affine census")
        task_sha = evidence.digest(task_name)
        binary, proof = case_dir / "candidates.bin", case_dir / "directions.bin"
        source = "code/space_native/" + KERNELS[case["target_length"]]
        evidence.digest(source)
        if create:
            summary = case_dir / "generation.json"
            subprocess.run([str(native / Path(source).stem), "--input", str(task), "--output", str(binary),
                            "--summary", str(summary)], check=True, capture_output=True)
            generation = json.loads(summary.read_bytes())
            check_counts(generation, case)
            count = generation["direction_marked_orbit_count"]
            summary = case_dir / "cover.json"
            residual = case_dir / "residual.bin"
            subprocess.run([str(native / "alternative_contraction_cover"), "--input", str(binary),
                            "--core-input", str(task), "--excluded", str(excluded_path),
                            "--length", str(case["target_length"]), "--proof", str(proof),
                            "--residual", str(residual), "--summary", str(summary)], check=True, capture_output=True)
            require(residual.stat().st_size == 0, "An alternative core has unresolved children")
            check_cover(json.loads(summary.read_bytes()), case, count)
        else:
            row = domain["cases"][case_index]
            require(row["id"] == identifier and row["core"] == case["core"]
                    and row["target_length"] == case["target_length"] and row["profile"] == case["profile"]
                    and row["task"] == task_name and row["task_sha256"] == task_sha,
                    "The proof manifest does not describe the reproduced domain")
            generation = row["generation"]
            check_counts(generation, case)
            count = generation["direction_marked_orbit_count"]
            for label, path in (("candidates", binary), ("directions", proof)):
                item = row[label]
                decompress(evidence, item["path"], item["sha256"], path)
                require(digest(path) == item["uncompressed_sha256"], "Changed uncompressed proof stream")
        summary = case_dir / "replay.json"
        subprocess.run([str(native / "alternative_contraction_cover"), "--input", str(binary),
                        "--core-input", str(task), "--excluded", str(excluded_path),
                        "--length", str(case["target_length"]), "--verify-proof", str(proof),
                        "--summary", str(summary)], check=True, capture_output=True)
        replay = json.loads(summary.read_bytes())
        check_cover(replay, case, count)
        require(replay["verification_mode"] is True, "Expected an independent proof replay")
        samples += sample_proof(binary, proof, case, supports[tuple(case["core"])], signatures, count)
        if create:
            row = {k: v for k, v in case.items() if k != "text"}
            row.update(task=task_name, task_sha256=task_sha,
                       generation={k: v for k, v in generation.items() if k != "elapsed_seconds"})
            for label, path in (("candidates", binary), ("directions", proof)):
                name = DATA + "/" + identifier + "." + label + ".bin.zst"
                compressed = evidence.path(name)
                compress(path, compressed)
                row[label] = dict(path=name, sha256=evidence.digest(name),
                                  uncompressed_sha256=digest(path), bytes=path.stat().st_size)
        rows.append(row)
        print(json.dumps(dict(case=identifier, marked_representatives=count, raw_pairs=case["profile"]["raw_marked_pairs"],
                              all_directions_replayed=True)), flush=True)
    if create:
        domain = dict(schema="alternative-contraction-domain-v1", affine_dimension=8,
                      excluded_cores=[list(key) for key in CORES], cases=rows)
        (directory / "domain.json").write_text(json.dumps(domain, indent=2)+"\n", encoding="ascii")
    evidence.digest(DATA + "/domain.json")
    totals = {str(length): dict(cases=sum(r["target_length"] == length for r in rows),
        marked_representatives=sum(r["generation"]["direction_marked_orbit_count"] for r in rows if r["target_length"] == length),
        raw_marked_pairs=sum(r["profile"]["raw_marked_pairs"] for r in rows if r["target_length"] == length)) for length in (50, 52, 54)}
    return dict(schema="alternative-contraction-cover-verification-v1", status="pass", totals=totals,
        marked_censuses_freshly_generated=create, every_direction_freshly_replayed=True,
        all_candidates_valid_and_contract_to_the_declared_source=True,
        python_reference_samples=samples, unresolved_cases=0, dependencies=evidence.bindings,
        conditional_premises=["The certified low-dimensional affine groups are complete.",
                              "The primary core-family child census is complete."],
        is_global_completeness_certificate=False, elapsed_seconds=round(perf_counter()-started, 6))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--native-directory", type=Path, required=True)
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument("--create", action="store_true", help="Generate the complete marked families and compressed witnesses")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.root, args.native_directory, args.work_directory, args.create)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps({k: v for k, v in result.items() if k != "dependencies"}, indent=2))


if __name__ == "__main__":
    main()
