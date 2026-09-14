"""Certify graph-lift sectors whose target classes have distinct affine invariants."""

import argparse
import hashlib
import itertools
import json
from pathlib import Path
import struct
import subprocess
from time import perf_counter

from affine_codeword_graph import canonicalise_many
from space_codec import binary_rank, read_exact, validate_unital_support
from space_lifts import evaluation_rows, reconstruct_lift, reduce_vector
from verify_protocol_cover import Evidence, require
from verify_zero_fibre_spaces import graph_labels, hyperplane_invariant, invariant, prepare_support, sector


def point_list(mask):
    result = []
    while mask:
        bit = mask & -mask
        result.append(bit.bit_length()-1)
        mask ^= bit
    return tuple(result)


def affine_pivots(points, dimension):
    pivots = {}
    for row in evaluation_rows(points, dimension, 1):
        row = reduce_vector(row, pivots)
        if row:
            pivots[row.bit_length()-1] = row
    require(len(pivots) == dimension+1, "The parent is not full rank")
    return pivots


def check_map(source, target, frame, dimension):
    require(len(frame) == dimension+1 and all(type(x) is int and 0 <= x < 1 << dimension for x in frame)
            and binary_rank(frame[1:]) == dimension, "Invalid affine map domain or rank")
    image = []
    for point in source:
        value = frame[0]
        for bit, column in enumerate(frame[1:]):
            if point >> bit & 1:
                value ^= column
        image.append(value)
    require(sorted(image) == list(target), "The affine map does not reproduce its target")


def graph_domain(parents, dimension):
    profiles, pivots, outstanding = [], [], []
    for i, points in enumerate(parents):
        rank, basis = graph_labels(points, dimension-1)
        affine = affine_pivots(points, dimension-1)
        labels, cosets = 0, set()
        for ordinal in range(1, 1 << len(basis)):
            labels ^= basis[(ordinal & -ordinal).bit_length()-1]
            cosets.add(reduce_vector(labels, affine))
        require(0 not in cosets and len(cosets) == (1 << len(basis))-1, "The graph complement repeats an affine coset")
        profiles.append(dict(source_index=i, quadratic_rank=rank, lift_dimension=len(basis), nonaffine_lifts=len(cosets)))
        pivots.append(affine)
        outstanding.append(cosets)
    return profiles, pivots, outstanding


def separated_targets(targets, dimension):
    keys = []
    for target in targets:
        validate_unital_support(target, dimension)
        prepared = prepare_support(target, dimension)
        keys.append((invariant(prepared), hyperplane_invariant(prepared)))
    require(len(set(keys)) == len(targets), "This sector requires an additional exact negative-decision audit")


def replay(root, data):
    """Reproduce the finite proof from compact parents and affine maps only."""
    started = perf_counter()
    dimension = data["m"]
    require(data["schema"] == "complete-separated-graph-space-cover-v1" and data["status"] == "pass"
            and data["c"] == 54 and dimension in (11, 12, 13, 14), "Invalid graph-cover certificate scope")
    evidence = Evidence(root)
    for name, digest in data["dependencies"].items():
        evidence.digest(name, digest)
    index = evidence.read("data/spaces/index.json")
    parents = sector(evidence, index, 54, dimension-1)
    targets = sector(evidence, index, 54, dimension)
    separated_targets(targets, dimension)
    profiles, pivots, outstanding = graph_domain(parents, dimension)
    require(profiles == data["source_profile"] and len(parents) == data["source_classes"]
            and len(targets) == data["target_classes"], "The graph source profile is not exact")
    reached, count = set(), 0
    for orbit in data["orbits"]:
        target_index = orbit["target_index"]
        require(type(target_index) is int and 0 <= target_index < len(targets)
                and target_index not in reached, "Repeated or invalid target orbit")
        reached.add(target_index)
        representative = tuple(orbit["representative_points"])
        validate_unital_support(representative, dimension)
        check_map(representative, targets[target_index], orbit["representative_to_target"], dimension)
        require(orbit["lifts"], "A target orbit has no graph lift")
        for row in orbit["lifts"]:
            i, labels = row["source_index"], row["labels"]
            require(type(i) is int and 0 <= i < len(parents) and type(labels) is int and 0 <= labels < 1 << 54,
                    "Invalid graph-lift source or label")
            coset = reduce_vector(labels, pivots[i])
            require(coset in outstanding[i], "Repeated or inadmissible graph-lift coset")
            outstanding[i].remove(coset)
            child = reconstruct_lift(parents[i], (), labels, dimension-1)
            check_map(representative, child, row["representative_to_lift"], dimension)
            count += 1
    require(reached == set(range(len(targets))) and not any(outstanding)
            and count == data["nonaffine_graph_lifts"], "The graph-lift proof omits part of the finite domain")
    return dict(schema="separated-graph-witness-replay-v1", status="pass", c=54, m=dimension,
        source_classes=len(parents), target_classes=len(targets), nonaffine_graph_lifts=count,
        every_affine_shear_coset_independently_enumerated=True, every_positive_affine_map_replayed=True,
        targets_separated_by_independent_exact_invariants=True, complete_given_predecessor_catalogue=True,
        predecessor_completeness_is_separate=True, native_executable_required=False,
        dependencies=evidence.bindings, elapsed_seconds=perf_counter()-started)


def graph_orbits(dimension, masks, candidates, targets, target_masks, executable):
    forms = canonicalise_many([executable], [(dimension, point_list(mask)) for mask in masks])
    by_mask = dict(zip(masks, forms, strict=True))
    by_key, orbits = {}, []
    for mask, target_index in sorted(target_masks.items(), key=lambda row: row[1]):
        form = by_mask[mask]
        key = tuple(form["canonical_points"])
        require(key not in by_key, "The graph backend merges independently distinct targets")
        by_key[key] = target_index
        orbits.append(dict(target_index=target_index, representative_points=list(key),
                           representative_to_target=form["frame"], lifts=[]))
    for mask, (source_index, labels) in candidates.items():
        form = by_mask[mask]
        key = tuple(form["canonical_points"])
        require(key in by_key, "A graph lift has no catalogue class under exact codeword-incidence equivalence")
        orbits[by_key[key]]["lifts"].append(dict(source_index=source_index, labels=labels,
                                               representative_to_lift=form["frame"]))
    require(all(row["lifts"] for row in orbits) and len(orbits) == len(targets), "The graph quotient omits a target")
    return orbits


def cover_report(evidence, dimension, parents, targets, candidates, candidates_path, profiles, orbits, started, backend):
    return dict(schema="complete-separated-graph-space-cover-v1", status="pass", c=54, m=dimension,
        source_classes=len(parents), target_classes=len(targets), nonaffine_graph_lifts=len(candidates),
        complete_zero_fibre_cover_by_pair_averaging=True, affine_backend=backend,
        every_affine_shear_coset_independently_enumerated=True, every_positive_affine_map_replayed=True,
        targets_separated_by_independent_exact_invariants=True, complete_given_predecessor_catalogue=True,
        predecessor_completeness_is_separate=True, source_profile=profiles, orbits=orbits,
        candidate_mask_sha256=hashlib.sha256(candidates_path.read_bytes()).hexdigest(),
        dependencies=evidence.bindings, elapsed_seconds=perf_counter()-started)


def verify(root, dimension, candidates_path, native, work, graph_backend=None):
    started = perf_counter()
    require(dimension in (11, 12, 13, 14), "Unsupported fixed-width affine kernel")
    evidence = Evidence(root)
    index = evidence.read("data/spaces/index.json")
    parents = sector(evidence, index, 54, dimension-1)
    targets = sector(evidence, index, 54, dimension)
    separated_targets(targets, dimension)
    parent_index = {points: i for i, points in enumerate(parents)}
    require(len(parent_index) == len(parents), "Repeated literal parent support")
    profiles, pivots, outstanding = graph_domain(parents, dimension)
    width = (1 << dimension)//8
    candidates = {}
    with candidates_path.open("rb") as stream:
        while raw := stream.read(width):
            require(len(raw) == width, "Truncated graph candidate")
            mask = int.from_bytes(raw, "little")
            points = point_list(mask)
            validate_unital_support(points, dimension)
            core = tuple(x >> 1 for x in points)
            require(core in parent_index, "A graph candidate has a different projection")
            i = parent_index[core]
            labels = sum((x & 1) << bit for bit, x in enumerate(points))
            coset = reduce_vector(labels, pivots[i])
            require(coset in outstanding[i], "Repeated or inadmissible graph-lift coset")
            outstanding[i].remove(coset)
            require(mask not in candidates, "Repeated literal graph candidate")
            candidates[mask] = (i, labels)
    require(not any(outstanding), "The graph-lift candidate stream omits affine cosets")
    target_masks = {sum(1 << x for x in points): i for i, points in enumerate(targets)}
    masks = sorted(set(candidates) | set(target_masks))
    if graph_backend is not None:
        orbits = graph_orbits(dimension, masks, candidates, targets, target_masks, graph_backend)
        evidence.digest("code/native/auxiliary/affine_codeword_graph.cpp")
        evidence.digest("code/native/vendor/bliss-0.77.zip")
        return cover_report(evidence, dimension, parents, targets, candidates, candidates_path, profiles, orbits, started,
                            "spanning-codeword-incidence-graph")
    work.mkdir(parents=True, exist_ok=True)
    paths = {name: work / (name+".bin") for name in ("input", "signatures", "buckets", "assignments", "classes", "negatives")}
    require(not any(path.exists() for path in paths.values()), "An affine quotient output already exists")
    with paths["input"].open("xb") as stream:
        for mask in masks:
            stream.write(mask.to_bytes(width, "little"))
    subprocess.run([str(native / "runtime_affine_signature"), "--input", str(paths["input"]),
                    "--output", str(paths["signatures"]), "--dimension", str(dimension), "--threads", "4"], check=True)
    records = []
    with paths["signatures"].open("rb") as stream:
        while raw := stream.read(width+32):
            require(len(raw) == width+32, "Truncated affine signature record")
            records.append((struct.unpack("<4Q", raw[:32]), int.from_bytes(raw[32:], "little")))
    require([r[1] for r in sorted(records, key=lambda r: r[1])] == masks, "The signature stream changes its input domain")
    require(records == sorted(records), "The signature ledger is not in exact bucket order")
    groups = [list(rows) for _, rows in itertools.groupby(records, key=lambda row: row[0])]
    require(len(groups) == len(targets), "The signature partition has an unassigned or unresolved bucket")
    offset = 0
    bucket_targets = []
    with paths["buckets"].open("xb") as stream:
        for group in groups:
            hit = [target_masks[mask] for _, mask in group if mask in target_masks]
            require(len(hit) == 1, "A graph bucket does not contain exactly one released target")
            bucket_targets.append(hit[0])
            stream.write(struct.pack("<6Q", *group[0][0], offset, len(group))
                         +group[0][1].to_bytes(width, "little")+group[-1][1].to_bytes(width, "little"))
            offset += len(group)
    subprocess.run([str(native / f"m{dimension:02d}_exact_quotient_kernel"),
                    "--ledger", str(paths["signatures"]), "--bucket-index", str(paths["buckets"]),
                    "--first-bucket", "0", "--last-bucket", str(len(groups)),
                    "--assignments", str(paths["assignments"]), "--classes", str(paths["classes"]),
                    "--negatives", str(paths["negatives"])], check=True)
    require(paths["negatives"].stat().st_size == 0, "A separated graph bucket contains inequivalent candidates")
    assignment_width = (4+2*(dimension+1)+3)//4*4
    assignments = paths["assignments"].open("rb")
    classes = paths["classes"].open("rb")
    orbits = []
    try:
        for bucket_index, (group, target_index) in enumerate(zip(groups, bucket_targets, strict=True)):
            row = read_exact(classes, width+24)
            bucket, local_class, representative_index = struct.unpack("<QII", row[:16])
            representative_mask = int.from_bytes(row[16:16+width], "little")
            members, = struct.unpack("<Q", row[-8:])
            require((bucket, local_class, representative_index, representative_mask, members)
                    == (bucket_index, 0, 0, group[0][1], len(group)), "The exact affine class record differs from its bucket")
            representative = point_list(representative_mask)
            orbit = dict(target_index=target_index, representative_points=list(representative), lifts=[])
            for _, mask in group:
                raw = read_exact(assignments, assignment_width)
                local, *frame = struct.unpack(f"<I{dimension+1}H", raw[:4+2*(dimension+1)])
                require(local == 0, "Invalid local affine assignment")
                points = point_list(mask)
                check_map(representative, points, frame, dimension)
                if mask in target_masks:
                    require(target_masks[mask] == target_index, "Wrong target representative in orbit")
                    orbit["representative_to_target"] = frame
                if mask in candidates:
                    source_index, labels = candidates[mask]
                    orbit["lifts"].append(dict(source_index=source_index, labels=labels, representative_to_lift=frame))
            require("representative_to_target" in orbit, "An orbit has no catalogue image")
            orbits.append(orbit)
        require(not classes.read(1) and not assignments.read(1), "Extraneous exact quotient records")
    finally:
        classes.close()
        assignments.close()
    for name in ("runtime_affine_signature", f"m{dimension:02d}_exact_quotient_kernel"):
        evidence.digest(f"code/space_native/length54/{name}.cpp")
    return cover_report(evidence, dimension, parents, targets, candidates, candidates_path, profiles, orbits, started,
                        "exact-basis-image-backtracking")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--dimension", type=int)
    parser.add_argument("--candidate-masks", type=Path)
    parser.add_argument("--native-directory", type=Path)
    parser.add_argument("--work-directory", type=Path)
    parser.add_argument("--graph-backend", type=Path, help="Optional spanning-codeword graph executable for affine maps")
    parser.add_argument("--certificate", type=Path, help="Replay a supplied finite witness proof without native code")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.certificate:
        raw = args.certificate.read_bytes()
        result = replay(args.root, json.loads(raw))
        result["certificate_sha256"] = hashlib.sha256(raw).hexdigest()
    else:
        require(all((args.dimension, args.candidate_masks, args.native_directory, args.work_directory)), "Missing enumeration arguments")
        result = verify(args.root, args.dimension, args.candidate_masks, args.native_directory.resolve(), args.work_directory.resolve(),
                        args.graph_backend.resolve() if args.graph_backend else None)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps({k: v for k, v in result.items() if k not in ("source_profile", "orbits", "dependencies")}, indent=2))


if __name__ == "__main__":
    main()
