"""Verify the finite source domain and logical census through length 48."""

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path

from logical_spaces import logical_label_space
from protocol_checks import LogicalTensor, check_matrix, distance_up_to, find_equivalence_basis
from protocol_domain import (open_domain, read_header, read_pointings, reconstruct,
                             source_key, source_points)
from selected_spaces import read_selected
from verify_protocol_case_census import normalized_keys
from verify_protocols import verify as verify_frontier
from verify_tensor_authorities import restriction, subspace_bases
from verify_finite_protocol_sources import source_keys as finite_source_keys, verify as verify_finite_sources


def expected_source_keys(root):
    data = json.loads((root / "data/spaces/index.json").read_bytes())
    return {(row["c"], sector["m"], i) for row in data["lengths"] if row["c"] <= 48
            for sector in row["sectors"] for i in range(sector["expected_classes"])}


def nondominated_metrics(metrics):
    metrics = set(metrics)
    return {row for row in metrics if not any(
        other[:2] == row[:2] and other[2] <= row[2] and other[3] <= row[3] and other != row
        for other in metrics)}


def verify(root, full_domain=False):
    root = Path(root).resolve()
    source_path = root / "data/protocol_sectors/through48_source_domain.json"
    census_path = root / "data/protocol_sectors/through48_logical_censuses.json"
    source_bytes, census_bytes = source_path.read_bytes(), census_path.read_bytes()
    domain, census = json.loads(source_bytes), json.loads(census_bytes)
    if (domain["schema"] != "finite-protocol-source-domain-v1"
            or census["schema"] != "finite-through48-protocol-censuses-v1"
            or domain["maximum_protocol_length"] != 48 or census["maximum_protocol_length"] != 48):
        raise ValueError("Incorrect through-48 scope")
    if (census["minimum_distance"] != 3 or census["matrix_scope"] != "full_projective"
            or census["output_equivalence"] != "CNOT+S"):
        raise ValueError("Incorrect protocol or output equivalence scope")
    sources = domain["sources"]
    if [r["index"] for r in sources] != list(range(len(sources))):
        raise ValueError("Incomplete source identities")
    keys = [source_key(row) for row in sources]
    separate = {(row["c"], row["m"], row["index"]) for row in domain["separate_sources"]}
    finite_path = root / "data/protocol_sectors/finite_shorter_source_blocks.json"
    finite_data = json.loads(finite_path.read_bytes())
    finite_keys = finite_source_keys(finite_data)
    finite_closure = verify_finite_sources(root)
    if (len(keys) != len(set(keys)) or set(keys) & separate
            or (set(keys) | separate) & finite_keys
            or set(keys) | separate | finite_keys != expected_source_keys(root)):
        raise ValueError("Source domains do not partition the complete catalogue")
    supports, bindings = read_selected(root, keys)
    bindings.update(finite_closure["source_shard_bindings"])
    for source in sources:
        source_points(source, supports)
    selected = census["selected_quotients"]
    identities, dimensions = {}, Counter()
    for row in selected:
        index = row["pointing_index"]
        record = row["source_index"], row["parity_case"], row["origin"]
        points, h = reconstruct(sources, supports, record)
        if index in identities or not 0 <= index < domain["pointing_count"]:
            raise ValueError("Overlapping selected quotient inputs")
        if points != row["points"] or h != row["ambient_dimension"]:
            raise ValueError("Incorrect selected input coordinates")
        dimension = logical_label_space(points, h, 3)["quotient_dimension"]
        if dimension != row["quotient_dimension"] or dimension < 12:
            raise ValueError("Incorrect independently recomputed quotient dimension")
        identities[index] = record
        dimensions[dimension] += 1
    if (sorted(dimensions.items()) != [tuple(row) for row in census["quotient_selected_counts"]]
            or census["small_quotient_interval"] != [0, 11]
            or len(selected)+census["small_quotient_supports"] != domain["pointing_count"]
            or census["pointing_count"] != domain["pointing_count"]):
        raise ValueError("Quotient partitions do not cover the finite domain")
    path = (root / domain["pointing_data"]).resolve()
    if not path.is_relative_to(root):
        raise ValueError("External finite domain")
    with path.open("rb") as raw:
        digest = hashlib.file_digest(raw, "sha256").hexdigest()
    if digest != domain["pointing_sha256"]:
        raise ValueError("Changed compact pointing domain")
    counts, parities, lengths = Counter(), Counter(), Counter()
    with open_domain(path) as stream:
        header = read_header(stream)
        if (header["source_count"], header["count"], header["maximum_protocol_length"]) != (
                len(sources), domain["pointing_count"], 48):
            raise ValueError("Wrong pointing header")
        for index, record in enumerate(read_pointings(stream, header)):
            source, parity, origin = record
            if index in identities and record != identities[index]:
                raise ValueError("Selected input has wrong global record identity")
            counts[source] += 1
            parities[parity] += 1
            if full_domain:
                points, _ = reconstruct(sources, supports, record)
                if not 0 < len(points) <= 48:
                    raise ValueError("Pointing outside protocol-length scope")
                lengths[len(points)] += 1
    if ([counts[i] for i in range(len(sources))] != domain["source_pointing_counts"]
            or dict(parities) != {int(k): v for k, v in domain["parity_counts"].items()}
            or full_domain and dict(lengths) != {int(k): v for k, v in domain["protocol_length_counts"].items()}):
        raise ValueError("Finite pointing statistics changed")

    frontier_path = root / "data/protocols/pareto_frontier.json"
    verify_frontier(frontier_path)
    frontier = json.loads(frontier_path.read_bytes())
    outputs = {(r["q"], tuple(int(word, 0) for word in r["tensor_key_words"])): r["output_id"]
               for r in frontier["outputs"]}
    tensors = {r["output_id"]: LogicalTensor(tuple(int(x, 2) for x in r["generator_matrix_rows"][:r["q"]]))
               for r in frontier["protocols"]}
    covered, q7_tensors, witness_checks = defaultdict(set), {}, []
    finite_metrics = set()
    for index, branch in enumerate(census["branches"]):
        q = branch["q"]
        if type(q) is not int or not 1 <= q <= 8:
            raise ValueError("Unexpected logical census dimension")
        low, high = branch["quotient_dimension_interval"]
        ds = {d for d in [0, 11, *dimensions] if low <= d <= high}
        if covered[q] & ds or branch["minimum_distance"] != 3:
            raise ValueError("Overlapping logical census intervals")
        covered[q] |= ds
        count = sum(n for d, n in dimensions.items() if low <= d <= high)
        if low == 0:
            count += census["small_quotient_supports"]
        stats = branch["statistics"]
        if (count != branch["selected_supports"]
                or stats["supports_processed"]-stats["quotient_dimension_filtered_supports"] != count
                or stats["raw_enumerated_supports"]+stats["marked_orbit_enumerated_supports"] != stats["eligible_supports"]):
            raise ValueError("Incomplete logical census domain")
        mode = branch["enumeration_mode"]
        if mode not in ("raw_isotropic_subspaces", "complete_marked_code_orbits", "hybrid_raw_and_marked_code_orbits"):
            raise ValueError("Unsupported complete enumeration family")
        if mode == "raw_isotropic_subspaces":
            if (stats["marked_orbit_enumerated_supports"] or not stats["isotropic_subspace_statistics_exact"]
                    or not stats["radical_statistics_exact"]
                    or sum(stats["radical_dimension_counts"].values())+stats["nondegenerate_subspaces"] != stats["isotropic_subspaces"]):
                raise ValueError("Inconsistent exact raw census")
        elif stats["isotropic_subspaces"] is not None or stats["isotropic_subspace_statistics_exact"]:
            raise ValueError("Orbit census is incorrectly labelled an exact raw count")
        if q == 8 and (stats["nondegenerate_subspaces"] or branch["witnesses"]):
            raise ValueError("The dimension-eight census is not empty")
        if q == 7 and high <= 14 and (stats["nondegenerate_subspaces"] or branch["witnesses"]):
            raise ValueError("The smaller-quotient dimension-seven census is not empty")
        observed = set()
        for position, witness in enumerate(branch["witnesses"]):
            masks, columns = check_matrix(witness["generator_matrix_rows"], q)
            tensor = LogicalTensor(masks[:q])
            key, = normalized_keys([witness["output_key_words"]])
            identifier = outputs.get((q, key))
            if identifier is None or not tensor.intrinsic():
                raise ValueError("Unknown or degenerate witness output")
            basis = find_equivalence_basis(tensor, tensors[identifier])
            distance, coefficient = distance_up_to(columns, q)
            if (basis is None or distance is None or distance < 3
                    or (witness["q"], witness["n"], witness["S"], witness["d_Z"], witness["error_coefficient"]) != (
                        q, len(columns), len(masks), distance, coefficient)):
                raise ValueError("Incorrect finite witness parameters or output")
            dominators = [r for r in frontier["protocols"] if r["output_id"] == identifier
                          and r["d_Z"] == distance and r["n"] <= witness["n"] and r["S"] <= witness["S"]]
            if not dominators:
                raise ValueError("A finite witness is not covered by the final frontier")
            observed.add(key)
            finite_metrics.add((identifier, distance, witness["n"], witness["S"]))
            if q == 7:
                q7_tensors[key] = tensor
            witness_checks.append(dict(branch=index, witness=position,
                dominator=min(r["index"] for r in dominators), output_basis=list(basis)))
        if q == 7 and len(observed) != stats["canonical_output_orbits"]:
            raise ValueError("Not all seven-dimensional output types are explicitly witnessed")
    for q in range(1, 8):
        if covered[q] != {0, 11, *dimensions}:
            raise ValueError("Incomplete logical dimension or quotient partition")
    if covered[8] != {d for d in dimensions if d >= 15}:
        raise ValueError("Incomplete dimension-eight successor domain")
    merged = nondominated_metrics(finite_metrics)
    released = {(r["output_id"], r["d_Z"], r["n"], r["S"]) for r in frontier["protocols"] if r["n"] <= 48}
    if merged != released:
        raise ValueError("The finite censuses do not merge to exactly the released through-48 frontier")
    hyperplanes = []
    for key, tensor in q7_tensors.items():
        found = [list(basis) for basis in subspace_bases(7, 6) if restriction(tensor, basis).intrinsic()]
        if not found:
            raise ValueError("A primitive seven-dimensional output invalidates higher-q closure")
        hyperplanes.append(dict(key_words=[f"0x{x:016x}" for x in key], basis=found[0], count=len(found)))
    return dict(schema="finite-through48-sector-verification-v1", status="pass",
        source_domain_sha256=hashlib.sha256(source_bytes).hexdigest(),
        logical_census_sha256=hashlib.sha256(census_bytes).hexdigest(),
        frontier_sha256=hashlib.sha256(frontier_path.read_bytes()).hexdigest(),
        main_source_spaces=len(sources), separate_source_spaces=len(separate),
        finite_block_source_spaces=len(finite_keys), finite_block_pointing_supports=finite_closure["pointing_supports"],
        finite_source_census_sha256=finite_closure["data_sha256"],
        finite_block_witnesses_strictly_dominated=True,
        pointing_supports=domain["pointing_count"], selected_quotients_recomputed=len(selected),
        logical_census_partitions=len(census["branches"]), witnesses_recomputed=len(witness_checks),
        merged_frontier_points=len(merged), exact_frontier_merge_agrees=True,
        all_pointings_reconstructed=full_domain,
        primitive_q7_and_all_q8_excluded_given_complete_finite_censuses=True,
        q7_nondegenerate_hyperplanes=hyperplanes, source_shard_bindings=bindings,
        witness_checks=witness_checks, origin_orbit_coverage_freshly_recomputed=False,
        logical_censuses_freshly_recomputed=False, separate_source_closure_is_additional_obligation=True,
        is_global_completeness_certificate=False)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=root)
    parser.add_argument("--full-domain", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.root, args.full_domain)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps({k: v for k, v in result.items() if k not in ("source_shard_bindings", "witness_checks")}, indent=2))


if __name__ == "__main__":
    main()
