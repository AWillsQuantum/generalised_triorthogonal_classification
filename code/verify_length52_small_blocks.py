"""Verify the length-52 source and small-quotient finite census partition."""

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path

from protocol_checks import LogicalTensor, check_matrix, distance_up_to, find_equivalence_basis
from verify_high_dimensional_sector import check_intervals
from verify_low_dimensional_blocks import histogram
from verify_low_q_pruning import verify as verify_low_q
from verify_protocol_case_census import normalized_keys
from verify_protocols import verify as verify_frontier


def check_raw_pointings(m, sources, points, parities):
    half = ((1 << m)+52)//2
    difference = sources*(2*half+1)-points
    if difference < 0 or difference % half or difference//half > sources:
        raise ValueError("Incorrect raw pointing count for length 52")
    periodic = difference//half
    expected = [(1 << m)*sources-(1 << (m-1))*periodic, 52*sources-26*periodic, 0, sources, 0]
    if parities != expected:
        raise ValueError("Pointing parity counts disagree with translation periods")
    return periodic


def check_branch(branch, profile, sources):
    q, floor = branch["q"], branch["minimum_distance"]
    if q not in range(1, 6) or floor != (3 if q == 5 else 4):
        raise ValueError("Wrong logical dimension or distance floor")
    minimum = 5 if q == 5 else 1
    selected = {d: n for d, n in profile.items() if d >= minimum}
    size = sum(selected.values())
    interval = branch["support_range"]
    start, count, end = interval["start"], interval["count"], interval["end_exclusive"]
    if any(type(v) is not int for v in (start, count, end)) or not 0 <= start <= end <= size or end != start+count:
        raise ValueError("Incorrect finite input interval")
    local = histogram(branch["input_quotient_profile"], count)
    if any(d < minimum or n > selected.get(d, 0) for d, n in local.items()):
        raise ValueError("Interval quotient distribution is outside its source domain")
    if start == 0 and count == size and local != selected:
        raise ValueError("Full interval quotient distribution is incomplete")
    low, high = branch["quotient_dimension_interval"]
    if type(low) is not int or type(high) is not int or not 0 <= low <= high <= 64:
        raise ValueError("Incorrect quotient stratum")
    eligible = sum(n for d, n in local.items() if max(q, low) <= d <= high)
    filtered = sum(n for d, n in local.items() if not low <= d <= high)
    stats = branch["statistics"]
    if branch["enumeration_mode"] == "alternating_tensor_exclusion":
        if q not in (1, 2, 4) or floor != 4 or stats is not None or branch["witnesses"] or branch["positive_supports"]:
            raise ValueError("Incorrect scope for the alternating tensor exclusion")
        return local
    if branch["enumeration_mode"] == "empty_by_quotient_dimension":
        if stats is not None or eligible or branch["witnesses"] or branch["positive_supports"]:
            raise ValueError("A nonempty finite domain is labelled empty")
        return local
    if branch["enumeration_mode"] not in ("raw_isotropic_subspaces", "hybrid_raw_and_marked_code_orbits"):
        raise ValueError("Unsupported finite census mode")
    if (stats["source_spaces"] != sources or stats["supports_processed"] != count
            or stats["eligible_supports"] != eligible or stats["quotient_dimension_filtered_supports"] != filtered
            or stats["raw_enumerated_supports"]+stats["marked_orbit_enumerated_supports"] != eligible):
        raise ValueError("Incomplete finite-subspace traversal")
    if stats["marked_orbit_enumerated_supports"]:
        if (q == 5 or branch["enumeration_mode"] != "hybrid_raw_and_marked_code_orbits"
                or stats["isotropic_subspace_statistics_exact"] is not False
                or stats["radical_statistics_exact"] is not False or stats["isotropic_subspaces"] is not None):
            raise ValueError("An orbit census is incorrectly presented as exact raw counts")
    elif (stats["marked_code_orbits"] or stats["isotropic_subspace_statistics_exact"] is not True
            or stats["radical_statistics_exact"] is not True
            or sum(stats["radical_dimension_counts"].values())+stats["nondegenerate_subspaces"] != stats["isotropic_subspaces"]):
        raise ValueError("Inconsistent exact raw subspace counts")
    if (any(not key.isdigit() or not 1 <= int(key) <= q or type(n) is not int or n < 0
            for key, n in stats["radical_dimension_counts"].items())
            or any(type(stats[key]) is not int or stats[key] < 0
                   for key in ("eligible_supports", "nondegenerate_subspaces", "canonical_output_orbits"))
            or stats["isotropic_subspaces"] is not None and (type(stats["isotropic_subspaces"]) is not int or stats["isotropic_subspaces"] < 0)):
        raise ValueError("Invalid exact finite count")
    positives = branch["positive_supports"]
    if bool(positives) != bool(stats["nondegenerate_subspaces"]):
        raise ValueError("Positive support set and nondegenerate count disagree")
    seen, outputs = set(), set()
    for positive in positives:
        index = positive["selected_support_index"]
        keys = normalized_keys(positive["output_keys"])
        if index in seen or not start <= index < end or not keys:
            raise ValueError("Positive support is duplicate or outside its input interval")
        seen.add(index)
        outputs.update(keys)
    witness_outputs = normalized_keys([row["output_key_words"] for row in branch["witnesses"]])
    if outputs != witness_outputs or len(outputs) != stats["canonical_output_orbits"]:
        raise ValueError("A finite output lacks a witness or is absent from the positive profile")
    return local


def check_quotient_cover(branches, profile, q):
    selected = {d: n for d, n in profile.items() if d >= (5 if q == 5 else 1)}
    for dimension in selected:
        if dimension < q or q == 5 and dimension > 14:
            continue
        active = [row for row in branches if row["q"] == q
                  and row["quotient_dimension_interval"][0] <= dimension <= row["quotient_dimension_interval"][1]]
        intervals, total = [], Counter()
        for branch in active:
            r = branch["support_range"]
            if r["count"]:
                intervals.append((r["start"], r["end_exclusive"]))
            total.update(dict(branch["input_quotient_profile"]))
        check_intervals(intervals, sum(selected.values()))
        if dict(total) != selected:
            raise ValueError("Interval quotient profiles do not sum to the full logical input domain")


def verify(root, path):
    raw = path.read_bytes()
    data = json.loads(raw)
    if (data["schema"] != "finite-length52-small-quotient-blocks-v1" or data["matrix_scope"] != "full_projective"
            or data["protocol_length_interval"] != [51, 52]):
        raise ValueError("Wrong finite protocol sector")
    domain_path = root / data["source_domain"]
    domain_raw = domain_path.read_bytes()
    if hashlib.sha256(domain_raw).hexdigest() != data["source_domain_sha256"]:
        raise ValueError("Source domain binding changed")
    domain = json.loads(domain_raw)
    if (domain["schema"] != "finite-length52-protocol-source-domain-v1"
            or domain["minimum_protocol_length"] != 51 or domain["maximum_protocol_length"] != 52
            or domain["matrix_scope"] != "full_projective"):
        raise ValueError("Wrong source interval scope")
    catalogue = json.loads((root / "data/spaces/index.json").read_bytes())
    length, = [row for row in catalogue["lengths"] if row["c"] == 52]
    expected = {row["m"]: row["expected_classes"] for row in length["sectors"] if row["expected_classes"]}
    intervals = defaultdict(list)
    for i, block in enumerate(domain["main_blocks"]):
        if block["index"] != i or block["c"] != 52:
            raise ValueError("Noncontiguous finite source blocks")
        intervals[block["m"]].append(tuple(block["source_interval"]))
    for source in domain["separate_sources"]:
        if source["c"] != 52:
            raise ValueError("Unexpected separate parent length")
        intervals[source["m"]].append((source["index"], source["index"]+1))
    if set(intervals) != set(expected) or sum(expected.values()) != domain["source_spaces"]:
        raise ValueError("Incomplete source dimension partition")
    for m, count in expected.items():
        check_intervals(intervals[m], count)
    profile_path = root / data["profile_dataset"]
    profiles = json.loads(profile_path.read_bytes())["profiles"]
    frontier_path = root / "data/protocols/pareto_frontier.json"
    verify_frontier(frontier_path)
    verify_low_q(root)
    released = json.loads(frontier_path.read_bytes())
    outputs = {(row["q"], tuple(int(word, 0) for word in row["tensor_key_words"])): row["output_id"]
               for row in released["outputs"]}
    tensors = {row["output_id"]: LogicalTensor(tuple(int(s, 2) for s in row["generator_matrix_rows"][:row["q"]]))
               for row in released["protocols"]}
    source_map = {row["index"]: row for row in domain["main_blocks"] if row["logical_cover"] == "raw_small_quotients"}
    seen_blocks, seen_profiles = set(), set()
    witnesses, matrices, totals, marked_branches = [], {}, Counter(), []
    counts = defaultdict(Counter)
    for block in data["blocks"]:
        index = block["index"]
        if index in seen_blocks or index not in source_map:
            raise ValueError("Unknown or repeated source block")
        seen_blocks.add(index)
        source = source_map[index]
        if any(block[key] != source[key] for key in ("c", "m", "source_interval")) or block["protocol_length_interval"] != [51, 52]:
            raise ValueError("Finite census and source interval differ")
        first, end = block["source_interval"]
        count, pointed = end-first, block["pointing_supports"]
        parities = block["parity_counts"]
        if (len(parities) != 5 or any(type(v) is not int or v < 0 for v in parities)
                or sum(parities) != pointed or parities[2] or parities[4] or parities[3] != count):
            raise ValueError("Incomplete pointing parity partition")
        if "pointed_supports" in source and source["pointed_supports"] != pointed:
            raise ValueError("Wrong source pointing cardinality")
        if block["pointing_mode"] == "all_origins":
            totals["translation_period_sources"] += check_raw_pointings(block["m"], count, pointed, parities)
        elif block["pointing_mode"] != "affine_automorphism_orbits":
            raise ValueError("Unknown origin quotient")
        distribution = {3: histogram(block["distance_three_quotient_profile"], pointed),
                        4: histogram(block["distance_four_quotient_profile"], pointed)}
        deferred = sum(n for d, n in distribution[3].items() if d >= 15)
        if deferred != source.get("larger_quotient_supports", 0):
            raise ValueError("Unaccounted larger quotient domain")
        totals.update(source_spaces=count, pointed_supports=pointed, larger_quotient_supports=deferred)
        for position, branch in enumerate(block["branches"]):
            q = branch["q"]
            if branch["enumeration_mode"] == "alternating_tensor_exclusion" and not block.get("constant_stabiliser_d4_closure"):
                raise ValueError("Alternating exclusion lacks a constant-stabiliser certificate")
            try:
                check_branch(branch, distribution[3 if q == 5 else 4], count)
            except ValueError as error:
                raise ValueError(f"Block {index}, branch {position}, q={q}: {error}") from error
            totals["branches"] += 1
            if branch["statistics"] is not None:
                counts[q].update({key: branch["statistics"][key] for key in (
                    "raw_enumerated_supports", "marked_orbit_enumerated_supports", "isotropic_subspaces", "nondegenerate_subspaces")
                    if branch["statistics"][key] is not None})
                if branch["statistics"]["marked_orbit_enumerated_supports"]:
                    marked_branches.append(dict(block=index, branch=position, q=q))
            if q == 5:
                for positive in branch["positive_supports"]:
                    i = positive["profile_index"]
                    if i in seen_profiles or not 0 <= i < len(profiles) or profiles[i]["index"] != i:
                        raise ValueError("Unknown or repeated complete q5 support profile")
                    if normalized_keys(positive["output_keys"]) != normalized_keys(profiles[i]["output_keys"]):
                        raise ValueError("Q5 successor profile differs from the finite census")
                    seen_profiles.add(i)
            for local, witness in enumerate(branch["witnesses"]):
                rows = tuple(witness["generator_matrix_rows"])
                key, = normalized_keys([witness["output_key_words"]])
                identifier = outputs.get((q, key))
                if identifier is None or witness["q"] != q:
                    raise ValueError("Unrecognised witness output")
                cache_key = q, rows, identifier
                if cache_key not in matrices:
                    masks, columns = check_matrix(rows, q)
                    tensor = LogicalTensor(masks[:q])
                    basis = find_equivalence_basis(tensor, tensors[identifier]) if tensor.intrinsic() else None
                    distance, coefficient = distance_up_to(columns, q)
                    if basis is None or distance is None:
                        raise ValueError("Invalid intrinsic witness output or unverified exact distance")
                    matrices[cache_key] = len(columns), len(rows), distance, coefficient, basis
                n, S, distance, coefficient, basis = matrices[cache_key]
                if ((witness["n"], witness["S"], witness["d_Z"], witness["error_coefficient"]) != (n, S, distance, coefficient)
                        or not 51 <= n <= 52 or distance < branch["minimum_distance"]):
                    raise ValueError("Incorrect exact witness parameters")
                dominators = [row["index"] for row in released["protocols"] if row["output_id"] == identifier
                              and row["d_Z"] == distance and row["n"] <= n and row["S"] <= S]
                if not dominators:
                    raise ValueError("Finite witness is not covered by the released frontier")
                witnesses.append(dict(block=index, branch=position, witness=local, dominator=min(dominators),
                                      output_basis=list(basis)))
        for q in range(1, 6):
            check_quotient_cover(block["branches"], distribution[3 if q == 5 else 4], q)
    if seen_blocks != set(source_map) or seen_profiles != set(range(len(profiles))):
        raise ValueError("Incomplete finite source or positive-output profile cover")
    from verify_cubic_distance_four import verify as verify_cubic
    cubic = verify_cubic(root, path, root / "certificates/length52_cubic_blocks")
    if marked_branches:
        raise ValueError("A higher-distance orbit branch still lacks a direct mathematical closure")
    return dict(schema="finite-length52-small-quotient-verification-v1", status="pass",
        source_domain_sha256=hashlib.sha256(domain_raw).hexdigest(), data_sha256=hashlib.sha256(raw).hexdigest(),
        profile_sha256=hashlib.sha256(profile_path.read_bytes()).hexdigest(),
        frontier_sha256=hashlib.sha256(frontier_path.read_bytes()).hexdigest(),
        source_blocks=len(seen_blocks), q5_positive_supports=len(seen_profiles), totals=dict(totals),
        census_counts={str(q): dict(row) for q, row in counts.items()},
        witnesses_checked=len(witnesses), distinct_witness_checks=len(matrices), witness_checks=witnesses,
        complete_source_intervals_bound=True, raw_q5_small_quotient_domain_complete=True,
        exact_d3_small_output_pareto_exclusions_verified=True,
        separately_certified_sources=domain["separate_sources"],
        selected_cubic_source_blocks=sum(row["logical_cover"] == "selected_cubic_supports" for row in domain["main_blocks"]),
        selected_cubic_censuses_separate_obligation=True, larger_quotient_censuses_separate_obligation=True,
        higher_logical_dimensions_separate_obligation=True,
        higher_distance_marked_branches=marked_branches,
        higher_distance_marked_census_replay_separate_obligation=bool(marked_branches),
        constant_stabiliser_d4_sources=len(cubic["sources"]), constant_stabiliser_d4_exclusions_verified=True,
        all_censuses_freshly_recomputed=False, is_global_completeness_certificate=False)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=root / "data/protocol_sectors/length52_small_quotient_blocks.json")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify(root, args.input)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps({k: v for k, v in result.items() if k != "witness_checks"}, indent=2))


if __name__ == "__main__":
    main()
