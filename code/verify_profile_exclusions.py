"""Apply a complete necessary-profile census to finite support-profile datasets."""

import argparse
import hashlib
import json
from pathlib import Path

from verify_protocol_case_census import normalized_keys
from verify_primitive_restrictions import verify as verify_primitive_restrictions


def necessary_profiles(census):
    if (census["status"] != "pass" or census["q6_tensors_tested"] != 10*(1 << 20)
            or census["q6_exact_tensor_count"] != census["q6_tensors_tested"]):
        raise ValueError("Incomplete six-dimensional tensor census")
    profiles = [normalized_keys([[key] for key in row["q5_output_keys"]])
                for row in census["profiles"]]
    if len(profiles) != census["minimal_output_profiles"] or any(not p for p in profiles):
        raise ValueError("Incomplete necessary-profile list")
    if any(left.issubset(right) for i, left in enumerate(profiles)
           for j, right in enumerate(profiles) if i != j):
        raise ValueError("Necessary profiles do not form an inclusion antichain")
    return profiles


def verify(data, profiles):
    if (data["predecessor_logical_dimension"] != 5 or data["minimum_distance"] != 3
            or data["profile_scope"] not in ("individual_support", "block_union")
            or len(data["profiles"]) != data["profile_count"]
            or [r["index"] for r in data["profiles"]] != list(range(data["profile_count"]))):
        raise ValueError("Invalid support-profile domain")
    candidates = []
    for row in data["profiles"]:
        available = normalized_keys(row["output_keys"])
        if any(len(key) != 1 or key[0] >= 1 << 25 for key in available):
            raise ValueError("Invalid five-dimensional tensor key")
        matches = [i for i, profile in enumerate(profiles) if profile.issubset(available)]
        if matches:
            candidates.append(dict(index=row["index"], matched_necessary_profiles=matches))
    if len(candidates) != data["expected_q6_candidates"]:
        raise ValueError("The finite profile exclusions disagree with the retained result")
    return dict(id=data["id"], profiles=data["profile_count"], profile_scope=data["profile_scope"],
                excluded=data["profile_count"]-len(candidates), candidates=candidates,
                all_q6_absent_on_listed_profiles=not candidates)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    census_path = args.root / "certificates/output_profiles/q6_necessary_profiles.json"
    census_bytes = census_path.read_bytes()
    profiles = necessary_profiles(json.loads(census_bytes))
    primitive_path = args.root / "certificates/tensors/q5_orbits.json"
    primitive_bytes = primitive_path.read_bytes()
    primitive = verify_primitive_restrictions(json.loads(primitive_bytes))
    primitive_keys = normalized_keys(row["output_key_words"] for row in primitive["restrictions"])
    results, bindings = [], {}
    for path in sorted((args.root / "data/protocol_profiles").glob("*.json")):
        raw = path.read_bytes()
        bindings[path.relative_to(args.root).as_posix()] = hashlib.sha256(raw).hexdigest()
        data = json.loads(raw)
        result = verify(data, profiles)
        result["primitive_q7_candidate_indices"] = [row["index"] for row in data["profiles"]
            if normalized_keys(row["output_keys"]) & primitive_keys]
        result["primitive_q7_absent_on_listed_profiles"] = not result["primitive_q7_candidate_indices"]
        results.append(result)
    if not results:
        raise ValueError("No finite profile datasets were provided")
    report = dict(schema="finite-output-profile-exclusions-v1", status="pass",
                  necessary_profile_census_sha256=hashlib.sha256(census_bytes).hexdigest(),
                  primitive_restriction_census_sha256=hashlib.sha256(primitive_bytes).hexdigest(),
                  profile_source_sha256=bindings, sectors=results,
                  premise="The listed output profiles are complete on their stated support domains",
                  support_enumerations_recomputed=False,
                  is_global_completeness_certificate=False)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="ascii")
    print(json.dumps([dict(id=r["id"], profiles=r["profiles"], excluded=r["excluded"],
                           candidates=len(r["candidates"])) for r in results], indent=2))
