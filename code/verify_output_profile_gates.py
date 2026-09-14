"""Recompute finite necessary output profiles and anchored dimension-seven exclusions."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

PROFILES = (
    (0xA08A400,),
    (0xA08A400, 0x400A08A400),
    (0x40200A000, 0x4002090400),
    (0x4002090400,),
    (0x400A08A400,),
)


def run(binary, arguments):
    completed = subprocess.run([str(binary), *arguments], text=True, capture_output=True, check=True)
    result = json.loads(completed.stdout)
    if result["status"] != "pass":
        raise ValueError("A finite profile computation did not cover its full domain")
    for name in ("elapsed_seconds", "workers"):
        result.pop(name, None)
    return result


def verify(binaries, output, workers):
    if workers < 1:
        raise ValueError("Workers must be positive")
    output.mkdir(parents=True, exist_ok=True)
    files = {}
    def publish(name, data):
        raw = (json.dumps(data, indent=2) + "\n").encode()
        (output / name).write_bytes(raw)
        files[name] = hashlib.sha256(raw).hexdigest()
        print(json.dumps(dict(completed=name, sha256=files[name])), flush=True)
    worker_args = ["--workers", str(workers)]
    all_profiles = run(binaries / "q6_hyperplane_output_profile_census", worker_args)
    if (all_profiles["q5_tensor_count"] != 1 << 25 or all_profiles["q5_orbit_count"] != 88
            or all_profiles["q5_nondegenerate_orbit_count"] != 65
            or all_profiles["q6_gram_class_count"] != 10
            or all_profiles["q6_exact_tensor_count"] != 10 * (1 << 20)
            or all_profiles["q6_tensors_tested"] != all_profiles["q6_exact_tensor_count"]
            or len(all_profiles["profiles"]) != all_profiles["minimal_output_profiles"]):
        raise ValueError("Incomplete q6 normal-form census")
    publish("q6_necessary_profiles.json", all_profiles)
    profile_args = [value for profile in PROFILES for value in
                    ("--profile", ",".join(f"0x{key:016x}" for key in profile))]
    binary = binaries / "q7_q6_output_profile_gate"
    uniqueness = run(binary, [*worker_args, *profile_args, "--q6-profile-census-only"])
    expected_keys = {key for profile in PROFILES for key in profile}
    if (uniqueness["q6_tensors_tested"] != 10 * (1 << 20)
            or uniqueness["q6_exact_tensor_count"] != 10 * (1 << 20)
            or uniqueness["q6_gram_class_count"] != 10
            or not uniqueness["all_target_profiles_distinct"]
            or not uniqueness["all_matching_normal_forms_have_claimed_key"]
            or {int(row["canonical_key"], 16) for row in uniqueness["targets"]} != expected_keys
            or any(row["canonical_key_mismatches"] != 0
                   or row["direct_canonical_checks"] != row["normal_form_matches"]
                   or row["normal_form_matches"] == 0 for row in uniqueness["targets"])):
        raise ValueError("The accelerated q6 recogniser is not fully certified")
    publish("q6_profile_recogniser.json", uniqueness)
    # This optimisation is enabled only after the complete recogniser census above.
    extension = run(binary, [*worker_args, *profile_args, "--trust-certified-q6-profiles"])
    actual_profiles = {tuple(int(key, 16) for key in row["q6_keys"])
                       for row in extension["support_profiles"]}
    if (extension["extensions_per_anchor"] != 1 << 22
            or extension["exact_extensions_per_anchor"] != 1 << 22
            or extension["candidates_tested"] != len(expected_keys) * (1 << 22)
            or extension["q6_target_key_count"] != len(expected_keys)
            or actual_profiles != set(PROFILES)
            or any(row["surviving_anchored_extensions"] != 0 for row in extension["support_profiles"])):
        raise ValueError("Anchored q7 extension exclusion is incomplete or has a survivor")
    publish("q7_anchored_exclusions.json", extension)
    report = dict(status="pass", q6_normal_form_tensors=10 * (1 << 20),
                  q7_anchored_tensors=len(expected_keys) * (1 << 22),
                  q6_necessary_profile_count=all_profiles["minimal_output_profiles"],
                  q7_excluded_profiles=[list(profile) for profile in PROFILES],
                  primitive_q7_sector_excluded_by_this_test=False,
                  is_global_protocol_completeness_certificate=False, files=files)
    publish("certificate.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary-directory", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    verify(args.binary_directory, args.output_directory, args.workers)
