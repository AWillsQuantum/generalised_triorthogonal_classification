"""Check block-to-support profile refinement and the exact pointed-support quotient."""

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import subprocess

from verify_profile_exclusions import necessary_profiles, verify as verify_profiles
from verify_protocol_case_census import normalized_keys, write_pointed_cases


def verify(root, native, work):
    root = Path(root).resolve()
    bindings = {}
    def load(name):
        path = (root / name).resolve()
        if not path.is_relative_to(root):
            raise ValueError("External refinement dependency")
        raw = path.read_bytes()
        bindings[name] = hashlib.sha256(raw).hexdigest()
        return json.loads(raw)
    census = load("certificates/output_profiles/q6_necessary_profiles.json")
    necessary = necessary_profiles(census)
    blocks = load("data/protocol_profiles/length54_small_quotient_blocks.json")
    block_result = verify_profiles(blocks, necessary)
    data = load("data/protocol_sectors/length54_q6_candidates.json")
    refined = load(data["refined_profile_dataset"])
    refined_result = verify_profiles(refined, necessary)
    survivors = {r["index"] for r in refined_result["candidates"]}
    canonical = load(data["canonical_case_dataset"])
    by_block = defaultdict(list)
    identities = {}
    for row in refined["profiles"]:
        identity = row["block_index"], row["block_support_index"]
        if identity in identities:
            raise ValueError("Repeated individual support profile")
        identities[identity] = row
        by_block[row["block_index"]].append(row)
    if set(by_block) != {r["index"] for r in block_result["candidates"]}:
        raise ValueError("The refinement does not cover exactly all unresolved blocks")
    for index, rows in by_block.items():
        block = blocks["profiles"][index]
        union = set().union(*(normalized_keys(row["output_keys"]) for row in rows))
        if len(rows) != block["positive_supports"] or union != normalized_keys(block["output_keys"]):
            raise ValueError("Block profile union or positive-support count does not reproduce")
    if [r["index"] for r in data["cases"]] != list(range(len(data["cases"]))):
        raise ValueError("Missing candidate indices")
    recovered = set()
    for row in data["cases"]:
        profile = identities[row["block_index"], row["block_support_index"]]
        if profile["index"] not in survivors or profile["index"] in recovered:
            raise ValueError("Candidate support absent from, or repeated in, the refinement")
        recovered.add(profile["index"])
        if normalized_keys(row["q5_output_keys"]) != normalized_keys(profile["output_keys"]):
            raise ValueError("Candidate q5 profile does not match its refined profile")
        target = row["canonical_case"]
        if type(target) is not int or not 0 <= target < len(canonical["cases"]):
            raise ValueError("Invalid canonical case reference")
        if normalized_keys(row["q5_output_keys"]) != normalized_keys(canonical["cases"][target]["q5_output_keys"]):
            raise ValueError("Equivalent pointed supports have different complete q5 profiles")
    if recovered != survivors:
        raise ValueError("Some refined candidates have no pointed geometry")
    combined = [dict(row, index=i) for i, row in enumerate(data["cases"] + canonical["cases"])]
    work = Path(work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    manifest, output = work / "supports.utsp", work / "canonical.json"
    with manifest.open("wb") as stream:
        write_pointed_cases(stream, dict(maximum_protocol_length=54, cases=combined))
    if output.exists():
        raise ValueError("Canonical output already exists")
    subprocess.run([str(native), "manifest-support-canonical-keys", "--input", str(manifest),
                    "--output", str(output), "--workers", "1"], check=True, capture_output=True)
    result = json.loads(output.read_bytes())
    rows = result["records"]
    if (result["status"] != "complete" or len(rows) != len(combined)
            or [r["support_index"] for r in rows] != list(range(len(combined)))):
        raise ValueError("Incomplete pointed-support canonicalisation")
    keys = [tuple(int(word, 16) for word in row["key_words"]) for row in rows]
    start = len(data["cases"])
    if len(set(keys[start:])) != len(canonical["cases"]):
        raise ValueError("Canonical case dataset repeats a pointed-support class")
    for row in data["cases"]:
        if keys[row["index"]] != keys[start + row["canonical_case"]]:
            raise ValueError("Candidate-to-census pointed equivalence fails")
    if {r["canonical_case"] for r in data["cases"]} != set(range(len(canonical["cases"]))):
        raise ValueError("The pointed quotient contains an unaccounted class")
    return dict(schema="output-profile-refinement-replay-v1", status="pass",
                blocks=blocks["profile_count"], excluded_blocks=block_result["excluded"],
                refined_blocks=len(by_block), refined_positive_supports=len(identities),
                candidate_supports=len(data["cases"]), canonical_support_classes=len(canonical["cases"]),
                class_multiplicity_histogram=dict(Counter(Counter(r["canonical_case"] for r in data["cases"]).values())),
                complete_block_unions_reproduced=True, exact_pointed_equivalences_recomputed=True,
                dependencies=bindings, is_global_completeness_certificate=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.root, args.native, args.work_directory)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="ascii")
    print(json.dumps({k: v for k, v in result.items() if k != "dependencies"}, indent=2))
