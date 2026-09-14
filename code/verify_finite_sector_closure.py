"""Compose exact-distance dominance and higher-logical-dimension closure on finite supports."""

import argparse
import hashlib
import json
from pathlib import Path

from logical_spaces import logical_label_space
from verify_primitive_restrictions import verify as verify_primitive_restrictions
from verify_protocol_case_census import normalized_keys
from verify_protocols import verify as verify_witnesses


def verify(root):
    root = Path(root).resolve()
    bindings = {}
    def load(name):
        path = (root / name).resolve()
        if not path.is_relative_to(root):
            raise ValueError("External closure dependency")
        raw = path.read_bytes()
        bindings[name] = hashlib.sha256(raw).hexdigest()
        return json.loads(raw)
    frontier_name = "data/protocols/pareto_frontier.json"
    frontier = load(frontier_name)
    verify_witnesses(root / frontier_name)
    primitive = verify_primitive_restrictions(load("certificates/tensors/q5_orbits.json"))
    forbidden_q5 = normalized_keys(r["output_key_words"] for r in primitive["restrictions"])
    gate_certificate = load("certificates/output_profiles/certificate.json")
    for name, digest in gate_certificate["files"].items():
        relative = "certificates/output_profiles/" + name
        load(relative)
        if bindings[relative] != digest:
            raise ValueError("The finite output-profile gate binding changed")
    gate = load("certificates/output_profiles/q7_anchored_exclusions.json")
    if (gate_certificate["status"] != "pass" or gate["status"] != "pass"
            or gate["extensions_per_anchor"] != 1 << 22
            or gate["candidates_tested"] != gate["q6_target_key_count"]*(1 << 22)):
        raise ValueError("Incomplete anchored extension gate")
    excluded_profiles = [normalized_keys([[key] for key in row["q6_keys"]])
                         for row in gate["support_profiles"] if row["surviving_anchored_extensions"] == 0]
    sectors = []
    configurations = (
        ("length54_q6_small_quotient", "q6_finite_census.json"),
        ("length52_q6_large_quotient", "q6_large_quotient_census.json"),
    )
    for name, certificate_name in configurations:
        data_name = "data/protocol_sectors/" + name + ".json"
        data = load(data_name)
        certificate = load("certificates/" + certificate_name)
        if (certificate["status"] != "pass" or certificate["input_sha256"] != bindings[data_name]
                or certificate["cases"] != len(data["cases"])
                or not certificate["every_support_output_profile_recomputed"]
                or certificate["logical_qubits"] != 6 or certificate["minimum_distance"] != 3):
            raise ValueError("A complete q6 support census is missing")
        comparisons = []
        for case in data["cases"]:
            points, h = case["points"], case["ambient_dimension"]
            if logical_label_space(points, h, 4)["quotient_dimension"] != 0:
                raise ValueError("A possible exact distance above three remains")
            if normalized_keys(case["q5_output_keys"]) & forbidden_q5:
                raise ValueError("A primitive seven-dimensional output remains possible")
            q6 = normalized_keys(case["expected"]["output_keys"])
            if not any(q6.issubset(profile) for profile in excluded_profiles):
                raise ValueError("Nonprimitive seven-dimensional outputs have not been excluded")
            for key in sorted(q6):
                identifier = "Q6_" + "_".join(f"{word:016x}" for word in key)
                witnesses = [row for row in frontier["protocols"]
                             if row["output_id"] == identifier and row["d_Z"] == 3
                             and row["n"] <= len(points) and row["S"] <= 6+h]
                if not witnesses:
                    raise ValueError("A finite-sector metric has no dominating catalogue point")
                witness = min(witnesses, key=lambda row: (row["n"], row["S"], row["index"]))
                comparisons.append(dict(case=case["index"], output_id=identifier,
                                        n=len(points), S=6+h, exact_distance=3,
                                        incumbent_protocol_index=witness["index"],
                                        relation="equal_metrics" if (witness["n"], witness["S"]) == (len(points), 6+h)
                                        else "strictly_dominated"))
        sectors.append(dict(id=name, supports=len(data["cases"]), comparisons=comparisons,
                            q6_has_no_unrepresented_pareto_metrics=True,
                            primitive_and_nonprimitive_q7_absent=True,
                            all_q_at_least_seven_absent=True))
    return dict(schema="finite-sector-exact-distance-closure-v1", status="pass", sectors=sectors,
                dependencies=bindings, exact_distance_four_and_higher_excluded_by_linear_equations=True,
                higher_dimension_rule="Primitive restriction theorem and descent through nondegenerate hyperplanes",
                scope="Only the listed finite supports and their complete output profiles",
                is_global_completeness_certificate=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.root)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="ascii")
    print(json.dumps(dict(status=result["status"], supports=sum(r["supports"] for r in result["sectors"]),
                         metric_comparisons=sum(len(r["comparisons"]) for r in result["sectors"]),
                         global_completeness=False)))
