"""Independently check distance-three and distance-four quotient dimensions."""

import argparse
import hashlib
import json
from pathlib import Path

from logical_spaces import logical_label_space


def verify(data):
    rows = []
    cases = data["cases"]
    if [r["index"] for r in cases] != list(range(len(cases))):
        raise ValueError("The finite case set is not contiguous")
    for case in cases:
        profile = {}
        for distance, name in ((3, "distance_three_quotient_dimension"),
                                (4, "distance_four_quotient_dimension")):
            space = logical_label_space(case["points"], case["ambient_dimension"], distance)
            if space["quotient_dimension"] != case[name]:
                raise ValueError(f"Incorrect distance-{distance} quotient in case {case['index']}")
            profile[str(distance)] = space["quotient_dimension"]
        rows.append(dict(index=case["index"], quotient_dimensions=profile))
    return dict(schema="independent-distance-filtration-v1", status="pass",
                cases=len(rows), all_distance_four_quotients_zero=bool(rows) and all(
                    row["quotient_dimensions"]["4"] == 0 for row in rows),
                is_global_completeness_certificate=False, profiles=rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    raw = args.input.read_bytes()
    result = verify(json.loads(raw))
    result["input_sha256"] = hashlib.sha256(raw).hexdigest()
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="ascii")
    print(json.dumps({k: v for k, v in result.items() if k != "profiles"}, indent=2))
