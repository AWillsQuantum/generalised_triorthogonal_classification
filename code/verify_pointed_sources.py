"""Bind explicitly pointed protocol supports to compact parent catalogue records."""

import argparse
import hashlib
import json
from pathlib import Path

from selected_spaces import read_selected
from verify_small_source_closure import all_pointings


def check_raw_block_position(row, block, parent_points):
    source = row["space"]
    first, end = block["source_interval"]
    if (block["index"] != row["block_index"] or (block["c"], block["m"]) != (source["c"], source["m"])
            or not first <= source["index"] < end or source["c"] != 54
            or block["pointing_mode"] != "all_origins" or block["protocol_length_interval"] != [53, 54]
            or any(d < 5 for d, count in block["distance_three_quotient_profile"] if count)):
        raise ValueError("The selected input does not have a full raw source-block index")
    cases = all_pointings(parent_points, source["m"], 54)
    cases.sort(key=lambda r: (r["ambient_dimension"], r["points"], r["parity_case"], r["origin"]))
    count = (1 << source["m"])+55
    if len(cases) != count or block["pointing_supports"] != (end-first)*count:
        raise ValueError("The raw block has a different per-source pointing cardinality")
    offset = row["block_support_index"]-(source["index"]-first)*count
    if not 0 <= offset < count:
        raise ValueError("Pointed input index belongs to a different parent")
    case = cases[offset]
    if (case["points"] != row["points"] or case["ambient_dimension"] != row["ambient_dimension"]
            or case["parity_case"] != 0 or case["origin"] != source["origin"]):
        raise ValueError("Pointed geometry differs from its exact source-block position")
    return offset


def verify(root, data, block_data=None):
    cases = data["cases"]
    keys = {(row["space"]["c"], row["space"]["m"], row["space"]["index"]) for row in cases}
    supports, bindings = read_selected(root, keys)
    positions = []
    for row in cases:
        source = row["space"]
        key = source["c"], source["m"], source["index"]
        origin = source["origin"]
        if type(origin) is not int or not 0 <= origin < 1 << source["m"]:
            raise ValueError("Invalid origin")
        expected = sorted(x ^ origin for x in supports[key])
        if expected != row["points"] or row["ambient_dimension"] != source["m"]:
            raise ValueError("Pointed support does not reproduce from its catalogue source")
        if block_data is not None:
            positions.append(check_raw_block_position(row, block_data["blocks"][row["block_index"]], supports[key]))
    return dict(schema="pointed-catalogue-source-replay-v1", status="pass", cases=len(cases),
                distinct_parent_spaces=len(keys), pointing_case="even translated support",
                source_catalogue_sha256=bindings, every_point_reproduced=True,
                exact_source_block_positions_checked=block_data is not None,
                within_parent_pointing_offsets=positions,
                is_global_completeness_certificate=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--block-input", type=Path)
    args = parser.parse_args()
    raw = args.input.read_bytes()
    block_raw = args.block_input.read_bytes() if args.block_input else None
    result = verify(args.root, json.loads(raw), json.loads(block_raw) if block_raw else None)
    if block_raw is not None:
        result["block_input_sha256"] = hashlib.sha256(block_raw).hexdigest()
    result["input_sha256"] = hashlib.sha256(raw).hexdigest()
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="ascii")
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ("source_catalogue_sha256", "within_parent_pointing_offsets")}, indent=2))


if __name__ == "__main__":
    main()
