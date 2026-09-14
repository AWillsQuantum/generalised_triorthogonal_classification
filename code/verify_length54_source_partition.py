"""Compose source-interval covers without confusing them with full protocol closure."""

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path

from verify_high_dimensional_sector import check_intervals


def verify(root):
    root = Path(root).resolve()
    index_path = root / "data/spaces/index.json"
    catalogue = json.loads(index_path.read_bytes())
    length, = [entry for entry in catalogue["lengths"] if entry["c"] == 54]
    if not length["data_complete"]:
        raise ValueError("Incomplete length-54 source catalogue")
    expected = {sector["m"]: sector["expected_classes"] for sector in length["sectors"]}
    intervals, bindings = defaultdict(list), {}
    source_blocks = 0
    for stem, certificate_name in (("length54_low_dimensional_blocks", "low_dimensional_protocol_blocks"),
                                   ("length54_high_dimensional", "high_dimensional_protocol_sector")):
        data_path = root / "data/protocol_sectors" / f"{stem}.json"
        certificate_path = root / "certificates" / f"{certificate_name}.json"
        data_raw, certificate_raw = data_path.read_bytes(), certificate_path.read_bytes()
        data, certificate = json.loads(data_raw), json.loads(certificate_raw)
        digest = hashlib.sha256(data_raw).hexdigest()
        if certificate["status"] != "pass" or certificate["data_sha256"] != digest:
            raise ValueError("Unverified or changed finite source partition")
        bindings[data_path.relative_to(root).as_posix()] = digest
        bindings[certificate_path.relative_to(root).as_posix()] = hashlib.sha256(certificate_raw).hexdigest()
        for block in data["blocks"]:
            if block["c"] == 54:
                if block["m"] not in expected:
                    raise ValueError("Unexpected source dimension")
                intervals[block["m"]].append(block["source_interval"])
                source_blocks += 1
    for m, count in expected.items():
        check_intervals(intervals[m], count)
    return dict(schema="length54-protocol-source-partition-v1", status="pass",
                source_length=54, source_spaces=sum(expected.values()), source_blocks=source_blocks,
                source_catalogue_index_sha256=hashlib.sha256(index_path.read_bytes()).hexdigest(),
                sources_by_affine_dimension={str(m): count for m, count in expected.items()},
                source_intervals_disjoint_and_exhaustive=True,
                protocol_length_53_parent_length52_excluded_by="exact zero-column domination",
                only_pareto_relevant_pointing_cases=["even support without zero", "remove origin", "affine hyperplane"],
                supporting_files=bindings,
                separate_obligations=["complete shorter-length protocol frontier",
                                      "complete larger-quotient logical-subspace enumeration",
                                      "complete higher-logical-dimension successor closure",
                                      "unital-source classification completeness and inequivalence"],
                is_global_completeness_certificate=False)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify(root)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="ascii")
    print(json.dumps({key: value for key, value in result.items() if key != "supporting_files"}, indent=2))


if __name__ == "__main__":
    main()
