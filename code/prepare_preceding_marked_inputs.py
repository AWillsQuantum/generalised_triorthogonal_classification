"""Write the complete marked-core inputs for length 50 or 52 from compact supports."""

import argparse
import hashlib
import json
from pathlib import Path

from prepare_space_contractions import task_text
from selected_spaces import read_selected
from verify_protocol_cover import Evidence, require


def prepare(root, length, directory):
    require(length in (50, 52), "Unsupported marked target length")
    evidence = Evidence(root)
    census = evidence.read(f"data/space_contractions/preceding_sectors/c{length}_m08.json")
    base = evidence.read("certificates/low_dimensional_spaces.json")
    require(base["status"] == "pass" and base["complete_and_pairwise_inequivalent"], "The complete affine base is required")
    forms = {(r["c"], r["m"], r["index"]): r for r in base["canonical_maps"]}
    cases = [(role, row) for role in ("primary", "alternative") for row in census[role+"_cores"]]
    selected = {(r["c"], r["m"], r["index"]) for _, r in cases}
    require(len(selected) == len(cases), "Duplicate marked core")
    supports, bindings = read_selected(evidence.root, selected)
    for path, digest in bindings.items():
        evidence.digest(path, digest)
    directory = Path(directory).resolve()
    require(not directory.exists(), "Use a new marked-input directory")
    directory.mkdir(parents=True)
    records = []
    for role, row in cases:
        c, m, i = row["c"], row["m"], row["index"]
        identifier = f"c{c}_m{m:02d}_{i:09d}"
        text, profile = task_text(identifier, supports[c, m, i], forms[c, m, i], length)
        require(role != "primary" or profile == row["profile"], "The primary native input profile differs")
        payload = text.encode("ascii")
        path = directory / (identifier+".task")
        path.write_bytes(payload)
        records.append(dict(c=c, m=m, index=i, role=role, task=path.name,
                            sha256=hashlib.sha256(payload).hexdigest(), **profile))
    result = dict(schema="preceding-marked-inputs-v1", status="pass", target_length=length,
        target_affine_dimension=8, primary_cores=sum(r["role"] == "primary" for r in records),
        alternative_cores=sum(r["role"] == "alternative" for r in records), records=records,
        dependencies=evidence.bindings, child_enumeration_performed=False)
    (directory / "domain.json").write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--length", type=int, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.root, args.length, args.output_directory)
    print(json.dumps({k: v for k, v in result.items() if k not in ("records", "dependencies")}, indent=2))


if __name__ == "__main__":
    main()
