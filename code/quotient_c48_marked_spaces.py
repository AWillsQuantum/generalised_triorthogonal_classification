"""Form the complete marked candidate union and its exact affine quotient."""

import argparse
import json
from pathlib import Path
import shutil

from finite_space_census import quotient_generated_domain
from prepare_zero_fibre_input import file_digest
from verify_protocol_cover import require


def run(root, generated, native, graph, work):
    certificate = generated / "complete.json"
    data = json.loads(certificate.read_bytes())
    require(data["status"] == "pass" and data["complete_source_domain"]
            and data["source_count"] == 152
            and [row["source_index"] for row in data["records"]] == list(range(152)),
            "Incomplete marked source cover")
    work.mkdir(parents=True, exist_ok=True)
    path = work / "candidates.bin"
    count = 0
    with path.open("xb") as stream:
        for row in data["records"]:
            source = generated / Path(row["task"]).stem / "full_rank.bin"
            require(source.stat().st_size == 32*row["full_rank_orbits"]
                    and file_digest(source) == row["full_rank_binary_sha256"], "Different marked candidate stream")
            with source.open("rb") as incoming:
                shutil.copyfileobj(incoming, stream, 1 << 20)
            count += row["full_rank_orbits"]
    domain = dict(candidates=count, candidate_sha256=file_digest(path),
                  marked_generation_sha256=file_digest(certificate),
                  source_domain_sha256=data["input_domain_sha256"],
                  raw_pairs=sum(row["raw_pairs"] for row in data["records"]),
                  full_rank_pairs=sum(row["full_rank_pairs"] for row in data["records"]),
                  deficient_pairs=sum(row["deficient_pairs"] for row in data["records"]))
    (work / "domain.json").write_text(json.dumps(domain, indent=2)+"\n", encoding="ascii")
    print(json.dumps(dict(stage="marked_union", **domain)), flush=True)
    return quotient_generated_domain(root, 48, 8, domain, native, graph, work)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--generated", type=Path, required=True)
    parser.add_argument("--native-directory", type=Path, required=True)
    parser.add_argument("--graph-native", type=Path, required=True)
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.root.resolve(), args.generated.resolve(), args.native_directory.resolve(),
                 args.graph_native.resolve(), args.work_directory.resolve())
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    print(json.dumps({k:v for k,v in result.items() if k in
                     ("status", "candidate_count", "affine_classes", "new_classes", "unreached_old_targets", "elapsed_seconds")}))


if __name__ == "__main__":
    main()
