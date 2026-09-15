"""Verify the presentation-only update against the original certified frontier."""

import argparse
import hashlib
import json
from pathlib import Path

from output_factorisation import verify_presentation_update


def verify(original_path, current_path):
    raw = original_path.read_bytes()
    current = json.loads(current_path.read_bytes())
    if hashlib.sha256(raw).hexdigest() != current["output_representative_presentation"]["previous_catalogue_sha256"]:
        raise ValueError("Incorrect original catalogue hash")
    return verify_presentation_update(json.loads(raw), current)


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original", type=Path, default=root / "provenance/pareto_frontier_before_factorisation.json")
    parser.add_argument("--catalogue", type=Path, default=root / "data/protocols/pareto_frontier.json")
    args = parser.parse_args()
    print(json.dumps(verify(args.original, args.catalogue), indent=2))
