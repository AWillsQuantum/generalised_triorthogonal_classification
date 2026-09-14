"""Build standalone fixed-width C++20 space-enumeration kernels."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, default=ROOT / "build")
    parser.add_argument("--cxx", default="g++")
    parser.add_argument("--select", nargs="*", help="Optional executable stems; default builds all")
    args = parser.parse_args()
    build = args.build_dir.resolve()
    build.mkdir(parents=True, exist_ok=True)
    sources = sorted(p for directory in ("length54", "contractions")
                     for p in (ROOT / directory).glob("*.cpp"))
    if len({p.stem for p in sources}) != len(sources):
        raise ValueError("Repeated native executable name")
    if args.select:
        unknown = set(args.select) - {p.stem for p in sources}
        if unknown:
            raise ValueError(f"Unknown kernels: {sorted(unknown)}")
        sources = [p for p in sources if p.stem in args.select]
    records = []
    for source in sources:
        output = build / source.stem
        subprocess.run([args.cxx, "-O3", "-DNDEBUG", "-std=c++20", "-fopenmp",
                        str(source), "-o", str(output)], check=True)
        records.append(dict(source=source.relative_to(ROOT).as_posix(),
                            sha256=hashlib.sha256(source.read_bytes()).hexdigest()))
        print(source.stem, flush=True)
    report = dict(status="pass", claim="all listed source units compile independently",
                  mathematical_validation_claimed=False, sources=records)
    (build / "build_verification.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
