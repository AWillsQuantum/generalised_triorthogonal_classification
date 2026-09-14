"""Build the same finite affine algorithms for a specified length and dimension."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--length", type=int, required=True)
    parser.add_argument("--dimension", type=int, required=True)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--cxx", default="g++")
    parser.add_argument("--include-lifts", action="store_true", help="Also build accelerated m=9,10,11 lift routines")
    args = parser.parse_args()
    c, m = args.length, args.dimension
    if not 16 <= c <= 54 or c % 2 or not 8 <= m <= 15:
        raise ValueError("Unsupported finite affine configuration")
    root = Path(__file__).resolve().parent
    build = args.build_dir.resolve()
    build.mkdir(parents=True, exist_ok=True)
    sources = {}
    names = [("signature", "signature"), ("exact_quotient", "exact_quotient")]
    if args.include_lifts:
        if m not in (9, 10, 11) or c not in (50, 52, 54):
            raise ValueError("The accelerated lift configurations are c=50,52,54 and m=9,10,11")
        names.append((f"extensions_m{m:02d}", "extensions"))
        if m == 11:
            names.append(("profile_m11", "profile"))
    for name, stem in names:
        source = root / (name+".cpp")
        executable = build / f"{stem}_c{c}_m{m:02d}"
        subprocess.run([args.cxx, "-O3", "-DNDEBUG", "-std=c++20", "-fopenmp",
                        f"-DUTSP_LENGTH={c}", f"-DUTSP_DIMENSION={m}",
                        f"-DUTSP_MINIMUM_FILTER={int(c != 50)}",
                        str(source), "-o", str(executable)], check=True)
        sources[source.name] = hashlib.sha256(source.read_bytes()).hexdigest()
        print(executable.name, flush=True)
    sources["configuration.hpp"] = hashlib.sha256((root / "configuration.hpp").read_bytes()).hexdigest()
    (build / f"build_c{c}_m{m:02d}.json").write_text(json.dumps(dict(status="pass", c=c, m=m, sources=sources), indent=2)+"\n")


if __name__ == "__main__":
    main()
