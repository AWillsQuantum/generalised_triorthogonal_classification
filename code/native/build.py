"""Build the C++20 mathematical kernels and their finite regression checks."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import zipfile

ROOT = Path(__file__).resolve().parent
BLISS_SHA256 = "acc8b98034f30fad24c897f365abd866c13d9f1bb207e398d0caf136875972a4"
BLISS_UNITS = ("abstractgraph", "bliss_C", "defs", "digraph", "graph", "orbit",
               "partition", "uintseqhash", "utils")


def run(command):
    subprocess.run([str(value) for value in command], check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, default=ROOT / "build")
    parser.add_argument("--cxx", default="g++")
    parser.add_argument("--ar", default="ar")
    parser.add_argument("--gmp-include", type=Path)
    parser.add_argument("--gmp-library", help="Library path; default is the system -lgmp")
    parser.add_argument("--test", action="store_true")
    args = parser.parse_args()
    build = args.build_dir.resolve()
    build.mkdir(parents=True, exist_ok=True)
    archive = ROOT / "vendor/bliss-0.77.zip"
    if hashlib.sha256(archive.read_bytes()).hexdigest() != BLISS_SHA256:
        raise ValueError("Unexpected Bliss source archive")
    dependency = build / "dependencies"
    dependency.mkdir(exist_ok=True)
    with zipfile.ZipFile(archive) as source:
        for member in source.infolist():
            destination = (dependency / member.filename).resolve()
            if not destination.is_relative_to(dependency.resolve()):
                raise ValueError("Unsafe third-party archive path")
            if (member.external_attr >> 16) & 0o170000 == 0o120000:
                raise ValueError("Third-party archive contains a symbolic link")
        source.extractall(dependency)
    bliss = dependency / "bliss-0.77/src"
    include = build / "include"
    (include / "bliss").mkdir(parents=True, exist_ok=True)
    for header in bliss.glob("*.hh"):
        shutil.copyfile(header, include / "bliss" / header.name)
    flags = ["-O3", "-DNDEBUG", "-std=c++20", "-pthread", "-DUTSP_HAVE_BLISS", "-DBLISS_USE_GMP",
             "-I" + str(ROOT / "include"), "-I" + str(include)]
    if args.gmp_include:
        flags.append("-I" + str(args.gmp_include.resolve()))
    library = args.gmp_library or "-lgmp"
    bliss_objects = []
    for name in BLISS_UNITS:
        obj = build / ("bliss_" + name + ".o")
        run([args.cxx, *flags, "-c", bliss / (name + ".cc"), "-o", obj])
        bliss_objects.append(obj)
    static = build / "libbliss.a"
    run([args.ar, "rcs", static, *bliss_objects])
    objects = {}
    for source in sorted((ROOT / "src").glob("*.cpp")):
        obj = build / (source.stem + ".o")
        run([args.cxx, *flags, "-Wall", "-Wextra", "-Wpedantic", "-c", source, "-o", obj])
        objects[source.stem] = obj
    binary = build / "utsp-native"
    run([args.cxx, *objects.values(), "-pthread", static, library, "-o", binary])
    for source in sorted((ROOT / "auxiliary").glob("*.cpp")):
        core = [obj for name, obj in objects.items() if name != "utsp_native"]
        run([args.cxx, *flags, source, *core, static, library, "-o", build / source.stem])
    source_files = [Path(__file__), archive]
    for pattern in ("src/*.cpp", "include/utsp/*.hpp", "auxiliary/*.cpp", "tests/*.cpp"):
        source_files.extend(sorted(ROOT.glob(pattern)))
    result = dict(status="built", source_files=len(objects),
                  auxiliary_programs=[p.stem for p in sorted((ROOT / "auxiliary").glob("*.cpp"))],
                  source_sha256={path.relative_to(ROOT).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
                                 for path in source_files}, tests={})
    if args.test:
        core = [obj for name, obj in objects.items() if name != "utsp_native"]
        programs = [(binary, ["self-test"])]
        for name in ("primitive_parent_regression", "parent_invariant_cache", "small_quotient_group"):
            target = build / name
            linked = core if name == "primitive_parent_regression" else [p for p in core if p.stem != "marked_code_canonical"]
            run([args.cxx, *flags, ROOT / "tests" / (name + ".cpp"), *linked,
                 static, library, "-o", target])
            programs.append((target, []))
        for executable, arguments in programs:
            output = subprocess.run([str(executable), *arguments], check=True, text=True, capture_output=True)
            result["tests"][executable.name] = dict(returncode=output.returncode, stdout=output.stdout)
            print(output.stdout, flush=True)
        result["status"] = "pass"
    (build / "verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
