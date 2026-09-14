"""Enumerate marked orbits, then retain the full-affine-rank support union."""

import argparse
import json
from pathlib import Path
import subprocess
from time import perf_counter

import numpy as np

from prepare_zero_fibre_input import file_digest
from space_batch import decode_bitmaps, validate_batch
from space_codec import binary_rank
from verify_protocol_cover import require


def check_contractions(bitmaps, core, multiplicity):
    require(bitmaps.ndim == 2 and bitmaps.shape[1] == 32 and bitmaps.dtype == np.uint8,
            "Invalid dimension-eight bitmap array")
    fibres = np.unpackbits(bitmaps, axis=1, bitorder="little").reshape(-1, 128, 2)
    expected = np.zeros(128, dtype=np.uint8)
    expected[list(core)] = 1
    require(np.all((fibres[:, :, 0] ^ fibres[:, :, 1]) == expected)
            and np.all(np.sum(fibres[:, :, 0] & fibres[:, :, 1], axis=1) == multiplicity),
            "An emitted support contracts to the wrong core or fibre multiplicity")


def generate(inputs, native, output, selected=None):
    started = perf_counter()
    domain = json.loads((inputs / "domain.json").read_bytes())
    require(domain["schema"] == "c48-marked-contraction-domain-v1", "Different marked input domain")
    output.mkdir(parents=True, exist_ok=True)
    rows = []
    for index, row in enumerate(domain["records"]):
        if selected is not None and index not in selected:
            continue
        task = inputs / row["task"]
        require(file_digest(task) == row["sha256"], "The native input has changed")
        core = tuple(map(int, task.read_text(encoding="ascii").splitlines()[3].split()))
        directory = output / task.stem
        directory.mkdir(exist_ok=True)
        binary, summary, trace = (directory / name for name in ("marked.bin", "census.json", "orbits.tsv"))
        receipt = directory / "binding.json"
        binding = dict(input_sha256=row["sha256"], native_sha256=file_digest(native))
        if not summary.exists():
            require(not binary.exists() and not trace.exists(), "Incomplete previous marked output")
            subprocess.run([str(native), "--input", str(task), "--output", str(binary),
                            "--summary", str(summary), "--debug", str(trace)], check=True,
                           capture_output=True, text=True)
            receipt.write_text(json.dumps(binding, indent=2)+"\n", encoding="ascii")
        require(json.loads(receipt.read_bytes()) == binding, "Different native or input binding")
        census = json.loads(summary.read_bytes())
        require(census["all_checks_pass"] and census["raw_marked_pair_count"] == row["raw_pairs"]
                and census["compatible_fiber_set_count"] == row["compatible_fibres"]
                and census["lift_quotient_dimension"] == row["ell"]
                and census["core_stabilizer_order"] == row["stabiliser_order"], "Marked census mismatch")
        require(binary.stat().st_size == 40*census["direction_marked_orbit_count"], "Truncated marked output")
        full_count = full_mass = deficient_count = deficient_mass = 0
        valid_path = directory / "full_rank.bin"
        with binary.open("rb") as stream, valid_path.open("wb") as valid:
            while raw := stream.read(4096*40):
                records = np.frombuffer(raw, dtype=np.uint8).reshape(-1, 40)
                check_contractions(records[:, :32], core, row["multiplicity"])
                supports = decode_bitmaps(records[:, :32].copy().tobytes(), 8, 48)
                weights = records[:, 32:].copy().view("<u8").reshape(-1)
                require(np.all(weights > 0), "Empty marked orbit")
                if row["proper_span_requires_output_rank_filter"]:
                    full = np.array([binary_rank(int(x) ^ int(points[0]) for x in points) == 8
                                     for points in supports], dtype=bool)
                else:
                    full = np.ones(len(supports), dtype=bool)
                validate_batch(supports[full], 8, 48)
                valid.write(records[full, :32].copy().tobytes())
                full_count += int(np.sum(full))
                full_mass += int(np.sum(weights[full]))
                deficient_count += int(np.sum(~full))
                deficient_mass += int(np.sum(weights[~full]))
        require(full_mass+deficient_mass == row["raw_pairs"]
                and full_count+deficient_count == census["direction_marked_orbit_count"],
                "Marked orbit masses do not exhaust the source")
        result = dict(source_index=index, task=row["task"], input_sha256=row["sha256"],
                      marked_orbits=full_count+deficient_count, full_rank_orbits=full_count,
                      full_rank_pairs=full_mass, deficient_orbits=deficient_count,
                      deficient_pairs=deficient_mass, raw_pairs=row["raw_pairs"],
                      full_rank_binary_sha256=file_digest(valid_path),
                      marked_binary_sha256=file_digest(binary), trace_sha256=file_digest(trace))
        (directory / "validation.json").write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
        rows.append(result)
        print(json.dumps(result), flush=True)
        (output / "progress.json").write_text(json.dumps(dict(records=rows, elapsed_seconds=perf_counter()-started), indent=2)+"\n")
    result = dict(schema="c48-marked-generation-v1", status="pass", source_count=len(rows),
                  complete_source_domain=selected is None, input_domain_sha256=file_digest(inputs / "domain.json"),
                  native_sha256=file_digest(native), records=rows, elapsed_seconds=perf_counter()-started)
    (output / ("complete.json" if selected is None else "selected.json")).write_text(json.dumps(result, indent=2)+"\n", encoding="ascii")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-index", type=int, action="append")
    args = parser.parse_args()
    generate(args.inputs.resolve(), args.native.resolve(), args.output.resolve(),
             None if args.source_index is None else set(args.source_index))


if __name__ == "__main__":
    main()
