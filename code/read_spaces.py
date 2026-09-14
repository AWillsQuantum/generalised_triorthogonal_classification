"""Stream compact space representatives, expanded matrices or polynomials."""

import argparse
from contextlib import contextmanager
import json
from pathlib import Path

from space_codec import (generator_rows, indicator_anf, polynomial_text,
                         read_header, read_records, validate_unital_support)


@contextmanager
def open_shard(path):
    path = Path(path)
    with path.open("rb") as raw:
        if path.suffix == ".zst":
            try:
                import zstandard
            except ImportError as error:
                raise RuntimeError("Reading .zst shards requires the zstandard package") from error
            with zstandard.ZstdDecompressor().stream_reader(raw) as decoded:
                yield decoded
        else:
            yield raw


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("shard", type=Path)
    parser.add_argument("--limit", type=int, default=10)
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--matrix", action="store_true")
    parser.add_argument("--polynomial", action="store_true")
    parser.add_argument("--validate", action="store_true")
    args = parser.parse_args()
    if args.limit < 0:
        parser.error("--limit must be nonnegative")
    with open_shard(args.shard) as stream:
        header = read_header(stream)
        records = read_records(stream, header)
        count = header.count if args.all else min(args.limit, header.count)
        for offset in range(count):
            points = next(records)
            result = dict(index=header.first_index + offset, c=header.c, m=header.m,
                          support=list(points))
            if args.validate:
                validate_unital_support(points, header.m)
                result["valid_unital_space"] = True
            if args.matrix:
                result["generator_matrix_rows"] = generator_rows(points, header.m)
            if args.polynomial:
                terms = indicator_anf(points, header.m)
                result["polynomial"] = polynomial_text(terms, header.m)
            print(json.dumps(result))
        if count == header.count:
            try:
                next(records)
            except StopIteration:
                pass


if __name__ == "__main__":
    main()
