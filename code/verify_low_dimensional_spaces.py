"""Check every space through affine dimension seven against complete seed orbits."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from affine_graph import canonicalise_many, rank
from read_spaces import open_shard
from space_codec import read_header, read_records, validate_unital_support
from verify_rm37 import certify, polynomial_support


def lower_degree_seeds(m):
    """Affine normal forms for RM(m-4,m), for 4 <= m <= 6."""
    if m not in (4, 5, 6):
        raise ValueError("Only constant, affine and quadratic sectors are included")
    yield ()
    yield (0,)
    if m >= 5:
        yield (1,)
    if m == 6:
        for half_rank in range(1, m // 2 + 1):
            quadratic = tuple(3 << (2*i) for i in range(half_rank))
            yield quadratic
            yield (0, *quadratic)
            if 2 * half_rank < m:
                yield (*quadratic, 1 << (2 * half_rank))


def expected_lower_spaces(command, maximum_length):
    inputs = []
    for m in range(4, 7):
        for monomials in lower_degree_seeds(m):
            mask = polynomial_support(monomials, m, m-4)
            points = tuple(x for x in range(1 << m) if mask >> x & 1)
            if (points and len(points) <= maximum_length
                    and rank(x ^ points[0] for x in points) == m):
                inputs.append((m, mask))
    forms = canonicalise_many(command, inputs)
    return {(mask.bit_count(), m, int(form["canonical_mask"], 16))
            for (m, mask), form in zip(inputs, forms, strict=True)}


def verify(root, command, maximum_length=54):
    if not 16 <= maximum_length <= 54:
        raise ValueError("Catalogue scope is lengths 16 through 54")
    root = root.resolve()
    seed_bytes = (root / "data/seeds/rm37.json").read_bytes()
    seed_summary, seed_forms = certify(json.loads(seed_bytes), command)
    expected = expected_lower_spaces(command, maximum_length)
    expected.update((row["weight"], 7, int(row["canonical_mask"], 16))
                    for row in seed_forms if row["affine_rank"] == 7
                    and 0 < row["weight"] <= maximum_length)
    rows = []
    bindings = {}
    for certificate_path in sorted((root / "certificates/support_encoding").glob("*.json")):
        certificate = json.loads(certificate_path.read_bytes())
        if certificate["m"] > 7 or certificate["c"] > maximum_length:
            continue
        for shard in certificate["shards"]:
            path = (root / shard["path"]).resolve()
            if not path.is_relative_to(root) or shard["path"] in bindings:
                raise ValueError("Repeated or external catalogue shard")
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            if digest != shard["sha256"]:
                raise ValueError("Catalogue shard digest mismatch")
            bindings[shard["path"]] = digest
            with open_shard(path) as stream:
                header = read_header(stream)
                if (header.c, header.m, header.count, header.first_index) != (
                        certificate["c"], certificate["m"], shard["count"], shard["first_index"]):
                    raise ValueError("Catalogue shard scope mismatch")
                for offset, support in enumerate(read_records(stream, header)):
                    validate_unital_support(support, header.m)
                    rows.append((header.c, header.m, header.first_index + offset,
                                 sum(1 << x for x in support)))
    forms = canonicalise_many(command, [(m, mask) for _, m, _, mask in rows])
    keys = [(c, m, int(form["canonical_mask"], 16))
            for (c, m, _, _), form in zip(rows, forms, strict=True)]
    if len(set(keys)) != len(keys):
        raise ValueError("The catalogue repeats an affine class")
    if set(keys) != expected:
        raise ValueError(f"Catalogue/seed mismatch: {len(expected-set(keys))} missing, "
                         f"{len(set(keys)-expected)} extra")
    counts = Counter((c, m) for c, m, _, _ in rows)
    return dict(schema="low-dimensional-space-completeness-v1", status="pass",
                maximum_length=maximum_length, maximum_affine_dimension=7,
                complete_and_pairwise_inequivalent=True, classes=len(rows),
                seed_orbit_mass=seed_summary["exact_orbit_mass"],
                seeds_sha256=hashlib.sha256(seed_bytes).hexdigest(),
                scope="nonempty projective unital spaces; constant, affine, quadratic normal forms and full cubic orbit exhaustion",
                counts=[dict(c=c, m=m, classes=n) for (c, m), n in sorted(counts.items())],
                catalogue_sha256=bindings,
                canonical_maps=[dict(form, c=c, index=index)
                                for (c, m, index, _), form in zip(rows, forms, strict=True)])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--maximum-length", type=int, default=54)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = verify(args.root, [args.native], args.maximum_length)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k not in ("canonical_maps", "catalogue_sha256")}, indent=2))
