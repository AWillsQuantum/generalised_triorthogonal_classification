"""Recompute finite tensor censuses and the q3-to-q5 covering lemma."""

import argparse
from collections import Counter
import hashlib
from itertools import combinations
import json
from math import prod
from pathlib import Path
import subprocess

from protocol_checks import LogicalTensor, find_equivalence_basis, gate_tensor


def standard_tensor(q, signature):
    """Signature order: singles, lexicographic pairs, lexicographic triples."""
    monomials = [indices for degree in (1, 2, 3) for indices in combinations(range(q), degree)]
    if not 0 <= signature < 1 << len(monomials):
        raise ValueError("Tensor signature outside its declared dimension")
    columns = []
    for bit, indices in enumerate(monomials):
        if signature & (1 << bit):
            for subset in range(1, 1 << len(indices)):
                columns.append(sum(1 << i for j, i in enumerate(indices) if subset & (1 << j)))
    return LogicalTensor(tuple(sum(((column >> i) & 1) << j for j, column in enumerate(columns))
                               for i in range(q)))


def subspace_bases(ambient, dimension):
    """Unique binary reduced-echelon bases, with least-significant pivots."""
    if not 0 <= dimension <= ambient:
        raise ValueError("Invalid subspace dimension")
    for pivots in combinations(range(ambient), dimension):
        free = [(row, column) for row, pivot in enumerate(pivots)
                for column in range(pivot + 1, ambient) if column not in pivots]
        for state in range(1 << len(free)):
            rows = [1 << pivot for pivot in pivots]
            for bit, (row, column) in enumerate(free):
                rows[row] |= ((state >> bit) & 1) << column
            yield tuple(rows)


def restriction(tensor, basis):
    return LogicalTensor(tuple(tensor.words[v] for v in basis))


def verify_q5_cover(census):
    if (census["status"] != "pass" or census["tensor_count"] != 1 << 25
            or census["generated_linear_group_order"] != prod((1 << 5) - (1 << i) for i in range(5))
            or census["direct_canonical_checks"] != 88
            or census["direct_canonical_mismatches"] != 0):
        raise ValueError("The exhaustive five-dimensional tensor census failed")
    rows = census["orbits"]
    if len(rows) != 88 or sum(row["orbit_size"] for row in rows) != 1 << 25:
        raise ValueError("Incorrect five-dimensional tensor orbit mass")
    if len({tuple(row["canonical_key_words"]) for row in rows}) != len(rows):
        raise ValueError("Duplicate five-dimensional tensor key")
    if any(row["orbit_size"] * row["stabilizer_size"] != census["generated_linear_group_order"] for row in rows):
        raise ValueError("Orbit-stabiliser identity failed")
    q4_bases, q3_bases = tuple(subspace_bases(5, 4)), tuple(subspace_bases(4, 3))
    if len(q4_bases) != 31 or len(q3_bases) != 15:
        raise ValueError("Incorrect hyperplane counts")
    targets = {name: gate_tensor(3, name) for name in ("CS12CS13", "T1CCZ123")}
    classes = []
    for row in rows:
        if row["radical_dimension"]:
            continue
        tensor = standard_tensor(5, int(row["canonical_standard_signature"], 16))
        if not tensor.intrinsic():
            raise ValueError("Advertised intrinsic tensor is degenerate")
        hyperplanes = 0
        witnesses = {}
        counts = Counter()
        for q4_basis in q4_bases:
            q4_tensor = restriction(tensor, q4_basis)
            if not q4_tensor.intrinsic():
                continue
            hyperplanes += 1
            for local_q3 in q3_bases:
                q3_tensor = restriction(q4_tensor, local_q3)
                for name, target in targets.items():
                    basis = find_equivalence_basis(q3_tensor, target)
                    if basis is None:
                        continue
                    counts[name] += 1
                    if name not in witnesses:
                        actual = [sum(1 << i for i in range(5)
                                      if sum(((q4_basis[j] >> i) & 1) for j in range(4) if v & (1 << j)) & 1)
                                  for v in local_q3]
                        witnesses[name] = dict(q4_basis=q4_basis, q3_basis=actual,
                                               target_basis_in_q3_coordinates=basis)
        if hyperplanes and not witnesses:
            raise ValueError("A nonprimitive output has no target chain")
        classes.append(dict(tensor_key_words=row["canonical_key_words"],
                            standard_signature=row["canonical_standard_signature"],
                            nondegenerate_hyperplanes=hyperplanes,
                            target_chain_counts=dict(counts), witnesses=witnesses))
    if len(classes) != 65 or sum(not row["nondegenerate_hyperplanes"] for row in classes) != 2:
        raise ValueError("Incorrect intrinsic or primitive tensor count")
    return dict(schema="triorthogonal-q5-target-chain-cover-v1", status="pass",
                intrinsic_classes=65, nonprimitive_classes=63, primitive_classes=2,
                q4_hyperplanes_per_q5=31, q3_hyperplanes_per_q4=15,
                target_gates=list(targets), classes=classes,
                protocol_enumeration_completeness_established=False)


def verify_alternating_mass(result, q, expected_orbits):
    order = prod((1 << q) - (1 << i) for i in range(q))
    rows = result["alternating_orbits"]
    if (result["general_linear_group_order"] != order or len(rows) != expected_orbits
            or sum(row["orbit_size"] for row in rows) != 1 << (q * (q-1) * (q-2) // 6)
            or any(row["orbit_size"] * row["stabilizer_order"] != order for row in rows)):
        raise ValueError("Alternating tensor orbit mass does not replay")
    if result["rank_one_pairs_tested"] != expected_orbits * ((1 << q) - 1):
        raise ValueError("Incomplete rank-one pair coverage")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parents[1] / "certificates/tensors")
    args = parser.parse_args()
    results = {}
    for name, command in (
        ("q5_orbits", ["five-qubit-tensor-orbits"]),
        ("q7_primitive", ["q7-primitive-tensor-authority", "--workers", "1"]),
        ("q8_primitive", ["q8-primitive-tensor-authority-screen"]),
    ):
        process = subprocess.run([str(args.native.resolve()), *command], capture_output=True, text=True, check=True)
        result = json.loads(process.stdout)
        result.pop("elapsed_seconds", None)
        result.pop("workers", None)
        results[name] = result
    results["q5_chain_cover"] = verify_q5_cover(results["q5_orbits"])
    verify_alternating_mass(results["q7_primitive"], 7, 12)
    verify_alternating_mass(results["q8_primitive"], 8, 32)
    if (results["q7_primitive"]["primitive_tensor_orbits"] != 2
            or results["q7_primitive"]["rank_one_primitive_pairs"] != 1
            or results["q8_primitive"]["primitive_candidates"]
            or results["q8_primitive"]["rank_one_primitive_pairs"]):
        raise ValueError("Primitive tensor check failed")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    hashes = {}
    for name, result in results.items():
        payload = (json.dumps(result, indent=2) + "\n").encode("ascii")
        (args.output_dir / (name + ".json")).write_bytes(payload)
        hashes[name] = hashlib.sha256(payload).hexdigest()
    print(json.dumps(dict(status="pass", certificates=hashes)))


if __name__ == "__main__":
    main()
