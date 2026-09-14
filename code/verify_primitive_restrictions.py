"""Certify explicit five-dimensional restrictions of the two primitive q7 types."""

import argparse
import hashlib
import json
from pathlib import Path

from protocol_checks import find_equivalence_basis, gate_tensor
from verify_tensor_authorities import restriction, standard_tensor, subspace_bases


def verify(census):
    if census["status"] != "pass" or census["tensor_count"] != 1 << 25:
        raise ValueError("A complete five-dimensional tensor census is required")
    certificates = []
    for epsilon in (0, 1):
        gate = ("T1" if epsilon else "") + "CCZ123CCZ145CCZ167"
        tensor = gate_tensor(7, gate)
        if not tensor.intrinsic() or any(restriction(tensor, basis).intrinsic()
                                         for basis in subspace_bases(7, 6)):
            raise ValueError("The stated seven-dimensional representative is not primitive")
        basis = (1, 2, 4, 8, 16)
        smaller = restriction(tensor, basis)
        if not smaller.intrinsic() or any(restriction(smaller, h).intrinsic()
                                          for h in subspace_bases(5, 4)):
            raise ValueError("The explicit restriction is not primitive and intrinsic")
        matches = []
        for row in census["orbits"]:
            if row["radical_dimension"]:
                continue
            canonical = standard_tensor(5, int(row["canonical_standard_signature"], 16))
            equivalence = find_equivalence_basis(smaller, canonical)
            if equivalence is not None:
                matches.append(dict(output_key_words=row["canonical_key_words"],
                                    canonical_standard_signature=row["canonical_standard_signature"],
                                    canonical_basis_in_restriction_coordinates=equivalence))
        if len(matches) != 1:
            raise ValueError("The restriction does not identify a unique output class")
        certificates.append(dict(epsilon=epsilon, q7_representative_gate=gate,
                                 restriction_basis=basis, **matches[0]))
    return dict(schema="primitive-q7-q5-restriction-v1", status="pass",
                q7_types=2, q7_hyperplanes_per_type=127, q5_hyperplanes_per_type=31,
                restrictions=certificates,
                necessary_rule="A realised primitive q7 output requires the corresponding q5 output on the same support and distance floor",
                is_global_completeness_certificate=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    raw = (args.root / "certificates/tensors/q5_orbits.json").read_bytes()
    result = verify(json.loads(raw))
    result["q5_census_sha256"] = hashlib.sha256(raw).hexdigest()
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="ascii")
    print(json.dumps(result, indent=2))
