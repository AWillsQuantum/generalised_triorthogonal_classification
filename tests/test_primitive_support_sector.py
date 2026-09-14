import copy
import hashlib
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from verify_primitive_support_sector import support_key, validate_zero_result, verify


class PrimitiveSupportTests(unittest.TestCase):
    def example(self):
        points = list(range(15, 0, -1))
        rows = [sum((point >> i & 1) << j for j, point in enumerate(points)) for i in range(4)]
        words = [format(15 | 4 << 8, "x"), *map(lambda x: format(x, "x"), rows*2)]
        key = hashlib.sha256(json.dumps([15, 4, words], separators=(",", ":")).encode("ascii")).hexdigest()
        return dict(cases=[dict(index=0, ambient_dimension=4, support_length=15, points=points,
            key_words=words, canonical_key_sha256=key, automorphism_group_order=20160,
            result=dict(quotient_dimension=1, subspace_orbits=0, weighted_subspace_count=0))],
            normalisation_inputs=[dict(index=0, case_index=0, ambient_dimension=4,
                points=sorted(points), family="larger_quotient_supports", deleted_zero_column=False)])

    def test_ordered_key_and_unordered_support_geometry(self):
        data = self.example()
        result = verify(data)
        self.assertEqual(result["distance_filtration"], [dict(d3=1, d4=0, count=1)])
        self.assertFalse(result["is_global_completeness_certificate"])
        self.assertFalse(result["normalisation_equivalences_freshly_recomputed"])

    def test_bad_key_or_domain_rejected(self):
        for mutation in ("key", "index", "point", "dimension", "count", "mapping"):
            data = self.example()
            row = data["cases"][0]
            if mutation == "key":
                row["key_words"][-1] = "0"
            elif mutation == "index":
                row["index"] = 1
            elif mutation == "point":
                row["points"][0] = row["points"][1]
            elif mutation == "dimension":
                row["result"]["quotient_dimension"] = 2
            elif mutation == "count":
                row["result"]["weighted_subspace_count"] = 1
            else:
                data["normalisation_inputs"] = []
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                verify(data)

    def test_only_complete_target_zeros_accepted(self):
        result = dict(logical_dimension=5, target_dimension=5, target_nondegenerate_seed=3,
            target_signature_count=2, used_complete_primitive_target_recognizer=True,
            sector="all", subspace_orbits=0, weighted_subspace_count=0)
        validate_zero_result(result)
        for name, value in (("sector", "even"), ("weighted_subspace_count", 1),
                            ("target_signature_count", 1), ("representatives", [[1, 2]]),
                            ("used_complete_primitive_target_recognizer", False)):
            bad = copy.deepcopy(result)
            bad[name] = value
            with self.subTest(name=name), self.assertRaises(ValueError):
                validate_zero_result(bad)


if __name__ == "__main__":
    unittest.main()
