import copy
import io
from pathlib import Path
import struct
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from compact_to_native import MANIFEST_HEADER
from verify_protocol_case_census import check_result, normalized_keys, write_pointed_cases


class FiniteCensusTests(unittest.TestCase):
    def data(self):
        return dict(logical_qubits=1, minimum_distance=3, maximum_protocol_length=54,
                    cases=[dict(index=0, ambient_dimension=4, points=list(range(16)),
                                expected=dict(isotropic_subspaces=1, nondegenerate_subspaces=1,
                                              output_keys=[["0x1"]]))])

    def result(self):
        return dict(status="complete", matrix_scope="full_projective", logical_qubits=1,
                    minimum_distance=3, enumeration_mode="raw_isotropic_subspaces",
                    support_range=dict(start=0, count=1, end_exclusive=1),
                    statistics=dict(supports_processed=1, isotropic_subspace_statistics_exact=True,
                                    radical_statistics_exact=True, isotropic_subspaces=1,
                                    nondegenerate_subspaces=1),
                    positive_supports=[dict(manifest_support_index=0, original_record_index=0,
                                            output_keys=[dict(logical_qubits=1, canonical_key_words=[1])])])

    def test_writer_flags_and_input_validation(self):
        stream = io.BytesIO()
        write_pointed_cases(stream, self.data())
        self.assertEqual(MANIFEST_HEADER.unpack(stream.getvalue()[:MANIFEST_HEADER.size]),
                         (b"UTSPTS1\0", 1, 1, 1, 54, 1))
        for field, value in (("index", 1), ("points", [0, 1, 2, 3]),
                             ("points", list(range(15))), ("points", [0, 0, *range(2, 16)])):
            data = self.data()
            data["cases"][0][field] = value
            with self.assertRaises(ValueError):
                write_pointed_cases(io.BytesIO(), data)

    def test_typed_output_profile_and_exact_count_checks(self):
        self.assertEqual(normalized_keys([[1], ["0x01"]]), {(1,)})
        for keys in ([[]], [[True]], [[-1]], [[1 << 64]]):
            with self.assertRaises(ValueError):
                normalized_keys(keys)
        self.assertEqual(check_result(self.data(), self.result())["status"], "pass")
        for field, value in (("status", "incomplete"), ("logical_qubits", 2),
                             ("minimum_distance", 2), ("positive_supports", [])):
            result = self.result()
            result[field] = value
            with self.assertRaises(ValueError):
                check_result(self.data(), result)
        result = self.result()
        result["statistics"]["isotropic_subspaces"] = 2
        with self.assertRaises(ValueError):
            check_result(self.data(), result)
        result = self.result()
        result["positive_supports"].append(copy.deepcopy(result["positive_supports"][0]))
        with self.assertRaises(ValueError):
            check_result(self.data(), result)


if __name__ == "__main__":
    unittest.main()
