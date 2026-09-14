import io
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))

from protocol_domain import (POINTING, read_header, read_native_header, read_native_records,
    read_pointings, source_points, write_base_native, write_header, write_native)
from verify_small_source_closure import all_pointings, check_zero_census
from verify_protocol_quotient_partition import check_distribution
from verify_through48_sector import nondominated_metrics


class ProtocolDomainTests(unittest.TestCase):
    def setUp(self):
        self.sources = [dict(index=0, space=dict(c=16, m=4, index=0))]
        self.supports = {(16, 4, 0): tuple(range(16))}

    def test_compact_record_round_trip(self):
        stream = io.BytesIO()
        write_header(stream, 1, 2, 16)
        stream.write(POINTING.pack(0, 0, 3))
        stream.write(POINTING.pack(0, 1, 3))
        stream.seek(0)
        header = read_header(stream)
        self.assertEqual(list(read_pointings(stream, header)), [(0, 0, 3), (0, 1, 3)])

    def test_truncation_and_trailing_bytes(self):
        stream = io.BytesIO()
        write_header(stream, 1, 1, 16)
        stream.write(POINTING.pack(0, 0, 3))
        raw = stream.getvalue()
        for bad in (raw[:-1], raw+b"x"):
            handle = io.BytesIO(bad)
            header = read_header(handle)
            with self.assertRaises(ValueError):
                list(read_pointings(handle, header))

    def test_bad_source_and_parity(self):
        for record in ((1, 0, 0), (0, 5, 0), (0, 1, -1)):
            stream = io.BytesIO()
            write_header(stream, 1, 1, 16)
            stream.write(POINTING.pack(*record))
            stream.seek(0)
            header = read_header(stream)
            with self.assertRaises(ValueError):
                list(read_pointings(stream, header))

    def test_native_expansion(self):
        stream = io.BytesIO()
        write_native(stream, dict(source_count=1, count=2, maximum_protocol_length=16),
                     self.sources, self.supports, [(0, 0, 0), (0, 1, 0)])
        stream.seek(0)
        header = read_native_header(stream)
        rows = list(read_native_records(stream, header))
        self.assertEqual(rows[0]["points"], tuple(range(16)))
        self.assertEqual(rows[1]["points"], tuple(range(1, 16)))
        self.assertEqual(header["flags"], 1)

    def test_coordinate_frame_and_base_expansion(self):
        source = dict(index=0, space=dict(c=4, m=3, index=0), coordinate_frame=[1, 4, 2, 1])
        supports = {(4, 3, 0): (0, 1, 2, 4)}
        self.assertEqual(source_points(source, supports), (0, 1, 3, 5))
        stream = io.BytesIO()
        write_base_native(stream, [source], supports, 4)
        stream.seek(0)
        header = read_native_header(stream)
        self.assertEqual(header["flags"], 3)
        self.assertEqual(next(read_native_records(stream, header))["points"], (0, 1, 3, 5))
        source["coordinate_frame"] = [0, 1, 1, 2]
        with self.assertRaises(ValueError):
            source_points(source, supports)

    def test_all_pointings_remove_identical_supports(self):
        cases = all_pointings(tuple(range(16)), 4, 16)
        self.assertEqual(len(cases), 3)
        self.assertEqual(sorted(len(row["points"]) for row in cases), [15, 16, 16])
        self.assertEqual({row["ambient_dimension"] for row in cases}, {4, 5})

    def test_positive_result_cannot_certify_absence(self):
        result = dict(status="complete", matrix_scope="full_projective", logical_qubits=5,
            minimum_distance=3, enumeration_mode="raw_isotropic_subspaces",
            support_range=dict(start=0, count=1, end_exclusive=1), pareto_protocols=[],
            statistics=dict(supports_processed=1, quotient_dimension_filtered_supports=0,
                marked_orbit_enumerated_supports=0, isotropic_subspace_statistics_exact=True,
                radical_statistics_exact=True, nondegenerate_subspaces=1,
                radical_dimension_counts={}, isotropic_subspaces=1))
        with self.assertRaises(ValueError):
            check_zero_census(result, 5, 1, 3)

    def test_exact_distance_pareto_partition(self):
        values = {("A", 3, 15, 5), ("A", 3, 31, 6), ("A", 4, 31, 6),
                  ("B", 3, 31, 6), ("A", 3, 14, 7)}
        self.assertEqual(nondominated_metrics(values), values-{("A", 3, 31, 6)})

    def test_quotient_partition_rejects_missing_large_case(self):
        domain = dict(pointing_count=3, maximum_protocol_length=48)
        census = dict(small_quotient_interval=[0, 11], small_quotient_supports=2,
                      quotient_selected_counts=[[12, 1]])
        result = dict(manifest_supports=3, records_processed=3, record_start=0,
            logical_dimension=64, minimum_distance=3, manifest_length_filter="protocol",
            length_limit=48, eligible_supports=0, isotropic_subspace_count=0,
            label_spaces_by_dimension=[[3, 2], [12, 1]])
        self.assertEqual(check_distribution(domain, census, result), [[3, 2], [12, 1]])
        result["label_spaces_by_dimension"] = [[3, 1], [12, 1], [15, 1]]
        with self.assertRaises(ValueError):
            check_distribution(domain, census, result)


if __name__ == "__main__":
    unittest.main()
