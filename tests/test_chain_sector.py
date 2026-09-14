from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from replay_chain_sector import cases_for_interval


class ChainDomainTests(unittest.TestCase):
    def test_interval_keeps_exact_input_identity(self):
        data = dict(q5_input_support_class_indices=[0], q6_candidates=[dict(support_class_index=0)])
        domain = dict(classes=[dict(ambient_dimension=1, support_length=2, key_words=["102", "2", "2"])])
        for q in (5, 6):
            cases, total = cases_for_interval(data, domain, q, 0, 1)
            self.assertEqual(total, 1)
            self.assertEqual(cases[0]["points"], [0, 1])
            self.assertEqual(cases[0]["support_class_index"], 0)
        for start, count in ((-1, 1), (1, 1), (0, 0)):
            with self.assertRaises(ValueError):
                cases_for_interval(data, domain, 5, start, count)


if __name__ == "__main__":
    unittest.main()
