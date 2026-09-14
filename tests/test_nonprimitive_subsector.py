from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from verify_nonprimitive_subsector import orbit_set


class MarkedOrbitTests(unittest.TestCase):
    def test_exact_words_and_positive_orbit_sizes(self):
        self.assertEqual(orbit_set([dict(marked_key_words=["a", "11"], orbit_size=7)]), {(10, 17): 7})
        for rows in ([dict(marked_key_words=["a"], orbit_size=0)],
                     [dict(marked_key_words=["a"], orbit_size=1),
                      dict(marked_key_words=["0a"], orbit_size=2)]):
            with self.subTest(rows=rows), self.assertRaises(ValueError):
                orbit_set(rows)


if __name__ == "__main__":
    unittest.main()
