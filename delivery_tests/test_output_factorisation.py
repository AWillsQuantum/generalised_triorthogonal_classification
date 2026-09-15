"""Exact factorisation, equivariance, full-phase and publication regressions."""

from copy import deepcopy
from pathlib import Path
import json
import random
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from output_factorisation import (clifford_correction, factor_registry, independent_blocks,
    normal_form, primitive_projectors, restrict, update_catalogue, verify_presentation_update)
from protocol_checks import gate_tensor, verify_basis


class FactorisationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        source = ROOT / "provenance/pareto_frontier_before_factorisation.json"
        if not source.exists():
            source = ROOT / "data/protocols/pareto_frontier.json"
        cls.original = json.loads(source.read_bytes())
        cls.registry = factor_registry(cls.original["outputs"])
        cls.updated = update_catalogue(cls.original)

    def test_all_classes_and_scientific_fields(self):
        report = verify_presentation_update(self.original, self.updated)
        self.assertTrue(report["scientific_fields_unchanged"])
        self.assertEqual(report["outputs"], 62)
        self.assertEqual(report["protocols"], 74)

    def test_hidden_six_T(self):
        row = next(r for r in self.updated["outputs"] if r["output_id"] == "Q6_000000400a08a400")
        self.assertEqual(row["representative_gate"], "T1T2T3T4T5T6")
        self.assertEqual(row["factorised_representative_gate"], "T^(tensor 6)")
        self.assertEqual(row["factorised_representative_latex"], r"\mathsf{T}^{\otimes 6}")

    def test_familiar_products(self):
        for q, gate, label in [(5, "CS12CCZ345", "CS tensor CCZ"),
                               (6, "CCZ123CCZ456", "CCZ^(tensor 2)"),
                               (7, "T1CS23CS45CS67", "T tensor CS^(tensor 3)")]:
            form, _ = normal_form(gate_tensor(q, gate), self.registry)
            self.assertEqual(form["factorised_representative_gate"], label)

    def test_basis_invariance_for_every_class(self):
        rng = random.Random(20260915)
        for old, new in zip(self.original["outputs"], self.updated["outputs"]):
            tensor = gate_tensor(old["q"], old["representative_gate"])
            for _ in range(2):
                basis = [1 << i for i in range(tensor.q)]
                if tensor.q > 1:
                    for _ in range(12*tensor.q):
                        i, j = rng.sample(range(tensor.q), 2)
                        basis[i] ^= basis[j]
                changed = restrict(tensor, basis)
                form, witness = normal_form(changed, self.registry)
                self.assertEqual(form["factorised_representative_gate"], new["factorised_representative_gate"])
                self.assertEqual(form["representative_gate"], new["representative_gate"])
                verify_basis(changed, gate_tensor(tensor.q, form["representative_gate"]), witness)

    def test_intrinsic_required(self):
        with self.assertRaisesRegex(ValueError, "intrinsic"):
            independent_blocks(gate_tensor(2, "T1"))

    def test_same_diagonal_does_not_mean_same_tensor(self):
        cs, ccz = gate_tensor(2, "CS12"), gate_tensor(3, "CCZ123")
        self.assertTrue(all(cs.value(v,v,v) == 0 for v in range(4)))
        self.assertTrue(all(ccz.value(v,v,v) == 0 for v in range(8)))
        self.assertEqual(len(primitive_projectors(cs)), 1)
        self.assertEqual(len(primitive_projectors(ccz)), 1)

    def test_unknown_factor_is_not_silently_named(self):
        with self.assertRaisesRegex(ValueError, "absent"):
            normal_form(gate_tensor(2, "CS12"), self.registry[:1])

    def test_tampering_rejected(self):
        changes = [lambda p: p["protocols"][0].__setitem__("n", 99),
                   lambda p: p["outputs"][1].__setitem__("factorised_representative_gate", "CS"),
                   lambda p: p["outputs"][1]["gate_basis_in_tensor_coordinates"].__setitem__(0, 0),
                   lambda p: p["outputs"][1]["factorisation_certificate"].__setitem__("diagonal_clifford_correction", {})]
        for change in changes:
            candidate = deepcopy(self.updated)
            change(candidate)
            with self.assertRaises(ValueError):
                verify_presentation_update(self.original, candidate)


if __name__ == "__main__":
    unittest.main()
