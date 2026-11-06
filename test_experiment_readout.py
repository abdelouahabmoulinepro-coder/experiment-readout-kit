import csv
import tempfile
import unittest
from pathlib import Path

from experiment_readout import analyze, load_assignments, wilson, write_report


def make_rows(cn, cc, tn, tc):
    return ([{"unit_id": f"c{i}", "variant": "control", "converted": int(i < cc)} for i in range(cn)] +
            [{"unit_id": f"t{i}", "variant": "treatment", "converted": int(i < tc)} for i in range(tn)])


class ReadoutTests(unittest.TestCase):
    def test_equal_rates_are_inconclusive(self):
        result = analyze(make_rows(100, 20, 100, 20))
        self.assertEqual(result["decision"], "inconclusive")
        self.assertEqual(result["absolute_lift"], 0)

    def test_clear_positive_requires_interval_above_zero(self):
        result = analyze(make_rows(1000, 100, 1000, 180))
        self.assertEqual(result["decision"], "positive")
        self.assertGreater(result["difference_interval"][0], 0)

    def test_negative_effect(self):
        result = analyze(make_rows(1000, 180, 1000, 100))
        self.assertEqual(result["decision"], "negative")

    def test_srm_invalid_blocks_positive_effect(self):
        result = analyze(make_rows(1800, 180, 200, 60))
        self.assertEqual(result["decision"], "invalid_srm")
        self.assertFalse(result["sample_ratio_check"]["passed"])

    def test_configured_allocation_ratio(self):
        result = analyze(make_rows(700, 70, 300, 36), expected_control_share=.7)
        self.assertTrue(result["sample_ratio_check"]["passed"])

    def test_zero_conversions_has_finite_interval(self):
        result = analyze(make_rows(50, 0, 50, 0))
        self.assertEqual(result["control"]["conversion_rate"], 0)
        self.assertGreater(result["control"]["wilson_interval"][1], 0)
        self.assertIsNone(result["relative_lift"])

    def test_wilson_interval_is_bounded(self):
        low, high = wilson(0, 20, 1.96)
        self.assertEqual(low, 0)
        self.assertGreater(high, 0)
        self.assertLessEqual(high, 1)

    def test_duplicate_unit_rejected_on_csv_load(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.csv"
            path.write_text("unit_id,variant,converted\nu1,control,0\nu1,treatment,1\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "dupliquée"):
                load_assignments(path, "control", "treatment")

    def test_invalid_conversion_and_unknown_arm(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.csv"
            path.write_text("unit_id,variant,converted\nu1,control,2\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "0 ou 1"):
                load_assignments(path, "control", "treatment")
            path.write_text("unit_id,variant,converted\nu1,other,1\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "inconnue"):
                load_assignments(path, "control", "treatment")

    def test_empty_arm_and_bad_config_rejected(self):
        with self.assertRaisesRegex(ValueError, "au moins une unité"):
            analyze(make_rows(2, 1, 0, 0))
        with self.assertRaises(ValueError):
            analyze(make_rows(2, 1, 2, 1), expected_control_share=1)
        with self.assertRaises(ValueError):
            analyze(make_rows(2, 1, 2, 1), confidence=.5)

    def test_report_artifacts(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            write_report(analyze(make_rows(100, 10, 100, 12)), out)
            self.assertTrue((out / "readout.json").exists())
            self.assertIn("Limites", (out / "report.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
