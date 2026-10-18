import unittest

import numpy as np

from scripts.evaluate_external_candidate import bootstrap_intervals


class EvaluateExternalCandidateTest(unittest.TestCase):
    def test_stratified_bootstrap_returns_ordered_intervals(self) -> None:
        intervals = bootstrap_intervals(
            np.array([0, 0, 0, 1, 1, 1]),
            np.array([0.1, 0.2, 0.4, 0.6, 0.8, 0.9]),
            seed=2027,
            replicates=100,
        )
        self.assertEqual(set(intervals), {"roc_auc", "average_precision"})
        for interval in intervals.values():
            self.assertLessEqual(interval["lower_95"], interval["upper_95"])
            self.assertGreaterEqual(interval["lower_95"], 0.0)
            self.assertLessEqual(interval["upper_95"], 1.0)


if __name__ == "__main__":
    unittest.main()
