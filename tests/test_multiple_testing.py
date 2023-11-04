import numpy as np

from goatools.multiple_testing import HolmBonferroni


def test_holmbonferroni():
    corrected_pvals = HolmBonferroni(
        [0.01, 0.01, 0.03, 0.05, 0.005], a=0.05
    ).corrected_pvals
    assert np.allclose(corrected_pvals, np.array([0.04, 0.04, 0.06, 0.05, 0.025]))


def test_holmbonferroni_pval_all_1():
    corrected_pvals = HolmBonferroni([1, 1, 1, 1, 1], a=0.05).corrected_pvals
    assert np.allclose(corrected_pvals, np.array([1, 1, 1, 1, 1]))


def test_holmbonferroni_pval_mixed():
    corrected_pvals = HolmBonferroni([1, 0.01, 0.03, 0.05, 1], a=0.05).corrected_pvals
    assert np.allclose(corrected_pvals, np.array([1.0, 0.05, 0.12, 0.15, 1.0]))
