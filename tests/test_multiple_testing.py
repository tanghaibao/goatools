import numpy as np
import pytest
from statsmodels.stats.multitest import multipletests

from goatools.multiple_testing import Bonferroni, HolmBonferroni, Sidak


def test_holmbonferroni():
    corrected_pvals = HolmBonferroni(
        [0.01, 0.01, 0.03, 0.05, 0.005], a=0.05
    ).corrected_pvals
    # 0.05 is the largest, scaled by 1; step-down raises it to the 0.06 before it
    assert np.allclose(corrected_pvals, np.array([0.04, 0.04, 0.06, 0.06, 0.025]))


def test_holmbonferroni_pval_all_1():
    corrected_pvals = HolmBonferroni([1, 1, 1, 1, 1], a=0.05).corrected_pvals
    assert np.allclose(corrected_pvals, np.array([1, 1, 1, 1, 1]))


def test_holmbonferroni_pval_mixed():
    corrected_pvals = HolmBonferroni([1, 0.01, 0.03, 0.05, 1], a=0.05).corrected_pvals
    assert np.allclose(corrected_pvals, np.array([1.0, 0.05, 0.12, 0.15, 1.0]))


@pytest.mark.parametrize("cls, method", [(HolmBonferroni, "holm"), (Sidak, "sidak"), (Bonferroni, "bonferroni")])
def test_local_corrections_match_statsmodels(cls, method):
    """Local corrections give the same adjusted p-values as statsmodels, ties included."""
    rng = np.random.default_rng(0)
    cases = [[0.03, 0.04], [1, 1, 1], [0.5]] + [
        list(np.round(rng.uniform(0, 0.2, n), 3)) for n in (5, 50, 500)]  # rounding gives ties
    for pvals in cases:
        expected = multipletests(pvals, 0.05, method)[1]
        assert np.allclose(cls(pvals, 0.05).corrected_pvals, expected), pvals
