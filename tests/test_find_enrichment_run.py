#!/usr/bin/env python
"""Test running of enrichment example in run.sh.

  $ python find_enrichment.py --pval=0.05 --indent data/study data/population data/association.
"""

from __future__ import print_function

__copyright__ = "Copyright (C) 2010-2018, DV Klopfenstein, H Tang. All rights reserved."

import collections as cx
from goatools.base import get_godag
from goatools.cli.find_enrichment import GoeaCliFnc
from goatools.test_data.cli.find_enrichment_dflts import ArgsDict


# This test is not on Travis because non-code file, data/study, is not found by Travis-CI
# This test can and should be run before any pull requests using 'make test'
def test_find_enrichment():
    """Recreate run in run.sh."""
    # Set params
    objtest = ArgsDict()
    get_godag(objtest.namespace['obo'], loading_bar=None)
    objtest.namespace['indent'] = True
    args = objtest.ntobj(**objtest.namespace)
    # Run test
    objcli = GoeaCliFnc(args)

    # Check results
    ## expected_cnts = {'fdr_bh': 17, 'sidak': 5, 'holm': 5, 'bonferroni': 5}
    expected_cnts = {'fdr_bh': 19, 'sidak': 9, 'holm': 9, 'bonferroni': 9}
    _chk_results(objcli.results_all, expected_cnts, objcli)
    print("TEST PASSED")


def _chk_results(results, expected_cnts, objcli):
    """Check results."""
    ctr = cx.Counter()
    alpha = objcli.args.alpha
    for ntres in results:
        for method in objcli.methods:
            ctr[method] += getattr(ntres, "p_{METHOD}".format(METHOD=method)) < alpha
    for method, num_sig in ctr.most_common():
        assert num_sig == expected_cnts[method], '{EXP} {ACT}'.format(
            EXP=expected_cnts, ACT=ctr.most_common())



if __name__ == '__main__':
    test_find_enrichment()

# Copyright (C) 2010-2018, DV Klopfenstein, H Tang. All rights reserved.
