#!/usr/bin/env python
"""Test running of enrichment example in run.sh.

  $ python find_enrichment.py --pval=0.05 --indent data/study data/population data/association.
"""

from __future__ import print_function

__copyright__ = "Copyright (C) 2010-2018, DV Klopfenstein, H Tang. All rights reserved."

import os
import collections as cx
from goatools.cli.find_enrichment import rd_files
from goatools.cli.find_enrichment import get_objgoea


REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def test_find_enrichment():
    """Recreate run in run.sh."""
    # Set params
    ntobj = cx.namedtuple("args_namespc", ("filenames obo "
                                           "pval alpha pvalcalc method no_propagate_counts "
                                           "compare ratio "
                                           "outfile indent min_overlap "))
    filenames = ['data/study', 'data/population', 'data/association']
    methods = ['bonferroni', 'sidak', 'holm', 'fdr_bh']
    alpha = 0.05
    args = ntobj(
        filenames=[os.path.join(REPO, f) for f in filenames],
        obo='go-basic.obo',
        pval=0.05,
        alpha=alpha,
        pvalcalc='fisher',
        method=",".join(methods),
        no_propagate_counts=False,
        compare=False,
        ratio=None,
        outfile=None,
        indent=True,
        min_overlap=0.7)

    # Run test
    study, pop, assoc = rd_files(args.filenames, args.compare)
    objgoea = get_objgoea(pop, assoc, args)
    results = objgoea.run_study(study)
    # Check results
    expected_cnts = {'fdr_bh': 17, 'sidak': 5, 'holm': 5, 'bonferroni': 5}
    _chk_results(results, expected_cnts, methods, alpha)
    print("TEST PASSED")


def _chk_results(results, expected_cnts, methods, alpha):
    """Check results."""
    ctr = cx.Counter()
    for ntres in results:
        for method in methods:
            ctr[method] += getattr(ntres, "p_{METHOD}".format(METHOD=method)) < alpha
    for method, num_sig in ctr.most_common():
        assert num_sig == expected_cnts[method]



if __name__ == '__main__':
    test_find_enrichment()

# Copyright (C) 2010-2018, DV Klopfenstein, H Tang. All rights reserved.
