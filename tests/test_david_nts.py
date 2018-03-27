#!/usr/bin/env python
"""Test working with DAVID results in DAVID chart files."""

__copyright__ = "Copyright (C) 2016-2017, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import os
import collections as cx
import timeit
from goatools.parsers.david_chart import DavidChartReader
from goatools.test_data.godag_timed import prt_hms


def test_david_chart():
    """Read in a small obo, print list of GO terms and plot."""
    repo = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    david_dir = "{REPO}/data/gjoneska".format(REPO=repo)
    ntobj = cx.namedtuple("david6p8", "TOTAL FDR Bonferroni Benjamini PValue")
    # pylint: disable=bad-whitespace
    fin2exp = {
        "gjoneska2015_chart6p8_Consistent_Decrease.txt": ntobj._make([  297,  24,  25,  48, 209]),
        "gjoneska2015_chart6p8_Transient_Decrease.txt":  ntobj._make([   79,   0,   2,   2,  53]),
        "gjoneska2015_chart6p8_Consistent_Increase.txt": ntobj._make([  418,  33,  33,  76, 307]),
        "gjoneska2015_chart6p8_Transient_Increase.txt":  ntobj._make([  473, 125, 129, 233, 384]),
        "gjoneska2015_chart6p8_Late_Decrease.txt":       ntobj._make([  595, 110, 106, 215, 462]),
        "gjoneska2015_chart6p8_Late_Increase.txt":       ntobj._make([  804,  81,  81, 168, 613]),
    }
    tic = timeit.default_timer()
    fin2obj = {f:DavidChartReader(os.path.join(david_dir, f)) for f in fin2exp.keys()}
    prt_hms(tic, "Created DavidChartReader objects")
    for fin, obj in fin2obj.items():
        ntexp = fin2exp[fin]
        assert ntexp.TOTAL == len(obj.nts)
        obj.prt_num_sig()
        ctr = obj.get_num_sig()
        for fld, cnt_actual in ctr.most_common():
            assert cnt_actual == getattr(ntexp, fld)


if __name__ == '__main__':
    test_david_chart()

# Copyright (C) 2016-2017, DV Klopfenstein, H Tang, All rights reserved.
