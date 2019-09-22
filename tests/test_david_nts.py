#!/usr/bin/env python
"""Test working with DAVID results in DAVID chart files."""

from __future__ import print_function

__copyright__ = "Copyright (C) 2016-2017, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import os
import collections as cx
import timeit
from goatools.parsers.david_chart import DavidChartReader
from goatools.godag.prttime import prt_hms


def test_david_chart():
    """Read in a small obo, print list of GO terms and plot."""
    repo = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    david_dir = "{REPO}/data/gjoneska_pfenning".format(REPO=repo)
    ntobj = cx.namedtuple("david6p8", "TOTAL FDR Bonferroni Benjamini PValue")
    # pylint: disable=bad-whitespace
    fin2exp = {
        "david_chart6p8_Consistent_Decrease.txt": ntobj._make([ 1773, 259, 249,  432, 1316]),
        "david_chart6p8_Transient_Decrease.txt":  ntobj._make([  423,   0,   2,    2,  246]),
        "david_chart6p8_Consistent_Increase.txt": ntobj._make([ 2359, 353, 308,  781, 1868]),
        "david_chart6p8_Transient_Increase.txt":  ntobj._make([ 2191, 658, 652, 1105, 1786]),
        "david_chart6p8_Late_Decrease.txt":       ntobj._make([ 2752, 591, 568, 1153, 2187]),
        "david_chart6p8_Late_Increase.txt":       ntobj._make([ 4597, 708, 616, 1715, 3603]),
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
            assert cnt_actual == getattr(ntexp, fld), "{FIN}: {FLD} Act({ACT}) Exp({EXP})".format(
                FIN=fin, FLD=fld, ACT=cnt_actual, EXP=getattr(ntexp, fld))


if __name__ == '__main__':
    test_david_chart()

# Copyright (C) 2016-2017, DV Klopfenstein, H Tang, All rights reserved.
