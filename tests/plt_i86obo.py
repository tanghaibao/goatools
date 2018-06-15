#!/usr/bin/env python
"""Read in a small obo, print list of GO terms and plot."""

__copyright__ = "Copyright (C) 2016-2017, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import os
import sys
from goatools.obo_parser import GODag
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.gosubdag.plot.gosubdag_plot import GoSubDagPlot

def test_rpt(prt=sys.stdout):
    """Read in a small obo, print list of GO terms and plot."""
    png = "i86_godag.png"
    gosubdag = _get_gosubdag()
    _prt_goids(gosubdag, prt)
    goobjplt = _get_goobjplt(gosubdag)
    goobjplt.plt_dag(png)

def _get_goobjplt(gosubdag):
    """STEP 3) Get a plotting object."""
    go_sources = set(["GO:0036476", "GO:0007516"])
    gopltdag = GoSubDag(go_sources, gosubdag.go2obj)
    return GoSubDagPlot(gopltdag)

def _prt_goids(gosubdag, prt):
    """STEP 2) Print a list of GO IDs in the GoSubDag."""
    sortby = lambda nt: [nt.NS, -1*nt.dcnt, nt.depth]
    gosubdag.prt_goids(gosubdag.go_sources, sortby=sortby, prt=prt)

def _get_gosubdag():
    """STEP 1) Get GoSubDag containing a small test obo."""
    repo = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../")
    fin_obo = os.path.join(repo, "data/i86.obo")
    go2obj = GODag(fin_obo)
    return GoSubDag(None, go2obj)

if __name__ == '__main__':
    test_rpt()

# Copyright (C) 2016-2017, DV Klopfenstein, H Tang, All rights reserved.
