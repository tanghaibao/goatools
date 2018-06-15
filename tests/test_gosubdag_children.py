#!/usr/bin/env python
"""Test creation of GoSubDag's rcntobj data member."""

from __future__ import print_function

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import os
from goatools.base import get_godag
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.gosubdag.plot.gosubdag_plot import GoSubDagPlot

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../")


def test_plotgosubdag():
    """Test creation of GoSubDag's rcntobj data member."""
    objrun = Run("../goatools/tests/data/mini_obo.obo")
    #objrun.prt_goids_all(prt)
    in_exp = [
        (["GO:0000002"], set(["GO:0000001", "GO:0000002"]), {}),
        (["GO:0000002"], set(["GO:0000001", "GO:0000002",
                              "GO:0000003", "GO:0000005",
                              "GO:0000006", "GO:0000008",
                              "GO:0000010",             ]), {'children':True}),
        (["GO:0000002"], set(["GO:0000001", "GO:0000002"]), {'children':0}),
        (["GO:0000002"], set(["GO:0000001", "GO:0000002",
                              "GO:0000003", "GO:0000005",
                              "GO:0000006", "GO:0000008",
                              "GO:0000010",             ]), {'children':1}),
        (["GO:0000008"], set(["GO:0000001", "GO:0000003",
                              "GO:0000006", "GO:0000008"]), {}),
        (["GO:0000008"], set(["GO:0000001", "GO:0000003",
                              "GO:0000006", "GO:0000008"]), {'children':False}),
        (["GO:0000008"], set(["GO:0000001", "GO:0000003",
                              "GO:0000006", "GO:0000008"]), {'children':None}),
        (["GO:0000008"], set(["GO:0000001", "GO:0000003",
                              "GO:0000006", "GO:0000008",
                              "GO:0000002", "GO:0000005",
                              "GO:0000004", "GO:0000007",
                              "GO:0000009", "GO:0000010"]), {'children':True}),
    ]
    for in_goids, exp_goids, kws in in_exp:
        objrun.run(in_goids, exp_goids, **kws)
    # objrun.plt_goids("mini_obo.png", None)


class Run(object):
    """Printing GO IDs and Plotting; GODag from obo using GoSubDag."""

    def __init__(self, obo):
        self.go2obj_all = get_godag(os.path.join(REPO, obo))
        self.gosubdag_all = GoSubDag(None, self.go2obj_all)
        self.prtfmt = self.gosubdag_all.prt_attr['fmta']

    def prt_goids_all(self, prt):
        """Print all GO IDs, including alternate GO IDs, in GODag."""
        self.gosubdag_all.prt_goids(prtfmt=self.prtfmt, prt=prt)

    def plt_goids(self, fout_img, go_sources):
        """Plot GO IDs."""
        # % src/bin/go_plot.py GOs --obo=../goatools/data/i86.obo --outfile=t00.jpg --mark_alt_id
        gosubdag = GoSubDag(go_sources, self.go2obj_all)
        objplt = GoSubDagPlot(gosubdag, mark_alt_id=True)
        objplt.plt_dag(os.path.join(REPO, fout_img))

    def run(self, go_sources, exp_gos, **kws):
        """Create GoSubDag using specified GO sources."""
        print("\nSRCS: {GOs}".format(GOs=go_sources))
        gosubdag = GoSubDag(go_sources, self.go2obj_all, **kws)
        gosubdag.prt_goids(gosubdag.go2nt)
        assert set(gosubdag.go2nt) == exp_gos, "ACT({}) != EXP({})\n{} {}".format(
            sorted(gosubdag.go2nt), sorted(exp_gos), go_sources, kws)


if __name__ == '__main__':
    test_plotgosubdag()

# Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved.
