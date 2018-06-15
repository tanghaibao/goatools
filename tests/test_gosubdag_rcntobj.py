#!/usr/bin/env python
"""Test creation of GoSubDag's rcntobj data member."""

__copyright__ = "Copyright (C) 2016-2017, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import os
import sys
from goatools.base import get_godag
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.gosubdag.godag_rcnt import CountRelatives
from goatools.gosubdag.plot.gosubdag_plot import GoSubDagPlot

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../")


def test_plotgosubdag(prt=sys.stdout):
    """Test creation of GoSubDag's rcntobj data member."""
    objrun = Run("data/i86.obo")
    objrun.prt_goids_all(prt)
    go_sources = set([
        'GO:0000004',  # a BP 15 L00 D00     biological_process
        'GO:0008151',  # a BP 10 L01 D01 B   cellular process
        'GO:0007516',  #   BP  0 L04 D05 ABC hemocyte development
        'GO:0036476']) #   BP  0 L06 D06 AB  neuron death in response to hydrogen peroxide
    objrun.plt_goids("test_gosubdag_rcntobj.png", go_sources)

    # TEST rcntobj creation
    #pylint: disable=bad-whitespace
    rcntobj = CountRelatives(objrun.go2obj_all)
    kws_exp = [
        ({},                  {'rcntobj':rcntobj}),  # New CountRelatives
        ({'rcntobj':None},    {'rcntobj':None}),
        ({'rcntobj':False},   {'rcntobj':None}),
        ({'rcntobj':True},    {'rcntobj':rcntobj}),  # New CountRelatives
        ({'rcntobj':rcntobj}, {'rcntobj':rcntobj}),  # Use CountRelatives in rcntobj
    ]
    for idx, (kws, expected_fields) in enumerate(kws_exp):
        gosubdag = GoSubDag(go_sources, objrun.go2obj_all, prt=prt, **kws)
        _chk_obj(idx, getattr(gosubdag, 'rcntobj'), expected_fields['rcntobj'], CountRelatives)


class Run(object):
    """Printing GO IDs and Plotting; GODag from obo using GoSubDag."""

    def __init__(self, obo):
        self.cwd = os.getcwd()
        self.go2obj_all = get_godag(os.path.join(REPO, "../goatools/", obo))
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
        objplt.plt_dag(os.path.join(self.cwd, fout_img))

def _chk_obj(idx, act_obj, exp_obj, cls):
    """Check that object creation agrees with expected results."""
    if exp_obj is None:
        assert act_obj is None, "{I}) EXP({E}) ACT({A})".format(I=idx, E=exp_obj, A=act_obj)
    else:
        assert isinstance(act_obj, cls)


if __name__ == '__main__':
    test_plotgosubdag()

# Copyright (C) 2016-2017, DV Klopfenstein, H Tang. All rights reserved.
