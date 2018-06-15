#!/usr/bin/env python
"""Test plotting of various GoSubDag options."""

__copyright__ = "Copyright (C) 2016-2017, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import os
import sys
from goatools.base import get_godag
from goatools.base import dnld_gaf
from goatools.associations import read_gaf
from goatools.semantic import TermCounts
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.gosubdag.plot.gosubdag_plot import GoSubDagPlot


REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../")

# TBD MOVE TO GOATOOLS TEST PKG
class Run(object):
    """Objects for running plotting test."""

    def __init__(self, obo, gaf, prt):
        self.prt = prt
        self.cwd = os.getcwd()
        # Gene Ontologies
        self.go2obj_all = get_godag(os.path.join(REPO, "../goatools/", obo))
        # Annotations
        #_file_gaf = dnld_gaf(os.path.join(REPO, gaf))
        _file_gaf = dnld_gaf(gaf)
        print("GAF: {GAF}\n".format(GAF=_file_gaf))
        self.gene2gos = read_gaf(_file_gaf)
        self.tcntobj = TermCounts(self.go2obj_all, self.gene2gos)
        # GoSubDag
        self.gosubdag_all = GoSubDag(None, self.go2obj_all, tcntobj=self.tcntobj, prt=prt)
        self.prtfmt = self.gosubdag_all.prt_attr['fmta']

    def prt_goids_all(self, prt):
        """Print all GO IDs, including alternate GO IDs, in GODag."""
        self.gosubdag_all.prt_goids(prtfmt=self.prtfmt, prt=prt)

    def plt_goids(self, fout_img, go_sources):
        """Plot GO IDs."""
        # % src/bin/go_plot.py GOs --obo=../goatools/data/i86.obo --outfile=t00.jpg --mark_alt_id
        gosubdag = GoSubDag(go_sources, self.gosubdag_all.go2obj, prt=self.prt,
                            # rcntobj=False,
                            rcntobj=self.gosubdag_all.rcntobj,
                            go2nt=self.gosubdag_all.go2nt)
        prtfmt = gosubdag.prt_attr['fmta']
        goids_plt = GoSubDagPlot(gosubdag).get_goids_plt()
        self.prt.write("\n{N} GO IDs\n".format(N=len(goids_plt)))
        gosubdag.prt_goids(goids_plt, prtfmt=prtfmt, prt=self.prt)
        objplt = GoSubDagPlot(gosubdag, mark_alt_id=True)
        objplt.plt_dag(os.path.join(self.cwd, fout_img))



def test_plotgosubdag(prt=sys.stdout):
    """Test plotting of various GoSubDag options."""
    objrun = Run("data/i86.obo", "goa_human", prt)
    # objrun.prt_goids_all(prt)
    go_sources = set([
        'GO:0000004',  # a BP 15 L00 D00     biological_process
        'GO:0008151',  # a BP 10 L01 D01 B   cellular process
        'GO:0007516',  #   BP  0 L04 D05 ABC hemocyte development
        'GO:0036476']) #   BP  0 L06 D06 AB  neuron death in response to hydrogen peroxide
    objrun.plt_goids("test_gosubdag_tcntobj.png", go_sources)

#     kws_exp = [
#         ({},                  {'rcntobj':rcntobj}),
#         ({'rcntobj':None},    {'rcntobj':None}),
#         ({'rcntobj':False},   {'rcntobj':None}),
#         ({'rcntobj':True},    {'rcntobj':rcntobj}),
#         ({'rcntobj':rcntobj}, {'rcntobj':rcntobj}),
#
#         #({},                  {'tcntobj':tcntobj}),
#         #({'tcntobj':None},    {'tcntobj':None}),
#         #({'tcntobj':False},   {'tcntobj':None}),
#         #({'tcntobj':True},    {'tcntobj':tcntobj}),
#         #({'tcntobj':tcntobj}, {'tcntobj':tcntobj}),
#     ]
#     for idx, (kws, expected_fields) in enumerate(kws_exp):
#         gosubdag = GoSubDag(go_sources, objrun.go2obj_all, prt=prt, **kws)
#         _chk_obj(getattr(gosubdag, 'rcntobj'), expected_fields['rcntobj'], CountRelatives)
#
# def _chk_obj(act_obj, exp_obj, cls):
#     """Check that object creation agrees with expected results."""
#     if exp_obj is None:
#         assert act_obj is None
#     else:
#         assert isinstance(act_obj, cls)
#
# # def _chk_rcntobj(idx, kws, gosubdag, expected_fields):
# #     """Check that an rcntobj was created or not created."""
# #     print idx, kws, expected_fields, gosubdag.rcntobj
#
#     #     # goids = kws['GO'] if 'GO' in kws else set(get_go2obj_unique(go2obj))
#     #     # print('CLI: FAST GoSubDag 0 -------------------')
#     #     # gosubdag = GoSubDag(goids, go2obj, rcntobj=False)
#     #     # print('CLI: RCNTOBJ({})'.format(gosubdag.rcntobj))
#     #     # gosubdag.prt_goids()
#     #     # print('CLI: FAST GoSubDag 1 -------------------')
#     #     # tcntobj = self._get_tcntobj(kws, gosubdag.go2obj)
#     #     # print('CLI: TermCounts INITed -------------------')
#     #     # self.gosubdag = GoSubDag(goids, go2obj, tcntobj=tcntobj)
#     #     # # self.gosubdag.prt_goids()
#     #     # print('CLI: FAST GoSubDag 2 -------------------')
#     #     # objcolor = Go2Color(self.gosubdag, None)
#     #     # objplt = GoSubDagPlot(self.gosubdag, Go2Color=objcolor, **kws)
#     #     # fout_img = self.get_outfile(goids, **kws)
#     #     # objplt.plt_dag(fout_img)


if __name__ == '__main__':
    test_plotgosubdag()

# Copyright (C) 2016-2017, DV Klopfenstein, H Tang. All rights reserved.
