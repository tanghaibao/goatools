#!/usr/bin/env python
"""Ancestors/Descendants."""

from __future__ import print_function

import os
import sys
import timeit
import numpy as np
from numpy.random import shuffle
from scipy import stats

from goatools.base import download_go_basic_obo
from goatools.obo_parser import GODag
from goatools.godag.prttime import prt_hms
from goatools.gosubdag.gosubdag import GoSubDag


def test_go_pools():
    """Print a comparison of GO terms from different species in two different comparisons."""
    objr = _Run()
    # Check that subset GoSubDags have the same ancestors/descendants as Full GoSubDag
    # pylint: disable=no-member
    for qty in np.random.randint(10, 100, size=10):
        print("")
        goids = objr.get_goids_rand(qty)
        # No relationships loaded; GoSubDag subset equivalent to Full subset?
        gosubdag_r0 = objr.get_gosubdag_r0(goids)
        for goid in gosubdag_r0.go2obj:
            r0_u = gosubdag_r0.rcntobj.go2ancestors.get(goid)
            assert r0_u == objr.gosubdag_r0.rcntobj.go2ancestors.get(goid)
            r0_d = gosubdag_r0.rcntobj.go2descendants.get(goid)
            assert r0_d == objr.gosubdag_r0.rcntobj.go2descendants.get(goid)
        # All relationships loaded; GoSubDag(r0) vs. GoSubDag(r1)
        gosubdag_r1 = objr.get_gosubdag_r1(goids)
        assert gosubdag_r0.go_sources == gosubdag_r1.go_sources
        assert set(gosubdag_r0.go2obj).issubset(gosubdag_r1.go2obj)
        cnts = {'r0_u':[], 'r1_u':[], 'r0_d':[], 'r1_d':[]}
        for goid in gosubdag_r0.go2obj:
            r1_d = gosubdag_r1.rcntobj.go2descendants.get(goid, set())
            r0_d = gosubdag_r0.rcntobj.go2descendants.get(goid, set())
            assert r0_d.issubset(r1_d), "R1({}) R0({})".format(len(r1_d), len(r0_d))
            r0_u = gosubdag_r0.rcntobj.go2ancestors.get(goid, set())
            r1_u = gosubdag_r1.rcntobj.go2ancestors.get(goid, set())
            assert r0_u.issubset(r1_u), "R1({}) R0({})".format(len(r1_u), len(r0_u))
            cnts['r0_u'].append(len(r0_u))
            cnts['r1_u'].append(len(r1_u))
            cnts['r0_d'].append(len(r0_d))
            cnts['r1_d'].append(len(r1_d))
        objr.prt_cnts(cnts)


class _Run:
    """Group entire go-basic.obo"""

    obo = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../go-basic.obo")

    def __init__(self):
        download_go_basic_obo(self.obo, sys.stdout, loading_bar=None)
        self.godag_r0 = GODag(self.obo)
        self.godag_r1 = GODag(self.obo, optional_attrs=set(['relationship']))
        self.goids = list(set(o.id for o in self.godag_r0.values()))
        # GoSubDag (plain)
        tic = timeit.default_timer()
        self.gosubdag_r0 = GoSubDag(self.goids, self.godag_r0, prt=None)
        prt_hms(tic, "GoSubDag r0 {N:4} GOs {S:3} srcs".format(
            N=len(self.gosubdag_r0.go2obj), S=len(self.gosubdag_r0.go_sources)))
        # GoSubDag with relationships
        self.gosubdag_r1 = GoSubDag(self.goids, self.godag_r1, prt=None, relationships=True)
        prt_hms(tic, "GoSubDag r1 {N:4} GOs {S:3} srcs".format(
            N=len(self.gosubdag_r1.go2obj), S=len(self.gosubdag_r1.go_sources)))

    def prt_cnts(self, cnts):
        """Compare ancestor/descendant counts with relatives=False/True."""
        k2v = {k:self.str_stats(v) for k, v in cnts.items()}
        print(k2v)

    @staticmethod
    def str_stats(vals):
        """Print statistics on values."""
        ntd = stats.describe(vals)
        std = int(round(np.sqrt(ntd.variance)))
        return "({m} {M}) STD={STD:,}".format(m=ntd.minmax[0], M=ntd.minmax[1], STD=std)

    def get_gosubdag_r0(self, goids):
        """Return a GoSubDag with N randomly chosen GO sources."""
        tic = timeit.default_timer()
        gosubdag = GoSubDag(goids, self.godag_r0, relationships=None,
                            #rcntobj=self.gosubdag_r0.rcntobj,
                            prt=None)
        prt_hms(tic, "GoSubDag r0 {N:4} GOs {S:3} srcs".format(
            N=len(gosubdag.go2obj), S=len(gosubdag.go_sources)))
        return gosubdag

    def get_gosubdag_r1(self, goids):
        """Return a GoSubDag with N randomly chosen GO sources."""
        tic = timeit.default_timer()
        gosubdag = GoSubDag(goids, self.godag_r1, relationships=True,
                            #rcntobj=self.gosubdag_r1.rcntobj,
                            prt=None)
        prt_hms(tic, "GoSubDag r1 {N:4} GOs {S:3} srcs".format(
            N=len(gosubdag.go2obj), S=len(gosubdag.go_sources)))
        return gosubdag

    def get_goids_rand(self, qty):
        """Return N randomly chosen GO IDs."""
        shuffle(self.goids)
        return self.goids[:qty]


if __name__ == '__main__':
    test_go_pools()
