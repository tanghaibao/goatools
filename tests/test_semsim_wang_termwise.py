#!/usr/bin/env python3
"""Test Wang_2007"""

__copyright__ = "Copyright (C) 2020-present, DV Klopfenstein. All rights reserved."
__author__ = "DV Klopfenstein"

from os.path import join
from os.path import splitext
import sys
from sys import stdout
import timeit
from numpy.random import shuffle
import networkx as nx

from goatools.randseed import RandomSeed32
from goatools.base import get_godag
from goatools.semsim.termwise.wang import SsWang
from goatools.godag.prttime import prt_hms
from goatools.godag.reldepth import get_go2reldepth

from tests.utils import REPO

# https://github.com/mojaie/pygosemsim
# pygosemsim requires networkx package
# ss for export PYTHONPATH
#    * pygosemsim does not handle alt GO IDs
from pygosemsim import graph
from pygosemsim import similarity


# pylint: disable=line-too-long
def test_semsim_wang(seed=None, num_calcs=1000, prt=stdout):
    """Wang Semantic Similarity tests"""
    # Log file
    logfile = join(REPO, 'test_semsim_wang_termwise.log')
    ## assert not exists(logfile), 'REMOVE TO RUN: {}'.format(logfile)
    # Check that all relationships seem in DAG are expected by SsWang
    fin_godag = join(REPO, 'go-basic.obo')
    # Run randoms
    edge_weights = {
        'is_a': 0.8,
        'part_of': 0.6,
        'regulates': 0.0,
        'negatively_regulates': 0.0,
        'positively_regulates': 0.0,
    }

    # Using only 'is_a' and 'part_of' is not the same as setting the 'regulates' weights to 0:
    # relationships = {'part_of'}
    # pygosemsim loads all relationships, then sets 'regulates' relationships to 0 edge_weight
    relationships = {'part_of', 'regulates', 'negatively_regulates', 'positively_regulates'}
    run = Run(fin_godag, num_calcs, relationships, edge_weights, seed, prt)
    run.randoms(logfile)

    # Test both main GO and alt GO, GO:0008150 GO:0000004, have the same comparison to GO:0008152
    goids = {'GO:0008150', 'GO:0000004', 'GO:0008152'}
    wang = SsWang(goids, run.godag, relationships, edge_weights)
    assert wang.get_sim('GO:0008150', 'GO:0008152') == wang.get_sim('GO:0000004', 'GO:0008152')


class Run:
    """Wang Semantic Similarity tests"""

    def __init__(self, fin_godag, num_calcs, relationships, w_e, seed, prt):
        tic = timeit.default_timer()
        self.godag = get_godag(fin_godag, optional_attrs=['relationship'], prt=prt)
        tic = prt_hms(tic, 'GOATOOLS read godag')
        # Needed because pysemsim not understand cygwin pathes
        self.graph = graph.from_resource(splitext(fin_godag)[0])
        tic = prt_hms(tic, 'pygosemsim read godag')
        self.seedobj = RandomSeed32(seed)
        self.goids = self._init_goids(num_calcs)
        tic = timeit.default_timer()
        self.wang = SsWang(self.goids, self.godag, relationships, w_e)
        self.go2reldepth = 
        tic = prt_hms(tic, 'GOATOOLS wang setup')

    def prt_ancestors(self, goid, prt_if_diff=False):
        """Print ancestors for both Wang and GOATOOLS"""
        a_w = nx.ancestors(self.graph, goid)
        a_g = self.wang.go2dag[goid].ancestors
        if prt_if_diff and a_w != a_g:
            print('{GO} {w:2} Wang  {g:2} GOATOOLS'.format(GO=goid, w=len(a_w), g=len(a_g)))

    def randoms(self, logfile):
        """Run random simulations. Compare SsWang in GOATOOLS to pygosemsim"""
        with open(logfile, 'w') as prt:
            self.seedobj.prt(prt)
            self.seedobj.prt(stdout)
            self._randoms(prt)
            print('  **WROTE: {LOG}'.format(LOG=logfile))

    def _randoms(self, prt):
        """Randomly select GO terms for semantic similarity calculations"""
        #pylint: disable=line-too-long
        goids = self.goids
        go_pairs = [(goids[i], goids[i+1]) for i in range(0, len(self.goids), 2)]
        tic = timeit.default_timer()
        # Information on Python's round, which is used in 2 spots in pygosemsim:
        #     https://stackoverflow.com/questions/13479163/round-float-to-x-decimals
        #     from decimal import Decimal
        #     >>> Decimal('66.66666666666').quantize(Decimal('1e-4'))
        #     Decimal('66.6667')
        #     >>> Decimal('1.29578293').quantize(Decimal('1e-6'))
        #     Decimal('1.295783')
        # In issue, https://github.com/micropython/micropython/issues/3516,
        # https://github.com/mdickinson dreams of deprecating the two-argument form of round in Python....
        #     https://github.com/micropython/micropython/issues/3516#issuecomment-625298591
        # Use the decimal type instead: https://docs.python.org/3.10/library/decimal.html
        acts = [self.wang.get_sim(a, b) for a, b in go_pairs]
        tic = prt_hms(tic, 'GOATOOLS wang calc')
        exps = [similarity.wang(self.graph, a, b) for a, b in go_pairs]
        tic = prt_hms(tic, 'pysemsim wang')
        assert len(acts) == len(exps)
        failures = 0
        for idx, (act, exp, (go_a, go_b)) in enumerate(zip(acts, exps, go_pairs)):
            assert act is not None, self._prt_ab(idx, go_a, go_b, act, exp, stdout)
            assert exp is not None, self._prt_ab(idx, go_a, go_b, act, exp, stdout)
            if abs(exp - act) > 0.02:
                for strm in [prt, stdout]:
                    go_a = goids[2*idx]
                    go_b = goids[2*idx+1]
                    self._prt_ab(idx, go_a, go_b, act, exp, strm)
                stdout.flush()
                prt.flush()
                failures += 1
                self.prt_ancestors(go_a, True)
                self.prt_ancestors(go_b, True)
            else:
                prt.write('{i} PASS {A} {B} pygosemsim={b:f} GOATOOLS={a:f}\n'.format(
                    i=idx, A=goids[2*idx], B=goids[2*idx+1], a=act, b=exp))

    def _prt_ab(self, idx, go_a, go_b, act, exp, strm):
        """Print GO IDs and similarities"""
        ## dif = abs(exp - act) if exp is not None and act is not None else 'XXX'
        strm.write('{i:3}) FAIL {x:2} {A}   {y:2} {B} pygosemsim={b:f} GOATOOLS={a:f} DIFF={ab:f}\n'.format(
            i=idx,
            A=go_a, x=self.godag[go_a].reldepth,
            B=go_b, y=self.godag[go_a].reldepth,
            a=act, b=exp, ab=abs(exp - act)))

    def _init_goids(self, num_calcs):
        """Pick random GO IDs"""
        goids = sorted([o.item_id for o in set(self.godag.values()) if o.namespace == 'biological_process'])
        shuffle(goids)
        return goids[:2*num_calcs]



if __name__ == '__main__':
    SEED = None if len(sys.argv) == 1 else int(sys.argv[1], 16)
    test_semsim_wang(SEED)

# Copyright (C) 2020-present DV Klopfenstein. All rights reserved.
