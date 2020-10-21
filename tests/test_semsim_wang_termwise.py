#!/usr/bin/env python3
"""Test Wang_2007"""

__copyright__ = "Copyright (C) 2020-present, DV Klopfenstein. All rights reserved."
__author__ = "DV Klopfenstein"

from os.path import join
from os.path import splitext
from os.path import exists
import sys
from sys import stdout
import timeit
from numpy.random import shuffle

from goatools.randseed import RandomSeed32
from goatools.base import get_godag
from goatools.semsim.termwise.wang import SsWang
from goatools.godag.consts import RELATIONSHIP_SET
from goatools.godag.prttime import prt_hms
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.godag.go_tasks import get_go2ancestors

from tests.utils import REPO

# https://github.com/mojaie/pygosemsim
# pygosemsim requires networkx package
# ss for export PYTHONPATH
#    * pygosemsim does not handle alt GO IDs
from pygosemsim import graph
from pygosemsim import similarity


def test_semsim_wang(seed=None, prt=stdout):
    """Wang Semantic Similarity tests"""
    fin_godag = join(REPO, 'go-basic.obo')
    run = Run(fin_godag, seed, prt)
    run.chk_relationships()
    relationships = {'part_of'}
    wang = SsWang(run.godag, relationships)
    wang.add_goids({"GO:0004340", "GO:0019158"})
    print(similarity.wang(run.graph, "GO:0004340", "GO:0019158"))
    print(wang.get_sim("GO:0004340", "GO:0019158"))
    run.randoms(100, relationships)


class Run:
    """Wang Semantic Similarity tests"""

    def __init__(self, fin_godag, seed, prt):
        self.godag = get_godag(fin_godag, optional_attrs=['relationship'], prt=prt)
        # Needed because pysemsim not understand cygwin pathes
        self.graph = graph.from_resource(splitext(fin_godag)[0])
        self.seedobj = RandomSeed32(seed)

    def randoms(self, num_calcs, relationships):
        """Run random simulations. Compare SsWang in GOATOOLS to pygosemsim"""
        logfile = join(REPO, 'test_semsim_wang_termwise.log')
        assert not exists(logfile), 'REMOVE TO RUN: {}'.format(logfile)
        with open(logfile, 'w') as prt:
            self.seedobj.prt(prt)
            self.seedobj.prt(stdout)
            self._randoms(num_calcs, relationships, prt)
            print('  **WROTE: {LOG}'.format(LOG=logfile))

    def _randoms(self, num_calcs, relationships, prt):
        """Randomly select GO terms for semantic similarity calculations"""
        #pylint: disable=line-too-long
        goids = sorted([o.item_id for o in set(self.godag.values()) if o.namespace == 'biological_process'])
        shuffle(goids)
        rng = range(0, num_calcs*2, 2)
        tic = timeit.default_timer()
        s_godag = self.godag
        goobjs = [s_godag[go] for go in goids[:2*num_calcs]]
        go2ancestors = get_go2ancestors(goobjs, relationships)
        tic = prt_hms(tic, 'goatools go2ancestors')
        gosubdag = GoSubDag(goids[:2*num_calcs], self.godag, relationships)
        tic = prt_hms(tic, 'goatools gosubdag')
        wang = SsWang(s_godag, relationships)
        wang.add_goids(goids[:2*num_calcs])
        tic = prt_hms(tic, 'goatools wang load')
        acts = [wang.get_sim(goids[i], goids[i+1]) for i in rng]
        tic = prt_hms(tic, 'goatools wang calc')
        exps = [similarity.wang(self.graph, goids[i], goids[i+1]) for i in rng]
        tic = prt_hms(tic, 'pysemsim wang')
        assert len(acts) == len(exps)
        assert len(acts) == num_calcs
        for idx, (act, exp) in enumerate(zip(acts, exps)):
            if abs(exp - act) < 0.0001:
                for strm in [prt, stdout]:
                    strm.write('{i} {A} {B} pygosemsim={b:f} GOATOOLS={a:f}'.format(
                        i=idx, A=goids[2*idx], B=goids[2*idx+1], a=act, b=exp))
                stdout.flush()
                prt.flush()
                sys.exit(1)

    def chk_relationships(self):
        """Check that actual relationships are expected"""
        rels_all = set()
        for goterm in self.godag.values():
            rels_cur = goterm.relationship.keys()
            if rels_cur:
                rels_all.update(rels_cur)
        assert rels_all == RELATIONSHIP_SET, 'UNEXPECTED RELATIONSHIPS'
        print('**PASSED: EXPECTED GODag  RELATIONSHIPS: {R}'.format(R=sorted(rels_all)))
        rels_all.add('is_a')
        rels_act = set(SsWang.dflt_rel2scf.keys())
        assert rels_all == rels_act, 'BAD SsWang RELATIONSHIPS: {Rs}'.format(Rs=rels_act)
        print('**PASSED: EXPECTED SsWang RELATIONSHIPS: {R}'.format(R=sorted(rels_act)))


if __name__ == '__main__':
    test_semsim_wang()

# Copyright (C) 2020-present DV Klopfenstein. All rights reserved.
