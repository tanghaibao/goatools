#!/usr/bin/env python3
"""Test S-value for Table 1 in Wang_2007"""

__copyright__ = "Copyright (C) 2020-present, DV Klopfenstein. All rights reserved."
__author__ = "DV Klopfenstein"

from os.path import join
from sys import stdout

from goatools.base import get_godag
from goatools.semsim.termwise.wang import SsWang
from goatools.godag.consts import RELATIONSHIP_SET

from tests.utils import REPO
from tests.data.ssWang.tbl1 import GO2SVALUE


def test_semsim_wang(prt=stdout):
    """Wang Semantic Similarity tests"""
    fin_godag = join(REPO, 'go-basic.obo')
    run = Run(fin_godag, prt)
    run.chk_relationships()


class Run:
    """Wang Semantic Similarity tests"""

    def __init__(self, fin_godag, prt):
        self.godag = get_godag(fin_godag, optional_attrs=['relationship'], prt=prt)

    @staticmethod
    def _chk_svalues_a(dag_a):
        """Check values against Table 1"""
        assert len(dag_a.go2svalue) == len(GO2SVALUE)
        for goid, svalue_act in dag_a.go2svalue.items():
            svalue_exp = GO2SVALUE[goid]
            assert abs(svalue_exp - svalue_act) < .001, 'MISMATCH EXP({}) != ACT({})'.format(
                svalue_exp, svalue_act)

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
