#!/usr/bin/env python3
"""Test S-value for Table 1 in Wang_2007"""

__copyright__ = "Copyright (C) 2020-present, DV Klopfenstein. All rights reserved."
__author__ = "DV Klopfenstein"

from os.path import join
from sys import stdout

from goatools.obo_parser import GODag
from goatools.semsim.termwise.wang import SsWang

from tests.utils import REPO
from tests.data.ssWang.tbl1 import GO2SVALUE


def test_semsim_wang(prt=stdout):
    """Wang Semantic Similarity tests"""
    fin_godag = join(REPO, 'tests/data/ssWang/fig1.obo')
    obj = Run(fin_godag, prt)
    obj.run_semsim_wang_tbl1('GO:0043231', 'GO:0043229')  # main main
    obj.run_semsim_wang_tbl1('GO:0043231', 'GO:0042222')  # main  alt
    obj.run_semsim_wang_tbl1('GO:0043333', 'GO:0043229')  #  alt main
    obj.run_semsim_wang_tbl1('GO:0043333', 'GO:0042222')  #  alt  alt


class Run:
    """Wang Semantic Similarity tests"""

    def __init__(self, fin_godag, prt):
        self.godag = GODag(fin_godag, optional_attrs=['relationship'], prt=prt)

    def run_semsim_wang_tbl1(self, go_a, go_b):
        """Test S-value for Table 1 in Wang_2007 (Alt ID)a65Gkk"""
        print('RUN: {A} {B}'.format(A=go_a, B=go_b))
        relationships = ['part_of']
        wang = SsWang(self.godag, relationships)
        wang.add_goid(go_a)
        dag_a = wang.go2subdag[go_a]
        self._chk_svalues_a(dag_a)

        # Wang 2.2 Test Semantic similarity of GO terms
        wang.add_goid(go_b)
        sim_ab = wang.get_semsim(go_a, go_b)
        assert abs(sim_ab - 0.7727) < .0001

    @staticmethod
    def _chk_svalues_a(dag_a):
        """Check values against Table 1"""
        assert len(dag_a.go2svalue) == len(GO2SVALUE)
        for goid, svalue_act in dag_a.go2svalue.items():
            svalue_exp = GO2SVALUE[goid]
            assert abs(svalue_exp - svalue_act) < .001, 'MISMATCH EXP({}) != ACT({})'.format(
                svalue_exp, svalue_act)


if __name__ == '__main__':
    test_semsim_wang()

# Copyright (C) 2020-present DV Klopfenstein. All rights reserved.
