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
    rels = ['part_of']
    # pylint: disable=line-too-long
    assert abs(obj.run_semsim_wang_tbl1('GO:0043231', 'GO:0043229', rels) - 0.7727) < .0001  # main main
    assert abs(obj.run_semsim_wang_tbl1('GO:0043231', 'GO:0042222', rels) - 0.7727) < .0001  # main  alt
    assert abs(obj.run_semsim_wang_tbl1('GO:0043333', 'GO:0043229', rels) - 0.7727) < .0001  #  alt main
    assert abs(obj.run_semsim_wang_tbl1('GO:0043333', 'GO:0042222', rels) - 0.7727) < .0001  #  alt  alt
    rels = []
    ss_val0 = obj.run_semsim_wang_tbl1('GO:0043231', 'GO:0043229', rels)  # main main
    ss_val1 = obj.run_semsim_wang_tbl1('GO:0043231', 'GO:0042222', rels)  # main  alt
    ss_val2 = obj.run_semsim_wang_tbl1('GO:0043333', 'GO:0043229', rels)  #  alt main
    ss_val3 = obj.run_semsim_wang_tbl1('GO:0043333', 'GO:0042222', rels)  #  alt  alt
    assert ss_val0 == ss_val1
    assert ss_val0 == ss_val2
    assert ss_val0 == ss_val3


# pylint: disable=too-few-public-methods
class Run:
    """Wang Semantic Similarity tests"""

    def __init__(self, fin_godag, prt):
        self.godag = GODag(fin_godag, optional_attrs=['relationship'], prt=prt)

    def run_semsim_wang_tbl1(self, go_a, go_b, relationships):
        """Test S-value for Table 1 in Wang_2007 (Alt ID)a65Gkk"""
        wang = SsWang(self.godag, relationships)
        wang.add_goid(go_a)
        dag_a = wang.go2subdag[go_a]
        if relationships == ['part_of']:
            self._chk_svalues_a(dag_a)

        # Wang 2.2 Test Semantic similarity of GO terms
        wang.add_goid(go_b)
        ss_rel = wang.get_semsim(go_a, go_b)
        print('RUN: {A} {B} rels={R} SS={S:6.4f}'.format(
            A=go_a, B=go_b, R=relationships, S=ss_rel))
        return ss_rel

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
