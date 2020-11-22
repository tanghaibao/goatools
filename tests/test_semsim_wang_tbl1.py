#!/usr/bin/env python3
"""Wang Semantic Similarity tests using Wang_2007 for expected results"""

__copyright__ = "Copyright (C) 2020-present, DV Klopfenstein. All rights reserved."
__author__ = "DV Klopfenstein"

from os.path import join
from os.path import splitext
from sys import stdout

from goatools.obo_parser import GODag
from goatools.semsim.termwise.wang import SsWang

from tests.utils import REPO
from tests.data.ssWang.tbl1 import GO2SVALUE

# https://github.com/mojaie/pygosemsim
# pygosemsim requires networkx package
# ss for export PYTHONPATH
#    * pygosemsim does not handle alt GO IDs
from pygosemsim import graph
from pygosemsim import similarity


def test_semsim_wang(prt=stdout):
    """Wang Semantic Similarity tests using Wang_2007 for expected results"""
    fin_godag = join(REPO, 'tests/data/ssWang/fig1.obo')
    obj = Run(fin_godag, prt)
    rels = ['part_of']

    # RUN WITH optional RELATIONSHIPS
    # pylint: disable=line-too-long
    assert abs(obj.run_semsim_wang_tbl1('GO:0043231', 'GO:0043229', rels) - 0.7727) < .0001  # main main
    assert abs(obj.run_semsim_wang_tbl1('GO:0043231', 'GO:0042222', rels) - 0.7727) < .0001  # main  alt
    assert abs(obj.run_semsim_wang_tbl1('GO:0043333', 'GO:0043229', rels) - 0.7727) < .0001  #  alt main
    assert abs(obj.run_semsim_wang_tbl1('GO:0043333', 'GO:0042222', rels) - 0.7727) < .0001  #  alt  alt
    # NOTE: pygosemsim rounds the semantic similarity from the original, 0.7727 to 0.773 and loads 'part_of' relationships
    print(similarity.wang(obj.graph, 'GO:0043231', 'GO:0043229'))

    # RUN without optional RELATIONSHIPS
    rels = []
    ss_val0 = obj.run_semsim_wang_tbl1('GO:0043231', 'GO:0043229', rels)  # main main
    ss_val1 = obj.run_semsim_wang_tbl1('GO:0043231', 'GO:0042222', rels)  # main  alt
    ss_val2 = obj.run_semsim_wang_tbl1('GO:0043333', 'GO:0043229', rels)  #  alt main
    ss_val3 = obj.run_semsim_wang_tbl1('GO:0043333', 'GO:0042222', rels)  #  alt  alt
    assert ss_val0 == ss_val1
    assert ss_val0 == ss_val2
    assert ss_val0 == ss_val3
    print('**PASSED')


# pylint: disable=too-few-public-methods
class Run:
    """Wang Semantic Similarity tests"""

    def __init__(self, fin_godag, prt):
        self.godag = GODag(fin_godag, optional_attrs=['relationship'], prt=prt)
        self.graph = graph.from_resource(splitext(fin_godag)[0])

    def run_semsim_wang_tbl1(self, go_a, go_b, relationships):
        """Test S-value for Table 1 in Wang_2007 (Alt ID)a65Gkk"""
        wang = SsWang({go_a, go_b}, self.godag, relationships)
        dag_a = wang.go2dag[go_a]
        if relationships == ['part_of']:
            self._chk_svalues_a(dag_a)

        # Wang 2.2 Test Semantic similarity of GO terms
        ss_rel = wang.get_sim(go_a, go_b)
        print('RUN: {A} {B} rels={R} SS={S:6.4f}'.format(
            A=go_a, B=go_b, R=relationships, S=ss_rel))
        return ss_rel

    @staticmethod
    def _chk_svalues_a(dag_a):
        """Check values against Table 1"""
        assert len(dag_a.go2svalue) == len(GO2SVALUE)
        for goid, svalue_act in dag_a.go2svalue.items():
            svalue_exp = GO2SVALUE[goid]
            assert abs(svalue_exp - svalue_act) < .0001, 'MISMATCH EXP({}) != ACT({})'.format(
                svalue_exp, svalue_act)


if __name__ == '__main__':
    test_semsim_wang()

# Copyright (C) 2020-present DV Klopfenstein. All rights reserved.
