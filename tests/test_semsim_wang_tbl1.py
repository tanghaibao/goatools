#!/usr/bin/env python3
"""Test S-value for Table 1 in Wang_2007"""

__copyright__ = "Copyright (C) 2020-present, DV Klopfenstein. All rights reserved."
__author__ = "DV Klopfenstein"

from os.path import join
from sys import stdout
from utils import REPO

from goatools.obo_parser import GODag
from goatools_alpha.semsim.termwise.wang import SsWang
from data.ssWang.tbl1 import GO2SVALUE


def test_semsim_wang_tbl1(prt=stdout):
    """Test S-value for Table 1 in Wang_2007"""
    fin_godag = join(REPO, 'tests/data/ssWang/fig1.obo')

    # Wang 2.1 Test svalues using Fig 1 and Table 1
    relationships = ['part_of']
    godag = GODag(fin_godag, optional_attrs=['relationship'], prt=prt)
    wang = SsWang(godag, relationships)
    go_a = 'GO:0043231'
    wang.add_goid(go_a)
    dag_a = wang.go2subdag[go_a]
    _chk_svalues_a(dag_a)

    # Wang 2.2 Test Semantic similarity of GO terms
    go_b = 'GO:0043229'
    wang.add_goid(go_b)
    sim_ab = wang.get_semsim(go_a, go_b)
    assert abs(sim_ab - 0.7727) < .0001


def _chk_svalues_a(dag_a):
    """Check values against Table 1"""
    assert len(dag_a.go2svalue) == len(GO2SVALUE)
    for goid, svalue_act in dag_a.go2svalue.items():
        svalue_exp = GO2SVALUE[goid]
        assert abs(svalue_exp - svalue_act) < .001, 'MISMATCH EXP({}) != ACT({})'.format(
            svalue_exp, svalue_act)


if __name__ == '__main__':
    test_semsim_wang_tbl1()

# Copyright (C) 2020-present DV Klopfenstein. All rights reserved.
