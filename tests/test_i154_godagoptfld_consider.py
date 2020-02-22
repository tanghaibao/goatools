#!/usr/bin/env python3
"""Test for issue 148, Lin Similarity if a term has no annotations"""

import sys
import timeit
from goatools.godag.prttime import prt_hms
from goatools.obo_parser import GODag
from goatools.base import download_go_basic_obo


def test_i154_semsim_lin():
    """Test for issue 148, Lin Similarity if a term has no annotations"""
    fin_dag = download_go_basic_obo()
    tic = timeit.default_timer()

    optional_attrs = {'consider', 'replaced_by'}
    load_obsolete = True
    prt = sys.stdout

    godag = GODag(fin_dag, optional_attrs, load_obsolete, prt)
    prt_hms(tic, 'Loaded GO DAG')
    assert godag['GO:0000067'].consider
    assert godag['GO:0003734'].replaced_by == 'GO:0030532'

    godag = GODag(fin_dag, 'consider', load_obsolete, prt)
    prt_hms(tic, 'Loaded GO DAG')
    assert godag['GO:0000067'].consider


if __name__ == '__main__':
    test_i154_semsim_lin()
