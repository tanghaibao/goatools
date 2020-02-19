#!/usr/bin/env python3
"""Test for issue 148, Lin Similarity if a term has no annotations"""

import os
import sys
from itertools import combinations_with_replacement as combo_w_rplc
from goatools.obo_parser import GODag
from goatools.anno.gaf_reader import GafReader
from goatools.semantic import TermCounts
from goatools.semantic import resnik_sim
from goatools.semantic import lin_sim
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.gosubdag.plot.gosubdag_plot import GoSubDagPlot

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def test_i148b_semsim_lin(do_plt=False):
    """Test for issue 148, Lin Similarity if a term has no annotations"""
    fin_gaf = os.path.join(REPO, 'tests/data/yangRWC/fig2a_nonleaf0.gaf')
    godag = GODag(os.path.join(REPO, "tests/data/yangRWC/fig2a.obo"))
    annoobj = GafReader(fin_gaf, godag=godag)

    associations = annoobj.get_id2gos('CC')
    tcntobj = TermCounts(godag, associations)

    if do_plt:
        _do_plt(tcntobj, godag)

    goids = list(godag.keys())

    ##print(lin_sim('GO:0000006', 'GO:0000002', godag, tcntobj, 1.0))
    ## print(lin_sim('GO:0005575', 'GO:0005575', godag, tcntobj, 1.0))
    ##return 

    # Calculate Resnik values
    p2r = {frozenset([a, b]): resnik_sim(a, b, godag, tcntobj) for a, b in combo_w_rplc(goids, 2)}
    _prt_values('Resnik', goids, p2r)

    # Calculate Lin values
    p2l = {frozenset([a, b]): lin_sim(a, b, godag, tcntobj) for a, b in combo_w_rplc(goids, 2)}
    _prt_values('Lin', goids, p2l)
    _chk_lin(p2l)
    return

    # Calculate Resnik values
    p2r = {frozenset([a, b]): resnik_sim(a, b, godag, tcntobj) for a, b in combo_w_rplc(goids, 2)}
    _prt_values('Resnik', goids, p2r)

    # Calculate Lin values
    p2l = {frozenset([a, b]): lin_sim(a, b, godag, tcntobj) for a, b in combo_w_rplc(goids, 2)}
    _prt_values('Lin', goids, p2l)
    _chk_lin(p2l)


def _prt_values(desc, goids, p2v, prt=sys.stdout):
    """Print values"""
    prt.write('\n{DESC}\n'.format(DESC=desc))
    prt.write('           {HDR}\n'.format(HDR=' '.join(goids)))
    none = 'None     '
    for go_row in goids:
        prt.write('{GO} '.format(GO=go_row))
        for go_col in goids:
            val = p2v[frozenset([go_row, go_col])]
            txt = '{L:<9.6} '.format(L=val) if val is not None else none
            prt.write('{T:10} '.format(T=txt))
        prt.write('\n')

def _chk_lin(p2l):
    """Check Lin values"""
    for go_pair, lin in p2l.items():
        # Check that a GO compared to iteself if 1.0
        if len(go_pair) == 1:
            assert lin == 1.0, 'Lin({L:10.6f}) {GO}'.format(L=lin, GO=list(go_pair)[0])
        # Check that a main GO compared against it's alt GO is 1.0
        elif go_pair == frozenset(['GO:0005575', 'GO:0000001']):
            assert lin == 1.0, 'Lin({L:10.6f}) GO:0005575 GO:0000001'.format(L=lin)
        if lin is not None:
            assert lin <= 1.0, 'Lin({L:10.6f}) {GOs}'.format(L=lin, GOs=' '.join(sorted(go_pair)))


def _do_plt(tcntobj, godag):
    """Plot the test GO-DAG"""
    gosubdag = GoSubDag(tcntobj.go2obj.keys(), godag, tcntobj=tcntobj)
    GoSubDagPlot(gosubdag).plt_dag('i148b.png')


if __name__ == '__main__':
    test_i148b_semsim_lin(len(sys.argv) != 1)
