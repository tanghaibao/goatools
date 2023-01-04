#!/usr/bin/env python3
"""Test for issue 148, Lin Similarity if a term has no annotations"""


from os.path import join
from os.path import dirname
from os.path import abspath
from sys import stdout

from itertools import combinations_with_replacement as combo_w_rplc

from goatools.base import get_godag
from goatools.associations import dnld_annofile
from goatools.anno.gpad_reader import GpadReader
#### from goatools.semantic import semantic_similarity
from goatools.semantic import TermCounts
#### from goatools.semantic import get_info_content
#### from goatools.semantic import deepest_common_ancestor
#### from goatools.semantic import resnik_sim
from goatools.semantic import lin_sim
#### from goatools.godag.consts import NS2GO

REPO = join(dirname(abspath(__file__)), "..")

def test_i148_semsim_lin(prt=stdout):
    """Test for issue 148, Lin Similarity if a term has no annotations"""
    fin_gpad = join(REPO, 'goa_human.gpad')
    dnld_annofile(fin_gpad, 'gpad')

    godag = get_godag(join(REPO, "go-basic.obo"), loading_bar=None)
    annoobj = GpadReader(fin_gpad, godag=godag)

    # Researcher-provided GO IDs
    goids = [
        'GO:0042581',
        'GO:0101002',
        'GO:0042582',
        'GO:0070820',
        'GO:0008021',
        'GO:0005766',
        'GO:0016591']

    associations = annoobj.get_id2gos('CC')
    termcounts = TermCounts(godag, associations)

    # Calculate Lin values
    p2v = {frozenset([a, b]): lin_sim(a, b, godag, termcounts) for a, b in combo_w_rplc(goids, 2)}
    _prt_lins_similarity(goids, p2v, prt)


def _prt_lins_similarity(goids, p2v, prt):
    """Print Lin's semantic similarity between all researcher provided GO IDs (goids)"""
    prt.write("\nLin's semantic similarity between researcher-provided GO IDs:\n")
    prt.write(f'           {" ".join(goids)}\n')
    none = 'None     '
    for go_row in goids:
        prt.write(f'{go_row} ')
        for go_col in goids:
            lins_val = p2v[frozenset([go_row, go_col])]
            lins_str = f'{lins_val:<9.6} ' if lins_val is not None else none
            prt.write(f'{lins_str:10} ')
        prt.write('\n')


if __name__ == '__main__':
    test_i148_semsim_lin()
