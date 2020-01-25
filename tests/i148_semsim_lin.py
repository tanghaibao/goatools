#!/usr/bin/env python
"""Test for issue 148, Lin Similarity Problems"""

from __future__ import print_function

# Computing basic semantic similarities between GO terms

# Adapted from book chapter written by _Alex Warwick Vesztrocy and Christophe Dessimoz_

# How to compute semantic similarity between GO terms.

# First we need to write a function that calculates the minimum number
# of branches connecting two GO terms.

import os
import sys
from itertools import combinations_with_replacement as combo_w_rplc
from goatools.base import get_godag
from goatools.associations import dnld_annofile
from goatools.anno.gpad_reader import GpadReader
#### from goatools.semantic import semantic_similarity
from goatools.semantic import TermCounts
#### from goatools.semantic import get_info_content
#### from goatools.semantic import deepest_common_ancestor
from goatools.semantic import resnik_sim
from goatools.semantic import lin_sim
#### from goatools.godag.consts import NS2GO

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def test_i148_semsim_lin(prt=sys.stdout):
    """Test for issue 148, Lin Similarity Problems"""
    fin_gpad = os.path.join(REPO, 'goa_human.gpad')
    dnld_annofile(fin_gpad, 'gpad')

    godag = get_godag(os.path.join(REPO, "go-basic.obo"), loading_bar=None)
    annoobj = GpadReader(fin_gpad, godag=godag)

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

    # Print Lin values
    prt.write('           {HDR}\n'.format(HDR=' '.join(goids)))
    for go_row in goids:
        prt.write('{GO} '.format(GO=go_row))
        for go_col in goids:
            prt.write('{L:10.6f} '.format(L=p2v[frozenset([go_row, go_col])]))
        prt.write('\n')


if __name__ == '__main__':
    test_i148_semsim_lin()
