#!/usr/bin/env python
"""Test that comparing two identical GO IDs returns true"""
# https://github.com/tanghaibao/goatools/issues/150

import os
## import sys
from goatools.obo_parser import GODag
## from goatools.anno.gaf_reader import GafReader
## from goatools.semantic import TermCounts
from goatools.semantic import semantic_similarity

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def test_semantic_i150():
    """Test that comparing two identical GO IDs returns true"""
    fin_dag = os.path.join(REPO, 'tests/data/yangRWC/fig1a.obo')
    ## fin_gaf = os.path.join(REPO, 'tests/data/yangRWC/fig2a_nonleaf0.gaf')
    # Read files
    godag = GODag(fin_dag)
    ## objanno = GafReader(fin_gaf)
    ## gene2gos = objanno.get_id2gos(namespace='CC')
    ## # Termcounts
    ## termcounts = TermCounts(godag, gene2gos, prt=sys.stdout)
    # Compare all GO terms with itself
    for goterm in set(godag.values()):
        goid = goterm.item_id
        assert semantic_similarity(goid, goid, godag) == 1.0


if __name__ == '__main__':
    test_semantic_i150()
