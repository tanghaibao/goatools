#!/usr/bin/env python
"""Test initializing TermCounts with annotations made to alternate GO ID"""


import os
import collections as cx
from goatools.obo_parser import GODag
from goatools.associations import get_b2aset
from goatools.anno.idtogos_reader import IdToGosReader
from goatools.semantic import TermCounts


REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

# The number of unique genes annotated to a GO Term
GOID2NUM = {
    'GO:000000A': 20,
    'GO:000000B': 10,
    'GO:000000C': 10,
    'GO:000000D': 10,
    'GO:000000E': 10,
    'GO:000000F': 30,
    'GO:000000G': 10,
}

def test_semantic_similarity():
    """Test initializing TermCounts with annotations made to alternate GO ID"""
    godag = GODag(os.path.join(REPO, '../goatools/tests/data/yangRWC/fig2a.obo'))
    assoc = _get_id2gos(os.path.join(REPO, '../goatools/tests/data/yangRWC/fig2a.anno'), godag)
    tcntobj = TermCounts(godag, assoc)
    # N_v: Test accuracy of Python equivalent to Java: getNumberOfAnnotations
    # Test number of unique genes annotated to a GO Term PLUS genes annotated to a descendant
    assert tcntobj.gocnts['GO:000000A'] == 100, tcntobj.gocnts
    assert tcntobj.gocnts['GO:000000B'] == 40, tcntobj.gocnts
    assert tcntobj.gocnts['GO:000000C'] == 50, tcntobj.gocnts
    assert tcntobj.gocnts['GO:000000D'] == 10, tcntobj.gocnts
    assert tcntobj.gocnts['GO:000000E'] == 10, tcntobj.gocnts
    assert tcntobj.gocnts['GO:000000F'] == 30, tcntobj.gocnts
    assert tcntobj.gocnts['GO:000000G'] == 10, tcntobj.gocnts

def _get_id2gos(file_id2gos, godag):
    """Get annotations"""
    if os.path.exists(file_id2gos):
        return IdToGosReader(file_id2gos, godag=godag).get_id2gos('CC')
    go2genes = cx.defaultdict(set)
    genenum = 0
    for goid, qty in GOID2NUM.items():
        for _ in range(qty):
            go2genes[goid].add(genenum)
            genenum += 1
    id2gos = get_b2aset(go2genes)
    IdToGosReader.wr_id2gos(file_id2gos, id2gos)
    return id2gos


if __name__ == '__main__':
    test_semantic_similarity()
