#!/usr/bin/env python
"""Test Yang's RWC measure added to other semantic similariy measures"""


import os
import collections as cx
from goatools.obo_parser import GODag
from goatools.associations import get_b2aset
from goatools.anno.idtogos_reader import IdToGosReader
from goatools.semantic import TermCounts


REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def test_semantic_similarity():
    """Test faster version of sematic similarity"""
    godag = GODag(os.path.join(REPO, 'tests/data/yangRWC/fig1b.obo'))
    assoc = _get_id2gos(os.path.join(REPO, 'tests/data/yangRWC/fig1b.anno'), godag)
    tcntobj = TermCounts(godag, assoc)
    assert tcntobj.gocnts['GO:000000I'] == 20
    assert tcntobj.gocnts['GO:000000L'] == 20
    assert tcntobj.gocnts['GO:000000M'] == 20
    assert tcntobj.gocnts['GO:000000N'] == 20

def _get_id2gos(file_id2gos, godag):
    """Get annotations"""
    if os.path.exists(file_id2gos):
        return IdToGosReader(file_id2gos, godag=godag).get_id2gos('CC')
    id2num = {
        'GO:000000A':  1,
        'GO:000000B':  1,
        'GO:000000C': 10,
        'GO:000000D': 10,
        'GO:000000E': 10,
        'GO:000000F': 10,
        'GO:000000G': 10,
        'GO:000000H': 10,
        'GO:000000I': 18,
    }
    go2genes = cx.defaultdict(set)
    genenum = 0
    for goid, qty in id2num.items():
        for _ in range(qty):
            go2genes[goid].add(genenum)
            genenum += 1
    id2gos = get_b2aset(go2genes)
    IdToGosReader.wr_id2gos(file_id2gos, id2gos)
    return id2gos



if __name__ == '__main__':
    test_semantic_similarity()
