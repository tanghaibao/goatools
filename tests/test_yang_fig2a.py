#!/usr/bin/env python
"""Test initializing TermCounts with annotations made to alternate GO ID"""


import os
import sys
import collections as cx
from goatools.obo_parser import GODag
from goatools.utils import get_b2aset
from goatools.anno.idtogos_reader import IdToGosReader
from goatools.semantic import TermCounts


REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

# The number of unique genes annotated to a GO Term
NAME2NUM = {
    'A': 20,
    'B': 10,
    'C': 10,
    'D': 10,
    'E': 10,
    'F': 10,
    'G': 30,
}

def test_semantic_similarity():
    """Test initializing TermCounts with annotations made to alternate GO ID"""
    godag = GODag(os.path.join(REPO, '../goatools/tests/data/yangRWC/fig2a.obo'))
    file_id2gos = os.path.join(REPO, '../goatools/tests/data/yangRWC/fig2a.anno')
    name2go = {o.name: o.item_id for o in godag.values()}
    assoc = _get_id2gos(file_id2gos, godag, name2go, NAME2NUM)
    tcntobj = TermCounts(godag, assoc)
    # N_v: Test accuracy of Python equivalent to Java: getNumberOfAnnotations
    # Test number of unique genes annotated to a GO Term PLUS genes annotated to a descendant
    assert tcntobj.gocnts[name2go['A']] == 100, tcntobj.gocnts
    assert tcntobj.gocnts[name2go['B']] == 40, tcntobj.gocnts
    assert tcntobj.gocnts[name2go['C']] == 50, tcntobj.gocnts
    assert tcntobj.gocnts[name2go['D']] == 10, tcntobj.gocnts
    assert tcntobj.gocnts[name2go['E']] == 10, tcntobj.gocnts
    assert tcntobj.gocnts[name2go['F']] == 10, tcntobj.gocnts
    assert tcntobj.gocnts[name2go['G']] == 30, tcntobj.gocnts

def _get_id2gos(file_id2gos, godag, name2go, name2num):
    """Get annotations"""
    if os.path.exists(file_id2gos):
        return IdToGosReader(file_id2gos, godag=godag).get_id2gos('CC')
    go2genes = cx.defaultdict(set)
    genenum = 0
    for name, qty in name2num.items():
        goid = name2go[name]
        for _ in range(qty):
            go2genes[goid].add(genenum)
            genenum += 1
    id2gos = get_b2aset(go2genes)
    IdToGosReader.wr_id2gos(file_id2gos, id2gos)
    return id2gos

def gen_anno_small():
    """Generate a maller nnotations containing 10% of the oringal genes"""
    godag = GODag(os.path.join(REPO, '../goatools/tests/data/yangRWC/fig2a.obo'))
    name2go = {o.name: o.item_id for o in godag.values()}
    file_id2gos = os.path.join(REPO, '../goatools/tests/data/yangRWC/fig2a_small.anno')
    name2num = {e:i/10 for e, i in NAME2NUM.items()}
    _get_id2gos(file_id2gos, godag, name2go, name2num)
    print(name2num)


if __name__ == '__main__':
    if len(sys.argv) == 1:
        test_semantic_similarity()
    else:
        gen_anno_small()
