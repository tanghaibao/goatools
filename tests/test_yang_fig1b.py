#!/usr/bin/env python
"""Test Yang's RWC measure added to other semantic similariy measures"""


import os
import collections as cx
from goatools.obo_parser import GODag
from goatools.utils import get_b2aset
from goatools.anno.idtogos_reader import IdToGosReader
from goatools.semantic import TermCounts


REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def test_semantic_similarity():
    """Test faster version of sematic similarity"""
    godag = GODag(os.path.join(REPO, 'tests/data/yangRWC/fig1b.obo'))
    name2go = {o.name: o.item_id for o in godag.values()}
    assoc = _get_id2gos(os.path.join(REPO, 'tests/data/yangRWC/fig1b.anno'), godag, name2go)
    tcntobj = TermCounts(godag, assoc)
    assert tcntobj.gocnts[name2go['I']] == 20
    assert tcntobj.gocnts[name2go['L']] == 20
    assert tcntobj.gocnts[name2go['M']] == 20
    assert tcntobj.gocnts[name2go['N']] == 20

def _get_id2gos(file_id2gos, godag, name2go):
    """Get annotations"""
    if os.path.exists(file_id2gos):
        return IdToGosReader(file_id2gos, godag=godag).get_id2gos('CC')
    id2num = {
        name2go['A']:  1,
        name2go['B']:  1,
        name2go['C']: 10,
        name2go['D']: 10,
        name2go['E']: 10,
        name2go['F']: 10,
        name2go['G']: 10,
        name2go['H']: 10,
        name2go['I']: 18,
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
