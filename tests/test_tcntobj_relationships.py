#!/usr/bin/env python
"""Test loading of relationships, like part_of, into TermCounts"""

import os
import sys
from goatools.base import download_go_basic_obo
from goatools.obo_parser import GODag
from goatools.associations import dnld_annotation
from goatools.anno.gpad_reader import GpadReader
from goatools.semantic import TermCounts


REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
NSS = ['BP', 'MF', 'CC', 'all']
RELS = {'part_of',}

def test_tcntobj_relationships(prt=sys.stdout):
    """Test loading of relationships, like part_of, into TermCounts"""
    fin_obo = os.path.join(REPO, "go-basic.obo")
    fin_anno = os.path.join(REPO, 'goa_human.gpad')

    download_go_basic_obo(fin_obo, prt, loading_bar=None)
    dnld_annotation(fin_anno)

    # Load ontologies
    go2obj_r0 = GODag(fin_obo)
    go2obj_r1 = GODag(fin_obo, optional_attrs=['relationship'])

    # Load annotations
    annoobj = GpadReader(fin_anno, godag=go2obj_r0)

    # Create TermCounts objects
    ns2tcntobj_r0 = {ns:TermCounts(go2obj_r0, annoobj.get_id2gos(ns)) for ns in NSS}
    ns2tcntobj_r1 = {ns:TermCounts(go2obj_r1, annoobj.get_id2gos(ns), RELS) for ns in NSS}
    _chk_pass_fail(ns2tcntobj_r0, ns2tcntobj_r1)


def _chk_pass_fail(ns2tcntobj_r0, ns2tcntobj_r1):
    """Check to see that term counts are different w and wo/relationships"""
    print('\nCOMPARE GO Counts wo/relationships and with: {Rs}'.format(Rs=' '.join(sorted(RELS))))
    for nspc in NSS:
        cnt = 0
        go2cnts_r1 = ns2tcntobj_r1[nspc].gocnts
        for goid, cnt_r0 in ns2tcntobj_r0[nspc].gocnts.most_common():
            cnt_r1 = go2cnts_r1[goid]
            assert cnt_r0 <= cnt_r1
            if cnt_r0 != cnt_r1:
                cnt += 1
        print('{NS:3} {N:5,} more GO ID counts using relationships'.format(NS=nspc, N=cnt))
        assert ns2tcntobj_r0[nspc].gocnts != ns2tcntobj_r1[nspc].gocnts


if __name__ == '__main__':
    test_tcntobj_relationships()
