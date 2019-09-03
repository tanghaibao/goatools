#!/usr/bin/env python
"""Test faster version of sematic similarity"""

from __future__ import print_function

# Computing basic semantic similarities between GO terms

# Adapted from book chapter written by _Alex Warwick Vesztrocy and Christophe Dessimoz_

# How to compute semantic similarity between GO terms.

# First we need to write a function that calculates the minimum number
# of branches connecting two GO terms.

import os
import timeit
from collections import Counter
## from goatools.base import get_godag
## from goatools.associations import dnld_assc
## from goatools.semantic import semantic_similarity
## from goatools.semantic import TermCounts
## from goatools.semantic import get_info_content
## from goatools.semantic import deepest_common_ancestor
## from goatools.semantic import resnik_sim
## from goatools.semantic import lin_sim
## from goatools.godag.consts import NS2GO

from goatools.anno.gpad_reader import GpadReader
from goatools.semantic import TermCounts

from tests.utils import get_godag
from tests.utils import get_anno_fullname
from tests.utils import prt_hms

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def test_semantic_similarity():
    """Test faster version of sematic similarity"""
    godag_r0 = get_godag('go-basic.obo')
    ## godag_r1 = get_godag('go-basic.obo', optional_attrs=['relationship'])
    annoobj = GpadReader(get_anno_fullname('goa_human.gpad'), godag=godag_r0)
    ns2assoc = annoobj.get_ns2assc()
    assoc = annoobj.get_id2gos('all')

    # Get TermCounts for each namespace and for all namespaces
    ns2tcnt = {ns:TermCounts(godag_r0, ns2assoc[ns]) for ns in ['BP', 'MF', 'CC']}
    tic = timeit.default_timer()
    tcntobj = TermCounts(godag_r0, assoc)
    prt_hms(tic, 'CUR ACTUAL   {N:,} TermCounts initialized'.format(N=len(tcntobj.gocnts)))
    # Compare various TermCount counts
    for nspc in ['BP', 'MF', 'CC']:
        for goid, cnt in ns2tcnt[nspc].gocnts.items():
            assert tcntobj.gocnts[goid] == cnt

    # Compare old and new count
    tic = timeit.default_timer()
    gocnts_old = _old_init_count_terms(godag_r0, assoc.values())
    assert gocnts_old
    prt_hms(tic, 'OLD EXPECTED {N:,} TermCounts initialized'.format(N=len(gocnts_old)))
    for goid, cnt_old in gocnts_old.items():
        assert cnt_old == tcntobj.gocnts[goid]


def _old_init_count_terms(go2obj, annots_values):
    '''
        Fills in the counts and overall aspect counts.
    '''
    gocnts = Counter()
    gonotindag = set()
    # Fill gocnts with GO IDs in annotations and their corresponding counts
    for terms in annots_values: # key is 'gene'
        # Make a union of all the terms for a gene, if term parents are
        # propagated but they won't get double-counted for the gene
        allterms = set()
        for go_id in terms:
            goobj = go2obj.get(go_id, None)
            if goobj is not None:
                allterms.add(go_id)
                allterms |= goobj.get_all_parents()
            else:
                gonotindag.add(go_id)
        # Add 1 for each GO annotated to this gene product
        for parent in allterms:
            gocnts[parent] += 1
    if gonotindag:
        print("{N} Assc. GO IDs not found in the GODag\n".format(N=len(gonotindag)))
    return gocnts

if __name__ == '__main__':
    test_semantic_similarity()
