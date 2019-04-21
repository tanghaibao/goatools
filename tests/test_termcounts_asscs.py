#!/usr/bin/env python
"""Test TermCounts object used in Resnik and Lin similarity calculations."""

from __future__ import print_function

import os
import sys
import timeit
import datetime
from goatools.base import get_godag
from goatools.associations import dnld_assc
from goatools.semantic import TermCounts
from goatools.semantic import get_info_content
from goatools.test_data.gafs import ASSOCIATIONS

TIC = timeit.default_timer()

def test_semantic_similarity(usr_assc=None):
    """Computing basic semantic similarities between GO terms."""
    not_these = {'goa_uniprot_all.gaf', 'goa_uniprot_all_noiea.gaf'}
    associations = sorted(ASSOCIATIONS.difference(not_these))
    go2obj = get_go2obj()
    # goids = go2obj.keys()
    # http://current.geneontology.org/annotations/
    if usr_assc is not None:
        associations = [usr_assc]
    cwd = os.getcwd()
    not_found = set()
    for assc_name in associations:  # Limit test numbers for speed
        tic = timeit.default_timer()
        # Get all the annotations from arabidopsis.
        assc_gene2gos = dnld_assc(os.path.join(cwd, assc_name), go2obj, prt=sys.stdout)
        if not assc_gene2gos:
            not_found.add(assc_name)
            continue

        # Calculate the information content of the single term, GO:0048364
        #       "Information content (GO:0048364) = 7.75481392334

        # First get the counts of each GO term.
        termcounts = TermCounts(go2obj, assc_gene2gos)
        go_cnt = termcounts.gocnts.most_common()
        #print termcounts.gocnts.most_common()

        if go_cnt:
            print("{ASSC}".format(ASSC=assc_name))
            print(sorted(termcounts.aspect_counts.most_common()))
            gocnt_max = go_cnt[0][1]
            prt_info(termcounts, go_cnt, None)
            prt_info(termcounts, go_cnt, gocnt_max/2.0)
            prt_info(termcounts, go_cnt, gocnt_max/10.0)
        print("{HMS} {hms} {ASSC}\n".format(ASSC=assc_name, HMS=_hms(TIC), hms=_hms(tic)))
    print('{HMS} {N} Associations'.format(HMS=_hms(TIC), N=len(associations)))
    if not_found:
        _prt_not_found(not_found)

def _prt_not_found(not_found):
    print('**WARNING: {N} EMPTY ASSOCIATIONS:'.format(N=len(not_found)))
    for idx, assc in enumerate(not_found):
        print('    {I}) {ASSC}'.format(I=idx, ASSC=assc))

def _hms(tic):
    """Get Timing."""
    return '{HMS}'.format(HMS=str(datetime.timedelta(seconds=(timeit.default_timer()-tic))))

def prt_info(termcounts, go_cnt, max_val):
    """Print the information content of a frequently used GO ID."""
    go_id, cnt = get_goid(go_cnt, max_val)
    infocontent = get_info_content(go_id, termcounts)
    msg = 'Information content ({GO} {CNT:7,}) = {INFO:8.6f} {NAME}'
    print(msg.format(GO=go_id, CNT=cnt, INFO=infocontent, NAME=termcounts.go2obj[go_id].name))

def get_goid(go_cnt, max_val):
    """Get frequently used GO ID."""
    if max_val is not None:
        for goid, cnt in go_cnt:
            if cnt < max_val:
                return goid, cnt
        return go_cnt[-1][0], go_cnt[-1][1]
    return go_cnt[0][0], go_cnt[0][1]

def get_go2obj():
    """Read GODag and return go2obj."""
    godag = get_godag(os.path.join(os.getcwd(), "go-basic.obo"), loading_bar=None)
    return {go:o for go, o in godag.items() if not o.is_obsolete}

if __name__ == '__main__':
    ASSC_NAME = None if len(sys.argv) == 1 else sys.argv[1]
    test_semantic_similarity(ASSC_NAME)
