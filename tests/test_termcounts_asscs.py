#!/usr/bin/env python
"""Test TermCounts object used in Resnik and Lin similarity calculations."""

from __future__ import print_function

import os
import sys
import timeit
import datetime
from goatools.base import get_godag
from goatools.semantic import TermCounts
from goatools.semantic import get_info_content
from goatools.test_data.gafs import ASSOCIATIONS
from goatools.associations import dnld_annotation
from goatools.anno.gaf_reader import GafReader
from goatools.godag.consts import NS2NAMESPACE

TIC = timeit.default_timer()
REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def test_semantic_similarity(usr_assc=None):
    """Computing basic semantic similarities between GO terms."""
    not_these = {'goa_uniprot_all.gaf', 'goa_uniprot_all_noiea.gaf'}
    associations = sorted(ASSOCIATIONS.difference(not_these))
    go2obj = get_go2obj()
    # goids = go2obj.keys()
    # http://current.geneontology.org/annotations/
    if usr_assc is not None:
        associations = [usr_assc]
    not_found = set()
    errs = []
    for assc_name in associations:  # Limit test numbers for speed
        tic = timeit.default_timer()
        # Get all the annotations from arabidopsis.
        fin_gaf = os.path.join(REPO, assc_name)
        if not os.path.exists(fin_gaf):
            dnld_annotation(fin_gaf)
        annoobj = GafReader(fin_gaf)
        #### for nspc in ['BP', 'MF', 'CC']:
        assc_gene2gos = annoobj.get_id2gos('all')
        if not assc_gene2gos:
            not_found.add(assc_name)
            continue

        # Calculate the information content of the single term, GO:0048364
        #       "Information content (GO:0048364) = 7.75481392334

        # Initialize the counts of each GO term.
        tcntobj = TermCounts(go2obj, assc_gene2gos)
        go_cnt = tcntobj.gocnts.most_common()

        #print tcntobj.gocnts.most_common()

        if go_cnt:
            print("{ASSC}".format(ASSC=assc_name))
            print(tcntobj.aspect_counts)
            gocnt_max = go_cnt[0][1]
            prt_info(tcntobj, go_cnt, None)
            prt_info(tcntobj, go_cnt, gocnt_max/2.0)
            prt_info(tcntobj, go_cnt, gocnt_max/10.0)
        print("{HMS} {hms} {ASSC}\n".format(ASSC=assc_name, HMS=_hms(TIC), hms=_hms(tic)))
    print('{HMS} {N} Associations'.format(HMS=_hms(TIC), N=len(associations)))
    if not_found:
        _prt_not_found(not_found)
    if errs:
        fout_err = 'namespace_errors.txt'
        with open(fout_err, 'w') as prt:
            for err in errs:
                prt.write(err)
            print('  {N} ERRORS WROTE: {TXT}'.format(N=len(errs), TXT=fout_err))


def _prt_not_found(not_found):
    print('**WARNING: {N} EMPTY ASSOCIATIONS:'.format(N=len(not_found)))
    for idx, assc in enumerate(not_found):
        print('    {I}) {ASSC}'.format(I=idx, ASSC=assc))

def _hms(tic):
    """Get Timing."""
    return '{HMS}'.format(HMS=str(datetime.timedelta(seconds=(timeit.default_timer()-tic))))

def prt_info(tcntobj, go_cnt, max_val):
    """Print the information content of a frequently used GO ID."""
    go_id, cnt = get_goid(go_cnt, max_val)
    infocontent = get_info_content(go_id, tcntobj)
    msg = 'Information content ({GO} {CNT:7,}) = {INFO:8.6f} {NAME}'
    print(msg.format(GO=go_id, CNT=cnt, INFO=infocontent, NAME=tcntobj.go2obj[go_id].name))

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
    godag = get_godag(os.path.join(REPO, "go-basic.obo"), loading_bar=None)
    return {go:o for go, o in godag.items() if not o.is_obsolete}

if __name__ == '__main__':
    ASSC_NAME = None if len(sys.argv) == 1 else sys.argv[1]
    test_semantic_similarity(ASSC_NAME)
