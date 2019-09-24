#!/usr/bin/env python
"""Test TermCounts object used in Resnik and Lin similarity calculations."""

from __future__ import print_function

import os
import sys
import timeit
import datetime
from goatools.test_data.gafs import ASSOCIATIONS
from goatools.associations import dnld_annofile
from goatools.anno.factory import get_objanno

TIC = timeit.default_timer()
REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def test_gafs(usr_assc=None):
    """Computing basic semantic similarities between GO terms."""
    not_these = {'goa_uniprot_all.gaf', 'goa_uniprot_all_noiea.gaf'}
    associations = sorted(ASSOCIATIONS.difference(not_these))
    # http://current.geneontology.org/annotations/
    if usr_assc is not None:
        associations = [usr_assc]

    failed = []
    not_found = set()
    ## https://github.com/geneontology/go-annotation/issues/2659
    ## To get the line numbers containing the errors:
    ##     assc_name example:
    ##         goa_human.gaf
    ##         goa_chicken.gaf
    ##         goa_cow.gaf
    ## only = {'gramene_oryza.gaf', 'pamgo_mgrisea.gaf', 'tair.gaf'}
    for assc_name in associations:  # Limit test numbers for speed
        ## if assc_name not in only:
        ##     continue
        obj = _get_objanno(assc_name, 'gaf')
        if not obj.chk_associations('{BASE}.err'.format(BASE=assc_name)):
            failed.append(assc_name)

    print('{HMS} {N} Associations'.format(HMS=_hms(TIC), N=len(associations)))
    if not_found:
        _prt_not_found(not_found)
    for assc in failed:
        print('**FAILED: {A}'.format(A=assc))
    # assert not failed

def _prt_not_found(not_found):
    print('**WARNING: {N} EMPTY ASSOCIATIONS:'.format(N=len(not_found)))
    for idx, assc in enumerate(not_found):
        print('    {I}) {ASSC}'.format(I=idx, ASSC=assc))

def _hms(tic):
    """Get Timing."""
    return '{HMS}'.format(HMS=str(datetime.timedelta(seconds=(timeit.default_timer()-tic))))

def _get_objanno(fin_anno, anno_type=None, **kws):
    """Get association object"""
    full_anno = os.path.join(REPO, fin_anno)
    dnld_annofile(full_anno, anno_type)
    obj = get_objanno(full_anno, anno_type, **kws)
    return obj

if __name__ == '__main__':
    ASSC_NAME = None if len(sys.argv) == 1 else sys.argv[1]
    test_gafs(ASSC_NAME)
