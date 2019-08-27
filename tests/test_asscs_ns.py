#!/usr/bin/env python
"""Test TermCounts object used in Resnik and Lin similarity calculations."""

from __future__ import print_function

import os
import sys
import timeit
import datetime
import collections as cx
from goatools.base import get_godag
from goatools.test_data.gafs import ASSOCIATIONS
from goatools.associations import dnld_annotation
from goatools.anno.gaf_reader import GafReader
from goatools.godag.consts import NS2NAMESPACE
from goatools.godag.consts import NAMESPACE2NS
# from goatools.godag.consts import NAMESPACE2GO

TIC = timeit.default_timer()
REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def test_semantic_similarity(usr_assc=None):
    """Computing basic semantic similarities between GO terms."""
    not_these = {'goa_uniprot_all.gaf', 'goa_uniprot_all_noiea.gaf'}
    assc_names = sorted(ASSOCIATIONS.difference(not_these))
    go2obj = get_go2obj()
    # http://current.geneontology.org/annotations/
    if usr_assc is not None:
        assc_names = [usr_assc]
    not_found = set()
    gaf2errs = cx.defaultdict(list)
    for assc_name in assc_names:  # Limit test numbers for speed
        tic = timeit.default_timer()
        # Get all the annotations from arabidopsis.
        fin_gaf = os.path.join(REPO, assc_name)
        if not os.path.exists(fin_gaf):
            dnld_annotation(fin_gaf)
        annoobj = GafReader(fin_gaf)
        for nta in annoobj.associations:
            if nta.GO_ID in go2obj:
                goterm = go2obj[nta.GO_ID]
                namespace_anno = NS2NAMESPACE.get(nta.NS)
                if namespace_anno != goterm.namespace:
                    gaf2errs[assc_name].append(nta)
            else:
                not_found.add(nta.GO_ID)
    print('{HMS} {N} Associations'.format(HMS=_hms(TIC), N=len(assc_names)))
    if not_found:
        _prt_not_found(not_found)
    if gaf2errs:
        _wr_errs('namespace_errors.txt', gaf2errs, go2obj)

def _wr_errs(fout_err, gaf2errs, go2obj):
    """Write errors in namespaces seen in annotation files"""
    with open(fout_err, 'w') as prt:
        err_cnt = 0
        gaf_errs = sorted(gaf2errs.items(), key=lambda t: len(t[1]))
        for gaf, errs in gaf_errs:
            err_cnt += len(errs)
            msg = '{N} mismarked namespaces in {GAF}'.format(GAF=gaf, N=len(errs))
            print(msg)
            prt.write('\n{TITLE}\n'.format(TITLE=msg))
            for nta in errs:
                prt.write('\n{GO} ACTUAL({ns}) EXPECTED({NS}) {GAF}:\n'.format(
                    GO=nta.GO_ID, ns=NS2NAMESPACE[nta.NS], NS=go2obj[nta.GO_ID].namespace, GAF=gaf))
                for fld, val in nta._asdict().items():
                    prt.write('{FLD:20}: {VAL}\n'.format(FLD=fld, VAL=val))

        print('  {N} GAFs WITH {E} TOTAL ERRORS WROTE: {TXT}'.format(
            N=len(gaf_errs), E=err_cnt, TXT=fout_err))


def _prt_not_found(not_found):
    print('**WARNING: {N} EMPTY ASSOCIATIONS:'.format(N=len(not_found)))
    for idx, assc in enumerate(not_found):
        print('    {I}) {ASSC}'.format(I=idx, ASSC=assc))

def _hms(tic):
    """Get Timing."""
    return '{HMS}'.format(HMS=str(datetime.timedelta(seconds=(timeit.default_timer()-tic))))


def get_go2obj():
    """Read GODag and return go2obj."""
    godag = get_godag(os.path.join(REPO, "go-basic.obo"), loading_bar=None)
    return {go:o for go, o in godag.items() if not o.is_obsolete}

if __name__ == '__main__':
    ASSC_NAME = None if len(sys.argv) == 1 else sys.argv[1]
    test_semantic_similarity(ASSC_NAME)
