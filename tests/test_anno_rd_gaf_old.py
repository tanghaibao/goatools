#!/usr/bin/env python
"""Ensure NEW results are equal to OLD results: read_ncbi_gene2go."""

from __future__ import print_function

import os
import collections as cx
from goatools.anno.gaf_reader import GafReader
from goatools.associations import dnld_annofile


REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def test_anno_read():
    """Test reading annotation file."""
    fin_anno = os.path.join(REPO, 'goa_human.gaf')
    dnld_annofile(fin_anno, 'gaf')

    print('\nTEST STORING ONLY ONE SPECIES')
    obj = GafReader(fin_anno)
    obj.prt_summary_anno2ev()
    ## new = obj.read_gaf()
    new = obj.get_id2gos_nss()
    old = read_gaf(obj)
    _prt_differences(new, old, obj)
    print('{N} NEW'.format(N=len(new)))
    print('{N} OLD'.format(N=len(old)))
    assert len(new) == len(old), 'new({N}) != old({O})'.format(N=len(new), O=len(old))

    print('\nTEST KWS: keep_ND and keep_NOT')
    # pylint: disable=bad-whitespace
    kws_lst = [
        {'keep_ND': False, 'keep_NOT': False},
        {'keep_ND': False, 'keep_NOT': True},
        {'keep_ND': True,  'keep_NOT': False},
        {'keep_ND': True,  'keep_NOT': True},
    ]
    for kws in kws_lst:
        print('\nTEST KWS:', kws)
        ## new = obj.read_gaf(namespace='BP', **kws)
        new = obj.get_id2gos_nss(**kws)
        old = read_gaf(obj, **kws)
        _prt_differences(new, old, obj)
        assert len(new) == len(old), 'new({N}) != old({O})'.format(N=len(new), O=len(old))

    print('\nTEST GETTING REVERSE ASSOCIATIONS: GO2GENES')
    ## new = obj.read_gaf(go2geneids=True)
    new = obj.get_id2gos_nss(go2geneids=True)
    old = read_gaf(obj, go2geneids=True)
    _prt_differences(new, old, obj)
    assert len(new) == len(old), 'new({N}) != old({O})'.format(N=len(new), O=len(old))

    print('\nTEST RETURNING ASSOCIATIONS FOR SELECTED EVIDENCE CODES')
    evcodes = set(['ISO', 'IKR'])
    print("\nTEST 9606 ev_include={CODES}".format(CODES=' '.join(evcodes)))
    ## new = obj.read_gaf(ev_include=evcodes)
    new = obj.get_id2gos_nss(ev_include=evcodes)
    old = read_gaf(obj, ev_include=evcodes)
    _prt_differences(new, old, obj)
    assert new == old

def _prt_differences(new, old, obj):
    """Print differences between the new and old associations"""
    ids_diff = set(new.keys()).symmetric_difference(old.keys())
    print('{N} NEW'.format(N=len(new)))
    print('{N} OLD'.format(N=len(old)))
    for db_id in ids_diff:
        print('NEW: ', new.get(db_id))
        print('OLD: ', old.get(db_id))
        for nta in obj.associations:
            if nta.DB_ID == db_id:
                print(nta)

# Formerly was in goatools/anno/gaf_reader.py
def read_gaf(obj, **kws):
    """Read Gene Association File (GAF). Return data."""
    # Simple associations
    id2gos = cx.defaultdict(set)
    # keyword arguments for choosing which GO IDs to keep
    # Optional detailed associations split by taxid and having both ID2GOs & GO2IDs
    taxid2asscs = kws.get('taxid2asscs', None)
    b_geneid2gos = not kws.get('go2geneids', False)
    evs = kws.get('ev_include', None)
    eval_nd = _get_nd(kws.get('keep_ND', False))
    eval_not = _get_not(kws.get('keep_NOT', False))
    # Optionally specify a subset of GOs based on their evidence.
    # By default, return id2gos. User can cause go2geneids to be returned by:
    #   >>> read_ncbi_gene2go(..., go2geneids=True
    for ntgaf in obj.associations:
        if eval_nd(ntgaf) and eval_not(ntgaf):
            if evs is None or ntgaf.Evidence_Code in evs:
                geneid = ntgaf.DB_ID
                go_id = ntgaf.GO_ID
                if b_geneid2gos:
                    id2gos[geneid].add(go_id)
                else:
                    id2gos[go_id].add(geneid)
                if taxid2asscs is not None:
                    if ntgaf.Taxon:
                        taxid = ntgaf.Taxon[0]
                        taxid2asscs[taxid]['ID2GOs'][geneid].add(go_id)
                        taxid2asscs[taxid]['GO2IDs'][go_id].add(geneid)
    obj.hms('OLD SLOW read_gaf(kws={KWS})'.format(KWS=str(kws)))
    return dict(id2gos) # return simple associations

def _get_nd(keep_nd):
    """Allow GAF values always or never."""
    if keep_nd:
        return lambda nt: True
    return lambda nt: nt.Evidence_Code != 'ND'

def _get_not(keep_not):
    """Allow GAF values always or never."""
    if keep_not:
        return lambda nt: True
    return lambda nt: nt.Qualifier.isdisjoint({'NOT', 'not'})


if __name__ == '__main__':
    test_anno_read()
