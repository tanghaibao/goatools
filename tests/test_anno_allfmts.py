#!/usr/bin/env python
"""Test all annotation formats"""

import os
from goatools.base import get_godag
from goatools.associations import dnld_annofile
from goatools.anno.opts import AnnoOptions
from goatools.anno.factory import get_objanno
from goatools.evidence_codes import EvidenceCodes

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
INC_GOOD = set(EvidenceCodes.code2nt.keys()).difference({'IEA'})


def test_anno_read():
    """Test all annotation formats"""
    godag = get_godag(os.path.join(REPO, 'go-basic.obo'))

    # pylint: disable=superfluous-parens
    print('- DOWNLOAD AND LOAD -----------------------------------------------')
    annoobjs = [
        _get_objanno('gene2go', taxid=10090),
        _get_objanno('goa_human.gaf'),
        _get_objanno('goa_human.gpad'),
        _get_objanno('goa_human.gpad', godag=godag),
        _get_objanno('data/association', 'id2gos'),
        _get_objanno('data/association', 'id2gos', godag=godag),
    ]

    print('- prt_summary_anno2ev ---------------------------------------------')
    for idx, obj in enumerate(annoobjs):
        print('>>>>> {I} >>>>> prt_summary_anno2ev {ANNO}'.format(I=idx, ANNO=obj.name))
        obj.prt_summary_anno2ev()
        obj.chk_associations()
        if obj.name in {'gaf', 'gpad'}:
            _prt_fld(obj, 'Extension')

    print('- get_id2gos ------------------------------------------------------')
    for idx, obj in enumerate(annoobjs):
        id2gos = obj.get_id2gos()
        assert obj.get_id2gos(ev_include=INC_GOOD) == obj.get_id2gos(ev_exclude={'IEA'}), \
            'INC ALL({A}) != EXC IEA({I}): {DIF}'.format(
                A=len(obj.get_id2gos(ev_incclude=INC_GOOD)),
                I=len(obj.get_id2gos(ev_exclude={'IEA'})),
                # DIF=set(obj.get_id2gos(ev_inc=INC_GOOD).keys()).symmetric_difference(obj.get_id2gos(ev_exclude={'IEA'})))
                DIF='')
        num_ids = len(id2gos)
        print('>>>>> {I} >>>>> get_id2gos {N:6,} {ANNO}'.format(I=idx, N=num_ids, ANNO=obj.name))
        assert next(iter(next(iter(id2gos.values()))))[:3] == "GO:"
        if obj.filename[-16:] == 'data/association':
            assert num_ids == 34284

    print('- get_ns2... ------------------------------------------------------')
    for idx, obj in enumerate(annoobjs):
        if obj.name in {'gpad', 'id2gos'} and obj.godag is None:
            print('{IDX}) SKIPPING(No NS): {C}:get_ns2ntsanno'.format(IDX=idx, C=obj.__class__.__name__))
            continue
        _tst_ns2(obj, idx)


def _tst_ns2(obj, idx):
    """Test functions which use ns2 functions."""
    # ALL annotations for a species
    nts_all = obj.get_associations()
    num_nts_all = len(nts_all)
    num_nts_act = 0

    # Separate ALL annotations into BP MF CC
    ns2ntsanno = obj.get_ns2ntsanno()
    assert set(ns2ntsanno.keys()) == {'BP', 'MF', 'CC'}

    # Reduce annotations to remove IEA
    ns2anno_exp = {}
    kws = {'ev_include': INC_GOOD}
    for nspc, nts_orig in ns2ntsanno.items():
        nts_redu = obj.reduce_annotations(nts_orig, AnnoOptions(obj.evobj, **kws))
        num_nts_orig = len(nts_orig)
        # Check that only current namespace is seen on namedtuples
        assert set(nt.NS for nt in nts_orig if nt.NS == nspc) == {nspc}
        num_nts_act += num_nts_orig
        num_nts_redu = len(nts_redu)
        print('{IDX}) ns2ntanno {NS} {ALL:7,} -> {N:7,} -> {R:7,} annos: {TYPE}'.format(
            IDX=idx, NS=nspc, ALL=num_nts_all, N=num_nts_orig, R=num_nts_redu, TYPE=obj.name))
        ns2anno_exp[nspc] = obj.get_dbid2goids(nts_redu)

        assert nts_all >= num_nts_orig
        if obj.name == 'id2gos':
            assert num_nts_orig == num_nts_redu
        else:
            assert num_nts_orig > num_nts_redu

    assert num_nts_all >= num_nts_act

    # Compare step-by-step transformation with top-level function, id2gos
    ns2anno_act = obj.get_ns2assc(**kws)
    for nspc, anno_exp in ns2anno_exp.items():
        anno_act = ns2anno_act[nspc]
        assert set(anno_exp.keys()) == set(anno_act.keys())
        for geneid, gos_exp in anno_exp.items():
            gos_act = anno_act[geneid]
            assert gos_exp == gos_act

def _prt_fld(obj, fld):
    """Print examples of field values"""
    idx = 0
    for ntd in obj.associations:
        val = getattr(ntd, fld)
        if val:
            print(obj.name, str(val))
            idx += 1
            if idx == 10:
                break

def _get_objanno(fin_anno, anno_type=None, **kws):
    """Get association object"""
    full_anno = os.path.join(REPO, fin_anno)
    dnld_annofile(full_anno, anno_type)
    obj = get_objanno(full_anno, anno_type, **kws)
    if obj.name in {'gene2go', 'gaf'} or 'godag' in kws:
        assert hasattr(obj.associations[0], 'NS')
    return obj


if __name__ == '__main__':
    test_anno_read()
