#!/usr/bin/env python
"""Test all annotation formats"""

import os
from goatools.base import get_godag
from goatools.associations import dnld_annofile
from goatools.utils import get_b2aset
from goatools.anno.opts import AnnoOptions
from goatools.anno.factory import get_objanno
from goatools.evidence_codes import EvidenceCodes

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
INC_GOOD = set(EvidenceCodes.code2nt.keys()).difference({'IEA'})


# pylint: disable=superfluous-parens
def test_anno_read():
    """Test all annotation formats"""
    godag = get_godag(os.path.join(REPO, 'go-basic.obo'))

    # pylint: disable=superfluous-parens
    print('- DOWNLOAD (if needed) AND LOAD -----------------------------------')
    annoobjs = [
        # gene2go
        _get_objanno('gene2go', taxid=10090),
        _get_objanno('gene2go', taxid=10090, namespaces={'BP'}),
        _get_objanno('gene2go', taxid=10090, namespaces={'MF'}),
        _get_objanno('gene2go', taxid=10090, namespaces={'CC'}),
        # gaf
        _get_objanno('goa_human.gaf'),
        _get_objanno('goa_human.gaf', namespaces={'BP'}),
        _get_objanno('goa_human.gaf', namespaces={'MF'}),
        _get_objanno('goa_human.gaf', namespaces={'CC'}),
        # gpad
        _get_objanno('goa_human.gpad', godag=godag),
        _get_objanno('goa_human.gpad', godag=godag, namespaces={'BP'}),
        _get_objanno('goa_human.gpad', godag=godag, namespaces={'MF'}),
        _get_objanno('goa_human.gpad', godag=godag, namespaces={'CC'}),
        _get_objanno('goa_human.gpad'),
        _get_objanno('goa_human.gpad', namespaces={'BP'}),
        _get_objanno('goa_human.gpad', namespaces={'MF'}),
        _get_objanno('goa_human.gpad', namespaces={'CC'}),
        # id2gos
        _get_objanno('data/association', 'id2gos'),
        _get_objanno('data/association', 'id2gos', namespaces={'BP'}),
        _get_objanno('data/association', 'id2gos', namespaces={'MF'}),
        _get_objanno('data/association', 'id2gos', namespaces={'CC'}),
        _get_objanno('data/association', 'id2gos', godag=godag),
        _get_objanno('data/association', 'id2gos', godag=godag, namespaces={'BP'}),
        _get_objanno('data/association', 'id2gos', godag=godag, namespaces={'MF'}),
        _get_objanno('data/association', 'id2gos', godag=godag, namespaces={'CC'}),
    ]

    print('- RUN get_id2gos --------------------------------------------------')
    _run_get_id2gos2(annoobjs)

    print('- RUN prt_summary_anno2ev -----------------------------------------')
    _run_prt_summary_anno2ev(annoobjs)

    print('- RUN print extension ---------------------------------------------')
    _run_print_extension(annoobjs)

    print('- RUN get_id2gos ALL ----------------------------------------------')
    _run_get_id2gos_nss(annoobjs)

    print('- RUN get_id2gos --------------------------------------------------')
    _run_get_id2gos(annoobjs)

    print('- RUN get_ns2... --------------------------------------------------')
    _run_get_ns(annoobjs)


def _run_prt_summary_anno2ev(annoobjs):
    """Test prt_summary_anno2ev"""
    for idx, obj in enumerate(annoobjs):
    # pylint: disable=superfluous-parens
        print('>>>>> {I} >>>>> prt_summary_anno2ev {ANNO}'.format(I=idx, ANNO=obj.get_desc()))
        obj.prt_summary_anno2ev()
        obj.chk_associations()

def _run_print_extension(annoobjs):
    """Test printing extensions"""
    for idx, obj in enumerate(annoobjs):
        print('>>>>> {I} >>>>> print Extension {ANNO}'.format(I=idx, ANNO=obj.get_desc()))
        if obj.name in {'gaf', 'gpad'}:
            _prt_fld(obj, 'Extension', 10)

def _run_get_id2gos_nss(annoobjs):
    """Test getting id2gos for all namespaces"""
    for idx, obj in enumerate(annoobjs):
        print('\n{I}) get_id2gos_nss {N:7,} annotations DESC({DESC}) NSs({NSs})'.format(
            I=idx, DESC=obj.get_desc(), NSs=obj.namespaces, N=len(obj.associations)))
        id2gos_nss = obj.get_id2gos_nss()  # Get all namespaces
        assert id2gos_nss

def _run_get_id2gos2(annoobjs):
    """Test get_id2gos"""
    idx = 0
    for obj in annoobjs:
        # Test only on annotations that loaded all: BP MF CC
        if not obj.namespaces and obj.has_ns():
            print('\n{I}) get_id2gos {DESC} {NSs} annoobj.namespaces({N:,}) annotations'.format(
                I=idx, DESC=obj.get_desc(), NSs=obj.namespaces, N=len(obj.associations)))
            id2gos = obj.get_id2gos('all')
            assert id2gos, 'NO ANNOTATIONS FOUND'
            assert id2gos == obj.get_id2gos_nss()
            # pylint: disable=line-too-long
            assert obj.get_id2gos('all', ev_include={'IEA'}) == obj.get_id2gos_nss(ev_include={'IEA'})
            idx += 1

def _run_get_id2gos(annoobjs):
    """Test get_id2gos"""
    for idx, obj in enumerate(annoobjs):
        print('\n{I}) get_id2gos {DESC} {NSs} {N:,} annotations'.format(
            I=idx, DESC=obj.get_desc(), NSs=obj.namespaces, N=len(obj.associations)))
        # If all namespaces are loaded, returns BP, else returns loaded NS
        print('Load all evidence codes')
        id2gos = obj.get_id2gos()
        assert id2gos, 'NO ANNOTATIONS FOUND'
        ## print(next(iter(obj.associations)))
        print('Load all evidence codes, except IEA')
        id2gos_inc = obj.get_id2gos(ev_include=INC_GOOD)
        id2gos_exc = obj.get_id2gos(ev_exclude={'IEA'})
        assert id2gos_exc, 'NO NON-IEA ANNOTATIONS FOUND'
        assert id2gos_inc == id2gos_exc, \
            'INC ALL({A}) != EXC IEA({I}): {DIF}'.format(
                A=len(id2gos_inc),
                I=len(id2gos_exc),
                # DIF=set(obj.get_id2gos(ev_inc=INC_GOOD).keys()).
                # symmetric_difference(obj.get_id2gos(ev_exclude={'IEA'})))
                DIF='')
        num_ids = len(id2gos)
        print('>>>>> {I} >>>>> get_id2gos {N:6,} go2ids[{B:6,}] {ANNO}'.format(
            I=idx, N=num_ids, ANNO=obj.get_desc(), B=len(get_b2aset(id2gos))))
        assert next(iter(next(iter(id2gos.values()))))[:3] == "GO:"
        if obj.filename[-16:] == 'data/association' and obj.godag is None:
            assert num_ids == 34284

def _run_get_ns(annoobjs):
    """Test getting namespace dict"""
    for idx, obj in enumerate(annoobjs):
        if obj.name in {'gpad', 'id2gos'} and obj.godag is None:
            print('{IDX}) SKIPPING(No NS): {C}:get_ns2ntsanno {ANNO}'.format(
                IDX=idx, C=obj.__class__.__name__, ANNO=obj.get_desc()))
            continue
        _tst_ns2(obj, idx)

# pylint: disable=too-many-locals
def _tst_ns2(obj, idx):
    """Test functions which use ns2 functions."""
    # ALL annotations for a species
    nts_all = obj.get_associations()
    num_nts_all = len(nts_all)
    num_nts_act = 0

    # Separate ALL annotations into BP MF CC
    ns2ntsanno = obj.get_ns2ntsanno()
    assert set(ns2ntsanno.keys()).issubset({'BP', 'MF', 'CC'}), ns2ntsanno.keys()

    # Reduce annotations to remove IEA
    ns2anno_exp = {}
    kws = {'ev_include': INC_GOOD}
    for nspc, nts_orig in sorted(ns2ntsanno.items()):
        opt = AnnoOptions(obj.evobj, **kws)
        nts_redu = obj.reduce_annotations(nts_orig, opt)
        num_nts_orig = len(nts_orig)
        # Check that only current namespace is seen on namedtuples
        assert set(nt.NS for nt in nts_orig if nt.NS == nspc) == {nspc}
        num_nts_act += num_nts_orig
        num_nts_redu = len(nts_redu)
        # pylint: disable=line-too-long
        print('{OPT} {IDX}) ns2ntanno {NS} {ALL:7,}=Loaded -> {N:7,} -> {R:7,} annos: {TYPE}'.format(
            OPT=opt, IDX=idx, NS=nspc, ALL=num_nts_all, N=num_nts_orig, R=num_nts_redu, TYPE=obj.get_desc()))
        ns2anno_exp[nspc] = obj.get_dbid2goids(nts_redu)

        assert num_nts_all >= num_nts_orig
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

def _prt_fld(obj, fld, max_num=10):
    """Print examples of field values"""
    idx = 0
    for ntd in obj.associations:
        val = getattr(ntd, fld)
        if val:
            print('PRINT THIS FIELD({F}) DESC({D}): VALUE({V})'.format(
                F=fld, D=obj.get_name(), V=str(val)))
            idx += 1
            if idx == max_num:
                break

def _get_objanno(fin_anno, anno_type=None, **kws):
    """Get association object"""
    full_anno = os.path.join(REPO, fin_anno)
    dnld_annofile(full_anno, anno_type)
    obj = get_objanno(full_anno, anno_type, **kws)
    if obj.name in {'gene2go', 'gaf'} or 'godag' in kws:
        assert hasattr(obj.associations[0], 'NS')
    # Check namespaces, if provided by the user
    if 'namespaces' in kws:
        _chk_namespaces(obj, kws['namespaces'])
    return obj

def _chk_namespaces(obj, namespaces):
    """Check namespaces"""
    # If namespaces is ignored due to godag=None, all namedtuples are loaded
    if obj.name in {'gpad', 'id2gos'} and obj.godag is None:
        assert not hasattr(next(iter(obj.associations)), 'NS')
        return
    for nta in obj.associations:
        assert nta.NS in namespaces, 'nta.NS({NS}) not in ({NSs})\n{NT}'.format(
            NSs=namespaces, NS=nta.NS, NT=nta)

def _get_num_annos(id2gos):
    """Return the number of annotations found in id2gos"""
    num = 0
    for gos in id2gos.values():
        num += len(gos)
    return num

if __name__ == '__main__':
    test_anno_read()
