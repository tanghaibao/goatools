#!/usr/bin/env python
"""Read GAF file and allow ND Evidence codes."""

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import os
import sys
import timeit
from goatools.anno.opts import AnnoOptions
from goatools.anno.factory import get_objanno
from goatools.associations import dnld_annofile

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def test_gaf_read():
    """Return GO associations from a GAF file. Download if necessary."""
    # On 2017/04/10, there were 3 GO IDs with ND Evidence Codes:
    #
    #    $ cut -f5,7 goa_human.gaf | grep ND | sort | uniq -c
    #        739 GO:0003674      ND
    #        484 GO:0005575      ND
    #        639 GO:0008150      ND

    # pylint: disable=superfluous-parens
    print('- DOWNLOAD AND LOAD -----------------------------------------------')
    annoobjs = [
        _Chk.get_objanno('gene2go', taxid=10090),
        _Chk.get_objanno('goa_human.gaf'),
        _Chk.get_objanno('goa_human.gpad'),
        _Chk.get_objanno('data/association', anno_type='id2gos'),
    ]
    for obj in annoobjs:
        _run(obj)

def _run(obj):
    """Test reducing annotations as specified by the user"""
    #pylint: disable=superfluous-parens
    print('\nTESTING: {T} --------------------------------------------------'.format(T=obj.name))
    chk = _Chk(obj)
    obj.chk_qualifiers()
    if obj.name != 'id2gos':
        assert chk.num_nots != 0, '{NAME} HAS 0 NOT Qualifiers'.format(NAME=obj.name)
    assc = obj.associations
    num_assc = len(assc)

    tic = timeit.default_timer()
    _ev = obj.evobj
    assc_dflt = obj.reduce_annotations(assc, AnnoOptions(_ev))
    tic = obj.hms('Default', tic)
    assc_nd0_not0 = obj.reduce_annotations(assc, AnnoOptions(_ev, keep_ND=False, keep_NOT=False))
    tic = obj.hms('{N:6,} assc_nd0_not0'.format(N=len(assc_nd0_not0)), tic)
    assc_nd0_not1 = obj.reduce_annotations(assc, AnnoOptions(_ev, keep_ND=False, keep_NOT=True))
    tic = obj.hms('{N:6,} assc_nd0_not1'.format(N=len(assc_nd0_not1)), tic)
    assc_nd1_not0 = obj.reduce_annotations(assc, AnnoOptions(_ev, keep_ND=True, keep_NOT=False))
    tic = obj.hms('{N:6,} assc_nd1_not0'.format(N=len(assc_nd1_not0)), tic)
    assc_nd1_not1 = obj.reduce_annotations(assc, AnnoOptions(_ev, keep_ND=True, keep_NOT=True))
    tic = obj.hms('{N:6,} assc_nd1_not1'.format(N=len(assc_nd1_not1)), tic)

    _red = obj.reduce_annotations
    inc = set(_ev.code2nt.keys()).difference({'IEA'})
    # pylint: disable=line-too-long
    assert _red(assc, AnnoOptions(obj.evobj, ev_exclude={'IEA'})) == _red(assc, AnnoOptions(obj.evobj, ev_include=inc))
    inc = _ev.grp2codes['High_Throughput']
    assert _red(assc, AnnoOptions(obj.evobj, ev_include={'High_Throughput'})) == _red(assc, AnnoOptions(obj.evobj, ev_include=inc))


    assert len(assc_dflt) == len(assc_nd0_not0)
    assert len(assc_nd1_not1) == num_assc
    # if obj.name in {'gaf', 'id2gos', 'gpad', 'id2go'}:
    assert len(assc_nd1_not0) == num_assc - chk.num_nots, '{N}: ACT({A:,}) != EXP({E:,})'.format(
        N=obj.name, A=len(assc_nd1_not0), E=num_assc - chk.num_nots)
    assert len(assc_nd0_not1) == num_assc - chk.num_nds, '{N}: ACT({A:,}) != EXP({E:,})'.format(
        N=obj.name, A=len(assc_nd0_not1), E=num_assc - chk.num_nds)

    # chk.prt_n_nts()


class _Chk(object):
    """Check results"""

    def __init__(self, obj):
        self.obj = obj
        self.nds = obj.nts_ev_nd()
        self.nots = obj.nts_qual_not()
        self.num_nds = len(self.nds)
        self.num_nots = len(self.nots)
        # pylint: disable=superfluous-parens
        print('{N} ND Evidence_Codes'.format(N=len(self.nds)))
        print('{N} NOT Qualifiers'.format(N=len(self.nots)))
        # self._prt_unusual_nots()

    def _prt_unusual_nots(self, prt=sys.stdout):
        """Print Qualifiers containint NOT"""
        for ntd in self.nots:
            for qual in ntd.Qualifier:
                if 'NOT' in qual.upper():
                    # if qual != 'not':  # and qual != 'NOT':
                    prt.write('{NT}\n'.format(NT=ntd))
                    # pass
                # if 'colocalizes_with' in qual:
                #     prt.write('{NT}\n'.format(NT=ntd))

    def prt_n_nts(self, num=10):
        """Print namedtuples"""
        cnt = 0
        for ntd in self.obj.associations:
            if len(ntd.Qualifier) > 2:
                # pylint: disable=superfluous-parens
                print(ntd)
                cnt += 1
                if cnt == num:
                    break

    @staticmethod
    def get_objanno(fin_anno, anno_type=None, **kws):
        """Get association object"""
        full_anno = os.path.join(REPO, fin_anno)
        dnld_annofile(full_anno, anno_type)
        obj = get_objanno(full_anno, anno_type, **kws)
        return obj


if __name__ == '__main__':
    test_gaf_read()

# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved.
