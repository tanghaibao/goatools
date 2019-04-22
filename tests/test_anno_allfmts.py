#!/usr/bin/env python
"""Test all annotation formats"""

import os
from goatools.associations import dnld_annofile
from goatools.anno.factory import get_objanno

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def test_anno_read():
    """Test all annotation formats"""

    # pylint: disable=superfluous-parens
    print('- DOWNLOAD AND LOAD -----------------------------------------------')
    annoobjs = [
        _get_objanno('gene2go', taxid=10090),
        _get_objanno('goa_human.gaf'),
        _get_objanno('goa_human.gpad'),
        _get_objanno('data/association', anno_type='id2gos'),
    ]

    print('- prt_summary_anno2ev ---------------------------------------------')
    for idx, obj in enumerate(annoobjs):
        print('>>>>> {I} >>>>> prt_summary_anno2ev {ANNO}'.format(I=idx, ANNO=obj.name))
        obj.prt_summary_anno2ev()

    print('- get_id2gos ------------------------------------------------------')
    for idx, obj in enumerate(annoobjs):
        id2gos = obj.get_id2gos()
        num_ids = len(id2gos)
        print('>>>>> {I} >>>>> get_id2gos {N:6,} {ANNO}'.format(I=idx, N=num_ids, ANNO=obj.name))
        assert next(iter(next(iter(id2gos.values()))))[:3] == "GO:"
        if obj.filename[-16:] == 'data/association':
            assert num_ids == 34284


def _get_objanno(fin_anno, anno_type=None, **kws):
    """Get association object"""
    full_anno = os.path.join(REPO, fin_anno)
    dnld_annofile(full_anno, anno_type)
    obj = get_objanno(full_anno, anno_type, **kws)
    return obj


if __name__ == '__main__':
    test_anno_read()
