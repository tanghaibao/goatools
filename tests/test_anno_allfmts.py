#!/usr/bin/env python
"""Test all annotation formats"""

import os
from goatools.associations import dnld_ncbi_gene_file
from goatools.associations import dnld_assc
from goatools.anno.factory import get_objanno
from goatools.anno.factory import get_anno_desc
from goatools.base import download_ncbi_associations

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def test_anno_read():
    """Test all annotation formats"""
    objanno = _Anno()

    # pylint: disable=superfluous-parens
    print('- DOWNLOAD AND LOAD -----------------------------------------------')
    annoobjs = [
        objanno.get_objanno('gene2go', taxid=10090),
        objanno.get_objanno('goa_human.gaf'),
        objanno.get_objanno('goa_human.gpad'),
        objanno.get_objanno('data/association', anno_type='id2gos'),
    ]

    print('- prt_summary_anno2ev ---------------------------------------------')
    for idx, obj in enumerate(annoobjs):
        print('>>>>> {I} >>>>> prt_summary_anno2ev {ANNO}'.format(I=idx, ANNO=obj.name))
        obj.prt_summary_anno2ev()


# pylint: disable=too-few-public-methods
class _Anno(object):
    """Holds annotation filenames, downloads files if necessary"""

    def get_objanno(self, fin_anno, anno_type=None, **kws):
        """Get association object"""
        full_anno = os.path.join(REPO, fin_anno)
        self._dnld(full_anno, anno_type)
        obj = get_objanno(full_anno, anno_type, **kws)
        return obj

    @staticmethod
    def _dnld(fin_anno, anno_type):
        """Download annotation file, if needed"""
        if os.path.exists(fin_anno):
            return
        anno_type = get_anno_desc(fin_anno, anno_type)
        if anno_type == 'gene2go':
            download_ncbi_associations(fin_anno)
        if anno_type in {'gaf', 'gpad'}:
            dnld_assc(fin_anno)

    @staticmethod
    def _dnld_gene2go(file_anno):
        """Download the annotation file, if needed."""
        if os.path.exists(file_anno):
            assert os.path.getsize(file_anno) > 1000000, "BAD ANNO({F})".format(F=file_anno)
            return
        dnld_ncbi_gene_file(file_anno, loading_bar=None)
        assert os.path.isfile(file_anno), "MISSING ANNO({F})".format(F=file_anno)
        assert os.path.getsize(file_anno) > 1000000, "BAD ANNO({F})".format(F=file_anno)


if __name__ == '__main__':
    test_anno_read()
