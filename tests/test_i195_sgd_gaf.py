#!/usr/bin/env python3
"""Investigate GAF reading error on saccharomyces"""

from os import system
from os.path import join
from os.path import exists
from os.path import basename
import wget

from tests.utils import REPO
from goatools.anno.gaf_reader import GafReader


def test_i195():
    """Investigate GAF reading error on saccharomyces"""
    fin_gaf1 = join(REPO, 'sgd.gaf')
    dnld_gaf1 = 'http://current.geneontology.org/annotations/sgd.gaf.gz'
    _dnld_gaf(fin_gaf1, dnld_gaf1)

    fin_gaf2 = join(REPO, 'gene_association.sgd.gaf')
    dnld_gaf2 = 'http://downloads.yeastgenome.org/curation/literature/gene_association.sgd.gaf.gz'
    _dnld_gaf(fin_gaf2, dnld_gaf2)

    # Read GAF
    print('READING: {GAF}'.format(GAF=basename(fin_gaf1)))
    ##DVK objanno_sgd1 = GafReader(fin_gaf1)

    print('READING: {GAF}'.format(GAF=basename(fin_gaf2)))
    objanno_sgd2 = GafReader(fin_gaf2)
    ##DVK assert objanno_sgd1 == objanno_sgd2


def _dnld_gaf(fin_gaf, dnld_gaf):
    """Download saccharomyces GAF"""
    if not exists(fin_gaf):
        system('rm -f {GAF}*'.format(GAF=fin_gaf))
        wget.download(dnld_gaf)
        fin_gaz = '{GAF}.gz'.format(GAF=fin_gaf)
        system('gunzip {ZIP}'.format(ZIP=fin_gaz))


if __name__ == '__main__':
    test_i195()
