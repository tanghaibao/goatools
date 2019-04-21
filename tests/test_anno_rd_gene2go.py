#!/usr/bin/env python
"""Ensure NEW results are equal to OLD results: read_ncbi_gene2go."""

from __future__ import print_function

import os
from goatools.associations import dnld_ncbi_gene_file
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.associations import read_ncbi_gene2go
from goatools.associations import read_ncbi_gene2go_old


REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def test_anno_read():
    """Test reading an NCBI gene2go annotation file."""
    fin_anno = os.path.join(REPO, 'gene2go')
    _dnld_anno(fin_anno)
    #godag = get_godag(os.path.join(REPO, 'go-basic.obo'), loading_bar=None)

    print('\nTEST STORING ONLY ONE SPECIES')
    obj = Gene2GoReader(fin_anno)
    assert len(obj.taxid2asscs) == 1
    obj.prt_summary_anno2ev()

    print('\nTEST STORING ALL SPECIES')
    obj = Gene2GoReader(fin_anno, taxids=True)
    assert len(obj.taxid2asscs) > 1, '**EXPECTED MORE: len(taxid2asscs) == {N}'.format(
        N=len(obj.taxid2asscs))
    obj.prt_summary_anno2ev()

    print('\nTEST GETTING ASSOCIATIONS FOR ONE SPECIES')
    print("\nTEST read_ncbi_gene2go_old: [9606]")
    old_g2go_hsa = read_ncbi_gene2go_old(fin_anno, [9606])
    assert old_g2go_hsa == read_ncbi_gene2go(fin_anno, [9606])
    print("\nTEST read_ncbi_gene2go_old: 9606")
    assert old_g2go_hsa == read_ncbi_gene2go(fin_anno, 9606)
    print("\nTEST read_ncbi_gene2go_old: None")
    assert old_g2go_hsa == read_ncbi_gene2go(fin_anno, None)

    print('\nTEST GETTING REVERSE ASSOCIATIONS: GO2GENES')
    go2geneids = True
    print("\nTEST read_ncbi_gene2go_old: 9606 go2geneids=True")
    old_go2gs_hsa = read_ncbi_gene2go_old(fin_anno, [9606], go2geneids=go2geneids)
    new_go2gs_hsa = read_ncbi_gene2go(fin_anno, 9606, go2geneids=go2geneids)
    print('OLD:', next(iter(old_go2gs_hsa.items())))
    print('NEW:', next(iter(new_go2gs_hsa.items())))
    assert old_go2gs_hsa == new_go2gs_hsa

    print('\nTEST RETURNING ASSOCIATIONS FOR SELECTED EVIDENCE CODES')
    evcodes = set(['ISO', 'IKR'])
    print("\nTEST read_ncbi_gene2go_old: 9606 evcodes=True")
    old_gene2gos_evc = read_ncbi_gene2go_old(fin_anno, [9606], evidence_set=evcodes)
    new_gene2gos_evc = read_ncbi_gene2go(fin_anno, 9606, evidence_set=evcodes)
    print('OLD:', next(iter(old_gene2gos_evc.items())))
    print('NEW:', next(iter(new_gene2gos_evc.items())))
    assert old_gene2gos_evc == new_gene2gos_evc


def _dnld_anno(file_anno):
    """Download the annotation file, if needed."""
    if os.path.exists(file_anno):
        assert os.path.getsize(file_anno) > 1000000, "BAD ANNO({F})".format(F=file_anno)
        return
    dnld_ncbi_gene_file(file_anno, loading_bar=None)
    assert os.path.isfile(file_anno), "MISSING ANNO({F})".format(F=file_anno)
    assert os.path.getsize(file_anno) > 1000000, "BAD ANNO({F})".format(F=file_anno)


if __name__ == '__main__':
    test_anno_read()
