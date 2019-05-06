#!/usr/bin/env python
"""Ensure NEW results are equal to OLD results: read_ncbi_gene2go."""

from __future__ import print_function

import os
import sys
from collections import defaultdict
from goatools.associations import dnld_ncbi_gene_file
from goatools.anno.genetogo_reader import Gene2GoReader


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
    ## new_g2go_hsa = read_ncbi_gene2go(fin_anno, [9606])
    new_g2go_hsa = obj.get_id2gos_nss(taxids=[9606])
    assert old_g2go_hsa == new_g2go_hsa, \
      'OLD({O}) != NEW({N})'.format(O=len(old_g2go_hsa), N=len(new_g2go_hsa))
    print("\nTEST read_ncbi_gene2go_old: 9606")
    ## assert old_g2go_hsa == read_ncbi_gene2go(fin_anno, 9606)
    assert old_g2go_hsa == obj.get_id2gos_nss(taxid=9606)

    print('\nTEST GETTING REVERSE ASSOCIATIONS: GO2GENES')
    go2geneids = True
    print("\nTEST read_ncbi_gene2go_old: 9606 go2geneids=True")
    old_go2gs_hsa = read_ncbi_gene2go_old(fin_anno, [9606], go2geneids=go2geneids)
    ## new_go2gs_hsa = read_ncbi_gene2go(fin_anno, 9606, go2geneids=go2geneids)
    new_go2gs_hsa = obj.get_id2gos_nss(taxid=9606, go2geneids=go2geneids)
    print('OLD:', next(iter(old_go2gs_hsa.items())))
    print('NEW:', next(iter(new_go2gs_hsa.items())))
    assert old_go2gs_hsa == new_go2gs_hsa, \
       'OLD({O}) != NEW({N})'.format(O=len(old_go2gs_hsa), N=len(new_go2gs_hsa))

    print('\nTEST RETURNING ASSOCIATIONS FOR SELECTED EVIDENCE CODES')
    evcodes = set(['ISO', 'IKR'])
    print("\nTEST read_ncbi_gene2go_old: 9606 evcodes=True")
    old_gene2gos_evc = read_ncbi_gene2go_old(fin_anno, taxids=[9606], ev_include=evcodes)
    ## new_gene2gos_evc = read_ncbi_gene2go(fin_anno, 9606, ev_include=evcodes)
    new_gene2gos_evc = obj.get_id2gos_nss(taxid=9606, ev_include=evcodes)
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

# Fomerly in goatools/associations.py file
def read_ncbi_gene2go_old(fin_gene2go, taxids=None, **kws):
    """Read NCBI's gene2go. Return gene2go data for user-specified taxids."""
    # kws: taxid2asscs ev_include
    # Simple associations
    id2gos = defaultdict(set)
    # Optional detailed associations split by taxid and having both ID2GOs & GO2IDs
    # e.g., taxid2asscs = defaultdict(lambda: defaultdict(lambda: defaultdict(set))
    taxid2asscs = kws.get('taxid2asscs', None)
    evs = kws.get('ev_include', None)
    # By default, return id2gos. User can cause go2geneids to be returned by:
    #   >>> read_ncbi_gene2go(..., go2geneids=True
    b_geneid2gos = not kws.get('go2geneids', False)
    if taxids is None: # Default taxid is Human
        taxids = [9606]
    with open(fin_gene2go) as ifstrm:
        # pylint: disable=too-many-nested-blocks
        for line in ifstrm:
            if line[0] != '#': # Line contains data. Not a comment
                line = line.rstrip() # chomp
                flds = line.split('\t')
                if len(flds) >= 5:
                    taxid_curr, geneid, go_id, evidence, qualifier = flds[:5]
                    taxid_curr = int(taxid_curr)
                    # NOT: Used when gene is expected to have function F, but does NOT.
                    # ND : GO function not seen after exhaustive annotation attempts to the gene.
                    ## if taxid_curr in taxids and qualifier != 'NOT' and evidence != 'ND':
                    if taxid_curr in taxids and 'NOT' not in qualifier and evidence != 'ND':
                        # Optionally specify a subset of GOs based on their evidence.
                        if evs is None or evidence in evs:
                            geneid = int(geneid)
                            if b_geneid2gos:
                                id2gos[geneid].add(go_id)
                            else:
                                id2gos[go_id].add(geneid)
                            if taxid2asscs is not None:
                                taxid2asscs[taxid_curr]['ID2GOs'][geneid].add(go_id)
                                taxid2asscs[taxid_curr]['GO2IDs'][go_id].add(geneid)
        sys.stdout.write("  {N:,} items READ: {ASSC}\n".format(N=len(id2gos), ASSC=fin_gene2go))
    return id2gos # return simple associations


if __name__ == '__main__':
    test_anno_read()
