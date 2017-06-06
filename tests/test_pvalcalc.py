#!/usr/bin/env python
"""Test that two different but equivalent fishers functions give the similar results."""

import sys
import os
import xlrd
import collections as cx

from goatools.go_enrichment import GOEnrichmentStudy
from goatools.pvalcalc import FisherFactory
from goatools.base import download_go_basic_obo
from goatools.obo_parser import GODag
from goatools.test_data.genes_NCBI_10090_ProteinCoding import GeneID2nt as GeneID2nt_mus
from goatools.associations import get_assoc_ncbi_taxids
from itertools import combinations

def test_pvalcalc(prt=None):
    pvalfnc_names = ['fisher', 'fisher_scipy_stats']
    fisher2pvals = _get_pvals(pvalfnc_names)
    _chk_pvals(fisher2pvals, prt)

def _chk_pvals(fisher2pvals, prt):
    for f1, f2 in combinations(fisher2pvals.keys(), 2):
        ctr = cx.Counter()
        pvals1 = cx.OrderedDict(sorted([(r.GO, r.p_uncorrected) for r in fisher2pvals[f1]]))
        pvals2 = cx.OrderedDict(sorted([(r.GO, r.p_uncorrected) for r in fisher2pvals[f2]]))
        assert len(pvals1) == len(pvals2)
        for go_id, pval1 in pvals1.items():
            pval2 = pvals2[go_id]
            ctr[pval1 == pval2] += 1
            if prt is not None:
                # Are values from 'fisher' and scipy stats 'fisher_exact' equivalent?
                if abs(pval1 - pval2) > 0.00001: 
                    prt.write("{GO} {N1}({P1:4.2f}) {N2}({P2:4.2f}) {N1}({p1}) {N2}({p2})\n".format(
                        GO=go_id, N1=f1, N2=f2, P1=pval1, P2=pval2, p1=pval1, p2=pval2))
        if prt is not None:
            for val, cnt in ctr.most_common():
                prt.write("{N:>5,} {RES:5}: {N1} == {N2}\n".format(N=cnt, RES=val, N1=f1, N2=f2))

def _get_pvals(pvalfnc_names, prt=sys.stdout):
    fisher2pvals = {}
    taxid = 10090 # Mouse study
    obo_dag = GODag(download_go_basic_obo(prt=prt))
    geneids_pop = GeneID2nt_mus.keys()
    assoc_geneid2gos = get_assoc_ncbi_taxids([taxid])
    geneids_study = _get_geneid2symbol("nbt.3102-S4_GeneIDs.xlsx", prt)
    for fisher in pvalfnc_names:
        goeaobj = GOEnrichmentStudy(
            geneids_pop,
            assoc_geneid2gos,
            obo_dag,
            propagate_counts = False,
            alpha = 0.05,
            methods = None,
            pvalcalc = fisher)
        fisher2pvals[fisher] = goeaobj._get_pval_uncorr(geneids_study, prt)
    return fisher2pvals

def _get_geneid2symbol(fin_xlsx, log):
    """Read xlsx file return dictionary with Entrez GeneID keys to Symbol data."""
    gene2symbol = {}
    data_dir = os.path.dirname(os.path.abspath(__file__)) + "/data/nbt_3102"
    tbl_genes = "{DIR}/{FIN}".format(DIR=data_dir, FIN=fin_xlsx)
    book = xlrd.open_workbook(tbl_genes)
    pg = book.sheet_by_index(0)
    for r in range(pg.nrows):
        symbol, geneid, pval = [pg.cell_value(r, c) for c in range(pg.ncols)]
        if geneid:
            gene2symbol[int(geneid)] = symbol
    return gene2symbol

if __name__ == '__main__':
    test_pvalcalc(sys.stdout)
