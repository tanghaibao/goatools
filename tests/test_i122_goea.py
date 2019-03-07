#!/usr/bin/env python3
"""Test to re-produce issue#122: Passes currently."""

from __future__ import print_function

import os
import xlrd
from goatools.base import get_godag
from goatools.associations import dnld_ncbi_gene_file
from goatools.associations import read_ncbi_gene2go
from goatools.go_enrichment import GOEnrichmentStudy

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def test_i122():
    """Test to re-produce issue#122: Passes currently."""
    obj = _Run(9606, 'gene2go', 'go-basic.obo')
    study_ids, population_ids = obj.get_genes_study_n_bg()

    # Result is the same whether fisher_scipy_stats of fisher
    pvalcalc = 'fisher_scipy_stats'
    goeaobj = GOEnrichmentStudy(population_ids, obj.gene2go, obj.godag, methods=['bonferroni', 'fdr_bh'], pvalcalc=pvalcalc)
    # Run GOEA Gene Ontology Enrichment Analysis
    results_goeas = goeaobj.run_study_nts(study_ids)
    print('NS GO         p stu_ratio pop_ratio    p-uncorr bonferro fdr_bh   stu  ')
    for ntd in results_goeas:
        if ntd.study_count == 0:
            doprt = False
            if ntd.p_bonferroni < 0.05:
                assert ntd.enrichment == 'p'
                doprt = True
            if ntd.p_fdr_bh < 0.05:
                assert ntd.enrichment == 'p'
                doprt = True
            if doprt:
                print(obj.str_nt(ntd))
    # print(next(iter(results_goeas))._fields)

class _Run():
    """Run test."""

    patrec = '{NS} {GO} {e} {RS} {RP:>12} {PVAL:8.2e} {BONF:8.2e} {BH:8.2e} {STU}'

    def __init__(self, taxid, fin_gene2go, fin_gobasic):
        _fin = os.path.join(REPO, fin_gene2go)
        dnld_ncbi_gene_file(_fin, loading_bar=None)
        self.gene2go = read_ncbi_gene2go(_fin, [taxid])

        _fin_obo = os.path.join(REPO, fin_gobasic)
        self.godag = get_godag(_fin_obo, loading_bar=None)

    def str_nt(self, ntd):
        return self.patrec.format(
            NS=ntd.NS, GO=ntd.GO,
            RS=ntd.ratio_in_study, RP=str(ntd.ratio_in_pop),
            e=ntd.enrichment,
            PVAL=ntd.p_uncorrected, BONF=ntd.p_bonferroni, BH=ntd.p_fdr_bh,
            STU=ntd.study_items)

    def get_genes_study_n_bg(self):
        """Get the study and background genes."""
        genes = {'stu':set(), 'pop':set()}
        fin_xlsx = 'data/i122/study_and_bg_genes.xlsx'
        book = xlrd.open_workbook(os.path.join(fin_xlsx))
        sheet = book.sheet_by_index(0)
        for rownum in range(sheet.nrows):
            gene_stu = self._get_gene(rownum, 0, sheet)
            gene_pop = self._get_gene(rownum, 1, sheet)
            if gene_stu:
                genes['stu'].add(gene_stu)
            if gene_pop:
                genes['pop'].add(gene_pop)
        print('{N} Study genes. {P} Background genes READ: {XLSX}'.format(
            N=len(genes['stu']), P=len(genes['pop']), XLSX=fin_xlsx))
        return genes['stu'], genes['pop']
   
    @staticmethod 
    def _get_gene(row, col, sheet):
        """Return the gene ID, if a gene is found."""
        gene_id = sheet.cell_value(row, col)
        if isinstance(gene_id, float):
            return int(gene_id)


if __name__ == '__main__':
    test_i122()
