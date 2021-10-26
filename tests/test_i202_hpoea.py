#!/usr/bin/env python3
"""Test to re-produce issue#202: Passes currently."""

from __future__ import print_function

import os
from os.path import join
from tests.utils import REPO

from goatools.base import get_godag
from goatools.associations import dnld_ncbi_gene_file
from goatools.associations import read_ncbi_gene2go
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.utils import read_geneset


def test_i202():
    """Test to re-produce issue#202: Passes currently."""
    fin_study = 'tests/data/i202_HPO_obo/genes.list'
    fin_pop = 'tests/data/i202_HPO_obo/gobackground.list'
    fin_obo = 'tests/data/i202_HPO_obo/hp.obo'

    study_ids = read_geneset(join(REPO, fin_study))
    population_ids = read_geneset(join(REPO, fin_pop))
    dag = GODag(
    return

    obj = _Run(9606, 'gene2go', 'go-basic.obo')

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
        _fin = join(REPO, fin_gene2go)
        dnld_ncbi_gene_file(_fin, loading_bar=None)
        self.gene2go = read_ncbi_gene2go(_fin, [taxid])

        _fin_obo = join(REPO, fin_gobasic)
        self.godag = get_godag(_fin_obo, loading_bar=None)

    def str_nt(self, ntd):
        return self.patrec.format(
            NS=ntd.NS, GO=ntd.GO,
            RS=ntd.ratio_in_study, RP=str(ntd.ratio_in_pop),
            e=ntd.enrichment,
            PVAL=ntd.p_uncorrected, BONF=ntd.p_bonferroni, BH=ntd.p_fdr_bh,
            STU=ntd.study_items)


if __name__ == '__main__':
    test_i202()
