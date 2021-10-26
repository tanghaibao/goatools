#!/usr/bin/env python3
"""Test fix for issue#202 Enrichment analyses w/human phenotype ontologies"""

from __future__ import print_function

from os.path import join
from tests.utils import REPO

from goatools.obo_parser import GODag
from goatools.anno.idtogos_reader import IdToGosReader
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.utils import read_geneset


def test_i202():
    """Test to re-produce issue#202: Passes currently."""
    fin_study = 'notebooks/data/hpo/genes.list'
    fin_pop   = 'notebooks/data/hpo/gobackground.list'
    fin_obo   = 'notebooks/data/hpo/hp.obo'
    fin_anno  = 'notebooks/data/hpo/hpo.annotation.tab'

    study_ids = read_geneset(join(REPO, fin_study))
    population_ids = read_geneset(join(REPO, fin_pop))
    godag = GODag(join(REPO, fin_obo))
    annoobj = IdToGosReader(join(REPO, fin_anno), godag=godag)

    ## obj = _Run(9606, 'gene2go', 'go-basic.obo')

    # Result is the same whether fisher_scipy_stats of fisher
    pvalcalc = 'fisher_scipy_stats'
    goeaobj = GOEnrichmentStudy(
        population_ids,
        annoobj.get_id2gos(),
        godag,
        methods=['bonferroni', 'fdr_bh'],
        pvalcalc=pvalcalc)
    # Run GOEA Gene Ontology Enrichment Analysis
    results_goeas = goeaobj.run_study_nts(study_ids)
    print('NS GO         p stu_ratio pop_ratio    p-uncorr bonferro fdr_bh   stu  ')
    for ntd in results_goeas:
        if ntd.study_count == 0:
            ## doprt = False
            if ntd.p_bonferroni < 0.05:
                assert ntd.enrichment == 'p'
                ## doprt = True
            if ntd.p_fdr_bh < 0.05:
                assert ntd.enrichment == 'p'
                ## doprt = True
            ## if doprt:
            ##     print(obj.str_nt(ntd))
    print(next(iter(results_goeas))._fields)


if __name__ == '__main__':
    test_i202()
