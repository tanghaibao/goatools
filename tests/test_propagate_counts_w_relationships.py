#!/usr/bin/env python
"""Test propagate_counts up relationships as well as parent-child links."""

import sys
import os
# from itertools import combinations
# import collections as cx

from goatools.go_enrichment import GOEnrichmentStudy
from goatools.base import get_godag
from goatools.test_data.genes_NCBI_10090_ProteinCoding import GENEID2NT as GeneID2nt_mus
from goatools.test_data.nature3102_goea import get_geneid2symbol
from goatools.associations import get_assoc_ncbi_taxids

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def test_pc_w_rels(prt=sys.stdout):
    """Test P-value calculations."""
    file_obo = os.path.join(REPO, "go-basic.obo")
    godag_r0 = get_godag(file_obo, prt, loading_bar=None)
    godag_r1 = get_godag(file_obo, prt, loading_bar=None, optional_attrs=['relationship'])

    # pylint: disable=line-too-long
    # Check that relationships are used in propgate counts
    results_r0 = _get_results(godag_r1, propagate_counts=True, relationships=False, prt=prt)
    results_r1 = _get_results(godag_r1, propagate_counts=True, relationships=True, prt=prt)
    results_rr = _get_results(godag_r1, propagate_counts=True, relationships={'regulates', 'negatively_regulates', 'positively_regulates'}, prt=prt)
    results_rp = _get_results(godag_r1, propagate_counts=True, relationships={'part_of'}, prt=prt)
    prt.write('{N} results with r0\n'.format(N=len(results_r0)))
    prt.write('{N} results with r1\n'.format(N=len(results_r1)))
    prt.write('{N} results with all regulates\n'.format(N=len(results_rr)))
    prt.write('{N} results with part_of\n'.format(N=len(results_rp)))
    assert len(results_r1) > len(results_rr)
    assert len(results_rr) > len(results_rp)
    assert len(results_rp) > len(results_r0)

    # TBD: Add warning message that relationships are ignored
    # Check that relationships are ignored in propagate counts if they were not loaded
    _get_results(godag_r0, propagate_counts=True, relationships=False, prt=prt)
    ## results_r0b = _get_results(godag_r0, propagate_counts=True, relationships=True, prt=prt)
    ## results_r0c = _get_results(godag_r0, propagate_counts=True, relationships={'regulates', 'negatively_regulates', 'positively_regulates'}, prt=prt)
    ## results_r0d = _get_results(godag_r0, propagate_counts=True, relationships={'part_of'}, prt=prt)
    ## assert len(results_r0a) == len(results_r0b)
    ## assert len(results_r0b) == len(results_r0c)
    ## assert len(results_r0c) == len(results_r0d)

def _get_results(godag, propagate_counts, relationships, prt=sys.stdout):
    """Run a GOEA. Return results"""
    taxid = 10090 # Mouse study
    geneids_pop = set(GeneID2nt_mus.keys())
    assoc_geneid2gos = get_assoc_ncbi_taxids([taxid], loading_bar=None)
    geneids_study = get_geneid2symbol("nbt.3102-S4_GeneIDs.xlsx")
    goeaobj = GOEnrichmentStudy(
        geneids_pop,
        assoc_geneid2gos,
        godag,
        propagate_counts=propagate_counts,
        relationships=relationships,
        alpha=0.05,
        methods=['fdr_bh'])
    return goeaobj.run_study(geneids_study, prt=prt)


if __name__ == '__main__':
    test_pc_w_rels()
