"""Tests GOEA using multiple-test corrections from statsmodels."""
# The tests in this file are intended to be run from the directory in which they reside.

__copyright__ = "Copyright (C) 2016, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import sys
from operator import attrgetter

sys.path.insert(0, '..') # Use local version of goatools during test
from goatools.obo_parser import GODag
from goatools.associations import read_associations

from PyBiocode.enrichanal.enrichanal_GO import GOEA

def test_fdr_bh(log):
    """Test Gene Ontology Enrichment Analysis using Benjamini/Hochberg multiple testing."""
    # ---------------------------------------------------------------------
    # Run Gene Ontology Analysis (GOEA)
    #
    # 1. Initialize
    obo_dag = GODag("go-basic.obo")
    assoc = read_associations("../data/association", no_top=True)
    popul_ids = [line.rstrip() for line in open("../data/population")]
    study_ids = [line.rstrip() for line in open("../data/study")]
    # 2. Run enrichment analysis
    goea = GOEA(obo_dag, assoc, log)
    goea.set_population(popul_ids)
    goea.set_params(alpha=0.05, method='fdr_bh')
    results_nt = goea.find_enrichment(study_ids)

    # ---------------------------------------------------------------------
    # Print results 3 ways: to screen, to tsv(tab-separated file), to xlsx(Excel spreadsheet)
    #
    field_names = ['NS', 'study_cnt', 'fdr_bh', 'level', 'depth', 'GO', 'name']
    # User customizable sort. Sort by: First: BP, MF, CC; Second: corrected pval, with smallest first.
    sort_by = lambda nt: [nt.NS, nt.fdr_bh]
    # 1. Print results to screen using format in prtfmt. For example:
    #
    #      BP 22 3.073e-03 L06 D07 GO:0006468 protein phosphorylation
    #      BP  9 1.023e-02 L07 D08 GO:0006511 ubiquitin-dependent protein catabolic process
    #      BP  2 1.023e-02 L05 D09 GO:0019877 diaminopimelate biosynthetic process
    #      BP  2 1.223e-02 L04 D08 GO:0006301 postreplication repair
    #      BP  2 1.223e-02 L05 D09 GO:0030418 nicotianamine biosynthetic process
    #      BP  2 1.492e-02 L04 D06 GO:0006909 phagocytosis
    #      BP  2 1.492e-02 L03 D03 GO:0051322 anaphase
    #      ...
    # Print format field names are the same names as in the 'field_names" variable.
    prtfmt = "{NS} {study_cnt:2} {fdr_bh:5.3e} L{level:02} D{depth:02} {GO} {name}\n"
    goea.prt_txt(sys.stdout, results_nt, field_names, prtfmt, sort_by=sort_by)

if __name__ == '__main__':
    test_fdr_bh(log=sys.stdout)

# Copyright (C) 2016, DV Klopfenstein, H Tang. All rights reserved.
