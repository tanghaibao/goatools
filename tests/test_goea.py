"""Tests GOEA using multiple-test corrections from statsmodels."""
# The tests in this file are intended to be run from the directory in which they reside.

__copyright__ = "Copyright (C) 2016, DV Klopfenstein, H Tang. All rigths reserved."
__author__ = "DV Klopfenstein"

import sys
from operator import attrgetter

sys.path.insert(0, '..') # Use local version of goatools during test
from goatools.obo_parser import GODag
from goatools.associations import read_associations

from PyBiocode.enrichanal.enrichanal_GO import GOEA

def test_fdr_bh(log):
    """Test Gene Ontology Enrichment Analysis using Benjamini/Hochberg multiple testing."""
    # Initialize
    obo_dag = GODag("go-basic.obo")
    assoc = read_associations("../data/association")
    popul_ids = [line.rstrip() for line in open("../data/population")]
    study_ids = [line.rstrip() for line in open("../data/study")]
    # Run enrichment analysis
    goea = GOEA(obo_dag, assoc, log)
    goea.set_population(popul_ids)
    goea.set_params(alpha=0.05, method='fdr_bh')
    results_nt = goea.find_enrichment(study_ids)
    # Print results 3 ways: to screen, to tsv(tab-separated file), to xlsx(Excel spreadsheet)
    field_names = ['study_cnt', 'fdr_bh', 'name']
    sort_by = lambda nt: nt.fdr_bh # Sort by corrected pval, with smallest first
    # 1. Print results to screen using format in prtfmt. For example:
    #
    #     22 1.353e-03 protein phosphorylation
    #      3 1.409e-03 palmitoyl-(protein) hydrolase activity
    #      6 1.554e-03 peptidase activity
    #      4 2.062e-03 proteasome core complex
    #      2 2.520e-03 CDP-alcohol phosphatidyltransferase activity
    #      2 2.520e-03 CDP-diacylglycerol-glycerol-3-phosphate 3-phosphatidyltransferase activity
    #      ...
    prtfmt = "{study_cnt:2} {fdr_bh:5.3e} {name}\n"
    goea.prt_txt(sys.stdout, results_nt, field_names, prtfmt, sort_by=sort_by)

if __name__ == '__main__':
    test_fdr_bh(log=sys.stdout)

# Copyright (C) 2016, DV Klopfenstein, H Tang. All rights reserved.
