"""Test Gene Ontology Enrichement Analysis using mutipletest methods in statsmodels."""

import sys
import os
sys.path.insert(0, "..") # Use local version of goatools during test
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.associations import read_associations

__copyright__ = "Copyright (C) 2010-2016, H Tang et al., All rights reserved."

def test_goea_statsmodels(log=sys.stdout):
    """Test GOEA with local multipletest correction methods."""
    goeaobj = get_goeaobj()
    study_ids = [line.rstrip() for line in open("../data/study")]
    prt_if = lambda nt: nt.p_uncorrected < 0.00005
    # These will specify to use the local methods
    #methods = ['bonferroni', 'sidak', 'holm']
    # These will specify to use the statsmodels methods
    methods = ['holm-sidak', 'simes-hochberg', 'hommel', 
               'fdr_bh', 'fdr_by', 'fdr_tsbh', 'fdr_tsbky']
    for method in methods:
        goea_results = goeaobj.run_study(study_ids, methods=[method])
        # Make format_string. Examples:
        # "{NS} {p_uncorrected:5.3e} {p_fdr_bh:5.3e} {name} ({study_count} gene(s))\n"
        fmtstr = "".join(["{NS} {p_uncorrected:5.3e} {",
                  "p_{M}:5.3e".format(M=method), 
                  "} {name} ({study_count} gene(s))\n"])
        log.write("\nSTATSMODELS METHOD: {M}\n".format(M=method))
        goeaobj.prt_txt(log, goea_results, fmtstr, prt_if=prt_if)
        fout_xlsx = "goea_statsmodels_{M}.xlsx".format(M=method)
        fout_tsv = "goea_statsmodels_{M}.tsv".format(M=method)
        goeaobj.prt_txt(log, goea_results, fmtstr, prt_if=prt_if)
        goeaobj.wr_xlsx(fout_xlsx, goea_results)
        goeaobj.wr_tsv(fout_tsv, goea_results)

def get_goeaobj(methods=None):
    """Test GOEA with method, fdr."""
    obo_dag = GODag("go-basic.obo")
    assoc = read_associations("../data/association", no_top=True)
    popul_ids = [line.rstrip() for line in open("../data/population")]
    goeaobj = GOEnrichmentStudy(popul_ids, assoc, obo_dag, methods=methods)
    return goeaobj

if __name__ == '__main__':
    test_goea_statsmodels()

# Copyright (C) 2010-2016, H Tang et al., All rights reserved.
