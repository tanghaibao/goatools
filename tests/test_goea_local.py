"""Test Gene Ontology Enrichement Analysis."""

import sys
import os
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.associations import read_associations

__copyright__ = "Copyright (C) 2010-2016, H Tang et al., All rights reserved."

def test_unknown_gos():
    """Ensure that a study with only unknown GO Terms will run gracefully."""
    os.system("python {SCR} --alpha=0.05 {STUDY} {POP} {ASSN} --obo={OBO}".format(
        SCR="../scripts/find_enrichment.py",
        OBO="../go-basic.obo",
        STUDY="data/study_unknown",
        POP="../data/population",
        ASSN="../data/association"))

def test_goea_fdr_dflt(log=sys.stdout):
    """Test GOEA with method, fdr. Print original summary"""
    goeaobj = get_goeaobj()
    study_ids = [line.rstrip() for line in open("../data/study")]
    goea_results = goeaobj.run_study(study_ids)
    goeaobj.print_summary(goea_results)

def test_goea_local(log=sys.stdout):
    """Test GOEA with local multipletest correction methods for local."""
    goeaobj = get_goeaobj()
    study_ids = [line.rstrip() for line in open("../data/study")]
    prt_if = lambda nt: nt.p_uncorrected < 0.00005
    for method in ("fdr", "bonferroni", "sidak", "holm"):
        goea_results = goeaobj.run_study(study_ids, methods=[method])
        # Make format_string. Examples:
        # "{NS} {p_uncorrected:5.3e} {p_fdr:5.3e} {name} ({study_count} gene(s))\n"
        # "{NS} {p_uncorrected:5.3e} {p_bonferroni:5.3e} {name} ({study_count} gene(s))\n"
        # "{NS} {p_uncorrected:5.3e} {p_sidak:5.3e} {name} ({study_count} gene(s))\n"
        fmtstr = "".join(["{NS} {p_uncorrected:5.3e} {",
                  "p_{M}:5.3e".format(M=method), 
                  "} {name} ({study_count} gene(s))\n"])
        goeaobj.prt_txt(log, goea_results, fmtstr, prt_if=prt_if)

def test_goea_bonferroni(log=sys.stdout):
    """Test GOEA with method, bonferroni."""
    goeaobj = get_goeaobj(['bonferroni'])
    study_ids = [line.rstrip() for line in open("../data/study")]
    goea_results = goeaobj.run_study(study_ids)
    # Only print if bonferonni value < 0.05
    prt_if = lambda nt: nt.p_bonferroni < 0.05
    # Print to tab-separated table and Excel spreadsheet
    goeaobj.wr_tsv("goea_bonferroni.tsv", goea_results, prt_if=prt_if)
    goeaobj.wr_xlsx("goea_bonferroni.xlsx", goea_results, prt_if=prt_if)
    # Print level in addition to all the regular fields
    # User can control which fields are printed and the order that they appear in the table
    prt_flds = "NS level GO enrichment name ratio_in_study ratio_in_pop p_uncorrected p_bonferroni".split()
    goeaobj.wr_xlsx("goea_bonferroni_lev.xlsx", goea_results, prt_if=prt_if, prt_flds=prt_flds)

def get_goeaobj(methods=None):
    """Test GOEA with method, fdr."""
    obo_dag = GODag("go-basic.obo")
    assoc = read_associations("../data/association", no_top=True)
    popul_ids = [line.rstrip() for line in open("../data/population")]
    goeaobj = GOEnrichmentStudy(popul_ids, assoc, obo_dag, methods=methods)
    return goeaobj

def run_all():
    """Run all local multiple tests."""
    test_unknown_gos()
    test_goea_fdr_dflt()
    test_goea_local()
    test_goea_bonferroni()

if __name__ == '__main__':
    run_all()

# Copyright (C) 2010-2016, H Tang et al., All rights reserved.
