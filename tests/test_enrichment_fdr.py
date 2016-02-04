"""Test Gene Ontology Enrichement Analysis."""

import sys
import os
sys.path.insert(0, "..") # Use local version of goatools during test
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.associations import read_associations

def test_unknown_gos():
    """Ensure that a study with only unknown GO Terms will run gracefully."""
    os.system("python {SCR} --alpha=0.05 {STUDY} {POP} {ASSN} --fdr --obo={OBO}".format(
        SCR="../scripts/find_enrichment.py",
        OBO="../go-basic.obo",
        STUDY="data/study_unknown",
        POP="../data/population",
        ASSN="../data/association"))

def test_goea_fdr():
    """Test GOEA with method, fdr."""
    run_goea(['fdr'])

def test_goea_bonferroni(log=sys.stdout):
    """Test GOEA with method, bonferroni."""
    goea_results, goeaobj = run_goea(['bonferroni'])
    # Test printing GOEA results in tables: text, tab-separated, Excel spreadsheet
    goeaobj.prt_tsv(log, goea_results) # Print to screen
    # Only print if bonferonni value < 0.05
    prt_if = lambda nt: nt.p_bonferroni < 0.05
    # Print to tab-separated table and Excel spreadsheet
    goeaobj.wr_tsv("goea_bonferroni.tsv", goea_results, prt_if=prt_if)
    goeaobj.wr_xlsx("goea_bonferroni.xlsx", goea_results, prt_if=prt_if)
    # Print level in addition to all the regular fields
    # User can control which fields are printed and the order that they appear in the table
    prt_flds = "NS level GO enrichment name ratio_in_study ratio_in_pop p_uncorrected p_bonferroni".split()
    goeaobj.wr_xlsx("goea_bonferroni_lev.xlsx", goea_results, prt_if=prt_if, prt_flds=prt_flds)

def run_goea(methods):
    """Test GOEA with method, fdr."""
    obo_dag = GODag("go-basic.obo")
    assoc = read_associations("../data/association", no_top=True)
    popul_ids = [line.rstrip() for line in open("../data/population")]
    study_ids = [line.rstrip() for line in open("../data/study")]
    goeaobj = GOEnrichmentStudy(popul_ids, assoc, obo_dag, methods=methods)
    goea_results = goeaobj.run_study(study_ids)
    goeaobj.print_summary(goea_results)
    return goea_results, goeaobj

def run_all():
    """Run all tests."""
    test_unknown_gos()
    test_goea_fdr()
    test_goea_bonferroni()

if __name__ == '__main__':
    #run_all()
    test_goea_bonferroni()

