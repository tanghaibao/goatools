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

def test_goea():
    """Test GOEA with method, fdr."""
    obo_dag = GODag("go-basic.obo")
    assoc = read_associations("../data/association", no_top=True)
    popul_ids = [line.rstrip() for line in open("../data/population")]
    study_ids = [line.rstrip() for line in open("../data/study")]
    goeaobj = GOEnrichmentStudy(popul_ids, assoc, obo_dag, methods=['fdr'])
    goea_results = goeaobj.run_study(study_ids)
    goeaobj.print_summary(goea_results)

def run_all():
    """Run all tests."""
    test_unknown_gos()
    test_goea()

if __name__ == '__main__':
  run_all()

