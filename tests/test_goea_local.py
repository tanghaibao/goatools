"""Test Gene Ontology Enrichement Analysis."""

import sys
import os
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.associations import read_associations
from goatools.godag.prtfncs import GoeaPrintFunctions
from goatools.base import get_godag

__copyright__ = "Copyright (C) 2010-2019, H Tang et al., All rights reserved."

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def test_unknown_gos():
    """Ensure that a study with only unknown GO Terms will run gracefully."""
    #pylint: disable=bad-whitespace
    code = os.system("python {SCR} --alpha=0.05 {STUDY} {POP} {ASSN} --obo={OBO}".format(
        SCR  ="{REPO}/scripts/find_enrichment.py".format(REPO=REPO),
        OBO  ="{REPO}/go-basic.obo".format(REPO=REPO),
        STUDY="{REPO}/tests/data/study_unknown".format(REPO=REPO),
        POP  ="{REPO}/tests/data/small_population".format(REPO=REPO),
        ASSN ="{REPO}/tests/data/small_association".format(REPO=REPO)))
    assert code != 0, "**FAILED: Simple find_enrichment test"

def test_goea_fdr_dflt():
    """Test GOEA with method, fdr. Print original summary"""
    goeaobj = get_goeaobj()
    study_fin = "{REPO}/tests/data/small_study".format(REPO=REPO)
    study_ids = [line.rstrip() for line in open(study_fin)]
    goea_results = goeaobj.run_study(study_ids)
    objprtres = GoeaPrintFunctions()
    objprtres.print_results(goea_results)
    objprtres.print_date()

def test_goea_local(log=sys.stdout):
    """Test GOEA with local multipletest correction methods for local."""
    goeaobj = get_goeaobj()
    study_fin = "{REPO}/tests/data/small_study".format(REPO=REPO)
    study_ids = [line.rstrip() for line in open(study_fin)]
    # prt_if = lambda nt: nt.p_uncorrected < 0.00005
    prt_if = None
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

def test_goea_bonferroni():
    """Test GOEA with method, bonferroni."""
    goeaobj = get_goeaobj(['bonferroni'])
    study_fin = "{REPO}/tests/data/small_study".format(REPO=REPO)
    study_ids = [line.rstrip() for line in open(study_fin)]

    fout_xlsx = "{REPO}/goea_bonferroni_usrflds.xlsx".format(REPO=REPO)
    goea_results = goeaobj.run_study(study_ids)
    prt_flds = ["GO", "NS", "enrichment", "name"]
    # Counts, ratios
    prt_flds.extend(["ratio_in_study", "ratio_in_pop"])
    # These fields have the same info as: ratio_in_study ratio_in_pop
    prt_flds.extend(["study_count", "study_n", "pop_count", "pop_n"])
    prt_flds.extend(["p_uncorrected", "depth", "p_bonferroni", "study_items"])
    goeaobj.wr_xlsx(fout_xlsx, goea_results, prt_flds=prt_flds)

    # Only print if bonferonni value < 0.05
    # prt_if = lambda nt: nt.p_bonferroni < 0.05
    prt_if = None
    # Print to tab-separated table and Excel spreadsheet
    goeaobj.wr_tsv("{REPO}/goea_bonferroni.tsv".format(REPO=REPO), goea_results, prt_if=prt_if)
    # Print level in addition to all the regular fields
    # User can control which fields are printed and the order that they appear in the table
    prt_flds = "NS level GO name ratio_in_study ratio_in_pop p_uncorrected p_bonferroni".split()
    fout_xlsx = "{REPO}/goea_bonferroni_lev.xlsx".format(REPO=REPO)
    goeaobj.wr_xlsx(fout_xlsx, goea_results, prt_if=prt_if, prt_flds=prt_flds)

def get_goeaobj(methods=None):
    """Test GOEA with method, fdr."""
    obo_fin = os.path.join(REPO, "go-basic.obo")
    obo_dag = get_godag(obo_fin, loading_bar=None)
    fin_assc = "{REPO}/tests/data/small_association".format(REPO=REPO)
    assoc = read_associations(fin_assc, 'id2gos', no_top=True)
    popul_fin = "{REPO}/tests/data/small_population".format(REPO=REPO)
    popul_ids = [line.rstrip() for line in open(popul_fin)]
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

# Copyright (C) 2010-2019, H Tang et al., All rights reserved.
