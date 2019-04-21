import sys
import os
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.associations import read_associations

"""Test Gene Ontology Enrichment Analysis using mutipletest methods in statsmodels."""

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/data/"

__copyright__ = "Copyright (C) 2010-2018, H Tang et al., All rights reserved."


def test_goea_statsmodels(log=sys.stdout):
    """Test GOEA with local multipletest correction methods for statsmodels."""
    goeaobj = get_goeaobj()
    study_ids = [line.rstrip() for line in open(ROOT + "small_study")]
    prt_if = lambda nt: nt.p_uncorrected < 0.0005
    ## These will specify to use the statsmodels methods
    methods_sm0 = ['holm-sidak', 'simes-hochberg', 'hommel',
                   'fdr_bh', 'fdr_by', 'fdr_tsbh', 'fdr_tsbky']
                   # 'fdr_gbs' generates a zerodivision warning
    # Prepend "sm_" or "statsmodels_" to a method to use that version
    methods_sm1 = ['sm_bonferroni', 'sm_sidak', 'sm_holm']
    methods = methods_sm0 + methods_sm1

    for method in methods:
        log.write("\nSTATSMODELS METHOD: {M}\n".format(M=method))
        goea_results = goeaobj.run_study(study_ids, methods=[method])
        # Make format_string. Examples:
        # "{NS} {p_uncorrected:5.3e} {p_fdr_bh:5.3e} {name} ({study_count} gene(s))\n"
        fmtstr = "".join(["{NS} {p_uncorrected:5.3e} {",
                  "p_{M}:5.3e".format(M=method),
                  "} {name} ({study_count} gene(s))\n"])
        fout_xlsx = "goea_statsmodels_{M}.xlsx".format(M=method)
        fout_tsv = "goea_statsmodels_{M}.tsv".format(M=method)
        goeaobj.prt_txt(log, goea_results, fmtstr, prt_if=prt_if)
        goeaobj.wr_xlsx(fout_xlsx, goea_results)
        goeaobj.wr_tsv(fout_tsv, goea_results)


def get_goeaobj(methods=None):
    """Test GOEA with method, fdr."""
    obo_dag = GODag(ROOT + "goslim_generic.obo")
    fin_assc = ROOT + "slim_association"
    assoc = read_associations(fin_assc, 'id2gos', no_top=True)
    popul_ids = [line.rstrip() for line in open(ROOT + "small_population")]
    goeaobj = GOEnrichmentStudy(popul_ids, assoc, obo_dag, methods=methods)
    return goeaobj


if __name__ == '__main__':
    test_goea_statsmodels()

# Copyright (C) 2010-2018, H Tang et al., All rights reserved.
