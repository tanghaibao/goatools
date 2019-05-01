#!/usr/bin/env python
"""Test GoeaPrintFunctions::print_results."""

import os
# from goatools.rpt.goea_nt_xfrm import MgrNtGOEAs  # get_goea_nts_all
from goatools.test_data.genes_NCBI_10090_ProteinCoding import GENEID2NT as GeneID2nt_mus
from goatools.test_data.nature3102_goea import get_geneid2symbol, get_goeaobj
from goatools.godag.prtfncs import GoeaPrintFunctions

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang, All rights reserved."

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../")


def test_example():
    """Test GoeaPrintFunctions::print_results."""
    # --------------------------------------------------------------------
    # --------------------------------------------------------------------
    # Gene Ontology Enrichment Analysis (GOEA)
    # --------------------------------------------------------------------
    # --------------------------------------------------------------------
    taxid = 10090 # Mouse study
    # Load ontologies, associations, and population ids
    geneids_pop = GeneID2nt_mus.keys()
    geneids2symbol_study = get_geneid2symbol("nbt.3102-S4_GeneIDs.xlsx")
    geneids_study = geneids2symbol_study.keys()
    goeaobj = get_goeaobj("fdr_bh", geneids_pop, taxid)
    # Run GOEA on study
    goea_results_all = goeaobj.run_study(geneids_study)
    goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]
    #goea_results_nt = MgrNtGOEAs(goea_results_sig).get_goea_nts_all()
    objprtres = GoeaPrintFunctions()
    objprtres.print_results(goea_results_sig)
    objprtres.print_date()


if __name__ == '__main__':
    test_example()

# Copyright (C) 2016-2019, DV Klopfenstein, H Tang, All rights reserved.
