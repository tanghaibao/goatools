#!/usr/bin/env python
"""Test GoeaResults in plotting package."""

import os
import sys
from goatools.rpt.goea_nt_xfrm import MgrNtGOEAs  # get_goea_nts_all
from goatools.test_data.genes_NCBI_10090_ProteinCoding import GENEID2NT as GeneID2nt_mus
from goatools.test_data.nature3102_goea import get_geneid2symbol, get_goeaobj
from goatools.gosubdag.plot.plot import plot_results
from goatools.gosubdag.plot.goea_results import GoeaResults

__copyright__ = "Copyright (C) 2016-2017, DV Klopfenstein, H Tang, All rights reserved."

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../")


def test_example(prt=sys.stdout):
    """Test GoeaResults in plotting package."""
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
    goea_results_nt = MgrNtGOEAs(goea_results_sig).get_goea_nts_all()
    # Test managing GOEA results
    objres = GoeaResults(goea_results_sig)
    run(objres, prt)
    objnts = GoeaResults(goea_results_nt)
    run(objnts, prt)
    # Plot GOEA results
    fout_img = os.path.join(REPO, "test_plot_objgoearesults_{NS}.png")
    plot_results(fout_img, goea_results_sig, id2symbol=geneids2symbol_study)

def run(obj, prt):
    """Run test for either GOEA objects or namedtuples."""
    go2colors = {}
    obj.prt_summary(prt)
    obj.set_goid2color_pval(go2colors)
    for goid, res in sorted(obj.go2res.items(), key=lambda t: t[1].p_fdr_bh):
        if res.enrichment == 'e':
            color = go2colors.get(goid, -1.0)
            if res.p_fdr_bh <= 0.005:
                assert color == obj.alpha2col[0.005], get_str(res, color)
            elif res.p_fdr_bh <= 0.010:
                assert color == obj.alpha2col[0.010], get_str(res, color)
            elif res.p_fdr_bh <= 0.050:
                assert color == obj.alpha2col[0.050], get_str(res, color)

def get_str(res, color):
    """Get short string with GOEA info relevant for this test."""
    return "{COL:15} {GO} {PVAL:8.2e} {E}\n".format(
        COL=color, GO=res.GO, PVAL=res.p_fdr_bh, E=res.enrichment)


if __name__ == '__main__':
    test_example()

# Copyright (C) 2016-2017, DV Klopfenstein, H Tang, All rights reserved.
