#!/usr/bin/env python
"""Test GoeaResults in plotting package."""

import collections as cx
from goatools.rpt.goea_nt_xfrm import MgrNtGOEAs  # get_goea_nts_all
from goatools.test_data.genes_NCBI_10090_ProteinCoding import GENEID2NT as GeneID2nt_mus
from goatools.test_data.nature3102_goea import get_geneid2symbol, get_goeaobj
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.gosubdag.plot.plot import plt_goids
from goatools.gosubdag.plot.plot import plot_gos
from goatools.gosubdag.plot.plot import plot_results

__copyright__ = "Copyright (C) 2016-2017, DV Klopfenstein, H Tang, All rights reserved."

def test_example():
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
    go2obj = goeaobj.obo_dag
    # Run GOEA on study
    goea_results_all = goeaobj.run_study(geneids_study)
    goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]
    goea_results_nt = MgrNtGOEAs(goea_results_sig).get_goea_nts_all()
    assert goea_results_nt
    ns2gos = get_ns2gos(goea_results_sig)
    # Test plotting GOEA results
    gosubdag = GoSubDag(set(r.GO for r in goea_results_sig), go2obj)
    plot_results("test_plot_goids_a_goea_{NS}.png", goea_results_sig,
                 id2symbol=geneids2symbol_study, parentcnt=True, childcnt=True)
    for nss, goids in ns2gos.items():
        plt_goids(gosubdag, "test_plot_goids_b_{NS}.png".format(NS=nss), goids)
        plot_gos("test_plot_goids_c_{NS}.png".format(NS=nss), goids, go2obj)

def get_ns2gos(goea_res):
    """Return dict with keys BP, MF, CC and values GO IDs."""
    ns2gos = cx.defaultdict(set)
    for res in goea_res:
        ns2gos[res.NS].add(res.GO)
    return ns2gos

def get_str(res, color):
    """Get short string with GOEA info relevant for this test."""
    return "{COL:15} {GO} {PVAL:8.2e} {E}\n".format(
        COL=color, GO=res.GO, PVAL=res.p_fdr_bh, E=res.enrichment)


if __name__ == '__main__':
    test_example()

# Copyright (C) 2016-2017, DV Klopfenstein, H Tang, All rights reserved.
