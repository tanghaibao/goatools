#!/usr/bin/env python3
"""Test GoeaResults can be converted to a Networkx Graph"""

import collections as cx
from goatools.rpt.goea_nt_xfrm import MgrNtGOEAs  # get_goea_nts_all
from goatools.test_data.genes_NCBI_10090_ProteinCoding import GENEID2NT as GeneID2nt_mus
from goatools.test_data.nature3102_goea import get_geneid2symbol, get_goeaobj
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.gosubdag.plot.plot import plt_goids
from goatools.gosubdag.plot.plot import plot_gos
from goatools.gosubdag.plot.plot import plot_results
from goatools.gosubdag.plot.go_to_nx import GOs_to_nx

__copyright__ = "Copyright (C) 2016-2017, DV Klopfenstein, H Tang, All rights reserved."

def test_i206():
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
    print('{N} GOEA results'.format(N=len(goea_results_sig)))

    grph = GOs_to_nx(goea_results_sig, sig=0.05)
    #print(grph)
    #### goea_results_nt = MgrNtGOEAs(goea_results_sig).get_goea_nts_all()
    #### assert goea_results_nt
    #### # Test plotting GOEA results
    #### gosubdag = GoSubDag(set(r.GO for r in goea_results_sig), go2obj)
    #### plot_results("test_plot_goids_a_goea_{NS}.png", goea_results_sig,
    ####              id2symbol=geneids2symbol_study, parentcnt=True, childcnt=True)
    #### # Plot b and c
    #### ns2gos = get_ns2gos(goea_results_sig)
    #### for nss, goids in ns2gos.items():
    ####     plt_goids(gosubdag, "test_plot_goids_b_{NS}.png".format(NS=nss), goids)
    ####     plot_gos("test_plot_goids_c_{NS}.png".format(NS=nss), goids, go2obj)

#### def get_ns2gos(goea_res):
####     """Return dict with keys BP, MF, CC and values GO IDs."""
####     ns2gos = cx.defaultdict(set)
####     for res in goea_res:
####         ns2gos[res.NS].add(res.GO)
####     return ns2gos


if __name__ == '__main__':
    test_i206()

# Copyright (C) 2016-2017, DV Klopfenstein, H Tang, All rights reserved.
