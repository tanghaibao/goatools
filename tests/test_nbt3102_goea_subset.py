#!/usr/bin/env python
"""Run a Gene Ontology Enrichment Analysis (GOEA), plots, etc.

    Nature 2014_0126;
			Computational analysis of cell-to-cell heterogeneity
      in single-cell RNA-sequencing data reveals hidden
      subpopulations of cells
    http://www.nature.com/nbt/journal/v33/n2/full/nbt.3102.html#methods
"""

from goatools.test_data.genes_NCBI_10090_ProteinCoding import GENEID2NT as GeneID2nt_mus
from goatools.test_data.nature3102_goea import get_geneid2symbol, get_goeaobj

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

def test_example():
    """Run Gene Ontology Enrichment Analysis (GOEA) on Nature data."""
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
    goea_results_sub = [r for r in goea_results_sig if r.study_count > r.study_n/10]
    # Print GOEA results to files: With study genes printed as geneids or symbols
    goeaobj.wr_xlsx("nbt3102_symbols.xlsx", goea_results_sub, itemid2name=geneids2symbol_study)
    goeaobj.wr_xlsx("nbt3102_geneids.xlsx", goea_results_sub)

if __name__ == '__main__':
    test_example()

# Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved.
