"""Run GOATOOLS Gene Ontology Analysis on Nature 3102 data. Return results."""

__copyright__ = "Copyright (C) 2016, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

from goatools.test_data.genes_NCBI_10090_ProteinCoding import GeneID2nt as GeneID2nt_mus
from test_nbt3102 import get_geneid2symbol, get_goeaobj

def get_goea_results(keep_if=None):
    """Demonstrate printing a subset of all available fields using two methods."""
    if keep_if is None:
        keep_if = lambda nt: getattr(nt, "p_fdr_bh") < 0.05 # keep if results are significant
    # 1. Gene Ontology Enrichment Analysis
    #    1a. Initialize: Load ontologies, associations, and population gene IDs
    taxid = 10090 # Mouse study
    geneids_pop = GeneID2nt_mus.keys() # Mouse protein-coding genes
    goeaobj = get_goeaobj("fdr_bh", geneids_pop, taxid)
    #    1b. Run GOEA
    geneids_study = get_geneid2symbol("nbt.3102-S4_GeneIDs.xlsx")
    return {
        'goea_results' : goeaobj.run_study(geneids_study, keep_if=keep_if),
        'goeaobj' : goeaobj,
        'geneids_study' : geneids_study,
        'geneids_pop' : geneids_pop}


# Copyright (C) 2016, DV Klopfenstein, H Tang. All rights reserved.
