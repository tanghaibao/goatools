"""Run GOATOOLS Gene Ontology Analysis on Nature 3102 data. Return results."""

__copyright__ = "Copyright (C) 2016-2017, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import os
import xlrd
from goatools.test_data.genes_NCBI_10090_ProteinCoding import GeneID2nt as GeneID2nt_mus
from goatools.base import get_godag
from goatools.associations import get_assoc_ncbi_taxids
from goatools.go_enrichment import GOEnrichmentStudy

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
        'geneids_pop' : geneids_pop,
        'obo_dag':goeaobj.obo_dag}

def get_geneid2symbol(fin_xlsx):
    """Read xlsx file return dictionary with Entrez GeneID keys to Symbol data."""
    gene2symbol = {}
    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../tests/data/nbt_3102")
    tbl_genes = "{DIR}/{FIN}".format(DIR=data_dir, FIN=fin_xlsx)
    book = xlrd.open_workbook(tbl_genes)
    sheet = book.sheet_by_index(0)
    for rownum in range(sheet.nrows):
        symbol, geneid, _ = [sheet.cell_value(rownum, c) for c in range(sheet.ncols)]
        if geneid:
            gene2symbol[int(geneid)] = symbol
    return gene2symbol

def get_goeaobj(method, geneids_pop, taxid):
    """Load: ontologies, associations, and population geneids."""
    obo_dag = get_godag()
    assoc_geneid2gos = get_assoc_ncbi_taxids([taxid])
    goeaobj = GOEnrichmentStudy(
        geneids_pop,
        assoc_geneid2gos,
        obo_dag,
        propagate_counts=False,
        alpha=0.05,
        methods=[method])
     # obo_dag is also found in goeaobj.obo_dag
    return goeaobj


# Copyright (C) 2016-2017, DV Klopfenstein, H Tang. All rights reserved.
