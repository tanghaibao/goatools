#!/usr/bin/env python
"""Run a Gene Ontology Enrichment Analysis (GOEA), plots, etc.

    Nature 2014_0126; 
			Computational analysis of cell-to-cell heterogeneity
      in single-cell RNA-sequencing data reveals hidden
      subpopulations of cells
    http://www.nature.com/nbt/journal/v33/n2/full/nbt.3102.html#methods

		     ... revealed a significant enrichment in the set
         of 401 genes that were differentially expressed
         between the identified clusters (P = 0.001
         Hypergeometric Test). Further, Gene Ontology (GO)
         enrichment analysis showed that the differentially
         expressed genes contained statistically
         significant enrichments of genes involved in:
             * glycolysis 
             * cellular response to IL-4 stimulation
               NOW: BP GO:0071353 1.668e-03 D06 cellular response to interleukin-4 (5 genes)
                  * BP GO:0070670: response to interleukin-4
                  * BP GO:0071353: cellular response to interleukin-4
             * positive regulation of B-cell proliferation 
               NOW: BP GO:0030890 2.706e-04 D09 positive regulation of B cell proliferation (7 genes)

         * 401 genes: Supplementary table 4
           http://www.nature.com/nbt/journal/v33/n2/extref/nbt.3102-S4.xlsx
         * GO enrichment results are in: Supplementary table 6
           http://www.nature.com/nbt/journal/v33/n2/extref/nbt.3102-S6.xlsx
"""
import os
import sys
import xlrd

from PyBiocode.dnld.NCBI.genes_NCBI_mus_ProteinCoding import GeneID2nt as GeneID2nt_mus
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.associations import get_assoc_ncbi_taxids
from goatools.godag_plot import plot_gos, plot_results, plot_goid2goobj

__copyright__ = "Copyright (C) 2016, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

def test_example(log=sys.stdout):
    """Run Gene Ontology Enrichment Analysis (GOEA) on Nature data."""
    # Mouse study
    taxid = 10090
    # Load ontologies, associations, and population ids
    geneids_pop = GeneID2nt_mus.keys()
    geneids_study = get_geneid2symbol("nbt.3102-S4_GeneIDs.xlsx", log)
    goeaobj = get_goeaobj("fdr_bh", geneids_pop, taxid)
    # Run GOEA on study
    keep_if = lambda nt: nt.p_fdr_bh < 0.05 # keep if results are significant
    goea_results = goeaobj.run_study(geneids_study, keep_if=keep_if)
    geneids = goeaobj.get_study_items(goea_results)
    # Print GOEA results to files
    goeaobj.wr_xlsx("nbt3102.xlsx", goea_results)
    goeaobj.wr_txt("nbt3102.txt", goea_results)
    # Plot significant GO terms w/annotated study info
    # with a variety of different plot options.
    plot_results("nbt3102_{NS}.png", goea_results)
    plot_results("nbt3102_{NS}_sym.png", goea_results, study_items=5, items_p_line=2, id2symbol=geneids_study)
    goid_subset = [
        'GO:0006096', # BP 4.24e-12 10 glycolytic process
        'GO:0071353', # BP 7.45e-06  5 cellular response to interleukin-4
        'GO:0030890', # BP 8.22e-07  7 positive regulation of B cell proliferation
    ]
    obo = goeaobj.obo_dag
    plot_gos("nbt3102_GOs.png", goid_subset, obo)
    plot_gos("nbt3102_GOs_genecnt.png", goid_subset, obo, goea_results=goea_results)
    plot_gos("nbt3102_GOs_genelst.png", goid_subset, obo, 
        study_items=True, goea_results=goea_results)
    plot_gos("nbt3102_GOs_symlst.png", goid_subset, obo, 
        study_items=True, id2symbol=geneids_study, goea_results=goea_results)
    plot_gos("nbt3102_GOs_symlst_small.png", goid_subset, obo, 
        study_items=5, id2symbol=geneids_study, goea_results=goea_results)
    plot_gos("nbt3102_GOs_GO0005743.png", ["GO:0005743"], obo, 
        items_p_line=2, study_items=6, 
        id2symbol=geneids_study, goea_results=goea_results)
    # One GO sub-plot per significant GO term from study
    for rec in goea_results:
        png = "nbt3102_{GO}.png".format(GO=rec.GO)
        goid2obo = {rec.GO:rec.goterm}
        plot_goid2goobj(png,
            goid2obo, # source GOs and their GOTerm object
            study_items=15, # Max number of gene symbols to print in each GO term
            id2symbol=geneids_study, # Contains GeneID-to-Symbol
            goea_results=goea_results) # pvals used for GO Term coloring
    # Are any significant geneids related to cell cycle?
    import test_genes_cell_cycle as CC
    genes_cell_cycle = CC.get_genes_cell_cycle(taxid, log=None)
    genes_cell_cycle_sig = genes_cell_cycle.intersection(geneids)
    CC.prt_genes("nbt3102_cell_cycle.txt", genes_cell_cycle_sig, taxid, log)


def get_goeaobj(method, geneids_pop, taxid):
    """Load: ontologies, associations, and population geneids."""
    fin_obo = "go-basic.obo"
    if not os.path.isfile(fin_obo):
        os.system("wget http://geneontology.org/ontology/go-basic.obo") 
    obo_dag = GODag(fin_obo)
    assoc_geneid2gos = get_assoc_ncbi_taxids([taxid])
    goeaobj = GOEnrichmentStudy(
        geneids_pop,
        assoc_geneid2gos,
        obo_dag,
        propagate_counts = False,
        alpha = 0.05,
        methods = [method])
    return goeaobj

def get_geneid2symbol(fin_xlsx, log):
    """Read xlsx file."""
    gene2symbol = {}
    data_dir = os.path.dirname(os.path.abspath(__file__)) + "/data/nbt_3102"
    tbl_genes = "{DIR}/{FIN}".format(DIR=data_dir, FIN=fin_xlsx)
    book = xlrd.open_workbook(tbl_genes)
    pg = book.sheet_by_index(0)
    for r in range(pg.nrows):
        symbol, geneid, pval = [pg.cell_value(r, c) for c in range(pg.ncols)]
        if geneid:
            gene2symbol[int(geneid)] = symbol
    return gene2symbol
    
if __name__ == '__main__':
    test_example()

# Copyright (C) 2016, DV Klopfenstein, H Tang, All rights reserved.
