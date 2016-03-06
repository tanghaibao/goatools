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
             * GO:0006096 "glycolysis" or "glycolytic process"
             * cellular response to IL-4 stimulation
               NOW: BP GO:0071353 1.668e-03 D06 cellular response to interleukin-4 (5 genes)
                  * BP GO:0070670: response to interleukin-4
                  * BP GO:0071353: cellular response to interleukin-4
             * positive regulation of B-cell proliferation 
               NOW: BP GO:0030890 2.706e-04 D09 positive regulation of B cell proliferation (7 genes)

         * 401 genes: Supplementary table 4
           note: Total gene count is 400 genes, not 401: Rpl41 is listed twice
           http://www.nature.com/nbt/journal/v33/n2/extref/nbt.3102-S4.xlsx
         * GO enrichment results are in: Supplementary table 6
           http://www.nature.com/nbt/journal/v33/n2/extref/nbt.3102-S6.xlsx
"""
import os
import sys
import xlrd
import wget
from collections import Counter, defaultdict, OrderedDict

from goatools.test_data.genes_NCBI_10090_ProteinCoding import GeneID2nt as GeneID2nt_mus
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.associations import get_assoc_ncbi_taxids
from goatools.godag_plot import plot_gos, plot_results, plot_goid2goobj

__copyright__ = "Copyright (C) 2016, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

def test_example(log=sys.stdout):
    """Run Gene Ontology Enrichment Analysis (GOEA) on Nature data."""
    # --------------------------------------------------------------------
    # --------------------------------------------------------------------
    # Gene Ontology Enrichment Analysis (GOEA)
    # --------------------------------------------------------------------
    # --------------------------------------------------------------------
    taxid = 10090 # Mouse study
    # Load ontologies, associations, and population ids
    geneids_pop = GeneID2nt_mus.keys()
    geneids_study = get_geneid2symbol("nbt.3102-S4_GeneIDs.xlsx", log)
    goeaobj = get_goeaobj("fdr_bh", geneids_pop, taxid)
    # Run GOEA on study
    keep_if = lambda nt: getattr(nt, "p_fdr_bh" ) < 0.05 # keep if results are significant
    goea_results = goeaobj.run_study(geneids_study, keep_if=keep_if)
    compare_results(goea_results)
    geneids = goeaobj.get_study_items(goea_results)
    # Print GOEA results to files
    goeaobj.wr_xlsx("nbt3102.xlsx", goea_results)
    goeaobj.wr_txt("nbt3102.txt", goea_results)
    # Plot all significant GO terms w/annotated study info (large plots)
    #plot_results("nbt3102_{NS}.png", goea_results)
    #plot_results("nbt3102_{NS}_sym.png", goea_results, study_items=5, items_p_line=2, id2symbol=geneids_study)



    # --------------------------------------------------------------------
    # --------------------------------------------------------------------
    # Further examination of GOEA results...
    # --------------------------------------------------------------------
    # --------------------------------------------------------------------
    obo = goeaobj.obo_dag
    dpi = 150 # For review: Figures can be saved in .jpg, .gif, .tif or .eps, at 150 dpi


    # --------------------------------------------------------------------
    # Item 1) Words in GO names associated with large numbers of study genes
    # --------------------------------------------------------------------
    # What GO term words are associated with the largest number of study genes?
    prt_word2genecnt("nbt3102_genecnt_GOword.txt", goea_results, log)
    # Curated selection of GO words associated with large numbers of study genes
    freq_seen = ['RNA', 'translation', 'mitochondrial', 'mitochondrion', 'ribosomal', 'ribosome',
        'adhesion', 'endoplasmic', 'nucleotide', 'apoptotic', 'myelin']
    # Collect the GOs which contains the chosen frequently seen words
    word2gos = get_word2gos(freq_seen, goea_results)
    go2res = {nt.GO:nt for nt in goea_results}
    # Print words of interest, the sig GO terms which contain that word, and study genes.
    prt_word_GO_genes("nbt3102_GO_word_genes.txt", word2gos, go2res, geneids_study, log)
    # Plot each set of GOs along w/study gene info 
    for word, gos in word2gos.items():
       goid2obo = {go:go2res[go].goterm for go in gos}
       plot_goid2goobj(
           "nbt3102_word_{WORD}.png".format(WORD=word),
           goid2obo, # source GOs and their GOTerm object
           study_items=15, # Max number of gene symbols to print in each GO term
           id2symbol=geneids_study, # Contains GeneID-to-Symbol
           goea_results=goea_results, # pvals used for GO Term coloring
           dpi=dpi)
      
    
    # --------------------------------------------------------------------
    # Item 2) Explore findings of Nature paper:
    #
    #     Gene Ontology (GO) enrichment analysis showed that the
    #     differentially expressed genes contained statistically
    #     significant enrichments of genes involved in 
    #         glycolysis,
    #         cellular response to IL-4 stimulation and 
    #         positive regulation of B-cell proliferation
    # --------------------------------------------------------------------
    goid_subset = [
        'GO:0006096', # BP 4.24e-12 10 glycolytic process
        'GO:0071353', # BP 7.45e-06  5 cellular response to interleukin-4
        'GO:0030890', # BP 8.22e-07  7 positive regulation of B cell proliferation
    ]
    plot_gos("nbt3102_GOs.png", goid_subset, obo, dpi=dpi)
    plot_gos("nbt3102_GOs_genecnt.png", goid_subset, obo, goea_results=goea_results, dpi=dpi)
    plot_gos("nbt3102_GOs_genelst.png", goid_subset, obo, 
        study_items=True, goea_results=goea_results, dpi=dpi)
    plot_gos("nbt3102_GOs_symlst.png", goid_subset, obo, 
        study_items=True, id2symbol=geneids_study, goea_results=goea_results, dpi=dpi)
    plot_gos("nbt3102_GOs_symlst_trunc.png", goid_subset, obo, 
        study_items=5, id2symbol=geneids_study, goea_results=goea_results, dpi=dpi)
    plot_gos("nbt3102_GOs_GO0005743.png", ["GO:0005743"], obo, 
        items_p_line=2, study_items=6, 
        id2symbol=geneids_study, goea_results=goea_results, dpi=dpi)


    # --------------------------------------------------------------------
    # Item 3) Create one GO sub-plot per significant GO term from study
    # --------------------------------------------------------------------
    for rec in goea_results:
        png = "nbt3102_{NS}_{GO}.png".format(GO=rec.GO.replace(':', '_'), NS=rec.NS)
        goid2obo = {rec.GO:rec.goterm}
        plot_goid2goobj(png,
            goid2obo, # source GOs and their GOTerm object
            study_items=15, # Max number of gene symbols to print in each GO term
            id2symbol=geneids_study, # Contains GeneID-to-Symbol
            goea_results=goea_results, # pvals used for GO Term coloring
            dpi=dpi)


    # --------------------------------------------------------------------
    # Item 4) Explore using manually curated lists of GO terms
    # --------------------------------------------------------------------
    goid_subset = [
      'GO:0030529', # CC D03 intracellular ribonucleoprotein complex (42 genes)
      'GO:0015934', # CC D05 large ribosomal subunit (4 genes)
      'GO:0015935', # CC D05 small ribosomal subunit (13 genes)
      'GO:0022625', # CC D06 cytosolic large ribosomal subunit (16 genes)
      'GO:0022627', # CC D06 cytosolic small ribosomal subunit (19 genes)
      'GO:0036464', # CC D06 cytoplasmic ribonucleoprotein granule (4 genes)
      'GO:0005840', # CC D05 ribosome (35 genes)
      'GO:0005844', # CC D04 polysome (6 genes)
    ]
    plot_gos("nbt3102_CC_ribosome.png", goid_subset, obo, 
        study_items=6, id2symbol=geneids_study, items_p_line=3,
        goea_results=goea_results, dpi=dpi)

    goid_subset = [
      'GO:0003723', # MF D04 RNA binding (32 genes)
      'GO:0044822', # MF D05 poly(A) RNA binding (86 genes)
      'GO:0003729', # MF D06 mRNA binding (11 genes)
      'GO:0019843', # MF D05 rRNA binding (6 genes)
      'GO:0003746', # MF D06 translation elongation factor activity (5 genes)
    ]
    for dpi in [150, 1200]: # 150 for review, 1200 for publication
        plot_gos("nbt3102_MF_RNA_dpi{DPI}.png".format(DPI=dpi), 
            goid_subset, 
            obo, 
            study_items=6, id2symbol=geneids_study, items_p_line=3,
            goea_results=goea_results, dpi=dpi)

    # --------------------------------------------------------------------
    # Item 5) Are any significant geneids related to cell cycle?
    # --------------------------------------------------------------------
    import test_genes_cell_cycle as CC
    genes_cell_cycle = CC.get_genes_cell_cycle(taxid, log=log)
    genes_cell_cycle_sig = genes_cell_cycle.intersection(geneids)
    CC.prt_genes("nbt3102_cell_cycle.txt", genes_cell_cycle_sig, taxid, log=None)



def get_goeaobj(method, geneids_pop, taxid):
    """Load: ontologies, associations, and population geneids."""
    fin_obo = "go-basic.obo"
    if not os.path.isfile(fin_obo):
        wget.download("wget http://geneontology.org/ontology/go-basic.obo") 
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
    """Read xlsx file return dictionary with Entrez GeneID keys to Symbol data."""
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

def compare_results(goea_results):
    """Compare GOATOOLS to results from Nature paper."""
    act_goids = [rec.GO for rec in goea_results]
    exp_goids = set(paper_top20())
    overlap = set(act_goids).intersection(exp_goids)
    fout_txt = "nbt3102_compare.txt"
    with open(fout_txt, 'w') as prt:
        prt.write("{N} GO terms overlapped with {M} top 20 GO terms in Nature paper\n".format(
            N = len(overlap), M = len(exp_goids)))
        idx = 1
        gos = set()
        for rec in goea_results:
            if rec.GO in exp_goids:
                gos.add(rec.GO)
                sig = '*' if rec.p_fdr_bh < 0.05 else ' '
                prt.write("{I:>2} {NS} {SIG} {GO} D{D:>02} {ALPHA:5.2e} {NAME}({N} genes)\n".format(
                    I=idx, NS=rec.NS, D=rec.goterm.depth, GO=rec.GO, NAME=rec.name, 
                    ALPHA=rec.p_fdr_bh, SIG=sig, N=rec.study_count))
                idx += 1
        nogo = exp_goids.difference(gos)
        prt.write("NOT LISTED: {GO}\n".format(GO=", ".join(nogo)))
   

def prt_word2genecnt(fout, goea_results, log):
    """Get words in GO term names and the number of study genes associated with GO words."""
    word2genes = defaultdict(set)
    for rec in goea_results:
        study_items = rec.study_items
        for word in rec.name.split():
            word2genes[word] |= study_items
    word2genecnt = Counter({w:len(gs) for w, gs in word2genes.items()})
    with open(fout, "w") as wordstrm:
        for word, cnt in word2genecnt.most_common():
            wordstrm.write("{CNT:>3} {WORD}\n".format(CNT=cnt, WORD=word))
    log.write("  WROTE: {F}\n".format(F=fout))
       
def get_word2gos(words, goea_results):
    """Get all GO terms which contain a word in 'words'."""        
    word2gos = defaultdict(set)
    for word in words: 
        for rec in goea_results: 
            if word in rec.name:
                word2gos[word].add(rec.GO)
    return OrderedDict([(w, word2gos[w]) for w in words])

def prt_word_GO_genes(fout, word2gos, go2res, geneids_study, log):
    """Print words in GO names that have large numbers of study genes."""
    with open(fout, "w") as prt:
      prt.write("""This file is generated by test_nbt3102.py and is intended to confirm
this statement in the GOATOOLS manuscript:

        We observed: 
            93 genes associated with RNA, 
            47 genes associated with translation, 
            38 genes associated with mitochondrial, and 
            37 genes associated with ribosomal, as reported by GOATOOLS.

""")
      for word, gos in word2gos.items():
        genes = set()
        for go in gos:
          genes |= go2res[go].study_items
        genes = sorted([geneids_study[g] for g in genes])
        prt.write("\n{WD}: {N} study genes, {M} GOs\n".format(WD=word, N=len(genes), M=len(gos)))
        for i, go in enumerate(gos):
          res = go2res[go]
          prt.write("{I}) {GO} {NAME} ({N} genes)\n".format(I=i, GO=go, NAME=res.name, N=res.study_count))
        prt.write("{N} study genes:\n".format(N=len(genes)))
        N = 10 # 10 genes per line
        mult = [genes[i:i+N] for i in range(0, len(genes), N)]
        prt.write("  {}\n".format("\n  ".join([", ".join(str(g) for g in sl) for sl in mult])))
      log.write("  WROTE: {F}\n".format(F=fout))

def paper_top20():
    """Return top 20 GO terms in Nature paper found using R's topGO with alogrithm=elim.

       Supplemental Table 6 description from supplemental information:
       ---------------------------------------------------------------

			     GO enrichment of differentially expressed genes
           between the sub-populations of cells revealed by non-linear
           PCA on the scLVM corrected expression levels. The R-package
           topGO was used with the "elim" algorithm and the top 20
           terms are shown.

           http://www.nature.com/nbt/journal/v33/n2/extref/nbt.3102-S1.pdf


       GOEA: How topGO tests and GOATOOLS test methods differ:
       -------------------------------------------------------

       The test statistics supported by topGO when using the topGO "elim" algorithm are
       [fisher, ks, t, globaltest, sum]. When using the "elim" algorithm, topGO
       does not automatically do multiple-test correction. 
       The documentation for topGO says:
       https://www.bioconductor.org/packages/3.3/bioc/vignettes/topGO/inst/doc/topGO.pdf

           For the methods that account for the GO topology like elim
           and weight, the problem of multiple testing is even more
           complicated. Here one computes the p-value of a GO term
           conditioned on the neighbouring terms. The tests are
           therefore not independent and the multiple testing theory
           does not directly apply. We like to interpret the p-values
           returned by these methods as corrected or not affected by
           multiple testing.

       For the GOATOOLS publication, we use the following:
           1.   GOATOOLS version 0.5.10
           2.   "Fisher's exact test" statistical analysis method
           2a.  "Benjamini/Hochberg" for multiple-test correction
           3a.  Ontologies from go-basic.obo version 1.2 release 2016-02-06
           3b.  Annotations from ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz modified 2/27/16
           4.   The population is 28,212 protein-coding genes for mouse
           4a.      18,396 of the 28,212 population genes contain GO annotations
           5.   The study size is the 400 genes in supplemental table 4 (Rpl41 is listed twice).
           5a.         372 of the 400 study genes contain GO annotations

     
    """
    return [
        # GO ID         idx Term                                  Annotated Sig Exp  result1
        # ---------     --- -----------------------------         --------- --- ---- -------
        "GO:0006412", #   1 translation                                 403 45  7.88 9.5e-12
        "GO:0006414", #   2 translational elongation                     44 11  0.86 5.3e-10
        "GO:0000028", #   3 ribosomal small subunit assembly              9  6  0.18 4.0e-09
        "GO:0006096", #   4 glycolysis                                   55 11  1.08 6.8e-09
        "GO:0071353", #   5 cellular response to interleukin-4           20  6  0.39 1.5e-06
        "GO:0030890", #   6 positive regulation of B cell proliferat...  37  7  0.72 6.0e-06
        "GO:0006172", #   7 ADP biosynthetic process                      8  4  0.16 9.1e-06
        "GO:0051099", #   8 positive regulation of binding               85  9  1.66 3.9e-05
        "GO:0008637", #   9 apoptotic mitochondrial changes              66  8  1.29 3.9e-05
        "GO:0051129", #  10 negative regulation of cellular componen... 308 21  6.02 0.00011
        "GO:0002474", #  11 antigen processing and presentation of p...  19  6  0.37 0.00012
        "GO:0046835", #  12 carbohydrate phosphorylation                 14  4  0.27 0.00012
        "GO:0042273", #  13 ribosomal large subunit biogenesis           14  4  0.27 0.00012
        "GO:0043066", #  14 negative regulation of apoptotic process    584 31 11.42 0.00016
        "GO:0043029", #  15 T cell homeostasis                           28  5  0.55 0.00018
        "GO:0015986", #  16 ATP synthesis coupled proton transport       16  4  0.31 0.00021
        "GO:0042274", #  17 ribosomal small subunit biogenesis           20 10  0.39 0.00022
        "GO:0030388", #  18 fructose 1,6-bisphosphate metabolic proc...   7  3  0.14 0.00024
        "GO:1901385", #  19 regulation of voltage-gated calcium chan...   7  3  0.14 0.00024
        "GO:0042102", #  20 positive regulation of T cell proliferat...  66  7  1.29 0.00028
    ]
    
if __name__ == '__main__':
    test_example()

# Copyright (C) 2016, DV Klopfenstein, H Tang, All rights reserved.
