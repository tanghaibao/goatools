"""GOEA and report generation w/bonferroni multiple test corrections from statsmodels.

        python test_goea_rpt_bonferroni.py
        python test_goea_rpt_bonferroni.py [LOG FILENAME]
"""

__copyright__ = "Copyright (C) 2016-2017, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import sys

from goatools.obo_parser import GODag
from goatools.associations import read_associations
from goatools.go_enrichment import GOEnrichmentStudy

def test_bonferroni(fout_log=None):
    """Do Gene Ontology Enrichment Analysis w/Bonferroni multipletest. Print results 3 ways."""
    # ---------------------------------------------------------------------
    # Run Gene Ontology Analysis (GOEA)
    #
    # 1. Initialize
    log = sys.stdout if fout_log is None else open(fout_log, 'w')
    results_nt, goea = run_bonferroni(log)

    # ---------------------------------------------------------------------
    # Print results 3 ways: to screen, to tsv (tab-separated file), to xlsx (Excel spreadsheet)
    fout_tsv = "goea_bonferroni.tsv"
    fout_xls = "goea_bonferroni.xlsx"
    
    # print these in tsv and xlsx
    print_fields = ['NS', 'study_count', 'p_uncorrected', 'p_bonferroni', 'level', 'depth', 'GO', 'name'] 
    # 1. Print results to screen using format in prtfmt. For example:
    #
    #      BP 22 3.073e-03 L06 D07 GO:0006468 protein phosphorylation
    #      BP  9 1.023e-02 L07 D08 GO:0006511 ubiquitin-dependent protein catabolic process
    #      BP  2 1.023e-02 L05 D09 GO:0019877 diaminopimelate biosynthetic process
    #      BP  2 1.223e-02 L04 D08 GO:0006301 postreplication repair
    #      BP  2 1.223e-02 L05 D09 GO:0030418 nicotianamine biosynthetic process
    #      BP  2 1.492e-02 L04 D06 GO:0006909 phagocytosis
    #      BP  2 1.492e-02 L03 D03 GO:0051322 anaphase
    #      ...
    prtfmt = " ".join(["{NS} {study_count:3} {p_uncorrected:5.3e}",
             "{p_bonferroni:5.3e} L{level:02} D{depth:02} {GO} {name}\n"])
    prt_if = lambda nt: nt.p_uncorrected < 0.05
    goea.prt_txt(log, results_nt, prtfmt, prt_if=prt_if)

    # 2. Write results to tsv file
    # Optional user defined formatting for specific fields
    fld2fmt = {'p_bonferroni':'{:8.2e}', 'p_uncorrected':'{:8.2e}'} 
    # Sort by: 1st) BP, MF, CC; 2nd) By GO depth, deepest GO first.
    sort_by = lambda nt: [nt.NS, -1*nt.depth] 
    goea.wr_tsv(fout_tsv, results_nt,
        prt_if=prt_if, sort_by=sort_by, fld2fmt=fld2fmt, prt_flds=print_fields)

    # 3. Write results to xlsx file, including specific study genes assc. w/significant GOs
    # Use these headers instead of the print_fields for the xlsx header
    hdrs = ['NS', 'pval', 'bonferroni', 'L', 'D', 'Term', 'Ontology Term Name', 'Cnt', 'Genes']
    print_fields = ['NS', 'p_uncorrected', 'p_bonferroni', 'level', 'depth', 'GO', 'name', 'study_count', 'study_items']
    goea.wr_xlsx(fout_xls, results_nt, 
        # optional key-word args (ie, kwargs, kws)
        prt_if=prt_if, sort_by=sort_by, hdrs=hdrs, fld2fmt=fld2fmt, prt_flds=print_fields) 
    if fout_log is not None:
        log.close()
        sys.stdout.write("  WROTE: {}\n".format(fout_log))

def run_bonferroni(log):
    """Do Gene Ontology Enrichment Analysis w/Bonferroni multipletest. Print results 3 ways."""
    # ---------------------------------------------------------------------
    # Run Gene Ontology Analysis (GOEA)
    #
    # 1. Initialize
    obo_dag = GODag("go-basic.obo")
    assoc = read_associations("../data/association", no_top=True)
    popul_ids = [line.rstrip() for line in open("../data/population")]
    study_ids = [line.rstrip() for line in open("../data/study")]
    # 2. Run enrichment analysis
    goea = GOEnrichmentStudy(popul_ids, assoc, obo_dag, alpha=0.05, methods=['bonferroni'])
    results_nt = goea.run_study(study_ids)
    return results_nt, goea

if __name__ == '__main__':
    test_bonferroni()

# Copyright (C) 2016-2017, DV Klopfenstein, H Tang. All rights reserved.
