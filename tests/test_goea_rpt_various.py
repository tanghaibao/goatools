"""GOEA and reports for various multipletest corrections from local and statsmodels.

        python test_goea_rpt_various.py
        python test_goea_rpt_various.py [LOG FILENAME]
"""

__copyright__ = "Copyright (C) 2016, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import sys

sys.path.insert(0, '..') # Use local version of goatools during test
from goatools.obo_parser import GODag
from goatools.associations import read_associations

from PyBiocode.enrichanal.enrichanal_GO import GOEA

def test_bonferronis(log, goeaobj):
    current_methods = ['local_bonferroni', 'statsmodels_bonferroni']
    goeaobj.set_params(methods=current_methods)
    study_ids = [line.rstrip() for line in open("../data/study")]
    run_goea(study_ids, "bonferronis", goeaobj, log)

def run_goea(study_ids, fout_base, goea, log):
    results_nt = goea.find_enrichment(study_ids)

    # ---------------------------------------------------------------------
    # Print results 3 ways: to screen, to tsv (tab-separated file), to xlsx (Excel spreadsheet)
    fout_tsv = "goea_{BASE}.tsv".format(BASE=fout_base)
    fout_xls = "goea_{BASE}.xlsx".format(BASE=fout_base)

    # Print these in tsv and xlsx in this order
    print_names = ['NS', 'study_cnt', 'pval_uncor', 
              'local_bonferroni_star',       'local_bonferroni', 
        'statsmodels_bonferroni_star', 'statsmodels_bonferroni', 
        'level', 'depth', 'GO', 'name'] 
    # Collect this. Used in prt_if to only print significant GO terms
    field_names = print_names + ['local_bonferroni_sig', 'statsmodels_bonferroni_sig'] 
    # *_sig is True or False: Print the GOEA GO Term only if it is significant
    prt_if = lambda nt: nt.local_bonferroni_sig or nt.statsmodels_bonferroni_sig

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
    # Print format field names are the same names as in the "field_names" variable.
    prtfmt = " ".join(["{NS} {study_cnt:2} {pval_uncor:5.3e}",
             "{local_bonferroni_star:1} {local_bonferroni:5}",
             "{statsmodels_bonferroni_star:1} {statsmodels_bonferroni:5.3e}",
             "L{level:02} D{depth:02} {GO} {name}\n"])
    goea.prt_txt(log, results_nt, field_names, prtfmt, prt_if=prt_if)

    ## 2. Write results to tsv file
    ## Optional user defined formatting for specific fields
    #fld2fmt = {'bonferroni':'{:8.2e}', 'pval_uncor':'{:8.2e}'} 
    ## Sort by: 1st) BP, MF, CC; 2nd) By GO depth, deepest GO first.
    #sort_by = lambda nt: [nt.NS, -1*nt.depth] 
    #goea.wr_tsv(fout_tsv, results_nt, field_names, 
    #    prt_if=prt_if, sort_by=sort_by, fld2fmt=fld2fmt, print_names=print_names)

    ## 3. Write results to xlsx file
    ## Use these headers instead of the print_names for the xlsx header
    ## TBD Check that header and size of fields printed match
    #goea.wr_xlsx(fout_xls, results_nt, field_names, 
    #    # optional key-word args (ie, kwargs, kws)
    #    prt_if=prt_if, sort_by=sort_by, fld2fmt=fld2fmt, print_names=print_names) 

def init_goea(log):
    """Read Ontologies and Annotations once."""
    # ---------------------------------------------------------------------
    # Run Gene Ontology Analysis (GOEA)
    #
    # 1. Initialize
    obo_dag = GODag("go-basic.obo")
    assoc = read_associations("../data/association", no_top=True)
    popul_ids = [line.rstrip() for line in open("../data/population")]
    # 2. Run enrichment analysis
    goeaobj = GOEA(obo_dag, assoc, log)
    goeaobj.set_population(popul_ids)
    return goeaobj

def close_log(fout_log, log):
    """Close log file and print name of log file."""
    if fout_log is not None:
        log.close()
        sys.stdout.write("  WROTE: {}\n".format(fout_log))

def run_all():
    """Run a variety of enrichment analysis."""
    fout_log = GOEA.get_fout_log()
    log = sys.stdout if fout_log is None else open(fout_log, 'w')
    goeaobj = init_goea(log) # Initialize once
    test_bonferronis(log, goeaobj)
    close_log(fout_log, log)

if __name__ == '__main__':
    run_all()

# Copyright (C) 2016, DV Klopfenstein, H Tang. All rights reserved.
