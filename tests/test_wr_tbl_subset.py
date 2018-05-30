"""Show two methods to write the same data into: xlsx file, tsv file, and text file."""

import os
import sys

import filecmp
from goatools.wr_tbl import get_fmtflds, wr_xlsx, wr_tsv
from goatools.test_data.nature3102_goea import get_goea_results
from goatools.rpt.goea_nt_xfrm import get_goea_nts_prt


def test_wr_methods(log=sys.stdout):
    """Demonstrate printing a subset of all available fields using two methods."""
    # 1. Gene Ontology Enrichment Analysis
    #    1a. Initialize: Load ontologies, associations, and population gene IDs
    nature_data = get_goea_results()
    goeaobj = nature_data['goeaobj']
    goea_results = nature_data['goea_results']
    # 2. Write results
    #    Write parameters:
    #    The format_string names below are the same names as in the namedtuple field_names.
    prtfmt = "{GO} {NS} {level:>2} {depth:>2} {p_fdr_bh:5.2e} {study_count:>5} {name}\n"
    wr_params = {
        # Format for printing in text format
        'prtfmt' : prtfmt,
        # Format for p-values in tsv and xlsx format
        'fld2fmt' : {'p_fdr_bh' : '{:8.2e}'},
        # Print a subset namedtuple fields, don't print all fields in namedtuple.
        'prt_flds' : get_fmtflds(prtfmt)
    }
    #    2a. Use the write functions inside the GOEnrichmentStudy class.
    cwddir = os.getcwd()
    tsv_obj = os.path.join(cwddir, 'nbt3102_subset_obj.tsv')
    tsv_nts = os.path.join(cwddir, 'nbt3102_subset_nt.tsv')
    _wr_3fmt_goeaobj(tsv_obj, goea_results, goeaobj, wr_params, log)
    #    2b. Use the write functions from the wr_tbl package to print a list of namedtuples.
    _wr_3fmt_wrtbl(tsv_nts, goea_results, wr_params, log)
    assert filecmp.cmp(tsv_obj, tsv_nts)



def _wr_3fmt_goeaobj(tsv_obj, goea_results, goeaobj, wr_params, log):
    """Demonstrate printing a subset of GOEA fields."""
    # List of all fields, printable or not, available from GOEA
    flds = goea_results[0].get_prtflds_all()
    log.write("\nALL GOEA FIELDS: {F}\n".format(F=" ".join(flds)))
    # Use the subset of namedtuple fields_names that are listed in the format string:
    # Same format: print to screen and print to file:
    goeaobj.prt_txt(log, goea_results, **wr_params) # Print to screen
    goeaobj.wr_txt("nbt3102_subset_obj.txt", goea_results, **wr_params)
    # Print to Excel Spreadsheet
    title = "Print subset of GOEA fields"
    goeaobj.wr_xlsx("nbt3102_subset_obj.xlsx", goea_results, title=title, **wr_params)
    # Print to tab-separated file
    goeaobj.wr_tsv(tsv_obj, goea_results, **wr_params)


def _wr_3fmt_wrtbl(tsv_nts, goea_results, wr_params, log):
    """Using functions in the wr_tbl pkg, demonstrate printing a subset of namedtuple fields."""
    goea_nts = get_goea_nts_prt(goea_results)
    # List of all fields, printable or not, in namedtuple (nt):
    log.write("\nnamedtuple FIELDS: {F}\n".format(F=" ".join(goea_nts[0]._fields)))
    # Print to Excel Spreadsheet
    title = "Print subset of namedtuple fields"
    wr_xlsx("nbt3102_subset_nt.xlsx", goea_nts, title=title, **wr_params)
    # Print to tab-separated file
    wr_tsv(tsv_nts, goea_nts, **wr_params)


if __name__ == '__main__':
    test_wr_methods(sys.stdout)
