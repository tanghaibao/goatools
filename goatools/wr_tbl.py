"""Contains generic table-writing functions. Data is expected to be a list of namedtuples.

       kwargs (kws):
           'prt_if': Only print a line if user-specfied test returns True.
               prt_if is a lambda function with the data item's namedtuple as input.
               Example: prt_if = lambda nt: nt.p_uncorrected < 0.05
           'sort_by' : User-customizable sort when printing.
               sortby is a lambda function with the data item's namedtuple as input.
               It is the 'key' used in the sorted function.
               Example: sort_by = lambda nt: [nt.NS, -1*nt.depth]
           'hdrs' : A list of column headers to use when printing the table.
               default: The fields in the data's namedtuple is used as the column headers.
           'sep': Separator used when printing the tab-separated table format.
               default: sep = '\t'
           'fld2col_widths: A dictionary of column widths used when writing xlsx files.
           'fld2fmt': Used in tsv files and xlsx files for formatting specific fields
           'print_names': Used when the user wants to print a subset of the nt fields.
               This may occur if a field is to be used in either sort_by or prt_if,
               but not to be printed.
"""

__copyright__ = "Copyright (C) 2016, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import re
import sys

def prt_txt(prt, data_nts, prtfmt, nt_fields=None, **kws):
    """Print list of namedtuples into a table using prtfmt."""
    if data_nts:
        if nt_fields is not None:
            _chk_flds_fmt(nt_fields, prtfmt)
        #_prt_txt_hdr(prt, prtfmt)
        if 'sort_by' in kws:
            data_nts = sorted(data_nts, key=kws['sort_by'])
        prt_if = kws['prt_if'] if 'prt_if' in kws else None
        for data_nt in data_nts:
            if prt_if is None or prt_if(data_nt):
                prt.write(prtfmt.format(**data_nt._asdict()))
    else:
        sys.stdout.write("      0 items. NOT WRITING w/format_string({F})\n".format(F=prtfmt))
        

def wr_xlsx(fout_xlsx, xlsx_data, **kws):
    """Write a spreadsheet into a xlsx file."""
    if xlsx_data:
        from xlsxwriter import Workbook
        workbook = Workbook(fout_xlsx)
        worksheet = workbook.add_worksheet()
        # Cell formatting
        fmt_hdr = workbook.add_format({'top':1, 'bottom':1, 'left':0, 'right':0, 'bold':True})
        fmt_txt = workbook.add_format({'border':0})
        flds = xlsx_data[0]._fields
        if 'fld2col_widths' in kws:
            _set_xlsx_colwidths(worksheet, kws['fld2col_widths'], flds)
        # Print header
        hdrs = flds if 'hdrs' not in kws else kws['hdrs']
        for col_idx, hdr in enumerate(hdrs):
            worksheet.write(0, col_idx, hdr, fmt_hdr)
        row_data_start = 1
        # Print data
        row_idx = _wrxlsxdata(xlsx_data, row_data_start, worksheet, fmt_txt, **kws)
        workbook.close()
        sys.stdout.write("  {:>5} items WROTE: {}\n".format(row_idx-row_data_start, fout_xlsx))
    else:
        sys.stdout.write("      0 items. NOT WRITING {}\n".format(fout_xlsx))

def _wrxlsxdata(xlsx_data, row_idx, worksheet, fmt_txt, **kws):
    """Write data into xlsx worksheet."""
    fld2fmt = None if 'fld2fmt' not in kws else kws['fld2fmt']
    prt_if = kws['prt_if'] if 'prt_if' in kws else None
    if 'sort_by' in kws:
        xlsx_data = sorted(xlsx_data, key=kws['sort_by'])
    for data_nt in xlsx_data:
        if prt_if is None or prt_if(data_nt):
            _wrxlsxrow(worksheet, row_idx, data_nt, fmt_txt, fld2fmt)
            row_idx += 1
    return row_idx

def _wrxlsxrow(worksheet, row_idx, data_nt, fmt_txt, fld2fmt):
    """Print an xlsx row by printing each column in order."""
    row_fld_vals = data_nt._asdict().items()
    for col_idx, (fld, val) in enumerate(row_fld_vals):
        # Optional user-formatting of specific fields, eg, pval: "{:8.2e}"
        if fld2fmt is not None and fld in fld2fmt:
            val = fld2fmt[fld].format(val)
        worksheet.write(row_idx, col_idx, val, fmt_txt)

def _set_xlsx_colwidths(worksheet, fld2col_widths, fldnames):
    """Set xlsx column widths using fld2col_widths."""
    for col_idx, fld in enumerate(fldnames):
        col_width = fld2col_widths.get(fld, None)
        if col_width is not None:
            worksheet.set_column(col_idx, col_idx, col_width)

def wr_tsv(fout_tsv, tsv_data, fldnames, **kws):
    """Write a file of tab-separated table data"""
    if tsv_data:
        ifstrm = sys.stdout if fout_tsv is None else open(fout_tsv, 'w')
        items = prt_tsv(ifstrm, tsv_data, fldnames, **kws)
        if fout_tsv is not None:
            sys.stdout.write("  {:>5} items WROTE: {}\n".format(items, fout_tsv))
            ifstrm.close()
    else:
        sys.stdout.write("      0 items. NOT WRITING {}\n".format(fout_tsv))

def prt_tsv(prt, data_nts, fldnames, **kws):
    """Print tab-separated table data"""
    # User-controlled printing options
    sep = "\t" if 'sep' not in kws else kws['sep']
    hdrs = fldnames if 'hdrs' not in kws else kws['hdrs']
    fld2fmt = None if 'fld2fmt' not in kws else kws['fld2fmt']
    if 'sort_by' in kws:
        data_nts = sorted(data_nts, key=kws['sort_by'])
    prt_if = kws['prt_if'] if 'prt_if' in kws else None
    # Write header
    prt.write("{}\n".format(sep.join(hdrs)))
    # Write data
    items = 0
    for nt_data_row in data_nts:
        if prt_if is None or prt_if(nt_data_row):
            if fld2fmt is not None:
                row_fld_vals = nt_data_row._asdict().items()
                row_vals = _fmt_fields(row_fld_vals, fld2fmt)
            else:
                row_vals = list(nt_data_row)
            prt.write("{}\n".format(sep.join(str(d) for d in row_vals)))
            items += 1
    return items

def _fmt_fields(fld_vals, fld2fmt):
    """Optional user-formatting of specific fields, eg, pval: '{:8.2e}'."""
    vals = []
    for fld, val in fld_vals:
        if fld in fld2fmt:
            val = fld2fmt[fld].format(val)
        vals.append(val)
    return vals

def _chk_flds_fmt(nt_fields, prtfmt):
    """Check that all fields in the prtfmt have corresponding data in the namedtuple."""
    fmtflds = get_fmtflds(prtfmt)
    missing_data = set(fmtflds).difference(set(nt_fields))
    # All data needed for print is present, return.
    if not missing_data:
        return
    #raise Exception('MISSING DATA({M}).'.format(M=" ".join(missing_data)))
    msg = ['CANNOT PRINT USING: "{PF}"'.format(PF=prtfmt.rstrip())]
    for fld in fmtflds:
        errmrk = "" if fld in nt_fields else "ERROR-->"
        msg.append("  {ERR:8} {FLD}".format(ERR=errmrk, FLD=fld)) 
    raise Exception('\n'.join(msg))

def get_fmtflds(prtfmt):
    """Return the fieldnames in the formatter text."""
    # Example prtfmt: "{NS} {study_cnt:2} {fdr_bh:5.3e} L{level:02} D{depth:02} {GO} {name}\n"
    return [f.split(':')[0] for f in re.findall(r'{(\S+)}', prtfmt)]

def get_fmtfldsdict(prtfmt):
    """Return the fieldnames in the formatter text."""
    # Example prtfmt: "{NS} {study_cnt:2} {fdr_bh:5.3e} L{level:02} D{depth:02} {GO} {name}\n"
    return {v:v for v in get_fmtflds(prtfmt)}

def _prt_txt_hdr(prt, prtfmt):
    """Print header for text report."""
    tblhdrs = get_fmtfldsdict(prtfmt)
    # If needed, reformat for format_string for header, which has strings, not floats.
    hdrfmt = re.sub(r':(\d+)\.\S+}', r':\1}', prtfmt)
    hdrfmt = re.sub(r':(0+)(\d+)}', r':\2}', hdrfmt)
    prt.write("#{}".format(hdrfmt.format(**tblhdrs)))

# Copyright (C) 2016, DV Klopfenstein, H Tang. All rights reserved.
