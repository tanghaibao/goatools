"""Contains generic table-writing functions. Data is expected to be a list of namedtuples.

       kwargs (kws):
           'title': First row will contain user-provided title string
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
           'prt_flds' : Used to print a subset of the fields in the namedtuple or
               to control the order of the print fields
           'fld2col_widths: A dictionary of column widths used when writing xlsx files.
           'fld2fmt': Used in tsv files and xlsx files for formatting specific fields

           For adding color or other formatting to a row based on value in a row:
               'ntfld_wbfmt': namedtuple field containing a value used as a key for a xlsx format
               'ntval2wbfmtdict': namedtuple value and corresponding xlsx format dict.
"""

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import re
import sys
from goatools.wr_tbl_class import get_hdrs

def prt_txt(prt, data_nts, prtfmt=None, nt_fields=None, **kws):
    """Print list of namedtuples into a table using prtfmt."""
    lines = get_lines(data_nts, prtfmt, nt_fields, **kws)
    if lines:
        for line in lines:
            prt.write(line)
    else:
        sys.stdout.write("      0 items. NOT WRITING\n")

def get_lines(data_nts, prtfmt=None, nt_fields=None, **kws):
    """Print list of namedtuples into a table using prtfmt."""
    lines = []
    # optional keyword args: prt_if sort_by
    if prtfmt is None:
        prtfmt = mk_fmtfld(data_nts[0], kws.get('joinchr', ' '), kws.get('eol', '\n'))
    # if nt_fields arg is None, use fields from prtfmt string.
    if nt_fields is not None:
        _chk_flds_fmt(nt_fields, prtfmt)
    if 'sort_by' in kws:
        data_nts = sorted(data_nts, key=kws['sort_by'])
    prt_if = kws.get('prt_if', None)
    for data_nt in data_nts:
        if prt_if is None or prt_if(data_nt):
            lines.append(prtfmt.format(**data_nt._asdict()))
    return lines

def prt_nts(data_nts, prtfmt=None, prt=sys.stdout, nt_fields=None, **kws):
    """Print list of namedtuples into a table using prtfmt."""
    prt_txt(prt, data_nts, prtfmt, nt_fields, **kws)

def wr_xlsx(fout_xlsx, data_xlsx, **kws):
    """Write a spreadsheet into a xlsx file."""
    from goatools.wr_tbl_class import WrXlsx
    # optional keyword args: fld2col_widths hdrs prt_if sort_by fld2fmt prt_flds
    items_str = kws.get("items", "items") if "items" not in kws else kws["items"]
    if data_xlsx:
        # Open xlsx file
        xlsxobj = WrXlsx(fout_xlsx, data_xlsx[0]._fields, **kws)
        worksheet = xlsxobj.add_worksheet()
        # Write title (optional) and headers.
        row_idx = xlsxobj.wr_title(worksheet)
        row_idx = xlsxobj.wr_hdrs(worksheet, row_idx)
        row_idx_data0 = row_idx
        # Write data
        row_idx = xlsxobj.wr_data(data_xlsx, row_idx, worksheet)
        # Close xlsx file
        xlsxobj.workbook.close()
        sys.stdout.write("  {N:>5} {ITEMS} WROTE: {FOUT}\n".format(
            N=row_idx-row_idx_data0, ITEMS=items_str, FOUT=fout_xlsx))
    else:
        sys.stdout.write("      0 {ITEMS}. NOT WRITING {FOUT}\n".format(
            ITEMS=items_str, FOUT=fout_xlsx))

def wr_xlsx_sections(fout_xlsx, xlsx_data, **kws):
    """Write xlsx file containing section names followed by lines of namedtuple data."""
    from goatools.wr_tbl_class import WrXlsx
    items_str = "items" if "items" not in kws else kws["items"]
    prt_hdr_min = 10
    num_items = 0
    if xlsx_data:
        # Basic data checks
        assert len(xlsx_data[0]) == 2, "wr_xlsx_sections EXPECTED: [(section, nts), ..."
        assert xlsx_data[0][1], \
            "wr_xlsx_sections EXPECTED SECTION({S}) LIST TO HAVE DATA".format(S=xlsx_data[0][0])
        # Open xlsx file and write title (optional) and headers.
        xlsxobj = WrXlsx(fout_xlsx, xlsx_data[0][1][0]._fields, **kws)
        worksheet = xlsxobj.add_worksheet()
        row_idx = xlsxobj.wr_title(worksheet)
        hdrs_wrote = False
        # Write data
        for section_text, data_nts in xlsx_data:
            num_items += len(data_nts)
            fmt = xlsxobj.wbfmtobj.get_fmt_section()
            row_idx = xlsxobj.wr_row_mergeall(worksheet, section_text, fmt, row_idx)
            if hdrs_wrote is False or len(data_nts) > prt_hdr_min:
                row_idx = xlsxobj.wr_hdrs(worksheet, row_idx)
                hdrs_wrote = True
            row_idx = xlsxobj.wr_data(data_nts, row_idx, worksheet)
        # Close xlsx file
        xlsxobj.workbook.close()
        sys.stdout.write("  {N:>5} {ITEMS} WROTE: {FOUT} ({S} sections)\n".format(
            N=num_items, ITEMS=items_str, FOUT=fout_xlsx, S=len(xlsx_data)))
    else:
        sys.stdout.write("      0 {ITEMS}. NOT WRITING {FOUT}\n".format(
            ITEMS=items_str, FOUT=fout_xlsx))

def wr_tsv(fout_tsv, tsv_data, **kws):
    """Write a file of tab-separated table data"""
    items_str = "items" if "items" not in kws else kws["items"]
    if tsv_data:
        ifstrm = sys.stdout if fout_tsv is None else open(fout_tsv, 'w')
        num_items = prt_tsv(ifstrm, tsv_data, **kws)
        if fout_tsv is not None:
            sys.stdout.write("  {N:>5} {ITEMS} WROTE: {FOUT}\n".format(
                N=num_items, ITEMS=items_str, FOUT=fout_tsv))
            ifstrm.close()
    else:
        sys.stdout.write("      0 {ITEMS}. NOT WRITING {FOUT}\n".format(
            ITEMS=items_str, FOUT=fout_tsv))

def prt_tsv_sections(prt, tsv_data, **kws):
    """Write tsv file containing section names followed by lines of namedtuple data."""
    prt_hdr_min = 10  # Print hdr on the 1st section and for any following 'large' sections
    num_items = 0
    if tsv_data:
        # Basic data checks
        assert len(tsv_data[0]) == 2, "wr_tsv_sections EXPECTED: [(section, nts), ..."
        assert tsv_data[0][1], \
            "wr_tsv_sections EXPECTED SECTION({S}) LIST TO HAVE DATA".format(S=tsv_data[0][0])
        hdrs_wrote = False
        sep = "\t" if 'sep' not in kws else kws['sep']
        prt_flds = kws['prt_flds'] if 'prt_flds' in kws else tsv_data[0]._fields
        fill = sep*(len(prt_flds) - 1)
        # Write data
        for section_text, data_nts in tsv_data:
            prt.write("{SEC}{FILL}\n".format(SEC=section_text, FILL=fill))
            if hdrs_wrote is False or len(data_nts) > prt_hdr_min:
                prt_tsv_hdr(prt, data_nts, **kws)
                hdrs_wrote = True
            num_items += prt_tsv_dat(prt, data_nts, **kws)
        return num_items
    else:
        return 0

def prt_tsv(prt, data_nts, **kws):
    """Print tab-separated table headers and data"""
    # User-controlled printing options
    prt_tsv_hdr(prt, data_nts, **kws)
    return prt_tsv_dat(prt, data_nts, **kws)

def prt_tsv_hdr(prt, data_nts, **kws):
    """Print tab-separated table headers"""
    sep = "\t" if 'sep' not in kws else kws['sep']
    flds_all = data_nts[0]._fields
    hdrs = get_hdrs(flds_all, **kws)
    prt.write("# {}\n".format(sep.join(hdrs)))

def prt_tsv_dat(prt, data_nts, **kws):
    """Print tab-separated table data"""
    sep = "\t" if 'sep' not in kws else kws['sep']
    fld2fmt = None if 'fld2fmt' not in kws else kws['fld2fmt']
    if 'sort_by' in kws:
        data_nts = sorted(data_nts, key=kws['sort_by'])
    prt_if = kws['prt_if'] if 'prt_if' in kws else None
    prt_flds = kws['prt_flds'] if 'prt_flds' in kws else data_nts[0]._fields
    items = 0
    for nt_data_row in data_nts:
        if prt_if is None or prt_if(nt_data_row):
            if fld2fmt is not None:
                row_fld_vals = [(fld, getattr(nt_data_row, fld)) for fld in prt_flds]
                row_vals = _fmt_fields(row_fld_vals, fld2fmt)
            else:
                row_vals = [getattr(nt_data_row, fld) for fld in prt_flds]
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

def mk_fmtfld(nt_item, joinchr=" ", eol="\n"):
    """Given a namedtuple, return a format_field string."""
    fldstrs = []
    # Default formats based on fieldname
    fld2fmt = {
        'hdrgo' : lambda f: "{{{FLD}:1,}}".format(FLD=f),
        'dcnt' : lambda f: "{{{FLD}:6,}}".format(FLD=f),
        'level' : lambda f: "L{{{FLD}:02,}}".format(FLD=f),
        'depth' : lambda f: "D{{{FLD}:02,}}".format(FLD=f),
    }
    for fld in nt_item._fields:
        if fld in fld2fmt:
            val = fld2fmt[fld](fld)
        else:
            val = "{{{FLD}}}".format(FLD=fld)
        fldstrs.append(val)
    return "{LINE}{EOL}".format(LINE=joinchr.join(fldstrs), EOL=eol)

# Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved.
