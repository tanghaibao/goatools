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
           'prt_flds' : Used to print a subset of the fields in the namedtuple or
               to control the order of the print fields
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

def prt_txt(prt, data_nts, prtfmt=None, nt_fields=None, **kws):
    """Print list of namedtuples into a table using prtfmt."""
    # optional keyword args: prt_if sort_by
    if data_nts:
        if prtfmt is None:
            prtfmt = mk_fmtfld(data_nts[0])
        # if nt_fields arg is None, use fields from prtfmt string.
        if nt_fields is not None:
            _chk_flds_fmt(nt_fields, prtfmt)
        if 'sort_by' in kws:
            data_nts = sorted(data_nts, key=kws['sort_by'])
        prt_if = kws['prt_if'] if 'prt_if' in kws else None
        for data_nt in data_nts:
            if prt_if is None or prt_if(data_nt):
                prt.write(prtfmt.format(**data_nt._asdict()))
    else:
        sys.stdout.write("      0 items. NOT WRITING w/format_string({F})\n".format(F=prtfmt))

def prt_nts(data_nts, prtfmt=None, prt=sys.stdout, nt_fields=None, **kws):
    """Print list of namedtuples into a table using prtfmt."""
    prt_txt(prt, data_nts, prtfmt, nt_fields, **kws)


class WrXlsx(object):
    """Class to store/manage Excel spreadsheet parameters."""

    dflt_fmt_hdr = {'top':1, 'bottom':1, 'left':0, 'right':0, 'bold':True}
    dflt_fmt_txt = [
        # Colors from "Grey Grey Color Palette" http://www.color-hex.com/color-palette/20477
        {'border':0},
        # single-line, weight=3, bg_color=light grey
        {'top':0, 'bottom':0, 'left':0, 'right':0, 'bold':True, 'bg_color':'#eeeded'},
        # double-line, weight=3
        {'top':0, 'bottom':6, 'left':0, 'right':0, 'bold':True, 'bg_color':'#c7c7c7'}]

    def __init__(self, fout_xlsx, xlsx_data, **kws):
        # KEYWORDS FOR WRITING DATA:
        self.fld2fmt = None if 'fld2fmt' not in kws else kws['fld2fmt']
        # User may specify to skip rows based on values in row
        self.prt_if = kws['prt_if'] if 'prt_if' in kws else None
        # User may specify a subset of columns to print or
        # a column ordering different from the _fields seen in the namedtuple
        self.b_format_txt = None
        self.prt_flds = None
        self._init_prt_flds(xlsx_data, **kws) # Init: b_format_txt, prtflds
        # Initialize sort_by
        self.sort_by = kws['sort_by'] if 'sort_by' in kws else None
        # Workbook
        from xlsxwriter import Workbook
        self.workbook = Workbook(fout_xlsx)
        self.fmt_hdr = self.workbook.add_format(self.dflt_fmt_hdr)
        self.fmt_txt = self._init_fmt_txt(**kws)

    def _init_prt_flds(self, xlsx_data, **kws):
        """Initialize prt_flds."""
        # Get user-defined fields or all namedtuple fields
        prt_flds = kws['prt_flds'] if 'prt_flds' in kws else xlsx_data[0]._fields
        b_format_txt = 'format_txt' in prt_flds
        # "format_txt" is used to set row colors, but values are not printed in columns.
        if b_format_txt:
            prt_flds = [f for f in prt_flds if f != "format_txt"]
        # Initialize data members
        self.b_format_txt = b_format_txt
        self.prt_flds = prt_flds

    def add_worksheet(self):
        """Add a worksheet to the workbook."""
        return self.workbook.add_worksheet()

    def get_fmt_txt(self, idx):
        """Return format for text cell."""
        assert idx < 3
        return self.fmt_txt[idx]

    def _init_fmt_txt(self, **kws):
        """Save user-provided cell format, if provided. If not, use default."""
        fmt_vals = [
            kws['format_txt0'] if 'format_txt0' in kws else self.dflt_fmt_txt[0],
            kws['format_txt1'] if 'format_txt1' in kws else self.dflt_fmt_txt[1],
            kws['format_txt2'] if 'format_txt2' in kws else self.dflt_fmt_txt[2]]
        return [
            self.workbook.add_format(fmt_vals[0]),
            self.workbook.add_format(fmt_vals[1]),
            self.workbook.add_format(fmt_vals[2])]

    @staticmethod
    def set_xlsx_colwidths(worksheet, fld2col_widths, fldnames):
        """Set xlsx column widths using fld2col_widths."""
        for col_idx, fld in enumerate(fldnames):
            col_width = fld2col_widths.get(fld, None)
            if col_width is not None:
                worksheet.set_column(col_idx, col_idx, col_width)

def wr_xlsx(fout_xlsx, xlsx_data, **kws):
    """Write a spreadsheet into a xlsx file."""
    # optional keyword args: fld2col_widths hdrs prt_if sort_by fld2fmt prt_flds
    if xlsx_data:
        xlsxobj = WrXlsx(fout_xlsx, xlsx_data, **kws)
        worksheet = xlsxobj.add_worksheet()
        flds_all = xlsxobj.prt_flds
        if 'fld2col_widths' in kws:
            xlsxobj.set_xlsx_colwidths(worksheet, kws['fld2col_widths'], flds_all)
        row_idx = 0
        # Print header
        hdrs = [h for h in _get_hdrs(flds_all, **kws) if h != "format_txt"]
        if 'title' in kws:
            worksheet.merge_range(row_idx, 0, row_idx, len(hdrs)-1, kws['title'], xlsxobj.fmt_hdr)
            row_idx += 1
        # Print header
        for col_idx, hdr in enumerate(hdrs):
            worksheet.write(row_idx, col_idx, hdr, xlsxobj.fmt_hdr)
        row_idx += 1
        row_idx_d0 = row_idx
        # Print data
        row_idx = _wrxlsxdata(xlsx_data, row_idx, worksheet, xlsxobj)
        xlsxobj.workbook.close()
        sys.stdout.write("  {:>5} items WROTE: {}\n".format(row_idx-row_idx_d0, fout_xlsx))
    else:
        sys.stdout.write("      0 items. NOT WRITING {}\n".format(fout_xlsx))

def _get_hdrs(flds_all, **kws):
    """Return headers, given user-specified key-word args."""
    # Return Headers if the user explicitly lists them.
    if 'hdrs' in kws:
        return kws['hdrs']
    # User may specify a subset of fields or a column order using prt_flds
    if 'prt_flds' in kws:
        return kws['prt_flds']
    # All fields in the namedtuple will be in the headers
    return flds_all

def _wrxlsxdata(xlsx_data, row_idx, worksheet, xlsxobj):
    """Write data into xlsx worksheet."""
    fld2fmt = xlsxobj.fld2fmt
    # User may specify to skip rows based on values in row
    prt_if = xlsxobj.prt_if
    # User may specify a subset of columns to print or
    # a column ordering different from the _fields seen in the namedtuple
    prt_flds = xlsxobj.prt_flds
    # User may specify a sort function or send the data list in already sorted
    b_format_txt = xlsxobj.b_format_txt
    if xlsxobj.sort_by is not None:
        xlsx_data = sorted(xlsx_data, key=xlsxobj.sort_by)
    for data_nt in xlsx_data:
        fmt_txt = xlsxobj.get_fmt_txt(0 if not b_format_txt else getattr(data_nt, "format_txt"))
        if prt_if is None or prt_if(data_nt):
            # Print an xlsx row by printing each column in order.
            for col_idx, fld in enumerate(prt_flds):
                # If field "format_txt" is present, use value for formatting, but don't print.
                val = getattr(data_nt, fld, "")
                # Optional user-formatting of specific fields, eg, pval: "{:8.2e}"
                if fld2fmt is not None and fld in fld2fmt:
                    val = fld2fmt[fld].format(val)
                worksheet.write(row_idx, col_idx, val, fmt_txt)
            row_idx += 1
    return row_idx

def wr_tsv(fout_tsv, tsv_data, **kws):
    """Write a file of tab-separated table data"""
    if tsv_data:
        ifstrm = sys.stdout if fout_tsv is None else open(fout_tsv, 'w')
        items = prt_tsv(ifstrm, tsv_data, **kws)
        if fout_tsv is not None:
            sys.stdout.write("  {:>5} items WROTE: {}\n".format(items, fout_tsv))
            ifstrm.close()
    else:
        sys.stdout.write("      0 items. NOT WRITING {}\n".format(fout_tsv))

def prt_tsv(prt, data_nts, **kws):
    """Print tab-separated table data"""
    # User-controlled printing options
    sep = "\t" if 'sep' not in kws else kws['sep']
    flds_all = data_nts[0]._fields
    hdrs = _get_hdrs(flds_all, **kws)
    fld2fmt = None if 'fld2fmt' not in kws else kws['fld2fmt']
    if 'sort_by' in kws:
        data_nts = sorted(data_nts, key=kws['sort_by'])
    prt_if = kws['prt_if'] if 'prt_if' in kws else None
    prt_flds = kws['prt_flds'] if 'prt_flds' in kws else data_nts[0]._fields
    # Write header
    prt.write("{}\n".format(sep.join(hdrs)))
    # Write data
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

def mk_fmtfld(nt_item):
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
    fldstrs.append("\n")
    return " ".join(fldstrs)

# Copyright (C) 2016, DV Klopfenstein, H Tang. All rights reserved.
