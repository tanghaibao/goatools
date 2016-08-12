"""Helper Class and functions for writing tables."""

class WrXlsxParams(object):
    """Class containing user-specified or default parameters."""

    dflt_fmt_txt = [
        # Colors from "Grey Grey Color Palette" http://www.color-hex.com/color-palette/20477
        # http://www.color-hex.com/color/808080
        # http://www.color-hex.com/color/d3d3d3
        {'border':0},
        {'top':0, 'bottom':0, 'left':0, 'right':0, 'bold':True, 'bg_color':'#eeeded'},
        {'top':0, 'bottom':0, 'left':0, 'right':0, 'bold':True, 'bg_color':'#d3d3d3'},
        {'border':1, 'bold':True}]

    def __init__(self, nt_flds, **kws):
        self.title = kws['title'] if 'title' in kws else None
        self.fld2fmt = kws['fld2fmt'] if 'fld2fmt' in kws else None
        self.fld2col_widths = kws['fld2col_widths'] if 'fld2col_widths' in kws else None
        # User may specify to skip rows based on values in row
        self.prt_if = kws['prt_if'] if 'prt_if' in kws else None
        # Initialize sort_by
        self.sort_by = kws['sort_by'] if 'sort_by' in kws else None
        # User may specify a subset of columns to print or
        # a column ordering different from the _fields seen in the namedtuple
        self.b_format_txt = None
        self.prt_flds = None
        self._init_prt_flds(kws['prt_flds'] if 'prt_flds' in kws else nt_flds)
        self.hdrs = self._init_hdrs(**kws)
        # Save user-provided cell format, if provided. If not, use default.
        self.fmt_txt_arg = [
            kws['format_txt0'] if 'format_txt0' in kws else self.dflt_fmt_txt[0],
            kws['format_txt1'] if 'format_txt1' in kws else self.dflt_fmt_txt[1],
            kws['format_txt2'] if 'format_txt2' in kws else self.dflt_fmt_txt[2],
            kws['format_txt3'] if 'format_txt3' in kws else self.dflt_fmt_txt[3]]

    def _init_prt_flds(self, prt_flds):
        """Initialize prt_flds."""
        b_format_txt = 'format_txt' in prt_flds
        # "format_txt" is used to set row colors, but values are not printed in columns.
        if b_format_txt:
            prt_flds = [f for f in prt_flds if f != "format_txt"]
        # Initialize data members
        self.b_format_txt = b_format_txt
        self.prt_flds = prt_flds

    def _init_hdrs(self, **kws):
        """Initialize column headers."""
        hdrs = get_hdrs(self.prt_flds, **kws)
        # Values in a "format_txt" "column" are used for formatting, not printing
        return [h for h in hdrs if h != "format_txt"]

class WrXlsx(object):
    """Class to store/manage Excel spreadsheet parameters."""

    dflt_fmt_hdr = {'top':0, 'bottom':0, 'left':0, 'right':0, 'bold':True}

    def __init__(self, fout_xlsx, nt_flds, **kws):
        # KEYWORDS FOR WRITING DATA:
        self.vars = WrXlsxParams(nt_flds, **kws)
        # Workbook
        from xlsxwriter import Workbook
        self.workbook = Workbook(fout_xlsx)
        self.fmt_hdr = self.workbook.add_format(self.dflt_fmt_hdr)
        self.fmt_txt = {
            'plain': self.workbook.add_format(self.vars.fmt_txt_arg[0]),
            'plain bold': self.workbook.add_format(self.vars.fmt_txt_arg[3]),
            'very light grey' : self.workbook.add_format(self.vars.fmt_txt_arg[1]),
            'light grey' :self.workbook.add_format(self.vars.fmt_txt_arg[2])}

    def wr_title(self, worksheet, row_idx=0):
        """Write title (optional)."""
        if self.vars.title is not None:
            return self.wr_row_mergeall(worksheet, self.vars.title, self.fmt_hdr, row_idx)
        return row_idx

    def wr_row_mergeall(self, worksheet, txtstr, fmt, row_idx):
        """Merge all columns and place text string in widened cell."""
        hdridxval = len(self.vars.hdrs) - 1
        worksheet.merge_range(row_idx, 0, row_idx, hdridxval, txtstr, fmt)
        return row_idx + 1

    def wr_hdrs(self, worksheet, row_idx):
        """Print row of column headers"""
        for col_idx, hdr in enumerate(self.vars.hdrs):
            worksheet.write(row_idx, col_idx, hdr, self.fmt_hdr)
        row_idx += 1
        return row_idx

    def wr_data(self, xlsx_data, row_idx, worksheet):
        """Write data into xlsx worksheet."""
        fld2fmt = self.vars.fld2fmt
        # User may specify to skip rows based on values in row
        prt_if = self.vars.prt_if
        # User may specify a subset of columns to print or
        # a column ordering different from the _fields seen in the namedtuple
        prt_flds = self.vars.prt_flds
        if self.vars.sort_by is not None:
            xlsx_data = sorted(xlsx_data, key=self.vars.sort_by)
        for data_nt in xlsx_data:
            #fmt_txt = self._get_fmt_txt(0 if not b_format_txt else getattr(data_nt, "format_txt"))
            fmt_txt = self._get_fmt_txt(data_nt)
            if prt_if is None or prt_if(data_nt):
                # Print an xlsx row by printing each column in order.
                for col_idx, fld in enumerate(prt_flds):
                    # If field "format_txt" is present, use value for formatting, but don't print.
                    val = getattr(data_nt, fld, "")
                    # Optional user-formatting of specific fields, eg, pval: "{:8.2e}"
                    # If field value is empty (""), don't use fld2fmt
                    if fld2fmt is not None and fld in fld2fmt and val != "":
                        val = fld2fmt[fld].format(val)
                    worksheet.write(row_idx, col_idx, val, fmt_txt)
                row_idx += 1
        return row_idx

    def add_worksheet(self):
        """Add a worksheet to the workbook."""
        worksheet = self.workbook.add_worksheet()
        if self.vars.fld2col_widths is not None:
            self.set_xlsx_colwidths(worksheet, self.vars.fld2col_widths, self.vars.prt_flds)
        return worksheet

    def _get_fmt_txt(self, data_nt=None):
        """Return format for text cell."""
        if data_nt is None or not self.vars.b_format_txt:
            return self.fmt_txt.get('plain')
        format_txt_val = getattr(data_nt, "format_txt")
        if format_txt_val == 1:
            return self.fmt_txt.get("very light grey")
        if format_txt_val == 2:
            return self.fmt_txt.get("light grey")
        fmt = self.fmt_txt.get(format_txt_val)
        if fmt is not None:
            return fmt
        return self.fmt_txt.get('plain')

    def get_fmt_section(self):
        """Grey if printing header GOs and plain if not printing header GOs."""
        if self.vars.b_format_txt:
            return self.fmt_txt.get("light grey")
        return self.fmt_txt.get("plain bold")

    @staticmethod
    def set_xlsx_colwidths(worksheet, fld2col_widths, fldnames):
        """Set xlsx column widths using fld2col_widths."""
        for col_idx, fld in enumerate(fldnames):
            col_width = fld2col_widths.get(fld, None)
            if col_width is not None:
                worksheet.set_column(col_idx, col_idx, col_width)

def get_hdrs(flds_all, **kws):
    """Return headers, given user-specified key-word args."""
    # Return Headers if the user explicitly lists them.
    if 'hdrs' in kws:
        return kws['hdrs']
    # User may specify a subset of fields or a column order using prt_flds
    if 'prt_flds' in kws:
        return kws['prt_flds']
    # All fields in the namedtuple will be in the headers
    return flds_all

