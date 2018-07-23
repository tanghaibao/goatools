"""Helper Class and functions for writing tables."""

import sys


class WrXlsxParams(object):
    """Class containing user-specified or default parameters."""

    def __init__(self, **kws):
        self.kws = kws
        self.title = kws.get('title', None)
        self.fld2fmt = kws.get('fld2fmt', None)
        self.fld2col_widths = kws.get('fld2col_widths', None)
        # User may specify to skip rows based on values in row
        self.prt_if = kws.get('prt_if', None)
        # Initialize sort_by
        self.sort_by = kws.get('sort_by', None)


class WrXlsx(object):
    """Class to store/manage Excel spreadsheet parameters."""

    dflt_fmt_hdr = {'top':0, 'bottom':0, 'left':0, 'right':0, 'bold':True}

    def __init__(self, fout_xlsx, nt_flds, **kws):
        # KEYWORDS FOR WRITING DATA:
        self.vars = WrXlsxParams(**kws)
        # Workbook
        from xlsxwriter import Workbook
        self.workbook = Workbook(fout_xlsx)
        self.wbfmtobj = WbFmt(nt_flds, self.workbook, **kws)
        self.hdrs = self.wbfmtobj.get_hdrs(**kws)
        self.fmt_hdr = self.workbook.add_format(self.dflt_fmt_hdr)

    def wr_title(self, worksheet, row_idx=0):
        """Write title (optional)."""
        if self.vars.title is not None:
            # Title is one line
            if isinstance(self.vars.title, str):
                return self.wr_row_mergeall(worksheet, self.vars.title, self.fmt_hdr, row_idx)
            # Title is multi-line
            else:
                ridx = row_idx
                for title_line in self.vars.title:
                    ridx = self.wr_row_mergeall(worksheet, title_line, self.fmt_hdr, ridx)
                return ridx
        return row_idx

    def wr_row_mergeall(self, worksheet, txtstr, fmt, row_idx):
        """Merge all columns and place text string in widened cell."""
        hdridxval = len(self.hdrs) - 1
        worksheet.merge_range(row_idx, 0, row_idx, hdridxval, txtstr, fmt)
        return row_idx + 1

    def wr_hdrs(self, worksheet, row_idx):
        """Print row of column headers"""
        for col_idx, hdr in enumerate(self.hdrs):
            # print("ROW({R}) COL({C}) HDR({H}) FMT({F})\n".format(
            #     R=row_idx, C=col_idx, H=hdr, F=self.fmt_hdr))
            worksheet.write(row_idx, col_idx, hdr, self.fmt_hdr)
        row_idx += 1
        return row_idx

    def wr_data(self, xlsx_data, row_i, worksheet):
        """Write data into xlsx worksheet."""
        fld2fmt = self.vars.fld2fmt
        # User may specify to skip rows based on values in row
        prt_if = self.vars.prt_if
        # User may specify a subset of columns to print or
        # a column ordering different from the _fields seen in the namedtuple
        prt_flds = self.wbfmtobj.get_prt_flds()
        get_wbfmt = self.wbfmtobj.get_wbfmt
        if self.vars.sort_by is not None:
            xlsx_data = sorted(xlsx_data, key=self.vars.sort_by)
        try:
            for data_nt in xlsx_data:
                if prt_if is None or prt_if(data_nt):
                    wbfmt = get_wbfmt(data_nt)  # xlsxwriter.format.Format created w/add_format
                    # Print an xlsx row by printing each column in order.
                    for col_i, fld in enumerate(prt_flds):
                        try:
                            # If fld "format_txt" present, use val for formatting, but don't print.
                            val = getattr(data_nt, fld, "")
                            # Optional user-formatting of specific fields, eg, pval: "{:8.2e}"
                            # If field value is empty (""), don't use fld2fmt
                            if fld2fmt is not None and fld in fld2fmt and val != "" and val != "*":
                                val = fld2fmt[fld].format(val)
                            worksheet.write(row_i, col_i, val, wbfmt)
                        except:
                            raise RuntimeError(self._get_err_msg(row_i, col_i, fld, val, prt_flds))
                    row_i += 1
        except RuntimeError as inst:
            import traceback
            traceback.print_exc()
            sys.stderr.write("\n  **FATAL in wr_data: {MSG}\n\n".format(MSG=str(inst)))
            sys.exit(1)
        return row_i

    @staticmethod
    def _get_err_msg(row, col, fld, val, prt_flds):
        """Return an informative message with details of xlsx write attempt."""
        import traceback
        traceback.print_exc()
        err_msg = (
            "ROW({R}) COL({C}) FIELD({F}) VAL({V})\n".format(R=row, C=col, F=fld, V=val),
            "PRINT FIELDS({N}): {F}".format(N=len(prt_flds), F=" ".join(prt_flds)))
        return "\n".join(err_msg)

    def add_worksheet(self):
        """Add a worksheet to the workbook."""
        wsh = self.workbook.add_worksheet()
        if self.vars.fld2col_widths is not None:
            self.set_xlsx_colwidths(wsh, self.vars.fld2col_widths, self.wbfmtobj.get_prt_flds())
        return wsh

    @staticmethod
    def set_xlsx_colwidths(worksheet, fld2col_widths, fldnames):
        """Set xlsx column widths using fld2col_widths."""
        for col_idx, fld in enumerate(fldnames):
            col_width = fld2col_widths.get(fld, None)
            if col_width is not None:
                worksheet.set_column(col_idx, col_idx, col_width)


class WbFmt(object):
    """Manages Workbook format objects."""

    dflt_wbfmtdict = [
        # Colors from "Grey Grey Color Palette" http://www.color-hex.com/color-palette/20477
        # http://www.color-hex.com/color/808080
        # http://www.color-hex.com/color/d3d3d3
        {'border':0},
        {'top':0, 'bottom':0, 'left':0, 'right':0, 'bold':True, 'bg_color':'#eeeded'},
        {'top':0, 'bottom':0, 'left':0, 'right':0, 'bold':True, 'bg_color':'#d3d3d3'},
        {'border':1, 'bold':True}]

        # Save user-provided cell format, if provided. If not, use default.
    def __init__(self, nt_flds, workbook, **kws):
        self.prt_flds = kws.get('prt_flds', nt_flds)
        self.b_format_txt = 'format_txt' in self.prt_flds
        self.ntfld_wbfmt = kws.get('ntfld_wbfmt', None)
        self.fmtname2wbfmtobj = self._init_fmtname2wbfmtobj(workbook, **kws)
        self.b_plain = not self.b_format_txt and self.ntfld_wbfmt is None

    def get_prt_flds(self):
        """Remove namedtuple fields used for formatting rows, but not printed."""
        return [f for f in self.prt_flds if f != "format_txt"]

    def get_hdrs(self, **kws):
        """Initialize column headers."""
        hdrs = get_hdrs(self.prt_flds, **kws)
        # Values in a "format_txt" "column" are used for formatting, not printing
        return [h for h in hdrs if h != "format_txt"]

    def _init_fmtname2wbfmtobj(self, workbook, **kws):
        """Initialize fmtname2wbfmtobj."""
        wbfmtdict = [
            kws.get('format_txt0', self.dflt_wbfmtdict[0]),
            kws.get('format_txt1', self.dflt_wbfmtdict[1]),
            kws.get('format_txt2', self.dflt_wbfmtdict[2]),
            kws.get('format_txt3', self.dflt_wbfmtdict[3])]
        fmtname2wbfmtobj = {
            'plain': workbook.add_format(wbfmtdict[0]),
            'plain bold': workbook.add_format(wbfmtdict[3]),
            'very light grey' : workbook.add_format(wbfmtdict[1]),
            'light grey' :workbook.add_format(wbfmtdict[2])}
        # Use a xlsx namedtuple field value to set row color
        ntval2wbfmtdict = kws.get('ntval2wbfmtdict', None)
        if ntval2wbfmtdict is not None:
            for ntval, wbfmtdict in ntval2wbfmtdict.items():
                fmtname2wbfmtobj[ntval] = workbook.add_format(wbfmtdict)
            if 'ntfld_wbfmt' not in kws:
                sys.stdout.write("**WARNING: 'ntfld_wbfmt' NOT PRESENT\n")
        return fmtname2wbfmtobj

    def get_wbfmt(self, data_nt=None):
        """Return format for text cell."""
        if data_nt is None or self.b_plain:
            return self.fmtname2wbfmtobj.get('plain')
        # User namedtuple field/value for color
        if self.ntfld_wbfmt is not None:
            return self.__get_wbfmt_usrfld(data_nt)
        # namedtuple format_txt for color/bold/border
        if self.b_format_txt:
            wbfmt = self.__get_wbfmt_format_txt(data_nt)
            if wbfmt is not None:
                return wbfmt
        # 'ntfld_wbfmt': namedtuple field which contains a value used as a key for a xlsx format
        # 'ntval2wbfmtdict': namedtuple value and corresponding xlsx format dict.
        return self.fmtname2wbfmtobj.get('plain')

    def __get_wbfmt_usrfld(self, data_nt):
        """Return format for text cell from namedtuple field specified by 'ntfld_wbfmt'"""
        if self.ntfld_wbfmt is not None:
            if isinstance(self.ntfld_wbfmt, str):
                ntval = getattr(data_nt, self.ntfld_wbfmt, None) # Ex: 'section'
                if ntval is not None:
                    return self.fmtname2wbfmtobj.get(ntval, None)
            #### elif isinstance(self.ntfld_wbfmt, dict):
            ####     print("DDDDDDDDDDDD IIIIIIIII CCCCCCCCCC TTTTTTTTTTTTTTTT")

    def __get_wbfmt_format_txt(self, data_nt):
        """Return format for text cell from namedtuple field, 'format_txt'."""
        format_txt_val = getattr(data_nt, "format_txt")
        if format_txt_val == 1:
            return self.fmtname2wbfmtobj.get("very light grey")
        if format_txt_val == 2:
            return self.fmtname2wbfmtobj.get("light grey")
        return self.fmtname2wbfmtobj.get(format_txt_val)

    def get_fmt_section(self):
        """Grey if printing header GOs and plain if not printing header GOs."""
        if self.b_format_txt:
            return self.fmtname2wbfmtobj.get("light grey")
        return self.fmtname2wbfmtobj.get("plain bold")


def get_hdrs(flds_all, **kws):
    """Return headers, given user-specified key-word args."""
    # Return Headers if the user explicitly lists them.
    hdrs = kws.get('hdrs', None)
    if hdrs is not None:
        return hdrs
    # User may specify a subset of fields or a column order using prt_flds
    if 'prt_flds' in kws:
        return kws['prt_flds']
    # All fields in the namedtuple will be in the headers
    return flds_all

#  Copyright (C) 2015-2018. All rights reserved.
