"""Manages print format for GOEA results."""

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

# import collections as cx


class PrtFmt(object):
    """Manages print format for GOEA results."""

    # Default Excel table column widths for GOEA results
    default_fld2col_widths = {
        'NS'        :  3,
        'GO'        : 12,
        'alt'       :  2,
        'level'     :  3,
        'depth'     :  3,
        'enrichment':  1,
        'name'      : 60,
        'ratio_in_study':  8,
        'ratio_in_pop'  : 12,
        'study_items'   : 15,
    }

    default_fld2fmt = {
        'NS'        : '{NS}',
        'GO'        : '{GO}',
        'alt'       : '{alt:1}',
        'level'     : 'L{level:02}',
        'depth'     : 'D{depth:02}',
        'reldepth'  : 'D{reldepth:02}',
        'dcnt'      : '{dcnt:5}',
        'tcnt'      : '{tcnt:7,}',
        'tfreq'     : '{tfreq:8.6f}',
        'tinfo'     : '{tinfo:5.2f}',
        'enrichment': '{enrichment:1}',
        'name'      : '{name:40}',
        'ratio_in_study': '{ratio_in_study:>7}',
        'ratio_in_pop'  : '{ratio_in_pop:>11}',
        'study_count'   : '{study_count:4}',
        'study_items'   : '{study_items}',
    }

    prtfmt_dflt = ("{GO} {NS} {p_uncorrected:5.2e} {ratio_in_study:>6} {ratio_in_pop:>9} "
                   "{depth:02} {name:40} {study_items}\n")

    def __init__(self):
        pass

    def get_prtfmt_str(self, flds, add_nl=True):
        fmts = self.get_prtfmt_list(flds, add_nl)
        return " ".join(fmts) 

    def get_prtfmt_list(self, flds, add_nl=True):
        """Get print format, given fields."""
        fmts = []
        for fld in flds:
            if fld[:2] == 'p_':
                fmts.append('{{{FLD}:8.2e}}'.format(FLD=fld))
            elif fld in self.default_fld2fmt:
                fmts.append(self.default_fld2fmt[fld])
            else:
                raise Exception("UNKNOWN FORMAT: {FLD}".format(FLD=fld))
        if add_nl:
            fmts.append("\n")
        return fmts

    @staticmethod
    def adjust_prtfmt(prtfmt):
        """Adjust format_strings for legal values."""
        prtfmt = prtfmt.replace("{p_holm-sidak", "{p_holm_sidak")
        prtfmt = prtfmt.replace("{p_simes-hochberg", "{p_simes_hochberg")
        return prtfmt

# Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved.
