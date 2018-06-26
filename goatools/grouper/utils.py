"""Utilities for grouper modules."""

__copyright__ = "Copyright (C) 2016-2017, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import collections as cx


def get_hdridx_flds():
    """Return "hdr idx" named tuple fields used in Excel spreadsheets."""
    return ["format_txt", "hdr_idx", "is_hdrgo", "is_usrgo", "num_usrgos", "hdr1usr01"]

def no_duplicates_sections2d(sections2d, prt=None):
    """Check for duplicate header GO IDs in the 2-D sections variable."""
    no_dups = True
    ctr = cx.Counter()
    for _, hdrgos in sections2d:
        for goid in hdrgos:
            ctr[goid] += 1
    for goid, cnt in ctr.most_common():
        if cnt == 1:
            break
        no_dups = False
        if prt is not None:
            prt.write("**SECTIONS WARNING FOUND: {N:3} {GO}\n".format(N=cnt, GO=goid))
    return no_dups


# Copyright (C) 2016-2017, DV Klopfenstein, H Tang, All rights reserved.
