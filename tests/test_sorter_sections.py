#!/usr/bin/env python
"""Test various options while sorting when using sections."""

import os
from goatools.base import get_godag
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.grouper.grprdflts import GrouperDflts
from goatools.grouper.hdrgos import HdrgosSections
from goatools.grouper.grprobj import Grouper
from goatools.grouper.sorter import Sorter
from goatools.grouper.wrxlsx import WrXlsxSortedGos
from goatools.test_data.sorter import USER_GOS
from goatools.test_data.sorter import SECTIONS

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def test_sort():
    """Test various options while sorting when using sections."""
    grprobj = _get_grprobj()
    # kws: hdrgo_prt section_prt use_sections

    #------------------------------------------------------------
    # TEST: hdrgo_prt and section_prt
    #------------------------------------------------------------
    # 1) Print sections:
    #   * Print headers: ph(print header) ch(color header)
    _wr_xlsx("t_a1_ps1_ph1", grprobj)
    #   * Omit headers (retain "a_ph1" sort)
    _wr_xlsx("t_a2_ps1_ph0", grprobj, hdrgo_prt=False)
    # 2) Print section as an additional column instead of a section row.
    #   * Use and print headers.
    _wr_xlsx("t_a3_ps0_ph1", grprobj, section_prt=False)
    #   * Use (but don't print) headers.
    _wr_xlsx("t_a4_ps0_ph0", grprobj, section_prt=False, hdrgo_prt=False)

    #------------------------------------------------------------
    # TEST: use_sections=False
    #       sections are ignored, but hdrgos defined in sections are used
    #------------------------------------------------------------
    _wr_xlsx("t_b1_us0_ph1", grprobj, use_sections=False)
    _wr_xlsx("t_b2_us0_ph0", grprobj, use_sections=False, hdrgo_prt=False)

    #------------------------------------------------------------
    # TEST: use_sections=True hdrgo_prt(T/F)
    #   These conditions force hdrgo_prt = False, if hdrgo_prt was not mentioned
    #       * section_sortby == True
    #       * section_sortby = user_sort
    #       * top_n == N
    #------------------------------------------------------------
    sortby = lambda nt: nt.depth
    _wr_xlsx("t_c1_ps1_ph0_hsortT", grprobj, section_sortby=True)
    _wr_xlsx("t_c2_ps1_ph0_hsortT_top3", grprobj, section_sortby=True, top_n=3)
    # Not commonly used: Uses order from a1_ps1_ph1
    _wr_xlsx("t_c3_ps1_ph0_hsortX_top5", grprobj, top_n=5)
    # Most commonly used; User provides the sort. Users often like to sort by pval, if exists:
    _wr_xlsx("t_c4_ps1_ph0_usort", grprobj, section_sortby=sortby)
    _wr_xlsx("t_c5_ps1_ph0_usort_top3", grprobj, section_sortby=sortby, top_n=3)
    _wr_xlsx("t_c6_ps0_ph0_usort_top3", grprobj, section_sortby=sortby, top_n=3, section_prt=False)


def _wr_xlsx(name, grprobj, **kws):
    """Group, sort, and print xlsx file."""
    # print('\nTEST {} kws_sortobj: {}'.format(name, kws))
    # KWS SORT OBJ
    kws_sort = {'sortby', 'hdrgo_sortby', 'section_sortby'}
    # KWS SORT FUNC: hdrgo_prt section_prt top_n use_sections prtfmt
    # Exclude ungrouped "Misc." section of sections var(sec_rd)
    fout_xlsx = "{NAME}.xlsx".format(NAME=name)
    # kws Sorter: hdrgo_prt section_prt top_n use_sections
    sortobj = Sorter(grprobj, **{k:v for k, v in kws.items() if k in kws_sort})
    desc2nts = sortobj.get_desc2nts(**kws)
    objwr = WrXlsxSortedGos(name, sortobj)
    # kws WrXlsxSortedGos wr_xlsx_nts: title hdrs
    objwr.wr_xlsx_nts(fout_xlsx, desc2nts, **kws)

def _get_grprobj():
    """Get object for grouping GO IDs."""
    fin_obo = os.path.join(REPO, "go-basic.obo")
    godag = get_godag(fin_obo, prt=None, loading_bar=False, optional_attrs=['relationship'])
    gosubdag = GoSubDag(USER_GOS, godag, relationships=True, tcntobj=None)
    grprdflt = GrouperDflts(gosubdag)
    hdrobj = HdrgosSections(gosubdag, grprdflt.hdrgos_dflt, SECTIONS)
    return Grouper("wrusrgos", USER_GOS, hdrobj, gosubdag)


if __name__ == '__main__':
    test_sort()
