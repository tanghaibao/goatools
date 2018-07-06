#!/usr/bin/env python
"""Test reading GO IDs from a file."""

import os
import sys
from goatools.grouper.grprdflts import GrouperDflts
from goatools.grouper.hdrgos import HdrgosSections
from goatools.grouper.grprobj import Grouper
from goatools.grouper.read_goids import read_sections
from goatools.grouper.wr_sections import WrSections
from goatools.test_data.gjoneska_goea_consistent_increase import goea_results

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def test_wr_sections_txt():
    """Group depth-02 GO terms under their most specific depth-01 GO parent(s)."""
    # Get GOs to be grouped
    usrgos = [getattr(nt, 'GO') for nt in goea_results]
    # Read OBO files once to save time
    grprdflt = GrouperDflts()

    # ------------------------------------------------------------------
    # Print usrgos in txt (Do not use sections containing hdrgos)
    # ------------------------------------------------------------------
    # Show GO grouping hdrgos and usrgos to show how usrgos are grouped
    _wr_sections_txt("a0_hdr1.txt", usrgos, sections_file=None, grprdflt=grprdflt)

    # ------------------------------------------------------------------
    # Print usrgos in txt using sections containing hdrgos
    # ------------------------------------------------------------------
    sec1 = _read_sections("./data/gjoneska/sections_in.txt")
    # Print usrgos in sections, showing how they were grouped under hdrgos
    _wr_sections_txt("sec0_hdr1.txt", usrgos, sec1, grprdflt=grprdflt)


def _read_sections(fin):
    """Get sections variable from file."""
    sec = os.path.join(REPO, fin)
    assert read_sections(sec), "EMPTY SECTIONS FILE({})".format(sec)
    return sec

def _wr_sections_txt(fout_txt, usrgos, sections_file, grprdflt):
    """Given a list of usrgos and sections, write text file."""
    try:
        sections = read_sections(sections_file)
        hdrobj = HdrgosSections(grprdflt.gosubdag, grprdflt.hdrgos_dflt, sections=sections)
        grprobj = Grouper(fout_txt, usrgos, hdrobj, grprdflt.gosubdag, go2nt=None)
        full_txt = os.path.join(REPO, fout_txt)
        WrSections(grprobj).wr_txt_section_hdrgos(full_txt, sortby=None, prt_section=True)
        assert os.path.exists(full_txt)
    except RuntimeError as inst:
        sys.stdout.write("\n  **FATAL: {MSG}\n\n".format(MSG=str(inst)))


if __name__ == '__main__':
    test_wr_sections_txt()
