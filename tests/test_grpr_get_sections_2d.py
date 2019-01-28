#!/usr/bin/env python
"""Test function, get_sections_2d, in the Grouper class.

The method, get_sections_2d should, return hdrgos used for
the user GO IDs in the same order as the 2-d sections list
provided when a Grouper object is initialized.

"""

import os
import sys
from goatools.base import get_godag
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.test_data.gjoneska_goea_consistent_increase import goea_results
from goatools.test_data.sections.gjoneska_pfenning import SECTIONS
from goatools.grouper.grprdflts import GrouperDflts
from goatools.grouper.hdrgos import HdrgosSections
from goatools.grouper.grprobj import Grouper

PRT = sys.stdout
REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def test_fnc():
    """Test function, get_sections_2d, in the Grouper class."""
    usrgo2nt = {getattr(nt, 'GO'):nt for nt in goea_results if getattr(nt, 'p_fdr_bh') < 0.05}
    usrgos = usrgo2nt.keys()
    grprdflt = _get_grprdflt()
    hdrobj = HdrgosSections(grprdflt.gosubdag, grprdflt.hdrgos_dflt, sections=SECTIONS, hdrgos=None)
    grprobj = Grouper("test", usrgos, hdrobj, grprdflt.gosubdag, go2nt=usrgo2nt)
    assert set(usrgos) == grprobj.usrgos
    sections_act = grprobj.get_sections_2d()
    chk_results(sections_act, grprobj)

def chk_results(sections_act, grprobj):
    """Get expected results."""
    hdrgos_act = grprobj.get_hdrgos()
    hdrgos_sec_act = set([g for _, gs in sections_act for g in gs])
    assert hdrgos_act == hdrgos_sec_act
    num_gos_orig = sum([len(gs) for _, gs in SECTIONS])
    PRT.write("{N} of {M} Sections Header GOs for {U} user GOs, {H} headers\n".format(
        N=len(hdrgos_sec_act), M=num_gos_orig, U=len(grprobj.usrgos), H=len(hdrgos_act)))
    sections_act_dict = {s:hs for s, hs in sections_act}
    # Check that order of actual header GOs is the same as found in the sections 2-d list
    for section_name, hdrgos_all in SECTIONS:
        hdrgos_act = sections_act_dict.get(section_name, None)
        if hdrgos_act is not None:
            h2i = {h:i for i, h in enumerate(hdrgos_act)}
            idx_act = None
            for hdrgo in hdrgos_all:
                idx = h2i.get(hdrgo, "")
                PRT.write("ORIG {I:2} {S} {H}\n".format(S=section_name, H=hdrgo, I=idx))
                if idx != "":
                    assert idx == 0 and idx_act is None or idx_act == idx - 1
                    idx_act = idx
            idx_act = None
            for hdrgo in hdrgos_act:
                idx = h2i.get(hdrgo, "")
                PRT.write("ACT  {I:2} {S} {H}\n".format(S=section_name, H=hdrgo, I=idx))
                assert idx == 0 and idx_act is None or idx_act == idx - 1
                idx_act = idx

def _get_gosubdag():
    """Get GO DAG."""
    fin = os.path.join(REPO, 'go-basic.obo')
    godag = get_godag(fin, prt=sys.stdout, loading_bar=False, optional_attrs=['relationship'])
    return GoSubDag(None, godag)

def _get_grprdflt():
    """Get Grouper defaults."""
    gosubdag = _get_gosubdag()
    fin_slim = os.path.join(REPO, 'goslim_generic.obo')
    return GrouperDflts(gosubdag, fin_slim)

if __name__ == '__main__':
    test_fnc()
