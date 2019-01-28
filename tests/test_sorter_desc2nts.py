#!/usr/bin/env python
"""Test various options while sorting when using sections."""

import os
from goatools.base import get_godag
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.grouper.grprdflts import GrouperDflts
from goatools.grouper.hdrgos import HdrgosSections
from goatools.grouper.grprobj import Grouper
from goatools.grouper.sorter import Sorter
from goatools.test_data.sorter import USER_GOS
from goatools.test_data.sorter import SECTIONS

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

GO_FLDS = ('format_txt', 'hdr_idx', 'is_hdrgo', 'is_usrgo', 'num_usrgos', 'hdr1usr01',
           'NS', 'level', 'depth', 'reldepth', 'GO', 'alt', 'GO_name',
           'dcnt', 'D1', 'childcnt', 'REL', 'REL_short', 'rel', 'id')

D2_FLDS = set(['flds', 'num_items', 'num_sections', 'hdrgo_prt', 'sortobj', 'sections'])
D1_FLDS = set(['flds', 'num_items', 'num_sections', 'hdrgo_prt', 'sortobj', 'flat'])


# pylint: disable=line-too-long
def test_desc2nts():
    """Test various options while sorting when using sections."""
    sortobj = _get_sortobj()

    desc2nts = sortobj.get_desc2nts_fnc(hdrgo_prt=True, section_prt=None, top_n=None, use_sections=True)
    assert desc2nts['hdrgo_prt']
    h1a = set(nt.GO for sec, nts in desc2nts['sections'] for nt in nts)
    num_sections = len(desc2nts['sections'])
    assert desc2nts['sections'] == sortobj.get_desc2nts_fnc()['sections'], "**FAILED: Default values"
    _chk_d2(desc2nts, num_sections)

    desc2nts = sortobj.get_desc2nts_fnc(hdrgo_prt=True, section_prt=True, top_n=None, use_sections=True)
    assert desc2nts['hdrgo_prt']
    assert set(nt.GO for sec, nts in desc2nts['sections'] for nt in nts) == h1a
    _chk_d2(desc2nts, num_sections)

    desc2nts = sortobj.get_desc2nts_fnc(hdrgo_prt=True, section_prt=False, top_n=None, use_sections=True)
    assert set(desc2nts.keys()) == D1_FLDS
    assert desc2nts['hdrgo_prt']
    assert set(nt.GO for nt in desc2nts['flat']) == h1a
    assert set(desc2nts['flds']).difference(set(GO_FLDS)) == set(['section'])
    assert desc2nts['num_sections'] == num_sections


    desc2nts = sortobj.get_desc2nts_fnc(hdrgo_prt=False, section_prt=None, top_n=None, use_sections=True)
    h0a = set(nt.GO for sec, nts in desc2nts['sections'] for nt in nts)
    assert not desc2nts['hdrgo_prt']
    assert len(h1a) > len(h0a), "**FAILED: MISSING HEADER GOs"
    _chk_d2(desc2nts, num_sections)

    desc2nts = sortobj.get_desc2nts_fnc(hdrgo_prt=False, section_prt=True, top_n=None, use_sections=True)
    assert not desc2nts['hdrgo_prt']
    _chk_d2(desc2nts, num_sections)

    desc2nts = sortobj.get_desc2nts_fnc(hdrgo_prt=False, section_prt=False, top_n=None, use_sections=True)
    assert set(desc2nts.keys()) == D1_FLDS
    assert not desc2nts['hdrgo_prt']
    assert set(nt.GO for nt in desc2nts['flat']) == h0a
    assert set(desc2nts['flds']).difference(set(GO_FLDS)) == set(['section'])
    assert desc2nts['num_sections'] == num_sections


    desc2nts = sortobj.get_desc2nts_fnc(hdrgo_prt=True, section_prt=None, top_n=None, use_sections=False)
    assert desc2nts['hdrgo_prt']
    _chk_flat(desc2nts, h1a)

    desc2nts = sortobj.get_desc2nts_fnc(hdrgo_prt=True, section_prt=True, top_n=None, use_sections=False)
    assert desc2nts['hdrgo_prt']
    _chk_flat(desc2nts, h1a)

    desc2nts = sortobj.get_desc2nts_fnc(hdrgo_prt=True, section_prt=False, top_n=None, use_sections=False)
    assert desc2nts['hdrgo_prt']
    _chk_flat(desc2nts, h1a)


    desc2nts = sortobj.get_desc2nts_fnc(hdrgo_prt=False, section_prt=None, top_n=None, use_sections=False)
    assert not desc2nts['hdrgo_prt']
    _chk_flat(desc2nts, h0a)

    desc2nts = sortobj.get_desc2nts_fnc(hdrgo_prt=False, section_prt=True, top_n=None, use_sections=False)
    assert not desc2nts['hdrgo_prt']
    _chk_flat(desc2nts, h0a)

    desc2nts = sortobj.get_desc2nts_fnc(hdrgo_prt=False, section_prt=False, top_n=None, use_sections=False)
    assert not desc2nts['hdrgo_prt']
    _chk_flat(desc2nts, h0a)


def _chk_d2(desc2nts, num_sections):
    """Check section fields."""
    assert set(desc2nts.keys()) == D2_FLDS
    assert desc2nts['flds'] == GO_FLDS
    assert desc2nts['num_sections'] == num_sections

def _chk_flat(desc2nts, h0a):
    """Check flat fields."""
    assert set(desc2nts.keys()) == D1_FLDS
    assert set(nt.GO for nt in desc2nts['flat']) == h0a
    assert desc2nts['flds'] == GO_FLDS

def _get_sortobj():
    """Get object for grouping GO IDs."""
    fin_godag = os.path.join(REPO, "go-basic.obo")
    godag = get_godag(fin_godag, prt=None, loading_bar=False, optional_attrs=['relationship'])
    gosubdag = GoSubDag(USER_GOS, godag, relationships=True, tcntobj=None)
    grprdflt = GrouperDflts(gosubdag)
    hdrobj = HdrgosSections(gosubdag, grprdflt.hdrgos_dflt, SECTIONS)
    grprobj = Grouper("wrusrgos", USER_GOS, hdrobj, gosubdag)
    return Sorter(grprobj)


if __name__ == '__main__':
    test_desc2nts()
