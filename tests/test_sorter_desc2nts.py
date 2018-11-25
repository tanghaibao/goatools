#!/usr/bin/env python
"""Test various options while sorting when using sections."""

from goatools.base import get_godag
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.grouper.grprdflts import GrouperDflts
from goatools.grouper.hdrgos import HdrgosSections
from goatools.grouper.grprobj import Grouper
from goatools.grouper.sorter import Sorter
from goatools.test_data.sorter import USER_GOS
from goatools.test_data.sorter import SECTIONS

FLDS = ('format_txt', 'hdr_idx', 'is_hdrgo', 'is_usrgo', 'num_usrgos', 'hdr1usr01',
        'NS', 'level', 'depth', 'reldepth', 'GO', 'alt', 'GO_name',
        'dcnt', 'D1', 'childcnt', 'REL', 'REL_short', 'rel', 'id')

# pylint: disable=line-too-long
def test_desc2nts():
    """Test various options while sorting when using sections."""
    sortobj = _get_sortobj()

    desc2nts = sortobj.get_desc2nts_fnc(hdrgo_prt=True, section_prt=None, top_n=None, use_sections=True)
    assert desc2nts.keys() == ['sections', 'hdrgo_prt', 'sortobj', 'flds']
    assert desc2nts['hdrgo_prt']
    h1a = set(nt.GO for sec, nts in desc2nts['sections'] for nt in nts)
    assert desc2nts['flds'] == FLDS
    assert desc2nts['sections'] == sortobj.get_desc2nts_fnc()['sections'], "**FAILED: Default values"

    desc2nts = sortobj.get_desc2nts_fnc(hdrgo_prt=True, section_prt=True, top_n=None, use_sections=True)
    assert desc2nts.keys() == ['sections', 'hdrgo_prt', 'sortobj', 'flds']
    assert desc2nts['hdrgo_prt']
    assert set(nt.GO for sec, nts in desc2nts['sections'] for nt in nts) == h1a
    assert desc2nts['flds'] == FLDS

    desc2nts = sortobj.get_desc2nts_fnc(hdrgo_prt=True, section_prt=False, top_n=None, use_sections=True)
    assert desc2nts.keys() == ['flat', 'hdrgo_prt', 'sortobj', 'flds']
    assert desc2nts['hdrgo_prt']
    assert set(nt.GO for nt in desc2nts['flat']) == h1a
    assert set(desc2nts['flds']).difference(set(FLDS)) == set(['section'])


    desc2nts = sortobj.get_desc2nts_fnc(hdrgo_prt=False, section_prt=None, top_n=None, use_sections=True)
    assert desc2nts.keys() == ['sections', 'hdrgo_prt', 'sortobj', 'flds']
    assert not desc2nts['hdrgo_prt']
    h0a = set(nt.GO for sec, nts in desc2nts['sections'] for nt in nts)
    assert len(h1a) > len(h0a), "**FAILED: MISSING HEADER GOs"
    assert desc2nts['flds'] == FLDS

    desc2nts = sortobj.get_desc2nts_fnc(hdrgo_prt=False, section_prt=True, top_n=None, use_sections=True)
    assert desc2nts.keys() == ['sections', 'hdrgo_prt', 'sortobj', 'flds']
    assert not desc2nts['hdrgo_prt']
    assert desc2nts['flds'] == FLDS

    desc2nts = sortobj.get_desc2nts_fnc(hdrgo_prt=False, section_prt=False, top_n=None, use_sections=True)
    assert desc2nts.keys() == ['flat', 'hdrgo_prt', 'sortobj', 'flds']
    assert not desc2nts['hdrgo_prt']
    assert set(nt.GO for nt in desc2nts['flat']) == h0a
    assert set(desc2nts['flds']).difference(set(FLDS)) == set(['section'])


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


def _chk_flat(desc2nts, h0a):
    """Check flat fields."""
    assert desc2nts.keys() == ['flat', 'hdrgo_prt', 'sortobj', 'flds']
    assert set(nt.GO for nt in desc2nts['flat']) == h0a
    assert desc2nts['flds'] == FLDS

def _get_sortobj():
    """Get object for grouping GO IDs."""
    godag = get_godag("go-basic.obo", prt=None, loading_bar=False, optional_attrs=['relationship'])
    gosubdag = GoSubDag(USER_GOS, godag, relationships=True, tcntobj=None)
    grprdflt = GrouperDflts(gosubdag)
    hdrobj = HdrgosSections(gosubdag, grprdflt.hdrgos_dflt, SECTIONS)
    grprobj = Grouper("wrusrgos", USER_GOS, hdrobj, gosubdag)
    return Sorter(grprobj)


if __name__ == '__main__':
    test_desc2nts()
