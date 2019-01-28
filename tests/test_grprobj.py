#!/usr/bin/env python
"""Test method, grouper, in class, CountRelatives."""

import os
import sys
from goatools.base import get_godag
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.grouper.grprdflts import GrouperDflts
from goatools.grouper.hdrgos import HdrgosSections
from goatools.grouper.grprobj import Grouper
from goatools.grouper.wr_sections import WrSectionsTxt

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def test_grouper_d2(do_plot=False):
    """Group depth-02 GO terms under their most specific depth-01 GO parent(s)."""
    print('CWD', os.getcwd())
    # Get GOs to be grouped
    # Since no "Grouping GOs" were provided, depth-01 GOs are used for grouping.
    grprdflt = _get_grprdflt()
    hdrobj = HdrgosSections(grprdflt.gosubdag, grprdflt.hdrgos_dflt, sections=None, hdrgos=None)
    grprobj = Grouper("Transient Increase", get_data0(), hdrobj, _get_gosubdag(), go2nt=None)
    objwr = WrSectionsTxt(grprobj)
    objwr.wr_txt_section_hdrgos("transient_increase_hdrgos.txt")
    objwr.wr_txt_grouping_gos()
    if do_plot:
        # Don't run in Travis-CI because it does not contain 'dot'
        from goatools.grouper.grprplt import GrouperPlot
        GrouperPlot(grprobj).plot_groups_unplaced()
    chk_hdrs(grprobj)

def chk_hdrs(grprobj, prt=sys.stdout):
    """Check GO group headers."""
    hdrgos_all = grprobj.get_hdrgos()
    hdrgos_u0 = grprobj.get_hdrgos_u0()
    hdrgos_u1 = grprobj.get_hdrgos_u1()
    prt.write("{N} hdrgos ({U} are also user GO IDs) used out of {M} available\n".format(
        N=len(hdrgos_all), U=len(hdrgos_u1), M=len(grprobj.hdrobj.hdrgos)))
    assert hdrgos_u0.union(hdrgos_u1) == hdrgos_all

def get_data0():
    """Nature GO ids."""
    return [
        "GO:0007049", "GO:0022402", "GO:0022403", "GO:0000279", "GO:0006259",
        "GO:0000278", "GO:0051301", "GO:0000087", "GO:0007067", "GO:0000280",
        "GO:0048285", "GO:0006996", "GO:0006260", "GO:0006974", "GO:0033554",
        "GO:0006281", "GO:0016043", "GO:0051716", "GO:0009987", "GO:0006323",
        "GO:0051276", "GO:0007059", "GO:0006139", "GO:0065004", "GO:0051726",
        "GO:0007017", "GO:0031497", "GO:0034728", "GO:0006950", "GO:0034641",
        "GO:0006334", "GO:0006333", "GO:0034621", "GO:0006807", "GO:0006261",
        "GO:0007126", "GO:0051327", "GO:0051321", "GO:0034622", "GO:0044260",
        "GO:0043933", "GO:0000226", "GO:0042770", "GO:0000075", "GO:0006270",
        "GO:0065003", "GO:0006310", "GO:0010564", "GO:0022607", "GO:0006325",
        "GO:0043170", "GO:0008629", "GO:0007346", "GO:0044085", "GO:0008630",
        "GO:0051052", "GO:0050896", "GO:0031570", "GO:0051053", "GO:0007018",
        "GO:0007051", "GO:0007093", "GO:0006275", "GO:0009411", "GO:0034645",
        "GO:0000910", "GO:0009059", "GO:0044237", "GO:0010212", "GO:0000077",
        "GO:0030261", "GO:0009615", "GO:0002376", "GO:0006955", "GO:0006952",
        "GO:0002682", "GO:0050896", "GO:0048518", "GO:0048002", "GO:0050776",
        "GO:0009605", "GO:0048583", "GO:0050778", "GO:0051707", "GO:0048584",
        "GO:0019882", "GO:0009607", "GO:0019884", "GO:0002684", "GO:0048522",
        "GO:0009611", "GO:0051704", "GO:0002252", "GO:0001819", "GO:0006954",
        "GO:0002478", "GO:0002474", "GO:0030029", "GO:0006950", "GO:0030036",
        "GO:0007155", "GO:0022610", "GO:0048856", "GO:0050793", "GO:0048731",
        "GO:0048518", "GO:0006629", "GO:0007275", "GO:0032502", "GO:0051239",
        "GO:0048513", "GO:0048522", "GO:0044255", "GO:0048869", "GO:0009611",
        "GO:0031589", "GO:0009605", "GO:0010646", "GO:0002376", "GO:0043436",
        "GO:0019752", "GO:0006082", "GO:0032787", "GO:0042127", "GO:0009987",
        "GO:0030154", "GO:0042180", "GO:0001944", "GO:0065008", "GO:0006631",
        "GO:0009966", "GO:0048583", "GO:0002682", "GO:0001568", "GO:0009653",
        "GO:0007399", "GO:0007160", "GO:0045321", "GO:0001775", "GO:0080134",
        "GO:0051093", "GO:0048519", "GO:0030155", "GO:0007167", "GO:0042221",
        "GO:0045595", "GO:0048514", "GO:0042060", "GO:0030029", "GO:0048523",
        "GO:0002684", "GO:0051234", "GO:0006810", "GO:0051179", "GO:0007268",
        "GO:0019226", "GO:0051649", "GO:0051641", "GO:0009987", "GO:0007267",
        "GO:0006811", "GO:0007154", "GO:0007611", "GO:0015672", "GO:0006812",
        "GO:0006813", "GO:0046907", "GO:0006793", "GO:0006796", "GO:0030001",
        "GO:0006091", "GO:0007612", "GO:0015031", "GO:0043632", "GO:0019941",
        "GO:0045184", "GO:0045333", "GO:0031175", "GO:0048812", "GO:0044057",
        "GO:0048858", "GO:0007399", "GO:0007409", "GO:0030030", "GO:0010646",
        "GO:0032990", "GO:0006816", "GO:0048666", "GO:0051179", "GO:0048667",
        "GO:0048699", "GO:0030001", "GO:0006811", "GO:0022008", "GO:0050804",
        "GO:0009987", "GO:0006812", "GO:0000904", "GO:0031644", "GO:0006796",
        "GO:0006793", "GO:0051969", "GO:0030182", "GO:0016310", "GO:0015674",
        "GO:0007242", "GO:0006468", "GO:0006810", "GO:0051234", "GO:0007268",
        "GO:0000902", "GO:0019226", "GO:0051056", "GO:0043687", "GO:0032989",
        "GO:0006464", "GO:0016192", "GO:0016043", "GO:0007411", "GO:0043412",
        "GO:0007610", "GO:0007267", "GO:0009966", "GO:0048468", "GO:0007154",
        "GO:0048731", "GO:0006928", "GO:0015672"]

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
    test_grouper_d2(True)
