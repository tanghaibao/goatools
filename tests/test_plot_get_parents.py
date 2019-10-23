#!/usr/bin/env python
"""Test the function, get_parents, in GoSubDag."""

from goatools.base import get_godag
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.gosubdag.plot.plot import plt_goids

def test_go_parents():
    """Run GO parent tests"""
    gosubdag_all = GoSubDag(None, get_godag("go-basic.obo", prt=None), rcntobj=True)
    run_1(gosubdag_all)
    run_2(gosubdag_all)

def run_1(gosubdag_all):
    """Test that using an Alt ID to color a parent in the middle does color."""
    # Main ID: GO:0035556 135 L04 D04 AB   intracellular signal transduction
    #  Alt ID: GO:0007242 135 L04 D04 AB   intracellular signal transduction
    alt_goid_middle = 'GO:0007242' # Alt GO ID. Key GO (GO:0035556) NOT in go2color
    alt_goid_color = '#f6cefc' # very light purple (Alt GO ID)
    goid_key = gosubdag_all.go2obj[alt_goid_middle].id

    #                  dcnt lev dep D1   Description
    goids = [ #        ---- --- --- ---- ----------------------------------------------
        "GO:0007165", # 716 L03 D03 AB   signal transduction (Header)
        "GO:0007166", # 336 L04 D04 AB   cell surface receptor signaling pathway
        "GO:0007186", #  99 L04 D04 AB   G-protein coupled receptor signaling pathway
        "GO:0097527", #   0 L04 D04 AB   necroptotic signaling pathway
        "GO:0007167", # 120 L05 D05 AB   enzyme linked receptor protein signaling pathway
        "GO:0042770", #  11 L05 D05 ABCG signal transduction in response to DNA damage
        "GO:0007229", #   0 L05 D05 AB   integrin-mediated signaling pathway
        "GO:0007205", #   0 L05 D05 AB   protein kinase C-activating GPCR signaling pathway
        "GO:0008630", #   1 L05 D06 ABCG intrinsic apoptotic signaling pw in rsp to DNA damage
        "GO:0070059", #   1 L05 D06 ABCG intrinsic apoptotic signaling pw in rsp to ER stress
        "GO:0035590", #   1 L06 D06 AB   purinergic nucleotide receptor signaling pathway
        "GO:0038063", #   0 L06 D07 AB   collagen-activated tyrosine kinase receptor signaling pw
    ]

    # If Alt ID is colored, then all equivalent GO IDs should be colored the same (unless overrode)
    go2color = {go:'#d6fffa' for go in goids} # klash ice
    go2color[alt_goid_middle] = alt_goid_color # very light purple (Alt GO ID)
    # Check that middle parent was NOT colored by user, even tho alt GO ID's color was set
    assert goid_key not in go2color
    # Plot
    godagplot = plt_goids(gosubdag_all, "test_get_parents1.png", goids, go2color=go2color)
    # Check that middle parent is colored properly, even if alt GO ID was used to set color
    assert godagplot.pydotnodego.go2color[alt_goid_middle] == alt_goid_color
    assert godagplot.pydotnodego.go2color[goid_key] == alt_goid_color
    # Check that original user data is NOT modified (User not expecting their data modified)
    assert goid_key not in go2color


def run_2(gosubdag):
    """Test GO colors at high and low levels of hierarchy."""
    goids = [
        'GO:0002682',  # GO:0002682 1,127 D03 A regulation of immune system process
        'GO:0002726']  # GO:0002726     2 D09 A +reg of T cell cytokine production
    gosubdag.prt_goids(goids)
    go2color = {
        'GO:0002682': '#b1fc99',  # pale light green
        'GO:0002726': '#f6cefc'}  # very light purple
    plt_goids(gosubdag, "test_get_parents2.png", goids, go2color=go2color, mark_alt_id=True)
    assert 'GO:0002682' in gosubdag.rcntobj.go2ancestors



if __name__ == '__main__':
    test_go_parents()
