#!/usr/bin/env python
"""Test method, sorter, in class, CountRelatives."""

from __future__ import print_function

import os
import sys

from goatools.base import get_godag
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.grouper.grprdflts import GrouperDflts
from goatools.grouper.hdrgos import HdrgosSections
from goatools.grouper.grprobj import Grouper
from goatools.grouper.sorter import Sorter

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


# pylint: disable=too-many-locals
def test_dflthdrs(prt=sys.stdout, do_plt=False):
    """Group depth-02 GO terms under their most specific depth-01 GO parent(s)."""
    # Initialize GoSubDag for grouping use once, since it takes a few seconds to initialize
    grprdflt = _get_grprdflt()

    # Get GOs to be grouped
    data = get_data0()

    # This may need to be updated if default hdrgos are changed
    exp_hdrs0 = set([
        "GO:0050789",  # BP 11,095 L01 D01 B     regulation of biological process
        "GO:0044848",  # BP     62 L01 D01 S     biological phase
        "GO:0050794",  # BP  8,031 L02 D02 AB    regulation of cellular process
        "GO:0019222",  # BP  3,227 L02 D02 AB    regulation of metabolic process
        "GO:0048583",  # BP  2,377 L02 D02 AB    regulation of response to stimulus
        "GO:0050793",  # BP  1,789 L02 D02 AB    regulation of developmental process
        "GO:0023051",  # BP  1,364 L02 D02 AB    regulation of signaling
        "GO:0002682",  # BP  1,183 L02 D02 AB    regulation of immune system process
        "GO:0007155",  # BP    165 L02 D02 P     cell adhesion
        "GO:0080134",  # BP    940 L03 D03 AB    regulation of response to stress
        "GO:0007165",  # BP    717 L03 D03 AB    signal transduction
        "GO:0050877",  # BP     96 L03 D03 K     neurological system process
        "GO:0007267"]) # BP     99 L03 D04 CDR   cell-cell signaling


    # Since no "GO group headers" (None) were provided, depth-01 GOs are used for grouping.
    hdrobj0 = HdrgosSections(grprdflt.gosubdag, grprdflt.hdrgos_dflt, sections=None, hdrgos=None)
    grprobj0 = Grouper("dflt", data, hdrobj0, grprdflt.gosubdag, go2nt=None)
    _, _, nts0_go, act_hdrs0 = run(grprobj0, hdrobj0, exp_hdrs0)

    # Grouping GOs are provided, these are added to the depth-01 defaults GOs are used for grouping.
    hdrgos = set([
        "GO:0099536", # BP     40 L04 D05 CDR   regulation of response to stimulus
        "GO:0051239", # BP  2,532 L02 D02 AB    regulation of multicellular organismal process
        "GO:0048519", # BP  3,293 L02 D02 AB    negative regulation of biological process
        "GO:0048518"])# BP  3,353 L02 D02 AB    positive regulation of biological process

    exp_hdrs1 = exp_hdrs0.union(hdrgos)
    name = "usrhdrs4"
    hdrobj1 = HdrgosSections(grprdflt.gosubdag, grprdflt.hdrgos_dflt, sections=None, hdrgos=hdrgos)
    grprobj1 = Grouper(name, data, hdrobj1, grprdflt.gosubdag, go2nt=None)
    sortobj1, _, nts1_go, act_hdrs1 = run(grprobj1, hdrobj1, exp_hdrs1)

    if do_plt:
        from goatools.grouper.grprplt import GrouperPlot
        prt.write("\nPLOT DAG\n")
        GrouperPlot(grprobj1).plot_grouped_gos()

    # GO:0099536 was a "middle" term (neither usrgo, not hdrgo) in run0, but is a hdrgo in run1

    # print "THIS"
    # grprdflt.gosubdag.prt_goids(nts1_go)
    # print "MINUS"
    # grprdflt.gosubdag.prt_goids(nts0_go)
    # print "EQUALS"
    # print nts1_go.difference(nts0_go)

    assert nts1_go.difference(nts0_go) == set(["GO:0099536"])
    assert act_hdrs1.difference(act_hdrs0) == set(hdrgos)

    hdrgo_prt = False
    sys.stdout.write("\n{NAME}: PRINT GOs hdrgo_prt({H}):\n".format(H=hdrgo_prt, NAME=name))
    sortobj1.prt_gos(hdrgo_prt=hdrgo_prt)
    nts2 = sortobj1.get_nts_flat(hdrgo_prt)
    nts2_go = set([nt.GO for nt in nts2])

    assert len(nts1_go) > len(nts2_go)
    assert nts1_go.intersection(data) == nts2_go
    assert nts2_go == data


def run(grprobj, hdrobj, exp_hdrs, hdrgo_prt=True):
    """Load sorter. Check results."""
    chk_hdrs(grprobj, hdrobj)
    act_hdrs = grprobj.get_hdrgos()

    print("ACTUAL")
    grprobj.gosubdag.prt_goids(sorted(act_hdrs))
    print("EXPECTED")
    grprobj.gosubdag.prt_goids(sorted(exp_hdrs))
    # assert act_hdrs == exp_hdrs

    sortobj = Sorter(grprobj, hdrgo_prt=hdrgo_prt)
    sys.stdout.write("\n{NAME} PRINT GOs hdrgo_prt({H}):\n".format(
        H=hdrgo_prt, NAME=grprobj.grpname))
    sortobj.prt_gos()
    nts = sortobj.get_nts_flat(hdrgo_prt)
    nts_go = set([nt.GO for nt in nts])
    usrgos = grprobj.usrgos
    assert nts_go.intersection(usrgos) == usrgos, \
        "ONLY {N} of {U} user gos found in grouped sorted GOs. MISSING: {GOs}".format(
            N=len(nts_go.intersection(usrgos)),
            GOs=" ".join(usrgos.difference(nts_go.intersection(usrgos))),
            U=len(usrgos))
    return sortobj, nts, nts_go, act_hdrs


def chk_hdrs(grprobj, hdrobj, prt=sys.stdout):
    """Check GO group headers."""
    hdrgos_all = grprobj.get_hdrgos()
    hdrgos_u0 = grprobj.get_hdrgos_u0()
    hdrgos_u1 = grprobj.get_hdrgos_u1()
    prt.write("{N} hdrgos ({U} are also user GO IDs) used out of {M} available\n".format(
        N=len(hdrgos_all), U=len(hdrgos_u1), M=len(hdrobj.hdrgos)))
    assert hdrgos_u0.union(hdrgos_u1) == hdrgos_all

def get_data0():
    """Nature GO ids."""
    return set([
        #"GO:0050789", # BP 1 11,101 L01 D01 B  reg. of biological process
        "GO:0051969", # BP       5 L03 D05 AB  reg. of transmission of nerve impulse
        "GO:0008629", # BP      13 L05 D05 AB  intrinsic apoptotic signaling pathway
        "GO:0051056", # BP      26 L05 D06 AB  reg. of small GTPase mediated signal transduction
        "GO:0031644", # BP      30 L04 D04 AB  reg. of neurological system process
        "GO:0006275", # BP      50 L05 D06 AB  reg. of DNA replication
        "GO:0051053", # BP   *  76 L05 D06 AB  negative reg. of DNA metabolic process
        "GO:0007167", # BP     121 L05 D05 AB  enzyme linked receptor protein signaling pathway
        "GO:0050804", # BP     120 L03 D04 AB  modulation of synaptic transmission
        "GO:0007242", # BP     135 L04 D04 AB  intracellular signal transduction
        "GO:0007346", # BP     157 L04 D04 AB  reg. of mitotic cell cycle
        "GO:0001819", # BP     154 L04 D04 AB  positive reg. of cytokine production
        "GO:0051052", # BP     225 L04 D05 AB  reg. of DNA metabolic process
        "GO:0050778", # BP     227 L04 D04 AB  positive reg. of immune response
        "GO:0030155", # BP     246 L02 D02 AB  reg. of cell adhesion
        "GO:0042127", # BP     268 L03 D03 AB  reg. of cell proliferation
        "GO:0010564", # BP     350 L04 D04 AB  reg. of cell cycle process
        "GO:0044057", # BP   * 392 L03 D03 AB  reg. of system process
        "GO:0051726", # BP     404 L03 D03 AB  reg. of cell cycle
        "GO:0002684", # BP   * 436 L03 D03 AB  positive reg. of immune system process
        "GO:0051093", # BP     549 L03 D03 AB  negative reg. of developmental process
        "GO:0050776", # BP     661 L03 D03 AB  reg. of immune response
        "GO:0048584", # BP     776 L03 D03 AB  positive reg. of response to stimulus
        "GO:0045595", # BP     828 L03 D03 AB  reg. of cell differentiation
        "GO:0080134", # BP     940 L03 D03 AB  reg. of response to stress
        "GO:0009966", # BP   1,108 L03 D04 AB  reg. of signal transduction
        "GO:0002682", # BP   1,183 L02 D02 AB  reg. of immune system process
        "GO:0010646", # BP   1,392 L03 D03 AB  reg. of cell communication
        "GO:0050793", # BP   1,789 L02 D02 AB  reg. of developmental process
        "GO:0048522", # BP   2,289 L03 D03 AB  positive reg. of cellular process
        "GO:0048523", # BP   2,372 L03 D03 AB  negative reg. of cellular process
        #"GO:0048583", # BP   2,377 L02 D02 AB  reg. of response to stimulus
        "GO:0051239", # BP   2,532 L02 D02 AB  reg. of multicellular organismal process
        "GO:0048519", # BP   3,293 L02 D02 AB  negative reg. of biological process
        "GO:0048518", # BP   3,353 L02 D02 AB  positive reg. of biological process
        #"GO:0044848", # BP 1    62 L01 D01 S   biological phase
        "GO:0000087", # BP 0     0 L04 D04 S   mitotic M phase
        "GO:0051327", # BP 0     0 L04 D04 S   meiotic M phase
        "GO:0000279", # BP 0     2 L03 D03 S   M phase
        "GO:0022403", # BP 0    46 L02 D02 S   cell cycle phase
        #"GO:0023052", # BP 1   116 L01 D01 R   signaling
        "GO:0019226", # BP 0     0 L04 D04 DKR transmission of nerve impulse
        "GO:0007268", # BP 0    12 L07 D08 CDR chemical synaptic transmission
        "GO:0007267", # BP 0    99 L03 D04 CDR cell-cell signaling
        #"GO:0022610", # BP 1   194 L01 D01 P   biological adhesion
        "GO:0007155", # BP 0   165 L02 D02 P   cell adhesion
        #"GO:0007610", # BP 1   219 L01 D01 O   behavior
        "GO:0007612", # BP 0    14 L04 D06 DKO learning
        "GO:0007611"])# BP 0    22 L03 D05 DKO learning or memory

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
    test_dflthdrs(do_plt=True)
