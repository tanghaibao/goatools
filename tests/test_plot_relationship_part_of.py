#!/usr/bin/env python
"""Plot both the standard 'is_a' field and the optional 'part_of' relationship."""

from __future__ import print_function

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved."

import os
import sys
from goatools.test_data.wr_subobo import WrSubObo
from goatools.base import download_go_basic_obo
from goatools.obo_parser import GODag
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.gosubdag.plot.gosubdag_plot import GoSubDagPlot



REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../")

NAME2GOIDS = {
    'smell': set([
        "GO:0007608",   # sensory perception of smell
        "GO:0050911"]), # detection of chemical stimulus involved in sensory perception of smell

    'secretory':set([
        "GO:0030141",   # CC 20 L06 D07 AC  P... p... secretory granule
        "GO:0034774"]), # CC  8 L05 D06 ABD P... .... secretory granule lumen

    # MISSING Edge source(GO:0007507); target(GO:0072359)
    # MISSING:  GO:0061371  # BP 0 L06 D06 B P p determination of heart left/right asymmetry
    # ERROR:    GO:0007507
    'heartjogging':set([
        "GO:0003304",  # BP 0 L06 D06 A P . myocardial epithelial involution in heart jogging
        "GO:0003146"]) # BP 0 L05 D07 A P p heart jogging
}

def test_plot_part_of():
    """Plot both the standard 'is_a' field and the 'part_of' relationship."""
    fout_log = "plot_relationship_part_of.log"
    obj = _Run()
    names = NAME2GOIDS
    # names = ["heartjogging"]
    with open(fout_log, 'w') as prt:
        for name in names:
            goids = NAME2GOIDS[name]
            obj.plot_all(goids, name, prt)
        print("  WROTE: {LOG}\n".format(LOG=fout_log))


# pylint: disable=too-few-public-methods
class _Run(object):
    """Holds data and steps for this test."""

    obopat = "../goatools/data/{NAME}.obo"

    def __init__(self):
        _fin_obo = os.path.join(REPO, "go-basic.obo")
        self.go2obj = GODag(_fin_obo, optional_attrs=['relationship'])

    def plot_all(self, goids, name, prt=sys.stdout):
        """Create plots with various numbers of relationships."""
        prt.write("\nCreate GoSubDag not loading any relationship")
        gosubdag_orig = GoSubDag(goids, self.go2obj, relationships=False, prt=prt)
        gosubdag_orig.prt_goids(gosubdag_orig.go2obj, prt=prt)
        prt.write("{N} GO IDS".format(N=len(gosubdag_orig.go2obj)))
        gopltdag = GoSubDagPlot(gosubdag_orig, mark_alt_id=True)
        gopltdag.plt_dag(os.path.join(REPO, "a_relationship_{NAME}_r0.png".format(NAME=name)))

        # goids.update(['GO:0007507'], ['GO:0072359'])
        prt.write("\nCreate GoSubDag while loading only the 'part_of' relationship")
        gosubdag = GoSubDag(goids, self.go2obj, relationships=['part_of'], prt=prt)
        gosubdag.prt_goids(gosubdag.go2obj, prt=prt)
        prt.write("{N} GO IDS".format(N=len(gosubdag.go2obj)))
        gopltdag = GoSubDagPlot(gosubdag, mark_alt_id=True)
        prt.write("GO SOURCES:")
        gosubdag.prt_goids(gosubdag.go_sources, prt=prt)
        gopltdag.plt_dag(os.path.join(REPO, "a_relationship_{NAME}_partof.png".format(NAME=name)))

        prt.write("\nCreate GoSubDag while loading all relationships")
        gosubdag = GoSubDag(goids, self.go2obj, relationships=True, prt=prt)
        prt.write("ALL {N} GO IDS:".format(N=len(gosubdag.go2obj)))
        gosubdag.prt_goids(gosubdag.go2obj, prt=prt)
        prt.write("2 GO SOURCES:")
        gosubdag.prt_goids(gosubdag.go_sources, prt=prt)
        goids_new = set(gosubdag.go2obj).difference(set(gosubdag_orig.go2obj))
        go2color = {go:'#d5ffff' for go in goids_new}
        prt.write("{N} NEW GO IDS:".format(N=len(goids_new)))
        gosubdag.prt_goids(goids_new, prt=prt)
        prt.write("{N} GO IDS".format(N=len(gosubdag.go2obj)))
        gopltdag = GoSubDagPlot(gosubdag, mark_alt_id=True, go2color=go2color)
        gopltdag.plt_dag(os.path.join(REPO, "a_relationship_{NAME}_r1.png".format(NAME=name)))

    def wr_subobo(self):
        """Write a subset obo to be used for testing."""
        # Load GO-DAG: Load optional 'relationship'
        for name, goids in NAME2GOIDS.items():
            fout_obo = self.get_obo_name(name)
            fin_obo = os.path.join(REPO, "go-basic.obo")
            download_go_basic_obo(fin_obo, prt=sys.stdout, loading_bar=None)
            obj = WrSubObo(fin_obo, optional_attrs=['relationship'])
            # obj = WrSubObo(fin_obo)
            obj.wrobo(fout_obo, goids)

    def get_obo_name(self, name):
        """Get the name of the obo for a small subset."""
        return os.path.join(REPO, self.obopat.format(NAME=name))


if __name__ == '__main__':
    #_Run().wr_subobo()
    test_plot_part_of()

# Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved.
