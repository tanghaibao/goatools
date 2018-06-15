#!/usr/bin/env python
"""Test writing an Excel spreadsheet given a set of GO IDs and go2obj."""

__copyright__ = "Copyright (C) 2016-2017, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import os
from numpy.random import shuffle
from goatools.obo_parser import GODag
from goatools.gosubdag.rpt.wr_xlsx import GoSubDagWr
# from goatools.test_data.gjoneska_2015_sections import SECTIONS


def main(num_goids=100):
    """Test writing an Excel spreadsheet given a set of GO IDs and go2obj."""
    objdat = Data("go-basic.obo")
    goids = objdat.get_goids(num_goids)
    objwr = GoSubDagWr(objdat.go2obj)
    fout_xlsx_flat = "a_goids_{N:05}_flat.xlsx".format(N=num_goids)
    objwr.wr_xlsx(fout_xlsx_flat, goids)
    #fout_xlsx_sect = "a_goids_{N:05}_sections.xlsx".format(N=num_goids)
    #objwr.wr_xlsx_sections(fout_xlsx_sect, goids, SECTIONS)


#pylint: disable=too-few-public-methods
class Data(object):
    """Holds data used in test."""

    def __init__(self, fin_obo):
        self.repo = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../")
        self.fin_obo = os.path.join(self.repo, fin_obo)
        self.dag = GODag(self.fin_obo)
        self.go2obj = {go:o for go, o in self.dag.items() if not o.is_obsolete}
        self.goids_all = self.go2obj.keys()

    def get_goids(self, num):
        """Return N randomly chosen GO IDs."""
        shuffle(self.goids_all)
        return set(self.goids_all[:num])


if __name__ == '__main__':
    main()

# Copyright (C) 2016-2017, DV Klopfenstein, H Tang, All rights reserved.
