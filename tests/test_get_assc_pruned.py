#!/usr/bin/env python
"""Test the function 'get_assc_pruned'."""

__copyright__ = "Copyright (C) 2010-2018, H Tang et al. All rights reserved."

import os
import sys
from goatools.associations import dnld_assc
from goatools.associations import get_assc_pruned
from goatools.utils import get_b2aset

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


# pylint: disable=too-few-public-methods
class Run(object):
    """Holds the data for the tests. Runs the tests. Checks the results."""

    def __init__(self):
        self.prt = sys.stdout
        _fin_assc = os.path.join(REPO, "goa_human.gaf")
        self.gene2gos_orig = dnld_assc(_fin_assc, go2obj=None, prt=self.prt)
        self.go2genes_orig = get_b2aset(self.gene2gos_orig)
        _num_genes = [len(gs) for gs in self.go2genes_orig.values()]
        self.min_genes = min(_num_genes)
        self.max_genes = max(_num_genes)
        assert self.gene2gos_orig == get_b2aset(self.go2genes_orig)

    def run(self, min_genecnt, max_genecnt, expchg):
        """Prune the association."""
        # Prune the association
        g2gos, rm_gos = get_assc_pruned(self.gene2gos_orig, min_genecnt, max_genecnt, self.prt)
        # Expect the assocation to be changed
        if expchg:
            self._chk_pruned(g2gos, rm_gos, min_genecnt, max_genecnt)
        # Expect no change in the association
        else:
            self._chk_no_change(g2gos, rm_gos)

    def _chk_pruned(self, g2gos, rm_gos, min_genecnt, max_genecnt):
        """Check that the association was pruned correctly."""
        assert self.gene2gos_orig != g2gos, "MIN({}) MAX({})".format(min_genecnt, max_genecnt)
        assert rm_gos
        if min_genecnt is None:
            min_genecnt = 1
        if max_genecnt is None:
            max_genecnt = self.max_genes
        for goid in rm_gos:
            num_genes = len(self.go2genes_orig[goid])
            assert num_genes <= min_genecnt or num_genes >= max_genecnt

    def _chk_no_change(self, g2gos, rm_gos):
        """Check that the 'pruned' assc. has not been changed."""
        assert self.gene2gos_orig == g2gos
        assert not rm_gos



def test_fnc():
    """Test the function 'get_assc_pruned'."""
    obj = Run()
    # pylint: disable=bad-whitespace
    #       min_genecnt           max_genecnt      Assc. expected to change
    #       -----------           -----------      ------------------------
    obj.run(None,                 None,              False)
    obj.run(0,                    None,              False)
    obj.run(1,                    None,              False)
    obj.run(obj.min_genes,        None,              False)
    obj.run(None,                 obj.max_genes,     False)
    obj.run(obj.min_genes+1, None,                    True)
    obj.run(None,                 obj.max_genes-1,    True)
    obj.run(2,                    None,               True)
    obj.run(3,                    None,               True)
    obj.run(4,                    None,               True)

if __name__ == '__main__':
    test_fnc()

# Copyright (C) 2010-2018, H Tang et al. All rights reserved."
