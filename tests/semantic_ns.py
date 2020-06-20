#!/usr/bin/env python
"""Full set of annotations can be used to set TermCounts. No need to break it up."""
# https://github.com/tanghaibao/goatools/issues/88

from __future__ import print_function

import os
import sys
from goatools.semantic import TermCounts
from goatools.base import get_godag
from goatools.associations import read_annotations

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def test_semantic_i88():
    """Full set of annotations can be used to set TermCounts. No need to break it up."""
    godag = get_godag("go-basic.obo")
    # Associations
    fin_gaf = os.path.join(REPO, "tair.gaf")
    gene2gos_all = read_annotations(gaf=fin_gaf, namespace='ALL')
    gene2gos_bp = read_annotations(gaf=fin_gaf, namespace='BP')
    gene2gos_mf = read_annotations(gaf=fin_gaf, namespace='MF')
    gene2gos_cc = read_annotations(gaf=fin_gaf, namespace='CC')
    # Termcounts
    prt = sys.stdout
    termcounts_all = TermCounts(godag, gene2gos_all, prt=prt)
    termcounts_bp = TermCounts(godag, gene2gos_bp, prt=prt)
    termcounts_mf = TermCounts(godag, gene2gos_mf, prt=prt)
    termcounts_cc = TermCounts(godag, gene2gos_cc, prt=prt)
    # Test content in subset is the same as in the full GO counts
    for goid, cnt in termcounts_bp.gocnts.items():
        assert termcounts_all.gocnts[goid] == cnt
    for goid, cnt in termcounts_mf.gocnts.items():
        assert termcounts_all.gocnts[goid] == cnt
    for goid, cnt in termcounts_cc.gocnts.items():
        assert termcounts_all.gocnts[goid] == cnt


if __name__ == '__main__':
    test_semantic_i88()
