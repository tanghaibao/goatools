#!/usr/bin/env python
"""Test get_all_children vs """

from __future__ import print_function

__copyright__ = "Copyright (C) 2010-2018, H Tang et al. All rights reserved."

import os
import sys
import timeit
from goatools.base import get_godag
from goatools.godag.go_tasks import get_go2children
from goatools.test_data.godag_timed import prt_hms
from goatools.test_data.checks import _chk_a2bset


def test_get_children(prt=sys.stdout):
    """Semantic Similarity test for Issue #86."""
    # Load GO-DAG
    fin_obo = "go-basic.obo"
    repo = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    godag = get_godag(os.path.join(repo, fin_obo))
    go2obj = {go:o for go, o in godag.items() if go == o.id}
    # Get all children for all GO IDs using get_all_children in GOTerm class
    tic = timeit.default_timer()
    go2children_orig = {}
    for goobj in go2obj.values():
        go2children_orig[goobj.id] = goobj.get_all_children()
    tic = prt_hms(prt, tic, "Get all goobj's children using GOTerm.get_all_children()")
    # Get all children for all GO IDs using GOTerm get_all_children
    go2children_fast = get_go2children(go2obj.values())
    prt_hms(prt, tic, "Get all goobj's children using go_tasks::get_go3children")
    # Compare children lists
    _chk_a2bset(go2children_orig, go2children_fast)


if __name__ == '__main__':
    PRT = None if len(sys.argv) != 1 else sys.stdout
    test_get_children(PRT)

# Copyright (C) 2010-2018, H Tang et al. All rights reserved.
