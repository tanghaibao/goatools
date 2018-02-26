#!/usr/bin/env python
"""Test get_all_parents vs """

from __future__ import print_function

__copyright__ = "Copyright (C) 2010-2018, H Tang et al. All rights reserved."

import os
import sys
import timeit
from goatools.base import get_godag
from goatools.godag.go_tasks import get_go2parents
from goatools.test_data.godag_timed import prt_hms


def test_get_parent(prt=sys.stdout):
    """Semantic Similarity test for Issue #86."""
    # Load GO-DAG
    fin_obo = "go-basic.obo"
    repo = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    godag = get_godag(os.path.join(repo, fin_obo))
    go2obj = {go:o for go, o in godag.items() if go == o.id}
    # Get all parents for all GO IDs using GOTerm get_all_parents
    tic = timeit.default_timer()
    go2parents_orig = {}
    for goobj in go2obj.values():
        go2parents_orig[goobj.id] = goobj.get_all_parents()
    tic = prt_hms(prt, tic, "Get all goobj's parents using GOTerm.get_all_parents()")
    # Get all parents for all GO IDs using GOTerm get_all_parents
    go2parents_task = get_go2parents(go2obj.values())
    prt_hms(prt, tic, "Get all goobj's parents using go_tasks::get_go2parents")
    # Compare parent lists
    assert set(go2parents_orig) == set(go2parents_task)
    for goid, parents_exp in go2parents_orig.items():
        parents_act = go2parents_task[goid]
        assert parents_act == parents_exp
    print("PASSED: get_parent test")


if __name__ == '__main__':
    PRT = None if len(sys.argv) != 1 else sys.stdout
    test_get_parent(PRT)

# Copyright (C) 2010-2018, H Tang et al. All rights reserved.
