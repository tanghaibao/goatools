#!/usr/bin/env python
"""Test get_all_parents vs """

from __future__ import print_function

__copyright__ = "Copyright (C) 2010-2019, DV Klopfenstein, H Tang et al. All rights reserved."

import os
import sys
import timeit
from goatools.base import get_godag
from goatools.godag.go_tasks import get_id2parents
from goatools.godag.prttime import prt_hms
from goatools.test_data.checks import CheckGOs


def test_get_parent(prt=sys.stdout):
    """Semantic Similarity test for Issue #86."""
    # Load GO-DAG
    fin_obo = "go-basic.obo"
    repo = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    godag = get_godag(os.path.join(repo, fin_obo))
    go2obj = {go:o for go, o in godag.items() if go == o.id}
    # ------------------------------------------------------------------------
    # Get all parents for all GO IDs using get_all_parents in GOTerm class
    tic = timeit.default_timer()
    go2parents_orig = {}
    ## go_noparents = set()
    for goterm in go2obj.values():
        parents = goterm.get_all_parents()
        #if parents:
        go2parents_orig[goterm.id] = parents
        #else:
        #    go_noparents.add(goterm.id)
    tic = prt_hms(tic, "Get all goobj's parents using GOTerm.get_all_parents()", prt)
    # ------------------------------------------------------------------------
    # Get all parents for all GO IDs using GOTerm get_all_parents
    go2parents_fast = get_id2parents(go2obj.values())
    tic = prt_hms(tic, "Get all goobj's parents using go_tasks::get_id2parents", prt)
    # ------------------------------------------------------------------------
    go2parents_fast2 = get_id2parents2(go2obj.values())
    tic = prt_hms(tic, "Get all goobj's parents using go_tasks::get_id2parents2", prt)
    # ------------------------------------------------------------------------
    # Compare parent lists
    chkr = CheckGOs('test_get_parents', go2obj)
    chkr.chk_a2bset_equiv(go2parents_orig, go2parents_fast)
    chkr.chk_a2bset_equiv(go2parents_orig, go2parents_fast2)
    print("PASSED: get_parent test")


# ------------------------------------------------------------------------------------
def get_id2parents2(objs):
    """Get all parent item IDs for each item in dict keys."""
    ofnc = IdToParents()
    ofnc.set_id2parents(objs)
    return ofnc.id2parents

# pylint:disable=too-few-public-methods
class IdToParents:
    """Test functions"""

    def __init__(self):
        self.id2parents = {}

    def set_id2parents(self, objs):
        """Get all parent item IDs for each item in dict keys."""
        for obj in objs:
            self._get_id2parents(obj.item_id, obj)

    def _get_id2parents(self, item_id, item_obj):
        """Add the parent item IDs for one item object and their parents."""
        id2parents = self.id2parents
        if item_id in id2parents:
            return id2parents[item_id]
        parent_ids = set()
        for parent_obj in item_obj.parents:
            parent_id = parent_obj.item_id
            parent_ids.add(parent_id)
            parent_ids |= self._get_id2parents(parent_id, parent_obj)
        id2parents[item_id] = parent_ids
        return parent_ids


if __name__ == '__main__':
    PRT = None if len(sys.argv) != 1 else sys.stdout
    test_get_parent(PRT)

# Copyright (C) 2010-2019, DV Klopfenstein, H Tang et al. All rights reserved.
