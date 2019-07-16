#!/usr/bin/env python
"""Test get_all_parents vs """

from __future__ import print_function

__copyright__ = "Copyright (C) 2010-2019, DV Klopfenstein, H Tang et al. All rights reserved."

import os
import sys
import timeit
from goatools.base import get_godag
from goatools.godag.go_tasks import get_id2upperselect
from goatools.test_data.godag_timed import prt_hms
from goatools.test_data.checks import CheckGOs
from goatools.godag.consts import RELATIONSHIP_SET


def test_get_upperselect(prt=sys.stdout):
    """Semantic Similarity test for Issue #86."""
    # Load GO-DAG
    fin_obo = "go-basic.obo"
    repo = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    godag = get_godag(os.path.join(repo, fin_obo), optional_attrs='relationship')
    go2obj = {go:o for go, o in godag.items() if go == o.id}
    chkr = CheckGOs('test_get_upper_select', go2obj)
    rels_all = _get_rels_all(go2obj)
    rels_combo = _get_rel_combos(rels_all)
    for relidx, rels_set in enumerate(rels_combo, 1):
        print('{I}) RELATIONSHIPS: {Rs}'.format(I=relidx, Rs=' '.join(sorted(rels_set))))
        # ------------------------------------------------------------------------
        # Get all parents for all GO IDs using get_all_parents in GOTerm class
        tic = timeit.default_timer()
        go2upperselect_orig = {o.item_id:get_all_upperselect(o, rels_set) for o in go2obj.values()}
        tic = prt_hms(tic, "Get all goobj's parents using get_all_upperselect(GOTerm)", prt)
        # ------------------------------------------------------------------------
        # Get all parents for all GO IDs using GOTerm get_all_parents
        go2upperselect_fast = get_id2upperselect(go2obj.values(), rels_set)
        tic = prt_hms(tic, "Get all goobj's parents using go_tasks::get_id2upperselect", prt)
        # ------------------------------------------------------------------------
        # Compare parent lists
        chkr.chk_a2bset(go2upperselect_orig, go2upperselect_fast)
        print("PASSED: get_upperselect test")

def _get_rel_combos(rels_all):
    """Get all combinations of all lengths of relationships"""
    rels_combo = []
    num_rels = len(rels_all)
    print('GODAG relationships[{N}]: {Rs}'.format(N=num_rels, Rs=rels_all))
    for cnt in range(2**num_rels):
        idxs = [i for i, v in enumerate('{N:0{L}b}'.format(N=cnt, L=num_rels)) if v == '1']
        if idxs:
            rels_cur = set(rels_all[i] for i in idxs)
            rels_combo.append(rels_cur)
            # print('{N:0{L}b}'.format(N=cnt, L=num_rels), idxs, rels_cur)
    return rels_combo

def _get_rels_all(go2obj):
    """Check that the list of relationships in consts is same as found in GODAG"""
    rels_all = set()
    for obj in go2obj.values():
        rels_all.update(obj.relationship.keys())
    assert rels_all == RELATIONSHIP_SET, rels_all.symmetric_difference(RELATIONSHIP_SET)
    return sorted(rels_all)

# ------------------------------------------------------------------------------------
def get_all_upperselect(goterm, relationship_set):
    """Return all parent GO IDs through both 'is_a' and all relationships."""
    # SLOW: GOTerm::get_all_upperselect
    all_upper = set()
    for upper in goterm.get_goterms_upper_rels(relationship_set):
        all_upper.add(upper.item_id)
        all_upper |= upper.get_all_upper()
    return all_upper


if __name__ == '__main__':
    PRT = None if len(sys.argv) != 1 else sys.stdout
    test_get_upperselect(PRT)

# Copyright (C) 2010-2019, DV Klopfenstein, H Tang et al. All rights reserved.
