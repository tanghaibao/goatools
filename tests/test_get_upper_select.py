#!/usr/bin/env python
"""Test get_all_parents vs """

from __future__ import print_function

__copyright__ = "Copyright (C) 2010-2019, DV Klopfenstein, H Tang et al. All rights reserved."

import os
import sys
import timeit
from goatools.base import get_godag
from goatools.godag.go_tasks import get_id2upperselect
from goatools.godag.prttime import prt_hms
from goatools.test_data.checks import CheckGOs
from goatools.godag.relationship_combos import RelationshipCombos


def test_get_upperselect(prt=sys.stdout):
    """Test getting parents and user-specfied ancestor relationships"""
    # Load GO-DAG
    repo = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    godag = get_godag(os.path.join(repo, 'go-basic.obo'), optional_attrs='relationship')
    run = RelationshipCombos(godag)
    run.chk_relationships_all()
    rels_combo = run.get_relationship_combos()
    print('{N} COMBINATIONS OF RELATIONSHIPS'.format(N=len(rels_combo)))

    for relidx, rels_set in enumerate(rels_combo, 1):
        print('{I}) RELATIONSHIPS[{N}]: {Rs}'.format(
            I=relidx, N=len(rels_set), Rs=' '.join(sorted(rels_set))))
        # ------------------------------------------------------------------------
        # Get all parents for all GO IDs using get_all_parents in GOTerm class
        tic = timeit.default_timer()
        # pylint: disable=line-too-long
        go2upperselect_orig = {o.item_id:get_all_upperselect(o, rels_set) for o in run.go2obj.values()}
        tic = prt_hms(tic, "Get all goobj's parents using get_all_upperselect(GOTerm)", prt)
        # ------------------------------------------------------------------------
        # Get all parents for all GO IDs using GOTerm get_all_parents
        go2upperselect_fast = get_id2upperselect(run.go2obj.values(), rels_set)
        tic = prt_hms(tic, "Get all goobj's parents using go_tasks::get_id2upperselect", prt)
        # ------------------------------------------------------------------------
        # Compare parent lists
        chkr = CheckGOs('test_get_upper_select', godag)
        chkr.chk_a2bset(go2upperselect_orig, go2upperselect_fast)  # EXPECTED, ACTUAL
        print("PASSED: get_upperselect RELATIONSHIPS[{N}]: {Rs}".format(
            N=len(rels_set), Rs=' '.join(sorted(rels_set))))

# ------------------------------------------------------------------------------------
def get_all_upperselect(goterm, relationship_set):
    """Return all parent GO IDs through both 'is_a' and all relationships."""
    # SLOW WHEN RUNNING MORE THAN ONE GO TERM: GOTerm::get_all_upperselect
    all_upper = set()
    for upper in goterm.get_goterms_upper_rels(relationship_set):
        all_upper.add(upper.item_id)
        all_upper |= get_all_upperselect(upper, relationship_set)
    return all_upper


if __name__ == '__main__':
    PRT = None if len(sys.argv) != 1 else sys.stdout
    test_get_upperselect(PRT)

# Copyright (C) 2010-2019, DV Klopfenstein, H Tang et al. All rights reserved.
