#!/usr/bin/env python
"""Test get_all_parents vs """

from __future__ import print_function

__copyright__ = "Copyright (C) 2010-2019, DV Klopfenstein, H Tang et al. All rights reserved."

import os
import sys
import timeit
from goatools.base import get_godag
from goatools.godag.go_tasks import get_id2lowerselect
from goatools.test_data.godag_timed import prt_hms
from goatools.test_data.checks import CheckGOs
from goatools.godag.consts import RELATIONSHIP_SET


def test_get_lowerselect(prt=sys.stdout):
    """Test getting parents and user-specfied ancestor relationships"""
    # Load GO-DAG
    run = _Run('go-basic.obo')
    rels_combo = run.get_rel_combos()

    for relidx, rels_set in enumerate(rels_combo, 1):
        print('{I}) RELATIONSHIPS[{N}]: {Rs}'.format(
            I=relidx, N=len(rels_set), Rs=' '.join(sorted(rels_set))))
        # ------------------------------------------------------------------------
        # Get all parents for all GO IDs using get_all_parents in GOTerm class
        tic = timeit.default_timer()
        # pylint: disable=line-too-long
        go2lowerselect_orig = {o.item_id:get_all_lowerselect(o, rels_set) for o in run.go2obj.values()}
        tic = prt_hms(tic, "Get all goobj's parents using get_all_lowerselect(GOTerm)", prt)
        # ------------------------------------------------------------------------
        # Get all parents for all GO IDs using GOTerm get_all_parents
        go2lowerselect_fast = get_id2lowerselect(run.go2obj.values(), rels_set)
        tic = prt_hms(tic, "Get all goobj's parents using go_tasks::get_id2lowerselect", prt)
        # ------------------------------------------------------------------------
        # Compare parent lists
        run.chkr.chk_a2bset(go2lowerselect_orig, go2lowerselect_fast)  # EXPECTED, ACTUAL
        print("PASSED: get_lowerselect RELATIONSHIPS[{N}]: {Rs}".format(
            N=len(rels_set), Rs=' '.join(sorted(rels_set))))

# pylint: disable=too-few-public-methods,old-style-class
class _Run:
    """Test getting parents and user-specfied ancestor relationships"""

    def __init__(self, fin_obo):
        _repo = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
        _godag = get_godag(os.path.join(_repo, fin_obo), optional_attrs='relationship')
        self.go2obj = {go:o for go, o in _godag.items() if go == o.id}
        self.chkr = CheckGOs('test_get_lower_select', self.go2obj)
        self.rels_all = self._get_rels_all()

    def get_rel_combos(self):
        """Get all combinations of all lengths of relationships"""
        rels_combo = []
        num_rels = len(self.rels_all)
        print('GODAG relationships[{N}]: {Rs}'.format(N=num_rels, Rs=self.rels_all))
        for cnt in range(2**num_rels):
            idxs = [i for i, v in enumerate('{N:0{L}b}'.format(N=cnt, L=num_rels)) if v == '1']
            if idxs:
                rels_cur = set(self.rels_all[i] for i in idxs)
                rels_combo.append(rels_cur)
                # print('{N:0{L}b}'.format(N=cnt, L=num_rels), idxs, rels_cur)
        return rels_combo

    def _get_rels_all(self):
        """Check that the list of relationships in consts is same as found in GODAG"""
        rels_all = set()
        for obj in self.go2obj.values():
            rels_all.update(obj.relationship.keys())
        assert rels_all == RELATIONSHIP_SET, rels_all.symmetric_difference(RELATIONSHIP_SET)
        return sorted(rels_all)

# ------------------------------------------------------------------------------------
def get_all_lowerselect(goterm, relationship_set):
    """Return all parent GO IDs through both 'is_a' and all relationships."""
    # SLOW WHEN RUNNING MORE THAN ONE GO TERM: GOTerm::get_all_lowerselect
    all_lower = set()
    for lower in goterm.get_goterms_lower_rels(relationship_set):
        all_lower.add(lower.item_id)
        all_lower |= get_all_lowerselect(lower, relationship_set)
    return all_lower


if __name__ == '__main__':
    PRT = None if len(sys.argv) != 1 else sys.stdout
    test_get_lowerselect(PRT)

# Copyright (C) 2010-2019, DV Klopfenstein, H Tang et al. All rights reserved.
