#!/usr/bin/env python3
"""Test setting edge weights for various relationships"""

__copyright__ = "Copyright (C) 2020-present, DV Klopfenstein. All rights reserved."
__author__ = "DV Klopfenstein"

from os.path import join
from sys import stdout

from goatools.base import get_godag
from goatools.semsim.termwise.wang import SsWang
from goatools.godag.consts import RELATIONSHIP_SET

from tests.utils import REPO


def test_semsim_wang(prt=stdout):
    """Test setting edge weights for various relationships"""
    # Log file
    # Check that all relationships seem in DAG are expected by SsWang
    fin_godag = join(REPO, 'go-basic.obo')
    godag_r0 = get_godag(fin_godag, prt=prt)

    passed = False
    try:
        wang = SsWang({}, godag_r0, {'part_of',})
    except RuntimeError as err:
        assert str(err) == '**ERROR: SsWang GODag not loaded with relationships', '({})'.format(err)
        passed = True
    assert passed

    wang = SsWang({}, godag_r0)
    assert wang.w_e == {'is_a': 0.8}

    wang = SsWang({}, godag_r0, rel2scf={'is_a': 0.9, 'part_of': 0.7})
    assert wang.w_e == {'is_a': 0.9}

    godag_r1 = get_godag(fin_godag, optional_attrs=['relationship'], prt=prt)
    _chk_relationships(godag_r1)
    # Run randoms
    relationships = {'part_of'}
    wang = SsWang({}, godag_r1, relationships, rel2scf={})
    assert wang.w_e == {'is_a': 0.8, 'part_of': 0.6}

    wang = SsWang({}, godag_r1, relationships, rel2scf={'is_a': 0.9, 'part_of': 0.7})
    assert wang.w_e == {'is_a': 0.9, 'part_of': 0.7}

    # pylint: disable=line-too-long
    wang = SsWang({}, godag_r1, relationships, rel2scf={'is_a': 0.9, 'part_of': 0.7, 'regulates':0.2})
    assert wang.w_e == {'is_a': 0.9, 'part_of': 0.7}

    wang = SsWang({}, godag_r1)
    assert wang.w_e == {'is_a': 0.8}

    wang = SsWang({}, godag_r1, rel2scf={'is_a': 0.9, 'part_of': 0.7})
    assert wang.w_e == {'is_a': 0.9}

    wang = SsWang({}, godag_r1, rel2scf={'is_a': 0.9, 'part_of': 0.7, 'regulates':0.2})
    assert wang.w_e == {'is_a': 0.9}

    wang = SsWang({}, godag_r1, relationships={'mock_rel'})
    assert wang.w_e == {'is_a': 0.8}
    print('**PASSED: Properly reported ERROR in relationship, mock_rel')

    wang = SsWang({}, godag_r1, rel2scf={'mock_rel':.7})
    assert wang.w_e == {'is_a': 0.8}


def _chk_relationships(godag):
    """Check that relationships in SsWang default match expected"""
    rels_all = set()
    for goterm in godag.values():
        rels_cur = goterm.relationship.keys()
        if rels_cur:
            rels_all.update(rels_cur)
    assert rels_all == RELATIONSHIP_SET, 'UNEXPECTED RELATIONSHIPS'
    print('**PASSED: EXPECTED GODag  RELATIONSHIPS: {R}'.format(R=sorted(rels_all)))
    rels_all.add('is_a')
    rels_act = set(SsWang.dflt_rel2scf.keys())
    assert rels_all == rels_act, 'BAD SsWang RELATIONSHIPS: {Rs}'.format(Rs=rels_act)
    print('**PASSED: EXPECTED SsWang RELATIONSHIPS: {R}'.format(R=sorted(rels_act)))


if __name__ == '__main__':
    test_semsim_wang()

# Copyright (C) 2020-present DV Klopfenstein. All rights reserved.
