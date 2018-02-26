"""Common checks in test data."""

__copyright__ = "Copyright (C) 2010-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"


def _chk_a2bset(exp_a2bset, act_a2bset):
    assert set(exp_a2bset) == set(act_a2bset)
    for goid, exp_goids in exp_a2bset.items():
        act_goids = act_a2bset[goid]
        assert act_goids == exp_goids


# Copyright (C) 2010-2018, DV Klopfenstein, H Tang, All rights reserved.
