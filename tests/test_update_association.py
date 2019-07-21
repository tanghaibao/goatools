#!/usr/bin/env python
"""Compare new propagate counts function with original function. Test assc results is same."""

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import os
import timeit
from goatools.obo_parser import GODag
from goatools.base import get_godag
from goatools.associations import dnld_assc
from goatools.anno.update_association import update_association
from goatools.test_data.godag_timed import prt_hms

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def test_update_association():
    """Compare new propagate counts function with original function. Test assc results is same."""
    assc_name = "goa_human.gaf" # gene_association.fb gene_association.mgi
    obo = os.path.join(REPO, "go-basic.obo")
    tic = timeit.default_timer()
    godag1 = get_godag(obo)
    godag2 = get_godag(obo)
    tic = prt_hms(tic, "Created two GODags: One for original and one for new propagate counts")
    assc_orig = dnld_assc(os.path.join(REPO, assc_name), godag1)
    tic = prt_hms(tic, "Associations Read")
    assc1 = {g:set(gos) for g, gos in assc_orig.items()}
    assc2 = {g:set(gos) for g, gos in assc_orig.items()}
    tic = prt_hms(tic, "Associations Copied: One for original and one for new")
    godag1.update_association(assc1)
    tic = prt_hms(tic, "ORIG: godag.update_association(assc)")
    update_association(assc2, godag2)
    tic = prt_hms(tic, "NEW:  update_association(go2obj, assc_orig)")
    _chk_assc(assc1, assc2)
    _chk_godag(godag1, obo)
    _chk_godag(godag2, obo)

def _chk_godag(go2obj_act, obo):
    """Check that the update_association function did not alter godag."""
    go2obj_exp = GODag(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../..", obo))
    assert len(go2obj_act) == len(go2obj_exp)
    assert set(go2obj_act) == set(go2obj_exp)
    for go_act, obj_act in go2obj_act.items():
        obj_exp = go2obj_exp[go_act]
        act_gos = set(o.id for o in obj_act.parents)
        exp_gos = set(o.id for o in obj_exp.parents)
        assert act_gos == exp_gos, "\nACT: {A}\nEXP: {E}".format(A=act_gos, E=exp_gos)

def _chk_assc(assc1, assc2):
    """Ensure the two associations are the same."""
    assert len(assc1) == len(assc2)
    assert set(assc1) == set(assc2)
    for gene, gos1 in assc1.items():
        gos2 = assc2[gene]
        assert gos1 == gos2


if __name__ == '__main__':
    test_update_association()

# Copyright (C) 2016-2019, DV Klopfenstein, H Tang, All rights reserved.
