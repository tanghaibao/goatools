#!/usr/bin/env python
"""Plot both the standard 'is_a' field and the optional 'part_of' relationship."""

from __future__ import print_function

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved."

import os
import sys
import timeit
import datetime
from goatools.base import download_go_basic_obo
from goatools.obo_parser import GODag
from goatools.gosubdag.gosubdag import GoSubDag


REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../..")

def test_gosubdag_relationships(prt=sys.stdout):
    """Plot both the standard 'is_a' field and the 'part_of' relationship."""
    goids = set([
        "GO:0032501",
        "GO:0044707",   # alt_id: GO:0032501  # BP  1011 L01 D01 B multicellular organismal process
        "GO:0050874",
        "GO:0007608",   # sensory perception of smell
        "GO:0050911"])  # detection of chemical stimulus involved in sensory perception of smell

    # Load GO-DAG: Load optional 'relationship'
    fin_obo = os.path.join(REPO, "go-basic.obo")
    download_go_basic_obo(fin_obo, prt, loading_bar=None)
    go2obj_plain = GODag(fin_obo)
    go2obj_relat = GODag(fin_obo, optional_attrs=['relationship'])

    print("\nCreate GoSubDag with GO DAG containing no relationships.")
    tic = timeit.default_timer()
    # Create Plot object; Plot both 'is_a' and optional 'part_of' relationship
    gosubdag = GoSubDag(goids, go2obj_plain, relationships=False, prt=prt)
    # gosubdag.prt_goids(gosubdag.go2obj)
    goids_plain = set(gosubdag.go2obj)
    tic = _rpt_hms(tic, len(gosubdag.go2obj))

    print("\nCreate GoSubDag while IGNORING relationships")
    # Create Plot object; Plot both 'is_a' and optional 'part_of' relationship
    gosubdag = GoSubDag(goids, go2obj_relat, relationships=False, prt=prt)
    # gosubdag.prt_goids(gosubdag.go2obj)
    goids_false = set(gosubdag.go2obj)
    tic = _rpt_hms(tic, len(gosubdag.go2obj))
    assert goids_plain == goids_false

    print("\nCreate GoSubDag while loading only the 'part_of' relationship")
    gosubdag = GoSubDag(goids, go2obj_relat, relationships=['part_of'], prt=prt)
    # gosubdag.prt_goids(gosubdag.go2obj)
    goids_part_of = set(gosubdag.go2obj)
    tic = _rpt_hms(tic, len(gosubdag.go2obj))
    assert goids_plain.intersection(goids_part_of) == goids_plain
    assert len(goids_part_of) > len(goids_plain)

    print("\nCreate GoSubDag while loading all relationships")
    gosubdag = GoSubDag(goids, go2obj_relat, relationships=True, prt=prt)
    # gosubdag.prt_goids(gosubdag.go2obj)
    goids_true = set(gosubdag.go2obj)
    tic = _rpt_hms(tic, len(gosubdag.go2obj))
    assert goids_part_of.intersection(goids_true) == goids_part_of
    assert len(goids_true) >= len(goids_part_of)

def _rpt_hms(tic, num_goids):
    """Report the elapsed time for particular events."""
    elapsed_time = str(datetime.timedelta(seconds=(timeit.default_timer()-tic)))
    print("Elapsed HMS: {HMS} {N} GO IDs".format(HMS=elapsed_time, N=num_goids))
    return timeit.default_timer()

if __name__ == '__main__':
    test_gosubdag_relationships()

# Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved.
