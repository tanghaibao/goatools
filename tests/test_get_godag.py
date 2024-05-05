#!/usr/bin/env python
"""Test downloading GO DAG."""

import os
import sys

from goatools.base import get_godag


def test_godag(prt=sys.stdout):
    """
    Test downloading GO DAG.
    """
    cwd = os.getcwd()
    for fin_obo in ["go-basic.obo", "goslim_generic.obo"]:
        fin_full = os.path.join(cwd, fin_obo)
        os.system("rm -f {OBO}".format(OBO=fin_obo))
        godag = get_godag(fin_full, prt)  # Get GODag object
        assert godag, "GO-DAG({OBO}) NOT PROPERLY LOADED".format(OBO=fin_obo)


def test_godag_ecs():
    """
    Test downloading GO DAG and extract ECs.
    """
    godag = get_godag("go-basic.obo", optional_attrs="xref")
    go = godag["GO:0000010"]
    ecs = [x for x in go.xref if x.startswith("EC:")]
    assert "EC:2.5.1.30" in ecs


if __name__ == "__main__":
    test_godag()
    test_godag_ecs()
