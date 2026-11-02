"""Issue #145: blank GO IDs in an association are not reported as missing GO IDs."""

import io
import os

from goatools.anno.update_association import _chk_goids_notfound, clean_anno, update_association
from goatools.obo_parser import GODag

MINI_OBO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "mini_obo.obo")


def test_i145_blank_goid_not_reported():
    """'' (e.g. from splitting "GO:1;GO:2;" on ";") is not a GO ID; real misses still are."""
    prt = io.StringIO()
    _chk_goids_notfound({"", "GO:0000001"}, {"GO:0000001"}, prt)
    assert prt.getvalue() == ""
    _chk_goids_notfound({"", "GO:9999999", "GO:9999998"}, {"GO:0000001"}, prt)
    assert prt.getvalue() == "2 GO IDs NOT FOUND IN THE GO DAG: GO:9999998 GO:9999999\n"


def test_i145_update_association_with_blank():
    """Counts still propagate, with no missing-GO warning for the blank entry."""
    godag = GODag(MINI_OBO, prt=None)
    assoc = {"g1": set("GO:0000004;".split(";"))}  # b, under A, under top
    prt = io.StringIO()
    update_association(assoc, godag, prt=prt)
    assert "NOT FOUND" not in prt.getvalue()
    assert {"GO:0000004", "GO:0000003", "GO:0000001"} <= assoc["g1"]


def test_i145_clean_anno_with_blank():
    """clean_anno does not count the blank entry as a GO ID missing from the DAG."""
    godag = GODag(MINI_OBO, prt=None)
    gene2gos, _, goids_notfound = clean_anno({"g1": {"", "GO:0000004"}}, godag, prt=None)
    assert gene2gos == {"g1": {"GO:0000004"}}
    assert goids_notfound == set()
