"""Test that GoSearch searches only GO term names by default (matching geneontology.org).

The issue reported that GoSearch returned more GO terms than geneontology.org because
it searched definition and comment fields in addition to names. The fix changes the
default to search only 'name', with optional extended search via goa_srch_hdrs.
"""

import os
import re

from goatools.go_search import GoSearch

__copyright__ = "Copyright (C) 2010-present, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

TESTDIR = os.path.dirname(os.path.abspath(__file__))
MINI_OBO = os.path.join(TESTDIR, "data", "mini_obo.obo")


def _all_goids(obo_dag):
    """Return set of all GO IDs in a GODag (used to search across all terms)."""
    return set(obo_dag.keys())


def test_gosearch_default_name_only():
    """Default GoSearch should only match GO term *names*, not definitions or synonyms.

    mini_obo.obo:
      GO:0000001 name="top"  def contains "mannose"  synonym "1,6-alpha-mannosyltransferase"
      GO:0000002 name="B"    def contains "maltose"  synonyms contain "maltose breakdown", etc.
    """
    srch = GoSearch(MINI_OBO, go2items={})
    all_gos = _all_goids(srch.obo_dag)

    # Searching a word that appears only in the *definition* must return nothing
    mannose_pat = re.compile(r"mannose", flags=re.IGNORECASE)
    assert srch.get_matching_gos(mannose_pat, gos=all_gos) == set(), \
        "Default search should NOT match definitions (mannose is in def, not name)"

    # Searching a word that appears only in *synonyms* must return nothing
    maltose_pat = re.compile(r"maltose", flags=re.IGNORECASE)
    assert srch.get_matching_gos(maltose_pat, gos=all_gos) == set(), \
        "Default search should NOT match synonyms (maltose is in synonym, not name)"

    # Searching a word that IS the term *name* must return that term
    top_pat = re.compile(r"^top$", flags=re.IGNORECASE)
    assert "GO:0000001" in srch.get_matching_gos(top_pat, gos=all_gos), \
        "Default search should match GO term names"


def test_gosearch_extended_defn():
    """Extended search with defn header should also match definitions."""
    srch = GoSearch(MINI_OBO, go2items={}, goa_srch_hdrs=['name', 'defn'])
    all_gos = _all_goids(srch.obo_dag)

    mannose_pat = re.compile(r"mannose", flags=re.IGNORECASE)
    result = srch.get_matching_gos(mannose_pat, gos=all_gos)
    assert "GO:0000001" in result, \
        "Extended search (name+defn) should match GO:0000001 via its definition"

    maltose_pat = re.compile(r"maltose", flags=re.IGNORECASE)
    result = srch.get_matching_gos(maltose_pat, gos=all_gos)
    assert "GO:0000002" in result, \
        "Extended search (name+defn) should match GO:0000002 via its definition"


def test_gosearch_extended_synonym():
    """Extended search with synonym header should also match synonym text."""
    srch = GoSearch(MINI_OBO, go2items={}, goa_srch_hdrs=['name', 'synonym'])
    all_gos = _all_goids(srch.obo_dag)

    # "maltose breakdown" is a synonym of GO:0000002 (name="B")
    maltose_pat = re.compile(r"maltose", flags=re.IGNORECASE)
    result = srch.get_matching_gos(maltose_pat, gos=all_gos)
    assert "GO:0000002" in result, \
        "Extended search (name+synonym) should match GO:0000002 via its synonym"

    # The name "B" itself must still match
    b_pat = re.compile(r"^B$")
    assert "GO:0000002" in srch.get_matching_gos(b_pat, gos=all_gos), \
        "Extended search should still match GO term names"


if __name__ == "__main__":
    test_gosearch_default_name_only()
    test_gosearch_extended_defn()
    test_gosearch_extended_synonym()
    print("All tests passed.")
