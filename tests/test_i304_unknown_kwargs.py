"""Issue #304: a misspelled keyword such as relationship= raises instead of being ignored."""

import os

import pytest

from goatools.go_enrichment import GOEnrichmentStudy
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.obo_parser import GODag

MINI_OBO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "mini_obo.obo")


@pytest.fixture(name="godag")
def _godag():
    return GODag(MINI_OBO, optional_attrs={"relationship"}, prt=None)


def test_i304_gosubdag(godag):
    """relationship= used to fall into **kws silently, dropping the relationships."""
    with pytest.raises(TypeError, match="'relationship'. Did you mean 'relationships'"):
        GoSubDag(["GO:0000004"], godag, relationship={"part_of"}, prt=None)
    assert GoSubDag(["GO:0000004"], godag, relationships=True, children=True, prt=None)


def test_i304_goea(godag):
    """Same for GOEnrichmentStudy; algorithm-named args such as elim_cutoff still pass."""
    args = (["g1"], {"g1": {"GO:0000004"}}, godag)
    with pytest.raises(TypeError, match="'relationship'. Did you mean 'relationships'"):
        GOEnrichmentStudy(*args, relationship={"part_of"}, log=None)
    with pytest.raises(TypeError, match="'method'. Did you mean 'methods'"):
        GOEnrichmentStudy(*args, method=["holm"], log=None)
    goea = GOEnrichmentStudy(*args, relationships={"part_of"}, algorithm="elim", elim_cutoff=0.05, log=None)
    assert (goea.relationships, goea.algorithm.cutoff) == ({"part_of"}, 0.05)
