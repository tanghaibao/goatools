"""Issue #304: a misspelled keyword such as relationship= raises instead of being ignored."""

import os

import pytest

from goatools.go_enrichment import GOEnrichmentStudy
from goatools.goea.algorithms import ClassicAlgorithm
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.gosubdag.plot.gosubdag_plot import GoSubDagPlot
from goatools.gosubdag.plot.plot import plt_goids
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


def test_i304_custom_algorithm_kwargs(godag):
    """Args named after a custom GoeaAlgorithm class still reach it."""

    class Custom(ClassicAlgorithm):
        """A user-supplied algorithm with its own option."""

        name = "custom"

        def __init__(self, threshold=0.5):
            self.threshold = threshold

    goea = GOEnrichmentStudy(["g1"], {"g1": {"GO:0000004"}}, godag,
                             algorithm=Custom, custom_threshold=0.1, log=None)
    assert goea.algorithm.threshold == 0.1


@pytest.fixture(name="no_draw")
def _no_draw(monkeypatch):
    """Build plots without running graphviz."""
    monkeypatch.setattr(GoSubDagPlot, "plt_dag", lambda self, fout_img, engine="pydot": None)


PART_OF_OBO = """format-version: 1.2

[Term]
id: GO:0000001
name: top
namespace: biological_process

[Term]
id: GO:0000002
name: whole
namespace: biological_process
is_a: GO:0000001

[Term]
id: GO:0000003
name: part
namespace: biological_process
is_a: GO:0000001
relationship: part_of GO:0000002

[Typedef]
id: part_of
name: part of
"""


def test_i304_plt_goids(tmp_path, no_draw):
    """plt_goids passes relationships to GoSubDag and plot-only options to the plot."""
    fin = tmp_path / "part_of.obo"
    fin.write_text(PART_OF_OBO)
    dag = GODag(str(fin), optional_attrs={"relationship"}, prt=None)
    plot = plt_goids(GoSubDag(None, dag, prt=None), str(tmp_path / "go.png"), ["GO:0000003"],
                     relationships={"part_of"}, mark_alt_id=True)
    assert "GO:0000002" in plot.gosubdag.go2obj  # reached only through part_of


def test_i304_grouper_prune(tmp_path, monkeypatch, no_draw):
    """Grouped plots with prune= used to forward dst_srcs_list to GoSubDag."""
    from goatools.grouper.grprobj import Grouper
    from goatools.grouper.grprplt import GrouperPlot
    from goatools.grouper.hdrgos import HdrgosSections
    from tests.test_grprobj import _get_gosubdag, _get_grprdflt, get_data0

    monkeypatch.chdir(tmp_path)
    grprdflt = _get_grprdflt()
    hdrobj = HdrgosSections(grprdflt.gosubdag, grprdflt.hdrgos_dflt, sections=None, hdrgos=None)
    grprobj = Grouper("prune", get_data0(), hdrobj, _get_gosubdag(), go2nt=None)
    hdrgo = next(iter(grprobj.get_hdrgos()))
    plot = GrouperPlot(grprobj).get_gosubdagplot(prune=True, hdrgo=hdrgo)
    assert plot.gosubdag.go2obj
