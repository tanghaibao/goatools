"""Test the elim GOEA algorithm (Alexa 2006, ported from topGO).

Uses a hand-built five-term DAG so that the exact set of eliminated genes is
known. The assertions check *which genes elim removed* by feeding the expected
counts back through the study's own Fisher function -- they deliberately do not
re-derive the hypergeometric p-value, which is scipy's job, not elim's.

        python test_goea_elim.py
"""

__copyright__ = "Copyright (C) 2010-present, H Tang et al., All rights reserved."

import sys

import pytest

from goatools.go_enrichment import GOEnrichmentStudy
from goatools.goea.algorithms import ElimAlgorithm, get_algorithm
from goatools.obo_parser import GODag

# root -- M --+-- A --.
#              |        >-- C -- D
#              +-- B --'
# C has two parents, and M is a grandparent that is NOT the root. The root is a
# poor test of ancestor-wide elimination: it holds every gene, so a one-sided
# Fisher test returns 1.0 whether or not anything was eliminated from it.
OBO = """format-version: 1.2

[Term]
id: GO:0000001
name: root
namespace: biological_process

[Term]
id: GO:0000006
name: term_m
namespace: biological_process
is_a: GO:0000001 ! root

[Term]
id: GO:0000002
name: term_a
namespace: biological_process
is_a: GO:0000006 ! term_m

[Term]
id: GO:0000003
name: term_b
namespace: biological_process
is_a: GO:0000006 ! term_m

[Term]
id: GO:0000004
name: term_c
namespace: biological_process
is_a: GO:0000002 ! term_a
is_a: GO:0000003 ! term_b

[Term]
id: GO:0000005
name: term_d
namespace: biological_process
is_a: GO:0000004 ! term_c
"""

ROOT, TERM_A, TERM_B, TERM_C, TERM_D, TERM_M = (
    "GO:0000001", "GO:0000002", "GO:0000003",
    "GO:0000004", "GO:0000005", "GO:0000006",
)

STUDY = ["s1", "s2", "s3", "s4", "s5"]
POP = STUDY + ["p{N}".format(N=i) for i in range(1, 16)]  # 20 genes total

# Direct annotations only; GOEnrichmentStudy propagates them up the DAG.
#   C: s1..s4 + p1     ->  5 genes, 4 of them study genes
#   A: adds p2..p6     -> 10 genes, 4 study
#   B: adds p7         ->  6 genes, 4 study
#   M: adds p8, p9     -> 13 genes, 4 study   (C's grandparent, not the root)
#   root: adds the rest-> 20 genes, 5 study
ASSOC = {
    "s1": {TERM_C}, "s2": {TERM_C}, "s3": {TERM_C}, "s4": {TERM_C},
    "p1": {TERM_C},
    "p2": {TERM_A}, "p3": {TERM_A}, "p4": {TERM_A}, "p5": {TERM_A}, "p6": {TERM_A},
    "p7": {TERM_B},
    "p8": {TERM_M}, "p9": {TERM_M},
    "s5": {ROOT},
    "p10": {ROOT}, "p11": {ROOT}, "p12": {ROOT},
    "p13": {ROOT}, "p14": {ROOT}, "p15": {ROOT},
}
C_GENES = {"s1", "s2", "s3", "s4", "p1"}  # what a significant C eliminates


@pytest.fixture(name="godag")
def _godag(tmp_path):
    fin = tmp_path / "tiny.obo"
    fin.write_text(OBO)
    return GODag(str(fin), prt=None)


def _run(godag, **kws):
    """GOEA over the tiny DAG; returns {GO: record} and the study object."""
    obj = GOEnrichmentStudy(
        POP, {g: set(v) for g, v in ASSOC.items()}, godag,
        methods=["bonferroni"], alternative="greater", log=None, **kws
    )
    return {r.GO: r for r in obj.run_study(STUDY, prt=None)}, obj


# --------------------------------------------------------------- basic wiring


def test_elim_is_registered():
    """elim is reachable by name, with topGO's defaults."""
    algo = get_algorithm("elim")
    assert isinstance(algo, ElimAlgorithm)
    assert algo.cutoff == 0.01  # topGO's default, applied raw
    assert algo.bonferroni is False


def test_elim_cutoff_kws_are_routed():
    """GOEnrichmentStudy kwargs prefixed 'elim_' configure the algorithm."""
    algo = get_algorithm("elim", elim_cutoff=0.05, elim_bonferroni=True)
    assert (algo.cutoff, algo.bonferroni) == (0.05, True)


def test_elim_rejects_bad_cutoff():
    """A cutoff outside (0, 1] is rejected up front."""
    for bad in (0.0, -0.1, 1.5):
        with pytest.raises(ValueError):
            ElimAlgorithm(cutoff=bad)


def test_elim_requires_propagated_counts(godag):
    """Without propagation a term's genes exclude its children's: elim is invalid."""
    with pytest.raises(ValueError) as exc:
        _run(godag, algorithm="elim", propagate_counts=False)
    assert "propagate_counts" in str(exc.value)


# ------------------------------------------------------- the algorithm itself


def test_propagation_sanity(godag):
    """The fixture means what the comments say it means."""
    res, _ = _run(godag)
    assert (res[TERM_C].pop_count, res[TERM_C].study_count) == (5, 4)
    assert (res[TERM_A].pop_count, res[TERM_A].study_count) == (10, 4)
    assert (res[TERM_B].pop_count, res[TERM_B].study_count) == (6, 4)
    assert (res[TERM_M].pop_count, res[TERM_M].study_count) == (13, 4)
    assert (res[ROOT].pop_count, res[ROOT].study_count) == (20, 5)
    assert TERM_D not in res or res[TERM_D].study_count == 0


def test_deepest_scored_term_is_unaffected(godag):
    """C is the deepest annotated term, so nothing has been eliminated from it."""
    classic, _ = _run(godag)
    elim, _ = _run(godag, algorithm="elim")
    assert elim[TERM_C].p_uncorrected == classic[TERM_C].p_uncorrected
    # ...and it is significant, which is what drives the rest of this test file
    assert elim[TERM_C].p_uncorrected <= 0.01


def test_significant_child_is_eliminated_from_both_parents(godag):
    """C's genes leave A and B -- and leave the population too, as topGO does."""
    elim, obj = _run(godag, algorithm="elim")
    calc = obj.pval_obj.calc_pvalue
    n_elim = len(C_GENES)  # 5
    n_elim_study = len(C_GENES.intersection(STUDY))  # 4

    # A: 10 genes / 4 study, minus C's 5 genes / 4 study
    assert elim[TERM_A].p_uncorrected == calc(
        4 - n_elim_study, len(STUDY) - n_elim_study, 10 - n_elim, len(POP) - n_elim
    )
    # B: 6 genes / 4 study, minus the same
    assert elim[TERM_B].p_uncorrected == calc(
        4 - n_elim_study, len(STUDY) - n_elim_study, 6 - n_elim, len(POP) - n_elim
    )


def test_elimination_reaches_all_ancestors_not_just_parents(godag):
    """M is C's grandparent and must still lose C's genes.

    This is the property that separates elim from a parent-only scheme
    (topGO propagates to nodesInInducedGraph, i.e. every ancestor).
    """
    classic, _ = _run(godag)
    elim, obj = _run(godag, algorithm="elim")
    calc = obj.pval_obj.calc_pvalue
    n_elim, n_elim_study = len(C_GENES), len(C_GENES.intersection(STUDY))

    exp_eliminated = calc(
        4 - n_elim_study, len(STUDY) - n_elim_study, 13 - n_elim, len(POP) - n_elim
    )
    # guard against a vacuous assertion: eliminating must actually change M
    assert exp_eliminated != classic[TERM_M].p_uncorrected
    assert elim[TERM_M].p_uncorrected == exp_eliminated
    assert elim[TERM_M].p_uncorrected > classic[TERM_M].p_uncorrected


def test_unreachable_cutoff_reduces_to_classic(godag):
    """If nothing can be significant, no genes are eliminated and elim == classic."""
    classic, _ = _run(godag)
    elim, _ = _run(godag, algorithm="elim", elim_cutoff=1e-300)
    assert set(elim) == set(classic)
    for goid, rec in classic.items():
        assert elim[goid].p_uncorrected == rec.p_uncorrected, goid


def test_elim_never_more_significant_than_classic(godag):
    """Removing genes can only ever weaken evidence for a term."""
    classic, _ = _run(godag)
    elim, _ = _run(godag, algorithm="elim")
    assert any(
        elim[g].p_uncorrected > classic[g].p_uncorrected for g in classic
    ), "fixture no longer exercises elimination"
    for goid, rec in classic.items():
        assert elim[goid].p_uncorrected >= rec.p_uncorrected, goid


def test_elim_reports_eliminated_gene_counts(godag):
    """Records carry how many genes elim removed, for transparency."""
    elim, _ = _run(godag, algorithm="elim")
    assert elim[TERM_C].elim_genes == 0
    assert elim[TERM_A].elim_genes == len(C_GENES)
    assert elim[TERM_B].elim_genes == len(C_GENES)
    assert elim[TERM_M].elim_genes == len(C_GENES)
    assert elim[ROOT].elim_genes == len(C_GENES)
    # terms that were never scored still carry the field, set to zero
    if TERM_D in elim:
        assert elim[TERM_D].elim_genes == 0


def test_bonferroni_cutoff_is_stricter(godag):
    """The paper's cutoff/N variant eliminates no more than the raw cutoff does."""
    raw, _ = _run(godag, algorithm="elim", elim_cutoff=0.05)
    adj, _ = _run(godag, algorithm="elim", elim_cutoff=0.05, elim_bonferroni=True)
    for goid in raw:
        assert adj[goid].p_uncorrected <= raw[goid].p_uncorrected, goid


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
