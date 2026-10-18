"""Test the pluggable GOEA term-scoring algorithm layer.

The 'classic' algorithm must remain a byte-for-byte no-op refactor: selecting it
explicitly, or not selecting anything at all, must give exactly the p-values
goatools produced before the algorithm layer existed.

        python test_goea_algorithms.py
"""

__copyright__ = "Copyright (C) 2010-present, H Tang et al., All rights reserved."

import os
import random
import sys
from unittest.mock import Mock, patch

import pytest

from goatools.associations import read_associations
from goatools.base import get_godag
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.goea.algorithms import (
    ALGORITHMS,
    ClassicAlgorithm,
    ElimAlgorithm,
    GoeaAlgorithm,
    get_algorithm,
)
from goatools.goea.algorithms.base import GoeaAlgoResult
from goatools.multiple_testing import calc_qval
from goatools.pvalcalc import FisherFactory, FisherScipyStats
from goatools.ratio import count_terms

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../data/"


# ------------------------------------------------------------------ registry


def test_get_algorithm_default_is_classic():
    """No algorithm requested -> classic, preserving historical behaviour."""
    assert isinstance(get_algorithm(), ClassicAlgorithm)
    assert isinstance(get_algorithm(None), ClassicAlgorithm)
    assert get_algorithm().name == "classic"
    assert "classic" in ALGORITHMS


def test_get_algorithm_accepts_name_instance_and_class():
    """A name, a ready instance, or the class itself are all acceptable."""
    obj = ClassicAlgorithm()
    assert isinstance(get_algorithm("classic"), ClassicAlgorithm)
    assert get_algorithm(obj) is obj
    assert isinstance(get_algorithm(ClassicAlgorithm), ClassicAlgorithm)


def test_get_algorithm_unknown_name():
    """An unknown algorithm names the valid choices rather than failing late."""
    with pytest.raises(ValueError) as exc:
        get_algorithm("no_such_algorithm")
    assert "classic" in str(exc.value)


def test_algorithm_base_is_abstract():
    """The base class refuses to score anything."""
    with pytest.raises(NotImplementedError):
        GoeaAlgorithm().run(None)


# ------------------------------------------------------------------ pvalcalc


def test_fisher_alternative_default_unchanged():
    """The historical default is a two-sided Fisher test."""
    pobj = FisherFactory(log=None).pval_obj
    assert pobj.alternative == "two-sided"
    # 8/10 study, 9/16 population -- see the worked example in pvalcalc.py
    assert pobj.calc_pvalue(8, 10, 9, 16) == pytest.approx(0.034965034965, abs=1e-9)


def test_fisher_alternative_greater():
    """'greater' tests over-representation only, as topGO's GOFisherTest does."""
    from scipy import stats

    pobj = FisherFactory(log=None, alternative="greater").pval_obj
    exp = stats.fisher_exact([[8, 2], [1, 5]], alternative="greater")[1]
    assert pobj.calc_pvalue(8, 10, 9, 16) == pytest.approx(exp, abs=1e-12)
    # one-sided is never larger than two-sided for an over-represented term
    two = FisherFactory(log=None).pval_obj.calc_pvalue(8, 10, 9, 16)
    assert pobj.calc_pvalue(8, 10, 9, 16) < two


def test_fisher_alternative_bad():
    """An unknown alternative is rejected at construction, not at first use."""
    with pytest.raises(ValueError) as exc:
        FisherFactory(log=None, alternative="sideways")
    assert "sideways" in str(exc.value)


# ------------------------------------------------------------- end-to-end


def _init_goea(optional_attrs=None, **kws):
    """GOEnrichmentStudy over the sample Arabidopsis data."""
    godag = get_godag(
        os.path.join(os.getcwd(), "go-basic.obo"),
        prt=None,
        optional_attrs=optional_attrs,
    )
    assoc = read_associations(ROOT + "association", "id2gos", no_top=True)
    popul_ids = [ln.rstrip() for ln in open(ROOT + "population")]
    study_ids = [ln.rstrip() for ln in open(ROOT + "study")]
    obj = GOEnrichmentStudy(
        popul_ids, assoc, godag, methods=["bonferroni"], log=None, **kws
    )
    return obj, study_ids


def _go2pval(results):
    return {r.GO: r.p_uncorrected for r in results}


def test_classic_algorithm_is_the_default():
    """Explicitly asking for 'classic' matches asking for nothing, exactly."""
    obj_def, study_ids = _init_goea()
    obj_cls, _ = _init_goea(algorithm="classic")
    assert isinstance(obj_def.algorithm, ClassicAlgorithm)

    pvals_def = _go2pval(obj_def.run_study(study_ids, prt=None))
    pvals_cls = _go2pval(obj_cls.run_study(study_ids, prt=None))
    assert pvals_def and pvals_def.keys() == pvals_cls.keys()
    for goid, pval in pvals_def.items():
        assert pvals_cls[goid] == pval, goid


def test_classic_algorithm_matches_direct_calculation():
    """The algorithm layer must not perturb the classic p-values at all."""
    obj, study_ids = _init_goea()
    results = obj.run_study(study_ids, prt=None)
    assert results
    calc_pvalue = obj.pval_obj.calc_pvalue
    for res in results:
        exp = calc_pvalue(res.study_count, res.study_n, res.pop_count, res.pop_n)
        assert res.p_uncorrected == exp, res.GO


def test_custom_algorithm_is_honoured():
    """A user-supplied algorithm object drives the reported p-values."""

    class _AllOnes(GoeaAlgorithm):
        """Scores every term 1.0 and tags it, to prove the seam is real."""

        name = "allones"

        def run(self, ctx):
            return GoeaAlgoResult(
                go2pval={go: 1.0 for go in ctx.goids},
                go2flds={go: {"tagged": True} for go in ctx.goids},
            )

    obj, study_ids = _init_goea(algorithm=_AllOnes())
    results = obj.run_study(study_ids, prt=None)
    assert results
    assert all(r.p_uncorrected == 1.0 for r in results)
    # extra fields returned by an algorithm land on the result records
    assert all(getattr(r, "tagged", False) for r in results)


def test_context_carries_relationships():
    """Ancestor traversal must be able to see the propagation relationships."""
    seen = {}

    class _Spy(GoeaAlgorithm):
        """Captures the context instead of scoring."""

        name = "spy"

        def run(self, ctx):
            seen["ctx"] = ctx
            return GoeaAlgoResult(go2pval={go: 1.0 for go in ctx.goids}, go2flds={})

    rels = {"part_of"}
    # propagating through a relationship needs the DAG loaded with them
    obj, study_ids = _init_goea(
        optional_attrs={"relationship"}, algorithm=_Spy(), relationships=rels
    )
    obj.run_study(study_ids, prt=None)
    ctx = seen["ctx"]
    assert ctx.relationships == rels
    assert ctx.study_n == len(ctx.study_ids)
    assert ctx.pop_n == obj.pop_n
    assert ctx.godag is obj.obo_dag
    assert ctx.goids


def _spy_pvalues(obj):
    """Record every table the study sends to its p-value backend."""
    spy = Mock(side_effect=obj.pval_obj.calc_pvalue)
    obj.pval_obj.calc_pvalue = spy
    return spy


@pytest.mark.parametrize("algorithm", ["classic", "elim"])
def test_identical_tables_computed_once_per_run(algorithm):
    """Every algorithm shares one table cache, which starts empty each run."""
    obj, study_ids = _init_goea(algorithm=algorithm)
    spy = _spy_pvalues(obj)
    runs = []
    for _ in range(2):
        spy.reset_mock()
        runs.append(_go2pval(obj.run_study(study_ids, prt=None)))
        tables = [c.args for c in spy.call_args_list]
        # a cache surviving the first run would leave the second with no calls
        assert tables and len(tables) == len(set(tables))
    assert runs[0] == runs[1]


def test_fdr_resampling_computes_each_table_once():
    """Count pairs recurring across resamples share one calculation."""
    obj, study_ids = _init_goea()
    term_pop = count_terms(obj.pop, obj.assoc, obj.obo_dag)
    random.seed(0)
    with patch.object(FisherScipyStats, "calc_pvalue", autospec=True,
                      side_effect=FisherScipyStats.calc_pvalue) as spy:
        distribution = calc_qval(len(study_ids), obj.pop_n, obj.pop, obj.assoc,
                                 term_pop, obj.obo_dag, T=20)
    assert len(distribution) == 20
    tables = [c.args[1:] for c in spy.call_args_list]  # drop self
    assert tables and len(tables) == len(set(tables))


@pytest.mark.parametrize("base, alternative", [(ClassicAlgorithm, "less"),
                                               (ElimAlgorithm, "greater")])
def test_fdr_null_scored_like_the_study(base, alternative):
    """Each random set uses the study's algorithm, alternative, terms and size."""
    contexts = []

    class _Spy(base):
        """Records every context it scores."""

        def run(self, ctx):
            contexts.append(ctx)
            return super().run(ctx)

    obj, study_ids = _init_goea(algorithm=_Spy(), alternative=alternative)
    spy = _spy_pvalues(obj)
    # genes outside the population must not inflate the random sets
    results = obj.run_study(study_ids + ["NOT_IN_POP_1", "NOT_IN_POP_2"],
                            prt=None, methods=["fdr"])
    assert results and all(0.0 <= r.p_fdr <= 1.0 for r in results)
    study, null = contexts[0], contexts[1:]
    assert null
    for ctx in null:
        # the same backend means the same Fisher alternative
        assert ctx.calc_pvalue.__wrapped__ is spy
        assert ctx.goids == study.goids  # terms with no study genes included
        assert ctx.study_n == len(ctx.study_ids) == study.study_n


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
