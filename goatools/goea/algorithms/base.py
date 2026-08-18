# -*- coding: UTF-8 -*-
"""Base class and shared context for GOEA term-scoring algorithms.

A GOEA algorithm turns the per-term gene sets of a study into one uncorrected
p-value per GO term. The default, `classic`, scores every term independently.
Algorithms that exploit the GO DAG topology (`elim`) instead score terms in a
topology-driven order, letting the result for one term influence another.
"""

__copyright__ = "Copyright (C) 2010-present, H Tang et al., All rights reserved."

import collections as cx


# Everything an algorithm may need to score the terms of one study.
#   goids         -- GO IDs to score
#   go2studyitems -- GO ID -> set of study genes annotated to it (propagated)
#   go2popitems   -- GO ID -> set of population genes annotated to it (propagated)
#   study_ids     -- the study genes found in the population
#   study_n       -- len(study_ids)
#   pop_n         -- size of the population
#   godag         -- GODag, for terms' levels and ancestors
#   relationships -- relationships used when propagating counts; ancestor
#                    traversal must use the same set or elimination and
#                    annotation would disagree
#   propagate_counts -- whether the annotations were propagated up the DAG.
#                    Topology-aware algorithms require this; without it a
#                    term's gene set is not a superset of its children's
#   calc_pvalue   -- fnc(study_count, study_n, pop_count, pop_n) -> pvalue
#   log           -- file-like or None
GoeaContext = cx.namedtuple(
    "GoeaContext",
    "goids go2studyitems go2popitems study_ids study_n pop_n "
    "godag relationships propagate_counts calc_pvalue log",
)

# go2pval  -- GO ID -> uncorrected p-value; must cover every ID in ctx.goids
# go2flds  -- GO ID -> dict of extra fields to set on the result record; may be
#             empty, and need not cover every GO ID
GoeaAlgoResult = cx.namedtuple("GoeaAlgoResult", "go2pval go2flds")


class GoeaAlgorithm:
    """Abstract base for GOEA term-scoring algorithms."""

    name = None

    # Fisher alternative this algorithm needs when the user has not chosen one.
    # None leaves goatools' default ("two-sided") alone.
    default_alternative = None

    # True if the algorithm's p-values already account for the dependence
    # between terms, so stacking a multiple-testing correction is questionable.
    pvals_precorrected = False

    def run(self, ctx):
        """Score ctx.goids. Returns a GoeaAlgoResult."""
        raise NotImplementedError(
            "{C} DOES NOT IMPLEMENT run()".format(C=type(self).__name__)
        )

    def __str__(self):
        return self.name
