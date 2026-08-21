# -*- coding: UTF-8 -*-
"""The classic GOEA algorithm: score every GO term independently."""

__copyright__ = "Copyright (C) 2010-present, H Tang et al., All rights reserved."

from goatools.goea.algorithms.base import GoeaAlgorithm, GoeaAlgoResult


class ClassicAlgorithm(GoeaAlgorithm):
    """Score each GO term on its own, ignoring the GO DAG topology.

    This is the historical goatools behaviour and remains the default.
    """

    name = "classic"

    def run(self, ctx):
        """One independent test per GO term."""
        calc_pvalue = ctx.calc_pvalue
        study_n, pop_n = ctx.study_n, ctx.pop_n
        go2studyitems, go2popitems = ctx.go2studyitems, ctx.go2popitems
        empty = set()
        go2pval = {
            goid: calc_pvalue(
                len(go2studyitems.get(goid, empty)),
                study_n,
                len(go2popitems.get(goid, empty)),
                pop_n,
            )
            for goid in ctx.goids
        }
        return GoeaAlgoResult(go2pval=go2pval, go2flds={})
