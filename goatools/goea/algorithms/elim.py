# -*- coding: UTF-8 -*-
"""The elim GOEA algorithm: decorrelate the GO graph by eliminating genes.

Alexa A, Rahnenfuehrer J, Lengauer T (2006). "Improved scoring of functional
groups from gene expression data by decorrelating GO graph structure."
Bioinformatics 22(13):1600-1607.

Ported from topGO's `.sigGroups.elim` (R/topGOalgo.R). The GO DAG is walked
bottom-up. When a term is found significant, the genes annotated to it are
marked "eliminated" in every one of its ancestors, so a broad parent term can
no longer look enriched merely by inheriting the genes of a specific child.
"""

__copyright__ = "Copyright (C) 2010-present, H Tang et al., All rights reserved."

import collections as cx

from goatools.godag.go_tasks import get_go2ancestors, get_go2parents
from goatools.goea.algorithms.base import GoeaAlgorithm, GoeaAlgoResult


class ElimAlgorithm(GoeaAlgorithm):
    """Score GO terms bottom-up, eliminating the genes of significant children.

    cutoff     -- a term is "significant", and so eliminates its genes from its
                  ancestors, when its p-value is <= this. topGO's default is
                  0.01, applied raw.
    bonferroni -- if True, use cutoff/N instead, where N is the number of terms
                  scored. The 2006 paper describes this Bonferroni-style
                  adjustment, but topGO ships with the line commented out and
                  uses the raw cutoff, so it is opt-in here.
    """

    name = "elim"

    # Elimination must be driven by over-representation only. Under a two-sided
    # test a significantly *depleted* term would eliminate its genes from its
    # ancestors, which inverts what the algorithm is for. topGO's GOFisherTest
    # is one-sided for the same reason.
    default_alternative = "greater"

    # topGO's vignette: for topology-aware algorithms "the tests are therefore
    # not independent and the multiple testing theory does not directly apply...
    # we like to interpret the p-values returned by these methods as corrected".
    pvals_precorrected = True

    def __init__(self, cutoff=0.01, bonferroni=False):
        if not 0.0 < cutoff <= 1.0:
            raise ValueError(
                "elim cutoff must fall in (0, 1]; GOT({C})".format(C=cutoff)
            )
        self.cutoff = cutoff
        self.bonferroni = bonferroni

    def run(self, ctx):
        """Score every GO ID in ctx.goids, bottom-up with gene elimination."""
        if not ctx.propagate_counts:
            raise ValueError(
                "THE elim ALGORITHM REQUIRES propagate_counts=True: it assumes a "
                "term's genes include those of its descendants"
            )
        # topGO sets aside terms with no significant genes, scores only the rest,
        # and reports p=1 for the remainder. Note this runs over ALL annotated
        # terms, not just ctx.goids: a caller-supplied `selected_goids` restricts
        # what is *reported*, but elimination still needs the whole DAG.
        scored = {go for go, items in ctx.go2studyitems.items() if items}
        # Terms with no significant genes keep p=1 and eliminate nothing; every
        # record carries elim_genes so the field is never sometimes-missing.
        go2pval = dict.fromkeys(ctx.goids, 1.0)
        go2flds = {go: {"elim_genes": 0} for go in ctx.goids}
        if not scored:
            return GoeaAlgoResult(go2pval=go2pval, go2flds=go2flds)

        cutoff = self.cutoff / len(scored) if self.bonferroni else self.cutoff
        go2ancestors, go2depth = self._get_topology(scored, ctx)

        study_ids, study_n, pop_n = ctx.study_ids, ctx.study_n, ctx.pop_n
        calc_pvalue = ctx.calc_pvalue
        # GO ID -> genes eliminated from it by its significant descendants
        go2elim = cx.defaultdict(set)

        # Deepest first. Ordering by longest-path-from-root reproduces topGO's
        # level order, and guarantees every descendant is scored before its
        # ancestors: if u is an ancestor of v then depth(v) > depth(u) strictly.
        # Terms sharing a depth are never ancestor-related, so their relative
        # order cannot matter -- the result is independent of the tie-break.
        for goid in sorted(scored, key=lambda go: (-go2depth[go], go)):
            study_items = ctx.go2studyitems[goid]
            pop_items = ctx.go2popitems.get(goid, set())
            elim = go2elim.get(goid)

            if elim:
                # Eliminated genes leave the term AND the population, exactly as
                # topGO's numMembers/numAllMembers do. elim is a subset of the
                # term's genes, so pop_count and pop_n shrink by the same count.
                n_elim = len(elim)
                pval = calc_pvalue(
                    len(study_items - elim),
                    study_n - len(elim & study_ids),
                    len(pop_items) - n_elim,
                    pop_n - n_elim,
                )
            else:
                n_elim = 0
                pval = calc_pvalue(
                    len(study_items), study_n, len(pop_items), pop_n
                )

            if goid in go2pval:
                go2pval[goid] = pval
                go2flds[goid]["elim_genes"] = n_elim

            if pval <= cutoff:
                # Mark this term's genes eliminated in ALL of its ancestors,
                # not merely its parents (topGO: nodesInInducedGraph).
                for ancestor in go2ancestors.get(goid, ()):
                    go2elim[ancestor].update(pop_items)

        return GoeaAlgoResult(go2pval=go2pval, go2flds=go2flds)

    @staticmethod
    def _get_topology(scored, ctx):
        """Ancestors and longest-path depths, over the terms being scored.

        Both must use the same relationships that propagated the annotations,
        or elimination would flow along edges the gene counts never did.
        """
        godag = ctx.godag
        terms = {godag[go] for go in scored if go in godag}
        go2ancestors = {
            go: ancestors & scored
            for go, ancestors in get_go2ancestors(terms, ctx.relationships).items()
        }
        go2parents = {
            go: parents & scored
            for go, parents in get_go2parents(
                {go: godag[go] for go in scored if go in godag}, ctx.relationships
            ).items()
        }
        return go2ancestors, ElimAlgorithm._get_depths(scored, go2parents)

    @staticmethod
    def _get_depths(scored, go2parents):
        """Longest distance from a root, computed iteratively (GO DAGs are deep)."""
        go2depth = {}
        for goid in scored:
            if goid in go2depth:
                continue
            stack = [goid]
            while stack:
                cur = stack[-1]
                if cur in go2depth:
                    stack.pop()
                    continue
                parents = go2parents.get(cur)
                if not parents:
                    go2depth[cur] = 0
                    stack.pop()
                    continue
                pending = [p for p in parents if p not in go2depth]
                if pending:
                    stack.extend(pending)
                else:
                    go2depth[cur] = 1 + max(go2depth[p] for p in parents)
                    stack.pop()
        return go2depth
