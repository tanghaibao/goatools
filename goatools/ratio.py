# -*- coding: UTF-8 -*-
"""Manage associations ratios."""

__copyright__ = "Copyright (C) 2010-2018, H Tang et al., All rights reserved."
__author__ = "various"

from collections import defaultdict, Counter


def count_terms(geneset, assoc, obo_dag):
    """count the number of terms in the study group
    """
    term_cnt = Counter()
    for gene in (g for g in geneset if g in assoc):
        for goid in assoc[gene]:
            if goid in obo_dag:
                term_cnt[obo_dag[goid].id] += 1

    return term_cnt

def get_terms(desc, geneset, assoc, obo_dag, log):
    """Get the terms in the study group
    """
    _chk_gene2go(assoc)
    term2itemids = defaultdict(set)
    genes = [g for g in geneset if g in assoc]
    for gene in genes:
        for goid in assoc[gene]:
            if goid in obo_dag:
                term2itemids[obo_dag[goid].id].add(gene)
    if log is not None:
        num_stu = len(genes)
        num_pop = len(geneset)
        perc = 100.0*num_stu/num_pop if num_pop != 0 else 0.0
        log.write("{P:3.0f}% {N:>6,} of {M:>6,} {DESC} items found in association\n".format(
            DESC=desc, N=num_stu, M=num_pop, P=perc))
    return term2itemids

def is_ratio_different(min_ratio, study_go, study_n, pop_go, pop_n):
    """
    check if the ratio go /n is different between the study group and
    the population
    """
    if min_ratio is None:
        return True
    stu_ratio = float(study_go) / study_n
    pop_ratio = float(pop_go) / pop_n
    if stu_ratio == 0.0:
        stu_ratio = 0.0000001
    if pop_ratio == 0.0:
        pop_ratio = 0.0000001
    if stu_ratio > pop_ratio:
        return stu_ratio / pop_ratio > min_ratio
    return pop_ratio / stu_ratio > min_ratio

def _chk_gene2go(assoc):
    """Check that associations is gene2go, not go2gene."""
    if not assoc:
        raise RuntimeError("NO ITEMS FOUND IN ASSOCIATIONS {A}".format(A=assoc))
    for key in assoc:
        if isinstance(key, str) and key[:3] == "GO:":
            raise Exception("ASSOCIATIONS EXPECTED TO BE gene2go, NOT go2gene: {EX}".format(
                EX=assoc.items()[:2]))
        return

# Copyright (C) 2010-2018, H Tang et al., All rights reserved.
