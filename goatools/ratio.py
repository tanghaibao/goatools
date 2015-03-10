#!/usr/bin/env python
# -*- coding: UTF-8 -*-


from collections import defaultdict


def count_terms(geneset, assoc, obo_dag):
    """count the number of terms in the study group
    """
    term_cnt = defaultdict(int)
    for gene in (g for g in geneset if g in assoc):
        for x in assoc[gene]:
            if x in obo_dag:
                term_cnt[obo_dag[x].id] += 1

    return term_cnt


def is_ratio_different(min_ratio, study_go, study_n, pop_go, pop_n):
    """
    check if the ratio go /n is different between the study group and
    the population
    """
    if min_ratio is None:
        return True
    s = float(study_go) / study_n
    p = float(pop_go) / pop_n
    if s > p:
        return s / p > min_ratio
    return p / s > min_ratio
