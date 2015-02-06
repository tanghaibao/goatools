#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
python %prog study.file population.file gene-association.file

This program returns P-values for functional enrichment in a cluster of
study genes using Fisher's exact test, and corrected for multiple testing
(including Bonferroni, Holm, Sidak, and false discovery rate)
"""
from __future__ import absolute_import
import sys
import collections
import os.path as op
import fisher
from .multiple_testing import Bonferroni, Sidak, HolmBonferroni, FDR, calc_qval


class GOEnrichmentRecord(object):
    """Represents one result (from a single GOTerm) in the GOEnrichmentStudy
    """
    _fields = "id enrichment description ratio_in_study ratio_in_pop"\
              " p_uncorrected p_bonferroni p_holm p_sidak p_fdr".split()

    def __init__(self, **kwargs):
        for f in self._fields:
            self.__setattr__(f, "n.a.")

        for k, v in kwargs.items():
            assert k in self._fields, "invalid field name %s" % k
            self.__setattr__(k, v)

        self.goterm = None  # the reference to the GOTerm

    def __setattr__(self, name, value):
        self.__dict__[name] = value

    def __str__(self, indent=False):
        field_data = [self.__dict__[f] for f in self._fields]
        field_formatter = ["%s"] * 3 + ["%d/%d"] * 2 + ["%.3g"] * 5
        assert len(field_data) == len(field_formatter)

        # default formatting only works for non-"n.a" data
        for i, f in enumerate(field_data):
            if f == "n.a.":
                field_formatter[i] = "%s"

        # print dots to show the level of the term
        dots = ""
        if self.goterm is not None and indent:
            dots = "." * self.goterm.level

        return dots + "\t".join(a % b for (a, b) in
                                zip(field_formatter, field_data))

    def __repr__(self):
        return "GOEnrichmentRecord(%s)" % self.id

    def find_goterm(self, go):
        if self.id in list(go.keys()):
            self.goterm = go[self.id]
            self.description = self.goterm.name

    def update_fields(self, **kwargs):
        for k, v in kwargs.items():
            assert k in self._fields, "invalid field name %s" % k
            self.__setattr__(k, v)

    def update_remaining_fields(self, min_ratio=None):
        study_count, study_n = self.ratio_in_study
        pop_count, pop_n = self.ratio_in_pop
        self.enrichment = 'e' if ((1.0 * study_count / study_n) >
                                  (1.0 * pop_count / pop_n)) else 'p'
        self.is_ratio_different = is_ratio_different(min_ratio, study_count,
                                                     study_n, pop_count, pop_n)


class GOEnrichmentStudy(object):
    """Runs Fisher's exact test, as well as multiple corrections
    """
    def __init__(self, pop, assoc, obo_dag, alpha=.05, study=None,
                 methods=["bonferroni", "sidak", "holm"]):

        self.pop = pop
        self.assoc = assoc
        self.obo_dag = obo_dag
        self.alpha = alpha
        self.methods = methods
        self.results = []

        obo_dag.update_association(assoc)
        self.term_pop = count_terms(pop, assoc, obo_dag)

        if study:
            self.run_study(study)

    def run_study(self, study):
        results = self.results

        term_study = count_terms(study, self.assoc, self.obo_dag)

        pop_n, study_n = len(self.pop), len(study)

        for term, study_count in list(term_study.items()):
            pop_count = self.term_pop[term]
            p = fisher.pvalue_population(study_count, study_n,
                                         pop_count, pop_n)

            one_record = GOEnrichmentRecord(
                id=term,
                p_uncorrected=p.two_tail,
                ratio_in_study=(study_count, study_n),
                ratio_in_pop=(pop_count, pop_n))

            results.append(one_record)

        # Calculate multiple corrections
        pvals = [r.p_uncorrected for r in results]
        all_methods = ("bonferroni", "sidak", "holm", "fdr")
        bonferroni, sidak, holm, fdr = None, None, None, None

        for method in self.methods:
            if method == "bonferroni":
                bonferroni = Bonferroni(pvals, self.alpha).corrected_pvals
            elif method == "sidak":
                sidak = Sidak(pvals, self.alpha).corrected_pvals
            elif method == "holm":
                holm = HolmBonferroni(pvals, self.alpha).corrected_pvals
            elif method == "fdr":
                # get the empirical p-value distributions for FDR
                p_val_distribution = calc_qval(study_count, study_n,
                                               pop_count, pop_n,
                                               self.pop, self.assoc,
                                               self.term_pop, self.obo_dag)
                fdr = FDR(p_val_distribution,
                          results, self.alpha).corrected_pvals
            else:
                raise Exception("multiple test correction methods must be "
                                "one of %s" % all_methods)

        all_corrections = (bonferroni, sidak, holm, fdr)

        for method, corrected_pvals in zip(all_methods, all_corrections):
            self.update_results(method, corrected_pvals)

        results.sort(key=lambda r: r.p_uncorrected)
        self.results = results

        for rec in results:
            # get go term for description and level
            rec.find_goterm(self.obo_dag)

        return results

    def update_results(self, method, corrected_pvals):
        if corrected_pvals is None:
            return
        for rec, val in zip(self.results, corrected_pvals):
            rec.__setattr__("p_"+method, val)

    def print_summary(self, min_ratio=None, indent=False, pval=0.05):
        # field names for output
        print("\t".join(GOEnrichmentRecord()._fields))

        for rec in self.results:
            # calculate some additional statistics
            # (over_under, is_ratio_different)
            rec.update_remaining_fields(min_ratio=min_ratio)

            if pval is not None and rec.p_bonferroni > pval:
                continue

            if rec.is_ratio_different:
                print(rec.__str__(indent=indent))


def count_terms(geneset, assoc, obo_dag):
    """count the number of terms in the study group
    """
    term_cnt = collections.defaultdict(int)
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
