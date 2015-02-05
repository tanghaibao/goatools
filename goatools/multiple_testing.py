#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
A list of commonly used multiple correction routines
"""
from __future__ import print_function
from __future__ import absolute_import
import sys
import random
import fisher
import numpy as np
import goatools.go_enrichment


class AbstractCorrection(object):

    def __init__(self, pvals, a=.05):
        self.pvals = self.corrected_pvals = np.array(pvals)
        self.n = len(self.pvals)    # number of multiple tests
        self.a = a                  # type-1 error cutoff for each test

        self.set_correction()

    def set_correction(self):
        # the purpose of multiple correction is to lower the alpha
        # instead of the canonical value (like .05)
        pass


class Bonferroni(AbstractCorrection):
    """
    >>> Bonferroni([0.01, 0.01, 0.03, 0.05, 0.005], a=0.05).corrected_pvals
    array([ 0.05 ,  0.05 ,  0.15 ,  0.25 ,  0.025])
    """
    def set_correction(self):
        self.corrected_pvals *= self.n


class Sidak(AbstractCorrection):
    """http://en.wikipedia.org/wiki/Bonferroni_correction
    >>> Sidak([0.01, 0.01, 0.03, 0.05, 0.005], a=0.05).corrected_pvals
    array([ 0.04898974,  0.04898974,  0.14696923,  0.24494871,  0.02449487])
    """

    def set_correction(self):
        if self.n != 0:
            correction = self.a * 1. / (1 - (1 - self.a) ** (1. / self.n))
        else:
            correction = 1
        self.corrected_pvals *= correction


class HolmBonferroni(AbstractCorrection):

    """http://en.wikipedia.org/wiki/Holm-Bonferroni_method
    given a list of pvals, perform the Holm-Bonferroni correction
    and return the indexes from original list that are significant.
    (cant use p-value as that may be repeated.)
    >>> HolmBonferroni([0.01, 0.01, 0.03, 0.05, 0.005], a=0.05).corrected_pvals
    array([ 0.04 ,  0.04 ,  0.06 ,  0.05 ,  0.025])
    """
    def set_correction(self):
        if len(self.pvals):
            idxs, correction = list(zip(*self.generate_significant()))
            idxs = list(idxs)
            self.corrected_pvals[idxs] *= correction

    def generate_significant(self):

        pvals = self.pvals
        pvals_idxs = list(zip(pvals, list(range(len(pvals)))))
        pvals_idxs.sort()

        lp = len(self.pvals)

        from itertools import groupby
        for pval, idxs in groupby(pvals_idxs, lambda x: x[0]):
            idxs = list(idxs)
            for p, i in idxs:
                if p * 1. / lp < self.a:
                    yield (i, lp)
            lp -= len(idxs)


class FDR(object):
    def __init__(self, p_val_distribution, results, a=.05):
        self.corrected_pvals = fdr = []
        for rec in results:
            q = (sum(1 for x in p_val_distribution if x < rec.p_uncorrected)
                 * 1.0 / len(p_val_distribution))
            fdr.append(q)


"""
Generate a p-value distribution based on re-sampling, as described in:
http://www.biomedcentral.com/1471-2105/6/168
"""


def calc_qval(study_count, study_n, pop_count, pop_n,
              pop, assoc, term_pop, obo_dag):
    print(("generating p-value distribution for FDR "
                         "calculation (this might take a while)"), file=sys.stderr)
    T = 1000    # number of samples
    distribution = []
    for i in range(T):
        new_study = random.sample(pop, study_n)
        new_term_study = go_enrichment.count_terms(new_study, assoc, obo_dag)

        smallest_p = 1
        for term, study_count in list(new_term_study.items()):
            pop_count = term_pop[term]
            p = fisher.pvalue_population(study_count,
                                         study_n,
                                         pop_count,
                                         pop_n)
            if p.two_tail < smallest_p:
                smallest_p = p.two_tail

        distribution.append(smallest_p)
        print(i, smallest_p, file=sys.stderr)
    return distribution


if __name__ == '__main__':
    import doctest
    doctest.testmod()
