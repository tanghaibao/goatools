#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
A list of commonly used multiple correction routines
"""
from __future__ import print_function
from __future__ import absolute_import
import sys
import random
import numpy as np
from .ratio import count_terms
import collections as cx

__copyright__ = "Copyright (C) 2010-2016, H Tang et al., All rights reserved."
__author__ = "various"

class Methods(object):
    """Class to manage multipletest methods from both local and remote sources."""

    all_methods = [
        ("local", ("bonferroni", "sidak", "holm", "fdr")),
        ("statsmodels", (
             'bonferroni',     # 0) bonferroni one-step correction
             'sidak',          # 1) sidak one-step correction
             'holm-sidak',     # 2) holm-sidak step-down method using Sidak adjustments
             'holm',           # 3) holm step-down method using Bonferroni adjustments
             'simes-hochberg', # 4) simes-hochberg step-up method  (independent)
             'hommel',         # 5) hommel closed method based on Simes tests (non-negative)
             'fdr_bh',         # 6) fdr_bh Benjamini/Hochberg  (non-negative)
             'fdr_by',         # 7) fdr_by Benjamini/Yekutieli (negative)
             'fdr_tsbh',       # 8) fdr_tsbh two stage fdr correction (non-negative)
             'fdr_tsbky',      # 9) fdr_tsbky two stage fdr correction (non-negative)
            )),
 
    ]
    prefixes = {'statsmodels':'sm_'}
    NtMethodInfo = cx.namedtuple("NtMethodInfo", "source method fieldname") 

    def __init__(self, usr_methods=None):
        self.statsmodels_multicomp = None
        if usr_methods is None:
            usr_methods = ['bonferroni']
        self._init_methods(usr_methods)

    def _init_methods(self, usr_methods):
       """From the methods list, set list of methods to be used during GOEA."""
       self.methods = []
       for usr_method in usr_methods:
           self._add_method(usr_method)

    def _add_method(self, method, method_source=None):
        """Determine method source if needed. Add method to list."""
        try:
            if method_source is not None:
                self._add_method_src(method_source, method)
            else:
                self._add_method_nosrc(method)
        except Exception as inst:
            raise Exception("{ERRMSG}".format(ERRMSG=inst))

    def _add_method_nosrc(self, usr_method):
        """Add method source, method, and fieldname to list of methods."""
        for method_source, available_methods in self.all_methods:
            if usr_method in available_methods:
                fieldname = self.get_fldnm_method(usr_method)
                nt = self.NtMethodInfo(method_source, usr_method, fieldname)
                self.methods.append(nt)
                return
        for src, prefix in self.prefixes.items():
          if usr_method.startswith(prefix):
              method_source = src
              method = usr_method[len(prefix):]
              nt = self.NtMethodInfo(method_source, method, usr_method)
              self.methods.append(nt)
              return
        raise self.rpt_invalid_method(usr_method)

    def getmsg_valid_methods(self):
        """Report the valid methods."""
        msg = []
        msg.append("    Available methods:")
        ctr = self._get_method_cnts()
        for midx, (method_source, methods) in enumerate(self.all_methods):
            msg.append("        {SRC}(".format(SRC=method_source))
            for method in methods:
                prefix = self.prefixes.get(method_source, "")
                prefix = prefix if ctr[method] != 1 else ""
                msg.append("            {P}{M}".format(P=prefix, M=method))
            msg.append("        )")
        return "\n".join(msg)

    def rpt_invalid_method(self, usr_method):
        """Report which methods are available."""
        msgerr = "FATAL: UNRECOGNIZED METHOD({M})".format(M=usr_method)
        msg = [msgerr, self.getmsg_valid_methods(), msgerr]
        raise Exception("\n".join(msg))

    def _get_method_cnts(self):
        """Count the number of times a method is seen."""
        ctr = cx.Counter()
        for method_source, methods in self.all_methods:
            for method in methods:
                ctr[method] += 1
        return ctr
 
    def rpt_valid_methods(self):
        """Report valid methods."""
        pass
            
    def _add_method_src(self, method_source, usr_method, fieldname=None):
        """Add method source and method to list of methods."""
        fieldname = self.get_fieldname(usr_method, method_source)
        available_methods = self.all_methods.get(method_source, None)
        if usr_method in available_methods:
            nt = self.NtMethodInfo(method_source, usr_method, fieldname)
            self.methods.append(nt)
        else: raise Exception("ERROR: FIELD({FN}) METHOD_SOURCE({MS}) AND METHOD({M})".format(
          FN=fieldname, MS=method_source, M=usr_method))

    @staticmethod
    def get_fldnm_method(method):
        """Given method and source, return fieldname for method."""
        fieldname = method.replace('-', '_')
        return fieldname

    def get_statsmodels_multipletests(self):
        """Only load statsmodels package if it is used."""
        if self.statsmodels_multicomp is not None:
            return self.statsmodels_multicomp
        from statsmodels.sandbox.stats.multicomp import multipletests
        self.statsmodels_multicomp = multipletests
        return self.statsmodels_multicomp

    def __iter__(self):
        return iter(self.methods)


class AbstractCorrection(object):

    def __init__(self, pvals, a=.05):
        self.pvals = self.corrected_pvals = np.array(pvals)
        self.n = len(self.pvals)    # number of multiple tests
        self.a = a                  # type-1 error cutoff for each test

        self.set_correction()
        # Reset all pvals > 1 to 1
        self.corrected_pvals[self.corrected_pvals > 1] = 1

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

def mcorrection_factory(pvals, alpha, method):
    """Return 'multiple correction' object of requested AbstractCorrection class."""
    correctioncls = globals().get(method, None)
    if correctioncls is not None:
        return correctioncls(pvals, alpha)
    


"""
Generate a p-value distribution based on re-sampling, as described in:
http://www.biomedcentral.com/1471-2105/6/168
"""


def calc_qval(study_n, pop_n,
              pop, assoc, term_pop, obo_dag, T=500):
    import fisher
    print(("Generate p-value distribution for FDR "
           "based on resampling (this might take a while)"), file=sys.stderr)
    distribution = []
    for i in range(T):
        new_study = random.sample(pop, study_n)
        new_term_study = count_terms(new_study, assoc, obo_dag)

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
        if i % 10  == 0:
            print("Sample {0} / {1}: p-value {2}".\
                        format(i, T, smallest_p), file=sys.stderr)
    return distribution


if __name__ == '__main__':
    import doctest
    doctest.testmod()

# Copyright (C) 2010-2016, H Tang et al., All rights reserved.
