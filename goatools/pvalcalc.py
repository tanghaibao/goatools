"""Options for calculating uncorrected p-values."""

from __future__ import print_function

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang et al., All rights reserved."
__author__ = "DV Klopfenstein"

import collections as cx
import sys

class PvalCalcBase(object):
    """Base class for initial p-value calculations."""

    def __init__(self, name, pval_fnc, log):
        self.log = log
        self.name = name
        self.pval_fnc = pval_fnc

    def calc_pvalue(self, study_count, study_n, pop_count, pop_n):
        """pvalues are calculated in derived classes."""
        fnc_call = "calc_pvalue({SCNT}, {STOT}, {PCNT} {PTOT})".format(
            SCNT=study_count, STOT=study_n, PCNT=pop_count, PTOT=pop_n)
        raise Exception("NOT IMPLEMENTED: {FNC_CALL} using {FNC}.".format(
            FNC_CALL=fnc_call, FNC=self.pval_fnc))


class FisherClass(PvalCalcBase):
    """From the 'fisher' package, use function, pvalue_population."""

    def __init__(self, name, log):
        import fisher
        super(FisherClass, self).__init__(name, fisher.pvalue_population, log)

    def calc_pvalue(self, study_count, study_n, pop_count, pop_n):
        """Calculate uncorrected p-values."""
        # k, n = study_true, study_tot,
        # K, N = population_true, population_tot
        # def pvalue_population(int k, int n, int K, int N): ...
        return self.pval_fnc(study_count, study_n, pop_count, pop_n).two_tail


class FisherScipyStats(PvalCalcBase):
    """From the scipy stats package, use function, fisher_exact."""

    fmterr = "STUDY={A}/{B} POP={C}/{D} scnt({scnt}) stot({stot}) pcnt({pcnt}) ptot({ptot})"

    def __init__(self, name, log):
        from scipy import stats
        super(FisherScipyStats, self).__init__(name, stats.fisher_exact, log)

    def calc_pvalue(self, study_count, study_n, pop_count, pop_n):
        """Calculate uncorrected p-values."""
        # http://docs.scipy.org/doc/scipy-0.17.0/reference/generated/scipy.stats.fisher_exact.html
        #
        #         Atlantic  Indian                              YES       NO
        # whales     8        2    | 10 whales    study_genes    8 scnt   2    | 10 = study_n
        # sharks     1        5    |  6 sharks    not s_genes    1        5    |  6
        #         --------  ------                            --------   -----
        #            9        7      16 = pop_n     pop_genes    9 pcnt   7      16 = pop_n
        #
        # We use the preceeding table to find the p-value for whales/sharks:
        #
        # >>> import scipy.stats as stats
        # >>> oddsratio, pvalue = stats.fisher_exact([[8, 2], [1, 5]])
        #                                              a  b    c  d
        avar = study_count
        bvar = study_n - study_count
        cvar = pop_count - study_count
        dvar = pop_n - pop_count - bvar
        assert cvar >= 0, self.fmterr.format(
            A=avar, B=bvar, C=cvar, D=dvar, scnt=study_count, stot=study_n, pcnt=pop_count, ptot=pop_n)
        # stats.fisher_exact returns oddsratio, pval_uncorrected
        _, p_uncorrected = self.pval_fnc([[avar, bvar], [cvar, dvar]])
        return p_uncorrected


class FisherFactory(object):
    """Factory for choosing a fisher function."""

    options = cx.OrderedDict([
        ('fisher', FisherClass),
        ('fisher_scipy_stats', FisherScipyStats),
    ])

    def __init__(self, **kws):
        self.log = kws['log'] if 'log' in kws else sys.stdout
        self.pval_fnc_name = kws["pvalcalc"] if "pvalcalc" in kws else "fisher"
        self.pval_obj = self._init_pval_obj()

    def _init_pval_obj(self):
        """Returns a Fisher object based on user-input."""
        if self.pval_fnc_name in self.options.keys():
            try:
                fisher_obj = self.options[self.pval_fnc_name](self.pval_fnc_name, self.log)
            except ImportError:
                print("fisher module not installed.  Falling back on scipy.stats.fisher_exact")
                fisher_obj = self.options['fisher_scipy_stats']('fisher_scipy_stats', self.log)

            return fisher_obj

        raise Exception("PVALUE FUNCTION({FNC}) NOT FOUND".format(FNC=self.pval_fnc_name))

    def __str__(self):
        return " ".join(self.options.keys())


# Copyright (C) 2016-2018, DV Klopfenstein, H Tang et al., All rights reserved.
