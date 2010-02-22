#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
python genemerge.py study.file population.file gene-association.file 

This program returns P-values for functional enrichment in a cluster of study genes using Fisher's exact test, and corrected for multiple testing (including Bonferroni, Holm, Sidak, and false discovery rate)
"""

import sys
import collections
import random

import fisher as f
from multiple_testing import Bonferroni, Sidak, HolmBonferroni 
from obo_parser import GODag 


class GOEnrichmentRecord(object):
    """Represents one result (from a single GOTerm) in the GOEnrichmentStudy
    """
    _fields = "id enrichment description ratio_in_study ratio_in_pop"\
            " p_uncorrected p_bonferroni p_holm p_sidak p_fdr".split()

    def __init__(self, **kwargs):
        for f in self._fields:
            self.__setattr__(f, "n.a.")

        for k, v in kwargs.iteritems():
            assert k in self._fields, "invalid field name %s" % k
            self.__setattr__(k, v)

        self.goterm = None # the reference to the GOTerm
    
    def __setattr__(self, name, value):
        self.__dict__[name] = value

    def __str__(self, indent=False):
        field_data = [self.__dict__[f] for f in self._fields]
        field_formatter = ["%s"] * 3 + ["%d/%d"] * 2 + ["%.3g"] * 5
        assert len(field_data)==len(field_formatter)
        
        # default formatting only works for non-"n.a" data
        for i, f in enumerate(field_data):
            if f=="n.a.":
                field_formatter[i] = "%s"

        # print dots to show the level of the term
        dots = ""
        if self.goterm is not None and indent:
            dots = "." * self.goterm.level

        return dots + "\t".join(a % b for (a, b) in \
                zip(field_formatter, field_data))

    def __repr__(self):
        return "GOEnrichmentRecord(%s)" % self.id

    def find_goterm(self, go):
        if self.id in go.keys():
            self.goterm = go[self.id]
            self.description = self.goterm.name

    def update_fields(self, **kwargs):
        for k, v in kwargs.iteritems():
            assert k in self._fields, "invalid field name %s" % k
            self.__setattr__(k, v)

    def update_remaining_fields(self, min_ratio=None):
        study_count, study_n = self.ratio_in_study
        pop_count, pop_n = self.ratio_in_pop
        self.enrichment = 'e' if 1.0* study_count/study_n > 1.0 * pop_count / pop_n else 'p'
        self.is_ratio_different = is_ratio_different(min_ratio, study_count, study_n, pop_count, pop_n)


class GOEnrichmentStudy(object):
    """Runs Fisher's exact test, as well as multiple corrections
    """
    def __init__(self, study_set, population_set, associations, alpha=None, \
            methods=["bonferroni", "sidak", "holm"]):
        self.results = results = []

        term_study = count_terms(study, assoc)
        term_pop = count_terms(pop, assoc)

        pop_n, study_n = len(pop), len(study)

        for term, study_count in term_study.items():
            pop_count = term_pop[term]
            left_p, right_p, p_val = f.pvalue(study_count, study_n, pop_count, pop_n)

            one_record = GOEnrichmentRecord(id=term, p_uncorrected=p_val,\
                    ratio_in_study=(study_count, study_n), 
                    ratio_in_pop=(pop_count, pop_n))

            results.append(one_record)

        # Calculate multiple corrections
        pvals = [r.p_uncorrected for r in results]
        all_methods = ("bonferroni", "sidak", "holm", "fdr")
        bonferroni, sidak, holm, fdr = None, None, None, None

        for method in methods:
            if method=="bonferroni":
                bonferroni = Bonferroni(pvals, alpha).corrected_pvals
            elif method=="sidak":
                sidak = Sidak(pvals, alpha).corrected_pvals
            elif method=="holm":
                holm = HolmBonferroni(pvals, alpha).corrected_pvals
            elif method=="fdr":
                # get the empirical p-value distributions for FDR
                print >>sys.stderr, "generating p-value distribution for FDR calculation " \
                    "(this might take a while)"
                p_val_distribution = calc_qval(study_count, study_n, pop_count, pop_n, \
                    pop, assoc, term_pop)
                fdr = []
                for rec in results:
                    q = sum(1 for x in p_val_distribution if x < rec.p_uncorrected) \
                            * 1./len(p_val_distribution)
                    fdr.append(q)
            else:
                raise Exception, "multiple test correction methods must be one of %s" % all_methods

        all_corrections = (bonferroni, sidak, holm, fdr)

        for method, corrected_pvals in zip(all_methods, all_corrections):
            self.update_results(method, corrected_pvals)

        results = list(self.filter_results(alpha))
        results.sort(key=lambda r: r.p_uncorrected) 
        self.results = results

    def update_results(self, method, corrected_pvals):
        if corrected_pvals is None: return
        for rec, val in zip(self.results, corrected_pvals):
            rec.__setattr__("p_"+method, val)

    def filter_results(self, alpha=None):
        if alpha is not None: alpha = float(alpha)
        for rec in self.results:
            if alpha is None: yield rec 
            elif rec.p_bonferroni < alpha: yield rec

    def populate_go(self, obo_file="gene_ontology.1_2.obo"):
        go = GODag(obo_file)
        for rec in self.results:
            # get go term for description and level
            rec.find_goterm(go)

    def print_summary(self, min_ratio=None):
        # field names for output
        print "\t".join(GOEnrichmentRecord()._fields)

        for rec in self.results:
            # calculate some additional statistics (over_under, is_ratio_different)
            rec.update_remaining_fields(min_ratio=min_ratio)

            if rec.is_ratio_different:
                print rec.__str__(indent=opts.indent)

"""
Generate a p-value distribution based on re-sampling, as described in:
http://www.biomedcentral.com/1471-2105/6/168
"""
#class FalseDiscoveryRate(AbstractCorrection):
def calc_qval(study_count, study_n, pop_count, pop_n, pop, assoc, term_pop):
    T = 1000 # number of samples
    distribution = []
    for i in xrange(T):
        new_study = random.sample(pop, study_n)
        new_term_study = count_terms(new_study, assoc)

        smallest_p = 1
        for term, study_count in new_term_study.items():
            pop_count = term_pop[term]
            left_p, right_p, p_val = f.pvalue(study_count, study_n, pop_count, pop_n)
            if p_val < smallest_p: smallest_p = p_val

        distribution.append(smallest_p)
        print >>sys.stderr, i, smallest_p
    return distribution


def read_geneset(study_fn, pop_fn, compare=False):
    pop = set(_.strip() for _ in open(pop_fn) if _.strip())
    study = frozenset(_.strip() for _ in open(study_fn) if _.strip())
    # some times the pop is a second group to compare, rather than the population
    # in that case, we need to make sure the overlapping terms are removed first
    if compare:
        common = pop & study
        pop |= study
        pop -= common
        study -= common
        print >>sys.stderr, "removed %d overlapping items" % (len(common), )

    return study, pop


def read_associations(assoc_fn):
    assoc = {}
    sep = " "
    for row in open(assoc_fn):
        if len(row.strip().split())<2: continue 
        try:
            # the accn may have a space. in which case get > 2 tokens.
            a, b = row.split(sep)
        except ValueError:
            sep = "\t"
            a, b = row.split(sep)
        b = set(b.replace(";"," ").split())
        assoc[a] = b

    return assoc


def check_bad_args(args):
    """check args. otherwise if one of the 3 args is bad
    it's hard to tell which one"""
    import os
    if not len(args) == 3: return "please send in 3 file names"
    for arg in args[:-1]:
        if not os.path.exists(arg):
            return "*%s* does not exist" % arg

    return False


def count_terms(geneset, assoc):
    """count the number of terms in the study group
    """
    term_cnt = collections.defaultdict(int)
    for gene in (g for g in geneset if g in assoc):
        for x in assoc[gene]:
            term_cnt[x] += 1

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



if __name__ == "__main__":

    import optparse
    p = optparse.OptionParser(__doc__)

    p.add_option('--alpha', dest='alpha', default=None, 
                 help="only print out the terms where the corrected p-value"
                 " is less than this value. [default: %default]")
    p.add_option('--compare', dest='compare', default=False, action='store_true',
                 help="the population file as a comparison group. if this flag is specified,"
                 " the population is used as the study plus the `population/comparison`")
    p.add_option('--ratio', dest='ratio', type='float', default=None,
                 help="only show values where the difference between study and population"
                 " ratios is greater than this. useful for excluding GO categories with"
                " small differences, but containing large numbers of genes. should be a "
                " value between 1 and 2. ")
    p.add_option('--fdr', dest='fdr', default=False,
                action='store_true',
                help="calculate the false discovery rate (alternative to the Bonferroni correction)")
    p.add_option('--indent', dest='indent', default=False,
                action='store_true', help="indent GO terms")

    (opts, args) = p.parse_args()
    bad = check_bad_args(args)
    if bad:
        print bad
        sys.exit(p.print_help())

    alpha = float(opts.alpha) if opts.alpha else 0.05

    min_ratio = opts.ratio
    if not min_ratio is None:
        assert 1 <= min_ratio <= 2

    study_fn, pop_fn, assoc_fn = args
    study, pop = read_geneset(study_fn, pop_fn, compare=opts.compare)
    assoc = read_associations(assoc_fn)

    methods=["bonferroni", "sidak", "holm"]
    if opts.fdr:
        methods.append("fdr")

    g = GOEnrichmentStudy(study, pop, assoc, alpha=alpha, methods=methods)
    g.populate_go(obo_file="gene_ontology.1_2.obo")
    g.print_summary(min_ratio)

