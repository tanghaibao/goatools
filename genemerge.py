#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
python genemerge.py gene-association.file description.file population.file study.file output

Results are written to stdout if output is `stdout`

This program returns P-values for functional enrichment in a cluster of study genes using fisher's exact test, and corrected for multiple testing (including Bonferroni, Holm, Sidak, and false discovery rate)
"""

import sys
import operator
import collections
import random

import fisher as f
from obo_parser import load_godag

#http://en.wikipedia.org/wiki/Bonferroni_correction
sidakp = lambda n, a: 1 - (1 - a) ** (1. /n )

def holm_bonferroni(pvals, a=0.05):
    """http://en.wikipedia.org/wiki/Holm-Bonferroni_method
    given a list of pvals, perform the Holm-Bonferroni correction
    and return the indexes from original list that are significant.
    (cant use p-value as that may be repeated.)
    >>> list(holm_bonferroni([0.01, 0.01, 0.03, 0.05, 0.005], a=0.05))
    [4, 0, 1]

    TODO: account for weights by number in study?
    """
    idxs = collections.defaultdict(list)
    for i, p in enumerate(pvals):
        idxs[p].append(i)

    lp = len(pvals)
    if lp == 0:
        raise StopIteration
    an = a / lp
    for p in sorted(pvals):
        if p < an:
            yield idxs[p].pop(0) # [1, 2, 3] becomes [2, 3] and yields 1
            lp -= 1
            an = a / lp
        else:
            break


def calc_qval(study_count, study_n, pop_count, pop_n, pop, assoc, term_cnt):
    """
    Generate a p-value distribution based on re-sampling, as described in:
    http://www.biomedcentral.com/1471-2105/6/168
    """
    T = 1000 # number of samples
    distribution = []
    for i in xrange(T):
        new_study = random.sample(pop, study_n)
        new_term_study = count_term_study(new_study, assoc)

        smallest_p = 1
        for term, study_count in new_term_study.items():
            pop_count = term_cnt[term]
            left_p, right_p, p_val = f.pvalue(study_count, study_n, pop_count, pop_n)
            if p_val < smallest_p: smallest_p = p_val

        distribution.append(smallest_p)
        print >>sys.stderr, i, smallest_p
    return distribution


def count_associations(assoc_fn):
    term_cnt = collections.defaultdict(int)
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
        if a not in pop: continue
        b = set(b.replace(";"," ").split())
        assoc[a] = b
        for term in b:
            term_cnt[term]+=1

    return assoc, term_cnt


def get_description(desc_fn):
    desc = {}
    for row in open(desc_fn):
        if row.strip()=="": continue 
        a,b = row.strip().split("\t")
        desc[a] = b
    return desc


def filter_results(results, alpha=None):
    if alpha is not None: alpha = float(alpha)
    for row in results:
        if "*" in row[-2:]: yield row
        else:
            if alpha is None: yield row
            elif row[-3] < alpha: yield row


def check_bad_args(args):
    """check args. otherwise if one of the 5 args is bad
    it's hard to tell which one"""
    import os
    if not len(args) == 5: return "please send in 5 file names"
    for arg in args[:-1]:
        if not os.path.exists(arg):
            return "*%s* does not exist" % arg
    return False


def count_term_study(study, assoc):
    """count the number of terms in the study group
    """
    term_study = collections.defaultdict(int)
    for gene in (g for g in study if g in assoc):
        for _ in assoc[gene]:
            term_study[_] += 1
    return term_study


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
    import doctest
    doctest.testmod()

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
    p.add_option('--hierarchy', dest='hierarchy', default=False,
                action='store_true', help="group overlapping GO terms and indent")

    opts, args = p.parse_args()
    bad = check_bad_args(args)
    if bad:
        print bad
        sys.exit(p.print_help())
    alpha = float(opts.alpha) if opts.alpha else 0.05

    min_ratio = opts.ratio
    if not min_ratio is None:
        assert 1 <= min_ratio <= 2

    assoc_fn, desc_fn, pop_fn, study_fn, out_fn = args

    # Calculations start here
    pop = set(_.strip() for _ in open(pop_fn) if _.strip())
    study = frozenset(_.strip() for _ in open(study_fn) if _.strip())
    if opts.compare:
        common = pop.intersection(study)
        pop = pop.union(study)
        pop = pop.difference(common)
        study = study.difference(common)
        print >>sys.stderr, "removed %d overlapping items" % (len(common), )

    assoc, term_cnt = count_associations(assoc_fn)
    desc = get_description(desc_fn)


    pop_n, study_n = len(pop), len(study)
    non_singletons = set(k for k, v in term_cnt.items() if v > 1)

    term_study = count_term_study(study, assoc)

    results = []
    for term, study_count in term_study.items():
        pop_count = term_cnt[term]
        left_p, right_p, p_val = f.pvalue(study_count, study_n, pop_count, pop_n)

        results.append([term, study_count, pop_count, p_val, left_p, right_p])


    # get the indexes of values that are significant by the holm-b test.
    holm_significant = list(holm_bonferroni([r[3] for r in results], a=alpha))

    correction = sum(1 for _ in term_study if _ in non_singletons)
    sidak = sidakp(correction, alpha) if correction > 0 else None

    # get the empirical p-value distributions for FDR
    if opts.fdr:
        print >>sys.stderr, "generating p-value distribution for FDR calculation " \
            "(this might take a while)"
        p_val_distribution = calc_qval(study_count, study_n, pop_count, pop_n, \
            pop, assoc, term_cnt)

    for i, r in enumerate(results):
        pval = float(r[3])
        results[i].append(pval * correction)
        results[i].append("*" if i in holm_significant else ".")
        results[i].append("*" if pval < sidak else ".")
        if opts.fdr:
            qval = sum(1 for x in p_val_distribution if x < pval) \
                    * 1./len(p_val_distribution)
        else:
            qval = 0
        results[i].append(qval)

    results = list(filter_results(results, opts.alpha))
    results.sort(key=operator.itemgetter(5)) # p_raw

    if opts.hierarchy:
        g = load_godag()

    fw = sys.stdout if out_fn=="stdout" else file(out_fn, "w") 
    # header for the output
    fw.write("go\tenriched_purified\tgo_desc\tgo/n_in_study\tgo/n_in_pop\tp_raw\tp_bonferroni\tp_holm\tp_sidak")
    if opts.fdr: fw.write("\tq_value")
    fw.write("\n")

    for term, study_count, pop_count, p_raw, left_p, right_p, p_corrected, p_holm, p_sidak, q_value in results:
        if opts.hierarchy:
            rec = g[term]
            term = "." * rec.level + term + " [%s]" % rec.name 

        D = desc.get(term, "No description")
        if is_ratio_different(min_ratio, study_count, study_n, pop_count, pop_n):
            over_under = 'e' if 1.0* study_count/study_n > 1.0 * pop_count / pop_n else 'p'
            fw.write("%s\t%s\t%s\t%d/%d\t%d/%d\t%.3g\t%.3g\t%s\t%s"%(term, over_under, D, study_count,\
                    study_n, pop_count, pop_n, p_raw, p_corrected, p_holm, p_sidak))
            if opts.fdr: fw.write("\t%.3f" % q_value)
            fw.write("\n")
    fw.close()

    print >>sys.stderr, "results written to %s" % out_fn
