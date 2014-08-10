#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
python %prog study.file population.file gene-association.file

This program returns P-values for functional enrichment in a cluster of study
genes using Fisher's exact test, and corrected for multiple testing (including
Bonferroni, Holm, Sidak, and false discovery rate).

About significance cutoff:
--alpha: test-wise alpha; for each GO term, what significance level to apply
        (most often you don't need to change this other than 0.05 or 0.01)
--pval: experiment-wise alpha; for the entire experiment, what significance
        level to apply after Bonferroni correction
"""

import sys
import os.path as op
sys.path.insert(0, op.join(op.dirname(__file__), ".."))
from goatools import GOEnrichmentStudy
from goatools.obo_parser import GODag


def read_geneset(study_fn, pop_fn, compare=False):
    pop = set(_.strip() for _ in open(pop_fn) if _.strip())
    study = frozenset(_.strip() for _ in open(study_fn) if _.strip())
    # some times the pop is a second group to compare, rather than the
    # population in that case, we need to make sure the overlapping terms
    # are removed first
    if compare:
        common = pop & study
        pop |= study
        pop -= common
        study -= common
        print >>sys.stderr, "removed %d overlapping items" % (len(common), )
        print >>sys.stderr, "Set 1: {0}, Set 2: {1}".\
            format(len(study), len(pop))

    return study, pop


def read_associations(assoc_fn):
    assoc = {}
    for row in open(assoc_fn):
        atoms = row.split()
        if len(atoms) == 2:
            a, b = atoms
        elif len(atoms) > 2 and row.count('\t') == 1:
            a, b = row.split("\t")
        else:
            continue
        b = set(b.split(";"))
        assoc[a] = b

    return assoc


def check_bad_args(args):
    """check args. otherwise if one of the 3 args is bad
    it's hard to tell which one"""
    import os
    if not len(args) == 3:
        return "please send in 3 file names"
    for arg in args[:-1]:
        if not os.path.exists(arg):
            return "*%s* does not exist" % arg

    return False


if __name__ == "__main__":

    import optparse
    p = optparse.OptionParser(__doc__)

    p.add_option('--alpha', default=0.05, type="float",
                 help="Test-wise alpha for multiple testing "
                 "[default: %default]")
    p.add_option('--pval', default=None, type="float",
                 help="Family-wise alpha (whole experiment), only print out "
                 "Bonferroni p-value is less than this value. "
                 "[default: %default]")
    p.add_option('--compare', dest='compare', default=False,
                 action='store_true',
                 help="the population file as a comparison group. if this "
                 "flag is specified, the population is used as the study "
                 "plus the `population/comparison`")
    p.add_option('--ratio', dest='ratio', type='float', default=None,
                 help="only show values where the difference between study "
                 "and population ratios is greater than this. useful for "
                 "excluding GO categories with small differences, but "
                 "containing large numbers of genes. should be a value "
                 "between 1 and 2. ")
    p.add_option('--fdr', dest='fdr', default=False,
                 action='store_true',
                 help="Calculate the false discovery rate (alt. to the "
                 "Bonferroni but slower)")
    p.add_option('--indent', dest='indent', default=False,
                 action='store_true', help="indent GO terms")

    (opts, args) = p.parse_args()
    bad = check_bad_args(args)
    if bad:
        print bad
        sys.exit(p.print_help())

    min_ratio = opts.ratio
    if min_ratio is not None:
        assert 1 <= min_ratio <= 2

    assert 0 < opts.alpha < 1, "Test-wise alpha must fall between (0, 1)"

    study_fn, pop_fn, assoc_fn = args
    study, pop = read_geneset(study_fn, pop_fn, compare=opts.compare)
    assoc = read_associations(assoc_fn)

    methods = ["bonferroni", "sidak", "holm"]
    if opts.fdr:
        methods.append("fdr")

    obo_dag = GODag(obo_file="gene_ontology.1_2.obo")
    g = GOEnrichmentStudy(pop, assoc, obo_dag, alpha=opts.alpha,
                          study=study, methods=methods)
    g.print_summary(min_ratio=min_ratio, indent=opts.indent, pval=opts.pval)
