#!/usr/bin/env python
# -*- coding: UTF-8 -*-
from __future__ import print_function

"""
python {} study.file population.file gene-association.file

This program returns P-values for functional enrichment in a cluster of study
genes using Fisher's exact test, and corrected for multiple testing (including
Bonferroni, Holm, Sidak, and false discovery rate).

About significance cutoff:
--alpha: test-wise alpha; for each GO term, what significance level to apply
        (most often you don't need to change this other than 0.05 or 0.01)
--pval: experiment-wise alpha; for the entire experiment, what significance
        level to apply after Bonferroni correction
""".format(__file__)

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
        print("removed %d overlapping items" % (len(common), ), file=sys.stderr)
        print("Set 1: {0}, Set 2: {1}".\
            format(len(study), len(pop)), file=sys.stderr)

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


def check_input_files(ns, p):
    """check filename args. otherwise if one of the 3 filenames is bad
    it's hard to tell which one"""
    import os
    if not len(ns.filenames) == 3:
        p.print_help()
        msg = """
  3 Expected files; Expected content: study population association",
  {} Actual   files: {}""".format(len(ns.filenames), ' '.join(ns.filenames))
        raise Exception(msg)
    for fin in ns.filenames:
        if not os.path.exists(fin):
            return "*{}* does not exist".format(fin)

    return False


if __name__ == "__main__":

    import argparse
    p = argparse.ArgumentParser(__doc__,
                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    p.add_argument('filenames', type=str, nargs='+',
                 help='data/study data/population data/association')
    p.add_argument('--alpha', default=0.05, type=float,
                 help="Test-wise alpha for multiple testing ")
    p.add_argument('--pval', default=None, type=float,
                 help="Family-wise alpha (whole experiment), only print out "
                 "Bonferroni p-value is less than this value. ")
    p.add_argument('--compare', dest='compare', default=False,
                 action='store_true',
                 help="the population file as a comparison group. if this "
                 "flag is specified, the population is used as the study "
                 "plus the `population/comparison`")
    p.add_argument('--ratio', dest='ratio', type=float, default=None,
                 help="only show values where the difference between study "
                 "and population ratios is greater than this. useful for "
                 "excluding GO categories with small differences, but "
                 "containing large numbers of genes. should be a value "
                 "between 1 and 2. ")
    p.add_argument('--fdr', dest='fdr', default=False,
                 action='store_true',
                 help="Calculate the false discovery rate (alt. to the "
                 "Bonferroni but slower)")
    p.add_argument('--indent', dest='indent', default=False,
                 action='store_true', help="indent GO terms")
    p.add_argument('--obo', default="go-basic.obo", type=str,
                 help="Specifies location and name of the obo file")

    args = p.parse_args()
    check_input_files(args, p)

    min_ratio = args.ratio
    if min_ratio is not None:
        assert 1 <= min_ratio <= 2

    assert 0 < args.alpha < 1, "Test-wise alpha must fall between (0, 1)"

    study_fn, pop_fn, assoc_fn = args.filenames
    study, pop = read_geneset(study_fn, pop_fn, compare=args.compare)
    assoc = read_associations(assoc_fn)

    methods = ["bonferroni", "sidak", "holm"]
    if args.fdr:
        methods.append("fdr")

    obo_dag = GODag(obo_file=args.obo)
    g = GOEnrichmentStudy(pop, assoc, obo_dag, alpha=args.alpha,
                          study=study, methods=methods)
    g.print_summary(min_ratio=min_ratio, indent=args.indent, pval=args.pval)
