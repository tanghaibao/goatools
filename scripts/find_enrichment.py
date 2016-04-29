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

__copyright__ = "Copyright (C) 2010-2016, H Tang et al., All rights reserved."
__author__ = "various"

import sys
import os.path as op
sys.path.insert(0, op.join(op.dirname(__file__), ".."))
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.obo_parser import GODag
from goatools.associations import read_associations
from goatools.multiple_testing import Methods


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

    p.add_argument('filenames', type=str, nargs=3,
                 help='data/study data/population data/association')
    p.add_argument('--alpha', default=0.05, type=float,
                 help="Test-wise alpha for multiple testing ")
    p.add_argument('--pval', default=.05, type=float,
                 help="Only print out when uncorrected p-value < this value.")
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
    p.add_argument('--indent', dest='indent', default=False,
                 action='store_true', help="indent GO terms")
    p.add_argument('--obo', default="go-basic.obo", type=str,
                 help="Specifies location and name of the obo file")
    p.add_argument('--no_propagate_counts', default=False, action='store_true',
                 help="Do not propagate counts to parent terms")
    p.add_argument('--outfile', default=None, type=str,
                 help="Write enrichment results into xlsx or tsv file")
    p.add_argument('--method', default="bonferroni,sidak,holm", type=str,
                 help=Methods().getmsg_valid_methods())

    if len(sys.argv) == 1:
        sys.exit(not p.print_help())

    args = p.parse_args()
    check_input_files(args, p)

    min_ratio = args.ratio
    if min_ratio is not None:
        assert 1 <= min_ratio <= 2

    study_fn, pop_fn, assoc_fn = args.filenames
    study, pop = read_geneset(study_fn, pop_fn, compare=args.compare)
    print("Study: {0} vs. Population {1}".format(len(study), len(pop)), file=sys.stderr)

    if not args.compare:  # sanity check
        if len(pop) < len(study):
            exit("\nERROR: The study file contains more elements than the population file. "
                 "Please check that the study file is a subset of the population file.\n")
        # check the fraction of genomic ids that overlap between study
        # and population
        overlap = float(len(study & pop)) / len(study)
        if 0.7 < overlap < 0.95:
            sys.stderr.write("\nWARNING: only {} fraction of genes/proteins in study are found in "
                             "the population  background.\n\n".format(overlap))
        if overlap <= 0.7:
            exit("\nERROR: only {} of genes/proteins in the study are found in the "
                 "background population. Please check.\n".format(overlap))

    assoc = read_associations(assoc_fn)

    methods = args.method.split(",")

    obo_dag = GODag(obo_file=args.obo)
    propagate_counts = not args.no_propagate_counts
    g = GOEnrichmentStudy(pop, assoc, obo_dag,
                          propagate_counts=propagate_counts,
                          alpha=args.alpha,
                          methods=methods)
    results = g.run_study(study)
    if args.outfile is None:
        g.print_summary(results, min_ratio=min_ratio, indent=args.indent, pval=args.pval)
    else:
        # Users can print to both tab-separated file and xlsx file in one run.
        outfiles = args.outfile.split(",")
        prt_if = None # Print all values
        if args.pval is not None:
            # Only print out when uncorrected p-value < this value.
            prt_if = lambda nt: nt.p_uncorrected < args.pval
        for outfile in outfiles:
            if outfile.endswith(".xlsx"):
                g.wr_xlsx(outfile, results, prt_if=prt_if)
            else:
                g.wr_tsv(outfile, results, prt_if=prt_if)

# Copyright (C) 2010-2016, H Tang et al., All rights reserved.
