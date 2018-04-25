# -*- coding: UTF-8 -*-
"""
python find_enrichment.py study.file population.file gene-association.file

This program returns P-values for functional enrichment in a cluster of study
genes using Fisher's exact test, and corrected for multiple testing (including
Bonferroni, Holm, Sidak, and false discovery rate).

About significance cutoff:
--alpha: test-wise alpha; for each GO term, what significance level to apply
        (most often you don't need to change this other than 0.05 or 0.01)
--pval: experiment-wise alpha; for the entire experiment, what significance
        level to apply after Bonferroni correction
"""

from __future__ import print_function

__copyright__ = "Copyright (C) 2010-2018, H Tang et al. All rights reserved."
__author__ = "various"

import os
import sys
import argparse
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.obo_parser import GODag
from goatools.associations import read_associations
from goatools.multiple_testing import Methods
from goatools.pvalcalc import FisherFactory


def read_geneset(study_fn, pop_fn, compare=False):
    """Open files containing genes. Return study genes and population genes."""
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
        sys.stderr.write("removed %d overlapping items\n" % (len(common)))
        sys.stderr.write("Set 1: {0}, Set 2: {1}\n".format(
            len(study), len(pop)))

    return study, pop


def _check_input_files(nspc, parser):
    """check filename args. otherwise if one of the 3 filenames is bad
    it's hard to tell which one"""
    if not len(nspc.filenames) == 3:
        parser.print_help()
        msg = """
  3 Expected files; Expected content: study population association",
  {} Actual   files: {}""".format(len(nspc.filenames), ' '.join(nspc.filenames))
        raise Exception(msg)
    for fin in nspc.filenames:
        if not os.path.exists(fin):
            return "*{}* does not exist".format(fin)

    return False


def get_arg_parser():
    """Get enrichment arg parser."""

    #pylint: disable=invalid-name
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
    p.add_argument('--method', default="bonferroni,sidak,holm,fdr_bh", type=str,
                   help=Methods().getmsg_valid_methods())
    p.add_argument('--pvalcalc', default="fisher", type=str,
                   help=str(FisherFactory()))
    p.add_argument('--min_overlap', default=0.7, type=float,
                   help="Check that a minimum amount of study genes are in the population")

    if len(sys.argv) == 1:
        sys.exit(not p.print_help())

    args = p.parse_args()  # Namespace object from argparse
    _check_input_files(args, p)
    return args

def rd_files(filenames, compare, prt=sys.stdout):
    """Read files and return study and population."""
    study_fn, pop_fn, assoc_fn = filenames
    assoc = read_associations(assoc_fn)
    study, pop = read_geneset(study_fn, pop_fn, compare=compare)
    if prt:
        prt.write("Study: {0} vs. Population {1}\n".format(len(study), len(pop)))
    return study, pop, assoc

def get_overlap(study, pop):
    """Get he ratio of study genes which are in the population."""
    return float(len(study & pop)) / len(study)

def chk_genes(study, pop, min_overlap):
    """Check gene sets."""
    if len(pop) < len(study):
        exit("\nERROR: The study file contains more elements than the population file. "
             "Please check that the study file is a subset of the population file.\n")
    # check the fraction of genomic ids that overlap between study
    # and population
    overlap = get_overlap(study, pop)
    if overlap < 0.95:
        sys.stderr.write("\nWARNING: only {} fraction of genes/proteins in study are found in "
                         "the population  background.\n\n".format(overlap))
    if overlap <= min_overlap:
        exit("\nERROR: only {} of genes/proteins in the study are found in the "
             "background population. Please check.\n".format(overlap))


def get_objgoea(pop, assoc, args):
    """Run gene ontology enrichment analysis (GOEA)."""
    obo_dag = GODag(obo_file=args.obo)
    methods = args.method.split(",")
    propagate_counts = not args.no_propagate_counts
    return GOEnrichmentStudy(pop, assoc, obo_dag,
                             propagate_counts=propagate_counts,
                             alpha=args.alpha,
                             pvalcalc=args.pvalcalc,
                             methods=methods)

def prt_results(results, objgoea, args):
    """Print GOEA results."""
    if args.outfile is None:
        min_ratio = args.ratio
        if min_ratio is not None:
            assert 1 <= min_ratio <= 2
        objgoea.print_summary(results, min_ratio=min_ratio, indent=args.indent, pval=args.pval)
    else:
        # Users can print to both tab-separated file and xlsx file in one run.
        outfiles = args.outfile.split(",")
        if args.pval is not None:
            # Only print results when uncorrected p-value < this value.A
            num_orig = len(results)
            results = [r for r in results if r.p_uncorrected <= args.pval]
            print("{N:7,} of {M:,} results have uncorrected P-values <= {PVAL}=pval\n".format(
                N=len(results), M=num_orig, PVAL=args.pval))
        for outfile in outfiles:
            if outfile.endswith(".xlsx"):
                objgoea.wr_xlsx(outfile, results, indent=args.indent)
            else:
                objgoea.wr_tsv(outfile, results, indent=args.indent)

# Copyright (C) 2010-2018, H Tang et al. All rights reserved.
