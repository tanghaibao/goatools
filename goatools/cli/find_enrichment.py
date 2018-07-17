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

from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.associations import read_associations
from goatools.multiple_testing import Methods
from goatools.pvalcalc import FisherFactory

from goatools.gosubdag.gosubdag import GoSubDag
from goatools.grouper.grprdflts import GrouperDflts
from goatools.grouper.hdrgos import HdrgosSections
from goatools.grouper.aart_geneproducts_all import AArtGeneProductSetsAll
from goatools.grouper.read_goids import read_sections


# pylint: disable=too-few-public-methods
class GoeaCliArgs(object):
    """Extracts arguments from the command-line."""

    def __init__(self):
        self.args = self._init_args()

    def _init_args(self):
        """Get enrichment arg parser."""

        #pylint: disable=invalid-name
        p = argparse.ArgumentParser(__doc__,
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

        p.add_argument('filenames', type=str, nargs=3,
                       help='data/study data/population data/association')
        p.add_argument('--alpha', default=0.05, type=float,
                       help='Test-wise alpha for multiple testing')
        p.add_argument('--pval', default=.05, type=float,
                       help='Only print results with uncorrected p-value < PVAL.')
        p.add_argument('--pval_field', default='uncorrected', type=str,
                       help='Only print results when PVAL_FIELD < PVAL.')
        p.add_argument('--outfile', default=None, type=str,
                       help='Write enrichment results into xlsx or tsv file')
        p.add_argument('--sections', default=None, type=str,
                       help=('Use sections file for printing grouped GOEA results. '
                             'Example SECTIONS values:\n'
                             'goatools.test_data.sections.gjoneska_pfenning \n'
                             'goatools/test_data/sections/gjoneska_pfenning.py \n'
                             'data/gjoneska_pfenning/sections_in.txt\n'))
        p.add_argument('--outfile_detail', type=str,
                       help=('Write enrichment results into a text file \n'
                             'containing the following information: \n'
                             '1) GOEA GO terms, grouped into sections \n\n'
                             '2) List of genes and ASCII art showing section membership \n'
                             '3) Detailed list of each gene and GO terms w/their P-values \n'))
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
        p.add_argument('--method', default="bonferroni,sidak,holm,fdr_bh", type=str,
                       help=Methods().getmsg_valid_methods())
        p.add_argument('--pvalcalc', default="fisher", type=str,
                       help=str(FisherFactory()))
        p.add_argument('--min_overlap', default=0.7, type=float,
                       help="Check that a minimum amount of study genes are in the population")
        p.add_argument('--goslim', default='goslim_generic.obo', type=str,
                       help="The GO slim file is used when grouping GO terms.")

        if len(sys.argv) == 1:
            sys.exit(not p.print_help())

        args = p.parse_args()  # Namespace object from argparse
        self._check_input_files(args, p)
        return args

    @staticmethod
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


class GoeaCliFnc(object):
    """For running a GOEA on the command-line."""

    def __init__(self, args):
        self.args = args
        self.sections = read_sections(self.args.sections) if self.args.sections else None
        _optional_attrs = ['relationship'] if self.sections else None
        self.godag = GODag(obo_file=self.args.obo, optional_attrs=_optional_attrs)
        # Get GOEnrichmentStudy
        _study, _pop, _assoc = self.rd_files()
        if not self.args.compare:  # sanity check
            self.chk_genes(_study, _pop)
        self.methods = self.args.method.split(",")
        self.objgoea = self._init_objgoea(_pop, _assoc)
        # Run GOEA
        self.results_all = self.objgoea.run_study(_study)

    def prt_results(self, results):
        """Print GOEA results."""
        if self.args.outfile is None:
            min_ratio = self.args.ratio
            if min_ratio is not None:
                assert 1 <= min_ratio <= 2
            self.objgoea.print_summary(
                results, min_ratio=min_ratio, indent=self.args.indent, pval=self.args.pval)
        else:
            # Users can print to both tab-separated file and xlsx file in one run.
            outfiles = self.args.outfile.split(",")
            for outfile in outfiles:
                if outfile.endswith(".xlsx"):
                    self.objgoea.wr_xlsx(outfile, results, indent=self.args.indent)
                else:
                    self.objgoea.wr_tsv(outfile, results, indent=self.args.indent)

    def get_results(self):
        """Given all GOEA results, return the significant results (< pval)."""
        return self.get_results_sig() if self.args.pval is not None else self.results_all

    def _init_objgoea(self, pop, assoc):
        """Run gene ontology enrichment analysis (GOEA)."""
        propagate_counts = not self.args.no_propagate_counts
        return GOEnrichmentStudy(pop, assoc, self.godag,
                                 propagate_counts=propagate_counts,
                                 alpha=self.args.alpha,
                                 pvalcalc=self.args.pvalcalc,
                                 methods=self.methods)

    def chk_genes(self, study, pop):
        """Check gene sets."""
        if len(pop) < len(study):
            exit("\nERROR: The study file contains more elements than the population file. "
                 "Please check that the study file is a subset of the population file.\n")
        # check the fraction of genomic ids that overlap between study
        # and population
        overlap = self.get_overlap(study, pop)
        if overlap < 0.95:
            sys.stderr.write("\nWARNING: only {} fraction of genes/proteins in study are found in "
                             "the population  background.\n\n".format(overlap))
        if overlap <= self.args.min_overlap:
            exit("\nERROR: only {} of genes/proteins in the study are found in the "
                 "background population. Please check.\n".format(overlap))

    def get_results_sig(self):
        """Get significant results."""
        # Only print results when uncorrected p-value < this value.
        pval_fld = self._get_pval_field()
        results = [r for r in self.results_all if getattr(r, pval_fld) <= self.args.pval]
        print("{N:7,} of {M:,} results have uncorrected P-values <= {PVAL}=pval\n".format(
            N=len(results), M=len(self.results_all), PVAL=self.args.pval))
        return results

    @staticmethod
    def get_overlap(study, pop):
        """Get he ratio of study genes which are in the population."""
        return float(len(study & pop)) / len(study)

    def _get_pval_field(self):
        """Get 'p_uncorrected' or the user-specified field for determining significant results."""
        pval_fld = self.args.pval_field
        if pval_fld[:2] != 'p_':
            pval_fld = 'p_' + pval_fld
        if self.results_all:
            assert hasattr(next(iter(self.results_all)), pval_fld), \
                'NO PVAL({P}). EXPECTED ONE OF: {E}'.format(
                    P=self.args.pval_field,
                    E=" ".join([k for k in dir(next(iter(self.results_all))) if k[:2] == 'p_']))
        return pval_fld

    def rd_files(self):
        """Read files and return study and population."""
        study_fn, pop_fn, assoc_fn = self.args.filenames
        assoc = read_associations(assoc_fn)
        study, pop = self._read_geneset(study_fn, pop_fn)
        print("Study: {0} vs. Population {1}\n".format(len(study), len(pop)))
        return study, pop, assoc

    def _read_geneset(self, study_fn, pop_fn):
        """Open files containing genes. Return study genes and population genes."""
        pop = set(_.strip() for _ in open(pop_fn) if _.strip())
        study = frozenset(_.strip() for _ in open(study_fn) if _.strip())
        # some times the pop is a second group to compare, rather than the
        # population in that case, we need to make sure the overlapping terms
        # are removed first
        if self.args.compare:
            common = pop & study
            pop |= study
            pop -= common
            study -= common
            sys.stderr.write("removed %d overlapping items\n" % (len(common)))
            sys.stderr.write("Set 1: {0}, Set 2: {1}\n".format(
                len(study), len(pop)))
        return study, pop

    def get_objaart(self):
        """Get background database info for making ASCII art."""
        print("AAAAAAAAAAAAAAAAAAAAAaa")
        goids = set(o.id for o in self.godag.values() if not o.children)
        # gosubanno = _init_gosubanno(10090)
        # _tobj = TermCounts(godag, gosubanno.gene2gos)
        gosubdag = GoSubDag(goids, self.godag, relationships=True) # , tcntobj=_tobj, prt=log)
        grprdflt = GrouperDflts(gosubdag, self.args.goslim)
        hdrobj = HdrgosSections(gosubdag, grprdflt.hdrgos_dflt, self.sections)
        kws = {
            'sortgo':lambda nt: [nt.NS, nt.dcnt],
            # fmtgo=('{p_fdr_bh:8.2e} {GO} '
            # Formatting for GO terms in grouped GO list
            'fmtgo':('{hdr1usr01:2} {NS} {GO} {s_fdr_bh:8} '
                     '{dcnt:5} {childcnt:3} R{reldepth:02} '
                     '{D1:5} {GO_name} ({study_count} study genes)\n'),
            # Formatting for GO terms listed under each gene
            'fmtgo2':('{hdr1usr01:2} {NS} {GO} {s_fdr_bh:8} '
                      '{dcnt:5} R{reldepth:02} '
                      '{GO_name} ({study_count} study genes)\n'),
            # itemid2name=ensmusg2symbol}
            }
        return AArtGeneProductSetsAll(grprdflt, hdrobj, **kws)




# Copyright (C) 2010-2018, H Tang et al. All rights reserved.
